library(tidyverse)
library(lubridate)
library(dplyr)
library(sf)
library(zoo)
#getwd() #current working directory
#setwd() #set working directory file path to folder with all files

#######################
### Reading in data ###
#######################

#groundwater levels:
usgs_data <- read.csv("WATERLEVEL_USGS_24_26.csv", colClasses = c("SiteNo" = "character"))
SJRWMD_data <- read.csv("WATERLEVEL_SJRWMD_24_26.csv", colClasses = c("SiteNo" = "character"))

#combining lake ids to groundwater ids
all_gws <- read.csv("all_sites_withlakeids_1000m4th.csv")
connected_ids <- all_gws[c("noSiteNo", "lake_id", "AgencyCd", "distance", "HorzDatum", "geologyCAT", "lkgeologyCAT", "WellDepth", "AquiferTyp")]
# All well depths are in feet # Aquifer type is either confined our unconfined
connected_ids$WellDepth_m <- connected_ids$WellDepth / 3.281 #converting to meters


#filtering by agency
con_usgs <- connected_ids%>%
  filter(AgencyCd == "USGS") 

#con_fldep <- connected_ids %>%
#  filter(AgencyCd == "FLDEP")
con_sjrwmd <- connected_ids %>%
  filter(AgencyCd == "SJRWMD") 

#SWOT DATA
swot <- read.csv("SWOT_levels.csv")

##########################################################################
#Testing of one gw site and lake to then implement in function below
# 282202081384601   7320211482
# Groundwater site 142m away
##########################################################################
#####################
#### Data Setup #####
#####################

#filtering gw to specific site id
gw_28171 <- usgs_data %>%
  filter(SiteNo == 282202081384601) %>%
  mutate(Time = ymd_hms(Time))

#calculating meters measurement to match SWOT units
gw_28171$m <- gw_28171$Water.level.in.feet.relative.to.NGVD29 / 3.281


#filtering to lake specific id only
swot_73202 <- swot %>% 
  filter(lake_id == 7320211482) %>%
  filter(quality_f <= 1) %>%
  mutate(time_str = ymd_hms(time_str))


#calculation of normalized water levels by subtracting mean
gw_28171$norm <- (gw_28171$m - mean(gw_28171$m, na.rm = TRUE))
swot_73202$norm <- (swot_73202$wse - mean(swot_73202$wse))

error_m <- 0.1 #meters
###################################
### cross correlation lag times ###
###################################

#shared daily timeline
all_dates <- data.frame(Time = seq(min(gw_28171$Time), max(gw_28171$Time), by = "1 day"))

#GW data (average to daily if it's hourly)
gw_daily <- gw_28171 %>%
  mutate(Time = as.Date(Time)) %>%
  group_by(Time) %>%
  summarize(gw_norm = mean(norm, na.rm = TRUE))

#SWOT data (assign to dates)
swot_daily <- swot_73202 %>%
  mutate(Time = as.Date(time_str)) %>%
  group_by(Time) %>%
  summarize(swot_norm = mean(norm, na.rm = TRUE))

#Merge and Interpolate the SWOT gaps
combined_ts <- all_dates %>%
  mutate(Time = as.Date(Time)) %>%
  left_join(gw_daily, by = "Time") %>%
  left_join(swot_daily, by = "Time") %>%
  mutate(
    # Linear interpolation for SWOT because of the gaps between passes
    swot_interp = na.approx(swot_norm, na.rm = FALSE),
    gw_interp = na.approx(gw_norm, na.rm = FALSE)
  ) %>%
  filter(!is.na(swot_interp) & !is.na(gw_interp))


#Run Cross-Correlation
cross_corr <- ccf(combined_ts$gw_interp, combined_ts$swot_interp, 
                  lag.max = 100, #Look for lags up to 100 days
                  plot = FALSE)

# Extract the results into a dataframe for plotting
ccf_df <- data.frame(lag = cross_corr$lag, acf = cross_corr$acf)

# Find the specific lag with the highest correlation
max_lag <- ccf_df[which.max(abs(ccf_df$acf)), ]

##################################################
### Plot of lag correlation coeffecient values ###
##################################################

ggplot(ccf_df, aes(x = lag, y = acf)) +
  geom_col(fill = "steelblue") +
  #geom_hline(color = "red", linetype = "dashed", 
  #           yintercept = c(-1.96/sqrt(nrow(combined_ts)), 1.96/sqrt(nrow(combined_ts)))) +
  labs(title = paste("Lag Correlation: GW vs SWOT"),
       subtitle = paste("Highest correlation at lag of", max_lag$lag, "days"),
       x = "Lag (Days)", y = "Correlation Coefficient") +
  theme_minimal()

###################
### normal plot ###
###################

ggplot() + 
  geom_ribbon(data = swot_73202, 
              aes(x = time_str, 
                  ymin = norm - error_m, 
                  ymax = norm + error_m), 
              color = "gray",
              alpha = 0.2) + # Very light transparency
  geom_line(data = gw_28171, 
            aes(x=Time, y = norm),
            color = "steelblue", linewidth = 1) +
  geom_line(data = swot_73202,
            aes(x=time_str,
                y= norm),
            color = "red2", linewidth = 1) + 
  geom_point(data = swot_73202,
             aes(x=time_str,
                 y = norm),
             color = "black", size = 2, shape = 21,fill = "red", stroke = 1) +
  scale_x_datetime(date_breaks = "3 month", 
                   date_labels = "%b") +
  labs(
    title = paste("Groundwater Levels at Lake Oliver \nUSGS Site 282202081384601"),
    x = "",
    y = "Water Level Normalized (m)") +
  theme_minimal()

##########################
### Plot with lag time ###
##########################

lag_days <- as.integer(max_lag)[1] #max lag calculated earlier

#shifting swot times by lag times
swot_73202 <- swot_73202 %>%
  mutate(time_shifted = time_str + days(lag_days))


ggplot() + 
  # Original Groundwater (The Reference)
  geom_line(data = gw_28171, 
            aes(x = Time, y = norm, color = "Groundwater"), 
            linewidth = 1) +
  
  # Shifted SWOT (To check the fit)
  geom_line(data = swot_73202,
            aes(x = time_shifted, y = norm, color = "SWOT (Shifted -72 Days)"), 
            linewidth = 1) + 
  
  geom_point(data = swot_73202,
             aes(x = time_shifted, y = norm),
             color = "black", size = 2, shape = 21, fill = "red", stroke = 1) +
  
  # Grid and Axis Labels
  scale_x_datetime(date_breaks = "3 month", 
                   date_labels = "%b %Y") +
  
  scale_color_manual(values = c("Groundwater" = "steelblue", 
                                "SWOT (Shifted -72 Days)" = "red2")) +
  
  labs(
    title = "Lake Oliver Lag Adjusted",
    subtitle = "SWOT shifted back 72 days to match Groundwater lead",
    x = "",
    y = "Water Level Normalized (m)",
    color = "") +
  
  theme_minimal() +
  theme(legend.position = "bottom", 
        axis.text.x = element_text(angle = 45, hjust = 1))

###################################################
# Combined Plot: Unadjusted vs. Lag-Adjusted SWOT #
###################################################

swot_long <- bind_rows(
  swot_73202 %>% mutate(PlotTime = time_str, Type = "Original (Unadjusted)"),
  swot_73202 %>% mutate(PlotTime = time_shifted, Type = paste("Adjusted (Lag:", lag_days, "Days)")))

ggplot() +
  geom_line(data = gw_28171, 
            aes(x = Time, y = norm, color = "Groundwater"), 
            linewidth = 0.8) +
  geom_line(data = swot_long, 
            aes(x = PlotTime, y = norm, color = "SWOT"), 
            linewidth = 0.8) +
  geom_point(data = swot_long, 
             aes(x = PlotTime, y = norm, fill = "SWOT"), 
             color = "black", size = 2, shape = 21, stroke = 0.5) +
  
  #creates side-by-side plots with locked axes
  facet_wrap(~Type) + 
  
  # Aesthetics and Labels
  scale_x_datetime(date_breaks = "4 months", date_labels = "%b %y") +
  scale_color_manual(values = c("Groundwater" = "steelblue", "SWOT" = "red2"),
                     breaks = "Groundwater") +
  scale_fill_manual(values = c("SWOT" = "red2")) +
  labs(
    title = "Lake Oliver SWOT Levels and Groundwater Measurement Plots", 
    y = "Normalized Water Level (m)") +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    axis.title.x = element_blank(),
    legend.title = element_blank(),   
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold", size = 11))
##############################################################################

#Function to auto compare groundwater data sites to swot data 
#Inputs of groundwater data and swot data sets set up above

#If statements to see what agency of data it is so it can use proper
#Swot ids that correlate to those site numbers and only use those site no.
#loops through associated agency gw sites and correlated lake Swot ids

#Setups for both data sets to be plottable
# #Swot quality flags of <= 1
# #mutate for time to be in ymd_hms for ggplot
#`#groundwater data converted to meters
# #Normalizing both by subtracting the mean of each

#Checks to see if more than 10 data points for both swot and gw sites 
#If true goes into plotting that gw site vs swot passes of related lake

#######################################################################
comparison <- function(gw_data, swot_data, mapping_table) {

  results_list <- list()
  agency <- gw_data$AgencyCd[1]
  
  # Determine the column name based on agency
  level_col <- if (agency == "USGS") {
    "Water.level.in.feet.relative.to.NGVD29"
  } else {
    "Water.level.in.feet.relative.to.NAVD88"
  }
  
  for (i in 1:nrow(mapping_table)) {
    current_site <- mapping_table$noSiteNo[i]
    current_lake <- mapping_table$lake_id[i]
    Lake_Geology <- mapping_table$lkgeologyCAT[i]
    Groundwater_Geology <- mapping_table$geologyCAT[i]
    
    # Filter data
    gw <- gw_data %>% filter(SiteNo == current_site) %>% mutate(Time = ymd_hms(Time))
    swot <- swot_data %>% filter(lake_id == current_lake, quality_f <= 1) %>% mutate(time_str = ymd_hms(time_str))
    
    # Trim SWOT with a 100-day Buffer
    if(nrow(gw) > 0) {
      gw_range <- range(gw$Time, na.rm = TRUE)
      
      # Expands the window by 100 days on both ends
      buffer_start <- gw_range[1] - days(100)
      buffer_end   <- gw_range[2] + days(100)
      
      swot <- swot %>% 
        filter(time_str >= buffer_start & time_str <= buffer_end)
    }

    # Tukey Filter (IQR) Outlier Removal
    if(nrow(swot) > 5){
      q1 <- quantile(swot$wse, 0.25, na.rm = TRUE)
      q3 <- quantile(swot$wse, 0.75, na.rm = TRUE)
      iqr <- q3 - q1
      swot <- swot %>% 
        filter(wse >= (q1 - 1.5 * iqr) & wse <= (q3 + 1.5 * iqr))
    }
    ###############################################################
    
    if(nrow(gw) <= 10 | nrow(swot) <= 10) next

    gw$m <- gw[[level_col]] / 3.281
    gw$norm <- gw$m - mean(gw$m, na.rm = TRUE)
    swot$norm <- swot$wse - mean(swot$wse)
    
    # 2. Cross-Correlation Calculation
    all_dates <- data.frame(Time = seq(min(gw$Time), max(gw$Time), by = "1 day"))
    gw_daily <- gw %>% mutate(Time = as.Date(Time)) %>% group_by(Time) %>% summarize(gw_norm = mean(norm, na.rm = TRUE))
    swot_daily <- swot %>% mutate(Time = as.Date(time_str)) %>% group_by(Time) %>% summarize(swot_norm = mean(norm, na.rm = TRUE))
    
    combined_ts <- all_dates %>%
      mutate(Time = as.Date(Time)) %>%
      left_join(gw_daily, by = "Time") %>%
      left_join(swot_daily, by = "Time") %>%
      mutate(swot_interp = na.approx(swot_norm, na.rm = FALSE),
             gw_interp = na.approx(gw_norm, na.rm = FALSE)) %>%
      filter(!is.na(swot_interp), !is.na(gw_interp))
    
    # SAFETY CHECK
    if (nrow(combined_ts) < 2) {
      message(paste("Skipping Site:", current_site, "- No temporal overlap between GW and SWOT."))
      next
    }
    cross_corr <- ccf(combined_ts$gw_interp, combined_ts$swot_interp, lag.max = 100, plot = FALSE)
    ccf_df <- data.frame(lag = cross_corr$lag, acf = cross_corr$acf)
    max_lag <- ccf_df[which.max(abs(ccf_df$acf)), ]
    lag_val <- max_lag$lag
    
    # PLOT A: LAG CORRELATION COEFFICIENTS
    p_corr <- ggplot(ccf_df, aes(x = lag, y = acf)) +

      geom_hline(yintercept = 0, color = "grey70") +
      
      # vertical lines
      geom_segment(aes(x = lag, xend = lag, y = 0, yend = acf, color = (lag == lag_val)), 
                   linewidth = 0.6) + # Adjust linewidth to your liking
      
      # a small dot at the top of each line for visibility
      geom_point(aes(color = (lag == lag_val)), size = 1) +
      
      scale_color_manual(values = c("steelblue", "orange"), guide = "none") +
      
      # Locks the Y-axis from -1 to 1
      scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
      scale_x_continuous(breaks = seq(-100, 100, 20)) + # Explicit breaks help readability
      
      labs(title = paste0("Lag Correlation Analysis: SWOT to Groundwater\n", agency, " Site ID: ", current_site),
           subtitle = paste("Highest correlation (r =", round(max_lag$acf, 2), ") at", lag_val, "days"),
           x = "Lag (Days)", y = "Correlation Coefficient") +
      theme_minimal()
    
    print(p_corr)
    
    # PLOT B: COMPARISON
    swot_long <- bind_rows(
      swot %>% mutate(PlotTime = time_str, Type = "Original (Unadjusted)"),
      swot %>% mutate(PlotTime = time_str + days(lag_val), Type = paste("Adjusted SWOT Lagged by:", lag_val, "Days"))
    )
    
    p_comp <- ggplot() +
      geom_line(data = gw, aes(x = Time, y = norm, color = "Groundwater"), linewidth = 1) +
      geom_line(data = swot_long, aes(x = PlotTime, y = norm, color = "SWOT"), linewidth = 0.5) +
      geom_point(data = swot_long, aes(x = PlotTime, y = norm, fill = "SWOT"), shape = 21, size = 1.25, stroke = 0.5) +
      facet_wrap(~Type) +
      scale_color_manual(values = c("Groundwater" = "steelblue", "SWOT" = "red2"), breaks = "Groundwater") +
      scale_fill_manual(values = c("SWOT" = "red2")) +
      labs(title = paste("SWOT Lake and ",agency,"Groundwater Site Elevations:\nSite ID:", current_site, "\nLake SWOT ID:", current_lake),
           subtitle = paste("Distance:", round(mapping_table$distance[i], 0), "m"),
           y = "Normalized Level (m)", x = NULL) +
      theme_minimal() +
      theme(legend.position = "bottom", legend.title = element_blank(), 
            axis.text.x = element_text(angle = 45, hjust = 1),
            strip.text = element_text(face = "bold"))
    
    print(p_comp)

    # Saves to folder
    #ggsave(paste0("Corr_Site_", current_site, "_Lake_", current_lake, ".png"), 
    #       plot = p_corr, width = 10, height = 6, dpi = 300)
    
    # Saves to folder
    #ggsave(paste0("Comp_Site_", current_site, "_Lake_", current_lake, ".png"), 
    #       plot = p_comp, width = 10, height = 6, dpi = 300)
    
    # Stores Data for the Final Table
    results_list[[length(results_list) + 1]] <- data.frame(
      Ag = agency,
      Site = current_site, 
      Lake = current_lake, 
      Lag = lag_val, 
      Corr = max_lag$acf, 
      Dist = mapping_table$distance[i],
      LKG = Lake_Geology,
      GWG = Groundwater_Geology)
  }
  return(bind_rows(results_list))
}

#############################################################################################################

usgstable <- comparison(usgs_data, swot, con_usgs)
sjrwmdtable <- comparison(SJRWMD_data, swot, con_sjrwmd)

names(usgstable) <- c("Agency", "Groundwater Site ID", "SWOT Lake ID", "Lag (days)", "Correlation r", "Distance (m)", "Lake Soil", "Groundwater Soil")
names(sjrwmdtable) <- c("Agency", "Groundwater Site ID", "SWOT Lake ID", "Lag (days)", "Correlation r", "Distance (m)", "Lake Soil", "Groundwater Soil")

final_combined_table <- rbind(usgstable, sjrwmdtable)
table_print <- final_combined_table

# Converts numeric columns to character before writing
table_print$`Groundwater Site ID` <- format(table_print$`Groundwater Site ID`, scientific = FALSE, trim = TRUE) 

write.csv(table_print, file = 'swot_vs_gw_tab.csv')
###############################################################################################################


######################################################
# Boxplotting Soil types and r lag correlation values#
######################################################
final_combined_table <- final_combined_table %>%
  mutate(Soil_Group = case_when(
    `Groundwater Soil` == "LIMESTONE" ~ "Limestone",
    `Groundwater Soil` == "PEAT" ~ "Peat",
    TRUE ~ "Sand, Silt, Clay"
  ))


ggplot(final_combined_table, aes(x = Soil_Group, y = `Correlation r`, fill = Soil_Group)) +
  # threshold lines
  geom_hline(yintercept = c(0.7, -0.7),
             linetype = "dotted",
             color = "red3",
             linewidth = 0.5) +
  
  # zero line
  geom_hline(yintercept = 0, color = "gray60", size = 0.5) +
  
  geom_boxplot(alpha = 0.75, outlier.shape = NA, width = 0.6) +
  
  geom_jitter(width = 0.15, alpha = 0.35, size = 1.5) +
  
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.5)) +
  
  # boxplots colors
  scale_fill_manual(values = c("Limestone" = "#56B4E9",
                               "Peat" = "#a11111",
                               "Sand, Silt, Clay" = "#F5F5a9")) +
  
  labs(
    title = "Lag Correlation By Soil Category",
    subtitle = "Dotted Lines indicate a strong correlation of |r| > 0.7",
    x = NULL,
    y = "Correlation Coefficient (r)"
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "none", # Hide legend since X-axis labels handle it
    axis.text.x = element_text(face = "bold", size = 11),
    panel.grid.minor = element_blank()
  )                   


############################
# Generate the Scatter Plot#
############################

# p values calculated by soil group 
stats_labels <- final_combined_table %>%
  group_by(Soil_Group) %>%
  filter(n() > 2) %>% # Need at least 3 points for a trend line/p-value
  group_modify(~ {
    fit <- cor.test(.x$`Distance (m)`, .x$`Correlation r`)
    data.frame(
      p_label = paste0("p = ", format.sep = "", round(fit$p.value, 3)),
      r_label = paste0("r = ", round(fit$estimate, 2))
    )})


# p value labels
custom_labels <- c(
  "Limestone"        = "Limestone\n(p = 0.342)",
  "Peat"             = "Peat\n(p = 0.104)",
  "Sand, Silt, Clay" = "Sand, Silt, Clay\n(p = 0.01)"
)

ggplot(final_combined_table, aes(x = `Distance (m)`, y = `Correlation r`, fill = Soil_Group)) +
  
  # Threshold lines
  geom_hline(yintercept = c(0.7, -0.7), linetype = "dotted", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "gray60", size = 0.5) +
  
  # Trend line 
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.75, linetype = "dashed", aes(color = Soil_Group)) +
  
  # Points
  geom_point(aes(fill = Soil_Group),   
             shape = 21,               
             color = "black",          
             size = 3,               
             stroke = 1,               
             alpha = 0.7) +
  
  # Formatting
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.1)))+ # Adds 10% extra space to the right 
  scale_y_continuous(limits = c(-1, 1), breaks = seq(-1, 1, 0.2)) +
  scale_fill_manual(values = c("Limestone" = "#56B4E9", 
                                "Peat"      = "#a11111", 
                                "Sand, Silt, Clay" = "#F5F5a9"), 
                    labels = custom_labels) +
  scale_color_manual(values = c("Limestone" = "#56aaee", 
                               "Peat"      = "#a33333", 
                               "Sand, Silt, Clay" = "#bfbf61"),
                     labels = custom_labels) +
  
  labs(
    title = "Correlation vs. Distance",
    subtitle = "Separate trend lines and p-values shown for each geological group",
    x = "Distance (meters)",
    y = "Correlation Coefficient (r)",
    fill = NULL,
    color = NULL,
  ) +
  
  theme_minimal() +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )



#### For cross-correlation a continuous/regular temporal grid was used. 
#### Linear interpolation of the SWOT data was used to create a daily time series, 
#### which allows the R ccf() function algorithm to test for lags at a 1-day resolution. 
#### The alternative was to go down to the resolution of SWOT data appearances however there
#### is such a low sampling size at *random intervals* and the resolution would need be changed
#### for every site this seemed like the better option. 
#### The linear interpolation is often used for matching hydraulic responses
#### between groundwater and surface water. 
