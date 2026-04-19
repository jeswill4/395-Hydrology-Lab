
# Load Required Packages
import os
import io
import time
import glob
import requests  
import pickle
import numpy as np
import pandas as pd
import geopandas as gpd
from geopandas import GeoDataFrame
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
import matplotlib.ticker as mticker
import matplotlib.cm as cm
import matplotlib.dates as mdates
from matplotlib.patches import Patch
from shapely.geometry import box
from shapely.geometry import mapping
from math import ceil
from concurrent.futures import ThreadPoolExecutor, as_completed 
print(f"✅ Required packages loaded")

# -----------------------------------------------------------------------------
# READ AOI & PLD SHAPEFILES !! USER INPUT !!
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

# Read AOI shapefile (example: hydrologic basins, country, state extent)
AOI_path = r"D:/395/Florida_poly/s_03mr26.shp" 
AOI_gdf = gpd.read_file(AOI_path)
print(f"🔵 AOI Basins: \n {AOI_gdf}")

# Read LakesATLAS (PLD lake information)
PLD_AOI_path = r"C:/Users/jlwilliams90/Desktop/395/PLD_NA/PLD_NA.shp"
PLD_AOI_gdf = gpd.read_file(PLD_AOI_path)
print(f"🔵 PLD Lakes: \n {PLD_AOI_gdf}")


# -----------------------------------------------------------------------------
# HYDROCRON SEARCH QUERY !! USER INPUT !!
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

# Set hydcron search query
feature="PriorLake"
start_time="2023-01-01T00:00:00Z"
end_time="2026-01-01T00:00:00Z"
new_end_time="2026-01-01T00:00:00Z"
fields="lake_id, reach_id, obs_id, overlap, n_overlap, time, time_tai, time_str, wse, wse_u, wse_r_u, wse_std, area_total, area_tot_u, area_detct, area_det_u, layovr_val, xtrk_dist, ds1_l, ds1_l_u, ds1_q, ds1_q_u, ds2_l, ds2_l_u, ds2_q, ds2_q_u, quality_f, dark_frac, ice_clim_f, ice_dyn_f, partial_f, xovr_cal_q, geoid_hght, solid_tide, load_tidef, load_tideg, pole_tide, dry_trop_c, wet_trop_c, iono_c, xovr_cal_c, lake_name, p_res_id, p_lon, p_lat, p_ref_wse, p_ref_area, p_date_t0, p_ds_t0, p_storage, cycle_id, pass_id, continent_id, range_start_time, range_end_time, crid, geometry, PLD_version, collection_shortname, collection_version, granuleUR, ingest_time"

# Define feature ids from AOI shapefile
feature_ids = PLD_AOI_gdf["lake_id"].tolist()

# Print for confirmation
print(f"🔵 LakeSP Product: {feature}")
print(f"🔵 LakeSP Fields: {fields}")
print(f"🔵 Start: {start_time}")
print(f"🔵 End: {end_time}")
print(f"🔵 Pulling data for n={len(feature_ids)} lakes")

# -----------------------------------------------------------------------------
# HYDROCRON FUNCTION
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

def fetch_hydrocron_list_parallel(
    feature="YourFeature", 
    feature_ids=None, 
    start_time="YYYY-MM-DDTHH:MM:SSZ", 
    end_time="YYYY-MM-DDTHH:MM:SSZ", 
    fields="lake_id", 
    output="csv", 
    max_workers=10 
):
    if feature_ids is None: 
        feature_ids = ["123456789"] 
    
    total = len(feature_ids) 
    all_responses = [] 
    failed_responses = [] 
    
    start_clock = time.time()  # ⏱️ Start timer 

    def fetch_one(index, feature_id): 
        url = (
            f"https://soto.podaac.earthdatacloud.nasa.gov/hydrocron/v1/timeseries"
            f"?feature={feature}&feature_id={feature_id}"
            f"&start_time={start_time}&end_time={end_time}"
            f"&fields={fields}&output={output}"
        )
        try:
            print(f"➡️ [{index}/{total}] Fetching: {feature_id}")
            response = requests.get(url, timeout=30).json()
            return {"lake_id": feature_id, "data": response, "index": index}
        except Exception as e:
            return {"lake_id": feature_id, "data": {"error": str(e)}, "index": index}

    with ThreadPoolExecutor(max_workers=max_workers) as executor: 
        futures = [executor.submit(fetch_one, i, fid) for i, fid in enumerate(feature_ids, start=1)] 

        for future in as_completed(futures): 
            result = future.result()
            i = result.get("index", "?")
            if "error" in result.get("data", {}): 
                print(f"❌ [{i}/{total}] Skipping lake {result['lake_id']} due to error: {result['data']['error']}")
                failed_responses.append(result)
                continue
            all_responses.append(result)
            print(f"✅ [{i}/{total}] Completed: {result['lake_id']}") 

    end_clock = time.time()  # ⏱️ End timer 
    duration = end_clock - start_clock
    print(f"\n🏁 Finished in {duration:.2f} seconds. ✅ {len(all_responses)} succeeded, ❌ {len(failed_responses)} failed.\n")

    return all_responses

print(f"✅ Hydrocron function ready!")
    
# -----------------------------------------------------------------------------
# GROUP LAKES BY PFAS-L4 BASINS
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

# Get unique PFAS_L4 values
unique_pfaf_l4 = PLD_AOI_gdf["PFAF_ID"].dropna().unique()
unique_pfaf_l4.sort()
print(f"🔢 Found {len(unique_pfaf_l4)} unique PFAF_L4 basins")

# Step 2: Create output folder
os.makedirs("00_HYDROCRON_BASINS_RAW", exist_ok=True)

# -----------------------------------------------------------------------------
# PULL HYDROCRON DATA PER BASIN
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

for pfaf_id in unique_pfaf_l4:
    output_path = f"00_HYDROCRON_BASINS_RAW/hydrocron_PFAF_{pfaf_id}.pkl"
    
    if os.path.exists(output_path):
        print(f"✅ Skipping PFAF_ID {pfaf_id}, already exists at: {output_path}")
        continue

    basin_lakes_gdf = PLD_AOI_gdf[PLD_AOI_gdf["PFAF_ID"] == pfaf_id]
    basin_lake_ids = basin_lakes_gdf["lake_id"].tolist()

    print(f"\n🌍 PFAF_ID Basin: {pfaf_id} → {len(basin_lake_ids)} lakes")

    result = fetch_hydrocron_list_parallel(
        feature=feature,
        feature_ids=basin_lake_ids,
        start_time=start_time,
        end_time=end_time,
        fields=fields,
        output="csv",
        max_workers=10
    )

    with open(output_path, "wb") as f:
        pickle.dump(result, f)

    print(f"💾 Saved results for PFAF_ID {pfaf_id} to: {output_path}")
    
# -----------------------------------------------------------------------------
# COMBINE ALL PFAF_ID BASIN PICKLES INTO ONE DATAFRAME
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

# Find all pickle files in the 00_HYDROCRON_BASINS_RAW folder
basin_pickles = glob.glob("00_HYDROCRON_BASINS_RAW/hydrocron_PFAF_*.pkl")

# Initialize a list to store all rows
all_data_rows = []

# Loop through each file and extract the observations
for pkl_path in basin_pickles:
    with open(pkl_path, "rb") as f:
        basin_results = pickle.load(f)
    
    # Each `basin_results` is a list of dictionaries like:
    # {"lake_id": ..., "data": [list of dicts]}
    for lake_result in basin_results:
        lake_id = lake_result["lake_id"]
        lake_data = lake_result["data"]["results"]

        if isinstance(lake_data, list):  # Valid result
            for row in lake_data:
                row["lake_id"] = lake_id
                all_data_rows.append(row)
        else:
            # Skip if data has error or is not a list
            continue

# Convert to DataFrame
hydrocron_csv_df = pd.DataFrame(all_data_rows)

# Save as a master pickle
combined_output_path = "00_HYDROCRON_BASINS_RAW/hydrocron_all_basins.pkl"
with open(combined_output_path, "wb") as f:
    pickle.dump(hydrocron_csv_df, f)
print(f"💾 Combined DataFrame saved to: {combined_output_path}")

print(f"💾 Saved results for PFAF_IDs")

# -----------------------------------------------------------------------------
# COMBINE ALL PFAF_ID BASIN PICKLES INTO ONE DATAFRAME
# -----------------------------------------------------------------------------
print("\n" + "-"*50 + "\n")

# Check if output already exists
combined_output_path = ("00_HYDROCRON_BASINS_RAW/hydrocron_all_basins.pkl")

if os.path.exists(combined_output_path):
    print(f"✅ Combined DataFrame already exists. Skipping recombination.")
    
else:
    print("\n" + "-"*50 + "\n")
    print("📦 Combining basin-level pickles into one DataFrame...\n")

    basin_pickles = glob.glob("00_HYDROCRON_BASINS_RAW/hydrocron_PFAF_*.pkl")
    all_dfs = []

    for pkl_path in basin_pickles:
        with open(pkl_path, "rb") as f:
            basin_results = pickle.load(f)

        for lake_result in basin_results:
            lake_id = lake_result.get("lake_id")
            lake_data = lake_result.get("data")

            if isinstance(lake_data, dict) and "results" in lake_data:
                csv_text = lake_data["results"]["csv"]

                try:
                    df = pd.read_csv(io.StringIO(csv_text))
                    df["lake_id"] = lake_id
                    all_dfs.append(df)
                except Exception as e:
                    print(f"⚠️ Failed to parse CSV for lake {lake_id}: {e}")
            else:
                print(f"⚠️ Skipping lake {lake_id}: missing 'results'")

    if all_dfs:
        hydrocron_df = pd.concat(all_dfs, ignore_index=True)
        with open(combined_output_path, "wb") as f:
            pickle.dump(hydrocron_df, f)
        print(f"💾 Combined DataFrame with {len(hydrocron_df)} rows saved to:\n- {combined_output_path}")
        print(hydrocron_df.head(3))
    else:
        print("❌ No valid lake data found. Combined DataFrame not created.")
    