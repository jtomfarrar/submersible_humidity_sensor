# %%
# plot variables from the Field20241118 dataset

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd

# %%
# set the directories
test = "Lab20250227"
shs_test_file = f"../data/{test}/233860_20250227_2224.nc"
# the real-time sentinel data file
sentinel_data_file = f"../data/Sentinel/WFIP_Sentinel2.csv"
datadir = f"../data/{test}"
img_dir = f"../img/{test}"
# lab test start and end time unknown
# use a wider time range to get the data
shs_start_time = pd.to_datetime("2025-02-27T20:06:00")
shs_end_time = pd.to_datetime("2025-02-27T22:23:00")

# make the directories
import os
os.makedirs(img_dir, exist_ok=True)

# %%
ds = xr.open_dataset(shs_test_file)
ds = ds.sel(time=slice(shs_start_time, shs_end_time))

# %%
## plot variables
fig, ax = plt.subplots(1, 1, figsize=(10, 5))   
plt.plot(ds.time, ds.temperature, label="temperature")
plt.plot(ds.time, ds.temperature1, label="temperature 1")
ax.set_title("Humditiy sensor tank test", fontsize=16)
plt.xlabel("Time", fontsize=14)
plt.ylabel("Temperature (degree C)", fontsize=14)
plt.xticks(rotation=45)
plt.grid(True)
plt.legend()
plt.tight_layout()
# plt.show()  # show the plot
# save the plot
plt.savefig(f"{img_dir}/temperature_{test}.png")

# %%
# read the sentinel 6 dataset
df = pd.read_csv(sentinel_data_file)


# %%
# Set the time column as index
df['time'] = pd.to_datetime(df['time'])
df = df.set_index('time')

# plot the air temperature against shs data

max_shs = max(ds.temperature.max(), ds.temperature1.max())
min_shs = min(ds.temperature.min(), ds.temperature1.min())
        
fig, ax = plt.subplots(1, 1, figsize=(10, 6))   
plt.plot(ds.time, ds.temperature, label="temperature")
plt.plot(ds.time, ds.temperature1, label="temperature 1")
plt.plot(df.index, df['temperature_WXT'], label="air temperature WXT")
plt.plot(df.index, df['temperature_HC2'], label="air temperature HC2")
ax.set_title("Humditiy sensor tank test 11/18/2024", fontsize=16)   
# Add footer
shs_duration = (shs_end_time - shs_start_time).total_seconds() / 60  # in minutes
# footer_text = f"Data period: {week_data.index.min().strftime('%Y-%m-%d %H:%M')} to {week_data.index.max().strftime('%Y-%m-%d %H:%M')}"
footer_text = f"SHS Period: {shs_start_time.strftime('%Y-%m-%d %H:%M')} to {shs_end_time.strftime('%Y-%m-%d %H:%M')} (Duration: {shs_duration:.1f} min)"
# footer_text += f"\nclosest WFIP Station is #5, occupied by Sentinel 6. Lat: {lat_sentinel:.4f} Lon: {lon_sentinel:.4f}"
fig.text(0.5, 0.01, footer_text, ha='center', fontsize=8, style='italic')
plt.xlim(shs_start_time , shs_end_time )
plt.ylim(min_shs, max_shs)
plt.xticks(rotation=45)
plt.ylabel("Temperature (degree C)", fontsize=14)
plt.legend()
plt.grid(True)
plt.tight_layout(rect=[0, 0.02, 1, 0.96])  # Adjusted to make room for two-line footer
# plt.show()  # show the plot
# save the plot
plt.savefig(f"{img_dir}/temperature_{test}.png")
