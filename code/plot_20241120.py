# %%
# plot variables from the Field20241120 dataset

import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd

# %%
# set the directories
test = "Field20241120"

datadir = f"../data/{test}"
img_dir = f"../img/{test}"
# test_start = "2024-11-18T17:18:52"
# test_end = "2024-11-18T18:20:53"
test_start = pd.to_datetime("2024-11-18T17:18:52")
test_end = pd.to_datetime("2024-11-18T18:21:53")

# make the directories
import os
os.makedirs(img_dir, exist_ok=True)

# %%
ds = xr.open_dataset(f"../data/{test}/233860_20241118_1826.nc")
# ds = xr.open_dataset(f"../data/Field20241120/233860_20241121_2006.nc")
ds = ds.sel(time=slice(test_start, test_end))

# %%
## plot variables
fig, ax = plt.subplots(1, 1, figsize=(10, 5))   
plt.plot(ds.time, ds.temperature, label="temperature")
plt.plot(ds.time, ds.temperature1, label="temperature 1")
ax.set_title("Humditiy sensor test 20241118", fontsize=16)
plt.xlabel("Time", fontsize=14)
plt.ylabel("Temperature (degree C)", fontsize=14)
plt.legend()
# plt.show()  # show the plot
# save the plot
plt.savefig(f"{img_dir}/temperature_{test}.png")

# %%
# read the sentinel 6 dataset
df = pd.read_csv(f"../data/Sentinel6/WFIP_Sentinel6.csv")
df = df.set_index("time")

# %%
# select the time range
# Convert index to datetime
df.index = pd.to_datetime(df.index)
df.index.name = 'time'

# Direct selection using loc
# Convert index to datetime if not already
df.index = pd.to_datetime(df.index)

# Get closest timestamps
start_idx = df.index.asof(test_start)
end_idx = df.index.asof(test_end)

# Select data
selected_df = df.loc[start_idx:end_idx]
selected_df.head()

# %%
# plot variables in sentinel 6 against temperature and temperature 1

fig, ax = plt.subplots(1, 1, figsize=(10, 5))
plt.plot(ds.time, ds.temperature, label="temperature")
plt.plot(ds.time, ds.temperature1, label="temperature 1")
# %%
