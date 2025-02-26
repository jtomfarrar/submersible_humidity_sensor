# %%
# plot variables from the Field20241118 dataset

import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import xarray as xr
import pandas as pd

# %%
# set the directories
test = "Field20241216"
shs_test_file = f"../data/{test}/233860_20250213_1924.rsk"
# shs_test_file = f"../data/{test}/233860_20250213_1924.nc"
datadir = f"../data/{test}"
img_dir = f"../img/{test}"

# closest WFIP Station is #5, occupied by Sentinel 6. 
# Location of test was directly adjacent to Sentinel 6 / WFIP Station 5: 41.15009, -71.18407  

shs_start_time = pd.to_datetime("2024-12-16T16:04:00")
shs_end_time = pd.to_datetime("2024-12-16T17:06:00")

# make the directories
import os
os.makedirs(img_dir, exist_ok=True)

# %%
# ds = xr.open_dataset(shs_test_file)
# ds = ds.sel(time=slice(shs_start_time, shs_end_time))

# # %%
# ## plot variables
# fig, ax = plt.subplots(1, 1, figsize=(10, 5))   
# plt.plot(ds.time, ds.temperature, label="temperature")
# plt.plot(ds.time, ds.temperature1, label="temperature 1")
# ax.set_title("Humditiy sensor test 20241118", fontsize=16)
# plt.xlabel("Time", fontsize=14)
# plt.ylabel("Temperature (degree C)", fontsize=14)
# plt.legend()
# # plt.show()  # show the plot
# # save the plot
# plt.savefig(f"{img_dir}/temperature_{test}.png")

# %%
# read the sentinel 6 dataset
df = pd.read_csv(f"../data/Sentinel6/WFIP_Sentinel6.csv")


# %%
# Set the time column as index
df['time'] = pd.to_datetime(df['time'])
df = df.set_index('time')

# %%
# sentinel 6 data
from datetime import timedelta
# Define time range for plotting
# Set a window several days wide to give context but center it on the SHS period
days_before = 3  # Show 3 days before SHS starts
days_after = 3   # Show 3 days after SHS ends

# Calculate the window
start_time = shs_start_time - timedelta(days=days_before)
end_time = shs_end_time + timedelta(days=days_after)

# Use floor/ceiling to normalize time boundaries for more robust selection
start_floored = start_time.floor('H')  # Round down to nearest hour
end_ceiling = end_time.ceil('H')       # Round up to nearest hour

print(f"Full data time range: {df.index.min()} to {df.index.max()}")
print(f"SHS period: {shs_start_time} to {shs_end_time}")
print(f"Selection window: {start_floored} to {end_ceiling}")

# Select data using the normalized boundaries
week_data = df[(df.index >= start_floored) & (df.index <= end_ceiling)].copy()


# Check if we have data
if week_data.empty:
    print(f"WARNING: No data found in the specified time range ({start_time} to {end_time})")
    print(f"Index range in original dataframe: {df.index.min()} to {df.index.max()}")
    # Fallback - take the last few days of data regardless of the time
    week_data = df.iloc[-96:].copy()  # Last 4 days (assuming hourly data)
    print(f"Using fallback: Last {len(week_data)} records from {week_data.index.min()} to {week_data.index.max()}")
else:
    print(f"Selected {len(week_data)} records from {week_data.index.min()} to {week_data.index.max()}")

# Get all columns to plot from the Sentinel 6 dataset
print("Available columns in Sentinel 6 dataset:", df.columns.tolist())

# Use all columns except any non-numeric ones
# First check the data types
dtypes = df.dtypes
numeric_cols = df.select_dtypes(include=['number']).columns.tolist()
print(f"Found {len(numeric_cols)} numeric columns to plot")

# These are the parameters we'll plot
parameters = numeric_cols

# Create a directory for individual plots
individual_plots_dir = f"{img_dir}/individual_plots"
os.makedirs(individual_plots_dir, exist_ok=True)
print(f"Created directory for individual plots: {individual_plots_dir}")

# Create individual plots for each parameter
print(f"Creating individual plots for {len(parameters)} variables...")

for param in parameters:
    # Create a new figure for each parameter
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Plot title with SHS period information
    plt.title(f'Sentinel 6: {param} with SHS period highlighted', fontsize=14)
    
    # Check for NaN values
    if week_data[param].isna().all():
        ax.text(0.5, 0.5, f"No data available for {param}", 
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12)
    else:
        # Plot the data
        ax.plot(week_data.index, week_data[param], 'b-', linewidth=2)
        
        # Calculate statistics
        mean_val = week_data[param].mean()
        max_val = week_data[param].max()
        min_val = week_data[param].min()
        std_val = week_data[param].std()
        
        # Add statistics box
        stats_text = (f"Statistics:\n"
                     f"Mean: {mean_val:.3f}\n"
                     f"Max: {max_val:.3f}\n"
                     f"Min: {min_val:.3f}\n"
                     f"Std: {std_val:.3f}")
        
        ax.text(0.02, 0.95, stats_text, transform=ax.transAxes, fontsize=10,
                bbox=dict(boxstyle="round,pad=0.5", fc="white", alpha=0.9))
    
    # Extend x-axis limits beyond the data to include SHS period
    data_min = week_data.index.min()
    data_max = week_data.index.max()
    
    # Ensure the x-axis includes the SHS period even if it's outside the data range
    plot_min = min(data_min, shs_start_time - timedelta(hours=1))
    plot_max = max(data_max, shs_end_time + timedelta(hours=1))
    
    ax.set_xlim(plot_min, plot_max)
    
    # Add SHS start and end time vertical lines
    ax.axvline(x=shs_start_time, color='r', linestyle='--', 
              linewidth=2, label='SHS Start')
    ax.axvline(x=shs_end_time, color='g', linestyle='--', 
              linewidth=2, label='SHS End')
    
    # Add shaded region between SHS start and end
    ax.axvspan(shs_start_time, shs_end_time, alpha=0.2, color='yellow', label='SHS Period')
    
    # If the SHS period is outside our data range, add a note
    if shs_start_time > data_max or shs_end_time < data_min:
        ax.text(0.5, 0.5, "Note: SHS period is outside the available data range", 
                horizontalalignment='center', verticalalignment='center',
                transform=ax.transAxes, fontsize=12, 
                bbox=dict(boxstyle="round,pad=0.5", fc="lightpink", alpha=0.8))
    
    # Labels and grid
    ax.set_xlabel('Time', fontsize=12)
    ax.set_ylabel(param, fontsize=12)
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    
    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
    fig.autofmt_xdate()
    
    # Add footer
    shs_duration = (shs_end_time - shs_start_time).total_seconds() / 60  # in minutes
    footer_text = f"Data period: {week_data.index.min().strftime('%Y-%m-%d %H:%M')} to {week_data.index.max().strftime('%Y-%m-%d %H:%M')}"
    footer_text += f"\nSHS Period: {shs_start_time.strftime('%Y-%m-%d %H:%M')} to {shs_end_time.strftime('%Y-%m-%d %H:%M')} (Duration: {shs_duration:.1f} min)"
    fig.text(0.5, 0.01, footer_text, ha='center', fontsize=8, style='italic')
    
    # Save individual plot
    # Clean parameter name for filename (replace special characters)
    clean_param = param.replace("/", "_").replace(" ", "_").replace("(", "").replace(")", "")
    individual_filename = f"{individual_plots_dir}/sentinel6_{clean_param}.png"
    plt.tight_layout(rect=[0, 0.02, 1, 0.96])  # Adjusted to make room for two-line footer
    plt.savefig(individual_filename, dpi=300, bbox_inches='tight')
    plt.close()

print(f"Created {len(parameters)} individual plots in {individual_plots_dir}")

# %%
