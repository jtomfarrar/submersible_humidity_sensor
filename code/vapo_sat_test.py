# %%
# test functions in vapo_sat.py
from vapo_sat import es, qs, dqsdT, LvK, Twet_autodiff, C, RdoRv, Cp

# open test files
import xarray as xr
import numpy as np  # Regular numpy for data analysis (not for autodiff)
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import pandas as pd

# Set test dataset path
test = "Lab20250227"
shs_test_file = f"../data/{test}/233860_20250227_2224.nc"

# Open the dataset
ds = xr.open_dataset(shs_test_file)
# ds = ds.sel(time=slice(shs_start_time, shs_end_time))
# Extract temperature data
temp_C = ds['temperature'].values
temp1_C = ds['temperature1'].values

# Get atmospheric pressure from dataset - convert from bar to Pa
# From the image, we can see the default_atmospheric_pressure is 10.13
if 'default_atmospheric_pressure' in ds and ds['default_atmospheric_pressure'].size > 0:
    # Value appears to be in bar (10.13), convert to Pa (1 bar = 100,000 Pa)
    pressure_Pa = float(ds['default_atmospheric_pressure'].values[0]) * 100000
    print(f"Using atmospheric pressure from dataset: {pressure_Pa/1000:.2f} kPa")
else:
    pressure_Pa = 101325.0  # Standard pressure in Pa
    print(f"Using standard atmospheric pressure: {pressure_Pa/1000:.2f} kPa")


# Convert time to pandas datetime for better plotting
time_pd = pd.to_datetime(ds.time.values)

# %%
# Test 1: Basic comparisons between the two temperature sensors
print("\nTest 1: Temperature sensor comparison")

# Calculate basic statistics
temp_diff = temp_C - temp1_C
mean_diff = np.mean(temp_diff)
std_diff = np.std(temp_diff)
max_diff = np.max(temp_diff)
min_diff = np.min(temp_diff)

print(f"Mean temperature difference (T1-T2): {mean_diff:.3f}°C")
print(f"Standard deviation of difference: {std_diff:.3f}°C")
print(f"Range of difference: {min_diff:.3f}°C to {max_diff:.3f}°C")

# %%
# Test 2: Calculate RH assuming one sensor is dry bulb and one is wet bulb
print("\nTest 2: RH calculation using dry/wet bulb pair")

# Function to calculate RH from dry and wet bulb temperatures
def calculate_RH_from_wet_dry(T_dry_C, T_wet_C, pressure_Pa):
    """Calculate relative humidity from dry and wet bulb temperatures"""
    # Skip invalid cases where wet bulb > dry bulb
    if T_wet_C > T_dry_C:
        return np.nan
        
    T_dry_K = T_dry_C + C
    T_wet_K = T_wet_C + C
    
    # Calculate saturation specific humidity at wet bulb temperature
    q_sat_wet = qs(pressure_Pa, T_wet_C)
    
    # Calculate actual specific humidity using psychrometric equation
    Lv = LvK((T_dry_K + T_wet_K) / 2)  # Mean temperature for latent heat
    q = q_sat_wet - (Cp/Lv) * (T_dry_C - T_wet_C)
    
    # Skip if calculated q is negative (physically impossible)
    if q <= 0:
        return np.nan
    
    # Calculate saturation specific humidity at dry bulb temperature
    q_sat_dry = qs(pressure_Pa, T_dry_C)
    
    # Calculate relative humidity (as a fraction)
    RH = q / q_sat_dry
    
    # Skip if RH > 100% (physically questionable)
    if RH > 1.0:
        return np.nan
        
    return RH * 100  # Return as percentage

# Calculate RH assuming temp is dry bulb and temp1 is wet bulb
rh_values = []
for i in range(len(temp_C)):
    rh = calculate_RH_from_wet_dry(temp_C[i], temp1_C[i], pressure_Pa)
    rh_values.append(rh)

rh_array = np.array(rh_values)
valid_rh = ~np.isnan(rh_array)

# Print statistical summary
if np.any(valid_rh):
    print(f"Valid RH measurements: {np.sum(valid_rh)} out of {len(rh_array)}")
    print(f"RH range: {np.min(rh_array[valid_rh]):.1f}% to {np.max(rh_array[valid_rh]):.1f}%")
    print(f"Mean RH: {np.mean(rh_array[valid_rh]):.1f}%")
else:
    print("No valid RH measurements (all wet bulb readings > dry bulb)")

# %%
# Test 3: Compare calculated RH with expected RH range for the experiment
print("\nTest 3: Expected RH vs. Calculated RH")

# create a visualization to help assess if the calculated values seem reasonable.

# Create figure with multiple subplots
fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)

subset = -16482  # To limit the number of points plotted for clarity

# Plot temperatures
axes[0].plot(time_pd[subset:], temp_C[subset:], 'r-', label='Temperature 1')
axes[0].plot(time_pd[subset:], temp1_C[subset:], 'b-', label='Temperature 2')
axes[0].set_ylabel('Temperature (°C)')
axes[0].set_title('Temperature Measurements')
axes[0].grid(True)
axes[0].legend()

# Plot temperature difference
axes[1].plot(time_pd[subset:], temp_diff[subset:], 'g-')
axes[1].set_ylabel('T1 - T2 (°C)')
axes[1].set_title('Temperature Difference')
axes[1].grid(True)

# Plot calculated RH where valid
if np.any(valid_rh):
    axes[2].plot(time_pd[valid_rh][subset:], rh_array[valid_rh][subset:], 'b-')
    axes[2].set_ylabel('RH (%)')
    axes[2].set_title('Calculated RH (assuming T1=dry, T2=wet)')
    axes[2].grid(True)
    axes[2].set_ylim(0, 100)
else:
    axes[2].text(0.5, 0.5, 'No valid RH values calculated', 
                 horizontalalignment='center', verticalalignment='center', 
                 transform=axes[2].transAxes)
    axes[2].set_ylabel('RH (%)')
    axes[2].set_title('Calculated RH (assuming T1=dry, T2=wet)')

axes[2].set_xlabel('Time')
axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
fig.autofmt_xdate()
plt.tight_layout()
plt.savefig(f'../img/{test}/temperature_rh_analysis.png', dpi=300)
plt.show()

# %%
# Test 4: Compare with theoretical psychrometric calculations
print("\nTest 4: Theoretical Psychrometric Calculations")

# Create a range of sample dry bulb temperatures based on our data
sample_T_dry = np.linspace(np.min(temp_C), np.max(temp_C), 5)
sample_RH = [30, 50, 70, 90]  # Sample RH values in percent

# Function to calculate wet bulb from dry bulb and RH
def calculate_wet_bulb(T_dry_C, RH_percent, pressure_Pa):
    """Calculate wet bulb temperature from dry bulb and RH"""
    T_dry_K = T_dry_C + C
    
    # Calculate saturation specific humidity at dry bulb
    q_sat = qs(pressure_Pa, T_dry_C)
    
    # Calculate actual specific humidity from RH
    q = (RH_percent/100) * q_sat
    
    # Calculate wet bulb temperature
    try:
        T_wet_K = Twet_autodiff(T_dry_K, q, pressure_Pa)
        return T_wet_K - C  # Return in Celsius
    except Exception as e:
        print(f"Error calculating wet bulb for T={T_dry_C}°C, RH={RH_percent}%: {str(e)}")
        return np.nan

# Create theoretical comparison table
print("\nDry Bulb (°C) | RH (%) | Theoretical Wet Bulb (°C) | Depression (°C)")
print("-" * 70)

theoretical_data = []

for T_dry in sample_T_dry:
    for RH in sample_RH:
        T_wet = calculate_wet_bulb(T_dry, RH, pressure_Pa)
        
        if not np.isnan(T_wet):
            depression = T_dry - T_wet
            print(f"{T_dry:12.1f} | {RH:6.1f} | {T_wet:25.2f} | {depression:14.2f}")
            
            theoretical_data.append({
                'T_dry': T_dry,
                'RH': RH,
                'T_wet': T_wet,
                'Depression': depression
            })

# %%
# Test 5: Compare measured data with theoretical calculations
print("\nTest 5: Comparison with Theoretical Models")

# Create a psychrometric chart with our measured data overlaid
plt.figure(figsize=(12, 8))

# Plot theoretical RH lines
colors = ['r', 'g', 'b', 'm']
RH_values = [30, 50, 70, 90]

for i, rh in enumerate(RH_values):
    dry_temps = np.linspace(np.min(temp_C)-2, np.max(temp_C)+2, 50)
    wet_temps = [calculate_wet_bulb(t, rh, pressure_Pa) for t in dry_temps]
    plt.plot(dry_temps, wet_temps, f'{colors[i]}-', label=f'RH={rh}%')

# Plot our temperature data
if np.any(valid_rh):
    # Find points where we had valid RH calculations
    plt.scatter(temp_C[valid_rh], temp1_C[valid_rh], c=rh_array[valid_rh], 
                cmap='viridis', s=30, alpha=0.7, label='Measured Data')
    cbar = plt.colorbar()
    cbar.set_label('Calculated RH (%)')
else:
    plt.scatter([], [], label='No valid RH measurements')

plt.grid(True)
plt.xlabel('Dry Bulb Temperature (°C)')
plt.ylabel('Wet Bulb Temperature (°C)')
plt.title('Psychrometric Chart with Measured Data')
plt.legend()
plt.savefig(f'../img/{test}psychrometric_comparison.png', dpi=300)
plt.show()

print("\nAnalysis complete!")