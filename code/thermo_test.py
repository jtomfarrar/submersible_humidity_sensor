"""
Direct Psychrometric Calculator
--------------------------------
A tool for analyzing psychrometric properties using dry/wet bulb measurements.
Features include:
- Single point analysis of psychrometric properties
- Time series visualization and processing
- Multiple psychrometric chart types
- Reference depression tables

Author: Yu Gao
Date: April 2025
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import xarray as xr
from datetime import datetime

# ======================================================================
# Configuration
# ======================================================================
TEST_NAME = "Lab20250227"
DATA_FILE = f"../data/{TEST_NAME}/233860_20250227_2224.nc"
IMG_DIR = f"../img/{TEST_NAME}/thermo_test/"
A = 0.000662  # Standard Ferrel coefficient (Pa⁻¹·°C⁻¹)

# Create output directory
os.makedirs(IMG_DIR, exist_ok=True)

# ======================================================================
# Data loading
# ======================================================================
print(f"Loading dataset: {DATA_FILE}")
try:
    ds = xr.open_dataset(DATA_FILE)
    
    # Extract temperature data (°C)
    T_dry_C_array = ds['temperature'].values
    T_wet_C_array = ds['temperature1'].values
    
    # Get atmospheric pressure (Pa)
    if 'default_atmospheric_pressure' in ds and ds['default_atmospheric_pressure'].size > 0:
        pressure_Pa = float(ds['default_atmospheric_pressure'].values[0]) * 10000  # kPa to Pa
        print(f"Using atmospheric pressure from dataset: {pressure_Pa} Pa")
    else:
        pressure_Pa = 101325  # Standard pressure in Pa
        print(f"Using standard atmospheric pressure: {pressure_Pa/1000:.2f} kPa")
    
    # Convert timestamp to pandas datetime
    timestamps = pd.to_datetime(ds.time.values)
    print(f"Successfully loaded {len(timestamps)} data points")
    
    # For analysis, let's use a single point first
    point_index = 0
    T_dry_C = T_dry_C_array[point_index]
    T_wet_C = T_wet_C_array[point_index]
    print(f"\nAnalyzing data point {point_index} at {timestamps[point_index]}")
except Exception as e:
    print(f"Error loading data file: {e}")
    exit(1)

# ======================================================================
# Single point psychrometric analysis
# ======================================================================
print(f"\nINPUT VALUES:")
print(f"  Dry bulb temperature: {T_dry_C}°C")
print(f"  Wet bulb temperature: {T_wet_C}°C")
print(f"  Atmospheric pressure: {pressure_Pa/1000:.2f} kPa")

# Calculate depression
depression = T_dry_C - T_wet_C

# Calculate saturation vapor pressures using August-Roche-Magnus formula
e_sat_dry = 610.94 * np.exp((17.625 * T_dry_C) / (T_dry_C + 243.04))
e_sat_wet = 610.94 * np.exp((17.625 * T_wet_C) / (T_wet_C + 243.04))

# Apply psychrometric equation with Ferrel coefficient
e = e_sat_wet - A * pressure_Pa * depression

# Calculate relative humidity
RH = 100.0 * e / e_sat_dry

# Calculate dew point from vapor pressure
alpha = np.log(e / 610.94)
dew_point = 243.04 * alpha / (17.625 - alpha)
    
# Calculate specific humidity
q = 0.622 * e / (pressure_Pa - 0.378 * e)

# Print calculation details
print(f"\nCALCULATION DETAILS:")
print(f"  Temperature depression: {depression:.3f}°C")
print(f"  Saturation vapor pressure at dry bulb: {e_sat_dry:.2f} Pa")
print(f"  Saturation vapor pressure at wet bulb: {e_sat_wet:.2f} Pa")
print(f"  Calculated actual vapor pressure: {e:.2f} Pa")

# Print final result summary
print("\nFINAL RESULTS:")
print(f"  Relative Humidity: {RH:.1f}%")
print(f"  Dew Point: {dew_point:.2f}°C")
print(f"  Specific Humidity: {q*1000:.2f} g/kg")

# ======================================================================
# Depression table for reference
# ======================================================================
print(f"\nTesting various depressions with dry bulb T = {T_dry_C}°C:")
print("Depression (°C) | RH (%) | Dew Point (°C) | Vapor Pressure (Pa)")
print("----------------------------------------------------------------")

for test_depression in [1.0, 2.0, 3.0, 4.0, 5.0]:
    T_wet_test = T_dry_C - test_depression
    
    # Calculate saturation vapor pressures
    e_sat_dry_test = 610.94 * np.exp((17.625 * T_dry_C) / (T_dry_C + 243.04))
    e_sat_wet_test = 610.94 * np.exp((17.625 * T_wet_test) / (T_wet_test + 243.04))
    
    # Calculate actual vapor pressure
    e_test = e_sat_wet_test - A * pressure_Pa * test_depression
    
    # Calculate RH
    RH_test = 100.0 * e_test / e_sat_dry_test
    
    # Calculate dew point
    alpha_test = np.log(e_test / 610.94)
    dew_point_test = 243.04 * alpha_test / (17.625 - alpha_test)
    
    print(f"{test_depression:13.1f} | {RH_test:6.1f} | {dew_point_test:14.2f} | {e_test:18.1f}")

# ======================================================================
# Depression curve generation
# ======================================================================
print("\nGenerating RH vs depression curve...")
depressions = np.linspace(0, 6, 61)
rh_values = []

for d in depressions:
    test_wet = T_dry_C - d
    e_sat_dry_curve = 610.94 * np.exp((17.625 * T_dry_C) / (T_dry_C + 243.04))
    e_sat_wet_curve = 610.94 * np.exp((17.625 * test_wet) / (test_wet + 243.04))
    e_curve = e_sat_wet_curve - A * pressure_Pa * d
    rh_curve = 100.0 * e_curve / e_sat_dry_curve
    rh_values.append(rh_curve)

# Plot RH vs depression curve
plt.figure(figsize=(10, 6))
plt.plot(depressions, rh_values, 'b-')
plt.xlabel('Temperature Depression (°C)')
plt.ylabel('Relative Humidity (%)')
plt.title(f'RH vs Depression at {T_dry_C}°C Dry Bulb')
plt.grid(True)
plt.xlim(0, 6)
plt.ylim(0, 100)
plt.savefig(f"{IMG_DIR}/rh_depression_curve_{TEST_NAME}.png", dpi=300)


# ======================================================================
# T-Wet vs T-Dry Psychrometric Chart
# ======================================================================
print("\nGenerating T-wet vs T-dry psychrometric chart...")

# Create figure
plt.figure(figsize=(10, 8))

# Define temperature ranges for the chart
dry_bulb_range = np.linspace(12, 24, 100)  # 12°C to 24°C
RH_levels = [30, 50, 70, 90]  # RH lines to show
colors = ['r', 'g', 'b', 'm']  # Colors for each RH line

# Calculate and plot RH lines
for i, rh in enumerate(RH_levels):
    wet_bulb_temps = []
    
    for t_dry in dry_bulb_range:
        # Calculate saturation vapor pressure at dry bulb
        e_sat_dry = 610.94 * np.exp((17.625 * t_dry) / (t_dry + 243.04))
        
        # Calculate actual vapor pressure for this RH
        e_actual = e_sat_dry * rh / 100.0
        
        # Find the wet bulb temperature that gives this vapor pressure
        # Use numerical approach - start with a guess and converge
        t_wet_guess = t_dry - 5  # Initial guess
        
        for _ in range(10):  # Usually converges within a few iterations
            e_sat_wet = 610.94 * np.exp((17.625 * t_wet_guess) / (t_wet_guess + 243.04))
            depression = t_dry - t_wet_guess
            e_from_wet = e_sat_wet - A * pressure_Pa * depression
            
            # Adjust guess based on difference
            t_wet_guess = t_wet_guess + (e_actual - e_from_wet) / 100
        
        wet_bulb_temps.append(t_wet_guess)
    
    # Plot RH line
    plt.plot(dry_bulb_range, wet_bulb_temps, '-', color=colors[i], label=f'RH={rh}%')

# Process time series if available
if len(T_dry_C_array) > 1:
    # Initialize array to store calculated RH values
    rh_values = np.zeros_like(T_dry_C_array)
    
    # Calculate RH for each data point
    for i in range(len(T_dry_C_array)):
        t_dry = T_dry_C_array[i]
        t_wet = T_wet_C_array[i]
        
        # Calculate depression
        depression = t_dry - t_wet
        
        # Calculate vapor pressures
        e_sat_dry = 610.94 * np.exp((17.625 * t_dry) / (t_dry + 243.04))
        e_sat_wet = 610.94 * np.exp((17.625 * t_wet) / (t_wet + 243.04))
        
        # Apply psychrometric equation
        e_actual = e_sat_wet - A * pressure_Pa * depression
        
        # Calculate RH
        rh = 100.0 * e_actual / e_sat_dry
        rh_values[i] = rh
    
    # Create scatter plot with color representing RH
    scatter = plt.scatter(T_dry_C_array, T_wet_C_array, c=rh_values, cmap='viridis', 
                         alpha=0.7, s=20, vmin=0, vmax=100)
    
    # Add colorbar
    cbar = plt.colorbar(scatter, label='Calculated RH (%)')

# Format chart
plt.grid(True)
plt.xlim(12, 24)
plt.ylim(10, 24)
plt.xlabel('Dry Bulb Temperature (°C)')
plt.ylabel('Wet Bulb Temperature (°C)')
plt.title('Psychrometric Chart (T1=dry, T2=wet)')
plt.legend()

# Save the chart
plt.tight_layout()
plt.savefig(f"{IMG_DIR}/psychrometric_twet_tdry_{TEST_NAME}.png", dpi=300)


print(f"T-wet vs T-dry psychrometric chart saved to {IMG_DIR}/psychrometric_twet_tdry_{TEST_NAME}.png")

# ======================================================================
# Time series processing
# ======================================================================
if len(T_dry_C_array) > 1:
    print("\nProcessing full time series data...")
    
    # Initialize arrays to store results
    rh_series = np.zeros_like(T_dry_C_array)
    dew_point_series = np.zeros_like(T_dry_C_array)
    
    # Process each data point
    for i in range(len(T_dry_C_array)):
        t_dry = T_dry_C_array[i]
        t_wet = T_wet_C_array[i]
            
        # Calculate depression
        dep = t_dry - t_wet
            
        # Calculate vapor pressures
        e_sat_dry_i = 610.94 * np.exp((17.625 * t_dry) / (t_dry + 243.04))
        e_sat_wet_i = 610.94 * np.exp((17.625 * t_wet) / (t_wet + 243.04))
        
        # Apply psychrometric equation
        e_i = e_sat_wet_i - A * pressure_Pa * dep
        
        # Calculate RH
        rh_i = 100.0 * e_i / e_sat_dry_i
        rh_series[i] = rh_i
        
        # Calculate dew point
        alpha_i = np.log(e_i / 610.94)
        dew_point_series[i] = 243.04 * alpha_i / (17.625 - alpha_i)
    
    # Plot time series results
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    # Temperature plot
    ax1.plot(timestamps, T_dry_C_array, 'r-', label='Dry Bulb')
    ax1.plot(timestamps, T_wet_C_array, 'b-', label='Wet Bulb')
    ax1.set_ylabel('Temperature (°C)')
    ax1.legend()
    ax1.grid(True)
    
    # RH plot
    ax2.plot(timestamps, rh_series, 'g-')
    ax2.set_ylabel('Relative Humidity (%)')
    ax2.grid(True)
    
    # Dew point plot
    ax3.plot(timestamps, dew_point_series, 'm-')
    ax3.plot(timestamps, T_dry_C_array, 'r--', alpha=0.3)  # Reference dry bulb
    ax3.set_ylabel('Dew Point (°C)')
    ax3.grid(True)
    
    plt.xlabel('Time')
    plt.tight_layout()
    
    # Save figure
    plt.savefig(f"{IMG_DIR}/time_series_{TEST_NAME}.png", dpi=300)
    

print(f"\nAnalysis complete. Images saved to {IMG_DIR}/")

