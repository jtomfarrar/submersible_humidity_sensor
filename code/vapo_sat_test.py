"""
Psychrometric analysis script for submersible humidity sensor
Processes wet-bulb/dry-bulb temperature data to calculate relative humidity
"""
import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from thermo import es, qs, Lv, Twet, CtoK, Rd, Rv, Cp

# ======================================================================
# Configuration
# ======================================================================
TEST_NAME = "Lab20250227"
DATA_FILE = f"../data/{TEST_NAME}/233860_20250227_2224.nc"
IMG_DIR = f"../img/{TEST_NAME}"

# Create output directory
os.makedirs(IMG_DIR, exist_ok=True)

# ======================================================================
# Data loading and preparation
# ======================================================================
def load_sensor_data(filepath):
    """Load temperature and pressure data from netCDF file"""
    print(f"Loading dataset: {filepath}")
    ds = xr.open_dataset(filepath)
    
    # Extract temperature data (°C)
    temp_C = ds['temperature'].values
    temp1_C = ds['temperature1'].values
    
    # Get atmospheric pressure (Pa)
    if 'default_atmospheric_pressure' in ds and ds['default_atmospheric_pressure'].size > 0:
        pressure_Pa = float(ds['default_atmospheric_pressure'].values[0]) * 100000  # bar to Pa
        print(f"Using atmospheric pressure from dataset: {pressure_Pa/1000:.2f} kPa")
    else:
        pressure_Pa = 101325.0  # Standard pressure in Pa
        print(f"Using standard atmospheric pressure: {pressure_Pa/1000:.2f} kPa")
    
    # Convert timestamp to pandas datetime
    time_pd = pd.to_datetime(ds.time.values)
    
    return temp_C, temp1_C, pressure_Pa, time_pd, ds

# ======================================================================
# Psychrometric calculation methods
# ======================================================================
def calculate_RH_using_energy_balance(T_dry_C, T_wet_C, pressure_Pa):
    """
    Calculate RH using energy balance approach
    
    This method uses the thermodynamic energy balance between dry and wet bulb
    temperatures to derive relative humidity.
    
    Args:
        T_dry_C: Dry bulb temperature (°C)
        T_wet_C: Wet bulb temperature (°C)
        pressure_Pa: Atmospheric pressure (Pa)
        
    Returns:
        Relative humidity (%) or np.nan if calculation is invalid
    """
    # Skip invalid cases
    if T_wet_C >= T_dry_C:
        return np.nan
    
    # Calculate ratio of specific heat to latent heat (temperature dependent)
    Lv_value = Lv(T_dry_C)  # J/kg
    CpoL = Cp / Lv_value    # Ratio of specific heat to latent heat
    
    # Calculate specific humidity at wet bulb
    q_sat_wet = qs(T_wet_C, pressure_Pa)
    
    # Calculate actual specific humidity using energy balance
    q_actual = -CpoL * (T_dry_C - T_wet_C) + q_sat_wet
    
    # Get saturation specific humidity at dry bulb
    q_sat_dry = qs(T_dry_C, pressure_Pa)
    
    # Calculate RH (with physical constraints)
    if q_actual <= 0 or q_actual > q_sat_dry:
        return np.nan
    
    RH = (q_actual / q_sat_dry) * 100
    
    # Validate reasonable range
    if RH > 100 or RH < 5:
        return np.nan
        
    return RH

def calculate_RH_from_wet_dry(T_dry_C, T_wet_C, pressure_Pa, psychrometric_coef=None):
    """
    Calculate relative humidity using psychrometric equation with dynamic coefficient
    
    Args:
        T_dry_C: Dry bulb temperature (°C)
        T_wet_C: Wet bulb temperature (°C)
        pressure_Pa: Atmospheric pressure (Pa)
        psychrometric_coef: Psychrometric coefficient (or None for auto-selection)
        
    Returns:
        Relative humidity (%) or np.nan if calculation is invalid
    """
    # Skip invalid cases
    if T_wet_C >= T_dry_C:
        return np.nan
    
    # Calculate depression
    depression = T_dry_C - T_wet_C
    
    # Skip extreme temperature differences
    if depression > 8.0:
        return np.nan
    
    # Dynamic coefficient based on depression and temperature
    if psychrometric_coef is None:
        # Coefficient decreases as depression increases for better accuracy
        if depression < 2.0:
            psychrometric_coef = 0.000662  # Standard coefficient
        elif depression < 3.0:
            psychrometric_coef = 0.000600
        elif depression < 4.0:
            psychrometric_coef = 0.000550
        else:
            psychrometric_coef = 0.000500  # Reduced for larger depressions
            
        # Adjust for higher temperatures
        if T_dry_C > 25.0:
            psychrometric_coef *= 0.95
    
    # Get saturation vapor pressures
    e_sat_wet = es(T_wet_C, pressure_Pa)
    e_sat_dry = es(T_dry_C, pressure_Pa)
    
    # Calculate actual vapor pressure using psychrometric equation
    e = e_sat_wet - psychrometric_coef * pressure_Pa * depression
    
    # Skip if calculated e is too small relative to saturation
    if e <= 0.001 * e_sat_dry:
        return np.nan
    
    # Calculate relative humidity
    RH = (e / e_sat_dry) * 100
    
    # Validate result is physically reasonable
    if RH > 100 or RH < 5:
        return np.nan
    
    return RH

def calculate_RH_approximate(T_dry_C, T_wet_C, pressure_Pa):
    """
    Calculate RH using simplified empirical approximation
    
    For cases where more rigorous methods fail, this provides a rough estimate.
    
    Args:
        T_dry_C: Dry bulb temperature (°C)
        T_wet_C: Wet bulb temperature (°C)
        pressure_Pa: Atmospheric pressure (Pa)
        
    Returns:
        Approximate relative humidity (%) or np.nan if invalid
    """
    # Skip invalid cases
    if T_wet_C >= T_dry_C:
        return np.nan
    
    # Calculate depression
    depression = T_dry_C - T_wet_C
    
    # Different approximation based on temperature range
    if T_dry_C < 25:
        # Simplified empirical formula for typical room temperatures
        RH = 100 - 5 * depression**1.5
    else:
        # Adjusted for higher temperatures
        RH = 100 - 4.5 * depression**1.6
    
    # Bound to valid range
    if RH < 5 or RH > 100:
        return np.nan
        
    return RH

def calculate_RH_empirical_robust(T_dry_C, T_wet_C, pressure_Pa=101325.0):
    """
    Calculate RH using a simplified but highly robust empirical approach
    that will work even when other methods fail.
    
    Based on simplified versions of standard psychrometric formulas.
    """
    # Basic validity check
    if T_wet_C >= T_dry_C:
        return np.nan
        
    depression = T_dry_C - T_wet_C
    
    # Simplified empirical formula with pressure adjustment
    # Constants optimized for typical indoor conditions (15-25°C)
    pressure_factor = np.sqrt(pressure_Pa / 101325.0)
    
    if depression <= 2.0:
        RH = 100 - 22.0 * depression * pressure_factor
    elif depression <= 4.0:
        RH = 100 - (20.0 + depression) * depression * pressure_factor * 0.97
    else:
        RH = 100 - (18.0 + 1.2 * depression) * depression * pressure_factor * 0.95
    
    # Cap at physical limits but be permissive
    return max(0.0, min(100.0, RH))

def get_best_RH_estimate(T_dry_C, T_wet_C, pressure_Pa):
    """
    Try multiple methods to calculate RH and return the most reliable result
    """
    # Method 1: Energy balance (most physically accurate)
    rh_energy = calculate_RH_using_energy_balance(T_dry_C, T_wet_C, pressure_Pa)
    if not np.isnan(rh_energy):
        return rh_energy
    
    # Method 2: Variable coefficient psychrometric equation
    rh_psychro = calculate_RH_from_wet_dry(T_dry_C, T_wet_C, pressure_Pa)
    if not np.isnan(rh_psychro):
        return rh_psychro
    
    # Method 3: Empirical approximation 
    rh_approx = calculate_RH_approximate(T_dry_C, T_wet_C, pressure_Pa)
    if not np.isnan(rh_approx):
        return rh_approx
    
    # Method 4: Final robust empirical fallback (almost never returns NaN)
    return calculate_RH_empirical_robust(T_dry_C, T_wet_C, pressure_Pa)

def calculate_wet_bulb(T_dry_C, RH_percent, pressure_Pa):
    """
    Calculate wet bulb temperature from dry bulb and RH
    
    Used for generating theoretical comparison lines.
    
    Args:
        T_dry_C: Dry bulb temperature (°C)
        RH_percent: Relative humidity (%)
        pressure_Pa: Atmospheric pressure (Pa)
        
    Returns:
        Wet bulb temperature (°C) or np.nan if calculation fails
    """
    # Calculate actual specific humidity from RH
    q_sat = qs(T_dry_C, pressure_Pa)
    q = (RH_percent/100) * q_sat
    
    # Calculate wet bulb temperature using thermo.py's Twet function
    try:
        T_wet_C = Twet(T_dry_C, q, pressure_Pa)
        return T_wet_C
    except Exception as e:
        print(f"Error calculating wet bulb for T={T_dry_C}°C, RH={RH_percent}%: {str(e)}")
        return np.nan

# ======================================================================
# Analysis and reporting functions
# ======================================================================
def analyze_temperature_stats(temp_C, temp1_C):
    """Calculate and report basic temperature statistics"""
    print("\nPart 1: Temperature Analysis")
    
    temp_diff = temp_C - temp1_C
    mean_diff = np.mean(temp_diff)
    std_diff = np.std(temp_diff)
    max_diff = np.max(temp_diff)
    min_diff = np.min(temp_diff)
    
    print(f"Mean temperature difference (T1-T2): {mean_diff:.3f}°C")
    print(f"Standard deviation of difference: {std_diff:.3f}°C")
    print(f"Range of difference: {min_diff:.3f}°C to {max_diff:.3f}°C")
    
    # Visualize temperature differences
    plt.figure(figsize=(10, 6))
    plt.hist(temp_diff, bins=50)
    plt.xlabel('Temperature difference (T1-T2) °C')
    plt.ylabel('Frequency')
    plt.grid(True)
    plt.title('Histogram of Temperature Differences')
    plt.savefig(f'{IMG_DIR}/temp_diff_histogram.png', dpi=300)
    plt.close()
    
    return temp_diff

def test_psychrometric_calculations(pressure_Pa):
    """Run test calculations with various parameters"""
    print("\nPart 2: Psychrometric Calculations")
    
    # Test with different coefficients
    print("\nTesting psychrometric coefficients for T_dry=19°C, T_wet=15°C:")
    T_dry_test = 19.0
    T_wet_test = 15.0
    
    for coef in [0.000662, 0.000700, 0.000800]:
        rh = calculate_RH_from_wet_dry(T_dry_test, T_wet_test, pressure_Pa, coef)
        print(f"  Psychrometric coef={coef:.6f} → RH={rh:.1f}%")
    
    # Test with different depressions
    print("\nVerifying psychrometric calculation at various temperature depressions:")
    for depression in [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]:
        T_dry = 19.0
        T_wet = T_dry - depression
        rh = calculate_RH_from_wet_dry(T_dry, T_wet, pressure_Pa)
        rh_energy = calculate_RH_using_energy_balance(T_dry, T_wet, pressure_Pa)
        print(f"  T_dry={T_dry}°C, T_wet={T_wet}°C (ΔT={depression}°C) → RH_psychro={rh:.1f}%, RH_energy={rh_energy:.1f}%")

def calculate_RH_for_both_configurations(temp_C, temp1_C, pressure_Pa):
    """Calculate RH assuming each sensor could be wet or dry bulb"""
    print("\nCalculating RH using both sensor configurations...")
    
    # Standard configuration: temp_C is dry bulb, temp1_C is wet bulb
    rh_standard = np.array([get_best_RH_estimate(temp_C[i], temp1_C[i], pressure_Pa) 
                            for i in range(len(temp_C))])
    valid_standard = ~np.isnan(rh_standard)
    
    # Reversed configuration: temp1_C is dry bulb, temp_C is wet bulb
    rh_reversed = np.array([get_best_RH_estimate(temp1_C[i], temp_C[i], pressure_Pa) 
                            for i in range(len(temp_C))])
    valid_reversed = ~np.isnan(rh_reversed)
    
    # Report statistics for both configurations
    if np.any(valid_standard):
        print("\nStandard sensor assignment (temp_C=dry, temp1_C=wet):")
        print(f"  Valid points: {np.sum(valid_standard)} of {len(rh_standard)} ({100*np.sum(valid_standard)/len(rh_standard):.1f}%)")
        print(f"  RH range: {np.min(rh_standard[valid_standard]):.1f}% to {np.max(rh_standard[valid_standard]):.1f}%")
        print(f"  Mean RH: {np.mean(rh_standard[valid_standard]):.1f}%")
    else:
        print("\nStandard sensor assignment: No valid RH values")
    
    if np.any(valid_reversed):
        print("\nReversed sensor assignment (temp1_C=dry, temp_C=wet):")
        print(f"  Valid points: {np.sum(valid_reversed)} of {len(rh_reversed)} ({100*np.sum(valid_reversed)/len(rh_reversed):.1f}%)")
        print(f"  RH range: {np.min(rh_reversed[valid_reversed]):.1f}% to {np.max(rh_reversed[valid_reversed]):.1f}%")
        print(f"  Mean RH: {np.mean(rh_reversed[valid_reversed]):.1f}%")
    else:
        print("\nReversed sensor assignment: No valid RH values")
    
    return rh_standard, valid_standard, rh_reversed, valid_reversed

def generate_theoretical_comparison(pressure_Pa):
    """Generate and display theoretical wet bulb temperatures for different RH levels"""
    print("\nPart 3: Theoretical Model Comparison")
    
    # Generate theoretical table
    print("\nTheoretical values for sample temperatures:")
    print("Dry Bulb (°C) | RH (%) | Theoretical Wet Bulb (°C) | Depression (°C)")
    print("-" * 70)
    
    # Sample temperatures and RH values
    sample_T_dry = np.linspace(15, 25, 3)
    sample_RH = [30, 50, 70, 90]
    
    for T_dry in sample_T_dry:
        for RH in sample_RH:
            T_wet = calculate_wet_bulb(T_dry, RH, pressure_Pa)
            
            if not np.isnan(T_wet):
                depression = T_dry - T_wet
                print(f"{T_dry:12.1f} | {RH:6.1f} | {T_wet:25.2f} | {depression:14.2f}")

def generate_visualizations(temp_C, temp1_C, temp_diff, time_pd, 
                          rh_standard, valid_standard, rh_reversed, valid_reversed,
                          pressure_Pa):
    """Generate and save data visualizations"""
    print("\nPart 4: Generating Visualizations")
    
    # Figure 1: Temperature and RH time series
    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    
    subset = slice(-16482, None)  # Last portion of data
    
    # Plot temperatures
    axes[0].plot(time_pd[subset], temp_C[subset], 'r-', label='Temperature 1')
    axes[0].plot(time_pd[subset], temp1_C[subset], 'b-', label='Temperature 2')
    axes[0].set_ylabel('Temperature (°C)')
    axes[0].set_title('Temperature Measurements')
    axes[0].grid(True)
    axes[0].legend()
    
    # Plot temperature difference
    axes[1].plot(time_pd[subset], temp_diff[subset], 'g-')
    axes[1].set_ylabel('T1 - T2 (°C)')
    axes[1].set_title('Temperature Difference')
    axes[1].grid(True)
    
    # Plot calculated RH - use configuration with more valid points
    if np.sum(valid_standard) >= np.sum(valid_reversed):
        mask = valid_standard[subset]
        rh_display = rh_standard
        rh_label = "Standard Configuration (T1=dry, T2=wet)"
    else:
        mask = valid_reversed[subset]
        rh_display = rh_reversed
        rh_label = "Reversed Configuration (T2=dry, T1=wet)"
    
    if np.any(mask):
        axes[2].plot(time_pd[subset][mask], rh_display[subset][mask], 'b-', 
                     label='Calculated RH')
        axes[2].set_ylabel('RH (%)')
        axes[2].set_title(f'Calculated RH ({rh_label})')
        axes[2].grid(True)
        axes[2].set_ylim(0, 100)
    else:
        axes[2].text(0.5, 0.5, 'No valid RH values calculated', 
                     horizontalalignment='center', verticalalignment='center', 
                     transform=axes[2].transAxes)
        axes[2].set_ylabel('RH (%)')
        axes[2].set_title('Calculated RH')
    
    axes[2].set_xlabel('Time')
    axes[2].xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    fig.autofmt_xdate()
    plt.tight_layout()
    plt.savefig(f'{IMG_DIR}/time_series_analysis.png', dpi=300)
    plt.close()
    
    # Figure 2: Psychrometric chart with measured data
    plt.figure(figsize=(12, 8))
    
    # Plot theoretical RH lines
    colors = ['r', 'g', 'b', 'm']
    RH_values = [30, 50, 70, 90]
    
    temp_min = min(np.min(temp_C), np.min(temp1_C)) - 2
    temp_max = max(np.max(temp_C), np.max(temp1_C)) + 2
    dry_temp_range = np.linspace(temp_min, temp_max, 50)
    
    for i, rh in enumerate(RH_values):
        wet_temps = [calculate_wet_bulb(t, rh, pressure_Pa) for t in dry_temp_range]
        valid = ~np.isnan(wet_temps)
        plt.plot(dry_temp_range[valid], np.array(wet_temps)[valid], 
                 f'{colors[i]}-', label=f'RH={rh}%')
    
    # Plot our temperature data using the configuration with more valid points
    if np.sum(valid_standard) >= np.sum(valid_reversed):
        plt.scatter(temp_C[valid_standard], temp1_C[valid_standard], 
                    c=rh_standard[valid_standard], cmap='viridis', s=30, alpha=0.7)
        plt.title('Psychrometric Chart (T1=dry, T2=wet)')
    else:
        plt.scatter(temp1_C[valid_reversed], temp_C[valid_reversed], 
                    c=rh_reversed[valid_reversed], cmap='viridis', s=30, alpha=0.7)
        plt.title('Psychrometric Chart (T2=dry, T1=wet)')
    
    plt.colorbar(label='Calculated RH (%)')
    plt.grid(True)
    plt.xlabel('Dry Bulb Temperature (°C)')
    plt.ylabel('Wet Bulb Temperature (°C)')
    plt.legend()
    plt.savefig(f'{IMG_DIR}/psychrometric_chart.png', dpi=300)
    plt.close()
    
    # Figure 3: RH vs Depression
    plt.figure(figsize=(10, 6))
    if np.sum(valid_standard) >= np.sum(valid_reversed):
        depression = temp_C[valid_standard] - temp1_C[valid_standard]
        plt.scatter(depression, rh_standard[valid_standard], c=temp_C[valid_standard], 
                    cmap='viridis', alpha=0.6)
        plt.colorbar(label='Dry Bulb Temperature (°C)')
        plt.title('RH vs Temperature Depression (T1=dry, T2=wet)')
    else:
        depression = temp1_C[valid_reversed] - temp_C[valid_reversed]
        plt.scatter(depression, rh_reversed[valid_reversed], c=temp1_C[valid_reversed], 
                    cmap='viridis', alpha=0.6)
        plt.colorbar(label='Dry Bulb Temperature (°C)')
        plt.title('RH vs Temperature Depression (T2=dry, T1=wet)')
    
    plt.xlabel('Temperature Depression (T_dry - T_wet) °C')
    plt.ylabel('Calculated RH (%)')
    plt.grid(True)
    plt.savefig(f'{IMG_DIR}/rh_vs_depression.png', dpi=300)
    plt.close()

# ======================================================================
# Main program execution
# ======================================================================
def main():
    """Main execution function"""
    # Load data
    temp_C, temp1_C, pressure_Pa, time_pd, ds = load_sensor_data(DATA_FILE)
    
    # Analyze temperature data
    temp_diff = analyze_temperature_stats(temp_C, temp1_C)
    
    # Run test calculations
    test_psychrometric_calculations(pressure_Pa)
    
    # Calculate RH for both possible sensor configurations
    rh_standard, valid_standard, rh_reversed, valid_reversed = calculate_RH_for_both_configurations(
        temp_C, temp1_C, pressure_Pa)
    
    # Generate theoretical comparison data
    generate_theoretical_comparison(pressure_Pa)
    
    # Generate visualizations
    generate_visualizations(temp_C, temp1_C, temp_diff, time_pd, 
                          rh_standard, valid_standard, rh_reversed, valid_reversed, 
                          pressure_Pa)
    
    print("\nAnalysis complete!")

# Execute if run as script
if __name__ == "__main__":
    main()