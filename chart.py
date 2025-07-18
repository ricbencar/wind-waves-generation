import math
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# SMB Wave Prediction Model (Deep Water)
# =============================================================================
# This function is adapted from the provided calculator.py to be self-contained
# within this script for chart generation.
def calculate_deep_water(wind_speed, fetch, gravity=9.81):
    """
    Calculates fetch-limited wave properties in deep water using the SMB method.

    Args:
        wind_speed (float): The wind speed at 10m height over water (m/s).
        fetch (float):      The effective fetch length (m).
        gravity (float):    The acceleration due to gravity (m/s^2).

    Returns:
        tuple: A tuple containing:
            - Hs (float): Predicted significant wave height (m).
            - Ts (float): Predicted significant wave period (s).
            - t_min (float): Minimum wind duration for fetch-limited state (hours).
    """
    # Dimensionless Fetch Calculation: F_hat = g * F / U^2
    dim_fetch = (gravity * fetch) / (wind_speed**2)

    # Significant Wave Height (Hs) Calculation
    # (g * Hs) / U^2 = 0.283 * tanh[0.0125 * (g * F / U^2)^0.42]
    gHs_U2 = 0.283 * math.tanh(0.0125 * (dim_fetch**0.42))
    Hs = gHs_U2 * (wind_speed**2 / gravity)

    # Significant Wave Period (Ts) Calculation
    # (g * Ts) / U = 7.54 * tanh[0.077 * (g * F / U^2)^0.25]
    gTs_U = 7.54 * math.tanh(0.077 * (dim_fetch**0.25))
    Ts = gTs_U * (wind_speed / gravity)

    # --- Minimum Wind Duration (t_min) Calculation ---
    # This complex empirical formula calculates the minimum time required for a
    # given wind speed and fetch to generate a fully fetch-limited sea. If the
    # actual storm duration is less than this value, the sea state would be
    # considered duration-limited.
    # Ref: Etemad-Shahidi et al. (2009), based on SPM data.
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours for practical use

    return Hs, Ts, t_min_hours

def generate_combined_chart(min_wind_speed, max_wind_speed, num_wind_steps,
                            min_fetch_km, max_fetch_km, num_fetch_steps):
    """
    Generates and displays a single contour chart for Hs, Ts, and t_min.

    The chart shows Hs, Ts, and t_min as functions of wind speed and fetch length
    in deep water conditions, using the SMB model. Different black line styles
    are used for each parameter.

    Args:
        min_wind_speed (float): Minimum wind speed for the chart (m/s).
        max_wind_speed (float): Maximum wind speed for the chart (m/s).
        num_wind_steps (int):   Number of steps for wind speed axis.
        min_fetch_km (float):   Minimum fetch length for the chart (km).
        max_fetch_km (float):   Maximum fetch length for the chart (km).
        num_fetch_steps (int):  Number of steps for fetch length axis.
    """
    print("Generating combined SMB parameters chart...")

    wind_speeds = np.linspace(min_wind_speed, max_wind_speed, num_wind_steps)
    fetch_kms = np.linspace(min_fetch_km, max_fetch_km, num_fetch_steps)

    WS_grid, F_grid_km = np.meshgrid(wind_speeds, fetch_kms)
    F_grid_m = F_grid_km * 1000 # Convert fetch from km to meters for calculation

    Hs_values = np.zeros_like(WS_grid)
    Ts_values = np.zeros_like(WS_grid)
    t_min_values = np.zeros_like(WS_grid)

    for i in range(num_fetch_steps):
        for j in range(num_wind_steps):
            Hs, Ts, t_min = calculate_deep_water(WS_grid[i, j], F_grid_m[i, j])
            Hs_values[i, j] = Hs
            Ts_values[i, j] = Ts
            t_min_values[i, j] = t_min

    # Set figure size to A3 landscape (420mm x 297mm converted to inches)
    plt.figure(figsize=(420/25.4, 297/25.4))

    # Plot Hs contours (solid black lines) with 0.25m step
    hs_levels = np.arange(np.floor(Hs_values.min() * 4) / 4, np.ceil(Hs_values.max() * 4) / 4 + 0.25, 0.25)
    hs_contour = plt.contour(WS_grid, F_grid_km, Hs_values,
                             levels=hs_levels, # Use specified levels
                             colors='black', linestyles='solid', linewidths=1.0)
    plt.clabel(hs_contour, inline=True, fontsize=10, fmt='Hs=%.2f') # Increased fontsize, 2 decimal places

    # Plot Ts contours (dashed black lines) with 0.5s step
    ts_levels = np.arange(np.floor(Ts_values.min() * 2) / 2, np.ceil(Ts_values.max() * 2) / 2 + 0.5, 0.5)
    ts_contour = plt.contour(WS_grid, F_grid_km, Ts_values,
                             levels=ts_levels, # Use specified levels
                             colors='black', linestyles='dashed', linewidths=1.0)
    plt.clabel(ts_contour, inline=True, fontsize=10, fmt='Ts=%.2f') # Increased fontsize, 2 decimal places

    # Plot t_min contours (dotted black lines) with 0.5h step
    # Filter out NaN values for t_min before plotting, if any
    t_min_finite_values = t_min_values[np.isfinite(t_min_values)]
    if t_min_finite_values.size > 0:
        dur_levels = np.arange(np.floor(t_min_finite_values.min() * 2) / 2, np.ceil(t_min_finite_values.max() * 2) / 2 + 0.5, 0.5)
        t_min_contour = plt.contour(WS_grid, F_grid_km, t_min_values,
                                    levels=dur_levels, # Use specified levels
                                    colors='black', linestyles='dotted', linewidths=1.0)
        plt.clabel(t_min_contour, inline=True, fontsize=10, fmt='Dur=%.2f') # Increased fontsize, 2 decimal places
    else:
        print("Warning: No finite t_min values to plot. Check input ranges.")

    plt.title('Combined SMB Wave Parameters: Hs, Ts, and Duration', fontsize=16) # Increased title fontsize
    plt.xlabel('Wind Speed (m/s)', fontsize=14) # Increased xlabel fontsize
    plt.ylabel('Fetch (km)', fontsize=14) # Added ylabel and increased fontsize

    plt.grid(True, linestyle='--', alpha=0.7)
    # Adjust subplot parameters to add minimal margin from corners
    plt.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.10) # Adjusted margins slightly for labels
    plt.tick_params(axis='both', which='major', labelsize=12) # Increased tick label fontsize
    plt.show()
    print("Combined chart generation complete.")


if __name__ == "__main__":
    # Define the parameters for the charts as per user's request
    WIND_SPEED_MIN = 5
    WIND_SPEED_MAX = 35
    NUM_WIND_STEPS = 50

    FETCH_KM_MIN = 1
    FETCH_KM_MAX = 50
    NUM_FETCH_STEPS = 50

    generate_combined_chart(WIND_SPEED_MIN, WIND_SPEED_MAX, NUM_WIND_STEPS,
                            FETCH_KM_MIN, FETCH_KM_MAX, NUM_FETCH_STEPS)
