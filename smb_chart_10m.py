import math
import numpy as np
import matplotlib.pyplot as plt

# =============================================================================
# SMB Wave Prediction Model (Depth-Limited Water)
# =============================================================================

def calculate_adjusted_wind_speed(U10):
    """
    Calculates the adjusted wind speed (Ua) from the 10-meter wind speed (U10).

    This adjustment accounts for the non-linear relationship between measured
    wind speed and the actual wind stress at the water surface, as specified
    in the Shore Protection Manual (SPM 1984).

    Args:
        U10 (float): The wind speed at 10m height over water (m/s).

    Returns:
        float: The adjusted wind speed (Ua) in m/s.
    """
    # Formula from SPM (1984): Ua = 0.71 * U10^1.23
    Ua = 0.71 * (U10**1.23)
    return Ua

def calculate_depth_limited(wind_speed_adjusted, fetch, depth, gravity=9.81):
    """
    Calculates wave properties for depth-limited conditions using the SMB method.

    This function uses SMB-based formulas adapted for shallow or transitional
    depths, where wave growth is influenced by the seabed. The Hs and Ts
    formulas are from the Shore Protection Manual. The duration is calculated
    using the direct formula from the Coastal Engineering Manual.

    Args:
        wind_speed_adjusted (float): Adjusted wind speed (Ua) (m/s).
        fetch (float):      Effective fetch length (m).
        depth (float):      Water depth (m).
        gravity (float):    Acceleration of gravity (m/s^2).

    Returns:
        tuple: A tuple containing Hs (m), Ts (s), and t_min (hours).
    """
    # --- Dimensionless Parameters ---
    # For depth-limited cases, both dimensionless fetch and dimensionless
    # depth are required to characterize the conditions.
    dim_fetch = (gravity * fetch) / wind_speed_adjusted**2
    dim_depth = (gravity * depth) / wind_speed_adjusted**2

    # --- Depth-Limited Significant Wave Height (Hs) ---
    # This formula combines the effects of fetch and depth. The outer tanh
    # term accounts for the depth limitation, while the inner tanh term
    # accounts for the fetch limitation.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    term_h = 0.00565 * (dim_fetch)**0.5
    tanh_depth_h = math.tanh(0.530 * (dim_depth)**0.75)
    Hs = (wind_speed_adjusted**2 / gravity) * 0.283 * tanh_depth_h * math.tanh(term_h / tanh_depth_h)

    # --- Depth-Limited Significant Wave Period (Ts) ---
    # Similar to the height calculation, this formula combines dimensionless
    # fetch and depth to determine the period.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    term_t = 0.0379 * (dim_fetch)**0.333
    tanh_depth_t = math.tanh(0.833 * (dim_depth)**0.375)
    Ts = (wind_speed_adjusted / gravity) * 7.54 * tanh_depth_t * math.tanh(term_t / tanh_depth_t)

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
    t_min_seconds = gt_min_U * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours for hours use

    return Hs, Ts, t_min_hours

def generate_combined_chart_depth_limited(min_wind_speed, max_wind_speed, num_wind_steps,
                                          min_fetch_km, max_fetch_km, num_fetch_steps,
                                          fixed_depth_m):
    """
    Generates and displays a single contour chart for Hs, Ts, and t_min
    for a fixed water depth.

    The chart shows Hs, Ts, and t_min as functions of wind speed and fetch length
    in depth-limited water conditions, using the SMB model. Different black line styles
    are used for each parameter.

    Args:
        min_wind_speed (float): Minimum wind speed for the chart (m/s).
        max_wind_speed (float): Maximum wind speed for the chart (m/s).
        num_wind_steps (int):   Number of steps for wind speed axis.
        min_fetch_km (float):   Minimum fetch length for the chart (km).
        max_fetch_km (float):   Maximum fetch length for the chart (km).
        num_fetch_steps (int):  Number of steps for fetch length axis.
        fixed_depth_m (float):  Fixed water depth for the calculations (m).
    """
    print(f"Generating combined SMB parameters chart for depth = {fixed_depth_m}m...")

    U10_grid, F_grid_km = np.meshgrid(
        np.linspace(min_wind_speed, max_wind_speed, num_wind_steps),
        np.linspace(min_fetch_km, max_fetch_km, num_fetch_steps)
    )
    F_grid_m = F_grid_km * 1000 # Convert fetch from km to meters for calculation

    # Calculate Ua for each U10 in the grid
    Ua_grid = calculate_adjusted_wind_speed(U10_grid)

    Hs_values = np.zeros_like(U10_grid)
    Ts_values = np.zeros_like(U10_grid)
    t_min_values = np.zeros_like(U10_grid)

    for i in range(num_fetch_steps):
        for j in range(num_wind_steps):
            Hs, Ts, t_min = calculate_depth_limited(Ua_grid[i, j], F_grid_m[i, j], fixed_depth_m)
            Hs_values[i, j] = Hs
            Ts_values[i, j] = Ts
            t_min_values[i, j] = t_min

    # Set figure size to A3 landscape (420mm x 297mm converted to inches)
    plt.figure(figsize=(420/25.4, 297/25.4))

    # Plot Hs contours (solid black lines) with 0.25m step
    # Filter out NaN values for Hs before plotting, if any
    Hs_finite_values = Hs_values[np.isfinite(Hs_values)]
    if Hs_finite_values.size > 0:
        hs_levels = np.arange(np.floor(Hs_finite_values.min() * 4) / 4, np.ceil(Hs_finite_values.max() * 4) / 4 + 0.25, 0.25)
        hs_contour = plt.contour(U10_grid, F_grid_km, Hs_values,
                                 levels=hs_levels,
                                 colors='black', linestyles='solid', linewidths=1.0)
        plt.clabel(hs_contour, inline=True, fontsize=10, fmt='Hs=%.2f')
    else:
        print("Warning: No finite Hs values to plot. Check input ranges.")


    # Plot Ts contours (dashed black lines) with 0.5s step
    # Filter out NaN values for Ts before plotting, if any
    Ts_finite_values = Ts_values[np.isfinite(Ts_values)]
    if Ts_finite_values.size > 0:
        ts_levels = np.arange(np.floor(Ts_finite_values.min() * 2) / 2, np.ceil(Ts_finite_values.max() * 2) / 2 + 0.5, 0.5)
        ts_contour = plt.contour(U10_grid, F_grid_km, Ts_values,
                                 levels=ts_levels,
                                 colors='black', linestyles='dashed', linewidths=1.0)
        plt.clabel(ts_contour, inline=True, fontsize=10, fmt='Ts=%.2f')
    else:
        print("Warning: No finite Ts values to plot. Check input ranges.")

    # Plot t_min contours (dotted black lines) with 0.5h step
    # Filter out NaN values for t_min before plotting, if any
    t_min_finite_values = t_min_values[np.isfinite(t_min_values)]
    if t_min_finite_values.size > 0:
        dur_levels = np.arange(np.floor(t_min_finite_values.min() * 2) / 2, np.ceil(t_min_finite_values.max() * 2) / 2 + 0.5, 0.5)
        t_min_contour = plt.contour(U10_grid, F_grid_km, t_min_values,
                                    levels=dur_levels,
                                    colors='black', linestyles='dotted', linewidths=1.0)
        plt.clabel(t_min_contour, inline=True, fontsize=10, fmt='Dur=%.2f')
    else:
        print("Warning: No finite t_min values to plot. Check input ranges.")

    plt.title(f'Combined SMB Wave Parameters (Depth = {fixed_depth_m}m): Hs, Ts, and Duration', fontsize=16)
    plt.xlabel('Wind Speed (U10) (m/s)', fontsize=14)
    plt.ylabel('Fetch (km)', fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.7)
    plt.subplots_adjust(left=0.05, right=0.98, top=0.90, bottom=0.10)
    plt.tick_params(axis='both', which='major', labelsize=12)
    plt.savefig('smb_chart_10m.pdf')
    print("Combined chart generation complete.")


if __name__ == "__main__":
    # Define the parameters for the charts as per user's request
    WIND_SPEED_MIN = 5
    WIND_SPEED_MAX = 35
    NUM_WIND_STEPS = 50

    FETCH_KM_MIN = 1
    FETCH_KM_MAX = 50
    NUM_FETCH_STEPS = 50

    FIXED_DEPTH = 10 # Fixed depth of 10 meters as requested

    generate_combined_chart_depth_limited(WIND_SPEED_MIN, WIND_SPEED_MAX, NUM_WIND_STEPS,
                                          FETCH_KM_MIN, FETCH_KM_MAX, NUM_FETCH_STEPS,
                                          FIXED_DEPTH)
