# =============================================================================
# SMB Wave Prediction Model (Deep Water)
# =============================================================================

import math
import numpy as np
import matplotlib.pyplot as plt

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

def calculate_deep_water(wind_speed_adjusted, fetch, gravity=9.8066):
    """
    Calculates fetch-limited wave properties in deep water using the Revised SMB method.

    This function uses the deep-water asymptotic limit of the Hurdle & Stive (1989)
    unified equations. This ensures that deep water calculations are perfectly
    consistent with depth-limited calculations as depth increases.

    Args:
        wind_speed_adjusted (float): The adjusted wind speed (Ua) (m/s).
        fetch (float):      The effective fetch length (m).
        gravity (float):    The acceleration due to gravity (m/s^2).

    Returns:
        tuple: A tuple containing:
            - Hs (float): Predicted significant wave height (m).
            - Ts (float): Predicted significant wave period (s).
            - t_min (float): Minimum wind duration for fetch-limited state (hours).
    """
    # --- Dimensionless Fetch Calculation ---
    # F_hat = g * F / Ua^2
    dim_fetch = (gravity * fetch) / (wind_speed_adjusted**2)

    # --- Significant Wave Height (Hs) Calculation (Hurdle & Stive, 1989) ---
    # Deep water limit of Eq 4.1: tanh(depth terms) -> 1
    # Revised Formula: (g * Hs) / Ua^2 = 0.25 * [tanh(4.3e-5 * F_hat)]^0.5
    # Note: tanh^0.5(x) denotes sqrt(tanh(x))
    term_fetch_h = 4.3e-5 * dim_fetch
    gHs_U2 = 0.25 * (math.tanh(term_fetch_h))**0.5
    Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

    # --- Significant Wave Period (Ts) Calculation (Hurdle & Stive, 1989) ---
    # Deep water limit of Eq 4.2: tanh(depth terms) -> 1
    # Revised Formula: (g * Ts) / Ua = 8.3 * [tanh(4.1e-5 * F_hat)]^(1/3)
    term_fetch_t = 4.1e-5 * dim_fetch
    gTs_U = 8.3 * (math.tanh(term_fetch_t))**(1/3)
    Ts = gTs_U * (wind_speed_adjusted / gravity)

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
            Hs, Ts, t_min = calculate_deep_water(Ua_grid[i, j], F_grid_m[i, j])
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
                                 levels=hs_levels, # Use specified levels
                                 colors='black', linestyles='solid', linewidths=1.0)
        plt.clabel(hs_contour, inline=True, fontsize=10, fmt='Hs=%.2f') # Increased fontsize, 2 decimal places
    else:
        print("Warning: No finite Hs values to plot. Check input ranges.")

    # Plot Ts contours (dashed black lines) with 0.5s step
    # Filter out NaN values for Ts before plotting, if any
    Ts_finite_values = Ts_values[np.isfinite(Ts_values)]
    if Ts_finite_values.size > 0:
        ts_levels = np.arange(np.floor(Ts_finite_values.min() * 2) / 2, np.ceil(Ts_finite_values.max() * 2) / 2 + 0.5, 0.5)
        ts_contour = plt.contour(U10_grid, F_grid_km, Ts_values,
                                 levels=ts_levels, # Use specified levels
                                 colors='black', linestyles='dashed', linewidths=1.0)
        plt.clabel(ts_contour, inline=True, fontsize=10, fmt='Ts=%.2f') # Increased fontsize, 2 decimal places
    else:
        print("Warning: No finite Ts values to plot. Check input ranges.")

    # Plot t_min contours (dotted black lines) with 0.5h step
    # Filter out NaN values for t_min before plotting, if any
    t_min_finite_values = t_min_values[np.isfinite(t_min_values)]
    if t_min_finite_values.size > 0:
        dur_levels = np.arange(np.floor(t_min_finite_values.min() * 2) / 2, np.ceil(t_min_finite_values.max() * 2) / 2 + 0.5, 0.5)
        t_min_contour = plt.contour(U10_grid, F_grid_km, t_min_values,
                                    levels=dur_levels, # Use specified levels
                                    colors='black', linestyles='dotted', linewidths=1.0)
        plt.clabel(t_min_contour, inline=True, fontsize=10, fmt='Dur=%.2f') # Increased fontsize, 2 decimal places
    else:
        print("Warning: No finite t_min values to plot. Check input ranges.")

    plt.title('Combined SMB Wave Parameters: Hs, Ts, and Duration', fontsize=16) # Increased title fontsize
    plt.xlabel('Wind Speed (U10) (m/s)', fontsize=14)
    plt.ylabel('Fetch (km)', fontsize=14)

    plt.grid(True, linestyle='--', alpha=0.7)
    # Adjust subplot parameters to add minimal margin from corners
    plt.subplots_adjust(left=0.05, right=0.98, top=0.92, bottom=0.10) # Adjusted margins slightly for labels
    plt.tick_params(axis='both', which='major', labelsize=12) # Increased tick label fontsize
    plt.savefig('smb_chart_deep.pdf')
    print("Combined chart generation complete.")


if __name__ == "__main__":
    # Define the parameters for the charts
    WIND_SPEED_MIN = 5
    WIND_SPEED_MAX = 35
    NUM_WIND_STEPS = 50

    FETCH_KM_MIN = 1
    FETCH_KM_MAX = 50
    NUM_FETCH_STEPS = 50

    generate_combined_chart(WIND_SPEED_MIN, WIND_SPEED_MAX, NUM_WIND_STEPS,
                            FETCH_KM_MIN, FETCH_KM_MAX, NUM_FETCH_STEPS)