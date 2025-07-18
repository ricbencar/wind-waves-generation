import math
import numpy as np

# =============================================================================
# SMB Wave Prediction Model
# =============================================================================
#
# DESCRIPTION:
# This script implements the Sverdrup-Munk-Bretschneider (SMB) method, a set
# of semi-empirical formulas for predicting wind-generated wave characteristics.
# It calculates the significant wave height (Hs), significant wave period (Ts),
# and the minimum wind duration (t_min) required to generate a fully-developed
# sea state for a given fetch.
#
# The script provides two distinct calculation modes based on water depth:
#   1. Deep Water: Assumes the water depth is large enough that the seabed
#      does not interact with or influence wave generation. This is typically
#      applicable for open-ocean and large, deep lake scenarios.
#   2. Depth-Limited: Applies to transitional or shallow water where the
#      water depth is a significant factor, limiting wave growth in
#      conjunction with wind speed and fetch. This is common in coastal
#      areas, estuaries, and shallower lakes.
#
# METHODOLOGY:
# The SMB method is a foundational tool in coastal engineering for making
# preliminary estimates of wave conditions where complex numerical models
# (like SWAN or WAM) are not required. The core principle is that wave growth
# is limited by either fetch (spatial constraint) or duration (temporal
# constraint).
#
#   - Fetch-Limited: Occurs when the wind has blown for a sufficient time over
#     a given fetch for waves to reach their maximum potential size.
#   - Duration-Limited: Occurs when the wind event is too short for waves to
#     become fully developed over the available fetch.
#
# This script calculates the fetch-limited wave parameters and the minimum
# duration required for that state to be valid for both deep and depth-limited
# conditions.
#
# CRITICAL ASSUMPTIONS & LIMITATIONS:
#   - Steady-State Wind: The model assumes that wind speed and direction are
#     constant over the entire fetch for the full duration. This idealized
#     condition is rarely met perfectly in nature.
#   - Input Quality: The model's accuracy is critically dependent on the
#     quality of the input data. Wind speed should be the standard 10-meter
#     overwater value, and the fetch should be an "effective fetch" that
#     accounts for the geometry of the water body.
#
# BIBLIOGRAPHY:
# The formulas and concepts herein are based on foundational work in coastal
# science and are documented in major engineering references, including:
# - U.S. Army. (2008). Coastal Engineering Manual (EM 1110-2-1100).
#   Washington, DC: U.S. Army Corps of Engineers.
# - U.S. Army. (1984). Shore Protection Manual. Vicksburg, MS: U.S. Army
#   Engineer Waterways Experiment Station.
# - Bretschneider, C. L. (1970). Wave forecasting relations for wave
#   generation. Look Lab, Hawaii, 1(3).
#
# =============================================================================

def calculate_deep_water(wind_speed, fetch, gravity=9.81):
    """
    Calculates fetch-limited wave properties in deep water using the SMB method.

    This function implements the widely-used deep-water SMB formulas, as revised
    by Bretschneider, for conditions where the seabed does not influence wave
    growth. It computes the significant wave height (Hs), significant wave
    period (Ts), and the minimum wind duration (t_min) required for the sea
    state to become fully developed for the given fetch.

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
    # --- Dimensionless Fetch Calculation ---
    # The SMB method is based on dimensional analysis. The dimensionless fetch
    # normalizes the fetch length (F) by the wind speed (U) and gravity (g),
    # forming the primary independent variable for the empirical formulas.
    # Formula: F_hat = g * F / U^2
    dim_fetch = (gravity * fetch) / (wind_speed**2)

    # --- Significant Wave Height (Hs) Calculation ---
    # This is the core Bretschneider formula for fetch-limited significant
    # wave height. It empirically relates the dimensionless wave height to the
    # dimensionless fetch. The constants are derived from field data.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    # Formula: (g * Hs) / U^2 = 0.283 * tanh[0.0125 * (g * F / U^2)^0.42]
    gHs_U2 = 0.283 * math.tanh(0.0125 * (dim_fetch**0.42))
    Hs = gHs_U2 * (wind_speed**2 / gravity)

    # --- Significant Wave Period (Ts) Calculation ---
    # This formula provides the corresponding dimensionless relationship for
    # the significant wave period.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    # Formula: (g * Ts) / U = 7.54 * tanh[0.077 * (g * F / U^2)^0.25]
    # The coefficient 7.54 is equivalent to 2.4 * pi.
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

def calculate_depth_limited(wind_speed, fetch, depth, gravity=9.81):
    """
    Calculates wave properties for depth-limited conditions.

    This function uses SMB-based formulas adapted for shallow or transitional
    depths, where wave growth is influenced by the seabed. The Hs and Ts
    formulas are from the Shore Protection Manual. The duration is calculated
    using the direct formula from the Coastal Engineering Manual.

    Args:
        wind_speed (float): Wind speed at 10m height (m/s).
        fetch (float):      Effective fetch length (m).
        depth (float):      Water depth (m).
        gravity (float):    Acceleration of gravity (m/s^2).

    Returns:
        tuple: A tuple containing Hs (m), Ts (s), and t_min (hours).
    """
    # --- Dimensionless Parameters ---
    # For depth-limited cases, both dimensionless fetch and dimensionless
    # depth are required to characterize the conditions.
    dim_fetch = (gravity * fetch) / wind_speed**2
    dim_depth = (gravity * depth) / wind_speed**2

    # --- Depth-Limited Significant Wave Height (Hs) ---
    # This formula combines the effects of fetch and depth. The outer tanh
    # term accounts for the depth limitation, while the inner tanh term
    # accounts for the fetch limitation.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    term_h = 0.00565 * (dim_fetch)**0.5
    tanh_depth_h = math.tanh(0.530 * (dim_depth)**0.75)
    Hs = (wind_speed**2 / gravity) * 0.283 * tanh_depth_h * math.tanh(term_h / tanh_depth_h)

    # --- Depth-Limited Significant Wave Period (Ts) ---
    # Similar to the height calculation, this formula combines dimensionless
    # fetch and depth to determine the period.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    term_t = 0.0379 * (dim_fetch)**0.333
    tanh_depth_t = math.tanh(0.833 * (dim_depth)**0.375)
    Ts = (wind_speed / gravity) * 7.54 * tanh_depth_t * math.tanh(term_t / tanh_depth_t)

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

    return Hs, Ts, t_min_hours # Corrected to return t_min_hours

def main():
    """
    Main function to drive the user interaction, calculations, and display results.
    """
    while True: # Loop for consecutive calculations
        print()
        print("======================================================")
        print("  Sverdrup-Munk-Bretschneider (SMB) Wave Calculator")
        print("======================================================")

        # --- Get User Inputs for Wind and Fetch ---
        while True:
            try:
                print()
                wind_speed = float(input("Enter wind speed at 10m height (m/s): "))
                fetch_km = float(input("Enter fetch length (km): "))
                fetch_m = fetch_km * 1000  # Convert fetch from km to meters
                break # Exit loop if inputs are valid numbers.
            except ValueError:
                print("Invalid input. Please enter numeric values.")

        # --- Select Water Depth Condition ---
        while True:
            print("\nSelect the water condition:")
            print()
            print("  1: Deep Water (Bottom does not affect waves)")
            print("  2: Depth-Limited (Transitional/Shallow, bottom affects waves)")
            print()
            choice = input("Enter your choice (1 or 2): ")
            if choice in ['1', '2']:
                break
            else:
                print("Invalid choice. Please enter 1 or 2.")

        # --- Perform Calculations based on User's Choice ---
        if choice == '1':
            # If the user chose 'Deep Water', call the corresponding function.
            height, period, duration = calculate_deep_water(wind_speed, fetch_m)
        else: # choice == '2'
            # If the user chose 'Depth-Limited', first get the depth input.
            while True:
                try:
                    depth = float(input("Enter water depth (m): "))
                    break # Exit loop if depth is a valid number.
                except ValueError:
                    print("Invalid input. Please enter a numeric value for depth.")
            # Then, call the depth-limited calculation function.
            height, period, duration = calculate_depth_limited(wind_speed, fetch_m, depth)

        # --- Display Final Results ---
        # The calculated results are printed to the console, formatted for clarity.
        print("\n------------------ RESULTS ------------------")
        print(f"Predicted Significant Wave Height (Hs): {height:.2f} meters")
        print(f"Predicted Significant Wave Period (Ts): {period:.2f} seconds")
        if not math.isnan(duration):
            print(f"Minimum Storm Duration Required:      {duration:.2f} hours")
        else:
            print("Minimum Storm Duration:             Not calculated for this case.")
        print("-------------------------------------------")

        # --- Ask user if they want to perform another calculation ---
        print()
        another_calculation = input("Do you want to perform another calculation? (yes/no): ").lower()
        if another_calculation not in ['yes', 'y', 'Y']:
            print("Exiting SMB Wave Calculator. Goodbye!")
            break # Exit the main loop

# This standard Python construct ensures that the `main()` function is called
# only when this script is executed directly. It prevents the code from running
# automatically if this script is imported as a module into another script.
if __name__ == "__main__":
    main()
