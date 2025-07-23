import math
import numpy as np
import sys # Import sys module to redirect stdout

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
# The script provides two distinct calculation modes based on water depth
# and limiting condition:
#   1. Fetch-Limited (Finite Depth or Deep Water): Assumes wave growth is
#      primarily limited by the available fetch. It can handle both deep water
#      conditions (where the seabed does not influence wave generation) and
#      finite depth conditions (where water depth is a significant factor).
#      Calculates Hs, Ts, and minimum duration for a given fetch and optional depth.
#   2. Duration-Limited (Finite Depth or Deep Water): Calculates wave parameters (Hs, Ts)
#      when the wind event is too short for waves to become fully developed
#      over the available fetch, and water depth also plays a role. Inputs are
#      wind speed, storm duration, and water depth. For this mode, it also
#      calculates the equivalent fetch that would produce the same wave conditions.
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
# conditions. It also provides a dedicated mode for duration-limited scenarios
# in both deep and finite water depth.
#
# CRITICAL ASSUMPTIONS & LIMITATIONS:
#   - Steady-State Wind: The model assumes that wind speed and direction are
#     constant over the entire fetch for the full duration. This idealized
#     condition is rarely met perfectly in nature.
#   - Input Quality: The model's accuracy is critically dependent on the
#     quality of the input data. Wind speed should be the standard 10-meter
#     overwater value, and the fetch should be an "effective fetch" that
#     accounts for the geometry of the water body.
#   - Heuristic for Duration-Limited (Finite Depth): The formulas used for
#     duration-limited conditions in finite depth (Option 2) are an adaptation
#     based on the structure of SMB equations for fetch-limited finite depth.
#     Direct empirical formulas for this specific combined scenario are less
#     common in basic SMB literature. While providing a reasonable estimate,
#     these should be used with awareness of their heuristic nature.
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

def calculate_deep_water(wind_speed_adjusted, fetch, gravity=9.81):
    """
    Calculates fetch-limited wave properties in deep water using the SMB method.

    This function implements the widely-used deep-water SMB formulas, as revised
    by Bretschneider, for conditions where the seabed does not influence wave
    growth. It computes the significant wave height (Hs), significant wave
    period (Ts), and the minimum wind duration (t_min) required for the sea
    state to become fully developed for the given fetch.

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
    # The SMB method is based on dimensional analysis. The dimensionless fetch
    # normalizes the fetch length (F) by the wind speed (U) and gravity (g),
    # forming the primary independent variable for the empirical formulas.
    # Formula: F_hat = g * F / U^2
    dim_fetch = (gravity * fetch) / (wind_speed_adjusted**2)

    # --- Significant Wave Height (Hs) Calculation ---
    # This is the core Bretschneider formula for fetch-limited significant
    # wave height. It empirically relates the dimensionless wave height to the
    # dimensionless fetch. The constants are derived from field data.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    # Formula: (g * Hs) / U^2 = 0.283 * tanh[0.0125 * (g * F / U^2)^0.42]
    gHs_U2 = 0.283 * math.tanh(0.0125 * (dim_fetch**0.42))
    Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

    # --- Significant Wave Period (Ts) Calculation ---
    # This formula provides the corresponding dimensionless relationship for
    # the significant wave period.
    # Ref: U.S. Army (1984), Shore Protection Manual.
    # Formula: (g * Ts) / U = 7.54 * tanh[0.077 * (g * F / U^2)^0.25]
    # The coefficient 7.54 is equivalent to 2.4 * pi.
    gTs_U = 7.54 * math.tanh(0.077 * (dim_fetch**0.25))
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

def calculate_depth_limited(wind_speed_adjusted, fetch, depth, gravity=9.81):
    """
    Calculates wave properties for depth-limited conditions.

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
    # Formula: (g * Ts) / U = 7.54 * tanh[0.0379 * (g * F / U^2)^0.333 / tanh(0.833 * (g * d / U^2)^0.375)]
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
    t_min_hours = t_min_seconds / 3600  # Convert to hours for practical use

    return Hs, Ts, t_min_hours

def calculate_duration_limited(wind_speed_adjusted, duration_hours, depth=None, gravity=9.81):
    """
    Calculates wave properties for duration-limited conditions, considering finite depth if provided.

    This function uses SMB-based formulas directly for duration-limited wave
    growth. If a depth is provided, it attempts to incorporate depth effects
    heuristically based on the structure of fetch-limited depth-limited formulas.
    It also calculates an "equivalent fetch" that would produce the same wave
    conditions if the system were fetch-limited.

    Args:
        wind_speed_adjusted (float): Adjusted wind speed (Ua) (m/s).
        duration_hours (float): Duration of the wind event (hours).
        depth (float, optional): Water depth (m). If None, deep water formulas are used.
        gravity (float):    Acceleration of gravity (m/s^2).

    Returns:
        tuple: A tuple containing:
            - Hs (float): Predicted significant wave height (m).
            - Ts (float): Predicted significant wave period (s).
            - equivalent_fetch_km (float): The equivalent fetch length (km) that would
                                           produce the same wave conditions if fetch-limited.
    """
    # Convert duration from hours to seconds
    duration_seconds = duration_hours * 3600

    # --- Dimensionless Duration Calculation ---
    # Dimensionless duration normalizes the duration (t) by wind speed (U) and gravity (g).
    # Formula: t_hat = g * t / U
    dim_duration = (gravity * duration_seconds) / wind_speed_adjusted

    if depth is None or depth == 0: # Treat as deep water if depth is not provided or zero
        # --- Significant Wave Height (Hs) Calculation for Duration-Limited (Deep Water) ---
        # This empirical formula relates dimensionless wave height to dimensionless duration.
        # Ref: U.S. Army (2008), Coastal Engineering Manual, Figure II-2-4-1.
        # Formula: (g * Hs) / U^2 = 0.283 * tanh[0.000528 * (g * t / U)^0.75]
        gHs_U2 = 0.283 * math.tanh(0.000528 * (dim_duration**0.75))
        Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

        # --- Significant Wave Period (Ts) Calculation for Duration-Limited (Deep Water) ---
        # This empirical formula relates dimensionless wave period to dimensionless duration.
        # Ref: U.S. Army (2008), Coastal Engineering Manual, Figure II-2-4-1.
        # Formula: (g * Ts) / U = 7.54 * tanh[0.00379 * (g * t / U)^0.41]
        gTs_U = 7.54 * math.tanh(0.00379 * (dim_duration**0.41))
        Ts = gTs_U * (wind_speed_adjusted / gravity)

        # --- Equivalent Fetch Calculation (Deep Water) ---
        # Equating the arguments of the tanh function for dimensionless Hs in duration-limited
        # and fetch-limited deep water conditions:
        # 0.000528 * (g * t / U)^0.75 = 0.0125 * (g * F / U^2)^0.42
        # Solving for (g * F / U^2) which is dim_fetch:
        dim_fetch_equivalent = ( (0.000528 / 0.0125) * (dim_duration**0.75) )**(1/0.42)
        equivalent_fetch_m = dim_fetch_equivalent * (wind_speed_adjusted**2 / gravity)
        equivalent_fetch_km = equivalent_fetch_m / 1000 # Convert to kilometers

    else:
        # --- Dimensionless Depth Calculation ---
        dim_depth = (gravity * depth) / wind_speed_adjusted**2

        # --- Heuristic for Duration-Limited Significant Wave Height (Finite Depth) ---
        # This is an adaptation. It combines the structure of the fetch-limited
        # depth-limited formula with the duration-limited deep water parameters.
        # It's based on the assumption that depth limits the maximum possible
        # wave height, similar to how it limits fetch-limited waves.
        # Constants for depth influence are taken from fetch-limited depth-limited Hs.
        term_h_duration = 0.000528 * (dim_duration)**0.75 # From deep water duration-limited Hs
        tanh_depth_h = math.tanh(0.530 * (dim_depth)**0.75) # From fetch-limited depth-limited Hs
        Hs = (wind_speed_adjusted**2 / gravity) * 0.283 * tanh_depth_h * math.tanh(term_h_duration / tanh_depth_h)

        # --- Heuristic for Duration-Limited Significant Wave Period (Finite Depth) ---
        # Similarly, this adapts the period formula.
        # Constants for depth influence are taken from fetch-limited depth-limited Ts.
        term_t_duration = 0.00379 * (dim_duration)**0.41 # From deep water duration-limited Ts
        tanh_depth_t = math.tanh(0.833 * (dim_depth)**0.375) # From fetch-limited depth-limited Ts
        Ts = (wind_speed_adjusted / gravity) * 7.54 * tanh_depth_t * math.tanh(term_t_duration / tanh_depth_t)

        # --- Equivalent Fetch Calculation (Finite Depth) ---
        # Equating the arguments of the inner tanh function for dimensionless Hs in duration-limited
        # and fetch-limited finite depth conditions:
        # 0.000528 * (g * t / U)^0.75 / tanh_depth_h = 0.00565 * (g * F / U^2)^0.5 / tanh_depth_h
        # This simplifies to:
        # 0.000528 * (g * t / U)^0.75 = 0.00565 * (g * F / U^2)^0.5
        # Solving for (g * F / U^2) which is dim_fetch:
        dim_fetch_equivalent = ( (0.000528 / 0.00565) * (dim_duration**0.75) )**2
        equivalent_fetch_m = dim_fetch_equivalent * (wind_speed_adjusted**2 / gravity)
        equivalent_fetch_km = equivalent_fetch_m / 1000 # Convert to kilometers

    return Hs, Ts, equivalent_fetch_km

def main():
    """
    Main function to drive the user interaction, calculations, and display results.
    """
    # Open the report file in write mode
    with open("report.txt", "w") as report_file:
        # Store original stdout
        original_stdout = sys.stdout
        # Redirect stdout to both the file and the console
        sys.stdout = Tee(sys.stdout, report_file)

        while True: # Loop for consecutive calculations
            print()
            print("======================================================")
            print("  Sverdrup-Munk-Bretschneider (SMB) Wave Calculator")
            print("======================================================")

            # --- Get User Inputs for Wind Speed (U10) ---
            while True:
                try:
                    print()
                    U10 = float(input("Enter wind speed at 10m height (U10) (m/s): "))
                    if U10 <= 0:
                        print("Wind speed (U10) must be a positive value.")
                    else:
                        break # Exit loop if input is a valid positive number.
                except ValueError:
                    print("Invalid input. Please enter a numeric value for wind speed (U10).")

            # --- Calculate Adjusted Wind Speed (Ua) ---
            Ua = calculate_adjusted_wind_speed(U10)
            print(f"Calculated Adjusted Wind Speed (Ua): {Ua:.2f} m/s")

            # --- Get User Inputs for Fetch and Duration ---
            # We will collect both fetch and duration to allow for comprehensive comparison.
            
            # Get Fetch Input
            while True:
                try:
                    fetch_km = float(input("Enter fetch length (km): "))
                    if fetch_km <= 0:
                        print("Fetch length must be a positive value.")
                    else:
                        fetch_m = fetch_km * 1000  # Convert fetch from km to meters
                        break
                except ValueError:
                    print("Invalid input. Please enter a numeric value for fetch.")

            # Get Duration Input
            while True:
                try:
                    duration_str = input("Enter storm duration (hours), or leave blank for effectively infinite duration: ").strip()
                    if duration_str == "":
                        duration_input_hours = float('inf') # Set to infinity if blank
                        break
                    duration_input_hours = float(duration_str)
                    if duration_input_hours <= 0:
                        print("Storm duration must be a positive value or left blank.")
                    else:
                        break
                except ValueError:
                    print("Invalid input. Please enter a numeric value for duration or leave blank.")

            # Get Depth Input (optional)
            while True:
                depth_str = input("Enter water depth (m) for finite depth, or leave blank for deep water: ").strip()
                if depth_str == "":
                    depth_input = None
                    break
                try:
                    depth_input = float(depth_str)
                    if depth_input < 0:
                        print("Depth cannot be negative. Please enter a positive value or leave blank.")
                    else:
                        break
                except ValueError:
                    print("Invalid input. Please enter a numeric value for depth or leave blank.")

            # --- Calculate for Fetch-Limited Scenario (based on input fetch) ---
            if depth_input is None:
                hs_fetch_potential, ts_fetch_potential, t_min_fetch_potential = calculate_deep_water(Ua, fetch_m)
            else:
                hs_fetch_potential, ts_fetch_potential, t_min_fetch_potential = calculate_depth_limited(Ua, fetch_m, depth_input)
            
            # --- Calculate for Duration-Limited Scenario (based on input duration) ---
            hs_duration_potential, ts_duration_potential, equivalent_fetch_duration_potential = calculate_duration_limited(Ua, duration_input_hours, depth_input)

            # --- Determine the Controlling Condition and Final Results ---
            # The actual wave height is the minimum of what fetch and duration allow.
            
            # First, determine the actual controlling Hs and Ts by taking the minimum of the two potentials
            # This ensures Hs never decreases with increased duration or fetch.
            if hs_fetch_potential <= hs_duration_potential:
                controlling_hs = hs_fetch_potential
                controlling_period = ts_fetch_potential
            else:
                controlling_hs = hs_duration_potential
                controlling_period = ts_duration_potential

            # Now, determine the controlling factor message and relevant values based on the physical constraints
            # This logic determines *why* the wave is limited (fetch or duration)
            if duration_input_hours < t_min_fetch_potential:
                # If the actual storm duration is less than the minimum required for the given fetch,
                # then duration is the primary physical limiting factor.
                controlling_factor_message = "Duration-Limited"
                relevant_duration_value = duration_input_hours
                relevant_duration_message = "Given Storm Duration"
                relevant_fetch_value = equivalent_fetch_duration_potential
                relevant_fetch_message = "Equivalent Fetch for Given Duration"
            else:
                # If the actual storm duration is greater than or equal to the minimum required,
                # then fetch is the primary physical limiting factor (or duration is sufficient).
                controlling_factor_message = "Fetch-Limited"
                relevant_duration_value = t_min_fetch_potential
                relevant_duration_message = "Minimum Duration for Given Fetch"
                relevant_fetch_value = fetch_km
                relevant_fetch_message = "Given Fetch"

            # --- Display Final Results ---
            print("\n------------------ RESULTS ------------------")
            print(f"Controlling Wave Growth Factor:     {controlling_factor_message}")
            print(f"Predicted Significant Wave Height (Hs): {controlling_hs:.2f} meters")
            print(f"Predicted Significant Wave Period (Ts): {controlling_period:.2f} seconds")
            
            if not math.isinf(relevant_duration_value) and not math.isnan(relevant_duration_value):
                print(f"{relevant_duration_message}:      {relevant_duration_value:.2f} hours")
            elif math.isinf(relevant_duration_value):
                print(f"{relevant_duration_message}:      Infinite hours (assumed)")
            
            if not math.isnan(relevant_fetch_value):
                print(f"{relevant_fetch_message}:                   {relevant_fetch_value:.2f} kilometers")
            
            print("-------------------------------------------")

            # --- Ask user if they want to perform another calculation ---
            print()
            another_calculation = input("Do you want to perform another calculation? (yes/no): ").lower()
            if another_calculation not in ['yes', 'y']:
                print("Exiting SMB Wave Calculator. Goodbye!")
                break # Exit the main loop
        
        # Restore original stdout
        sys.stdout = original_stdout

# Helper class to redirect print output to multiple destinations
class Tee(object):
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() # Ensure immediate writing
    def flush(self):
        for f in self.files:
            f.flush()

# This standard Python construct ensures that the `main()` function is called
# only when this script is executed directly. It prevents the code from running
# automatically if this script is imported as a module into another script.
if __name__ == "__main__":
    main()