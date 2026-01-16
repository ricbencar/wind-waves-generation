# =============================================================================
# SMB Wave Prediction Model
# =============================================================================
#
# DESCRIPTION:
# This script implements the Sverdrup-Munk-Bretschneider (SMB) method for
# predicting wind-generated wave characteristics.
#
# UPDATED FORMULATION (Hurdle & Stive, 1989):
# The core equations for Significant Wave Height (Hs) and Period (Ts) have been
# updated to the unified formulations proposed by Hurdle & Stive (1989). These
# equations replace the disjointed deep/shallow water formulas from the original
# SPM (1984), eliminating step-changes (discontinuities) at transition points.
#
# KEY FEATURES:
#   1. Revised Coefficients: Uses 0.25 (for Hs) and 8.3 (for Ts) instead of
#      the traditional 0.283 and 7.54, utilizing Hurdle & Stive's adjusted
#      hyperbolic structures for depth attenuation.
#   2. Unified Logic: The same equation structure applies to both deep and
#      finite-depth water, ensuring mathematical consistency across domains.
#   3. Asymptotic Consistency: Duration-limited formulas have been adapted with
#      the revised coefficients to ensure they asymptote to the same
#      fully-developed limits as the fetch-limited equations.
#   4. Minimum Duration: Retains the modern approximation by Etemad-Shahidi
#      et al. (2009) for calculating t_min.
#
# CALCULATION MODES:
#   1. Fetch-Limited: Uses Hurdle & Stive (1989) unified equations.
#   2. Duration-Limited: Uses CEM/SPM time-growth equations adapted with
#      Hurdle & Stive coefficients (0.25/8.3) to maintain model consistency.
#
# BIBLIOGRAPHY:
#
# Primary Engineering Manuals (Operational Standards)
# * U.S. Army. (2008). Coastal Engineering Manual (EM 1110-2-1100). Washington,
#   DC: U.S. Army Corps of Engineers.
#   - The current primary reference for USACE coastal projects, superseding
#     the Shore Protection Manual.
#
# * U.S. Army. (1984). Shore Protection Manual (Vol. 1 & 2). Vicksburg, MS:
#   U.S. Army Engineer Waterways Experiment Station.
#   - The classic reference that standardized the SMB equations for decades.
#
# * World Meteorological Organization (WMO). (2018). Guide to Wave Analysis and
#   Forecasting (WMO-No. 702). Geneva: Secretariat of the WMO.
#   - The international standard for meteorological wave forecasting.
#
# Foundational Papers (The "SMB" Method)
# * Sverdrup, H. U., & Munk, W. H. (1947). Wind, Sea, and Swell: Theory of
#   Relations for Forecasting. H.O. Pub. No. 601, U.S. Navy Hydrographic Office,
#   Washington, D.C.
#   - The original wartime research that established the "Significant Wave" concept.
#
# * Bretschneider, C. L. (1958). Revisions in Wave Forecasting: Deep and Shallow
#   Water. Proceedings of the 6th Conference on Coastal Engineering, pp. 30–67.
#   - The paper that revised the 1947 Sverdrup-Munk curves, adding the "B" to SMB.
#
# * Bretschneider, C. L. (1970). Wave forecasting relations for wave generation.
#   Look Lab, Hawaii, 1(3).
#   - Further refinement of the 1958 curves.
#
# Revisions, Critiques & Finite Depth Extensions
# * Hurdle, D. P., & Stive, R. (1989). Revision of SPM 1984 wave hindcast model
#   to avoid inconsistencies in engineering applications. Coastal Engineering,
#   12(4), 339–351.
#   - Corrects inconsistencies in the original 1984 SPM formulas.
#
# * Bishop, C. T., Donelan, M. A., & Kahma, K. K. (1992). Shore protection
#   manual's wave prediction reviewed. Coastal Engineering, 17(1-2), 25-48.
#   - Comprehensive comparison of SPM predictions against measured data.
#
# * Etemad-Shahidi, A., Kazeminezhad, M. H., & Mousavi, S. J. (2009). On the
#   prediction of wave parameters using simplified methods. Journal of Coastal
#   Research, SI 56, 505-509.
#   - Validates the SMB equations against modern methods.
#
# =============================================================================
import math
import numpy as np
import sys # Import sys module to redirect stdout

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
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours for practical use

    return Hs, Ts, t_min_hours

def calculate_depth_limited(wind_speed_adjusted, fetch, depth, gravity=9.8066):
    """
    Calculates wave properties for depth-limited conditions using the Revised SMB method.

    This function implements the unified equations from Hurdle & Stive (1989),
    which correct the inconsistencies found in the original SPM (1984) formulations
    at the transition between deep and shallow water.

    Args:
        wind_speed_adjusted (float): Adjusted wind speed (Ua) (m/s).
        fetch (float):      Effective fetch length (m).
        depth (float):      Water depth (m).
        gravity (float):    Acceleration of gravity (m/s^2).

    Returns:
        tuple: A tuple containing Hs (m), Ts (s), and t_min (hours).
    """
    # --- Dimensionless Parameters ---
    dim_fetch = (gravity * fetch) / wind_speed_adjusted**2
    dim_depth = (gravity * depth) / wind_speed_adjusted**2

    # --- Revised Significant Wave Height (Hs) ---
    # Hurdle & Stive (1989), Eq 4.1
    # Coefficient is 0.25 (vs 0.283 in SPM)
    # Depth term: tanh(0.6 * d_hat^0.75)
    # Fetch term inner: 4.3e-5 * F_hat / (depth_term^2)
    # Exponent: 0.5
    depth_term_h = math.tanh(0.6 * dim_depth**0.75)
    fetch_term_inner_h = (4.3e-5 * dim_fetch) / (depth_term_h**2)
    
    gHs_U2 = 0.25 * depth_term_h * (math.tanh(fetch_term_inner_h))**0.5
    Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

    # --- Revised Significant Wave Period (Ts) ---
    # Hurdle & Stive (1989), Eq 4.2
    # Coefficient is 8.3 (vs 7.54 in SPM)
    # Depth term: tanh(0.76 * d_hat^0.375)
    # Fetch term inner: 4.1e-5 * F_hat / (depth_term^3)
    # Exponent: 1/3
    depth_term_t = math.tanh(0.76 * dim_depth**0.375)
    fetch_term_inner_t = (4.1e-5 * dim_fetch) / (depth_term_t**3)
    
    gTs_U = 8.3 * depth_term_t * (math.tanh(fetch_term_inner_t))**(1/3)
    Ts = gTs_U * (wind_speed_adjusted / gravity)

    # --- Minimum Wind Duration (t_min) Calculation ---
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours for practical use

    return Hs, Ts, t_min_hours

def calculate_duration_limited(wind_speed_adjusted, duration_hours, depth=None, gravity=9.8066):
    """
    Calculates wave properties for duration-limited conditions, considering finite depth if provided.

    UPDATED COEFFICIENTS (Hurdle & Stive, 1989 Compatibility):
    This function has been updated to use the leading coefficients (0.25, 8.3)
    and depth-damping factors (0.6, 0.76) from the Hurdle & Stive (1989) revision.
    This ensures that as time -> infinity, the wave parameters asymptote to the
    same fully-developed limits as the revised fetch-limited model, preventing
    inconsistencies.

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
    # Formula: t_hat = g * t / U
    dim_duration = (gravity * duration_seconds) / wind_speed_adjusted

    if depth is None or depth == 0: # Treat as deep water if depth is not provided or zero
        # --- Significant Wave Height (Hs) - Deep Water ---
        # Coefficient updated to 0.25 (Hurdle & Stive asymptote)
        # Note: Retaining CEM time-growth exponent 0.75
        gHs_U2 = 0.25 * math.tanh(0.000528 * (dim_duration**0.75))
        Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

        # --- Significant Wave Period (Ts) - Deep Water ---
        # Coefficient updated to 8.3 (Hurdle & Stive asymptote)
        # Note: Retaining CEM time-growth exponent 0.41
        gTs_U = 8.3 * math.tanh(0.00379 * (dim_duration**0.41))
        Ts = gTs_U * (wind_speed_adjusted / gravity)

        # --- Equivalent Fetch Calculation (Deep Water) ---
        # Calculated using the Hurdle & Stive deep water fetch limit: 0.25 * (tanh(4.3e-5 * F_hat))^0.5
        # Equating H_hat(duration) = H_hat(fetch):
        # tanh(0.000528 * t_hat^0.75) = (tanh(4.3e-5 * F_hat))^0.5
        # F_hat = atanh( tanh(time_term)^2 ) / 4.3e-5
        # Note: This is an approximation to map the new coefficients.
        try:
            time_term = 0.000528 * (dim_duration**0.75)
            # Prevent domain error if time_term is large (tanh -> 1)
            if time_term > 10: 
                tanh_sq = 1.0
            else:
                tanh_sq = math.tanh(time_term)**2
            
            # Avoid singularity at fully developed state
            if tanh_sq >= 1.0:
                dim_fetch_equivalent = 200000 # Large number representing fully developed
            else:
                dim_fetch_equivalent = math.atanh(tanh_sq) / 4.3e-5
                
            equivalent_fetch_m = dim_fetch_equivalent * (wind_speed_adjusted**2 / gravity)
            equivalent_fetch_km = equivalent_fetch_m / 1000
        except ValueError:
            equivalent_fetch_km = 0.0

    else:
        # --- Dimensionless Depth Calculation ---
        dim_depth = (gravity * depth) / wind_speed_adjusted**2

        # --- Heuristic for Duration-Limited Hs (Finite Depth) ---
        # Updated Depth Coefficient: 0.530 -> 0.6 (Hurdle & Stive Eq 4.1)
        # Updated Leading Coefficient: 0.283 -> 0.25
        term_h_duration = 0.000528 * (dim_duration)**0.75
        tanh_depth_h = math.tanh(0.6 * (dim_depth)**0.75) 
        Hs = (wind_speed_adjusted**2 / gravity) * 0.25 * tanh_depth_h * math.tanh(term_h_duration / tanh_depth_h)

        # --- Heuristic for Duration-Limited Ts (Finite Depth) ---
        # Updated Depth Coefficient: 0.833 -> 0.76 (Hurdle & Stive Eq 4.2)
        # Updated Leading Coefficient: 7.54 -> 8.3
        term_t_duration = 0.00379 * (dim_duration)**0.41
        tanh_depth_t = math.tanh(0.76 * (dim_depth)**0.375)
        Ts = (wind_speed_adjusted / gravity) * 8.3 * tanh_depth_t * math.tanh(term_t_duration / tanh_depth_t)

        # --- Equivalent Fetch Calculation (Finite Depth) ---
        # Simplified mapping based on the updated coefficients
        # 0.000528 * t_hat^0.75 ~ 4.3e-5 * F_hat (Approximation of the inner arguments)
        dim_fetch_equivalent = (0.000528 / 4.3e-5) * (dim_duration**0.75)
        equivalent_fetch_m = dim_fetch_equivalent * (wind_speed_adjusted**2 / gravity)
        equivalent_fetch_km = equivalent_fetch_m / 1000

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
            print("  Revised Model (Hurdle & Stive, 1989)")
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
                    duration_str = input("Enter storm duration (hours), or leave blank for infinite duration: ").strip()
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