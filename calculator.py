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
#   3. Asymptotic Consistency: The model enforces strictly consistent asymptotic 
#      limits for both fetch and duration-limited growth.
#   4. Minimum Duration: Uses the Hurdle & Stive (1989) power law (Eq. 4.12)
#      to determine the minimum time required for fetch-limited conditions.
#
# CALCULATION MODES:
#   1. Fetch-Limited: Uses Hurdle & Stive (1989) unified equations.
#   2. Duration-Limited: Uses the "Effective Fetch" method (Hurdle & Stive, 
#      1989, Eq 4.13). This inverts the minimum duration power law to find
#      an equivalent fetch, ensuring kinematic consistency.
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

import math
import sys
from typing import Tuple, Optional

# --- Constants ---
G_STD = 9.8066
PI = 3.14159265358979323846

# =============================================================================
#  SECTION 1: PHYSICS ENGINE
# =============================================================================

def calculate_adjusted_wind_speed(U10: float) -> float:
    """
    Calculates the adjusted wind speed (Ua) from the 10-meter wind speed (U10).
    Formula: Ua = 0.71 * U10^1.23 (SPM 1984)
    """
    return 0.71 * (U10**1.23)

def solve_dispersion(T: float, d: float, tol: float = 1e-15, max_iter: int = 100) -> float:
    """
    Solves the transcendental Linear Dispersion Relation using Newton-Raphson.
    
    Equation: L = (g * T^2 / 2pi) * tanh( 2pi * d / L )
    We solve for dimensionless wavenumber kh = k * d.
    
    Args:
        T (float): Wave Period (s)
        d (float): Water Depth (m)
    
    Returns:
        float: The dimensionless wavenumber kh = k*d
    """
    if d <= 0:
        return 0.0
    
    omega = 2.0 * PI / T
    k0 = (omega * omega) / G_STD
    k0h = k0 * d
    
    # Initial Guess (Carvalho 2006 explicit approximation)
    # This provides a very close starting point for Newton-Raphson
    if k0h > 0:
        kh = k0h / math.tanh((6.0/5.0)**k0h * math.sqrt(k0h))
    else:
        return 0.0

    # Newton-Raphson Iteration
    for _ in range(max_iter):
        t_kh = math.tanh(kh)
        f = k0h - kh * t_kh
        sech = 1.0 / math.cosh(kh)
        df = -t_kh - kh * (sech * sech)
        
        if abs(df) < 1e-15:
            break
        
        dkh = f / df
        kh_new = kh - dkh
        
        if abs(dkh / kh) < tol:
            return kh_new
        kh = kh_new
        
    return kh

def calculate_deep_water(Ua: float, fetch_m: float) -> Tuple[float, float, float, float]:
    """
    Calculates fetch-limited wave properties in deep water using Revised SMB.
    
    Returns:
        (Hs, Ts, t_min_hours, equivalent_fetch_km)
    """
    # Dimensionless Fetch: F_hat = g * F / Ua^2
    dim_fetch = (G_STD * fetch_m) / (Ua**2)

    # Hs Calculation (Deep Water Asymptote)
    term_fetch_h = 4.3e-5 * dim_fetch
    gHs_U2 = 0.25 * (math.tanh(term_fetch_h))**0.5
    Hs = gHs_U2 * (Ua**2 / G_STD)

    # Ts Calculation (Deep Water Asymptote)
    term_fetch_t = 4.1e-5 * dim_fetch
    gTs_U = 8.3 * (math.tanh(term_fetch_t))**(1/3)
    Ts = gTs_U * (Ua / G_STD)

    # Minimum Duration (t_min) - Hurdle & Stive (1989) Eq 4.12
    dim_duration = 65.9 * (dim_fetch**(2.0/3.0))
    t_min_seconds = dim_duration * Ua / G_STD
    t_min_hours = t_min_seconds / 3600.0

    return Hs, Ts, t_min_hours, 0.0

def calculate_depth_limited(Ua: float, fetch_m: float, depth: float) -> Tuple[float, float, float, float]:
    """
    Calculates fetch-limited wave properties for finite depth using Revised SMB.
    
    Returns:
        (Hs, Ts, t_min_hours, equivalent_fetch_km)
    """
    dim_fetch = (G_STD * fetch_m) / (Ua**2)
    dim_depth = (G_STD * depth) / (Ua**2)

    # Hs Calculation - Hurdle & Stive (1989) Eq 4.1
    depth_term_h = math.tanh(0.6 * dim_depth**0.75)
    if depth_term_h < 1e-6:
        depth_term_h = 1e-6

    fetch_term_inner_h = (4.3e-5 * dim_fetch) / (depth_term_h**2)
    gHs_U2 = 0.25 * depth_term_h * (math.tanh(fetch_term_inner_h))**0.5
    Hs = gHs_U2 * (Ua**2 / G_STD)

    # Ts Calculation - Hurdle & Stive (1989) Eq 4.2
    depth_term_t = math.tanh(0.76 * dim_depth**0.375)
    if depth_term_t < 1e-6:
        depth_term_t = 1e-6

    fetch_term_inner_t = (4.1e-5 * dim_fetch) / (depth_term_t**3)
    gTs_U = 8.3 * depth_term_t * (math.tanh(fetch_term_inner_t))**(1/3)
    Ts = gTs_U * (Ua / G_STD)

    # Minimum Duration (t_min) - Recommends deep water power law for consistency
    dim_duration = 65.9 * (dim_fetch**(2.0/3.0))
    t_min_seconds = dim_duration * Ua / G_STD
    t_min_hours = t_min_seconds / 3600.0

    return Hs, Ts, t_min_hours, 0.0

def calculate_duration_limited(Ua: float, duration_hours: float, depth: Optional[float] = None) -> Tuple[float, float, float, float]:
    """
    Calculates wave properties for duration-limited conditions using Effective Fetch.
    
    Returns:
        (Hs, Ts, t_min_hours, equivalent_fetch_km)
    """
    duration_seconds = duration_hours * 3600.0
    
    # 1. Dimensionless Duration: t_hat = g * t / Ua
    dim_duration = (G_STD * duration_seconds) / Ua
    
    # 2. Effective Fetch (Eq 4.13): F_hat' = (t_hat / 65.9)^(1.5)
    dim_fetch_equivalent = (dim_duration / 65.9)**1.5
    
    # Convert back to meters
    equivalent_fetch_m = dim_fetch_equivalent * (Ua**2 / G_STD)
    equivalent_fetch_km = equivalent_fetch_m / 1000.0

    # 3. Calculate Waves using Effective Fetch
    if depth is None:
        Hs, Ts, _, _ = calculate_deep_water(Ua, equivalent_fetch_m)
    else:
        Hs, Ts, _, _ = calculate_depth_limited(Ua, equivalent_fetch_m, depth)

    # Return calculated values, setting t_min to the input duration for this branch
    return Hs, Ts, duration_hours, equivalent_fetch_km

# =============================================================================
#  SECTION 2: REPORTING & UTILITIES
# =============================================================================

class Tee(object):
    """Helper to redirect output to both console and file."""
    def __init__(self, *files):
        self.files = files
    def write(self, obj):
        for f in self.files:
            f.write(obj)
            f.flush() 
    def flush(self):
        for f in self.files:
            f.flush()

def format_row(label: str, value, unit: str = "", width_lbl: int = 25, width_val: int = 12) -> str:
    if isinstance(value, float):
        val_str = f"{value:.4f}" if abs(value) > 1e-9 else "0.00"
    else:
        val_str = str(value)
    
    return f"  {label:<{width_lbl}} : {val_str:<{width_val}}{unit}"

def run_simulation(U10: float, fetch_km: float, duration_hr: Optional[float], depth_m: Optional[float]):
    """
    Orchestrates the physics calculation and prints the detailed report.
    """
    # --- 1. Physics Calculations ---
    Ua = calculate_adjusted_wind_speed(U10)
    fetch_m = fetch_km * 1000.0
    
    is_deep_water = (depth_m is None)
    is_inf_duration = (duration_hr is None or math.isinf(duration_hr))

    # A. Calculate Potential A: Fetch Limited
    if is_deep_water:
        r_fetch = calculate_deep_water(Ua, fetch_m) # (Hs, Ts, t_min, eq_f)
    else:
        r_fetch = calculate_depth_limited(Ua, fetch_m, depth_m)

    # B. Calculate Potential B: Duration Limited
    dur_val = 1e9 if is_inf_duration else duration_hr
    r_dur = calculate_duration_limited(Ua, dur_val, depth_m)

    # C. Determine Controlling Factor
    # Comparing Hs to find the limiting state (Minimum of the two potentials)
    if r_fetch[0] <= r_dur[0]:
        # Fetch Limited
        final_Hs, final_Ts, final_t_min, _ = r_fetch
        control_msg = "FETCH-LIMITED"
        
        relevant_dur_val = final_t_min
        relevant_dur_label = "Min Duration Req."
        relevant_fetch_val = fetch_km
        relevant_fetch_label = "Given Fetch"
    else:
        # Duration Limited
        final_Hs, final_Ts, _, final_eq_fetch = r_dur
        control_msg = "DURATION-LIMITED"
        
        relevant_dur_val = duration_hr if not is_inf_duration else "Infinite"
        relevant_dur_label = "Given Duration"
        relevant_fetch_val = final_eq_fetch
        relevant_fetch_label = "Equivalent Fetch"

    # D. Precise Hydrodynamics (Dispersion & Celerity)
    # This section replicates the "SMBEngine::solve_dispersion" logic
    if is_deep_water:
        # Deep water analytic: L0 = gT^2 / 2pi
        L = (G_STD * final_Ts**2) / (2 * PI)
        k = 2 * PI / L
        kh = 100.0 # Effectively infinite
        d_actual = 1e9
    else:
        d_actual = depth_m
        kh = solve_dispersion(final_Ts, d_actual)
        k = kh / d_actual
        L = 2 * PI / k
        
    C = L / final_Ts # Celerity

    # --- 2. Report Generation ---
    
    print("============================================================")
    print("  SMB WAVE PREDICTION REPORT")
    print("============================================================")
    print("  Methodology: Revised SMB")
    print("  Ref: Hurdle & Stive (1989) Unified Formulations")
    
    print("\n------------------------------------------------------------")
    print("  1. INPUT PARAMETERS")
    print("------------------------------------------------------------")
    print(format_row("Wind Speed (U10)", U10, "m/s"))
    print(format_row("Adjusted Speed (Ua)", Ua, "m/s"))
    print(format_row("Fetch Length (F)", fetch_km, "km"))
    
    if is_inf_duration:
        print(format_row("Storm Duration (t_act)", "Infinite", "hours"))
    else:
        print(format_row("Storm Duration (t_act)", duration_hr, "hours"))
        
    if is_deep_water:
        print(format_row("Water Depth (d)", "Deep Water", "-"))
    else:
        print(format_row("Water Depth (d)", depth_m, "m"))

    print("\n------------------------------------------------------------")
    print("  2. LIMITING CONDITIONS")
    print("------------------------------------------------------------")
    print(format_row(">> Controlling Factor", control_msg))
    print("  >> Determining Parameters:")
    
    print(format_row("   " + relevant_dur_label, relevant_dur_val, "hours"))
    print(format_row("   " + relevant_fetch_label, relevant_fetch_val, "km"))

    print("\n------------------------------------------------------------")
    print("  3. PREDICTION RESULTS")
    print("------------------------------------------------------------")
    print(format_row("Sig. Wave Height (Hs)", final_Hs, "m"))
    print(format_row("Sig. Wave Period (Ts)", final_Ts, "s"))
    print(format_row("Wave Celerity (C)", C, "m/s"))
    print(format_row("Wave Length (L)", L, "m"))
    print(format_row("Wave Number (k)", k, "rad/m"))

    print("\n------------------------------------------------------------")
    print("  4. INSIGHTS & OBSERVATIONS")
    print("------------------------------------------------------------")
    
    # Insight A: Growth State
    if control_msg == "FETCH-LIMITED":
        print("  >> State: Fully Developed for Fetch.")
        print("     Waves have reached maximum size for this distance.")
    else:
        print("  >> State: Growing Sea State (Duration Limited).")
        print("     Waves are still growing over time.")

    # Insight B: Regime & Stability
    if is_deep_water:
        regime = "Deep Water (d/L > 0.5)"
    else:
        d_L = d_actual / L
        if d_L >= 0.5: regime = "Deep Water (d/L > 0.5)"
        elif d_L < 0.05: regime = "Shallow Water (d/L < 0.05)"
        else: regime = "Transitional/Intermediate"
    
    print(format_row(">> Flow Regime", regime))

    # Breaking (Miche Criterion)
    # H_max / L = 0.142 * tanh(kd)
    steepness = final_Hs / L
    limit_steepness = 0.142 * math.tanh(kh)
    
    print(format_row("   Calculated Steepness", steepness, "-"))
    print(format_row("   Miche Limit (H/L)", limit_steepness, "-"))
    
    if steepness > limit_steepness:
        print(format_row(">> Breaking Status", "UNSTABLE / BREAKING [!]", ""))
        print("     Wave height exceeds the Miche stability limit.")
    else:
        print(format_row(">> Breaking Status", "STABLE", ""))
        safety_margin = (1.0 - steepness/limit_steepness) * 100.0
        print(format_row("   Stability Margin", safety_margin, "%"))

    if not is_deep_water:
        ratio = final_Hs / d_actual
        print(format_row("   Depth Ratio (Hs/d)", ratio, "-"))

    print("\n------------------------------------------------------------")
    print("  5. SCENARIO ANALYSIS")
    print("------------------------------------------------------------")
    print("  [A] Fetch-Limited Potential")
    print(format_row("    - Potential Hs", r_fetch[0], "m"))
    print(format_row("    - Potential Ts", r_fetch[1], "s"))
    print(format_row("    - Time to Develop", r_fetch[2], "hours"))
    
    print()
    print("  [B] Duration-Limited Potential")
    if is_inf_duration:
        print("      - Not applicable (Duration is infinite)")
    else:
        print(format_row("    - Potential Hs", r_dur[0], "m"))
        print(format_row("    - Potential Ts", r_dur[1], "s"))
        print(format_row("    - Equiv. Fetch", r_dur[3], "km"))
    print()

def main():
    """
    Main entry point handling user input and single execution.
    """
    # --- 1. CLI Header & Input Phase (Console Only) ---
    print("======================================================")
    print("  SMB WAVE PREDICTION (CLI VERSION)")
    print("======================================================")

    try:
        # 1. Wind Speed
        while True:
            u_str = input("Enter Wind Speed (U10) [m/s]: ").strip()
            if u_str: 
                u10 = float(u_str)
                if u10 > 0: break
            print(">> Error: Wind speed must be > 0.")

        # 2. Fetch
        while True:
            f_str = input("Enter Fetch Length (F) [km]:  ").strip()
            if f_str:
                fetch = float(f_str)
                if fetch > 0: break
            print(">> Error: Fetch must be > 0.")

        # 3. Duration
        d_str = input("Enter Duration (hours) [Blank=Inf]: ").strip()
        duration = float(d_str) if d_str else None
        if duration is not None and duration <= 0:
            print(">> Note: Negative/Zero duration treated as infinite.")
            duration = None

        # 4. Depth
        while True:
            dp_str = input("Enter Depth (meters) [Blank=Deep]:  ").strip()
            if not dp_str:
                depth = None
                break
            try:
                depth = float(dp_str)
                if depth > 0: break
                print(">> Error: Depth must be > 0.")
            except ValueError:
                print(">> Error: Please enter a valid number.")

    except ValueError:
        print("\n>> Input Error: Please enter valid numeric values.")
        return
    except KeyboardInterrupt:
        print("\nOperation cancelled by user.")
        return

    # --- 2. Execution & Reporting Phase (Console + File) ---
    with open("report.txt", "w") as report_file:
        original_stdout = sys.stdout
        try:
            sys.stdout = Tee(original_stdout, report_file)
            run_simulation(u10, fetch, duration, depth)
        finally:
            sys.stdout = original_stdout

if __name__ == "__main__":
    main()