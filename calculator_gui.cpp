/**
 * ==============================================================================
 * SMB WAVE PREDICTION CALCULATOR & GUI
 * ==============================================================================
 * MODULE:   calculator_gui.cpp
 * TYPE:     Empirical Wave Prediction Tool (Win32 GUI)
 * METHOD:   Sverdrup-Munk-Bretschneider (SMB) + Exact Linear Dispersion
 * REF:      U.S. Army Coastal Engineering Manual (CEM) & Shore Protection Manual (SPM)
 * LICENSE:  MIT / Academic Open Source
 * ==============================================================================
 *
 * ------------------------------------------------------------------------------
 * 1. THEORETICAL MANUAL & PHYSICS ENGINE DOCUMENTATION
 * ------------------------------------------------------------------------------
 *
 * 1.1 PHYSICAL OVERVIEW
 * This software predicts the characteristics of wind-generated waves (Significant
 * Wave Height Hs and Period Ts) based on the Sverdrup-Munk-Bretschneider (SMB)
 * method.
 *
 * 1.2 ATMOSPHERIC FORCING (WIND STRESS)
 * The driving force for wave generation is the wind stress on the water surface.
 * Standard meteorological wind speed is measured at 10m elevation (U10).
 * However, the empirical curves in the SPM (1984) were calibrated using an
 * "Adjusted Wind Speed" (Ua), also known as the wind stress factor.
 *
 * FORMULA: Ua = 0.71 * (U10)^1.23
 * - U10: Wind speed at 10m elevation [m/s]
 * - Ua:  Adjusted wind speed [m/s]
 *
 * 1.3 DIMENSIONAL ANALYSIS & PARAMETERS
 * The governing physics are described by the following dimensionless groups:
 *
 * - Dimensionless Fetch (F_hat):    F_hat = (g * F) / Ua^2
 * - Dimensionless Depth (d_hat):    d_hat = (g * d) / Ua^2
 * - Dimensionless Duration (t_hat): t_hat = (g * t) / Ua
 *
 * Where:
 * - g = Gravitational acceleration (9.8066 m/s^2)
 * - F = Fetch Length [m]
 * - d = Water Depth [m]
 * - t = Duration [s]
 *
 * 1.4 WAVE GROWTH EQUATIONS (FETCH-LIMITED - Hurdle & Stive, 1989)
 * The revised equations provide a single unified formulation for both deep and
 * finite-depth water, avoiding inconsistencies found in the original SPM.
 *
 * A. Significant Wave Height (Hs):
 * (g * Hs) / Ua^2 = 0.25 * tanh[K_d1] * tanh^0.5 [ (4.3e-5 * F_hat) / tanh^2(K_d1) ]
 * Where K_d1 = 0.6 * (d_hat)^0.75
 *
 * B. Significant Wave Period (Ts):
 * (g * Ts) / Ua   = 8.3 * tanh[K_d2] * tanh^1/3 [ (4.1e-5 * F_hat) / tanh^3(K_d2) ]
 * Where K_d2 = 0.76 * (d_hat)^0.375
 *
 * Note: In Deep Water (d_hat -> inf), the depth terms tanh(K_d) approach 1.0.
 *
 * 1.5 WAVE GROWTH EQUATIONS (DURATION-LIMITED)
 * When the storm duration (t) is the limiting factor, the equations are adapted
 * with the updated coefficients (0.25 and 8.3) to ensure asymptotic consistency.
 *
 * (g * Hs) / Ua^2 = 0.25 * tanh[ 0.000528 * (t_hat)^0.75 ]
 * (g * Ts) / Ua   = 8.30 * tanh[ 0.003790 * (t_hat)^0.41 ]
 *
 * 1.6 EXACT WAVE KINEMATICS (DISPERSION RELATION)
 * To convert the Period (Ts) into Wavelength (L), we solve the Linear Dispersion
 * Relation for water waves. This is a transcendental equation that cannot be
 * solved analytically for intermediate depths.
 *
 * Equation: L = (g * T^2 / 2pi) * tanh( 2pi * d / L )
 *
 * Numerical Solution (Newton-Raphson):
 * We solve for the dimensionless wavenumber 'kh' = k * d.
 * Function: f(kh) = (w^2 * d / g) - kh * tanh(kh) = 0
 * Derivative: f'(kh) = -tanh(kh) - kh * sech^2(kh)
 *
 * 1.7 STABILITY CRITERIA (BREAKING)
 * The model checks if the predicted waves are physically sustainable.
 *
 * A. Miche Criterion (Steepness Breaking):
 * Waves break when the particle velocity at the crest exceeds the phase speed.
 * Limit: H / L <= 0.142 * tanh(k * d)
 * In deep water, this simplifies to H/L <= 1/7.
 *
 * ------------------------------------------------------------------------------
 * COMPILATION INSTRUCTIONS:
 * g++ -O3 -std=c++17 -static -static-libgcc -static-libstdc++ 
 * -o calculator_gui.exe calculator_gui.cpp -mwindows -lgdi32
 * ==============================================================================
 *
 * BIBLIOGRAPHY:
 *
 * Primary Engineering Manuals (Operational Standards)
 * * U.S. Army. (2008). Coastal Engineering Manual (EM 1110-2-1100). Washington,
 * DC: U.S. Army Corps of Engineers.
 * - The current primary reference for USACE coastal projects, superseding
 * the Shore Protection Manual.
 *
 * * U.S. Army. (1984). Shore Protection Manual (Vol. 1 & 2). Vicksburg, MS:
 * U.S. Army Engineer Waterways Experiment Station.
 * - The classic reference that standardized the SMB equations for decades.
 *
 * * World Meteorological Organization (WMO). (2018). Guide to Wave Analysis and
 * Forecasting (WMO-No. 702). Geneva: Secretariat of the WMO.
 * - The international standard for meteorological wave forecasting.
 *
 * Foundational Papers (The "SMB" Method)
 * * Sverdrup, H. U., & Munk, W. H. (1947). Wind, Sea, and Swell: Theory of
 * Relations for Forecasting. H.O. Pub. No. 601, U.S. Navy Hydrographic Office,
 * Washington, D.C.
 * - The original wartime research that established the "Significant Wave" concept.
 *
 * * Bretschneider, C. L. (1958). Revisions in Wave Forecasting: Deep and Shallow
 * Water. Proceedings of the 6th Conference on Coastal Engineering, pp. 30–67.
 * - The paper that revised the 1947 Sverdrup-Munk curves, adding the "B" to SMB.
 *
 * * Bretschneider, C. L. (1970). Wave forecasting relations for wave generation.
 * Look Lab, Hawaii, 1(3).
 * - Further refinement of the 1958 curves.
 *
 * Revisions, Critiques & Finite Depth Extensions
 * * Hurdle, D. P., & Stive, R. (1989). Revision of SPM 1984 wave hindcast model
 * to avoid inconsistencies in engineering applications. Coastal Engineering,
 * 12(4), 339–351.
 * - Corrects inconsistencies in the original 1984 SPM formulas.
 *
 * * Bishop, C. T., Donelan, M. A., & Kahma, K. K. (1992). Shore protection
 * manual's wave prediction reviewed. Coastal Engineering, 17(1-2), 25-48.
 * - Comprehensive comparison of SPM predictions against measured data.
 *
 * * Etemad-Shahidi, A., Kazeminezhad, M. H., & Mousavi, S. J. (2009). On the
 * prediction of wave parameters using simplified methods. Journal of Coastal
 * Research, SI 56, 505-509.
 * - Validates the SMB equations against modern methods.
 */

#define _USE_MATH_DEFINES 
#define NOMINMAX

#include <windows.h>
#include <commctrl.h>
#include <string>
#include <sstream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fstream>
#include <limits>

// ==============================================================================
//  SECTION 1: CORE TYPES & CONSTANTS
// ==============================================================================

namespace Core {
    using Real    = double;
    using String  = std::string;
}

using namespace Core;

namespace Phys {
    constexpr Real G_STD = 9.8066 ; // Standard Gravity [m/s^2]
    constexpr Real PI    = 3.14159265358979323846;
}

// ==============================================================================
//  SECTION 2: LOGGING UTILITIES
//  Handles the formatting of the text report displayed in the output window.
// ==============================================================================

class StringLogger {
    std::stringstream buffer;
    const int width_label = 25;
    const int width_val   = 12;
    const int width_unit  = 8;
    const int line_len    = 60;

public:
    // Prints a major section header with double borders
    void print_header(const String& title) {
        buffer << "\n" << std::string(line_len, '=') << "\n";
        buffer << "  " << title << "\n";
        buffer << std::string(line_len, '=') << "\n";
    }

    // Prints a sub-section header with single borders
    void print_subheader(const String& title) {
        buffer << "\n" << std::string(line_len, '-') << "\n";
        buffer << "  " << title << "\n";
        buffer << std::string(line_len, '-') << "\n";
    }
    
    // Prints a standard data row: "Label : Value Unit"
    void print_row(const String& desc, Real val, const String& unit) {
        std::stringstream ss; 
        if (std::abs(val) < 1e-9) ss << "0.00"; 
        else ss << std::fixed << std::setprecision(4) << val;
        
        buffer << "  " 
               << std::left << std::setw(width_label) << desc << " : " 
               << std::left << std::setw(width_val)   << ss.str()
               << std::left << std::setw(width_unit)  << unit << "\n";
    }

    // Prints a string data row: "Label : Value"
    void print_str(const String& desc, const String& val, const String& unit = "") {
        buffer << "  " 
               << std::left << std::setw(width_label) << desc << " : " 
               << std::left << std::setw(width_val)   << val
               << std::left << std::setw(width_unit)  << unit << "\n";
    }

    // Prints a plain text line (indented)
    void print_text(const String& text) {
        buffer << "  " << text << "\n";
    }

    // Inserts a blank line
    void newline() { buffer << "\n"; }

    // Returns the complete report string
    std::string get_content() const { return buffer.str(); }
};

// ==============================================================================
//  SECTION 3: PHYSICS ENGINE (SMB & DISPERSION)
// ==============================================================================

struct SMBResult {
    Real Hs;          // Significant Wave Height [m]
    Real Ts;          // Significant Wave Period [s]
    Real t_min;       // Minimum Duration required for fetch-limited state [Hours]
    Real equiv_fetch; // Equivalent Fetch length [km] (for duration-limited cases)
};

class SMBEngine {
public:
    // --- 3.1 Empirical Formulas (Revised SMB - Hurdle & Stive, 1989) ---

    /**
     * Calculates Adjusted Wind Speed (Ua).
     * Correction for the non-linear relationship between wind speed and stress.
     * Formula: Ua = 0.71 * U10^1.23
     */
    static Real calculate_adjusted_wind_speed(Real U10) {
        return 0.71 * std::pow(U10, 1.23);
    }

    /**
     * Calculates Fetch-Limited parameters for Deep Water.
     * Uses the asymptotic limits of the Hurdle & Stive (1989) unified equations.
     *
     * @param Ua: Adjusted wind speed [m/s]
     * @param fetch_m: Fetch length [m]
     */
    static SMBResult calculate_deep_water(Real Ua, Real fetch_m) {
        // Dimensionless Fetch: F_hat = g * F / Ua^2
        Real dim_fetch = (Phys::G_STD * fetch_m) / (Ua * Ua);

        // Revised Significant Wave Height (Hs) - Deep Water Asymptote
        // Formula: (g * Hs) / Ua^2 = 0.25 * [ tanh(4.3e-5 * F_hat) ]^0.5
        Real term_fetch_h = 4.3e-5 * dim_fetch;
        Real gHs_U2 = 0.25 * std::pow(std::tanh(term_fetch_h), 0.5);
        Real Hs = gHs_U2 * (Ua * Ua / Phys::G_STD);

        // Revised Significant Wave Period (Ts) - Deep Water Asymptote
        // Formula: (g * Ts) / Ua = 8.3 * [ tanh(4.1e-5 * F_hat) ]^(1/3)
        Real term_fetch_t = 4.1e-5 * dim_fetch;
        Real gTs_U = 8.3 * std::pow(std::tanh(term_fetch_t), 1.0/3.0);
        Real Ts = gTs_U * (Ua / Phys::G_STD);

        // Minimum Duration (t_min)
        // Solved via approximation from Etemad-Shahidi et al. (2009)
        Real log_dim_fetch = std::log(dim_fetch);
        Real A = 0.0161, B = 0.3692, C = 2.2024, D = 0.8798;
        Real exponent_term = std::pow(A * std::pow(log_dim_fetch, 2) - B * log_dim_fetch + C, 0.5) + D * log_dim_fetch;
        Real gt_min_U = 6.5882 * std::exp(exponent_term);
        Real t_min_hr = (gt_min_U * Ua / Phys::G_STD) / 3600.0;

        return {Hs, Ts, t_min_hr, 0.0};
    }

    /**
     * Calculates Fetch-Limited parameters for Finite Depth.
     * Uses the unified equations from Hurdle & Stive (1989).
     *
     * @param Ua: Adjusted wind speed [m/s]
     * @param fetch_m: Fetch length [m]
     * @param depth: Water depth [m]
     */
    static SMBResult calculate_depth_limited(Real Ua, Real fetch_m, Real depth) {
        Real dim_fetch = (Phys::G_STD * fetch_m) / (Ua * Ua);
        Real dim_depth = (Phys::G_STD * depth) / (Ua * Ua);

        // Revised Significant Wave Height (Hs) - Hurdle & Stive (1989) Eq 4.1
        // Depth-damping: K_d1 = tanh(0.6 * d_hat^0.75)
        Real depth_term_h = std::tanh(0.6 * std::pow(dim_depth, 0.75));
        Real fetch_term_inner_h = (4.3e-5 * dim_fetch) / (depth_term_h * depth_term_h);
        
        Real gHs_U2 = 0.25 * depth_term_h * std::pow(std::tanh(fetch_term_inner_h), 0.5);
        Real Hs = gHs_U2 * (Ua * Ua / Phys::G_STD);

        // Revised Significant Wave Period (Ts) - Hurdle & Stive (1989) Eq 4.2
        // Depth-damping: K_d2 = tanh(0.76 * d_hat^0.375)
        Real depth_term_t = std::tanh(0.76 * std::pow(dim_depth, 0.375));
        Real fetch_term_inner_t = (4.1e-5 * dim_fetch) / std::pow(depth_term_t, 3.0);
        
        Real gTs_U = 8.3 * depth_term_t * std::pow(std::tanh(fetch_term_inner_t), 1.0/3.0);
        Real Ts = gTs_U * (Ua / Phys::G_STD);

        // Minimum Duration (t_min) - Etemad-Shahidi et al. (2009)
        Real log_dim_fetch = std::log(dim_fetch);
        Real A = 0.0161, B = 0.3692, C = 2.2024, D = 0.8798;
        Real exponent_term = std::pow(A * std::pow(log_dim_fetch, 2) - B * log_dim_fetch + C, 0.5) + D * log_dim_fetch;
        Real gt_min_U = 6.5882 * std::exp(exponent_term);
        Real t_min_hr = (gt_min_U * Ua / Phys::G_STD) / 3600.0;

        return {Hs, Ts, t_min_hr, 0.0};
    }

    /**
     * Calculates Duration-Limited parameters.
     * Uses CEM/SPM formulations adapted with Hurdle & Stive coefficients
     * to ensure asymptotic consistency.
     *
     * @param Ua: Adjusted wind speed [m/s]
     * @param dur_hr: Storm duration [hours]
     * @param depth: Water depth [m]
     * @param is_deep: Flag for deep water regime
     */
    static SMBResult calculate_duration_limited(Real Ua, Real dur_hr, Real depth, bool is_deep) {
        Real dur_sec = dur_hr * 3600.0;
        Real dim_dur = (Phys::G_STD * dur_sec) / Ua;
        Real Hs, Ts, equiv_fetch_m;

        if (is_deep) {
            // Deep Water Duration Limited
            // Adjusted coefficient: 0.283 -> 0.25 (Hurdle & Stive limit)
            // (g*Hs)/Ua^2 = 0.25 * tanh[ 0.000528 * t_hat^0.75 ]
            Real gHs_U2 = 0.25 * std::tanh(0.000528 * std::pow(dim_dur, 0.75));
            Hs = gHs_U2 * (Ua * Ua / Phys::G_STD);

            // Adjusted coefficient: 7.54 -> 8.3 (Hurdle & Stive limit)
            Real gTs_U = 8.3 * std::tanh(0.00379 * std::pow(dim_dur, 0.41));
            Ts = gTs_U * (Ua / Phys::G_STD);

            // Equivalent Fetch: Inverting the Hurdle & Stive Deep Water relationship
            // Approximation for mapping coefficients
            Real time_term = 0.000528 * std::pow(dim_dur, 0.75);
            Real tanh_sq = std::tanh(time_term) * std::tanh(time_term);
            Real dim_fetch_eq;
            
            // Prevent singularity at fully developed state
            if (time_term > 10.0 || tanh_sq >= 1.0) {
                dim_fetch_eq = 200000.0; // Large number representing fully developed
            } else {
                dim_fetch_eq = std::atanh(tanh_sq) / 4.3e-5;
            }
            equiv_fetch_m = dim_fetch_eq * (Ua * Ua / Phys::G_STD);

        } else {
            // Finite Depth Duration Limited (Heuristic Adaptation)
            Real dim_depth = (Phys::G_STD * depth) / (Ua * Ua);

            // Coefficient updates: 0.530 -> 0.6, 0.283 -> 0.25
            Real term_h = 0.000528 * std::pow(dim_dur, 0.75);
            Real tanh_depth_h = std::tanh(0.6 * std::pow(dim_depth, 0.75));
            Hs = (Ua * Ua / Phys::G_STD) * 0.25 * tanh_depth_h * std::tanh(term_h / tanh_depth_h);

            // Coefficient updates: 0.833 -> 0.76, 7.54 -> 8.3
            Real term_t = 0.00379 * std::pow(dim_dur, 0.41);
            Real tanh_depth_t = std::tanh(0.76 * std::pow(dim_depth, 0.375));
            Ts = (Ua / Phys::G_STD) * 8.3 * tanh_depth_t * std::tanh(term_t / tanh_depth_t);

            // Equivalent Fetch Calculation (Simplified Mapping)
            // 0.000528 * t_hat^0.75 ~ 4.3e-5 * F_hat
            Real dim_fetch_eq = (0.000528 / 4.3e-5) * std::pow(dim_dur, 0.75);
            equiv_fetch_m = dim_fetch_eq * (Ua * Ua / Phys::G_STD);
        }

        return {Hs, Ts, 0.0, equiv_fetch_m / 1000.0};
    }

    // --- 3.2 Exact Dispersion Solver (Newton-Raphson) ---
    /**
     * Solves the transcendental Linear Dispersion Relation.
     * Equation: L = (g * T^2 / 2pi) * tanh( 2pi * d / L )
     *
     * The problem is cast in terms of the dimensionless wavenumber kh = k * d.
     * Function to zero: f(kh) = (w^2 * d / g) - kh * tanh(kh) = 0
     *
     * @param T: Wave Period [s]
     * @param d: Water Depth [m]
     * @return kh: The dimensionless wavenumber k*d
     */
    static Real solve_dispersion(Real T, Real d, Real tol=1e-15, int max_iter=100) {
        if (d <= 0) return 0.0;
        Real omega = 2.0 * Phys::PI / T;       // Angular frequency
        Real k0 = (omega * omega) / Phys::G_STD; // Deep water wavenumber (approx)
        Real k0h = k0 * d;

        // Initial Guess using Carvalho (2006) explicit approximation
        // This provides a starting point very close to the true root.
        Real kh = k0h;
        if (k0h > 0) {
             kh = k0h / std::tanh(std::pow(6.0/5.0, k0h) * std::sqrt(k0h));
        } else {
             return 0.0;
        }

        // Newton-Raphson Iteration loop
        for (int i = 0; i < max_iter; ++i) {
            Real t_kh = std::tanh(kh);
            Real f = k0h - kh * t_kh; // Residual
            
            Real sech = 1.0 / std::cosh(kh);
            Real df = -t_kh - kh * (sech * sech); // Derivative
            
            if (std::abs(df) < 1e-15) break; // Avoid division by zero
            
            Real dkh = f / df;
            Real kh_new = kh - dkh;
            
            // Convergence check
            if (std::abs(dkh / kh) < tol) return kh_new;
            kh = kh_new;
        }
        return kh;
    }
};

// ==============================================================================
//  SECTION 4: GUI IMPLEMENTATION & LOGIC
// ==============================================================================

#define IDC_EDIT_U10    101
#define IDC_EDIT_FETCH  102
#define IDC_EDIT_DUR    103
#define IDC_EDIT_DEPTH  104
#define IDC_BTN_CALC    105
#define IDC_OUTPUT      106

HWND hEditU10, hEditFetch, hEditDur, hEditDepth, hOutput, hBtnCalc;
HFONT hUIFont, hMonoFont;

struct AppInputs {
    double U10;
    double fetch_km;
    double duration_hr; // -1 if infinite
    double depth_m;     // -1 if deep water
};

std::string RunSimulation(AppInputs inp) {
    StringLogger log;
    log.print_header("SMB WAVE PREDICTION REPORT");

    // --- 1. Physics Calculations (SMB) ---
    Real Ua = SMBEngine::calculate_adjusted_wind_speed(inp.U10);
    Real fetch_m = inp.fetch_km * 1000.0;
    bool deep_water = (inp.depth_m < 0);
    bool inf_duration = (inp.duration_hr < 0);

    // Calculate Potential A: Fetch Limited
    // (Assumes infinite duration, limited only by distance)
    SMBResult r_fetch;
    if (deep_water) r_fetch = SMBEngine::calculate_deep_water(Ua, fetch_m);
    else r_fetch = SMBEngine::calculate_depth_limited(Ua, fetch_m, inp.depth_m);

    // Calculate Potential B: Duration Limited
    // (Assumes infinite fetch, limited only by time)
    SMBResult r_dur;
    Real dur_val = inf_duration ? 1e9 : inp.duration_hr; 
    r_dur = SMBEngine::calculate_duration_limited(Ua, dur_val, inp.depth_m, deep_water);

    // Controlling Factor Logic
    // The actual sea state is the smaller of the two potentials (Fetch vs Duration).
    Real final_Hs, final_Ts;
    String control_msg;
    Real rel_dur, rel_fetch;
    String dur_label, fetch_label;

    if (r_fetch.Hs <= r_dur.Hs) {
        final_Hs = r_fetch.Hs;
        final_Ts = r_fetch.Ts;
        control_msg = "FETCH-LIMITED";
        rel_dur = r_fetch.t_min;
        dur_label = "Min Duration Req.";
        rel_fetch = inp.fetch_km;
        fetch_label = "Given Fetch";
    } else {
        final_Hs = r_dur.Hs;
        final_Ts = r_dur.Ts;
        control_msg = "DURATION-LIMITED";
        rel_dur = inp.duration_hr;
        dur_label = "Given Duration";
        rel_fetch = r_dur.equiv_fetch;
        fetch_label = "Equivalent Fetch";
    }

    // --- 2. Precise Hydrodynamics (Dispersion) ---
    Real k, L, C, kh, d_actual;

    if (deep_water) {
        // Deep water dispersion (Analytic): L0 = gT^2 / 2pi
        L = (Phys::G_STD * final_Ts * final_Ts) / (2 * Phys::PI);
        k = 2 * Phys::PI / L;
        kh = 100.0; // Effectively infinite
        d_actual = 1e9; // Large number for display logic
    } else {
        d_actual = inp.depth_m;
        // Finite depth dispersion (Numerical): Solve Newton-Raphson
        kh = SMBEngine::solve_dispersion(final_Ts, d_actual);
        k = kh / d_actual;
        L = 2 * Phys::PI / k;
    }
    C = L / final_Ts; // Phase Speed

    // --- 3. Report Generation ---

    // INPUTS
    log.print_subheader("1. INPUT PARAMETERS");
    log.print_row("Wind Speed (U10)", inp.U10, "m/s");
    log.print_row("Adjusted Speed (Ua)", Ua, "m/s");
    log.print_row("Fetch Length", inp.fetch_km, "km");
    
    if (inf_duration) log.print_str("Storm Duration", "Infinite", "hours");
    else log.print_row("Storm Duration", inp.duration_hr, "hours");

    if (deep_water) log.print_str("Water Depth", "Deep Water", "-");
    else log.print_row("Water Depth", inp.depth_m, "m");

    // CONDITIONS
    log.print_subheader("2. LIMITING CONDITIONS");
    log.print_str("Controlling Factor", control_msg);
    log.newline();
    log.print_text("Determining Parameters:");
    if (inf_duration && control_msg == "FETCH-LIMITED") {
        log.print_row(dur_label, rel_dur, "hours");
    } else if (inf_duration) {
         log.print_str(dur_label, "Infinite", "hours");
    } else {
        log.print_row(dur_label, rel_dur, "hours");
    }
    log.print_row(fetch_label, rel_fetch, "km");

    // RESULTS
    log.print_subheader("3. PREDICTION RESULTS");
    log.print_row("Sig. Wave Height (Hs)", final_Hs, "m");
    log.print_row("Sig. Wave Period (Ts)", final_Ts, "s");
    log.print_row("Wave Celerity (C)", C, "m/s");
    log.print_row("Wave Length (L)", L, "m");
    log.print_row("Wave Number (k)", k, "rad/m");

    // INSIGHTS
    log.print_subheader("4. INSIGHTS & OBSERVATIONS");
    
    // Insight A: Growth State
    if (control_msg == "FETCH-LIMITED") {
        log.print_text(">> STATE: Fully Developed for Fetch.");
        log.print_text("   Waves have reached maximum size for this distance.");
    } else {
        log.print_text(">> STATE: Growing Sea State (Duration Limited).");
        log.print_text("   Waves are still growing over time.");
    }
    log.newline();

    // Insight B: Regime & Stability
    String regime;
    Real d_L = (deep_water) ? 1.0 : (d_actual / L);
    if (d_L >= 0.5) regime = "Deep Water (d/L > 0.5)";
    else if (d_L < 0.05) regime = "Shallow Water (d/L < 0.05)";
    else regime = "Transitional/Intermediate";
    
    log.print_str(">> Flow Regime", regime);
    
    // Breaking (Miche Criterion)
    // Limits wave steepness (H/L) in any water depth.
    // H_max / L = 0.142 * tanh(kd)
    Real steepness = final_Hs / L;
    Real limit_steepness = 0.142 * std::tanh(kh);
    
    log.print_row("   Calculated Steepness", steepness, "-");
    log.print_row("   Miche Limit (H/L)", limit_steepness, "-");
    
    if (steepness > limit_steepness) {
        log.print_str("   >> BREAKING STATUS", "UNSTABLE / BREAKING [!]", "");
        log.print_text("      Wave height exceeds the Miche stability limit.");
    } else {
        log.print_str("   >> BREAKING STATUS", "STABLE", "");
        Real safety_margin = (1.0 - steepness/limit_steepness) * 100.0;
        log.print_row("      Stability Margin", safety_margin, "%");
    }

    if (!deep_water) {
        Real ratio = final_Hs / d_actual;
        log.print_row("   Depth Ratio (Hs/d)", ratio, "-");
    }

    // SCENARIOS
    log.print_subheader("5. SCENARIO ANALYSIS");
    
    log.print_text("[A] Fetch-Limited Potential");
    log.print_row("    - Potential Hs", r_fetch.Hs, "m");
    log.print_row("    - Potential Ts", r_fetch.Ts, "s");
    log.print_row("    - Time to Develop", r_fetch.t_min, "hours");
    
    log.newline();
    log.print_text("[B] Duration-Limited Potential");
    if (inf_duration) {
        log.print_text("    - Not applicable (Duration is infinite)");
    } else {
        log.print_row("    - Potential Hs", r_dur.Hs, "m");
        log.print_row("    - Potential Ts", r_dur.Ts, "s");
        log.print_row("    - Equiv. Fetch", r_dur.equiv_fetch, "km");
    }

    return log.get_content();
}

// ==============================================================================
//  SECTION 5: WINDOWS GUI
// ==============================================================================

// Helper: Convert ANSI string to Wide String (for Windows API)
std::wstring to_wstring(const std::string& str) {
    if (str.empty()) return std::wstring();
    int size = MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), NULL, 0);
    std::wstring wstr(size, 0);
    MultiByteToWideChar(CP_UTF8, 0, &str[0], (int)str.size(), &wstr[0], size);
    return wstr;
}

// Helper: Normalize newlines for Windows Edit Control (\r\n)
std::wstring FixNewlines(const std::wstring& in) {
    std::wstring out;
    for (wchar_t c : in) {
        if (c == L'\n') { out += L"\r\n"; } else { out += c; }
    }
    return out;
}

void CreateGUI(HWND hwnd) {
    // Fonts
    hUIFont = CreateFontW(22, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, 
                         OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, DEFAULT_PITCH | FF_SWISS, L"Segoe UI");
    hMonoFont = CreateFontW(20, 0, 0, 0, FW_NORMAL, FALSE, FALSE, FALSE, DEFAULT_CHARSET, 
                           OUT_DEFAULT_PRECIS, CLIP_DEFAULT_PRECIS, DEFAULT_QUALITY, FIXED_PITCH | FF_DONTCARE, L"Consolas");

    // Layout Constants
    int y = 40;
    int x_lbl = 20;
    int x_in = 180;
    int w_in = 130;
    int h_ctl = 32;
    int step = 55;

    auto AddLabel = [&](const wchar_t* txt, int yy) {
        SendMessageW(CreateWindowW(L"STATIC", txt, WS_CHILD|WS_VISIBLE, x_lbl, yy, 155, h_ctl, hwnd, NULL, NULL, NULL), 
                    WM_SETFONT, (WPARAM)hUIFont, TRUE);
    };
    auto AddEdit = [&](int id, const wchar_t* def, int yy) {
        HWND h = CreateWindowW(L"EDIT", def, WS_CHILD|WS_VISIBLE|WS_BORDER|ES_AUTOHSCROLL, x_in, yy, w_in, h_ctl, hwnd, (HMENU)(INT_PTR)id, NULL, NULL);
        SendMessageW(h, WM_SETFONT, (WPARAM)hUIFont, TRUE);
        return h;
    };
    auto AddNote = [&](const wchar_t* txt, int yy) {
        SendMessageW(CreateWindowW(L"STATIC", txt, WS_CHILD|WS_VISIBLE, x_in + w_in + 10, yy, 120, h_ctl, hwnd, NULL, NULL, NULL), 
                    WM_SETFONT, (WPARAM)hUIFont, TRUE);
    };

    // Labels
    AddLabel(L"Wind Speed (U10):", y); hEditU10 = AddEdit(IDC_EDIT_U10, L"30.0", y); AddNote(L"m/s", y); y += step;
    AddLabel(L"Fetch Length:", y);     hEditFetch = AddEdit(IDC_EDIT_FETCH, L"10.0", y); AddNote(L"km", y); y += step;
    AddLabel(L"Duration:", y);         hEditDur = AddEdit(IDC_EDIT_DUR, L"", y); AddNote(L"hrs", y); y += step;
    AddLabel(L"Water Depth:", y);      hEditDepth = AddEdit(IDC_EDIT_DEPTH, L"", y); AddNote(L"m", y); y += step;

    // Calculate Button
    y += 20;
    hBtnCalc = CreateWindowW(L"BUTTON", L"CALCULATE", WS_TABSTOP|WS_VISIBLE|WS_CHILD|BS_DEFPUSHBUTTON, 
                             x_lbl, y, 320, 50, hwnd, (HMENU)IDC_BTN_CALC, NULL, NULL);
    SendMessageW(hBtnCalc, WM_SETFONT, (WPARAM)hUIFont, TRUE);

    // Output Box
    // Height reduced to 460 to match main window size of 540
    hOutput = CreateWindowW(L"EDIT", L"", WS_CHILD|WS_VISIBLE|WS_BORDER|ES_MULTILINE|ES_AUTOVSCROLL|WS_VSCROLL|ES_READONLY|WS_HSCROLL|ES_AUTOHSCROLL, 
                            380, 20, 600, 640, hwnd, (HMENU)IDC_OUTPUT, NULL, NULL);
    SendMessageW(hOutput, WM_SETFONT, (WPARAM)hMonoFont, TRUE);
}

struct ThreadData { AppInputs inp; HWND hBtn; HWND hOut; };

DWORD WINAPI Worker(LPVOID lpParam) {
    ThreadData* p = (ThreadData*)lpParam;
    std::string report = RunSimulation(p->inp);
    
    std::ofstream f("report.txt");
    if(f.is_open()) { f << report; f.close(); }

    std::wstring wtxt = FixNewlines(to_wstring(report));
    SetWindowTextW(p->hOut, wtxt.c_str());
    
    EnableWindow(p->hBtn, TRUE);
    SetWindowTextW(p->hBtn, L"CALCULATE");
    
    delete p;
    return 0;
}

LRESULT CALLBACK WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp) {
    if (msg == WM_CREATE) { CreateGUI(hwnd); SetWindowPos(hwnd, 0, 0, 0, 1024, 700, SWP_NOMOVE|SWP_NOZORDER); return 0; }
    if (msg == WM_DESTROY) { PostQuitMessage(0); return 0; }
    
    if (msg == WM_COMMAND && LOWORD(wp) == IDC_BTN_CALC) {
        auto GetVal = [&](HWND h, bool& empty) {
            wchar_t b[64]; GetWindowTextW(h, b, 63);
            if(b[0] == 0) { empty=true; return 0.0; }
            empty=false; return wcstod(b, 0);
        };

        AppInputs inp; bool e;
        inp.U10 = GetVal(hEditU10, e);
        if(e || inp.U10 <= 0) { MessageBoxW(hwnd, L"Wind Speed must be > 0", L"Error", MB_ICONERROR); return 0; }
        
        inp.fetch_km = GetVal(hEditFetch, e);
        if(e || inp.fetch_km <= 0) { MessageBoxW(hwnd, L"Fetch must be > 0", L"Error", MB_ICONERROR); return 0; }
        
        inp.duration_hr = GetVal(hEditDur, e);
        if(e) inp.duration_hr = -1;
        
        inp.depth_m = GetVal(hEditDepth, e);
        if(e) inp.depth_m = -1;

        EnableWindow(hBtnCalc, FALSE);
        SetWindowTextW(hBtnCalc, L"Computing...");

        ThreadData* td = new ThreadData{inp, hBtnCalc, hOutput};
        CloseHandle(CreateThread(0,0,Worker,td,0,0));
        return 0;
    }
    return DefWindowProcW(hwnd, msg, wp, lp);
}

int WINAPI WinMain(HINSTANCE h, HINSTANCE, LPSTR, int n) {
    WNDCLASSEXW wc = {sizeof(WNDCLASSEXW), 0, WndProc, 0, 0, h, LoadIcon(0, IDI_APPLICATION), 
                      LoadCursor(0, IDC_ARROW), (HBRUSH)(COLOR_WINDOW+1), 0, L"SMBClass", 0};
    RegisterClassExW(&wc);
    
    // Height set to 700 to accommodate the taller text output window
    HWND w = CreateWindowExW(0, L"SMBClass", L"SMB Wave Calculator", 
                             WS_OVERLAPPED | WS_CAPTION | WS_SYSMENU | WS_MINIMIZEBOX, 
                             CW_USEDEFAULT, CW_USEDEFAULT, 1024, 700, 0, 0, h, 0);
    ShowWindow(w, n);
    MSG m; while(GetMessageW(&m,0,0,0)) { TranslateMessage(&m); DispatchMessageW(&m); }
    return m.wParam;
}