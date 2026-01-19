#!/usr/bin/env python3
"""
smb-nomogram-shallow.py

This script generates a set of three nomograms for the Sverdrup-Munk-Bretschneider
(SMB) method for depth-limited (shallow water) wave prediction. It is configured
for a fixed water depth of 10 meters.

It creates a single, multi-page PDF file named 'smb-nomogram-shallow.pdf'
with separate pages for:
1. Significant Wave Height (Hs)
2. Significant Wave Period (Ts)
3. Minimum required wind duration (t_min)

Each nomogram includes an isopleth for the Lake Garda example
(Wind Speed = 25 m/s, Fetch = 45 km) calculated for a 10m depth.

UPDATED FORMULATION (Hurdle & Stive, 1989):
The calculations for Hs and Ts have been updated to use the revised unified
equations from Hurdle & Stive (1989). These formulas correct the inconsistencies
found in the original SPM (1984) shallow water equations, providing a smooth
transition between depth-limited and deep water conditions.
"""

import sys
import math
import numpy as np
import scipy
import os

# --- FIX FOR SCIPY ATTRIBUTE ERROR ---
# Older libraries like pynomo expect scipy.arange, which was removed in new versions.
# This line restores it by pointing it to numpy.arange.
scipy.arange = np.arange

# --- Optional parameter to control isopleth plotting ---
# Set this to True to plot the isopleths (example lines) or False to hide them.
PLOT_ISOPLETHS = True
# -------------------------------------------------------------

# --- Configuration for Shallow Water ---
WATER_DEPTH = 10.0  # Fixed water depth in meters
# -------------------------------------------------------------


# --- Ensure nomogen and pynomo are in the path ---
# This assumes nomogen.py is in the same or parent directory.
try:
    from nomogen import Nomogen
    from pynomo.nomographer import Nomographer
except ImportError:
    # If nomogen.py is in the parent directory:
    sys.path.insert(0, "..")
    from nomogen import Nomogen
    from pynomo.nomographer import Nomographer

from pypdf import PdfWriter
from pyx import text

# This line helps ensure PyX uses the LaTeX engine for high-quality text rendering.
text.set(text.LatexEngine)


# =============================================================================
# SMB WAVE PREDICTION MODEL (DEPTH-LIMITED / SHALLOW WATER)
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

def calculate_shallow_water(wind_speed_adjusted, fetch, depth, gravity=9.8066):
    """
    Calculates wave properties for depth-limited conditions using the Revised SMB method.
    """
    if wind_speed_adjusted <= 0 or fetch <= 0 or depth <= 0:
        return 0, 0, 0

    # --- Dimensionless Parameters ---
    dim_fetch = (gravity * fetch) / wind_speed_adjusted**2
    dim_depth = (gravity * depth) / wind_speed_adjusted**2

    # --- Revised Significant Wave Height (Hs) ---
    # Hurdle & Stive (1989), Eq 4.1
    depth_term_h = math.tanh(0.6 * dim_depth**0.75)
    # Protection against division by zero if depth is extremely small
    if depth_term_h == 0: depth_term_h = 1e-6
    
    fetch_term_inner_h = (4.3e-5 * dim_fetch) / (depth_term_h**2)
    gHs_U2 = 0.25 * depth_term_h * (math.tanh(fetch_term_inner_h))**0.5
    Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

    # --- Revised Significant Wave Period (Ts) ---
    # Hurdle & Stive (1989), Eq 4.2
    depth_term_t = math.tanh(0.76 * dim_depth**0.375)
    if depth_term_t == 0: depth_term_t = 1e-6

    fetch_term_inner_t = (4.1e-5 * dim_fetch) / (depth_term_t**3)
    gTs_U = 8.3 * depth_term_t * (math.tanh(fetch_term_inner_t))**(1/3)
    Ts = gTs_U * (wind_speed_adjusted / gravity)

    # --- Minimum Wind Duration (t_min) Calculation ---
    # REVISED (Hurdle & Stive, 1989, Eq 4.12 & Section 4.4)
    # The paper recommends using the deep water power law for limiting duration
    # in all circumstances to ensure consistent "Effective Fetch" inversion.
    dim_duration = 65.9 * (dim_fetch**(2.0/3.0))
    
    t_min_seconds = dim_duration * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600

    return Hs, Ts, t_min_hours

# =============================================================================
# WRAPPER FUNCTIONS FOR NOMOGEN
# =============================================================================
def get_Hs(U10, fetch_km):
    """Returns only the Significant Wave Height for shallow water."""
    Ua = calculate_adjusted_wind_speed(U10)
    fetch_m = fetch_km * 1000
    Hs, _, _ = calculate_shallow_water(Ua, fetch_m, WATER_DEPTH)
    return Hs

def get_Ts(U10, fetch_km):
    """Returns only the Significant Wave Period for shallow water."""
    Ua = calculate_adjusted_wind_speed(U10)
    fetch_m = fetch_km * 1000
    _, Ts, _ = calculate_shallow_water(Ua, fetch_m, WATER_DEPTH)
    return Ts

def get_Duration(U10, fetch_km):
    """Returns only the Minimum Duration for shallow water."""
    Ua = calculate_adjusted_wind_speed(U10)
    fetch_m = fetch_km * 1000
    _, _, t_min = calculate_shallow_water(Ua, fetch_m, WATER_DEPTH)
    return t_min

# =============================================================================
# NOMOGRAM GENERATION
# =============================================================================
def generate_nomograms():
    """
    Main function to define and generate the three nomograms into a single PDF.
    """
    # --- Define Input Ranges ---
    wind_min, wind_max = 5.0, 35.0
    fetch_min, fetch_max = 1, 50.0

    # --- Lake Garda Isopleth Example (re-calculated for d=10m) ---
    isopleth_U10 = 25.0
    isopleth_fetch = 45.0
    isopleth_Ua = calculate_adjusted_wind_speed(isopleth_U10) # Calculate Ua for isopleth
    hs_iso, ts_iso, dur_iso = calculate_shallow_water(isopleth_Ua, isopleth_fetch * 1000, WATER_DEPTH)
    print(f"--- Lake Garda Example (U10={isopleth_U10} m/s, F={isopleth_fetch} km, d={WATER_DEPTH}m) ---")
    print(f"Calculated Hs: {hs_iso:.2f} m, Ts: {ts_iso:.2f} s, Duration: {dur_iso:.2f} hours")
    print("------------------------------------------------")

    # --- Define parameters for the nomograms to be created ---
    nomograms_to_build = [
        {"nomo_func": get_Hs, "title": r"Significant Wave Height, $H_s$ (m)", "scale_type": "linear smart", "temp_file": "temp-hs-shallow.pdf"},
        {"nomo_func": get_Ts, "title": r"Significant Wave Period, $T_s$ (s)", "scale_type": "linear smart", "temp_file": "temp-ts-shallow.pdf"},
        {"nomo_func": get_Duration, "title": r"Minimum Duration, $t_{min}$ (hours)", "scale_type": "log", "temp_file": "temp-duration-shallow.pdf"},
    ]
    
    temp_files = []

    # --- Generate each nomogram as a separate, temporary PDF ---
    for nomo in nomograms_to_build:
        print(f"\nGenerating temporary nomogram for: {nomo['title']}")
        temp_files.append(nomo["temp_file"])
    
        # --- Find true min/max of the output function by searching a grid ---
        grid_points = 200 
        wind_range = np.linspace(wind_min, wind_max, grid_points)
        fetch_range = np.linspace(fetch_min, fetch_max, grid_points)
        
        output_values = [nomo["nomo_func"](w, f) for w in wind_range for f in fetch_range]
        
        out_min = min(output_values)
        out_max = max(output_values)

        # Ensure the isopleth value is within the scale range
        if "Height" in nomo["title"]:
            out_min = min(out_min, hs_iso)
            out_max = max(out_max, hs_iso)
        elif "Period" in nomo["title"]:
            out_min = min(out_min, ts_iso)
            out_max = max(out_max, ts_iso)
        elif "Duration" in nomo["title"]:
            out_min = min(out_min, dur_iso)
            out_max = max(out_max, dur_iso)

        wind_axis_params = {
            'u_min': wind_min, 'u_max': wind_max, 'title': r'Wind Speed, $U_{10}$ (m/s)', # Changed label to U10
            'scale_type': 'linear smart', 'tick_levels': 3, 'tick_text_levels': 1,
            'text_format': r'\Large{%g}',
        }
        
        fetch_axis_params = {
            'u_min': fetch_min, 'u_max': fetch_max, 'title': r'Fetch, F (km)',
            'scale_type': 'log', 'tick_levels': 3, 'tick_text_levels': 2,
            'text_format': r'\Large{%g}',
        }
        
        output_axis_params = {
            'u_min': out_min, 'u_max': out_max, 'title': nomo['title'],
            'scale_type': nomo['scale_type'], 'tick_levels': 3,
            'tick_text_levels': 2,
            'text_format': r'\Large\textbf{%g}',
            'line_params': {'linewidth': 'thick'},
        }
        
        block_params = {
            'block_type': 'type_9', 'f1_params': wind_axis_params,
            'f2_params': output_axis_params, 'f3_params': fetch_axis_params,
        }
        
        if PLOT_ISOPLETHS:
            block_params['isopleth_values'] = [[isopleth_U10, 'x', isopleth_fetch]] # Use U10 for isopleth
        
        main_params = {
            'filename': nomo["temp_file"], 'paper_height': 29.7, 'paper_width': 21,
            'block_params': [block_params], 'transformations': [('scale paper',)],
            'title_str': f"\\Large SMB Shallow Water Nomogram (d={WATER_DEPTH}m): {nomo['title']}",
            'isopleth_params': [{'color': 'Blue', 'linewidth': 'thick', 'circle_size': 0.1}],
            'npoints': 16,
        }

        Nomogen(nomo['nomo_func'], main_params)
        Nomographer(main_params)

    # --- Merge the temporary PDFs into a single file ---
    merger = PdfWriter()
    final_filename = "smb-nomogram-shallow.pdf"
    
    print(f"\nMerging temporary files into {final_filename}...")
    for pdf in temp_files:
        merger.append(pdf)
    
    # Standardize page sizes to ensure uniformity
    if len(merger.pages) > 0:
        print("Standardizing page sizes...")
        reference_mediabox = merger.pages[0].mediabox
        for page in merger.pages:
            page.mediabox = reference_mediabox

    merger.write(final_filename)
    merger.close()
    print(f"Successfully created {final_filename}")

    # --- Clean up temporary files ---
    print("Cleaning up temporary files...")
    for pdf in temp_files:
        try:
            os.remove(pdf)
        except OSError as e:
            print(f"Error removing file {pdf}: {e}")

    print("Done.")

if __name__ == "__main__":
    generate_nomograms()