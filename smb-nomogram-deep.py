#!/usr/bin/env python3
"""
smb-nomogram-deep.py

This script generates a set of three nomograms for the Sverdrup-Munk-Bretschneider
(SMB) method for deep water wave prediction. It creates a single, multi-page PDF
file named 'smb-nomogram-deep.pdf' with separate pages for:
1. Significant Wave Height (Hs)
2. Significant Wave Period (Ts)
3. Minimum required wind duration (t_min)

Each nomogram includes an isopleth for the Lake Garda example
(Wind Speed = 25 m/s, Fetch = 45 km).

The calculations are based on the deep-water formulas documented in the
U.S. Army Corps of Engineers' Coastal Engineering Manual and Shore Protection Manual.
"""

import sys
import math
import numpy as np
import scipy
import io
import re
import copy
import os

# --- FIX FOR SCIPY ATTRIBUTE ERROR ---
# Older libraries like pynomo expect scipy.arange, which was removed in new versions.
# This line restores it by pointing it to numpy.arange.
scipy.arange = np.arange

# --- Optional parameter to control isopleth plotting ---
# Set this to True to plot the isopleths (example lines) or False to hide them.
PLOT_ISOPLETHS = False
# -------------------------------------------------------------

# --- Ensure nomogen and pynomo are in the path ---
# If nomogen.py is in the parent directory:
sys.path.insert(0, "..")

from nomogen import Nomogen
from pynomo.nomographer import Nomographer
from pyx import color, text
from pypdf import PdfWriter, PdfReader # Used for merging PDFs

# This line helps ensure PyX uses the LaTeX engine for high-quality text rendering.
text.set(text.LatexEngine)

# =============================================================================
# SMB WAVE PREDICTION MODEL (DEEP WATER)
# =============================================================================
def calculate_deep_water(wind_speed, fetch, gravity=9.81):
    """
    Calculates fetch-limited wave properties in deep water using the SMB method.
    """
    if wind_speed <= 0 or fetch <= 0:
        return 0, 0, 0
        
    dim_fetch = (gravity * fetch) / (wind_speed**2)
    gHs_U2 = 0.283 * math.tanh(0.0125 * (dim_fetch**0.42))
    Hs = gHs_U2 * (wind_speed**2 / gravity)
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

# =============================================================================
# WRAPPER FUNCTIONS FOR NOMOGEN
# =============================================================================
def get_Hs(wind_speed, fetch_km):
    """Returns only the Significant Wave Height."""
    fetch_m = fetch_km * 1000
    Hs, _, _ = calculate_deep_water(wind_speed, fetch_m)
    return Hs

def get_Ts(wind_speed, fetch_km):
    """Returns only the Significant Wave Period."""
    fetch_m = fetch_km * 1000
    _, Ts, _ = calculate_deep_water(wind_speed, fetch_m)
    return Ts

def get_Duration(wind_speed, fetch_km):
    """Returns only the Minimum Duration."""
    fetch_m = fetch_km * 1000
    _, _, t_min = calculate_deep_water(wind_speed, fetch_m)
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

    # --- Lake Garda Isopleth Example ---
    isopleth_wind = 25.0
    isopleth_fetch = 45.0
    hs_iso, ts_iso, dur_iso = calculate_deep_water(isopleth_wind, isopleth_fetch * 1000)
    print("--- Lake Garda Example (U=25 m/s, F=45 km) ---")
    print(f"Calculated Hs: {hs_iso:.2f} m, Ts: {ts_iso:.2f} s, Duration: {dur_iso:.2f} hours")
    print("------------------------------------------------")

    # --- Define parameters for the nomograms to be created ---
    nomograms_to_build = [
        {"nomo_func": get_Hs, "title": r"Significant Wave Height, $H_s$ (m)", "scale_type": "linear smart", "temp_file": "temp-hs-deep.pdf"},
        {"nomo_func": get_Ts, "title": r"Significant Wave Period, $T_s$ (s)", "scale_type": "linear smart", "temp_file": "temp-ts-deep.pdf"},
        {"nomo_func": get_Duration, "title": r"Minimum Duration, $t_{min}$ (hours)", "scale_type": "log", "temp_file": "temp-duration-deep.pdf"},
    ]
    
    temp_files = []

    # --- Generate each nomogram as a separate, temporary PDF ---
    for nomo in nomograms_to_build:
        print(f"\nGenerating temporary nomogram for: {nomo['title']}")
        temp_files.append(nomo["temp_file"])
        
        # --- Find true min/max of the output function by searching a grid ---
        grid_points = 20
        wind_range = np.linspace(wind_min, wind_max, grid_points)
        fetch_range = np.linspace(fetch_min, fetch_max, grid_points)
        output_values = [nomo["nomo_func"](w, f) for w in wind_range for f in fetch_range]
        
        out_min = min(output_values)
        out_max = max(output_values)
        
        # Ensure the isopleth value is within the scale range, just in case
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
            'u_min': wind_min, 'u_max': wind_max, 'title': r'Wind Speed, U (m/s)',
            'scale_type': 'linear smart', 'tick_levels': 3, 'tick_text_levels': 1,
            'text_format': r'\Large{%g}',
        }
        
        fetch_axis_params = {
            'u_min': fetch_min, 'u_max': fetch_max, 'title': r'Fetch, F (km)',
            'scale_type': 'log', 'tick_levels': 3, 'tick_text_levels': 2,
            'text_format': r'\Large{%g}',
        }
        
        # --- Parameters for the calculated output axis ---
        output_axis_params = {
            'u_min': out_min, 'u_max': out_max, 'title': nomo['title'],
            'scale_type': nomo['scale_type'], 'tick_levels': 3,
            'tick_text_levels': 2,
            # Use LaTeX \textbf for bold numbers on the scale
            'text_format': r'\Large\textbf{%g}',
            # Use line_params to make the axis line thicker (bold)
            'line_params': {'linewidth': 'thick'},
        }
        
        block_params = {
            'block_type': 'type_9', 'f1_params': wind_axis_params,
            'f2_params': output_axis_params, 'f3_params': fetch_axis_params,
        }
        
        if PLOT_ISOPLETHS:
            block_params['isopleth_values'] = [[isopleth_wind, 'x', isopleth_fetch]]
        
        main_params = {
            'filename': nomo["temp_file"], 'paper_height': 29.7, 'paper_width': 21,
            'block_params': [block_params], 'transformations': [('scale paper',)],
            'title_str': f"\\Large SMB Deep Water Nomogram: {nomo['title']}",
            'isopleth_params': [{'color': 'Blue', 'linewidth': 'thick', 'circle_size': 0.1}],
            'npoints': 11,
        }

        Nomogen(nomo['nomo_func'], main_params)
        Nomographer(main_params)

    # --- Merge the temporary PDFs into a single file ---
    merger = PdfWriter()
    final_filename = "smb-nomogram-deep.pdf"
    
    print(f"\nMerging temporary files into {final_filename}...")
    # First, append all the generated PDF pages to the writer object
    for pdf in temp_files:
        merger.append(pdf)
    
    # To solve the issue of the last page being a different size, we will
    # iterate through the pages in the merged document and enforce the same
    # page size (MediaBox) for all pages, using the first page as the reference.
    if len(merger.pages) > 0:
        print("Standardizing page sizes to ensure uniformity...")
        # Get the MediaBox of the first page, which is assumed to be the correct A4 size.
        reference_mediabox = merger.pages[0].mediabox
        
        # Loop through all pages and set their MediaBox to the reference size.
        # This ensures every page in the final PDF has identical dimensions.
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
