# SMB Wave Prediction Model

A collection of Python scripts for predicting wind-generated wave characteristics using the Sverdrup-Munk-Bretschneider (SMB) method. These tools calculate significant wave height (Hs), significant wave period (Ts), and the minimum required storm duration based on meteorological inputs.

## Overview

This repository provides a practical implementation of the foundational SMB wave prediction model, a cornerstone of coastal and ocean engineering for decades. It is designed for engineers, scientists, and students who need to make preliminary estimates of wave conditions for design, planning, or research purposes.

The core model operates in two distinct modes, based on the primary limiting condition for wave growth:

1. **Fetch-Limited (Finite Depth or Deep Water):** For scenarios where wave growth is primarily limited by the available fetch. This mode can handle both deep water conditions (where the seabed does not influence wave generation) and finite depth conditions (where water depth is a significant factor). It calculates Hs, Ts, and the minimum duration required for a fully developed sea state over the given fetch.

2. **Duration-Limited (Finite Depth or Deep Water):** Calculates wave parameters when the wind event is too short for waves to become fully developed over the available fetch. This mode considers wind speed, storm duration, and water depth, and also calculates the equivalent fetch that would produce the same wave conditions.

The `calculator.py` script comprehensively evaluates both fetch and duration limits for given inputs to determine the actual controlling wave growth factor.

## Methodology

The scripts are based on the semi-empirical Sverdrup-Munk-Bretschneider (SMB) method, which relates dimensionless wave parameters to wind conditions. The core principle is that wave growth is limited by either fetch (spatial constraint) or duration (temporal constraint).

### Deep Water Formulas (Fetch-Limited)

For deep water, the script uses the revised Bretschneider (1970) equations:

* **Dimensionless Fetch:**

  $$\hat{F} = \frac{gF}{U^{2}}$$

* **Significant Wave Height (**$H\_s$**):**

  $$\frac{gH_{s}}{U^{2}} = 0.283 \tanh\left[0.0125\left(\hat{F}\right)^{0.42}\right]$$

* **Significant Wave Period (**$T\_s$**):**

  $$\frac{gT_{s}}{U} = 7.54 \tanh\left[0.077\left(\hat{F}\right)^{0.25}\right]$$

* **Minimum Duration (**$t\_{min}$**):**

  $$t_{min} = \frac{U}{g} \cdot 6.5882 \cdot \exp\left(\sqrt{0.0161 \cdot (\ln \hat{F})^2 - 0.3692 \cdot \ln \hat{F} + 2.2024} + 0.8798 \cdot \ln \hat{F}\right)$$

### Depth-Limited Formulas (Fetch-Limited)

For shallower water, the script uses formulas from the Shore Protection Manual that incorporate a dimensionless depth parameter ($\\hat{d} = gd/U^2$):

* **Significant Wave Height (**$H\_s$**):**

  $$H_s = \frac{U^2}{g} \cdot 0.283 \cdot \tanh(0.530 \hat{d}^{0.75}) \cdot \tanh\left[\frac{0.00565 \hat{F}^{0.5}}{\tanh(0.530 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period (**$T\_s$**):**

  $$T_s = \frac{U}{g} \cdot 7.54 \cdot \tanh(0.833 \hat{d}^{0.375}) \cdot \tanh\left[\frac{0.0379 \hat{F}^{0.333}}{\tanh(0.833 \hat{d}^{0.375})}\right]$$

* **Minimum Duration (**$t\_{min}$**):**

  $$t_{min} = \frac{U}{g} \cdot 6.5882 \cdot \exp\left(\sqrt{0.0161 \cdot (\ln \hat{F})^2 - 0.3692 \cdot \ln \hat{F} + 2.2024} + 0.8798 \cdot \ln \hat{F}\right)$$

### Duration-Limited Formulas (Finite Depth or Deep Water)

When the wind event duration ($t$) is the limiting factor, the following formulas are used.
Dimensionless duration: $\\hat{t} = \\frac{gt}{U}$.

* **Significant Wave Height (**$H\_s$**) - Deep Water:**

  $$\frac{gH_{s}}{U^{2}} = 0.283 \tanh\left[0.000528 \left(\hat{t}\right)^{0.75}\right]$$

* **Significant Wave Period (**$T\_s$**) - Deep Water:**

  $$\frac{gT_{s}}{U} = 7.54 \tanh\left[0.00379 \left(\hat{t}\right)^{0.41}\right]$$

* **Significant Wave Height (**$H\_s$**) - Finite Depth (Heuristic):**

  $$H_s = \frac{U^2}{g} \cdot 0.283 \cdot \tanh(0.530 \hat{d}^{0.75}) \cdot \tanh\left[\frac{0.000528 \hat{t}^{0.75}}{\tanh(0.530 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period (**$T\_s$**) - Finite Depth (Heuristic):**

  $$T_s = \frac{U}{g} \cdot 7.54 \cdot \tanh(0.833 \hat{d}^{0.375}) \cdot \tanh\left[\frac{0.00379 \hat{t}^{0.41}}{\tanh(0.833 \hat{d}^{0.375})}\right]$$

## Features

* **Dual-Mode Calculation:** Accurately applies formulas for fetch-limited (deep or finite depth) and duration-limited (deep or finite depth) conditions, and determines the controlling factor.

* **Comprehensive Outputs:** Calculates Significant Wave Height (Hs), Significant Wave Period (Ts), Minimum Storm Duration (t_min) where applicable, and Equivalent Fetch for duration-limited cases.

* **Interactive Interface:** A simple command-line interface guides the user through the input process, allowing for flexible input of fetch, duration, and optional depth.

* **Validated Formulas:** The implemented equations are based on authoritative sources, including the U.S. Army's Coastal Engineering Manual and Shore Protection Manual.

* **Data Visualization:** Generates contour charts for both deep and depth-limited conditions.

* **Nomogram Generation:** Creates printable nomograms for quick graphical estimations.

* **Tabular Output:** Generates a PDF table summarizing wave calculations for various parameters.

## Scripts Description

This repository contains the following Python scripts:

### `calculator.py`

This is the main interactive script for performing individual SMB wave calculations.

* **Functionality:**

  * Prompts the user for wind speed, fetch length, storm duration, and optionally water depth.

  * Calculates wave parameters for both fetch-limited and duration-limited scenarios based on the provided inputs.

  * Determines the **controlling wave growth factor** (whether fetch or duration is the primary limit) and displays the corresponding significant wave height (Hs), significant wave period (Ts), and relevant duration/fetch values.

  * Outputs all calculations to a `report.txt` file, mirroring the command line output.

* **Usage:** Run directly from the command line and follow the prompts.

  ```
  python calculator.py
  ```

### `chart.py`

Generates a combined contour chart for SMB wave parameters in **deep water** conditions.

* **Functionality:**

  * Uses `matplotlib` to create a single plot showing contours of Hs, Ts, and t_min.

  * Displays wave parameters as functions of wind speed and fetch length for deep water.

  * Utilizes different black line styles (solid for Hs, dashed for Ts, dotted for t_min) for clarity.

  * The chart is generated in A3 landscape format for better readability.

* **Usage:** Run directly to generate and display the chart.

  ```
  python chart.py
  ```

### `chart_10m.py`

Generates a combined contour chart for SMB wave parameters in **depth-limited water** conditions, specifically for a fixed water depth of 10 meters.

* **Functionality:**

  * Similar to `chart.py`, but tailored for depth-limited scenarios.

  * Calculates and plots contours of Hs, Ts, and t_min for a constant water depth (defaulting to 10m).

  * Provides a visual representation of how wave parameters change with wind speed and fetch in shallow water.

  * The chart is generated in A3 landscape format.

* **Usage:** Run directly to generate and display the chart. The `FIXED_DEPTH` variable can be modified within the script.

  ```
  python chart_10m.py
  ```

### `smb-nomogram-deep.py`

Generates a multi-page PDF containing three nomograms for **deep water** wave prediction.

* **Functionality:**

  * Creates separate nomograms for Significant Wave Height (Hs), Significant Wave Period (Ts), and Minimum required wind duration (t_min).

  * Outputs a single PDF file named `smb-nomogram-deep.pdf`.

  * Requires `pynomo` and `nomogen` libraries for nomogram generation.

* **Usage:** Run directly to generate the PDF.

  ```
  python smb-nomogram-deep.py
  ```

### `smb-nomogram-shallow.py`

Generates a multi-page PDF containing three nomograms for **depth-limited (shallow water)** wave prediction, configured for a fixed water depth of 10 meters.

* **Functionality:**

  * Similar to `smb-nomogram-deep.py`, but specifically for shallow water conditions.

  * Generates nomograms for Hs, Ts, and t_min at a constant water depth (defaulting to 10m).

  * Outputs a single PDF file named `smb-nomogram-shallow.pdf`.

  * Requires `pynomo` and `nomogen` libraries.

* **Usage:** Run directly to generate the PDF. The `WATER_DEPTH` variable can be modified within the script.

  ```
  python smb-nomogram-shallow.py
  ```

### `tables.py`

Generates a comprehensive PDF table summarizing SMB wave calculations for various combinations of wind speed, fetch, and depth.

* **Functionality:**

  * Calculates Hs, Ts, and t_min for predefined ranges of wind speeds (5-35 m/s), fetches (0-50 km), and depths (Deep Water, 100m, 50m, 25m, 10m, 5m, 1m).

  * Organizes the results into a well-formatted table within a PDF document.

  * Uses `reportlab` for PDF generation, ensuring a professional and readable output.

* **Usage:** Run directly to generate the PDF file `comprehensive_wave_calculations.pdf`.

  ```
  python tables.py
  ```

## How to Use

1. **Prerequisites:** Ensure you have Python 3 installed. You will also need to install the following libraries:

   * `numpy`

   * `scipy` (for `smb-nomogram-deep.py` and `smb-nomogram-shallow.py`, specifically for `scipy.arange` compatibility fix)

   * `matplotlib` (for `chart.py` and `chart_10m.py`)

   * `reportlab` (for `tables.py`)

   * `pynomo` (for `smb-nomogram-deep.py` and `smb-nomogram-shallow.py`)

   * `PyPDF2` (for `smb-nomogram-deep.py` and `smb-nomogram-shallow.py` to merge PDFs)

   * `PyX` (for `smb-nomogram-deep.py` and `smb-nomogram-shallow.py` for LaTeX rendering)

   You can install most of them using pip (`nomogen.py` is expected to be in same directory as script):

  ```
  pip install numpy matplotlib reportlab PyPDF2 PyX pynomo
  ```

2. **Running Scripts:**

   * For interactive calculations, run `python calculator.py`.

   * For deep water charts, run `python chart.py`.

   * For 10m depth-limited charts, run `python chart_10m.py`.

   * For deep water nomograms, run `python smb-nomogram-deep.py`.

   * For shallow water nomograms, run `python smb-nomogram-shallow.py`.

   * For the comprehensive table, run `python tables.py`.

## Assumptions and Limitations

* **Steady-State Wind:** The model assumes that the wind speed and direction are uniform and constant across the entire fetch for the specified duration. This is an idealization not always met in nature.

* **Input Data Quality:** The accuracy of the results is highly dependent on the quality of the inputs. For best results:

  * **Wind Speed:** Should be the standard 10-meter overwater wind speed, adjusted for atmospheric stability if possible.

  * **Fetch Length:** Should be the "effective fetch," which accounts for the geometry of the water body, not just a straight-line distance.

* **Heuristic for Duration-Limited (Finite Depth):** The formulas used for duration-limited conditions in finite depth (Option 2 in `calculator.py`) are an adaptation based on the structure of SMB equations for fetch-limited finite depth. Direct empirical formulas for this specific combined scenario are less common in basic SMB literature. While providing a reasonable estimate, these should be used with awareness of their heuristic nature.

## Bibliography

* U.S. Army. (2008). *Coastal Engineering Manual* (EM 1110-2-1100). Washington, DC: U.S. Army Corps of Engineers.

* U.S. Army. (1984). *Shore Protection Manual*. Vicksburg, MS: U.S. Army Engineer Waterways Experiment Station.

* Bretschneider, C. L. (1970). *Wave forecasting relations for wave generation*. Look Lab, Hawaii, 1(3).

* Etemad-Shahidi, A., Kazeminezhad, M. H., & Mousavi, S. J. (2009). On the prediction of wave parameters using simplified methods. *Journal of Coastal Research, SI 56*, 505-509.
