# SMB Wave Prediction Model

A collection of Python scripts for predicting wind-generated wave characteristics using the Sverdrup-Munk-Bretschneider (SMB) method. These tools calculate significant wave height ($H\_s$), significant wave period ($T\_s$), and the minimum required storm duration based on meteorological inputs.

## Overview

This repository provides a practical implementation of the foundational SMB wave prediction model, a cornerstone of coastal and ocean engineering for decades. It is designed for engineers, scientists, and students who need to make preliminary estimates of wave conditions for design, planning, or research purposes.

The core model operates in two distinct modes, based on the primary limiting condition for wave growth:

1. **Fetch-Limited (Finite Depth or Deep Water):** For scenarios where wave growth is primarily limited by the available fetch. This mode can handle both deep water conditions (where the seabed does not influence wave generation) and finite depth conditions (where water depth is a significant factor). It calculates $H\_s$, $T\_s$, and the minimum duration required for a fully developed sea state over the given fetch.

2. **Duration-Limited (Finite Depth or Deep Water):** Calculates wave parameters when the wind event is too short for waves to become fully developed over the available fetch. This mode considers wind speed, storm duration, and water depth, and also calculates the equivalent fetch that would produce the same wave conditions.

The `calculator.py` script now comprehensively evaluates both fetch and duration limits for given inputs to determine the actual controlling wave growth factor and provides consistent results.

## Wind Speed ($U\_{10}$) and Adjusted Wind Speed ($U\_a$)

* $U\_{10}$ represents the average wind speed measured at a standard height of 10 meters above the mean water level. This is the primary raw wind input data for the model.

* The SMB calculations, particularly those recommended in the Shore Protection Manual (SPM 1984), require an "adjusted wind speed," $U\_a$. This adjustment is an empirical correction designed to account for the non-linear relationship between the measured wind speed and the actual wind stress (or friction velocity) at the water surface, which is the fundamental driver of wave growth.

* The specific formula for $U\_a$ from SPM (1984) is:
  $$U_a = 0.71 \cdot U_{10}^{1.23}$$
  where $U\_{10}$ is in meters per second (m/s). This $U\_a$ value is consistently employed in all subsequent wave prediction equations within these scripts.

## Methodology

The scripts are based on the semi-empirical Sverdrup-Munk-Bretschneider (SMB) method, which relates dimensionless wave parameters to wind conditions. The core principle is that wave growth is limited by either fetch (spatial constraint) or duration (temporal constraint). All formulas below utilize the **adjusted wind speed ($U\_a$)**.

### Deep Water Formulas (Fetch-Limited)

For deep water, the script uses the revised Bretschneider (1970) equations:

* **Dimensionless Fetch:**
  $$\hat{F} = \frac{gF}{U_a^{2}}$$

* **Significant Wave Height ($H\_s$):**
  $$\frac{gH_{s}}{U_a^{2}} = 0.283 \tanh\left[0.0125\left(\hat{F}\right)^{0.42}\right]$$

* **Significant Wave Period ($T\_s$):**
  $$\frac{gT_{s}}{U_a} = 7.54 \tanh\left[0.077\left(\hat{F}\right)^{0.25}\right]$$

* **Minimum Duration ($t\_{min}$):**
  $$t_{min} = \frac{U_a}{g} \cdot 6.5882 \cdot \exp\left(\sqrt{0.0161 \cdot (\ln \hat{F})^2 - 0.3692 \cdot \ln \hat{F} + 2.2024} + 0.8798 \cdot \ln \hat{F}\right)$$

### Depth-Limited Formulas (Fetch-Limited)

For shallower water, the script uses formulas from the Shore Protection Manual that incorporate a dimensionless depth parameter ($\hat{d} = gd/U\_a^2$):

* **Significant Wave Height ($H\_s$):**
  $$H_s = \frac{U_a^2}{g} \cdot 0.283 \cdot \tanh(0.530 \hat{d}^{0.75}) \cdot \tanh\left[\frac{0.00565 \hat{F}^{0.5}}{\tanh(0.530 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period ($T\_s$):**
  $$T_s = \frac{U_a}{g} \cdot 7.54 \cdot \tanh(0.833 \hat{d}^{0.375}) \cdot \tanh\left[\frac{0.0379 \hat{F}^{0.333}}{\tanh(0.833 \hat{d}^{0.375})}\right]$$

* **Minimum Duration ($t_{min}$):**
  $$t_{min} = \frac{U_a}{g} \cdot 6.5882 \cdot \exp\left(\sqrt{0.0161 \cdot (\ln \hat{F})^2 - 0.3692 \cdot \ln \hat{F} + 2.2024} + 0.8798 \cdot \ln \hat{F}\right)$$

### Duration-Limited Formulas (Finite Depth or Deep Water)

When the wind event duration ($t$) is the limiting factor, the following formulas are used.
Dimensionless duration: $\hat{t} = \frac{gt}{U\_a}$.

* **Significant Wave Height ($H\_s$) - Deep Water:**
  $$\frac{gH_{s}}{U_a^{2}} = 0.283 \tanh\left[0.000528 \left(\hat{t}\right)^{0.75}\right]$$

* **Significant Wave Period ($T\_s$) - Deep Water:**
  $$\frac{gT_{s}}{U_a} = 7.54 \tanh\left[0.00379 \left(\hat{t}\right)^{0.41}\right]$$

* **Significant Wave Height ($H\_s$) - Finite Depth (Heuristic):**
  $$H_s = \frac{U_a^2}{g} \cdot 0.283 \cdot \tanh(0.530 \hat{d}^{0.75}) \cdot \tanh\left[\frac{0.000528 \hat{t}^{0.75}}{\tanh(0.530 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period ($T\_s$) - Finite Depth (Heuristic):**
  $$T_s = \frac{U_a}{g} \cdot 7.54 \cdot \tanh(0.833 \hat{d}^{0.375}) \cdot \tanh\left[\frac{0.00379 \hat{t}^{0.41}}{\tanh(0.833 \hat{d}^{0.375})}\right]$$

## Features

* **Adjusted Wind Speed ($U\_a$) Integration:** All calculations now consistently use the adjusted wind speed ($U\_a$) derived from the 10-meter wind speed ($U\_{10}$), following SPM (1984) guidelines.

* **Dual-Mode Calculation:** Accurately applies formulas for fetch-limited (deep or finite depth) and duration-limited (deep or finite depth) conditions, and determines the controlling factor.

* **Comprehensive Outputs:** Calculates Significant Wave Height ($H\_s$), Significant Wave Period ($T\_s$), Minimum Storm Duration ($t\_{min}$) where applicable, and Equivalent Fetch for duration-limited cases.

* **Interactive Interface:** A simple command-line interface guides the user through the input process, allowing for flexible input of fetch, duration, and optional depth.

* **Validated Formulas:** The implemented equations are based on authoritative sources, including the U.S. Army's Coastal Engineering Manual and Shore Protection Manual.

* **Data Visualization:** Generates contour charts for both deep and depth-limited conditions, outputting to PDF.

* **Nomogram Generation:** Creates printable nomograms for quick graphical estimations, outputting to PDF.

* **Tabular Output:** Generates a comprehensive PDF table summarizing wave calculations for various parameters.

## Scripts Description

This repository contains the following Python scripts:

### `calculator.py`

This is the main interactive script for performing individual SMB wave calculations.

* **Functionality:**
    * Prompts the user for **10-meter wind speed ($U\_{10}$)**, fetch length, storm duration, and optionally water depth.
    * Calculates the **adjusted wind speed ($U\_a$)** from $U\_{10}$ and uses $U\_a$ in all subsequent wave parameter calculations.
    * Calculates wave parameters for both fetch-limited and duration-limited scenarios based on the provided inputs.
    * Determines the **controlling wave growth factor** (whether fetch or duration is the primary limit) and displays the corresponding significant wave height ($H\_s$), significant wave period ($T\_s$), and relevant duration/fetch values.
    * Outputs all calculations to a `report.txt` file, mirroring the command line output.

* **Usage:** Run directly from the command line and follow the prompts.
    ```bash
    python calculator.py
    ```

### `chart.py`

Generates a combined contour chart for SMB wave parameters in **deep water** conditions.

* **Functionality:**
    * Takes **10-meter wind speed ($U\_{10}$)** as input for the wind axis and calculates the corresponding adjusted wind speed ($U\_a$) for wave calculations.
    * Uses `matplotlib` to create a single plot showing contours of $H\_s$, $T\_s$, and $t_{min}$.
    * Displays wave parameters as functions of $U\_{10}$ and fetch length for deep water.
    * Utilizes different black line styles (solid for $H\_s$, dashed for $T\_s$, dotted for $t_{min}$) for clarity.
    * The chart is generated in A3 landscape format and saved as `combined_smb_chart_Ua.pdf`.

* **Usage:** Run directly to generate the PDF chart.
    ```bash
    python chart.py
    ```

### `chart_10m.py`

Generates a combined contour chart for SMB wave parameters in **depth-limited water** conditions, specifically for a fixed water depth of 10 meters.

* **Functionality:**
    * Similar to `chart.py`, but tailored for depth-limited scenarios.
    * Takes **10-meter wind speed ($U\_{10}$)** as input for the wind axis and calculates the corresponding adjusted wind speed ($U\_a$) for wave calculations.
    * Calculates and plots contours of $H\_s$, $T\_s$, and $t_{min}$ for a constant water depth (defaulting to 10m).
    * Provides a visual representation of how wave parameters change with $U\_{10}$ and fetch in shallow water.
    * The chart is generated in A3 landscape format and saved as `smb_chart_10m.pdf`.

* **Usage:** Run directly to generate the PDF chart. The `FIXED_DEPTH` variable can be modified within the script.
    ```bash
    python chart_10m.py
    ```

### `smb-nomogram-deep.py`

Generates a multi-page PDF containing three nomograms for **deep water** wave prediction.

* **Functionality:**
    * Takes **10-meter wind speed ($U\_{10}$)** as input for the wind axis and calculates the corresponding adjusted wind speed ($U\_a$) for wave calculations on the nomogram.
    * Creates separate nomograms for Significant Wave Height ($H\_s$), Significant Wave Period ($T\_s$), and Minimum required wind duration ($t_{min}$).
    * Outputs a single PDF file named `smb-nomogram-deep.pdf`.
    * Requires `pynomo` and `nomogen` libraries for nomogram generation.

* **Usage:** Run directly to generate the PDF.
    ```bash
    python smb-nomogram-deep.py
    ```

### `smb-nomogram-shallow.py`

Generates a multi-page PDF containing three nomograms for **depth-limited (shallow water)** wave prediction, configured for a fixed water depth of 10 meters.

* **Functionality:**
    * Similar to `smb-nomogram-deep.py`, but specifically for shallow water conditions.
    * Takes **10-meter wind speed ($U\_{10}$)** as input for the wind axis and calculates the corresponding adjusted wind speed ($U\_a$) for wave calculations on the nomogram.
    * Generates nomograms for $H\_s$, $T\_s$, and $t_{min}$ at a constant water depth (defaulting to 10m).
    * Outputs a single PDF file named `smb-nomogram-shallow.pdf`.
    * Requires `pynomo` and `nomogen` libraries.

* **Usage:** Run directly to generate the PDF. The `WATER_DEPTH` variable can be modified within the script.
    ```bash
    python smb-nomogram-shallow.py
    ```

### `tables.py`

Generates a comprehensive PDF table summarizing SMB wave calculations for various combinations of wind speed, fetch, and depth.

* **Functionality:**
    * Calculates $H\_s$, $T\_s$, and $t_{min}$ for predefined ranges of **10-meter wind speeds ($U\_{10}$)** (10-30 m/s), fetches (5-50 km), and depths (Deep Water, 50m, 25m, 10m, 5m).
    * Calculates the **adjusted wind speed ($U\_a$)** from $U\_{10}$ for all wave calculations within the table.
    * Organizes the results into a well-formatted table within a PDF document.
    * Uses `reportlab` for PDF generation, ensuring a professional and readable output.

* **Usage:** Run directly to generate the PDF file `comprehensive_wave_calculations.pdf`.
    ```bash
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
    * `nomogen.py` is expected to be in same directory as nomogram scripts

    You can install most of them using pip:
    ```bash
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
    * **Wind Speed:** Should be the standard 10-meter overwater wind speed ($U\_{10}$). The model internally adjusts this to $U\_a$ for calculations.
    * **Fetch Length:** Should be the "effective fetch," which accounts for the geometry of the water body, not just a straight-line distance.

* **Heuristic for Duration-Limited (Finite Depth):** The formulas used for duration-limited conditions in finite depth (Option 2 in `calculator.py`) are an adaptation based on the structure of SMB equations for fetch-limited finite depth. Direct empirical formulas for this specific combined scenario are less common in basic SMB literature. While providing a reasonable estimate, these should be used with awareness of their heuristic nature.

## Bibliography

* U.S. Army. (2008). *Coastal Engineering Manual* (EM 1110-2-1100). Washington, DC: U.S. Army Corps of Engineers. URL: https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/u43544q/636F617374616C20656E67696E656572696E67206D616E75616C/
* U.S. Army. (1984). *Shore Protection Manual*. Vicksburg, MS: U.S. Army Engineer Waterways Experiment Station. URL: https://archive.org/details/shoreprotectionm01unit & https://archive.org/details/shoreprotectionm02unit/page/n387/mode/2up
* Bishop et al. (1992). Shore protection manual's wave prediction reviewed. Coastal Engineering, Volume 17, Issues 1–2, 1992, Pages 25-48, ISSN 0378-3839, https://doi.org/10.1016/0378-3839(92)90012-J.
* Etemad-Shahidi, A., Kazeminezhad, M. H., & Mousavi, S. J. (2009). On the prediction of wave parameters using simplified methods. *Journal of Coastal Research, SI 56*, 505-509. URL: https://www.jstor.org/stable/25737628
* Bretschneider, C. L. (1970). *Wave forecasting relations for wave generation*. Look Lab, Hawaii, 1(3).
