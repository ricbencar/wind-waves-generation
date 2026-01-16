# SMB Wave Prediction Model

A collection of Python scripts for predicting wind-generated wave characteristics using the Sverdrup-Munk-Bretschneider (SMB) method. These tools calculate significant wave height ($H_s$), significant wave period ($T_s$), and the minimum required storm duration based on meteorological inputs.

## Overview

This repository provides a practical implementation of the foundational SMB wave prediction model, a cornerstone of coastal and ocean engineering for decades. It is designed for engineers, scientists, and students who need to make preliminary estimates of wave conditions for design, planning, or research purposes.

The core model operates in two distinct modes, based on the primary limiting condition for wave growth:

1. **Fetch-Limited (Finite Depth or Deep Water):** For scenarios where wave growth is primarily limited by the available fetch. This mode can handle both deep water conditions (where the seabed does not influence wave generation) and finite depth conditions (where water depth is a significant factor). It calculates $H_s$, $T_s$, and the minimum duration required for a fully developed sea state over the given fetch.

2. **Duration-Limited (Finite Depth or Deep Water):** Calculates wave parameters when the wind event is too short for waves to become fully developed over the available fetch. This mode considers wind speed, storm duration, and water depth, and also calculates the equivalent fetch that would produce the same wave conditions.

The `calculator.py` script comprehensively evaluates both fetch and duration limits for given inputs to determine the actual controlling wave growth factor and provides consistent results.

## Wind Speed ($U_{10}$) and Adjusted Wind Speed ($U_a$)

* $U_{10}$ represents the average wind speed measured at a standard height of 10 meters above the mean water level. This is the primary raw wind input data for the model.

* The SMB calculations, particularly those recommended in the Shore Protection Manual (SPM 1984), require an "adjusted wind speed," $U_a$. This adjustment is an empirical correction designed to account for the non-linear relationship between the measured wind speed and the actual wind stress (or friction velocity) at the water surface, which is the fundamental driver of wave growth.

* The specific formula for $U_a$ from SPM (1984) is:
  $$U_a = 0.71 \cdot U_{10}^{1.23}$$
  where $U_{10}$ is in meters per second (m/s). This $U_a$ value is consistently employed in all subsequent wave prediction equations within these scripts.

## Methodology

The scripts are based on the semi-empirical Sverdrup-Munk-Bretschneider (SMB) method. Specifically, this implementation uses the **Revised SMB Model proposed by Hurdle & Stive (1989)**. This formulation provides a single set of unified equations that work smoothly across deep and shallow water, correcting inconsistencies found in the original Shore Protection Manual (SPM 1984) equations.

All formulas below utilize the **adjusted wind speed ($U_a$)**.

### Unified Fetch-Limited Formulas (Hurdle & Stive, 1989)

These equations apply to both deep water and finite depth conditions. The depth terms naturally approach 1.0 as depth increases, converging to the deep-water asymptotes.

* **Dimensionless Fetch:**
  $$\hat{F} = \frac{gF}{U_a^{2}}$$

* **Dimensionless Depth:**
  $$\hat{d} = \frac{gd}{U_a^{2}}$$

* **Significant Wave Height ($H_s$):**
  $$\frac{gH_{s}}{U_a^{2}} = 0.25 \tanh(0.6 \hat{d}^{0.75}) \tanh^{0.5}\left[\frac{4.3 \times 10^{-5} \hat{F}}{\tanh^2(0.6 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period ($T_s$):**
  $$\frac{gT_{s}}{U_a} = 8.3 \tanh(0.76 \hat{d}^{0.375}) \tanh^{1/3}\left[\frac{4.1 \times 10^{-5} \hat{F}}{\tanh^3(0.76 \hat{d}^{0.375})}\right]$$

* **Minimum Duration ($t_{min}$):**
  $$t_{min} = \frac{U_a}{g} \cdot 6.5882 \cdot \exp\left(\sqrt{0.0161 \cdot (\ln \hat{F})^2 - 0.3692 \cdot \ln \hat{F} + 2.2024} + 0.8798 \cdot \ln \hat{F}\right)$$
  *(Note: $t_{min}$ calculation retains the approximation by Etemad-Shahidi et al., 2009)*

### Duration-Limited Formulas (Finite Depth or Deep Water)

When the wind event duration ($t$) is the limiting factor, the coefficients have been adapted to ensure asymptotic consistency with the Hurdle & Stive fetch-limited equations.
Dimensionless duration: $\hat{t} = \frac{gt}{U_a}$.

* **Significant Wave Height ($H_s$) - Deep Water:**
  $$\frac{gH_{s}}{U_a^{2}} = 0.25 \tanh\left[0.000528 \left(\hat{t}\right)^{0.75}\right]$$

* **Significant Wave Period ($T_s$) - Deep Water:**
  $$\frac{gT_{s}}{U_a} = 8.3 \tanh\left[0.00379 \left(\hat{t}\right)^{0.41}\right]$$

* **Significant Wave Height ($H_s$) - Finite Depth (Heuristic):**
  $$H_s = \frac{U_a^2}{g} \cdot 0.25 \cdot \tanh(0.6 \hat{d}^{0.75}) \cdot \tanh\left[\frac{0.000528 \hat{t}^{0.75}}{\tanh(0.6 \hat{d}^{0.75})}\right]$$

* **Significant Wave Period ($T_s$) - Finite Depth (Heuristic):**
  $$T_s = \frac{U_a}{g} \cdot 8.3 \cdot \tanh(0.76 \hat{d}^{0.375}) \cdot \tanh\left[\frac{0.00379 \hat{t}^{0.41}}{\tanh(0.76 \hat{d}^{0.375})}\right]$$

## Features

* **Adjusted Wind Speed ($U_a$) Integration:** All calculations consistent across the suite (Python and C++) use the **Adjusted Wind Speed ($U_a$)** derived from the 10-meter wind speed ($U_{10}$), strictly following Shore Protection Manual (SPM 1984) guidelines ($U_a = 0.71 \cdot U_{10}^{1.23}$).

* **Dual-Mode Physics Engine:** The core logic accurately distinguishes between **Fetch-Limited** and **Duration-Limited** growth regimes. It computes potentials for both scenarios (checking against Minimum Duration $t_{min}$) and automatically determines the **Controlling Factor** to output the realistic sea state.

* **Multi-Regime Support:** The model supports calculations for both **Deep Water** and **Finite Depth**. The unified Hurdle & Stive (1989) equations ensure smooth mathematical transitions between regimes without the "step changes" found in older models.

* **Advanced Hydrodynamics (GUI Only):** The C++ implementation includes a **Newton-Raphson solver** for the transcendental Linear Dispersion Relation ($L = \frac{gT^2}{2\pi}\tanh(\frac{2\pi d}{L})$), providing exact outputs for Wavelength ($L$), Celerity ($C$), and Wave Number ($k$) rather than deep-water approximations.

* **Stability Analysis (GUI Only):** The GUI evaluates the **Miche Criterion** ($H/L \le 0.142 \tanh(kd)$) to determine the breaking status of the waves, warning if the sea state is unstable due to steepness.

* **Interactive Interfaces:**
    * **Python CLI:** A robust command-line tool for iterative text-based calculations.
    * **C++ GUI:** A high-performance, multi-threaded Windows application for responsive analysis.

* **Professional Visualization:**
    * **Contour Charts:** Generates high-quality A3 PDF contour plots for $H_s$, $T_s$, and $t_{min}$.
    * **Nomograms:** Creates precision alignment charts (N-charts) using numerical optimization for graphical estimation.
    * **Tabular Reports:** Generates comprehensive PDF lookup tables for wide ranges of wind/fetch/depth.

## Scripts Description

This repository contains the following computational modules:

### `calculator.py`

The primary interactive Python script for single-point SMB wave calculations.

* **Functionality:**
    * Prompts for **10-meter wind speed ($U_{10}$)**, fetch length, storm duration, and optional water depth.
    * Internally converts inputs to **Adjusted Wind Speed ($U_a$)**.
    * Calculates potentials for both **Fetch-Limited** and **Duration-Limited** scenarios using the updated Hurdle & Stive logic.
    * Logic: Determines the **Controlling Factor** by taking the minimum of the fetch-limited and duration-limited potentials.
    * For duration-limited cases, calculates the **Equivalent Fetch** (the fetch required to produce the same energy).
    * Outputs a detailed transcript to the console and logs it to `report.txt`.
* **Usage:**
    ```bash
    python calculator.py
    ```

### `calculator_gui.cpp`

A high-performance **Windows GUI application** written in C++ (Win32 API).

* **Functionality:**
    * **Real-time Interface:** graphical inputs for $U_{10}$, Fetch, Duration, and Depth.
    * **Advanced Physics:** Unlike the Python scripts, this module solves the **Linear Dispersion Relation** using a Newton-Raphson iterative solver to compute exact Wavelength ($L$) and Phase Speed ($C$).
    * **Regime Analysis:** Automatically detects and displays the flow regime (Deep, Shallow, or Transitional) based on $d/L$.
    * **Breaking Check:** Calculates wave steepness ($H/L$) and checks against the Miche Limit to report if waves are "STABLE" or "BREAKING".
    * **Architecture:** Uses a worker thread to perform calculations without freezing the UI.
* **Usage:** Must be compiled (see "How to Use").

### `smb_chart_deep.py`

Generates a combined contour chart for SMB wave parameters in **deep water** conditions.

* **Functionality:**
    * Plotting Axis: **10-meter wind speed ($U_{10}$)** vs. Fetch Length.
    * **Contours:** Plots three distinct parameter sets on one chart:
        * **Significant Wave Height ($H_s$):** Solid black lines.
        * **Significant Wave Period ($T_s$):** Dashed black lines.
        * **Minimum Duration ($t_{min}$):** Dotted black lines.
    * **Output:** Generates an A3 landscape PDF named `smb_chart_deep.pdf`.
* **Usage:**
    ```bash
    python smb_chart_deep.py
    ```

### `smb_chart_10m.py`

Generates a combined contour chart for SMB wave parameters in **depth-limited water**.

* **Functionality:**
    * **Fixed Depth:** Configured for a constant water depth (default: **10 meters**).
    * Applies the unified Hurdle & Stive depth-damping ($\tanh$) functions to all wave growth curves.
    * Visualizes how shallow water friction limits maximum wave height compared to the deep water chart.
    * **Output:** Generates an A3 landscape PDF named `smb_chart_10m.pdf`.
* **Usage:**
    ```bash
    python smb_chart_10m.py
    ```

### `smb-nomogram-deep.py`

Generates a professional-grade nomogram (alignment chart) for **deep water** prediction.

* **Functionality:**
    * Uses the `nomogen` library to numerically optimize scale positions.
    * Produces a multi-page PDF containing separate nomograms for $H_s$, $T_s$, and $t_{min}$.
    * **Isopleth:** Automatically plots a reference example line (Lake Garda scenario: $U_{10}=25$ m/s, Fetch=45 km).
    * **Output:** `smb-nomogram-deep.pdf`.
* **Usage:**
    ```bash
    python smb-nomogram-deep.py
    ```

### `smb-nomogram-shallow.py`

Generates nomograms for **depth-limited** conditions.

* **Functionality:**
    * Configured for a fixed depth of **10 meters**.
    * Solves the complex implicit relationships of the unified shallow water equations to build linear scales.
    * Includes the Lake Garda isopleth example adjusted for shallow depth.
    * **Output:** `smb-nomogram-shallow.pdf`.
* **Usage:**
    ```bash
    python smb-nomogram-shallow.py
    ```

### `tables.py`

Generates a comprehensive look-up table PDF.

* **Functionality:**
    * Iterates through a matrix of:
        * **Wind Speeds ($U_{10}$):** 10 to 30 m/s.
        * **Fetches:** 5 to 50 km.
        * **Depths:** 5m, 10m, 25m, 50m, and Deep Water.
    * Uses `reportlab` to create a formatted, readable PDF report.
    * **Output:** `comprehensive_wave_calculations.pdf`.
* **Usage:**
    ```bash
    python tables.py
    ```

### `nomogen.py`

A utility library required by the nomogram scripts.

* **Functionality:**
    * Implements numerical optimization (using `scipy` and `numpy`) to construct nomogram scales (Type 9 / N-charts) for arbitrary functions.
    * Handles transformations, scale generation, and LaTeX text rendering via `PyX`.

## How to Use

1. **Python Prerequisites:** Ensure you have Python 3 installed. You will need the following scientific libraries:

    * `numpy`
    * `scipy` (Note: Scripts include a fix for `scipy.arange` compatibility)
    * `matplotlib` (for charts)
    * `reportlab` (for tables)
    * `pynomo` (for nomograms)
    * `pypdf` (for merging PDFs)
    * `pyx` (for nomogram rendering)

    Install via pip:
    ```bash
    pip install numpy scipy matplotlib reportlab pypdf pyx pynomo
    ```

2. **C++ Compilation (Optional GUI):**
    To use the GUI, compile `calculator_gui.cpp` using a compiler like MinGW (G++). The following command ensures static linking and includes necessary Windows libraries:

    ```bash
    g++ -O3 -std=c++17 -static -static-libgcc -static-libstdc++ -o calculator_gui.exe calculator_gui.cpp -mwindows -lgdi32
    ```

3. **Running Scripts:**

    * **Interactive Mode:** `python calculator.py`
    * **Charts:** `python smb_chart_deep.py` or `python smb_chart_10m.py`
    * **Nomograms:** `python smb-nomogram-deep.py`
    * **Tables:** `python tables.py`

## Assumptions and Limitations

* **Steady-State Wind:** The model assumes that the wind speed and direction are uniform and constant across the entire fetch for the specified duration. This is an idealization not always met in nature.

* **Input Data Quality:** The accuracy of the results is highly dependent on the quality of the inputs. For best results:
    * **Wind Speed:** Should be the standard 10-meter overwater wind speed ($U_{10}$). The model internally adjusts this to $U_a$ for calculations.
    * **Fetch Length:** Should be the "effective fetch," which accounts for the geometry of the water body, not just a straight-line distance.

* **Heuristic for Duration-Limited (Finite Depth):** The formulas used for duration-limited conditions in finite depth (Option 2 in `calculator.py`) are an adaptation based on the structure of SMB equations for fetch-limited finite depth. Direct empirical formulas for this specific combined scenario are less common in basic SMB literature. While providing a reasonable estimate, these should be used with awareness of their heuristic nature.

## Bibliography

### Primary Engineering Manuals (Operational Standards)
* **U.S. Army.** (2008). *Coastal Engineering Manual* (EM 1110-2-1100). Washington, DC: U.S. Army Corps of Engineers.
    * The current primary reference for USACE coastal projects, superseding the Shore Protection Manual.
    * **URL:** https://www.publications.usace.army.mil/USACE-Publications/Engineer-Manuals/

* **U.S. Army.** (1984). *Shore Protection Manual* (Vol. 1 & 2). Vicksburg, MS: U.S. Army Engineer Waterways Experiment Station.
    * The classic reference that standardized the SMB equations for decades.
    * **URL (Vol 1):** https://usace.contentdm.oclc.org/digital/collection/p16021coll11/id/1934/

* **World Meteorological Organization (WMO).** (2018). *Guide to Wave Analysis and Forecasting* (WMO-No. 702). Geneva: Secretariat of the WMO.
    * The international standard for meteorological wave forecasting.
    * **URL:** https://library.wmo.int/records/item/31871-guide-to-wave-analysis-and-forecasting

### Foundational Papers (The "SMB" Method)
* **Sverdrup, H. U., & Munk, W. H.** (1947). *Wind, Sea, and Swell: Theory of Relations for Forecasting*. H.O. Pub. No. 601, U.S. Navy Hydrographic Office, Washington, D.C.
    * The original wartime research that established the "Significant Wave" concept.
    * **URL:** https://ia801302.us.archive.org/32/items/windseaswelltheo00sver/windseaswelltheo00sver.pdf

* **Bretschneider, C. L.** (1958). *Revisions in Wave Forecasting: Deep and Shallow Water*. Proceedings of the 6th Conference on Coastal Engineering, pp. 30–67.
    * The paper that revised the 1947 Sverdrup-Munk curves, adding the "B" to SMB.
    * **URL:** https://icce-ojs-tamu.tdl.org/icce/article/view/1868

* **Bretschneider, C. L.** (1970). *Wave forecasting relations for wave generation*. Look Lab, Hawaii, 1(3).
    * Further refinement of the 1958 curves.

### Revisions, Critiques & Finite Depth Extensions
* **Hurdle, D. P., & Stive, R.** (1989). *Revision of SPM 1984 wave hindcast model to avoid inconsistencies in engineering applications*. Coastal Engineering, 12(4), 339–351.
    * Corrects inconsistencies in the original 1984 SPM formulas.
    * **DOI:** https://doi.org/10.1016/0378-3839(89)90011-2

* **Bishop, C. T., Donelan, M. A., & Kahma, K. K.** (1992). *Shore protection manual's wave prediction reviewed*. Coastal Engineering, 17(1-2), 25-48.
    * Comprehensive comparison of SPM predictions against measured data.
    * **DOI:** https://doi.org/10.1016/0378-3839(92)90012-J

* **Etemad-Shahidi, A., Kazeminezhad, M. H., & Mousavi, S. J.** (2009). *On the prediction of wave parameters using simplified methods*. Journal of Coastal Research, SI 56, 505-509.
    * Validates the SMB equations against modern methods.
    * **URL:** https://www.jstor.org/stable/25737628