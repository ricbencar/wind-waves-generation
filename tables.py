# =============================================================================
# SMB Wave Prediction Model Functions
# =============================================================================

import math
from reportlab.lib.pagesizes import A4 # Import A4 for portrait
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet

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
    # Handle the case where fetch is zero
    if fetch <= 0:
        return 0.0, 0.0, 0.0

    # --- Dimensionless Fetch Calculation ---
    # F_hat = g * F / Ua^2
    dim_fetch = (gravity * fetch) / (wind_speed_adjusted**2)

    # --- Significant Wave Height (Hs) Calculation (Hurdle & Stive, 1989) ---
    # Deep water limit of Eq 4.1: tanh(depth terms) -> 1
    # Revised Formula: (g * Hs) / Ua^2 = 0.25 * [tanh(4.3e-5 * F_hat)]^0.5
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
    # Ref: Etemad-Shahidi et al. (2009), based on SPM data.
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours

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
    # Handle the case where fetch is zero
    if fetch <= 0:
        return 0.0, 0.0, 0.0

    # --- Dimensionless Parameters ---
    dim_fetch = (gravity * fetch) / wind_speed_adjusted**2
    dim_depth = (gravity * depth) / wind_speed_adjusted**2

    # --- Revised Significant Wave Height (Hs) ---
    # Hurdle & Stive (1989), Eq 4.1
    # Coefficient is 0.25 (vs 0.283 in SPM)
    # Depth term: tanh(0.6 * d_hat^0.75)
    # Fetch term inner: 4.3e-5 * F_hat / (depth_term^2)
    depth_term_h = math.tanh(0.6 * dim_depth**0.75)
    fetch_term_inner_h = (4.3e-5 * dim_fetch) / (depth_term_h**2)
    
    gHs_U2 = 0.25 * depth_term_h * (math.tanh(fetch_term_inner_h))**0.5
    Hs = gHs_U2 * (wind_speed_adjusted**2 / gravity)

    # --- Revised Significant Wave Period (Ts) ---
    # Hurdle & Stive (1989), Eq 4.2
    # Coefficient is 8.3 (vs 7.54 in SPM)
    # Depth term: tanh(0.76 * d_hat^0.375)
    # Fetch term inner: 4.1e-5 * F_hat / (depth_term^3)
    depth_term_t = math.tanh(0.76 * dim_depth**0.375)
    fetch_term_inner_t = (4.1e-5 * dim_fetch) / (depth_term_t**3)
    
    gTs_U = 8.3 * depth_term_t * (math.tanh(fetch_term_inner_t))**(1/3)
    Ts = gTs_U * (wind_speed_adjusted / gravity)

    # --- Minimum Wind Duration (t_min) Calculation ---
    # Ref: Etemad-Shahidi et al. (2009), based on SPM data.
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed_adjusted / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours

    return Hs, Ts, t_min_hours

def create_comprehensive_wave_table_pdf(output_filename="comprehensive_wave_calculations.pdf"):
    """
    Generates a PDF document with a comprehensive table of wave calculations
    for various wind speeds, fetches, and depths.
    """
    # Use A4 for portrait orientation
    doc = SimpleDocTemplate(output_filename, pagesize=A4)
    styles = getSampleStyleSheet()
    story = []

    # Add a title and introduction to the document
    title_text = "SMB (Sverdrup-Munk-Bretschneider)<br/>Wave Prediction Model"
    story.append(Paragraph(title_text, styles['h1']))
    story.append(Spacer(1, 12)) # Add some space

    intro_text = "This table presents significant wave height (Hs), significant wave period (Ts), and minimum storm duration for various combinations of wind speed, fetch, and water depth, calculated using the SMB model."
    story.append(Paragraph(intro_text, styles['Normal']))
    story.append(Spacer(1, 24)) # Add more space before the table

    # Define the table headers
    headers = ["U10 (m/s)", "Fetch (km)", "Depth (m)", "Hs (m)", "Ts (s)", "Dur (h)"] # Changed header to U10

    # Define the ranges for calculations
    fetches_km = list(range(5, 51, 5))  # 5 to 50 km, step 5 km
    U10_speeds_mps = list(range(10, 31, 5)) # 10 to 30 m/s, step 5 m/s
    depths_m = [5, 10, 25, 50, 1000] # Specified depths

    data_rows = []

    # Iterate through all combinations to generate data
    for U10_speed in U10_speeds_mps:
        Ua_speed = calculate_adjusted_wind_speed(float(U10_speed)) # Calculate Ua
        for fetch_km in fetches_km:
            fetch_m = float(fetch_km * 1000)  # Convert fetch from km to meters, ensure float
            for depth in depths_m:
                hs, ts, duration = 0.0, 0.0, 0.0 # Initialize with floats

                if depth == 999: # Treat 999m as deep water for calculation
                    hs, ts, duration = calculate_deep_water(Ua_speed, fetch_m)
                else: # All other specified depths are depth-limited
                    hs, ts, duration = calculate_depth_limited(Ua_speed, fetch_m, float(depth)) # Ensure depth is float

                data_rows.append([
                    f"{U10_speed:.1f}", # Use U10 for display
                    f"{fetch_km:.0f}",
                    f"{depth:.0f}",
                    f"{hs:.2f}",
                    f"{ts:.2f}",
                    f"{duration:.2f}" if not math.isnan(duration) else "N/A"
                ])

    # Combine headers and data for the table
    table_data = [headers] + data_rows

    # Create the table object, with repeatRows=1 to repeat the first row (headers) on each page
    table = Table(table_data, repeatRows=1)

    # Define table style for good readability and a "fancy" look
    table_style = TableStyle([
        # Header Styling
        ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#4A4A4A')), # Dark grey header
        ('TEXTCOLOR', (0, 0), (-1, 0), colors.whitesmoke),          # Header text color
        ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),            # Header font
        ('FONTSIZE', (0, 0), (-1, 0), 14),                          # Larger font size for header
        ('BOTTOMPADDING', (0, 0), (-1, 0), 14),                     # More padding below header

        # General Cell Styling
        ('ALIGN', (0, 0), (-1, -1), 'CENTER'),                      # Center align all content
        ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),                # Font for data cells
        ('FONTSIZE', (0, 1), (-1, -1), 14),                         # Increased font size for data cells (numbers)
        ('LEFTPADDING', (0, 0), (-1, -1), 8),                       # Increased left padding
        ('RIGHTPADDING', (0, 0), (-1, -1), 8),                      # Increased right padding
        ('TOPPADDING', (0, 0), (-1, -1), 8),                        # Increased top padding
        ('BOTTOMPADDING', (0, 0), (-1, -1), 8),                     # Increased bottom padding

        # Grid Lines
        ('GRID', (0, 0), (-1, -1), 1, colors.HexColor('#CCCCCC')),  # Lighter grey grid lines
        ('BOX', (0, 0), (-1, -1), 1.5, colors.HexColor('#4A4A4A')), # Thicker box border around the entire table (dark grey)
    ])

    # Apply alternating row colors for better readability with more contrast
    for i in range(1, len(table_data)):
        if i % 2 == 0:
            table_style.add('BACKGROUND', (0, i), (-1, i), colors.white) # White background
        else:
            table_style.add('BACKGROUND', (0, i), (-1, i), colors.lightgrey) # Light grey background

    table.setStyle(table_style)

    # Add the table to the story
    story.append(table)

    # Build the PDF document
    try:
        doc.build(story)
        print(f"PDF '{output_filename}' created successfully with comprehensive wave calculation table.")
    except Exception as e:
        print(f"Error creating PDF: {e}")

# Main execution block
if __name__ == "__main__":
    create_comprehensive_wave_table_pdf("tables.pdf")