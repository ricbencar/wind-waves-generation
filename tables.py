import math
from reportlab.lib.pagesizes import A4 # Import A4 for portrait
from reportlab.platypus import SimpleDocTemplate, Table, TableStyle, Paragraph, Spacer
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet

# =============================================================================
# SMB Wave Prediction Model Functions
# =============================================================================
# These functions are extracted from the provided calculator.py to ensure
# consistency in wave parameter calculations for both deep and depth-limited conditions.

def calculate_deep_water(wind_speed, fetch, gravity=9.81):
    """
    Calculates fetch-limited wave properties in deep water using the SMB method.

    Args:
        wind_speed (float): The wind speed at 10m height over water (m/s).
        fetch (float):      The effective fetch length (m).
        gravity (float):    The acceleration due to gravity (m/s^2).

    Returns:
        tuple: A tuple containing:
            - Hs (float): Predicted significant wave height (m).
            - Ts (float): Predicted significant wave period (s).
            - t_min (float): Minimum wind duration for fetch-limited state (hours).
    """
    # Handle the case where fetch is zero, as log(0) is undefined.
    # For zero fetch, no waves are generated, so Hs, Ts, and Duration are 0.
    if fetch <= 0:
        return 0.0, 0.0, 0.0

    # Dimensionless Fetch Calculation: F_hat = g * F / U^2
    dim_fetch = (gravity * fetch) / (wind_speed**2)

    # Significant Wave Height (Hs) Calculation: (g * Hs) / U^2 = 0.283 * tanh[0.0125 * (g * F / U^2)^0.42]
    gHs_U2 = 0.283 * math.tanh(0.0125 * (dim_fetch**0.42))
    Hs = gHs_U2 * (wind_speed**2 / gravity)

    # Significant Wave Period (Ts) Calculation: (g * Ts) / U = 7.54 * tanh[0.077 * (g * F / U^2)^0.25]
    gTs_U = 7.54 * math.tanh(0.077 * (dim_fetch**0.25))
    Ts = gTs_U * (wind_speed / gravity)

    # Minimum Wind Duration (t_min) Calculation: Empirical formula
    log_dim_fetch = math.log(dim_fetch) # This line caused the error when dim_fetch was 0
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed / gravity
    t_min_hours = t_min_seconds / 3600  # Convert to hours

    return Hs, Ts, t_min_hours

def calculate_depth_limited(wind_speed, fetch, depth, gravity=9.81):
    """
    Calculates wave properties for depth-limited conditions.

    This function uses SMB-based formulas adapted for shallow or transitional
    depths, where wave growth is influenced by the seabed.

    Args:
        wind_speed (float): Wind speed at 10m height (m/s).
        fetch (float):      Effective fetch length (m).
        depth (float):      Water depth (m).
        gravity (float):    Acceleration of gravity (m/s^2).

    Returns:
        tuple: A tuple containing Hs (m), Ts (s), and t_min (hours).
    """
    # Handle the case where fetch is zero, as log(0) is undefined in duration calculation.
    # For zero fetch, no waves are generated, so Hs, Ts, and Duration are 0.
    if fetch <= 0:
        return 0.0, 0.0, 0.0

    # Dimensionless Parameters
    dim_fetch = (gravity * fetch) / wind_speed**2
    dim_depth = (gravity * depth) / wind_speed**2

    # Depth-Limited Significant Wave Height (Hs)
    term_h = 0.00565 * (dim_fetch)**0.5
    tanh_depth_h = math.tanh(0.530 * (dim_depth)**0.75)
    Hs = (wind_speed**2 / gravity) * 0.283 * tanh_depth_h * math.tanh(term_h / tanh_depth_h)

    # Depth-Limited Significant Wave Period (Ts)
    term_t = 0.0379 * (dim_fetch)**0.333
    tanh_depth_t = math.tanh(0.833 * (dim_depth)**0.375)
    Ts = (wind_speed / gravity) * 7.54 * tanh_depth_t * math.tanh(term_t / tanh_depth_t)

    # Minimum Wind Duration (t_min) Calculation (same as deep water for consistency with original script)
    # This part was the source of the error if dim_fetch was 0
    log_dim_fetch = math.log(dim_fetch)
    A, B, C, D = 0.0161, 0.3692, 2.2024, 0.8798
    exponent_term = (A * log_dim_fetch**2 - B * log_dim_fetch + C)**0.5 + D * log_dim_fetch
    gt_min_U = 6.5882 * math.exp(exponent_term)
    t_min_seconds = gt_min_U * wind_speed / gravity
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
    headers = ["Wind (m/s)", "Fetch (km)", "Depth (m)", "Hs (m)", "Ts (s)", "Dur (h)"]

    # Define the ranges for calculations
    fetches_km = list(range(0, 51, 5))  # 0 to 50 km, step 5 km
    wind_speeds_mps = list(range(5, 36, 5)) # 5 to 35 m/s, step 5 m/s
    depths_m = [999, 100, 50, 25, 10, 5, 1] # Specified depths

    data_rows = []

    # Iterate through all combinations to generate data
    for wind_speed in wind_speeds_mps:
        for fetch_km in fetches_km:
            fetch_m = float(fetch_km * 1000)  # Convert fetch from km to meters, ensure float
            for depth in depths_m:
                hs, ts, duration = 0.0, 0.0, 0.0 # Initialize with floats

                if depth == 999: # Treat 999m as deep water for calculation
                    hs, ts, duration = calculate_deep_water(wind_speed, fetch_m)
                else: # All other specified depths are depth-limited
                    hs, ts, duration = calculate_depth_limited(wind_speed, fetch_m, float(depth)) # Ensure depth is float

                data_rows.append([
                    f"{wind_speed:.1f}",
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
