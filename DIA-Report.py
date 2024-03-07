# python3 Report.py _features.csv .stats.tsv .features.tsv  _result.tsv .tsv .pdf
# This script is Features.py + DiaNN-Parser.py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import io
import csv
import sys
import os
import glob

# Input file 
FeaturesFile = sys.argv[1]
StatsFile = sys.argv[2]
AdditionalDataFile = sys.argv[3]
result_tsv = sys.argv[4]
tsv = sys.argv[5]
Report = sys.argv[6]

# Initialize PDF
c = canvas.Canvas(Report, pagesize=letter)
width, height = letter

################ Processing Data for Plots from Script One ####################
base_name = os.path.splitext(os.path.basename(FeaturesFile))[0]

# Initialize lists for plotting
RetentionTime_values_high_charge = []  # For charge >= 2
MZ_values_high_charge = []  # For charge >= 2
RetentionTime_values = []  # For all features
MZ_values = []  # For all features

# Open and read the additional TSV file for mz and retention time values
with open(AdditionalDataFile, 'r') as file:
    tsv_reader = csv.DictReader(file, delimiter='\t')
    for row in tsv_reader:
        mz = float(row['mz'])
        rtApex = float(row['rtApex'])
        charge = int(row['charge'])
         # Store all features
        MZ_values.append(mz)
        RetentionTime_values.append(rtApex)

        if charge >= 2:
            MZ_values_high_charge.append(mz)
            RetentionTime_values_high_charge.append(rtApex)

lower_bounds = []
upper_bounds = []
window_sizes = []
injection_times = []
unique_targets = set()
scan_number_MSMS = 0
scan_number_MS1 = 0

with open(FeaturesFile) as file:
    csv_reader = csv.DictReader(file)  # Use DictReader to access columns by header
    for row in csv_reader:
        # Check and count non-'N/A' values for 'Precursor Ion m/z'
        if row.get('Precursor Ion m/z') and row['Precursor Ion m/z'] != 'N/A':
            scan_number_MSMS += 1
        # Check and count non-'N/A' values for 'MS1 Base Peak m/z'
        if row.get('MS1 Base Peak m/z') and row['MS1 Base Peak m/z'] != 'N/A':
            scan_number_MS1 += 1

        if row['Injection Time'] != 'N/A':
            injection_times.append(float(row['Injection Time']))
        # Calculate window size
        if row['Isolation Window Target'] != 'N/A' and row['Isolation window upper offset'] != 'N/A' and row['Isolation window lower offset'] != 'N/A':
            target = float(row['Isolation Window Target'])
            upper_offset = float(row['Isolation window upper offset'])
            lower_offset = float(row['Isolation window lower offset'])
            
            lower_bound = target - lower_offset
            upper_bound = target + upper_offset
            lower_bounds.append(lower_bound)
            upper_bounds.append(upper_bound)
            window_size = upper_offset + lower_offset
            window_sizes.append(window_size)
            unique_targets.add(target)

# Calculate average and median
average_injection_time = np.mean(injection_times)
median_injection_time = np.median(injection_times)
average_lower_bound = np.mean(lower_bounds)
average_upper_bound = np.mean(upper_bounds)
average_window_size = np.mean(window_sizes)

# Initialize variables for the second file
precursors_identified = 0
# Read data from the second TSV file
with open(StatsFile, newline='') as file:
    tsv_reader = csv.DictReader(file, delimiter='\t')  # TSV files use tab delimiters
    for row in tsv_reader:
        precursors_identified = row['Precursors.Identified']

################ Processing Data for Plots for DIA files ####################
df1 = pd.read_csv(result_tsv, sep='\t')
df2 = pd.read_csv(tsv, sep='\t')
# Merging two files
merged_df = pd.merge(df1, df2, left_on='ModifiedPeptide', right_on='Modified.Sequence', how='inner')
df_filtered = merged_df[merged_df['PrecursorCharge'] >= 2]

######################## Generate four Plots  ############################
bufs = []
for i in range(4):
    buf = io.BytesIO()
    plt.figure()

    if i == 0:  # Plot 1: RetentionTime vs M/Z (All Features)
        plt.scatter(RetentionTime_values, MZ_values, s=0.1, color='blue')
        plt.xlabel('RetentionTime')
        plt.ylabel('M/Z')
        plt.title('RetentionTime and M/Z (All Features)')
    elif i == 1:  # Plot 2: RetentionTime vs M/Z (Charge >= 2)
        plt.scatter(RetentionTime_values_high_charge, MZ_values_high_charge, s=0.1, color='blue')
        plt.xlabel('Retention Time')
        plt.ylabel('M/Z')
        plt.title('Retention Time and M/Z (Charge >= 2)')
    elif i == 2:  # Plot 3: All Charges from merged data
        plt.scatter(merged_df['RT'], merged_df['PrecursorMz'], alpha=0.5, s=0.5)
        plt.title('All Charges (DiaNN)')
        plt.xlabel('Retention Time (RT)')
        plt.ylabel('Precursor M/Z')
        plt.grid(True)
    elif i == 3:  # Plot 4: Charge ≥ 2 from filtered data
        plt.scatter(df_filtered['RT'], df_filtered['PrecursorMz'], alpha=0.5, s=0.5)
        plt.title('Charge ≥ 2 (DiaNN) ')
        plt.xlabel('Retention Time (RT)')
        plt.ylabel('Precursor M/Z')
        plt.grid(True)

    plt.savefig(buf, format='png')
    buf.seek(0)
    bufs.append(buf)
    plt.close()

################### Add Text Content to PDF ###################
report_content = f"""
Number of distinct isolation windows used: {len(unique_targets)}
Window Size: {average_upper_bound - average_lower_bound:.2f}

Average Injection Time: {average_injection_time:.2f} ms
Median Injection Time: {median_injection_time:.2f} ms
"""
# Add the new metrics from the second file
report_content += f"""
Precursors Identified: {precursors_identified}
"""
report_content += f"""
Scan number MS1: {scan_number_MS1}
Scan number MS/MS: {scan_number_MSMS}
"""
# Add report content to PDF
# Define starting Y position for text
current_y = height - 72  # Starting Y position at the top of the page
# Set a fixed amount to decrement the Y position for each line
line_height = 14
text_lines = report_content.strip().split('\n')
for line in text_lines:
    c.drawString(72, current_y, line)
    current_y -= line_height

current_y -= 30  # Additional space before starting images

# Assuming bufs is a list of BytesIO objects for your images
image_positions = []  # To store calculated positions and sizes for each image
for i, buf in enumerate(bufs):
    image = ImageReader(buf)
    aspect_ratio = image.getSize()[0] / float(image.getSize()[1])
    
    # Calculate image width and height based on desired width or max height
    image_width = (width - 144) / 2  # Two images per row, with margins
    image_height = image_width / aspect_ratio
    
    # Calculate X, Y positions for the current image
    x = 72 + (i % 2) * (width / 2)
    y = current_y - (i // 2 + 1) * image_height - 10 * (i // 2)  # Adjust Y based on row
    
    image_positions.append((image, x, y, image_width, image_height))
    
    # Update current_y for the first image only, to avoid reducing it for every image
    if i == 1 or i == 3:  # After drawing two images in a row
        current_y = y + 150  # Adjust space for next row of images

# Draw images on the PDF using calculated positions and sizes
for image, x, y, image_width, image_height in image_positions:
    c.drawImage(image, x, y, width=image_width, height=image_height)

# Save the PDF
c.save()
