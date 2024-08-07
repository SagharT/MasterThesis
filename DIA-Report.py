# python3 Report.py /MzML/_features.csv .stats.tsv .features.tsv  _result.tsv .tsv .pdf

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
TableFile = sys.argv[7]

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

window_sizes_by_target = {}
scan_number_MSMS = 0
scan_number_MS1 = 0
injection_times = []
demultiplexing_applied = False
target_details = {}

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

        if all(key in row for key in ['Isolation Window Target', 'Isolation window upper offset', 'Isolation window lower offset']) and all(row[key] != 'N/A' for key in ['Isolation Window Target', 'Isolation window upper offset', 'Isolation window lower offset']):
            target = float(row['Isolation Window Target'])
            upper_offset = float(row['Isolation window upper offset'])
            lower_offset = float(row['Isolation window lower offset'])
            window_size = round(upper_offset + lower_offset, 1)

            # Check for demultiplexing signal in the data and adjust window size
            if row.get('Demultiplexing') == 'Yes':  # This needs to be based on your specific data marking for demultiplexing
                demultiplexing_applied = True
                window_size *= 2  # Double the size if demultiplexed
                # Store details for each target
            target_details[target] = {
                'upper_offset': upper_offset,
                'lower_offset': lower_offset,
                'window_size': window_size
            }
                
            window_sizes_by_target[target] = window_size  # Store window size for each target

window_size_counts = {}
for size in window_sizes_by_target.values():
    if size not in window_size_counts:
        window_size_counts[size] = 0
    window_size_counts[size] += 1
    
num_isolation_windows = len(window_sizes_by_target)
if demultiplexing_applied:
    # Adjust the count of each window size based on the new total number of targets
    adjustment_factor = (num_isolation_windows - 2) // 2
    adjusted_window_size_counts = {}
    for size, count in window_size_counts.items():
        # Reduce the original count proportionally based on the adjustment factor
        adjusted_count = max(1, count * adjustment_factor // num_isolation_windows)
        adjusted_window_size_counts[size] = adjusted_count
    window_size_counts = adjusted_window_size_counts  # Use the adjusted counts for summary
    num_isolation_windows = adjustment_factor  # Update the total number of isolation windows

# Generate the summary of window sizes with corrected target counts
window_size_summary = "; ".join(f"{size} for {count} windows" for size, count in window_size_counts.items())


# Calculate average and median
average_injection_time = np.mean(injection_times)
median_injection_time = np.median(injection_times)


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
# Set the default font size for all plot elements
plt.rcParams.update({'font.size': 12})
bufs = []
for i in range(4):
    buf = io.BytesIO()
    plt.figure()

    if i == 0:  # Plot 1: RetentionTime vs m/z (All Features)
        plt.scatter(RetentionTime_values, MZ_values, s=0.1, color='blue')
        plt.xlabel('RetentionTime')
        plt.ylabel('m/z')
        plt.title('RetentionTime and m/z (All Features)')
    elif i == 1:  # Plot 2: RetentionTime vs m/z (Charge >= 2)
        plt.scatter(RetentionTime_values_high_charge, MZ_values_high_charge, s=0.1, color='blue')
        plt.xlabel('Retention Time')
        plt.ylabel('m/z')
        plt.title('Retention Time and m/z (Charge >= 2)')
    elif i == 2:  # Plot 3: All Charges from merged data
        plt.scatter(merged_df['RT'], merged_df['PrecursorMz'], alpha=0.5, s=0.5)
        plt.title('All Charges (Identified features, Dia-NN)')
        plt.xlabel('Retention Time (RT)')
        plt.ylabel('Precursor m/z')
        plt.grid(True)
    elif i == 3:  # Plot 4: Charge ≥ 2 from filtered data
        plt.scatter(df_filtered['RT'], df_filtered['PrecursorMz'], alpha=0.5, s=0.5)
        plt.title('Charge ≥ 2 (Identified features, Dia-NN) ')
        plt.xlabel('Retention Time (RT)')
        plt.ylabel('Precursor m/z')
        plt.grid(True)

    plt.savefig(buf, format='png')
    buf.seek(0)
    bufs.append(buf)
    plt.close()

################### Add Text Content to PDF ###################
demultiplexing_status = "DIA sample using overlapping windows and demultiplexing after" if demultiplexing_applied else "DIA sample"
report_content = f"""
{demultiplexing_status}.
Number of distinct isolation windows used: {num_isolation_windows}
Detailed Window Size: {window_size_summary}

Average Injection Time: {average_injection_time:.2f} ms
Median Injection Time: {median_injection_time:.2f} ms
"""
# Add the new metrics from the second file
report_content += f"""
Precursors Identified: {precursors_identified}
"""
report_content += f"""
Number of MS1 scans: {scan_number_MS1}
Number of MS/MS scans: {scan_number_MSMS}
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


# Prepare the header and data rows for the table
data = [['Target', 'Upper Window', 'Lower Window', 'Window Size']]
for target, details in target_details.items():
    data.append([
        target,
        details['upper_offset'],
        details['lower_offset'],
        details['window_size']
    ])

# Write data to CSV file
with open(TableFile, 'w', newline='') as csvfile:
    csv_writer = csv.writer(csvfile)
    
    # Write rows
    csv_writer.writerows(data)