# python3 DDA-Report.py /home/ubuntu/snake-test/20240226-DDAtest/MzML/231010_50ngHela31min-1_features.csv /home/ubuntu/snake-test/20240226-DDAtest/231010_50ngHela31min-1.features.tsv /home/ubuntu/snake-test/20240226-DDAtest/231010_50ngHela31min-1.mzid .pdf
# This script is Features.py + MZid-parser.py
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader
import xml.etree.ElementTree as ET
import io
import csv
import sys
import os
import glob

# Input file 
FeaturesFile = sys.argv[1]
AdditionalDataFile = sys.argv[2]
XmlFile = sys.argv[3]
Report = sys.argv[4]

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
# Parse the XML file
tree = ET.parse(XmlFile)
root = tree.getroot()

# Namespace to access elements correctly
namespaces = {'ns': 'http://psidev.info/psi/pi/mzIdentML/1.1'}

# Initialize lists for retention times, charges, and m/z values
retention_times = []
charges = []
mz_values = []

# Iterate through each SpectrumIdentificationResult
for spectrum_result in root.findall('.//ns:SpectrumIdentificationResult', namespaces):
    # Find and store the retention time for the current SpectrumIdentificationResult
    retention_time = None
    for cv_param in spectrum_result.findall('.//ns:cvParam', namespaces):
        if cv_param.get('name') == 'scan start time':
            retention_time = float(cv_param.get('value'))
            break  # Stop looking for retention time once found

    # If retention time is found, proceed to check each SpectrumIdentificationItem in the result
    if retention_time is not None:
        for spectrum_item in spectrum_result.findall('.//ns:SpectrumIdentificationItem', namespaces):
            q_value = None
            for cv_param in spectrum_item.findall('.//ns:cvParam', namespaces):
                if cv_param.get('accession') == 'MS:1002054':  # MS-GF:QValue
                    q_value = float(cv_param.get('value'))
                    break  # Stop once the QValue is found
            
            # If QValue meets the condition, collect the necessary information
            if q_value is not None and q_value < 0.01:  # Assuming the condition is QValue < 0.01
                charges.append(int(spectrum_item.get('chargeState')))
                mz_values.append(float(spectrum_item.get('experimentalMassToCharge')))
                retention_times.append(retention_time)

precursors_identified = len(retention_times)
filtered_retention_times = [retention_times[i] for i in range(len(charges)) if charges[i] >= 2]
filtered_mz_values = [mz_values[i] for i in range(len(charges)) if charges[i] >= 2]
################### Generate Plots and Save to BytesIO Buffers #################
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
        plt.scatter(retention_times, mz_values, alpha=0.5, s=0.5)
        plt.title('All Charges (DiaNN)')
        plt.xlabel('Retention Time (RT)')
        plt.ylabel('Precursor M/Z')
        plt.grid(True)
    elif i == 3:  # Plot 4: Charge ≥ 2 from filtered data
        plt.scatter(filtered_retention_times, filtered_mz_values, alpha=0.5, s=0.5)
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
