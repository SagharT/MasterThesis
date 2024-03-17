#python3 Compare.py 20240220/221020/MzML/*_features.csv 20240220/221020/DiaNN/*.stats.tsv 20240220/221020/*.mzidsummary.txt AnalysisReport.pdf

import os
import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

# Common setup
script_dir = os.path.dirname(os.path.abspath(__file__))
plot_filename = os.path.join(script_dir, 'plot1.png')
plot2_filename = os.path.join(script_dir, 'plot2.png')
Report = sys.argv[-1]

scannumbers = {}
precursors_identified_counts = {}
ng_regex = re.compile(r"(\d[\d_]*\d*\.?\d*)ng")

# File grouping with added flexibility for both file types
group1_files = [f for f in sys.argv[1:] if '_features.csv' in f]
group2_files = [f for f in sys.argv[1:] if '.stats.tsv' in f or '.mzidsummary.txt' in f]

######################## Process Group 1 files (CSV)
for inputfile in group1_files:
    scan_number_MSMS = 0
    scan_number_MS1 = 0

    with open(inputfile) as file:
        print(inputfile)
        csv_reader = csv.DictReader(file)  # Use DictReader to access columns by header
        for row in csv_reader:
            # Check and count non-'N/A' values for 'Precursor Ion m/z'
            if row.get('Precursor Ion m/z') and row['Precursor Ion m/z'] != 'N/A':
                scan_number_MSMS += 1
            # Check and count non-'N/A' values for 'MS1 Base Peak m/z'
            if row.get('MS1 Base Peak m/z') and row['MS1 Base Peak m/z'] != 'N/A':
                scan_number_MS1 += 1
        ng_match = ng_regex.search(inputfile)
        ng = ng_match.group(0)

        if '_' in ng:
            formatted_ng = ng.replace('_', ',')
        elif ng.startswith("0"):
            formatted_ng = '0,' + ng[1:]
        else:
            formatted_ng = ng

        # Use formatted_ng for dictionary operations
        if formatted_ng not in scannumbers:
            scannumbers[formatted_ng] = ([], [])
        scannumbers[formatted_ng][0].append(scan_number_MSMS)
        scannumbers[formatted_ng][1].append(scan_number_MS1)

# Print each key and value in the desired format
for key, (msms, ms1) in scannumbers.items():
    msms = np.average(msms)
    ms1 = np.average(ms1)
    print(f"{key}: {msms} {ms1}")

# Prepare data for plotting
formatted_ng_labels = list(scannumbers.keys())


# Define a sorting key function
def ng_value_sort_key(ng_value):
    # Replace commas with dots and remove 'ng' suffix, then convert to float
    numeric_part = ng_value.replace(',', '.').replace('ng', '')
    try:
        return float(numeric_part)
    except ValueError:
        # In case of a conversion error, return zero or handle it as appropriate
        return 0
    
# Assuming ng_value_sort_key function is already defined
# Sort NG values for consistent order in plotting
width = 0.35
sorted_ng_labels = sorted(scannumbers.keys(), key=ng_value_sort_key)
sorted_msms_counts = [np.mean(scannumbers[ng][0]) for ng in sorted_ng_labels]
sorted_ms1_counts = [np.mean(scannumbers[ng][1]) for ng in sorted_ng_labels]

x = np.arange(len(sorted_ng_labels))  # the label locations, now sorted

fig, ax = plt.subplots()
bars1 = ax.bar(x - width/2, sorted_ms1_counts, width, label='MS1', color='yellow')
bars2 = ax.bar(x + width/2, sorted_msms_counts, width, label='MS2', color='blue')

# Update text for labels, title, and custom x-axis tick labels to use sorted NG values
ax.set_xlabel('NG Value')
ax.set_ylabel('Average Count')
ax.set_title('Average MS1 and MSMS Spectra by Concentration')
ax.set_xticks(x)
ax.set_xticklabels(sorted_ng_labels, rotation=45, ha='right')  # Use sorted labels
ax.legend()


plt.tight_layout()
fig.savefig(plot_filename, bbox_inches='tight')


######################### Process Group 2 files (.stats.tsv or .mzidsummary.txt)
for inputfile in group2_files:
    precursors_identified = 0  # Default to 0
    formatted_ng = ''

    # Check file type and process accordingly
    if inputfile.endswith('.stats.tsv'):
        with open(inputfile, newline='') as file:
            tsv_reader = csv.DictReader(file, delimiter='\t')
            for row in tsv_reader:
                try:
                    precursors_identified += int(row['Precursors.Identified'])
                except ValueError:
                    pass  # Handle non-integer values gracefully
    elif inputfile.endswith('.mzidsummary.txt'):
        with open(inputfile, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 2:
                precursors_identified = int(lines[1].split()[0])
    
    # Common code for extracting and formatting ng value from filename
    ng_match = ng_regex.search(inputfile)
    if ng_match:
        ng = ng_match.group(0)
        if '_' in ng:
            formatted_ng = ng.replace('_', ',')
        elif ng.startswith("0"):
            formatted_ng = '0,' + ng[1:]
        else:
            formatted_ng = ng

    # Update the dictionary for precursors identified counts
    if formatted_ng not in precursors_identified_counts:
        precursors_identified_counts[formatted_ng] = []
    precursors_identified_counts[formatted_ng].append(precursors_identified)

# Common code for plotting and PDF generation
# Print each key and value in the desired format for the second group
for key, precursors_counts in precursors_identified_counts.items():
    avg_precursors_identified = np.average(precursors_counts)
    print(f"{key}: {avg_precursors_identified}")


formatted_ng_labels_group2 = list(precursors_identified_counts.keys())

all_formatted_ng_labels = sorted(set(formatted_ng_labels + formatted_ng_labels_group2))

avg_precursors_counts = [np.mean(counts) for counts in precursors_identified_counts.values()]

# Define a sorting key function
def ng_value_sort_key(ng_value):
    # Replace commas with dots and remove 'ng' suffix, then convert to float
    numeric_part = ng_value.replace(',', '.').replace('ng', '')
    try:
        return float(numeric_part)
    except ValueError:
        # In case of a conversion error, return zero or handle it as appropriate
        return 0
    
# Sort NG values and ensure consistent order for plotting
sorted_ng_keys = sorted(precursors_identified_counts.keys(), key=ng_value_sort_key)
sorted_avg_precursors_counts = [np.mean(precursors_identified_counts[ng]) for ng in sorted_ng_keys]


x = np.arange(len(sorted_ng_keys))  # the label locations
width = 0.35  # the width of the bars

fig2, ax = plt.subplots()
bars1 = ax.bar(x - width/2, sorted_avg_precursors_counts, width, label='precursors counts', color='green')

# Add some text for labels, title, and custom x-axis tick labels, etc.
ax.set_xlabel('NG Value')
ax.set_ylabel('Precursors Count')
ax.set_title('Precursors Counts by Concentration')
ax.set_xticks(x)
ax.set_xticklabels(sorted_ng_keys, rotation=45, ha='right')
ax.legend()


plt.tight_layout()

fig.savefig(plot_filename, bbox_inches='tight')
fig2.savefig(plot2_filename, bbox_inches='tight')


#Create the Report
c = canvas.Canvas(Report, pagesize=letter)
width, height = letter 
# Load the plot image
plot_image = ImageReader(plot_filename)
plot2_filename = ImageReader(plot2_filename)

# Draw the images onto the PDF
c.drawImage(plot_image, x=72, y=height-8*72, width=6*72, preserveAspectRatio=True, anchor='c')
c.showPage()  # Finish the first page
c.drawImage(plot2_filename, x=72, y=height-8*72, width=6*72, preserveAspectRatio=True, anchor='c')

# Save the PDF
c.save()
