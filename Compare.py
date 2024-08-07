#python3 Compare.py /MzML/*_features.csv /DiaNN/*.stats.tsv /*.mzidsummary.txt AnalysisReport.pdf
import pandas as pd
import os
import sys
import csv
import re
import numpy as np
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4, landscape
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

# Common setup
script_dir = os.path.dirname(os.path.abspath(__file__))
plot_filename = os.path.join(script_dir, 'plot1.png')
plot2_filename = os.path.join(script_dir, 'plot2.png')
Report = sys.argv[-1]

scannumbers = {}
precursors_identified_counts = {}
ng_regex = re.compile(r"\/([^\/]*?)(\d+_)?([\d]+ng)(.*?)(_\d)?\.")

# File grouping with added flexibility for both file types
group1_files = [f for f in sys.argv[1:] if '.features.csv' in f]
group2_files = [f for f in sys.argv[1:] if '.stats.tsv' in f or '.mzidsummary.txt' in f]

#################### Find the type of inputfile DDA or DIA
def analyze_isolation_window_target(file):
    # Load the CSV file
    df = pd.read_csv(file)
    
    # Select the "Isolation Window Target" column, replace 'N/A' with np.nan, and convert to float
    column = df['Isolation Window Target'].replace('N/A', np.nan).astype(float)
     
    # Calculate the number of unique numeric values
    unique_numeric_values = column.dropna().unique()
    
    # It is DIA if the pattern repeated at least 50 times
    if len(unique_numeric_values) < len(column.dropna()) / 50:
        return "DIA"
    else:
        return "DDA"

def create_keyname(file, no_match_name, filetype):
    ng_match = ng_regex.search(file)
    if ng_match is None:
        print("No match", file)
        keyname = no_match_name
    else:
        print(ng_match.group(1), ng_match.group(2), ng_match.group(3), ng_match.group(4))
        ng_part1 = ng_match.group(2)
        ng_part2 = ng_match.group(3)

        if ng_part1 is not None and 1 < len(ng_part1) <= 4:
            formatted_ng = ng_part1.replace('_', ',') + ng_part2
        elif ng_part2.startswith("0"):
            formatted_ng = '0,' + ng_part2[1:]
        elif ng_part1 is not None:
            formatted_ng = ng_part1 + "_" + ng_part2
        else:
            formatted_ng = ng_part2

        keyname = ng_match.group(1) + formatted_ng + ng_match.group(4)

    # Remove 'DIA' or 'DDA' from keyname (if not part of the intentional filetype addition)
    keyname = re.sub(r'(dia|dda)', '', keyname, flags=re.IGNORECASE)

    if len(unique_filetypes) > 1:
        keyname += " " + filetype
    return keyname

unique_filetypes = set([analyze_isolation_window_target(inputfile) for inputfile in group1_files])

######################## Process Group 1 files (CSV)
for inputfile in group1_files:
    filetype = analyze_isolation_window_target(inputfile)
    scan_number_MSMS = 0
    scan_number_MS1 = 0

    with open(inputfile) as file:
        csv_reader = csv.DictReader(file)  # Use DictReader to access columns by header
        for row in csv_reader:
            # Check and count non-'N/A' values for 'Precursor Ion m/z'
            if row.get('Precursor Ion m/z') and row['Precursor Ion m/z'] != 'N/A':
                scan_number_MSMS += 1
            # Check and count non-'N/A' values for 'MS1 Base Peak m/z'
            if row.get('MS1 Base Peak m/z') and row['MS1 Base Peak m/z'] != 'N/A':
                scan_number_MS1 += 1

        keyname = create_keyname(inputfile, os.path.basename(inputfile).replace(".features.csv", ""), filetype)

        # Use keyname for dictionary operations
        if keyname not in scannumbers:
            scannumbers[keyname] = ([], [])
        scannumbers[keyname][0].append(scan_number_MSMS)
        scannumbers[keyname][1].append(scan_number_MS1)

# Print each key and value in the desired format
for key, (msms, ms1) in scannumbers.items():
    msms = np.average(msms)
    ms1 = np.average(ms1)
    print(f"{key}: {msms} {ms1}")

# Prepare data for plotting
formatted_ng_labels = list(scannumbers.keys())


# Define a sorting key function
def ng_value_sort_key(ng_item):
    ng_key, ng_value = ng_item
    return np.mean(ng_value[0])
    
# Assuming ng_value_sort_key function is already defined
# Sort NG values for consistent order in plotting
width = 0.35
sorted_ng_items = sorted(scannumbers.items(), key=ng_value_sort_key)
sorted_ng_labels = [item[0] for item in sorted_ng_items]
sorted_msms_counts = [np.mean(item[1][0]) for item in sorted_ng_items]
sorted_ms1_counts = [np.mean(item[1][1]) for item in sorted_ng_items]

x = np.arange(len(sorted_ng_labels))  # the label locations, now sorted

fig, ax = plt.subplots()
bars1 = ax.bar(x - width/2, sorted_ms1_counts, width, label='MS1', color='yellow')
bars2 = ax.bar(x + width/2, sorted_msms_counts, width, label='MS2', color='blue')

# Update text for labels, title, and custom x-axis tick labels to use sorted NG values
ax.set_xlabel('Samples')
ax.set_ylabel('Average Count')
ax.set_title('Average MS1 and MSMS Spectra for each sample')
ax.set_xticks(x)
ax.set_xticklabels(sorted_ng_labels, rotation=45, ha='right')  # Use sorted labels
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))


plt.tight_layout()
fig.savefig(plot_filename, bbox_inches='tight')


######################### Process Group 2 files (.stats.tsv or .mzidsummary.txt)
for inputfile in group2_files:
    precursors_identified = -1  # Default to -1
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
        filetype = 'DIA'
    elif inputfile.endswith('.mzidsummary.txt'):
        with open(inputfile, 'r') as file:
            lines = file.readlines()
            if len(lines) >= 2:
                precursors_identified = int(lines[1].split()[0])
        filetype = 'DDA'
    
    keyname = create_keyname(inputfile, os.path.basename(inputfile).replace('.stats.tsv', "").replace('.mzidsummary.txt', ""), filetype)

    # Update the dictionary for precursors identified counts
    if precursors_identified != -1:
        if keyname not in precursors_identified_counts:
            precursors_identified_counts[keyname] = []
        precursors_identified_counts[keyname].append(precursors_identified)

# Print each key and value in the desired format for the second group
for key, precursors_counts in precursors_identified_counts.items():
    avg_precursors_identified = np.average(precursors_counts)
    print(f"{key}: {avg_precursors_identified}")


formatted_ng_labels_group2 = list(precursors_identified_counts.keys())

all_formatted_ng_labels = sorted(set(formatted_ng_labels + formatted_ng_labels_group2))

avg_precursors_counts = [np.mean(counts) for counts in precursors_identified_counts.values()]


# Define a sorting key function
def ng_value_sort_key(ng_item):
    ng_key, ng_value = ng_item
    return np.mean(ng_value[0])
    
# Sort NG values and ensure consistent order for plotting
sorted_ng_items = sorted(precursors_identified_counts.items(), key=ng_value_sort_key)
sorted_ng_keys = [item[0] for item in sorted_ng_items]
sorted_avg_precursors_counts = [np.mean(item[1]) for item in sorted_ng_items]


x = np.arange(len(sorted_ng_keys))  # the label locations
width = 0.35  # the width of the bars

fig2, ax = plt.subplots()
bars1 = ax.bar(x - width/2, sorted_avg_precursors_counts, width, label='precursor\ncounts', color='green')

# Add some text for labels, title, and custom x-axis tick labels, etc.
ax.set_xlabel('Samples')
ax.set_ylabel('Identified Precursors')
ax.set_title('Identified Precursors comparison between samples')
ax.set_xticks(x)
ax.set_xticklabels(sorted_ng_keys, rotation=45, ha='right')
ax.legend(loc='upper left', bbox_to_anchor=(1, 1))


plt.tight_layout()

fig.savefig(plot_filename, bbox_inches='tight')
fig2.savefig(plot2_filename, bbox_inches='tight')


#Create the Report
c = canvas.Canvas(Report, pagesize=landscape(A4))
width, height = landscape(A4) 
# Load the plot image
plot_image = ImageReader(plot_filename)
plot2_image = ImageReader(plot2_filename)

# Draw the images onto the PDF
c.drawImage(plot_image, x=0, y=0, width=width, preserveAspectRatio=True, anchor='c')
c.showPage()  # Finish the first page
c.drawImage(plot2_image, x=0, y=0, width=width, preserveAspectRatio=True, anchor='c')

# Save the PDF
c.save()
