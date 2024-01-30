'''
open different files stats,.tsv with different ng 
extract number of FWHM from 7th column: FWHM.Scans 
extract the number before ng from the first column: File.Name
for example 50 here: ./MzML/221025_HeLa_50ng_DIA.mzML
221025_HeLa_50ng_DIA.stats.tsv
file name will be the key and ng and FWHM values
make a graph for FWHM for each ng
'''
import glob
import sys
import os
import csv
import re
import matplotlib.pyplot as plt

directory = sys.argv[1]
os.chdir(directory)
stats_files = glob.glob('*.stats.tsv')
data = {}
for i in stats_files:
    with open(i ,"r", encoding="utf8") as file:
        tsv_reader = csv.reader(file, delimiter="\t")
        # Skip the first row, which is the header
        next(tsv_reader)
        row = next(tsv_reader)
        FWHM = float(row[6])
        match = re.search(r'_(\d+)ng', row[0])
        ng = int(match.group(1))
        data[ng] = FWHM

# Plotting the data
ngs = list(data.keys())
ngs.sort()
plt.plot(ngs, [data[ng] for ng in ngs], '-o')
for ng in ngs:
	plt.text(ng, data[ng], f'{ng}', fontsize=10)
plt.xlabel('ng')
plt.ylabel('FWHM')
plt.title('FWHM for each ng')
plt.legend()
plt.show()


