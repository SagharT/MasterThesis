'''
This script is designed to parse an mzML file, to extract specific features from the data. 
The pyteomics library,is used for parsing the mzML file. 
The script is extracting: precursor ion m/z, retention time, injection time, 
and details about the isolation window (target, upper offset, and lower offset). 
The extracted data is structured into a list of dictionaries.
Finally, this data is written to a CSV file.
Before running this script, ensure that 'pyteomics' and 'lxml' libraries are installed. 
These can be installed via pip.
'''
import csv
import sys
#pip install pyteomics
#pip install lxml (if import pyteomics then we need to write pyteomics.mzml)
from pyteomics import mzml

# Function to parse mzML and extract desired features
def parse_mzml(file_path):
    results = []

    with mzml.read(file_path) as reader:
        for spectrum in reader:
            if 'precursorList' in spectrum and spectrum['precursorList']['count'] > 0:
                precursor = spectrum['precursorList']['precursor'][0]
                selected_ion = precursor['selectedIonList']['selectedIon'][0]

                precursor_mz = selected_ion['selected ion m/z']
                retention_time = spectrum['scanList']['scan'][0]['scan start time']
                injection_time = spectrum['scanList']['scan'][0].get('ion injection time', 'N/A')
                isolation_window_target = precursor['isolationWindow'].get('isolation window target m/z', 'N/A')
                isolation_window_upper = precursor['isolationWindow'].get('isolation window upper offset', 'N/A')
                isolation_window_lower = precursor['isolationWindow'].get('isolation window lower offset', 'N/A')
#create a dictionary and add it to result
                results.append({
                    'Precursor Ion m/z': precursor_mz,
                    'Retention Time': retention_time,
                    'Injection Time': injection_time,
                    'Isolation Window Target': isolation_window_target,
                    'Isolation window upper offset':isolation_window_upper,
                    'Isolation window lower offset': isolation_window_lower
                })

    return results

# Write results to a CSV file
def write_to_csv(results, output_file):
    keys = results[0].keys()
    with open(output_file, 'w', newline='') as file:
        dict_writer = csv.DictWriter(file, keys)
        dict_writer.writeheader()
        dict_writer.writerows(results)

# Main execution
mzml_file = sys.argv[1]
output_csv = sys.argv[2]
# call the two functions
parsed_data = parse_mzml(mzml_file)
write_to_csv(parsed_data, output_csv)

print("Parsing completed:", output_csv)
