import xml.etree.ElementTree as ET
import sys
import matplotlib.pyplot as plt

file_path = sys.argv[1]

# Parse the XML file
#The parse function reads the XML document from the file and creates a tree structure that represents the XML document.
tree = ET.parse(file_path)
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


# Graph 1: Retention Time vs. M/Z Values for All Items
plt.figure(figsize=(10, 6))
plt.scatter(retention_times, mz_values, alpha=0.5, s=0.5)
plt.title('Retention Time vs. M/Z Values for All Items')
plt.xlabel('Retention Time (min)')
plt.ylabel('M/Z')
plt.grid(True)
plt.show()

# Preparing data for Graph 2: Filter for charge at least 2
filtered_retention_times = [retention_times[i] for i in range(len(charges)) if charges[i] >= 2]
filtered_mz_values = [mz_values[i] for i in range(len(charges)) if charges[i] >= 2]

# Graph 2: Retention Time vs. M/Z Values for Items with Charge >= 2
plt.figure(figsize=(10, 6))
plt.scatter(filtered_retention_times, filtered_mz_values, alpha=0.5, s=0.5)
plt.title('Retention Time vs. M/Z Values for Items with Charge >= 2')
plt.xlabel('Retention Time (min)')
plt.ylabel('M/Z')
plt.grid(True)
plt.show()
