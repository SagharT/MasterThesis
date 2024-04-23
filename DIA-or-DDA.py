import pandas as pd
import numpy as np
import sys

def analyze_isolation_window_target(file_path):
    # Load the CSV file
    df = pd.read_csv(file_path)
    
    # Select the "Isolation Window Target" column, replace 'N/A' with np.nan, and convert to float
    column = df['Isolation Window Target'].replace('N/A', np.nan).astype(float)
    
    # Check if there are enough numeric values
    if column.isna().all():
        print("The column contains only 'N/A' values or is empty.")
        return
    
    # Calculate the number of unique numeric values
    unique_numeric_values = column.dropna().unique()
    
    # It is DIA if the pattern repeated at least 50 times
    if len(unique_numeric_values) < len(column.dropna()) / 50:
        print("DIA file.")
        return "DIA file."
    else:
        print("DDA file.")
        return "DDA file."


if __name__ == "__main__":
    files = sys.argv[1:]  # Get all provided file paths
    dia_files_count = 0

    for file_path in files:
        result = analyze_isolation_window_target(file_path)
        if result == "DIA file.":
            dia_files_count += 1

    # Check if all provided files are DIA
    if dia_files_count == len(files):
        print("DIA file.")
        sys.exit(0)  # Exit with status 0 to indicate all files are DIA
    else:
        print("DDA file.")
        sys.exit(1)
    
