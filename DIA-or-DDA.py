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
    else:
        print("DDA file.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script_name.py path_to_your_file.csv")
        sys.exit(1)
    csvfile = sys.argv[1]
    analyze_isolation_window_target(csvfile)
