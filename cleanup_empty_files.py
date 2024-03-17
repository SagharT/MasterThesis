# cleanup_empty_files.py
import os
import sys

def remove_empty_files(directory):
    for root, _, files in os.walk(directory):
        for file in files:
            file_path = os.path.join(root, file)
            if os.path.getsize(file_path) == 0:
                os.remove(file_path)
                print(f"Removed empty file: {file_path}")

if __name__ == "__main__":
    for dir_path in sys.argv[1:]:
        remove_empty_files(dir_path)
