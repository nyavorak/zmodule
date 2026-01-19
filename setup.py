import subprocess
import os

folders = ["pyGRBaglow/", "pyGRBz/"]

def install_folder(folder):
    """Install pyGRBaglow and required modules
       Install pyGRBz and required modules
    """    
    try:
        print(f"Installing dependencies in {folder}...")
        subprocess.check_call(["pip", "install", "-e", "."], cwd=folder)
        print(f"Successfully installed {folder}.")
    except subprocess.CalledProcessError as e:
        print(f"Failed to install {folder}. Error: {e}")
        exit(1)

if __name__ == "__main__":
    # Iterate over each folder and install
    for folder in folders:
        if os.path.isdir(folder):
            install_folder(folder)
        else:
            print(f"Folder {folder} does not exist. Skipping...")
