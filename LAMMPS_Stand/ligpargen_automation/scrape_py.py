import os
import time
import shutil
from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.by import By

# === List of SMILES ===
import pandas as pd

# Replace with the path to your Excel file
df = pd.read_excel('Electrolytes.xlsx')

# Display the first few rows
smiles_list = [smi for i, smi in df['rdkit_smiles'].items()]

print(len(smiles_list), "SMILES found in the file.")

# === Paths ===
download_dir = os.getcwd()  # Current working directory for downloads
output_dir = os.path.join(os.getcwd(), "output_lmp")
chromedriver_path = os.path.join(os.getcwd(), "chromedriver")  # Ensure chromedriver is in the current directory

# === Create output_lmp directory if it doesn't exist ===
os.makedirs(output_dir, exist_ok=True)

# === Chrome Driver Setup ===
service = Service(chromedriver_path)
options = webdriver.ChromeOptions()
prefs = {
    "download.default_directory": download_dir,
    "download.prompt_for_download": False,
    "safebrowsing.enabled": True
}
options.add_experimental_option("prefs", prefs)
driver = webdriver.Chrome(service=service, options=options)

for smiles in smiles_list:
    try:
        print(f"Processing: {smiles}")
        # === Step 1: Open LigParGen ===
        driver.get("https://traken.chem.yale.edu/ligpargen/")

        # === Step 2: Enter SMILES ===
        text_box = driver.find_element(By.ID, "smiles")
        text_box.clear()
        text_box.send_keys(smiles)

        # === Step 3: Click Submit Molecule ===
        submit_button = driver.find_element(By.XPATH, "//button[text()='Submit Molecule']")
        submit_button.click()

        # === Step 4: Wait for new page to load ===
        time.sleep(5)

        # === Step 5: Click LAMMPS ===
        lammps_button = driver.find_element(By.XPATH, "//input[@type='submit' and @value='LAMMPS']")
        lammps_button.click()

        # === Step 6: Wait for file download to complete ===
        downloaded_file = None
        timeout = 30  # seconds
        for _ in range(timeout):
            files = [f for f in os.listdir(download_dir) if f.endswith(".lmp")]
            if files:
                downloaded_file = files[0]
                break
            time.sleep(1)

        # === Step 7: Move and Rename File ===
        if downloaded_file:
            src = os.path.join(download_dir, downloaded_file)
            dst = os.path.join(output_dir, f"{smiles.replace('/', '_')}.lmp")
            shutil.move(src, dst)
            print(f"Downloaded .lmp saved as: {dst}")
        else:
            print(f"Download failed for {smiles}.")

    except Exception as e:
        print(f"Error processing {smiles}: {e}")

# === Step 8: Close Browser ===
driver.quit()