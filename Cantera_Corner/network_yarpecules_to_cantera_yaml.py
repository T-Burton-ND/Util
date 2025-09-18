import cantera as ct
import pandas as pd
import numpy as np
import yarp as yp
from rdkit.Chem import MolFromSmiles as smi2mol
from rdkit.Chem import MolToSmiles as mol2smi
import io

def get_elements(smi):
    temp_dict = {}
    elements_list = yp.yarpecule(smi).elements
    for i in elements_list:
        if i.capitalize() not in temp_dict:
            temp_dict[i.capitalize()] = 1
        else:
            temp_dict[i.capitalize()] += 1
    return temp_dict

# filtered_h2o_min_dg.csv

def get_rate_constant(idx,rsmi_list,Pressure):
    # uses the erying equation, does not automatically pull the imaginary frequency for HAT correction (once that is not recorded by YARP it will be possible)
    molecularity = len(rsmi_list)
    K_HAT = 1
    kB = 1.380649*10**-23
    h = 6.62607015*10**-34
    R = 0.000082057366080960
    return (K_HAT*kB/h*(R/Pressure)**(molecularity-1)),(molecularity)

def build_yaml(yaml_name):
    df = pd.read_csv(yaml_name,header=0)
    Pressure = 1

    all_strucs = []
    initial_strucs = ['O',mol2smi(smi2mol('OCC1OC(O)C(O)C(O)C1O'))]

    for idx,row in df.iterrows():
        Rsmi = row['Rsmi_canon'].split('.')
        for i in Rsmi:
            x = mol2smi(smi2mol(i))
            if x not in all_strucs:
                all_strucs.append(i)
        Psmi = row['Psmi_canon'].split('.')
        for i in Psmi:
            x = mol2smi(smi2mol(i))
            if x not in all_strucs:
                all_strucs.append(i)
                
    f = io.StringIO()
    f.write("units: {time: s, quantity: mol, activation-energy: kcal/mol, pressure: atm, energy: kcal")
    f.write("}\n\n")
    f.write("phases:\n")
    f.write("- name: Sim\n")
    f.write("  elements: ['C','O','H']\n")
    f.write("  thermo: ideal-gas\n")
    f.write("  species: [")
    for i in all_strucs:
        if i != all_strucs[-1]:
            f.write(f"{i}, ")
        else:
            f.write(f"{i}")
    f.write("]\n")
    f.write("  kinetics: gas\n")
    f.write("  reactions: all\n")
    f.write("  state: {T: 298, P: ")
    f.write(f"{Pressure}, X: ")
    f.write("{'O':1, 'OCC1OC(O)C(O)C(O)C1O':1")
    for i in all_strucs:
        if i not in initial_strucs:
            f.write(f", {i}:0")
    f.write("}}\n\n")
    f.write("species:\n")
    for i in all_strucs:
        ele = get_elements(i)
        f.write(f"- name: '{i}'\n")
        f.write(f"  composition: {ele}\n")
        f.write("  thermo:\n")
        f.write("    model: constant-cp\n")
        f.write("  equation-of-state: ideal-gas\n\n")
    f.write("reactions:\n")
    for idx,row in df.iterrows():
        Rsmi = [mol2smi(smi2mol(i)) for i in row['Rsmi_canon'].split('.')]
        Psmi = [mol2smi(smi2mol(i)) for i in row['Psmi_canon'].split('.')]
        dG = row['DG']
        A,b = get_rate_constant(idx,Rsmi,Pressure)
        f.write("- equation: ")
        for smi in Rsmi:
            if smi == Rsmi[-1]:
                f.write(f"{i}")
            else:
                f.write(f"{i} + ")
        f.write(" => ")
        for smi in Psmi:
            if smi == Psmi[-1]:
                f.write(f"{i}\n")
            else:
                f.write(f"{i} + ")
        f.write("  rate-constant: {A: ")
        f.write(f"{A}, b: {b}, Ea: {dG}")
        f.write("}\n")

f.seek(0)
for idx,line in enumerate(f):
    if idx < 50:
        print(line)
            
            
