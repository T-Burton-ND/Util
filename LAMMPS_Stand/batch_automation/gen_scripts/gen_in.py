# Python Script to Generate LAMMPS Input File from YAML Config
# Usage: python generate_input.py config.yaml

import yaml
import sys

# ============================
# Load Configuration Variables
# ============================

with open(sys.argv[1], 'r') as file:
    config = yaml.safe_load(file)

# -------- File and Naming --------
system_file = config['system_file']
settings_name = config['settings_name']
prefix = config['prefix']
output_file = config['lammps_input']

# -------- Random Seed --------
vseed = config['vseed']

# -------- Simulation Parameters --------
units = config['units']
atom_style = config['atom_style']
boundary = config['boundary']
bond_style = config['bond_style']
angle_style = config['angle_style']
dihedral_style = config['dihedral_style']
improper_style = config['improper_style']
kspace_style = config['kspace_style']
pair_style = config['pair_style']
pair_modify = config['pair_modify']
special_bonds = config['special_bonds']

# -------- Neighbor Settings --------
neighbor_bin = config['neighbor_bin']
neigh_modify = config['neigh_modify']
neighbor_bin_npt = config['neighbor_bin_npt']

#-------- Simulation Size --------
rescale = config['rescale']
replicate = config['replicate']

# -------- Thermo and Dump Frequencies --------
thermo_freq = config['thermo_freq']
avg_freq = config['avg_freq']
dump4avg = config['dump4avg']
coords_freq = config['coords_freq']

# -------- MSD Parameters --------
msd_freq = config['msd_freq']
msd_ave = config['msd_ave']
msd_tot = config['msd_tot']
msd_flush = config.get('msd_flush', True)

# -------- RDF Parameters --------
rdf_bins = config['rdf_bins']
rdf_freq = config['rdf_freq']
rdf_ave = config['rdf_ave']
rdf_tot = config['rdf_tot']
rdf_type_1 = config['rdf_type_1']
rdf_type_2 = config['rdf_type_2']

# -------- Relaxation Parameters --------
nSteps_relax = config['nSteps_relax']
relax_timestep = config['relax_timestep']

# -------- Energy Minimization Parameters --------
min_etol = config['min_etol']
min_ftol = config['min_ftol']
min_maxiter = config['min_maxiter']
min_maxeval = config['min_maxeval']

# -------- Thawing Parameters --------
thawing_steps = config['thawing_steps']
thawing_tstart = config['thawing_tstart']
thawing_tstop = config['thawing_tstop']
thawing_tdamp = config['thawing_tdamp']
thawing_tdamp_steps = config.get('thawing_tdamp_steps', 100)
Td_thaw = thawing_tdamp_steps*relax_timestep

# -------- Equilibration Parameters --------
press0_equil = config['press0_equil']
pressf_equil = config['pressf_equil']
temp0_equil = config['temp0_equil']
tempf_equil = config['tempf_equil']
nSteps_equil = config['nSteps_equil']
equil_timestep = config['equil_timestep']
equil_tdamp_steps   = config.get('equil_tdamp_steps', 100)
equil_pdamp_steps   = config.get('equil_pdamp_steps', 1000)
barostat_drag       = config.get('barostat_drag', 1.0)

# -------- Production Paramet --------
nSteps_prod = config['nSteps_prod']
prod_timestep = config['prod_timestep']
linear_momentum = config['linear_momentum']
production_neighbor_bin = config['production_neighbor_bin']



# ============================
# Generate LAMMPS Input File
# ============================



output = f"""
# LAMMPS Input File (Generated from YAML)

#===========================================================
# Variables
#===========================================================

# System Inputs:
variable    system_name     index    {system_file}
variable    settings_name    index    {settings_name}
variable    prefix           index    {prefix}
variable    replicate        index    {replicate}
variable    rescale         index    {rescale}

# Output Variables:
variable    thermo_freq      index    {thermo_freq}
variable    avg_freq         index    {avg_freq}
variable    dump4avg         index    {dump4avg}
variable    coords_freq      index    {coords_freq}

# MSD Variables:
variable    msd_freq         index    {msd_freq}
variable    msd_ave          index    {msd_ave}
variable    msd_tot          index    {msd_tot}

# RDF Parameters:
variable rdf_bins equal {rdf_bins}             # Number of histogram bins
variable rdf_freq equal {rdf_freq}              # Output every N steps
variable rdf_ave  equal {rdf_ave}               # Number of values to average
variable rdf_tot  equal {rdf_tot}             # Total steps per output window
variable rdf_type1 equal {rdf_type_1}          # Type 1 for RDF
variable rdf_type2 equal {rdf_type_2}          # Type 2 for RDF

# Random Seed:
variable    vseed            index    {vseed}

# Relaxation Variables:
variable    nSteps_relax    index    {nSteps_relax}
variable    relax_timestep index  {relax_timestep}

# Equilibration Variables:
variable    nSteps_equil    index    {nSteps_equil}
variable    press0_equil    index    {press0_equil}
variable    pressf_equil    index    {pressf_equil}
variable    temp0_equil     index    {temp0_equil}
variable    tempf_equil    index    {tempf_equil}

# Production Variables:
variable    nSteps_prod    index    {nSteps_prod}
variable    production_timestep index  {prod_timestep}

#===========================================================
# Initialize System
#===========================================================

units {units}
atom_style {atom_style}
boundary {boundary}
bond_style {bond_style}
angle_style {angle_style}
dihedral_style {dihedral_style}
improper_style {improper_style}
kspace_style   {kspace_style}
pair_style {pair_style}
pair_modify {pair_modify}
special_bonds lj/coul 0.0 0.0 0.5

#===========================================================
# Setup System
#===========================================================

read_data {system_file}
include {settings_name}

change_box all x scale {rescale} y scale {rescale} z scale {rescale} remap
replicate {replicate} {replicate} {replicate} bond/periodic

group all_atoms type > 0
group cl type 1
group li type 2

neighbor        {neighbor_bin} bin
neigh_modify    {neigh_modify}

thermo_style custom step temp press ke pe etotal density
thermo {thermo_freq}

#===========================================================
# Relaxation (unoverlap + minimization)
#===========================================================

timestep {relax_timestep}

# --- Disable long-range during relaxation ---
kspace_style none
pair_style    lj/cut 10.0
include      {settings_name}           

# Gentle contact relaxation
fix nve_relax all nve/limit 0.05          # WAS 0.05; a bit looser is safer

print "======================"
print "======================"
print "Starting relaxation..."
print "======================"
print "======================"
run {nSteps_relax}
unfix nve_relax

# -------------------- Minimization -----------------------
print "======================"
print "======================"
print "Starting minimization..."
print "======================"
print "======================"

min_style cg
minimize {min_etol} {min_ftol} {min_maxiter} {min_maxeval}

# Restore long-range & settings (idempotent)
kspace_style {kspace_style}
pair_style   {pair_style}
include  {settings_name}

#===========================================================
# Thawing Phase (NVT Ramp to Relieve Pressure)
#===========================================================
velocity all_atoms create {thawing_tstart} {vseed} dist gaussian

# Convert step-based damping to time using dt

fix nvt_thaw all nvt temp {thawing_tstart} {thawing_tstop} {Td_thaw}

print "======================"
print "======================"
print "Starting NVT thawing..."
print "======================"
print "======================"

run {thawing_steps}
unfix nvt_thaw

#===========================================================
# Equilibration (NPT, Nose-Hoover)
#===========================================================

timestep {equil_timestep}

# Slightly larger neighbor skin during cell changes
neighbor        {neighbor_bin_npt} bin
neigh_modify    delay 0 every 1 check yes
comm_modify     cutoff 15.0 vel yes

variable Td_eq equal {equil_tdamp_steps}*{equil_timestep}
variable Pd_eq equal {equil_pdamp_steps}*{equil_timestep}

fix npt_equil all npt temp {temp0_equil} {tempf_equil} ${{Td_eq}} iso  {press0_equil} {pressf_equil} ${{Pd_eq}} drag {barostat_drag}

print "======================"
print "======================"
print "Starting equilibration..."
print "======================"
print "======================"

run {nSteps_equil}
unfix npt_equil

write_data equilibrated.data

#===========================================================
# Production (NVE)
#===========================================================


# --------- MSD Output for Types 1 and 2 ---------
compute msd_cl cl msd
fix msd_tot_cl  cl ave/time {msd_freq} {msd_ave} {msd_tot} c_msd_cl[4] file msd_cl.txt mode scalar

compute msd_li li msd
fix msd_tot_li  li ave/time {msd_freq} {msd_ave} {msd_tot} c_msd_li[4] file msd_li.txt mode scalar

# --------- RDF Output ---------
compute         rdf all rdf {rdf_bins} {rdf_type_1} {rdf_type_2}
fix             rdf_out all ave/time {rdf_freq} {rdf_ave} {rdf_tot} c_rdf[*] file rdf_{rdf_type_1}_{rdf_type_2}.txt.gz mode vector

# Trajectory dump (unwrapped)
dump production all custom {coords_freq} production.lammpstrj id mol type xu yu zu
dump_modify production format line "%d %d %d %.3f %.3f %.3f"
dump_modify production sort id

# Thermo (compact)
thermo_style custom step temp pe ke etotal press density
thermo_modify format float %10.3f
thermo 1000

# --------- Production Fix and Momentum Removal ---------
fix nve all nve
fix remove_momentum all momentum {linear_momentum} linear 1 1 1


# --------- Production Neighbor Update
neighbor        {production_neighbor_bin} bin
neigh_modify    delay 0 every 4 check yes

timestep {prod_timestep}

print "======================"
print "======================"
print "Starting production..."
print "======================"
print "======================"

run {nSteps_prod}

#===========================================================
# Post-Processing
#===========================================================

unfix nve
unfix remove_momentum
undump production
unfix rdf_out
unfix msd_tot_cl
unfix msd_tot_li
write_data postprod.data
"""

# ============================
# Write Output File
# ============================

with open(output_file, 'w') as f:
    f.write(output)

print(f"\nLAMMPS input file '{output_file}' created successfully.\n")
