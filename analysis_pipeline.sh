#!/bin/bash

# Author: [Your Name]
# Date: $(date +%Y-%m-%d)
# Description: GROMACS analysis pipeline for protein MD trajectory

# --- Input Files ---
TPR_FILE="md.tpr"      # Input GROMACS run input file
XTC_FILE="md.xtc"      # Input trajectory file

# --- Trajectory Processing ---
echo "Processing trajectory..."
gmx trjconv -s "${TPR_FILE}" -f "${XTC_FILE}" -o md_new.xtc -center -pbc mol -ur compact << EOF
Protein
System
EOF

# --- RMSD Analysis ---
echo "Calculating RMSD..."
gmx rms -s "${TPR_FILE}" -f md_new.xtc -o RMSD_Protein.xvg -tu ns << EOF
Backbone
Backbone
EOF

# --- RMSF Analysis ---
echo "Calculating RMSF..."
gmx rmsf -s "${TPR_FILE}" -f md_new.xtc -o RMSF.xvg -res << EOF
Protein
EOF

# --- Radius of Gyration Analysis ---
echo "Calculating Radius of Gyration..."
gmx gyrate -s "${TPR_FILE}" -f md_new.xtc -o Rg.xvg << EOF
Protein
EOF

# --- Hydrogen Bond Analysis ---
echo "Calculating Hydrogen Bonds..."
echo "q" | gmx make_ndx -f "${TPR_FILE}" -o index.ndx
gmx select -s "${TPR_FILE}" -on index.ndx -select "Protein or resname SOL" 
gmx hbond -f md_new.xtc -n index.ndx -s "${TPR_FILE}" -num Hbond.xvg -tu ns << EOF
Protein
SOL
EOF

# --- Solvent Accessible Surface Area Analysis ---
echo "Calculating SASA..."
gmx sasa -f md_new.xtc -s "${TPR_FILE}" -o SASA.xvg -tu ns << EOF
Protein
EOF

echo "Analysis completed!"