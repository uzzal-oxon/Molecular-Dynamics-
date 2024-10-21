#!/bin/bash

# Author: [Your Name]
# Date: $(date +%Y-%m-%d)
# Description: GROMACS MD simulation pipeline for protein in water

# --- Input Files ---
PDB_FILE="protein.pdb"  # Input protein structure file
IONS_MDP="ions.mdp"    # MDP file for ion generation
MINIM_MDP="minim.mdp"  # MDP file for energy minimization
NVT_MDP="nvt.mdp"      # MDP file for NVT equilibration
NPT_MDP="npt.mdp"      # MDP file for NPT equilibration
MD_MDP="md.mdp"        # MDP file for final MD simulation

# --- Forcefield and Topology generation ---
echo "Generating topology..."
gmx pdb2gmx -f "${PDB_FILE}" -ignh -o molecule_ff.pdb -p topol.top << EOF
8
1
EOF

# --- Box creation ---
echo "Creating simulation box..."
gmx editconf -f molecule_ff.pdb -o molecule_box.pdb -d 0.5 -bt triclinic

# --- System solvation ---
echo "Solvating the system..."
gmx solvate -cp molecule_box.pdb -cs spc216.gro -p topol.top -o molecule_solv.pdb

# --- Ion addition ---
echo "Adding ions..."
gmx grompp -f "${IONS_MDP}" -c molecule_solv.pdb -o ions.tpr -p topol.top
gmx genion -s ions.tpr -o molecule_ions.pdb -p topol.top -pname NA -nname CL -neutral -conc 0.1 << EOF
13
EOF

# --- Energy minimization ---
echo "Performing energy minimization..."
gmx grompp -f "${MINIM_MDP}" -c molecule_ions.pdb -p topol.top -o em.tpr
gmx mdrun -v -deffnm em

# --- NVT equilibration ---
echo "Running NVT equilibration..."
gmx grompp -f "${NVT_MDP}" -c em.gro -r em.gro -p topol.top -o nvt.tpr
gmx mdrun -v -deffnm nvt

# --- NPT equilibration ---
echo "Running NPT equilibration..."
gmx grompp -f "${NPT_MDP}" -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
gmx mdrun -v -deffnm npt

# --- Final MD simulation ---
echo "Running final MD simulation..."
gmx grompp -f "${MD_MDP}" -c npt.gro -t npt.cpt -p topol.top -o md.tpr
gmx mdrun -v -deffnm md

echo "MD simulation completed!"