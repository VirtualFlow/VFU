#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 21:56:20 2022

@author: akshat
"""
import os 
import subprocess
from lig_process import process_ligand

# smi = 'C1CC(CCC1NC(=O)COC2=CC=C(C=C2)Cl)NC(=O)COC3=CC=C(C=C3)Cl'
smi = 'NC1CCC(CC1)NC(=O)COC1=CC=C(Cl)C=C1'

protein_file = ''

# Process the ligand; convert smi->3d
process_ligand(smi, 'sdf')


# Generate a list of all 3d tautomers: 
lig_locations = os.listdir('./ligands/')

results = {}

for lig_ in lig_locations:
    
    # Run smina with autobox parameter: 
    # command = [ './executables/smina',                      # Docking executable
    #            '-r', './config/prot_1.pdb',                 # Defining the receptor
    #            '-l', './ligands/{}'.format(lig_),           # Defining the ligand to be docked
    #            '--exhaustiveness', '10',                    # Define exhaustiveness of search
    #            '--autobox_ligand', './config/lig_1.pdb',    # Define the reference ligand
    #            '--autobox_add', '3',                        # Space around ligand for docking
    #            '-o', './outputs/pose_{}.pdb'.format(lig_.split('.')[0]),     # Path of output poses 
    #            '--log', './outputs/log_{}.txt'.format(lig_.split('.')[0])    # Path of output log file 
    #            ]
    
    # Run smina with custom box parameters    
    # command = [ './executables/smina',                      # Docking executable
    #            '-r', './config/prot_1.pdb',                 # Defining the receptor
    #            '-l', './ligands/{}'.format(lig_),           # Defining the ligand to be docked
    #            '--exhaustiveness', '10',                    # Define exhaustiveness of search
               
    #            '--center_x', '-16',                         # Define center for search space (x-axis)
    #            '--center_y', '145',                         # Define center for search space (y-axis)
    #            '--center_z', '27',                          # Define center for search space (z-axis)
    #            '--size_x', '10',                            # Define the length of the search space box (x-axis)
    #            '--size_y', '10',                            # Define the length of the search space box (y-axis)
    #            '--size_z', '10',                            # Define the length of the search space box (z-axis)

    #            '-o', './outputs/pose_{}.pdb'.format(lig_.split('.')[0]),     # Path of output poses 
    #            '--log', './outputs/log_{}.txt'.format(lig_.split('.')[0])    # Path of output log file 
    #            ]
    
    
    # Run smina with flexible side chains:     
    command = [ './executables/smina',                      # Docking executable
               '-r', './config/prot_1.pdb',                 # Defining the receptor
               '-l', './ligands/{}'.format(lig_),           # Defining the ligand to be docked
               '--exhaustiveness', '10',                    # Define exhaustiveness of search
               
               '--center_x', '-16',                         # Define center for search space (x-axis)
               '--center_y', '145',                         # Define center for search space (y-axis)
               '--center_z', '27',                          # Define center for search space (z-axis)
               '--size_x', '10',                            # Define the length of the search space box (x-axis)
               '--size_y', '10',                            # Define the length of the search space box (y-axis)
               '--size_z', '10',                            # Define the length of the search space box (z-axis)

               '-o', './outputs/pose_{}.pdb'.format(lig_.split('.')[0]),     # Path of output poses 
               '--log', './outputs/log_{}.txt'.format(lig_.split('.')[0])    # Path of output log file 
               ]
    
    
    command_run = subprocess.run(command, capture_output=True)
    
    
    # Check the quality of generated structure (some post-processing quality control): 
    command = ['obenergy', './outputs/pose_{}.pdb'.format(lig_.split('.')[0])]
    command_obabel_check = subprocess.run(command, capture_output=True)
    command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
    total_energy         = float(command_obabel_check.split(' ')[-2])


    if total_energy < 10000: 
        # Collect the docking output: 
        docking_out = command_run.stdout.decode("utf-8")
        A = docking_out.split('\n')
        A = [x for x in A if '       ' in x]
        docking_score = [] # Results will be stored here!
        for item in A: 
            line_ = item.split('       ')
            if len(line_) == 3 and '      ' in line_[-1]: 
                docking_score.append(float(line_[1]))
        
                        # Docking Scores for all poses, Pose file, Log File               
        results[lig_] = [docking_score, './outputs/pose_{}.pdb'.format(lig_.split('.')[0]), './outputs/log_{}.txt'.format(lig_.split('.')[0])]
    else: 
        results[lig_] = 'Extremely high pose energy encountered.'
        
    raise Exception('T')