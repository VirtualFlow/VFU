#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep 24 17:17:13 2022

@author: akshat
"""
import os 
import sys
import subprocess
from lig_process import process_ligand
from utils import run_plants_docking, run_autodock_gpu_docking, run_EquiBind, run_rDock, run_leDock, process_idock_output
command = []

# Parameters:  
is_selfies     = False 
program_choice = 'idock' # smina/qvina/qvina-w/vina/vina_carb/vina_xb/gwovina/PLANTS/autodock_gpu/autodock_cpu/EquiBind/rDock/gnina/ledock/idock
receptor       = './config/prot_1.pdb'
# smi            = 'C1CC(CCC1NC(=O)COC2=CC=C(C=C2)Cl)NC(=O)COC3=CC=C(C=C3)Cl'
smi            = 'BrC=CC1OC(C2)(F)C2(Cl)C1.CC.[Cl][Cl]'



if is_selfies == True: 
    import selfies 
    try: 
        smi = selfies.decoder(smi)
    except: 
        raise Exception('Invalid SELFIES provided. Please make sure that the varibale smi contains a valid selfies string. ')

# Check if receptor file exists inside the config directory: 
files_ls = os.listdir('./config/')
if receptor.split('/')[-1] not in files_ls: 
    raise Exception('Receptor file not found inside the config directory. Please makes sure that the designated receptor file ({}) exists inside config'.format(receptor.split('/')[-1]))


# Docking search paramters: 
exhaustiveness = 10
center_x       = -16                       # Define center for search space (x-axis)
center_y       = 145                       # Define center for search space (y-axis)
center_z       = 27                        # Define center for search space (z-axis)
size_x         = 10                        # Define the length of the search space box (x-axis)
size_y         = 10                        # Define the length of the search space box (y-axis)
size_z         = 10                        # Define the length of the search space box (z-axis)

results = {}  # Storage for results





if program_choice == 'PLANTS': 
    results = run_plants_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
    sys.exit()
if program_choice == 'autodock_gpu': 
    results = run_autodock_gpu_docking(receptor, smi, program_choice)
    sys.exit()
if program_choice == 'autodock_cpu': 
    results = run_autodock_gpu_docking(receptor, smi, program_choice)
    sys.exit()
if program_choice == 'EquiBind': 
    results = run_EquiBind(receptor, smi)
    sys.exit()
if program_choice == 'rDock':   
    run_rDock(receptor, smi)
    sys.exit()
if program_choice == 'ledock': 
    run_leDock(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
    sys.exit()

    

# Assign the right program for docking:  
command.append('./executables/{}'.format(program_choice))


# Assign the receptor for docking:  
file_type_check = receptor.split('.')[-1]
if not(file_type_check == 'pdb' or file_type_check == 'pdbqt'): 
    raise Exception('invalid receptor file type provided!')
if program_choice == 'smina' or program_choice == 'gnina': 
    command.append('-r')
    command.append(receptor)
elif program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina': 
    command.append('--receptor')
    command.append(receptor)


# Assign the right ligand for docking
if program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'smina' or program_choice == 'gnina' or program_choice == 'qvina-w' or program_choice == 'qvina-w' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina':
    process_ligand(smi, 'pdbqt')
lig_locations = os.listdir('./ligands/')


for lig_ in lig_locations: 
    
    # Add in the ligand file and the exhaustiveness setting
    if program_choice == 'idock' or  program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina': 
        if program_choice == 'idock': 
            cmd = command + ['--ligand', './ligands/{}'.format(lig_)]
        else: 
            cmd = command + ['--ligand', './ligands/{}'.format(lig_), '--exhaustiveness', str(exhaustiveness)]
    elif program_choice == 'smina' or program_choice == 'gnina': 
        cmd = command + ['-l', './ligands/{}'.format(lig_), '--exhaustiveness', str(exhaustiveness)]
        
    # Add in the docking parameters: 
    cmd = cmd + ['--center_x', str(center_x)]
    cmd = cmd + ['--center_y', str(center_y)]
    cmd = cmd + ['--center_z', str(center_z)]
    cmd = cmd + ['--size_x', str(size_x)]
    cmd = cmd + ['--size_y', str(size_y)]
    cmd = cmd + ['--size_z', str(size_z)]

    # Add in parameters for generating output files: 
    if  program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina': 
        cmd = cmd + ['--out', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
    elif program_choice == 'smina' or program_choice == 'gnina': 
        cmd = cmd + ['-o', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
    
    if program_choice != 'idock': 
        cmd = cmd + ['--log', './outputs/log_{}.txt'.format(lig_.split('.')[0])]

    # Run the command: 
    command_run = subprocess.run(cmd, capture_output=True)
    
    
    if program_choice == 'idock': 
        process_idock_output(results)
        sys.exit()

    # Check the quality of generated structure (some post-processing quality control):
    # TODO: Make a function out of this! 
    try: 
        ob_cmd = ['obenergy', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 

    if total_energy < 10000: 
        # Collect the docking output: 
        docking_out = command_run.stdout.decode("utf-8")
        A = docking_out.split('\n')
        docking_score = []
        for item in A: 
            line_split = item.split(' ')
            line_split = [x for x in line_split if x != '']
            if len(line_split) == 4: 
                try: 
                    vr_1 = float(line_split[0])
                    vr_2 = float(line_split[1])
                    vr_3 = float(line_split[2])
                    vr_4 = float(line_split[3])
                    docking_score.append(vr_2)
                except: 
                    continue
        
        
                        # Docking Scores for all poses, Pose file, Log File               
        results[lig_] = [docking_score, './outputs/pose_{}.pdb'.format(lig_.split('.')[0]), './outputs/log_{}.txt'.format(lig_.split('.')[0])]
    else: 
        results[lig_] = 'Extremely high pose energy encountered.'

    