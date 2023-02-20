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
from utils import run_plants_docking, run_autodock_gpu_docking, run_EquiBind, run_rDock, run_leDock, process_idock_output, run_adfr_docking, run_flexx_docking
from utils import check_energy, run_ligand_fit, run_mcdock, run_AutodockZN, run_GalaxyDock3, run_dock6, run_fred_docking, run_iGemDock, perform_gold_docking
from utils import run_glide_docking, run_rosetta_docking, run_mdock_docking, run_seed_docking, run_molegro_docking, run_fitdock_docking
from utils import  run_lightdock_docking, run_RLDock_docking, run_MpSDockZN_docking


def run_pose_prediction_program(program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, receptor): 
    command = []
    
    # Check if receptor file exists inside the config directory: 
    files_ls = os.listdir('./config/')
    if receptor.split('/')[-1] not in files_ls: 
        raise Exception('Receptor file not found inside the config directory. Please makes sure that the designated receptor file ({}) exists inside /config/'.format(receptor.split('/')[-1]))
        
    if os.path.exists('./ligands') == False: 
        subprocess.run(['mkdir', './ligands'])
    
    if program_choice == 'flexx': 
        results = run_flexx_docking(receptor, smi)
        return results
    if program_choice == 'MCDock': 
        results = run_mcdock(receptor, smi)
        return results
    if program_choice == 'dock6': 
        results = run_dock6(receptor, smi)
        return results
    if program_choice == 'iGemDock': 
        results = run_iGemDock(receptor, smi, exhaustiveness)
        return results
    if program_choice == 'M-Dock': 
        results = run_mdock_docking(receptor, smi)
        return results
    if program_choice == 'SEED': 
        results = run_seed_docking(receptor, smi)
        return results
    if program_choice == 'molegro': 
        results = run_molegro_docking(receptor, smi)
        return results
    if program_choice == 'FitDock': 
        results = run_fitdock_docking(receptor, smi)
        return results
    if program_choice == 'LightDock': 
        results = run_lightdock_docking(receptor, smi, exhaustiveness)
        return results
    if program_choice == 'RLDock': 
        results = run_RLDock_docking(receptor, smi, exhaustiveness)
        return results
    if program_choice == 'MpSDockZN': 
        results = run_MpSDockZN_docking(receptor, smi)
        return results
    if program_choice == 'LigandFit': 
        results = run_ligand_fit(receptor, smi, center_x, center_y, center_z)
        return results
    if program_choice == 'AutodockZN': 
        results = run_AutodockZN(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness)
        return results
    if program_choice == 'GalaxyDock3': 
        results = run_GalaxyDock3(receptor, smi, center_x, center_y, center_z, exhaustiveness)
        return results
    if program_choice == 'FRED': 
        results = run_fred_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness)
        return results
    if program_choice == 'gold': 
        results = perform_gold_docking(receptor, smi, size_x, size_y, size_z, center_x, center_y, center_z)
        return results
    if program_choice == 'glide':
        results = run_glide_docking(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi)
        return results
    if program_choice == 'rosetta-ligand':
        results = run_rosetta_docking(receptor, smi, center_x, center_y, center_z, exhaustiveness)
        return results
    if program_choice == 'PLANTS': 
        results = run_plants_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
        return results
    if program_choice == 'autodock_gpu': 
        results = run_autodock_gpu_docking(receptor, smi, program_choice)
        return results
    if program_choice == 'autodock_cpu': 
        results = run_autodock_gpu_docking(receptor, smi, program_choice)
        return results
    if program_choice == 'EquiBind': 
        results = run_EquiBind(receptor, smi)
        return results
    if program_choice == 'rDock':   
        results = run_rDock(receptor, smi)
        return results
    if program_choice == 'ledock': 
        run_leDock(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
        return results
    if program_choice == 'adfr': 
        results = run_adfr_docking(receptor, smi)
        return results
        
    results = {}  # Storage for results
        
    # Assign the right program for docking:  
    command.append('./executables/{}'.format(program_choice))
    if program_choice == 'AutoDock-Koto': 
        program_choice = 'AutodockVina_1.2' # Run everything in the same setting as AutodockVina1.2
    
    # Assign the receptor for docking:  
    file_type_check = receptor.split('.')[-1]
    if not(file_type_check == 'pdb' or file_type_check == 'pdbqt'): 
        raise Exception('invalid receptor file type provided!')
    if program_choice == 'smina' or program_choice == 'gnina': 
        command.append('-r')
        command.append(receptor)
    elif program_choice == 'AutodockVina_1.2' or program_choice == 'autodock_vina' or program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina'  or program_choice == 'PSOVina': 
        command.append('--receptor')
        command.append(receptor)
    
    
    # Assign the right ligand for docking
    if  program_choice == 'AutodockVina_1.2' or program_choice == 'autodock_vina' or program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'smina' or program_choice == 'gnina' or program_choice == 'qvina-w' or program_choice == 'qvina-w' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina':
        process_ligand(smi, 'pdbqt')
    lig_locations = os.listdir('./ligands/')
    
    
    for lig_ in lig_locations: 
        
        # Add in the ligand file and the exhaustiveness setting
        if program_choice == 'AutodockVina_1.2' or program_choice == 'autodock_vina' or program_choice == 'idock' or  program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina': 
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
        if program_choice == 'AutodockVina_1.2' or program_choice == 'autodock_vina' or program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina': 
            cmd = cmd + ['--out', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
        elif program_choice == 'smina' or program_choice == 'gnina': 
            cmd = cmd + ['-o', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]

        # Run the command:     
        command_run = subprocess.run(cmd, capture_output=True)
        
        if program_choice == 'idock': 
            process_idock_output(results)
            sys.exit()
                
        # Check the quality of generated structure (some post-processing quality control):
        total_energy = check_energy(lig_)
    
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
                        _ = float(line_split[0])
                        vr_2 = float(line_split[1])
                        _ = float(line_split[2])
                        _ = float(line_split[3])
                        docking_score.append(vr_2)
                    except: 
                        continue
                if program_choice == 'vina_carb': 
                    line_split = item.split(' ')
                    line_split = [x for x in line_split if x != '']
                    if len(line_split) == 6: 
                        try: 
                            line_split = item.split(' ')
                            line_split = [x for x in line_split if x != '']
                            docking_score.append(float(line_split[1]))
                        except: 
                            continue 
                            # Docking Scores for all poses, Pose file, Log File               
            results[lig_] = [docking_score, './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
        else: 
            results[lig_] = 'Extremely high pose energy encountered/Unsuccessfull execution.'
        
        if 'VC_log.txt' in os.listdir(): 
            os.system('rm VC_log.txt')
    
    return results
