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

from pose_prediction import run_plants_docking, run_autodock_gpu_docking, run_EquiBind, run_rDock, run_leDock, process_idock_output, run_adfr_docking, run_flexx_docking
from pose_prediction import check_energy, run_ligand_fit, run_mcdock, run_AutodockZN, run_GalaxyDock3, run_dock6, run_fred_docking, run_iGemDock, perform_gold_docking
from pose_prediction import run_rosetta_docking, run_mdock_docking, run_seed_docking, run_molegro_docking, run_fitdock_docking
from pose_prediction import  run_lightdock_docking, run_RLDock_docking, run_MpSDockZN_docking, run_CovDock_docking, run_Glide_HTVS, run_Glide_SP, run_Glide_XP
from pose_prediction import perform_HDock_docking

from scoring_functions import run_nnscore2, run_rf_scoring, run_smina_scoring, run_ad4_scoring, run_vinandro_scoring, run_vina_scoring, run_gnina_scoring
from scoring_functions import run_PLANTS_chemplp_scoring, run_PLANTS_plp_scoring, run_PLANTS_plp95_scoring, contact_score, continuous_score, grid_score
from scoring_functions import gold_asp_scoring, gold_chemscore_scoring, gold_goldscore_scoring, gold_plp_scoring
from scoring_functions import run_mm_gbsa, Hawkins_gbsa

def run_pose_prediction_program(program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, receptor, chimera_path = '', dock6_path = '', ref_lig = '', seed_path = '', mdock_path = '', molegro_path = '', covalent_bond_constraints='', receptor_2=''): 
    '''
    This function runs docking simulations on a given ligand for a specified receptor using the chosen docking program.

    Parameters:    
        program_choice (str): The name of the docking program to use.
        center_x (int): The x-coordinate of the center of the search space.
        center_y (int): The y-coordinate of the center of the search space.
        center_z (int): The z-coordinate of the center of the search space.
        size_x (int): The size of the search space along the x-axis.
        size_y (int): The size of the search space along the y-axis.
        size_z (int): The size of the search space along the z-axis.
        exhaustiveness (int): The exhaustiveness parameter to use.
        smi (str): The SMILES representation of the ligand to be docked.
        receptor (str): The path to the receptor file for docking.
        chimera_path (str): Path in system for Chimera application
        dock6_path (str): Path in system for dock6 application
        ref_lig (str): Reference ligand file, required by some docking programs
        seed_path (str): Path in system for SEED application
        mdock_path (str): Path in system for MDOCK software
        molegro_path (str): Path in system for Molegro software
        
    Returns:
        A dictionary containing the results of the docking simulations, with keys being the names of the ligand files and values being lists of docking scores and the path to the docked pose file.
        
    Raises:    
        Exception: If the receptor file is not found in the config/ directory.
        Exception: If the receptor file type is invalid.

    Example usage:
        results = run_pose_prediction_program('AutoDock-Koto', 10, 10, 10, 20, 20, 20, 8, 'C1=CC(=CC=C1CSCC2C(C(C(O2)N3C=NC4=C(N=CN=C43)N)O)O)Cl', './config/receptor.pdbqt')
    '''
    command = []
    
    # Check if receptor file exists inside the config directory: 
    files_ls = os.listdir('./config/')
    if receptor.split('/')[-1] not in files_ls: 
        raise Exception('Receptor file not found inside the config directory. Please makes sure that the designated receptor file ({}) exists inside /config/'.format(receptor.split('/')[-1]))
        
    if os.path.exists('./ligands') == False: 
        subprocess.run(['mkdir', './ligands'])
    if os.path.exists('./outputs') == False: 
        subprocess.run(['mkdir', './outputs'])
        
    if program_choice == 'flexx': 
        results = run_flexx_docking(receptor, smi, ref_lig)
        return results
    if program_choice == 'MCDock': 
        results = run_mcdock(receptor, smi)
        return results
    if program_choice == 'dock6': 
        results = run_dock6(receptor, smi, chimera_path, dock6_path, ref_lig)
        return results
    if program_choice == 'iGemDock': 
        results = run_iGemDock(receptor, smi, exhaustiveness)
        return results
    if program_choice == 'M-Dock': 
        results = run_mdock_docking(receptor, smi, mdock_path, ref_lig)
        return results
    if program_choice == 'SEED': 
        results = run_seed_docking(receptor, smi, chimera_path, seed_path)
        return results
    if program_choice == 'molegro': 
        results = run_molegro_docking(receptor, smi, ref_lig, molegro_path)
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
        results = run_rDock(receptor, smi, ref_lig)
        return results
    if program_choice == 'ledock': 
        run_leDock(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
        return results
    if program_choice == 'adfr': 
        results = run_adfr_docking(receptor, smi)
        return results
    if program_choice == 'CovDock': 
        results =  run_CovDock_docking(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi, covalent_bond_constraints)
        return results
    if program_choice == 'GlideHTVS': 
        results =  run_Glide_HTVS(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi)
        return results
    if program_choice == 'GlideXP': 
        results =  run_Glide_XP(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi)
        return results
    if program_choice == 'GlideSP': 
        results =  run_Glide_SP(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi)
        return results
    if program_choice == 'HDock': 
        results =  perform_HDock_docking(receptor, receptor_2)
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
    elif program_choice == 'vina_gpu' or program_choice == 'vina_gpu_2.0' or program_choice == 'qvina_w_gpu' or program_choice == 'qvina_gpu' or program_choice == 'AutodockVina_1.2' or program_choice == 'AutodockVina_1.1.2' or program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina'  or program_choice == 'PSOVina': 
        command.append('--receptor')
        command.append(receptor)
    
    
    # Assign the right ligand for docking
    if  program_choice == 'vina_gpu' or program_choice == 'vina_gpu_2.0' or program_choice == 'qvina_w_gpu' or program_choice == 'qvina_gpu' or program_choice == 'AutodockVina_1.2' or program_choice == 'AutodockVina_1.1.2' or program_choice == 'idock' or program_choice == 'qvina' or program_choice == 'smina' or program_choice == 'gnina' or program_choice == 'qvina-w' or program_choice == 'qvina-w' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina':
        process_ligand(smi, 'pdbqt')
    lig_locations = os.listdir('./ligands/')
    
    
    for lig_ in lig_locations: 
        
        # Add in the ligand file and the exhaustiveness setting
        if program_choice == 'vina_gpu' or program_choice == 'vina_gpu_2.0' or program_choice == 'qvina_w_gpu' or program_choice == 'qvina_gpu' or program_choice == 'AutodockVina_1.2' or program_choice == 'AutodockVina_1.1.2' or program_choice == 'idock' or  program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina': 
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
        if program_choice == 'vina_gpu' or program_choice == 'vina_gpu_2.0' or program_choice == 'qvina_w_gpu' or program_choice == 'qvina_gpu' or program_choice == 'AutodockVina_1.2' or program_choice == 'AutodockVina_1.1.2' or program_choice == 'qvina' or program_choice == 'qvina-w' or program_choice == 'vina' or program_choice == 'vina_carb' or program_choice == 'vina_xb' or program_choice == 'gwovina' or program_choice == 'PSOVina': 
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

def run_scoring_prediction_program(scoring_function, ligand_path, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, receptor, chimera_path = '', dock6_path = ''): 
    """Runs a molecular scoring program on a docked ligand-receptor pair and returns the resulting score(s).

    Args:
        scoring_function (str): The name of the scoring program to run. 
        ligand_path (str): The file path to the ligand to be docked or scored.
        center_x (float): The x-coordinate of the center of the binding site.
        center_y (float): The y-coordinate of the center of the binding site.
        center_z (float): The z-coordinate of the center of the binding site.
        size_x (float): The x-dimension of the box used for the docking or scoring program.
        size_y (float): The y-dimension of the box used for the docking or scoring program.
        size_z (float): The z-dimension of the box used for the docking or scoring program.
        exhaustiveness (int): The exhaustiveness parameter used for some docking programs.
        smi (str): The SMILES string of the ligand to be docked or scored.
        receptor (str): The file path to the receptor to which the ligand should be docked or scored.
        chimera_path (str): Path in system for Chimera application
        dock6_path (str): Path in system for dock6 application
    Returns:
        scores (float or strings): The score(s) resulting from the molecular docking or scoring program. 
    """
    if scoring_function == 'nnscore2': 
        scores = run_nnscore2(receptor, ligand_path)
    if scoring_function == 'rf-score': 
        scores = run_rf_scoring(receptor, ligand_path) 
    if scoring_function == 'smina-scoring': 
        scores = run_smina_scoring(receptor, ligand_path)
    if scoring_function == 'ad4_scoring': 
        scores = run_ad4_scoring(receptor, ligand_path)
    if scoring_function == 'vinandro_scoring': 
        scores = run_vinandro_scoring(receptor, ligand_path)
    if scoring_function == 'vina_scoring': 
        scores = run_vina_scoring(receptor, ligand_path)
    if scoring_function == 'gnina_scoring': 
        scores = run_gnina_scoring(receptor, ligand_path)
    if scoring_function == 'chemplp_scoring': 
        scores = run_PLANTS_chemplp_scoring(receptor, ligand_path)
    if scoring_function == 'PLP_scoring': 
        scores = run_PLANTS_plp_scoring(receptor, ligand_path)
    if scoring_function == 'PLP95_scoring': 
        scores = run_PLANTS_plp95_scoring(receptor, ligand_path)
    if scoring_function == 'contact_scoring': 
        scores = contact_score(receptor, chimera_path, dock6_path, ligand_path, center_x, center_y, center_z, size_x, size_y, size_z)
    if scoring_function == 'continuous_scoring': 
        scores = continuous_score(receptor, chimera_path, dock6_path, ligand_path)
    if scoring_function == 'grid_scoring': 
        scores = grid_score(receptor, chimera_path, dock6_path, ligand_path, center_x, center_y, center_z, size_x, size_y, size_z)
    if scoring_function == 'asp':
        scores = gold_asp_scoring(receptor_filepath=receptor, ligand_filepath=ligand_path)
    if scoring_function == 'chemscore':
        scores = gold_chemscore_scoring(receptor_filepath=receptor, ligand_filepath=ligand_path)
    if scoring_function == 'goldscore':
        scores = gold_goldscore_scoring(receptor_filepath=receptor, ligand_filepath=ligand_path)
    if scoring_function == 'plp':
        scores = gold_plp_scoring(receptor_filepath=receptor, ligand_filepath=ligand_path)
    if scoring_function == 'mm_gbsa_scoring': 
        scores = run_mm_gbsa(chimera_path, ligand_path, receptor)
    if scoring_function == "Hawkins_gbsa":
        scores = Hawkins_gbsa(receptor, chimera_path, dock6_path, ligand_path, center_x, center_y, center_z, size_x, size_y, size_z)

                     
    return scores
    
    



