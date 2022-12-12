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
from utils import check_energy, run_mm_gbsa, run_ligand_fit, run_mcdock, run_AutodockZN, run_GalaxyDock3, run_dock6, run_fred_docking, run_iGemDock, perform_gold_docking
from utils import run_glide_docking, run_rosetta_docking, run_mdock_docking, run_seed_docking, run_nnscore2, run_rf_scoring, run_molegro_docking, run_fitdock_docking
from utils import run_smina_scoring, run_gnina_scoring, run_lightdock_docking, run_RLDock_docking

command = []

# Parameters:  
is_selfies     = False 
is_peptide     = False
program_choice = 'PSOVina' # smina/qvina/qvina-w/vina/vina_carb/vina_xb/gwovina/PLANTS/autodock_gpu/autodock_cpu/EquiBind/rDock/gnina/ledock/idock
                                 # /autodock_vina/adfr/AutodockVina_1.2/AutodockZN/flexx/MM-GBSA/MCDock/LigandFit/GalaxyDock3/dock6/FRED/iGemDock/gold
                                 # /glide/rosetta-ligand/M-Dock/SEED/nnscore2/rf-score/molegro/FitDock/PSOVina/smina-scoring/gnina-scoring/AutoDock-Koto
                                 # /LightDock/RLDock
                           
receptor       = './config/5wiu_test.pdbqt'
if program_choice == 'nnscore2': 
    run_nnscore2(receptor)
if program_choice == 'rf-score': 
    rescored_pose, rf_scores = run_rf_scoring(receptor)
if program_choice == 'smina-scoring': 
    score = run_smina_scoring(receptor)
if program_choice == 'gnina-scoring': 
    score = run_gnina_scoring(receptor)

smi            = 'C=C=C=C' # 'BrC=CC1OC(C2)(F)C2(Cl)C1.CC.[Cl][Cl]', 'C=C=C=C', 'SQETFSDLWKLLPEN'
exhaustiveness = 10


if is_selfies == True: 
    import selfies 
    try: 
        smi = selfies.decoder(smi)
    except: 
        raise Exception('Invalid SELFIES provided. Please make sure that the varibale smi contains a valid selfies string. ')
if is_peptide == True: 
    
    amino_acids = {'G', 'A', 'L', 'M', 'F', 'W', 'K', 'Q', 'E', 'S', 'P', 'V', 'I', 'C', 'Y', 'H', 'R', 'N', 'D', 'T'}    
    for aa in smi:
        if aa not in amino_acids:
            raise Exception('Non-standard/missing amino acid:', aa)
    import rdkit 
    from rdkit import Chem
    mol = Chem.MolFromSequence(smi)
    if mol == None: 
        raise Exception('RDKit unable to recognize amino acid sequence:', smi)
    smi = Chem.MolToSmiles(mol)
    
if 'Si' in smi or 'Sn' in smi or 'B' in smi: 
    raise Exception('Si/Sn/B found. Please be careful when using these for docking. Namely, the original AutoDock programs do not suport these atom types. Please make sure that the selected docking program support these atoms. Please remove this exception if everything is fine :)')

if program_choice == 'flexx': 
    results = run_flexx_docking(receptor, smi)
    sys.exit()
if program_choice == 'MM-GBSA': 
    results = run_mm_gbsa()
    sys.exit()
if program_choice == 'MCDock': 
    results = run_mcdock(receptor, smi)
    sys.exit()
if program_choice == 'dock6': 
    results = run_dock6(receptor, smi)
    sys.exit()
if program_choice == 'iGemDock': 
    results = run_iGemDock(receptor, smi, exhaustiveness)
    sys.exit()
if program_choice == 'M-Dock': 
    results = run_mdock_docking(receptor, smi)
    sys.exit()
if program_choice == 'SEED': 
    results = run_seed_docking(receptor, smi)
    sys.exit()
if program_choice == 'molegro': 
    results = run_molegro_docking(receptor, smi)
    sys.exit()
if program_choice == 'FitDock': 
    results = run_fitdock_docking(receptor, smi)
    sys.exit()
if program_choice == 'LightDock': 
    results = run_lightdock_docking(receptor, smi, exhaustiveness)
    sys.exit()
if program_choice == 'RLDock': 
    results = run_RLDock_docking(receptor, smi, exhaustiveness)
    sys.exit()


# Docking search paramters: 
center_x       = -17.820                      # Define center for search space (x-axis)
center_y       = 16.140                       # Define center for search space (y-axis)
center_z       = -18.643                      # Define center for search space (z-axis)
size_x         = 18.75                        # Define the length of the search space box (x-axis)
size_y         = 17.25                        # Define the length of the search space box (y-axis)
size_z         = 18.75                        # Define the length of the search space box (z-axis)

if program_choice == 'LigandFit': 
    results = run_ligand_fit(receptor, smi, center_x, center_y, center_z)
    sys.exit()
        
if program_choice == 'AutodockZN': 
    results = run_AutodockZN(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness)
    sys.exit()

if program_choice == 'GalaxyDock3': 
    results = run_GalaxyDock3(receptor, smi, center_x, center_y, center_z, exhaustiveness)
    sys.exit()

if program_choice == 'FRED': 
    results = run_fred_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness)
    sys.exit()

if program_choice == 'gold': 
    results = perform_gold_docking(receptor, smi, size_x, size_y, size_z, center_x, center_y, center_z)
    sys.exit()

if program_choice == 'glide':
    results = run_glide_docking(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi)
    sys.exit()
if program_choice == 'rosetta-ligand':
    results = run_rosetta_docking(receptor, smi, center_x, center_y, center_z, exhaustiveness)
    sys.exit()

results = {}  # Storage for results




# Check if receptor file exists inside the config directory: 
files_ls = os.listdir('./config/')
if receptor.split('/')[-1] not in files_ls: 
    raise Exception('Receptor file not found inside the config directory. Please makes sure that the designated receptor file ({}) exists inside config'.format(receptor.split('/')[-1]))

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
    results = run_rDock(receptor, smi)
    sys.exit()
if program_choice == 'ledock': 
    run_leDock(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z)
    sys.exit()
if program_choice == 'adfr': 
    results = run_adfr_docking(receptor, smi)
    sys.exit()



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
    
    if program_choice != 'idock' and program_choice != 'autodock_vina' and program_choice != 'AutodockVina_1.2': 
        cmd = cmd + ['--log', './outputs/log_{}.txt'.format(lig_.split('.')[0])]


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
                    vr_1 = float(line_split[0])
                    vr_2 = float(line_split[1])
                    vr_3 = float(line_split[2])
                    vr_4 = float(line_split[3])
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
        results[lig_] = [docking_score, './outputs/pose_{}.pdb'.format(lig_.split('.')[0]), './outputs/log_{}.txt'.format(lig_.split('.')[0])]
    else: 
        results[lig_] = 'Extremely high pose energy encountered/Unsuccessfull execution.'
    
    if 'VC_log.txt' in os.listdir(): 
        os.system('rm VC_log.txt')


