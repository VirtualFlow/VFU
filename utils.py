#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 20:01:20 2022

@author: akshat
"""
import os 
import subprocess
from lig_process import process_ligand

def run_plants_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z): 
    
    results = {}
    print('Note: I am defaulting to setting cluster_structures to 10; please change me in function run_plants_docking, line 40')
    print('Note: I am defaulting to setting cluster_rmsd to 2; please change me in function run_plants_docking, line 41')
    
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
        
    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')

    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        
        # Prepare the config file: 
        with open('./config.txt', 'w') as f: 
            f.writelines('protein_file {}\n'.format(receptor))
            f.writelines('ligand_file {}\n'.format(lig_path))
            f.writelines('output_dir ./result_{}\n'.format(lig_.split('.')[0])) 
            f.writelines('write_multi_mol2 0\n') 
            f.writelines('bindingsite_center {} {} {}\n'.format(center_x, center_y, center_z))
            f.writelines('bindingsite_radius {}\n'.format( max([size_x, size_y, size_z]) )) 
            f.writelines('cluster_structures {}\n'.format( 10 ) )
            f.writelines('cluster_rmsd {}\n'.format( 2.0 ) ) 
            
        # Run the command:         
        plants_cmd = ['./executables/PLANTS', '--mode', 'screen', 'config.txt']
        plants_cmd = subprocess.run(plants_cmd, capture_output=True)
        plants_cmd = plants_cmd.stdout.decode("utf-8").split('\n')[-6]
        score_     = float(plants_cmd.split(' ')[-1])
        
        # Copy the direcory in the right place: 
        os.system('cp -a {} ./outputs/'.format('./result_{}'.format(lig_.split('.')[0]) ))
        os.system('rm -rf result_{}'.format(lig_.split('.')[0]))
        os.system('rm config.txt')
        
        results[lig_] = [score_, './results/result_{}'.format(lig_.split('.')[0])]

    return results        

def run_autodock_gpu_docking(receptor, smi, program_choice): 

    print('Note: For use of vina gpu, the receptor needs to be prepared in a specif way. Have a look at the examples provided in https://github.com/ccsb-scripps/AutoDock-GPU & the example dir we provided within executables/vf_gpu_example.zip')
    command = []
    
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')
    if recetor_format[-1] != 'fld' and recetor_format[-2] != 'maps': 
        raise Exception('Receptor needs to be of file type .maps.fld (example: 1stp_protein.maps.fld). Please try again, after incorporating this correction.')
    
    # check for the existence of the executable: 
    if 'gpu' in program_choice: 
        executable = [x for x in os.listdir('./executables') if 'autodock_gpu' in x][0]
    elif 'cpu' in program_choice: 
        executable = [x for x in os.listdir('./executables') if 'autodock_cpu' in x][0]
    else: 
        raise Exception('Executable must be of format autodock_cpu/gpu')
    
    if len(executable) == 0: 
        raise Exception('Executable not found. Executable needs to have autodock_gpu in name (example: autodock_gpu_1wi/autodock_cpu_1wi)')
                        
    # Assign the right program for docking:  
    command.append('./executables/{}'.format(program_choice))
    
    # Prepare the ligands: 
    process_ligand(smi, 'pdbqt')
    lig_locations = os.listdir('./ligands/')

    # Get ready for running the files: 
        
    results = {}
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        
        vina_gpu_cmd = command + ['--ffile', '{}'.format(receptor)]
        vina_gpu_cmd = vina_gpu_cmd + ['--lfile', '{}'.format(lig_path)]

        vina_gpu_cmd = subprocess.run(vina_gpu_cmd, capture_output=True)
        vina_gpu_cmd = vina_gpu_cmd.stdout.decode("utf-8").split('\n')[-6]
            
        if vina_gpu_cmd[-1] != 'All jobs ran without errors.\n': 
            print('An error was encountered when executing docking for: ', lig_path)
            results[lig_] = ['FAIL', vina_gpu_cmd]
        else: 
            lines = [x.strip() for x in vina_gpu_cmd if 'best energy' in x][0]
            docking_score = float(lines.split(',')[1].split(' ')[-2])
            results[lig_] = [docking_score, vina_gpu_cmd]
    return results



def run_EquiBind(receptor, smi): 
    files_ls = os.listdir('./')
    if not('data' in files_ls and 'inference.py' in files_ls):
        raise Exception('Please make sure process EquiBind based on the instructions provided in the readme. I could not find the key EquiBind files.')
    
    if receptor.split('.')[2] != 'pdb': 
        raise Exception('For EquiBind, protein file needs to be in pdb file type. Please incorporate this correction')
    
    # Process the ligands
    process_ligand(smi, 'sdf') # mol2 ligand format is supported in plants

    # Make a direcotry containing all the tasks to be performed: 
    os.sytem('mkdir ./data/to_predict')
    
    results = {}
    lig_locations = os.listdir('./ligands/')
    for i,lig_ in enumerate(lig_locations): 
        lig_path = 'ligands/{}'.format(lig_)
        os.system('./data/to_predict/test{}'.format(i))
        os.system('cp {} {}'.format(receptor, './data/to_predict/test{}/rec.pdb'.format(i)))    # Copy the protein file: 
        os.system('cp {} {}'.format(lig_path, './data/to_predict/test{}/ligand.sdf'.format(i))) # Copy the ligand file: 
        results[lig_path] = './data/results/output/test{}'.format(i)
            
    os.system('python inference.py --config=configs_clean/inference.yml')
    print('Results are saved in: ./data/results/output')

    return results


def run_rDock(receptor, smi): 
    
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    
    
    # Create ligands as '.sd' file type: 
    process_ligand(smi, 'sd') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')
    
    ref_lig = '' # TODO!!!
    raise Exception('Note: a reference ligand file needs to be filled in the line above. Please do so and comment this line!')
    
    # Creation of the prm file: 
    print('Please have a look at the prm parameters. Inside [TODO]; we have assigned some default values.')
    with open('config.prm', 'a+') as f: 
        f.writelines('RBT_PARAMETER_FILE_V1.00')
        f.writelines('TITLE gart_DUD')
        f.writelines('RECEPTOR_FILE {}'.format(receptor))
        f.writelines('SECTION MAPPER')
        f.writelines('    SITE_MAPPER RbtLigandSiteMapper')
        f.writelines('    REF_MOL {}'.format(ref_lig))
        f.writelines('    RADIUS 6.0')
        f.writelines('    SMALL_SPHERE 1.0')
        f.writelines('    MIN_VOLUME 100')
        f.writelines('    MAX_CAVITIES 1')
        f.writelines('    VOL_INCR 0.0')
        f.writelines('   GRIDSTEP 0.5')
        f.writelines('END_SECTION')
        f.writelines('SECTION CAVITY')
        f.writelines('    SCORING_FUNCTION RbtCavityGridSF')
        f.writelines('    WEIGHT 1.0')
        f.writelines('END_SECTION')

    # Cavity generation: 
    os.system('rbcavity -was -d -r config.prm')
    
    # Perform docking: 
    os.system('mkdir rDock_outputs')
    for i,lig_ in enumerate(lig_locations): 
        os.system('mkdir rDock_outputs/{}'.format(i))
        lig_path = 'ligands/{}'.format(lig_)
        os.system('rbdock -i {} -o {} -r config.prm -p dock.prm -n 50'.format(lig_path, 'rDock_outputs/{}'.format(i)))
    
    return 







