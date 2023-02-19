#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 22:41:22 2023

@author: akshat
"""
import os 
import subprocess

def run_nnscore2(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/ligand.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    # Perform the calculation: 
    vina_loc = os.getcwd() + '/executables/vina'
    os.system('export VINA_EXEC={}; python ./config/NNScore2.py -receptor {} -ligand {} -vina_executable $VINA_EXEC > output.txt'.format(vina_loc, receptor, lig_path))

    os.system('cp ./output.txt ./outputs/{}.txt'.format(lig_path.split('/')[-1].split('.')[0]))
    os.system('rm output.txt')
    return 

def run_rf_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/ligand.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    # Perform the calculation: 
    os.system('./executables/rf-score-vs --receptor {} {} -O ./outputs/ligands_rescored.pdbqt'.format(receptor, lig_path))
    os.system('./executables/rf-score-vs --receptor {} {} -ocsv > temp.csv'.format(receptor, lig_path))

    with open('./temp.csv', 'r') as f: 
        lines = f.readlines()
    rf_scores = []
    for item in lines[1: ]: 
        rf_scores.append( float(item.split(',')[-1]) )
    
    os.system('rm temp.csv')
    return ['./outputs/ligands_rescored.pdbqt', rf_scores]

def run_smina_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/pose_ligand_1.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_ad4_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/pose_ligand_1.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'ad4_scoring']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_vinandro_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/pose_ligand_1.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'vinardo']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_vina_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/pose_ligand_1.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'vina']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_gnina_scoring(receptor): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
        
    lig_path = './config/pose_ligand_1.pdbqt'
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        raise Exception('Ligand needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
    if os.path.exists('./executables/gnina') == False: 
        raise Exception('Gnina executable {} not found.'.format('./executables/gnina'))

    cmd = ['./executables/gnina', '--receptor', receptor, '-l', lig_path, '--score_only']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    score = command_out[-8:-4]
    key_ = {}
    for item in score: 
        A = item.split(':')
        key_[A[0]] = float(A[1].split(' ')[1])

    return key_

def run_PLANTS_chemplp_scoring(receptor, ligand_file):
    
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        raise Exception('Ligand needs to be in mol2 format. Please try again, after incorporating this correction.')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		chemplp\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))
    
    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    
    return float(command_out.split(' ')[-1])

def run_PLANTS_plp_scoring(receptor, ligand_file):
    
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        raise Exception('Ligand needs to be in mol2 format. Please try again, after incorporating this correction.')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		plp\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))

    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    
    return float(command_out.split(' ')[-1])

def run_PLANTS_plp95_scoring(receptor, ligand_file):
    
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        raise Exception('Ligand needs to be in mol2 format. Please try again, after incorporating this correction.')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		plp95\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))

    
    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    
    return float(command_out.split(' ')[-1])


def contact_score(receptor_file, chimera_path, dock6_path, ligand_file): 

    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    
    with open('./dock6_sript.sh', 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(chimera_path)])
        f.writelines(['export DOCK6={}\n'.format(dock6_path)])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(receptor_file)])
        f.writelines(['$DOCK6/bin/sphgen INSPH\n']) # Del: rm OUTSPH rec.sph temp1.ms temp3.atc
        f.writelines(['$DOCK6/bin/sphere_selector rec.sph {} 12.0\n'.format(ligand_file)])
        f.writelines(['$DOCK6/bin/showbox < box.in\n'])
        f.writelines(['$DOCK6/bin/grid -i grid.in\n'])
        f.writelines(['$DOCK6/bin/dock6 -i Contact_Score.in\n'])
    
    os.system('chmod 777 dock6_sript.sh')
    os.system('cp config/dockprep.py ./dockprep.py')
    
    # Create INSPH File: 
    with open('./INSPH', 'w') as f: 
        f.writelines('rec.ms\n')
        f.writelines('R\n')
        f.writelines('X\n')
        f.writelines('0.0\n')
        f.writelines('4.0\n')
        f.writelines('1.4\n')
        f.writelines('rec.sph\n')
    
    # Create box.in File: 
    with open('./box.in', 'w') as f: 
        f.writelines('N\n')
        f.writelines('U\n')
        f.writelines('{}   {}    {}\n'.format(9, 14, 7))
        f.writelines('25 25 25\n')
        f.writelines('rec_box.pdb\n')
    
    with open('./grid.in', 'w') as f: 
        f.writelines('compute_grids                  yes\n')
        f.writelines('energy_score                   yes\n')
        f.writelines('energy_cutoff_distance         9999\n')
        f.writelines('atom_model                     a\n')
        f.writelines('bump_filter                    yes\n')
        f.writelines('receptor_file                  {}\n'.format(receptor_file))
        f.writelines('box_file                       rec_box.pdb\n')
        f.writelines('vdw_definition_file            {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path))
        f.writelines('score_grid_prefix              grid\n')
        f.writelines('grid_spacing                   0.3\n')
        f.writelines('output_molecule                no\n')
        f.writelines('contact_score                  yes\n')
        f.writelines('attractive_exponent            6\n')
        f.writelines('repulsive_exponent             12\n')
        f.writelines('distance_dielectric            yes\n')
        f.writelines('dielectric_factor              4\n')
        f.writelines('bump_overlap                   0.75\n')
        f.writelines('contact_cutoff_distance        4.5\n')
    
    with open('./Contact_Score.in', 'w') as f: 
        f.writelines(['conformer_search_type                                        rigid\n'])
        f.writelines(['use_internal_energy                                          yes\n'])
        f.writelines(['internal_energy_rep_exp                                      12\n'])
        f.writelines(['internal_energy_cutoff                                       100.0\n'])
        f.writelines(['ligand_atom_file                                             {}\n'.format(ligand_file)])
        f.writelines(['limit_max_ligands                                            no\n'])
        f.writelines(['skip_molecule                                                no\n'])
        f.writelines(['read_mol_solvation                                           no\n'])
        f.writelines(['calculate_rmsd                                               no\n'])
        f.writelines(['use_database_filter                                          no\n'])
        f.writelines(['orient_ligand                                                no\n'])
        f.writelines(['bump_filter                                                  no\n'])
        f.writelines(['score_molecules                                              yes\n'])
        f.writelines(['contact_score_primary                                        yes\n'])
        f.writelines(['contact_score_secondary                                      no\n'])
        f.writelines(['contact_score_cutoff_distance                                4.5\n'])
        f.writelines(['contact_score_clash_overlap                                  0.75\n'])
        f.writelines(['contact_score_clash_penalty                                  50\n'])
        f.writelines(['contact_score_grid_prefix                                    grid\n'])
        f.writelines(['grid_score_secondary                                         no\n'])
        f.writelines(['multigrid_score_secondary                                    no\n'])
        f.writelines(['dock3.5_score_secondary                                      no\n'])
        f.writelines(['continuous_score_secondary                                   no\n'])
        f.writelines(['footprint_similarity_score_secondary                         no\n'])
        f.writelines(['pharmacophore_score_secondary                                no\n'])
        f.writelines(['descriptor_score_secondary                                   no\n'])
        f.writelines(['gbsa_zou_score_secondary                                     no\n'])
        f.writelines(['gbsa_hawkins_score_secondary                                 no\n'])
        f.writelines(['SASA_score_secondary                                         no\n'])
        f.writelines(['amber_score_secondary                                        no\n'])
        f.writelines(['minimize_ligand                                              yes\n'])
        f.writelines(['simplex_max_iterations                                       1000\n'])
        f.writelines(['simplex_tors_premin_iterations                               0\n'])
        f.writelines(['simplex_max_cycles                                           1\n'])
        f.writelines(['simplex_score_converge                                       0.1\n'])
        f.writelines(['simplex_cycle_converge                                       1.0\n'])
        f.writelines(['simplex_trans_step                                           1.0\n'])
        f.writelines(['simplex_rot_step                                             0.1\n'])
        f.writelines(['simplex_tors_step                                            10.0\n'])
        f.writelines(['simplex_random_seed                                          0\n'])
        f.writelines(['simplex_restraint_min                                        no\n'])
        f.writelines(['atom_model                                                   all\n'])
        f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path)])
        f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(dock6_path)])
        f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(dock6_path)])
        f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
        f.writelines(['write_orientations                                           no\n'])
        f.writelines(['num_scored_conformers                                        1\n'])
        f.writelines(['rank_ligands                                                 no\n'])
    
    os.system('./dock6_sript.sh')
    os.system('rm box.in Contact_Score.in dock6_sript.sh dockprep.py dockprep.pyc grid.bmp grid.cnt grid.in grid.nrg INSPH OUTSPH rec.ms rec.sph rec_box.pdb rec_charged.mol2 rec_noH.pdb selected_spheres.sph')
    with open('./ligand_out_scored.mol2', 'r') as f: 
        lines = f.readlines()
    score = float([x for x in lines[2].split(' ') if x!= ''][-1])
    
    return score


def continuous_score(receptor_file, chimera_path, dock6_path, ligand_file): 

    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    
    with open('./dock6_sript.sh', 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(chimera_path)])
        f.writelines(['export DOCK6={}\n'.format(dock6_path)])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(receptor_file)])
        f.writelines(['$DOCK6/bin/dock6 -i Continuous_Score.in\n'])
        
    os.system('cp config/dockprep.py ./dockprep.py')
    os.system('chmod 777 dock6_sript.sh')
    
    with open('./Continuous_Score.in', 'w') as f: 
        f.writelines(['conformer_search_type                                        rigid\n'])
        f.writelines(['use_internal_energy                                          yes\n'])
        f.writelines(['internal_energy_rep_exp                                      12\n'])
        f.writelines(['internal_energy_cutoff                                       100.0\n'])
        f.writelines(['ligand_atom_file                                             {}\n'.format(ligand_file)])
        f.writelines(['limit_max_ligands                                            no\n'])
        f.writelines(['skip_molecule                                                no\n'])
        f.writelines(['read_mol_solvation                                           no\n'])
        f.writelines(['calculate_rmsd                                               no\n'])
        f.writelines(['use_database_filter                                          no\n'])
        f.writelines(['orient_ligand                                                no\n'])
        f.writelines(['bump_filter                                                  no\n'])
        f.writelines(['score_molecules                                              yes\n'])
        f.writelines(['contact_score_primary                                        no\n'])
        f.writelines(['contact_score_secondary                                      no\n'])
        f.writelines(['grid_score_primary                                           no\n'])
        f.writelines(['grid_score_secondary                                         no\n'])
        f.writelines(['multigrid_score_primary                                      no\n'])
        f.writelines(['multigrid_score_secondary                                    no\n'])
        f.writelines(['dock3.5_score_primary                                        no\n'])
        f.writelines(['dock3.5_score_secondary                                      no\n'])
        f.writelines(['continuous_score_primary                                     yes\n'])
        f.writelines(['continuous_score_secondary                                   no\n'])
        f.writelines(['cont_score_rec_filename                                      rec_charged.mol2\n'])
        f.writelines(['cont_score_att_exp                                           6\n'])
        f.writelines(['cont_score_rep_exp                                           12\n'])
        f.writelines(['cont_score_rep_rad_scale                                     1.0\n'])
        f.writelines(['cont_score_use_dist_dep_dielectric                           yes\n'])
        f.writelines(['cont_score_dielectric                                        4.0\n'])
        f.writelines(['cont_score_vdw_scale                                         yes\n'])
        f.writelines(['cont_score_turn_off_vdw                                      yes\n'])
        f.writelines(['cont_score_es_scale                                          1.0\n'])
        f.writelines(['footprint_similarity_score_secondary                         no\n'])
        f.writelines(['pharmacophore_score_secondary                                no\n'])
        f.writelines(['descriptor_score_secondary                                   no\n'])
        f.writelines(['gbsa_zou_score_secondary                                     no\n'])
        f.writelines(['gbsa_hawkins_score_secondary                                 no\n'])
        f.writelines(['SASA_score_secondary                                         no\n'])
        f.writelines(['amber_score_secondary                                        no\n'])
        f.writelines(['minimize_ligand                                              yes\n'])
        f.writelines(['simplex_max_iterations                                       1000\n'])
        f.writelines(['simplex_tors_premin_iterations                               0\n'])
        f.writelines(['simplex_max_cycles                                           1\n'])
        f.writelines(['simplex_score_converge                                       0.1\n'])
        f.writelines(['simplex_cycle_converge                                       1.0\n'])
        f.writelines(['simplex_trans_step                                           1.0\n'])
        f.writelines(['simplex_rot_step                                             0.1\n'])
        f.writelines(['simplex_tors_step                                            10.0\n'])
        f.writelines(['simplex_random_seed                                          0\n'])
        f.writelines(['simplex_restraint_min                                        no\n'])
        f.writelines(['atom_model                                                   all\n'])
        f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path)])
        f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(dock6_path)])
        f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(dock6_path)])
        f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
        f.writelines(['write_orientations                                           no\n'])
        f.writelines(['num_scored_conformers                                        1\n'])
        f.writelines(['rank_ligands                                                 no\n'])
    
    os.system('./dock6_sript.sh')
    
    os.system('rm Continuous_Score.in dock6_sript.sh dockprep.py dockprep.pyc rec_charged.mol2 rec_noH.pdb')
    
    with open('./ligand_out_scored.mol2', 'r') as f: 
        lines = f.readlines()
    score = float([x for x in lines[2].split(' ') if x!= ''][-1])
    
    return score 

def grid_score(receptor_file, chimera_path, dock6_path, ligand_file): 

    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    
    with open('./dock6_sript.sh', 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(chimera_path)])
        f.writelines(['export DOCK6={}\n'.format(dock6_path)])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(receptor_file)])
        f.writelines(['$DOCK6/bin/sphgen INSPH\n']) 
        f.writelines(['$DOCK6/bin/sphere_selector rec.sph {} 12.0\n'.format(ligand_file)])
        f.writelines(['$DOCK6/bin/showbox < box.in\n'])
        f.writelines(['$DOCK6/bin/grid -i grid.in\n'])
        f.writelines(['$DOCK6/bin/dock6 -i Grid_Score.in\n'])
    
    os.system('chmod 777 dock6_sript.sh')
    
    # Create INSPH File: 
    with open('./INSPH', 'w') as f: 
        f.writelines('rec.ms\n')
        f.writelines('R\n')
        f.writelines('X\n')
        f.writelines('0.0\n')
        f.writelines('4.0\n')
        f.writelines('1.4\n')
        f.writelines('rec.sph\n')
    
    # Create box.in File: 
    with open('./box.in', 'w') as f: 
        f.writelines('N\n')
        f.writelines('U\n')
        f.writelines('{}   {}    {}\n'.format(9, 14, 7))
        f.writelines('25 25 25\n')
        f.writelines('rec_box.pdb\n')
    
    with open('./grid.in', 'w') as f: 
        f.writelines('compute_grids                  yes\n')
        f.writelines('energy_score                   yes\n')
        f.writelines('energy_cutoff_distance         9999\n')
        f.writelines('atom_model                     a\n')
        f.writelines('bump_filter                    yes\n')
        f.writelines('receptor_file                  {}\n'.format(receptor_file))
        f.writelines('box_file                       rec_box.pdb\n')
        f.writelines('vdw_definition_file            {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path))
        f.writelines('score_grid_prefix              grid\n')
        f.writelines('grid_spacing                   0.3\n')
        f.writelines('output_molecule                no\n')
        f.writelines('contact_score                  yes\n')
        f.writelines('attractive_exponent            6\n')
        f.writelines('repulsive_exponent             12\n')
        f.writelines('distance_dielectric            yes\n')
        f.writelines('dielectric_factor              4\n')
        f.writelines('bump_overlap                   0.75\n')
        f.writelines('contact_cutoff_distance        4.5\n')
        
        
    with open('./Grid_Score.in', 'w') as f: 
        f.writelines('conformer_search_type                                        rigid\n')
        f.writelines('use_internal_energy                                          yes\n')
        f.writelines('internal_energy_rep_exp                                      12\n')
        f.writelines('internal_energy_cutoff                                       100.0\n')
        f.writelines('ligand_atom_file                                             {}\n'.format(ligand_file))
        f.writelines('limit_max_ligands                                            no\n')
        f.writelines('skip_molecule                                                no\n')
        f.writelines('read_mol_solvation                                           no\n')
        f.writelines('calculate_rmsd                                               no\n')
        f.writelines('use_database_filter                                          no\n')
        f.writelines('orient_ligand                                                no\n')
        f.writelines('bump_filter                                                  no\n')
        f.writelines('score_molecules                                              yes\n')
        f.writelines('contact_score_primary                                        no\n')
        f.writelines('contact_score_secondary                                      no\n')
        f.writelines('grid_score_primary                                           yes\n')
        f.writelines('grid_score_secondary                                         no\n')
        f.writelines('grid_score_rep_rad_scale                                     1.0\n')
        f.writelines('grid_score_vdw_scale                                         1\n')
        f.writelines('grid_score_es_scale                                          1\n')
        f.writelines('grid_score_grid_prefix                                       grid\n')
        f.writelines('multigrid_score_secondary                                    no\n')
        f.writelines('dock3.5_score_secondary                                      no\n')
        f.writelines('continuous_score_secondary                                   no\n')
        f.writelines('footprint_similarity_score_secondary                         no\n')
        f.writelines('pharmacophore_score_secondary                                no\n')
        f.writelines('descriptor_score_secondary                                   no\n')
        f.writelines('gbsa_zou_score_secondary                                     no\n')
        f.writelines('gbsa_hawkins_score_secondary                                 no\n')
        f.writelines('SASA_score_secondary                                         no\n')
        f.writelines('amber_score_secondary                                        no\n')
        f.writelines('minimize_ligand                                              yes\n')
        f.writelines('simplex_max_iterations                                       1000\n')
        f.writelines('simplex_tors_premin_iterations                               0\n')
        f.writelines('simplex_max_cycles                                           1\n')
        f.writelines('simplex_score_converge                                       0.1\n')
        f.writelines('simplex_cycle_converge                                       1.0\n')
        f.writelines('simplex_trans_step                                           1.0\n')
        f.writelines('simplex_rot_step                                             0.1\n')
        f.writelines('simplex_tors_step                                            10.0\n')
        f.writelines('simplex_random_seed                                          0\n')
        f.writelines('simplex_restraint_min                                        no\n')
        f.writelines('atom_model                                                   all\n')
        f.writelines('vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path))
        f.writelines('flex_defn_file                                               {}/parameters/flex.defn\n'.format(dock6_path))
        f.writelines('flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(dock6_path))
        f.writelines('ligand_outfile_prefix                                        ligand_out\n')
        f.writelines('write_orientations                                           no\n')
        f.writelines('num_scored_conformers                                        1\n')
        f.writelines('rank_ligands                                                 no\n')
    
        os.system('./dock6_sript.sh')
        
        os.system('rm box.in Contact_Score.in dock6_sript.sh dockprep.py dockprep.pyc grid.bmp grid.cnt grid.in grid.nrg INSPH OUTSPH rec.ms rec.sph rec_box.pdb rec_charged.mol2 rec_noH.pdb selected_spheres.sph')
        with open('./ligand_out_scored.mol2', 'r') as f: 
            lines = f.readlines()
        score = float([x for x in lines[2].split(' ') if x!= ''][-1])
        
        return score 

def run_mm_gbsa(): 
    chimera_path  = '/home/akshat/chimera' # Please update the Chimera path 
    ligand_file   = 'ligand.mol2' 
    receptor_file = 'receptor.pdb'
    
    output = {}
    
    if os.path.exists(chimera_path) == False: 
        raise Exception('Location of Chemira not found (used location from variable chimera_path) when trying to initiate mm_fbsa calculation.')
    
    # Check if paths exists'
    if os.path.exists(ligand_file) == False: 
        raise Exception('Ligand file ligand.mol2 not found. Please add the file to the current working directory')
    if os.path.exists(receptor_file) == False: 
        raise Exception('Receptor file receptor.pdb not found. Please add the file to the current working directory')

    # Check to make sure ligand is in mol2 format: 
    lig_format = ligand_file.split('.')[1]
    if lig_format != 'mol2': 
        raise Exception('Please ensure ligand is in mol2 file')

    with open('./GBSA.sh', 'w') as f: 
        # Getting Ligand Parameters: 
        f.writelines('export Chimera={}\n'.format(chimera_path))
        f.writelines('charge=`$Chimera/bin/chimera --nogui --silent ligand.mol2 ./config/charges.py`\n')
        f.writelines('antechamber -i ligand.mol2 -fi mol2 -o ligand_bcc.mol2 -fo mol2 -at gaff2 -c gas -rn LIG -nc $charge -pf y\n')
        f.writelines('parmchk2 -i ligand_bcc.mol2 -f mol2 -o ligand.frcmod\n')

        # Building Topology Files:
        f.writelines('tleap -f ./config/tleap_r.in\n')
        f.writelines('tleap -f ./config/tleap_c.in\n')
        
        # Run MD: 
        f.writelines('sander -O -i ./config/min.in -p complex.prmtop -c complex.inpcrd -r min.rst -ref complex.inpcrd -o minim.out\n')
        
        # Running MMPBSA.py
        f.writelines('MMPBSA.py -O -i ./config/gbsa.in -cp complex.prmtop -rp receptor.prmtop -lp ligand.prmtop -y  min.rst\n')

    os.system('chmod 777 GBSA.sh')
    os.system('./GBSA.sh') # Run the calculation
    
    # Remove auxillary files: 
    os.system('rm complex.inpcrd complex.prmtop leap.log ligand.frcmod ligand.inpcrd ligand.prmtop ligand_bcc.mol2 mdcrd mdinfo min.rst minim.out receptor.inpcrd receptor.prmtop reference.frc')
    
    # Read in the result: 
    try: 
        with open('./FINAL_RESULTS_MMPBSA.dat', 'r') as f: 
            result = f.readlines()
        result = [x for x in result if 'DELTA TOTAL' in x][0]
        result = float([x for x in result.split(' ') if x != ''][2])
        
        output[ligand_file] = result
    except: 
        output[ligand_file] = 'FAIL'
    
    os.system('rm FINAL_RESULTS_MMPBSA.dat')
    
    return output