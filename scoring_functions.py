#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 22:41:22 2023

@author: akshat
"""
import os 
import subprocess
import shutil


def convert_ligand_format(ligand_, new_format): 
    """Converts a ligand file to a different file format using the Open Babel tool.

        Args:
            ligand_ (str): The path to the input ligand file.
            new_format (str): The desired output format for the ligand file.
    
        Returns:
            None
    
        Raises:
            Exception: If the input file does not exist, or if the Open Babel tool is not installed.
    
        Examples:
            To convert a ligand file from mol2 format to pdbqt format:
            >>> convert_ligand_format('./ligands/ligand1.mol2', 'pdbqt')
    """
    input_format = ligand_.split('.')[-1]
    os.system('obabel {} -O {}'.format(ligand_, ligand_.replace(input_format, new_format)))

def run_nnscore2(receptor, lig_path): 
    """
    Perform scoring for docking ligands using the NNScore2 method.
    
    Args:
    - receptor: string, path to the receptor file in pdbqt format.
    - lig_path: string, path to the ligand file in pdbqt format.
    
    Returns:
    - scores: list of strings, the best score for the docking run.
    
    Raises:
    - Exception: If receptor is not in pdbqt format or not found.
    - Exception: If lig_path is not in pdbqt format or not found.
    
    The function uses the Vina executable and NNScore2.py script in the executables directory to perform the docking run. The output of the run is saved in the output.txt file, which is then read to extract the best score.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

    # Perform the calculation: 
    vina_loc = os.getcwd() + '/executables/vina'
    os.system('export VINA_EXEC={}; python ./config/NNScore2.py -receptor {} -ligand {} -vina_executable $VINA_EXEC > output.txt'.format(vina_loc, receptor, lig_path))

    # os.system('cp ./output.txt ./outputs/{}.txt'.format(lig_path.split('/')[-1].split('.')[0]))
    with open('./output.txt', 'r') as f: 
        lines = f.readlines()
    scores = [x for x in lines if 'Best Score:' in x]
    scores = [A.split('(')[-1].split(')')[0] for A in scores]
    
    os.system('rm output.txt')
    return scores

def run_rf_scoring(receptor, lig_path): 
    """
    Runs RF-Score calculation for given receptor-ligand complex.

    Parameters:
    -----------
    receptor : str
        File path of the receptor in PDB format.
    lig_path : str
        File path of the ligand in PDB or PDBQT format.

    Returns:
    --------
    list
        A list containing the file path of the rescored ligand in PDBQT format and a list of RF-Score scores for the ligand.

    Raises:
    -------
    Exception
        If the receptor file is not in PDB format or not found.
        If the ligand file is not found.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Reception path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

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

def run_smina_scoring(receptor, lig_path): 
    """
    Runs smina scoring on a receptor-ligand complex and returns the score.

    Args:
    receptor (str): path to the receptor PDBQT file
    lig_path (str): path to the ligand PDBQT file

    Returns:
    float: the smina score for the complex
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))        
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_ad4_scoring(receptor, ligand_path): 
    """
    Runs AutoDock4 scoring on a receptor-ligand complex and returns the score.

    Args:
    receptor (str): path to the receptor PDBQT file
    lig_path (str): path to the ligand PDBQT file
    
    Returns: 
    float: binding affinity (in kcal/mol) as a floating point number
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(ligand_path) == False: 
        raise Exception('Ligand path {} not found.'.format(ligand_path))
    
    lig_format = ligand_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_path, 'pdbqt')
        ligand_path = ligand_path.replace(lig_format, 'pdbqt')
        
    cmd = ['./executables/smina', '--receptor', receptor, '-l', ligand_path, '--score_only', '--scoring', 'ad4_scoring']    
    
    command_run = subprocess.run(cmd, capture_output=True)

    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_vinandro_scoring(receptor, lig_path): 
    """
    Runs Vinardo scoring on the given receptor-ligand pair and returns the binding affinity score.

    Parameters:
    -----------
    receptor: str
        Path to the receptor file in pdbqt format.
    lig_path: str
        Path to the ligand file in pdbqt format.
    
    Returns:
    --------
    float
        Vinardo binding affinity score for the receptor-ligand pair.

    Raises:
    -------
    Exception
        If receptor file format is not pdbqt.
    Exception
        If receptor file is not found in the provided path.
    Exception
        If ligand file is not found in the provided path.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
   
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'vinardo']    
    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_vina_scoring(receptor, lig_path): 
    """
    Runs Vina scoring on the given receptor-ligand pair and returns the binding affinity score.

    Parameters:
    -----------
    receptor: str
        Path to the receptor file in pdbqt format.
    lig_path: str
        Path to the ligand file in pdbqt format.
    
    Returns:
    --------
    float
        Vina binding affinity score for the receptor-ligand pair.

    Raises:
    -------
    Exception
        If receptor file format is not pdbqt.
    Exception
        If receptor file is not found in the provided path.
    Exception
        If ligand file is not found in the provided path.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

    cmd = ['./executables/smina', '--receptor', receptor, '-l', lig_path, '--score_only', '--scoring', 'vina']    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = float([x for x in command_out if 'Affinity' in x][0].split(' ')[1])

    return command_out

def run_gnina_scoring(receptor, lig_path): 
    """
    Runs Gnina scoring on the given receptor-ligand pair and returns the binding affinity score.

    Parameters:
    -----------
    receptor: str
        Path to the receptor file in pdbqt format.
    lig_path: str
        Path to the ligand file in pdbqt format.
    
    Returns:
    --------
    float
        Gnina binding affinity score for the receptor-ligand pair.

    Raises:
    -------
    Exception
        If receptor file format is not pdbqt.
    Exception
        If receptor file is not found in the provided path.
    Exception
        If ligand file is not found in the provided path.
    Exception
        If Gnina executable is not found in /executables/.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    if os.path.exists('./executables/gnina') == False: 
        raise Exception('Gnina executable {} not found.'.format('./executables/gnina'))
    if os.path.exists(lig_path) == False: 
        raise Exception('Ligand path {} not found.'.format(lig_path))
        
    lig_format = lig_path.split('.')[-1]
    if lig_format != 'pdbqt': 
        print('Ligand needs to be in pdbqt format. Converting ligand format using obabel.')
        convert_ligand_format(lig_path, 'pdbqt')
        lig_path = lig_path.replace(lig_format, 'pdbqt')

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
    """
    Runs the PLANTS molecular docking program with the ChemPLP scoring function to predict
    the binding affinity of a ligand to a receptor. 
    
    Args:
    - receptor (str): path to the receptor file in mol2 format.
    - ligand_file (str): path to the ligand file in pdb or mol2 format.
    
    Returns:
    - score (float): predicted binding affinity score.
    """
    
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		chemplp\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))
    
    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    os.system('rm plants_config')
    return float(command_out.split(' ')[-1])

def run_PLANTS_plp_scoring(receptor, ligand_file):
    """
    Runs the PLANTS molecular docking program with the PLP scoring function to predict
    the binding affinity of a ligand to a receptor. 
    
    Args:
    - receptor (str): path to the receptor file in mol2 format.
    - ligand_file (str): path to the ligand file in pdb or mol2 format.
    
    Returns:
    - score (float): predicted binding affinity score.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		plp\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))

    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    os.system('rm plants_config')
    return float(command_out.split(' ')[-1])

def run_PLANTS_plp95_scoring(receptor, ligand_file):
    """
    Runs the PLANTS molecular docking program with the PLP95 scoring function to predict
    the binding affinity of a ligand to a receptor. 
    
    Args:
    - receptor (str): path to the receptor file in mol2 format.
    - ligand_file (str): path to the ligand file in pdb or mol2 format.
    
    Returns:
    - score (float): predicted binding affinity score.
    """
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')
    
    with open('plants_config', 'w') as f: 
        f.writelines('scoring_function 		plp95\n')
        f.writelines('protein_file 			{}\n'.format(receptor))
        f.writelines('ligand_file 			{}\n'.format(ligand_file))

    cmd = ['./executables/PLANTS', '--mode', 'rescore', './plants_config']    
    command_run = subprocess.run(cmd, capture_output=True)
    command_out = command_run.stdout.decode("utf-8").split('\n')
    command_out = [x for x in command_out if 'best score:' in x][-1]    
    os.system('rm plants_config')
    return float(command_out.split(' ')[-1])


def contact_score(receptor_file, chimera_path, dock6_path, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z): 
    """
    Calculate the Contact Score between a receptor and a ligand using DOCK6.

    Args:
        - receptor_file (str): Path to the receptor file (in PDB format).
        - chimera_path (str): Path to UCSF Chimera executable.
        - dock6_path (str): Path to DOCK6 program.
        - ligand_file (str): Path to the ligand file (in MOL2 format).
        - center_x (float): x-coordinate of the center of the box.
        - center_y (float): y-coordinate of the center of the box.
        - center_z (float): z-coordinate of the center of the box.
        - size_x (float): Size of the box in the x direction.
        - size_y (float): Size of the box in the y direction.
        - size_z (float): Size of the box in the z direction.

    Returns:
        - score (float): Contact Score between the receptor and the ligand.
    """

    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')

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
        f.writelines('{}   {}    {}\n'.format(center_x, center_y, center_z))
        f.writelines('{} {} {}\n'.format(size_x, size_y, size_z))
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
    """
    This function performs docking of a given ligand to a receptor using DOCK6 and returns the continuous score of the docked complex.

    Args:
        - receptor_file (str): path to the receptor file in PDB format.
        - chimera_path (str): path to the Chimera installation directory.
        - dock6_path (str): path to the DOCK6 installation directory.
        - ligand_file (str): path to the ligand file in MOL2 format.

    Returns:
        - score (float): the continuous score of the docked complex.

    Raises:
        - Exception: if the receptor file is not in PDB format.
    """

    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')

    
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

def grid_score(receptor_file, chimera_path, dock6_path, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z): 
    """
    This function calculates the grid score of a given ligand on the receptor protein. The function performs the following steps:
        1. Checks if the receptor and ligand files are in the correct format. If not, it converts the ligand file to the correct format.
        2. Generates a box around the receptor protein.
        3. Generates spheres around the ligand and selects spheres that are within 12 Angstrom of the ligand.
        4. Calculates the grid scores of the ligand in the box.
    
    Args:
        - receptor_file (str): The path to the receptor file in PDB format.
        - chimera_path (str): The path to the Chimera installation directory.
        - dock6_path (str): The path to the Dock6 installation directory.
        - ligand_file (str): The path to the ligand file in MOL2 format.
        - center_x (float): The x-coordinate of the center of the grid box.
        - center_y (float): The y-coordinate of the center of the grid box.
        - center_z (float): The z-coordinate of the center of the grid box.
        - size_x (float): The size of the grid box in the x-dimension.
        - size_y (float): The size of the grid box in the y-dimension.
        - size_z (float): The size of the grid box in the z-dimension.
        
    Returns:
        - score (float): The grid score of the ligand on the receptor protein.
    """
    
    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')
    
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
        f.writelines('{}   {}    {}\n'.format(center_x, center_y, center_z))
        f.writelines('{} {} {}\n'.format(size_x, size_y, size_z))
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

def run_mm_gbsa(chimera_path, ligand_file, receptor_file): 
    
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

def Zou_GBSA(receptor_file, chimera_path, dock6_path, ligand_file, center_x, center_y, center_z, size_x, size_y, size_z):
    """
    Args:
    - receptor_file (str): The path to the receptor file in PDB format.
    - chimera_path (str): The path to the Chimera installation directory.
    - dock6_path (str): The path to the Dock6 installation directory.
    - ligand_file (str): The path to the ligand file in MOL2 format.
    - center_x (float): The x-coordinate of the center of the grid box.
    - center_y (float): The y-coordinate of the center of the grid box.
    - center_z (float): The z-coordinate of the center of the grid box.
    - size_x (float): The size of the grid box in the x-dimension.
    - size_y (float): The size of the grid box in the y-dimension.
    - size_z (float): The size of the grid box in the z-dimension.
    """
    recetor_format = receptor_file.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    lig_format = ligand_file.split('.')[-1]
    if lig_format != 'mol2': 
        print('Ligand needs to be in mol2 format. Converting ligand format using obabel.')
        convert_ligand_format(ligand_file, 'mol2')
        ligand_file = ligand_file.replace(lig_format, 'mol2')


    with open('./GBSA.sh', 'w') as f: 
        f.writelines(['export Chimera={}\n'.format(chimera_path)])
        f.writelines(['export DOCK6={}\n'.format(dock6_path)])
        f.writelines(['$Chimera/bin/chimera --nogui {} dockprep.py\n'.format(receptor_file)])
        f.writelines(['$DOCK6/bin/sphgen INSPH\n']) 
        f.writelines(['$DOCK6/bin/sphere_selector rec.sph {} 12.0 \n'.format(ligand_file)])
        f.writelines(['$DOCK6/bin/showbox < box.in\n'])
        f.writelines(['$DOCK6/bin/grid -i grid.in\n'])
        f.writelines(['cd nchemgrid_GB\n'])
        f.writelines(['$DOCK6/bin/nchemgrid_GB\n'])
        f.writelines(['cd ../nchemgrid_SA\n'])
        f.writelines(['$DOCK6/bin/nchemgrid_SA\n'])

        f.writelines(['$DOCK6/bin/dock6 -i Hawkins_GBSA_Score.in\n'])


    os.system('cp config/dockprep.py ./dockprep.py')
    os.system('chmod 777 GBSA.sh')  

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
        f.writelines('{}   {}    {}\n'.format(center_x, center_y, center_z))
        f.writelines('{} {} {}\n'.format(size_x, size_y, size_z))
        f.writelines('rec_box.pdb\n')

    with open('./grid.in', 'w') as f: 
        f.writelines('compute_grids                  yes\n')
        f.writelines('grid_spacing                   0.3\n')
        f.writelines('output_molecule                no\n')
        f.writelines('contact_score                  no\n')
        f.writelines('energy_score                   yes\n')
        f.writelines('energy_cutoff_distance         9999\n')
        f.writelines('atom_model                     a\n')
        f.writelines('attractive_exponent            6\n')
        f.writelines('repulsive_exponent             12\n')
        f.writelines('distance_dielectric            no\n')
        f.writelines('dielectric_factor              1\n')
        f.writelines('bump_filter                    yes\n')
        f.writelines('bump_overlap                   0.75\n')
        f.writelines('receptor_file                  {}\n'.format(receptor_file))
        f.writelines('box_file                       rec_box.pdb\n')
        f.writelines('vdw_definition_file            {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path))
        f.writelines('score_grid_prefix              solvent_grid\n')

    with open('./Hawkins_GBSA_Score.in', 'w') as f: 
        f.writelines('conformer_search_type                                        rigid\n')
        f.writelines('use_internal_energy                                          no\n')
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
        f.writelines('grid_score_primary                                           no\n')
        f.writelines('grid_score_secondary                                         no\n')
        f.writelines('multigrid_score_primary                                      no\n')
        f.writelines('multigrid_score_secondary                                    no\n')
        f.writelines('dock3.5_score_primary                                        no\n')
        f.writelines('dock3.5_score_secondary                                      no\n')
        f.writelines('continuous_score_primary                                     no\n')
        f.writelines('continuous_score_secondary                                   no\n')
        f.writelines('footprint_similarity_score_primary                           no\n')
        f.writelines('footprint_similarity_score_secondary                         no\n')
        f.writelines('pharmacophore_score_primary                                  no\n')
        f.writelines('pharmacophore_score_secondary                                no\n')
        f.writelines('descriptor_score_primary                                     no\n')
        f.writelines('descriptor_score_secondary                                   no\n')
        f.writelines('gbsa_zou_score_primary                                       no\n')
        f.writelines('gbsa_zou_score_secondary                                     no\n')
        f.writelines('gbsa_hawkins_score_primary                                   yes\n')
        f.writelines('gbsa_hawkins_score_secondary                                 no\n')
        f.writelines('gbsa_hawkins_score_rec_filename                              rec_charged.mol2')
        f.writelines('gbsa_hawkins_score_solvent_dielectric                        78.5\n')
        f.writelines('gbsa_hawkins_use_salt_screen                                 no\n')
        f.writelines('gbsa_hawkins_score_gb_offset                                 0.09\n')
        f.writelines('gbsa_hawkins_score_cont_vdw_and_es                           no\n')
        f.writelines('gbsa_hawkins_score_grid_prefix                               solvent_grid\n')
        f.writelines('SASA_score_secondary                                         no\n')
        f.writelines('amber_score_secondary                                        no\n')
        f.writelines('minimize_ligand                                              no\n')
        f.writelines('atom_model                                                   all\n')
        f.writelines('vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path))
        f.writelines('flex_defn_file                                               {}/parameters/flex.defn'.format(dock6_path))
        f.writelines('flex_drive_file                                              {}/parameters/flex_drive.tbl'.format(dock6_path))
        f.writelines('ligand_outfile_prefix                                        gbsa_hawkins\n')
        f.writelines('write_orientations                                           no\n')
        f.writelines('num_scored_conformers                                        1\n')
        f.writelines('rank_ligands                                                 no\n')

    os.system('./GBSA.sh')

    # Check that this OS system remove removes everything correctly
    os.system('rm box.in grid.in INSPH OUTSPH rec_box.pdb dockprep.py dockprep.pyc rec_charged.mol2 solvent_grid.bmp solvent_grid.nrg rec.ms rec.sph rec.pdb rec_noH.pdb selected_spheres.sph')
    # need to import shutil function
    shutil.rmtree('nchemgrid_GB')
    shutil.rmtree('nchemgrid_SA')
    with open('./gbsa_hawkins_scored.mol2', 'r') as f: 
        lines = f.readlines()
    score = float([x for x in lines[2].split(' ') if x!= ''][-1])
    
    return score