#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 17:23:20 2023

@author: akshat
"""
from pose_prediction import run_pose_prediction_program


def read_config_file(): 
    """
    Reads a configuration file and extracts the values of certain variables.

    The configuration file should be named 'config.txt' and should be located in the same directory as this function.
    The configuration file should contain a series of key-value pairs, one per line, in the following format:
        key=value
    Lines starting with '#' or containing only whitespace are ignored.

    Returns:
        program_choice (str): The name of the program to use.
        center_x (float): The x-coordinate of the center of the search space.
        center_y (float): The y-coordinate of the center of the search space.
        center_z (float): The z-coordinate of the center of the search space.
        size_x (float): The size of the search space along the x-axis.
        size_y (float): The size of the search space along the y-axis.
        size_z (float): The size of the search space along the z-axis.
        exhaustiveness (int): The exhaustiveness parameter to use.
        smi (str): The SMILES string representing the ligand to dock.
        is_selfies (bool): Whether the input ligand is represented using SELFIES encoding.
        is_peptide (bool): Whether the input ligand is a peptide.
        receptor (str): The path to the receptor file.

    Raises:
        FileNotFoundError: If the configuration file cannot be found.

    Example usage:
    program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor = read_config_file()
    """
    # process the config file: 
    with open('./config.txt', 'r') as f: 
        lines = f.readlines()
    lines = [x for x in lines if x[0] != '#' and x[0] != '\n'] # Remove all lines that correspond to comments
        
    variables = {}
    for string in lines:
        key, value = string.strip().split('=')
        variables[key] = value    
        
    program_choice = variables['program_choice']
    center_x = float(variables['center_x'])
    center_y = float(variables['center_y'])
    center_z = float(variables['center_z'])
    size_x = float(variables['size_x'])
    size_y = float(variables['size_y'])
    size_z = float(variables['size_z'])
    exhaustiveness = int(variables['exhaustiveness'])
    smi = variables['smi']
    is_selfies = variables['is_selfies']
    is_peptide = variables['is_peptide']
    receptor = variables['receptor']
    
    return program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor

if __name__ == "__main__": 
    
    # Read in parameters from config.txt: 
    program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor = read_config_file()
    if '+' in program_choice: 
        program_choice, scoring_function = program_choice.split('+')[0], program_choice.split('+')[1]
        
    # Convert molecules from selfies->smiles/AA_string->smiles if option is specified within config.txt:
    if is_selfies == 'True': 
        import selfies 
        try: 
            smi = selfies.decoder(smi)
        except: 
            raise Exception('Invalid SELFIES provided. Please make sure that the varibale smi contains a valid selfies string. ')
    if is_peptide == 'True': 
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
        raise print('WARNING: Si/Sn/B found. Please be careful when using these for docking. Namely, the original AutoDock programs do not suport these atom types. Please make sure that the selected docking program support these atoms. Please remove this exception if everything is fine :)')


    # Perform pose prediction:
    pose_pred_out = run_pose_prediction_program(program_choice, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, receptor)

    # receptor       = './config/5wiu_test.pdbqt'
    # if program_choice == 'nnscore2': 
    #     run_nnscore2(receptor)
    #     sys.exit()
    # if program_choice == 'rf-score': 
    #     rescored_pose, rf_scores = run_rf_scoring(receptor)
    # if program_choice == 'smina-scoring': 
    #     score = run_smina_scoring(receptor)
    #     sys.exit()
    # if program_choice == 'gnina-scoring': 
    #     score = run_gnina_scoring(receptor)
    #     sys.exit()
    # if program_choice == 'ad4_scoring': 
    #     score = run_ad4_scoring(receptor)
    #     sys.exit()
    # if program_choice == 'vinandro_scoring': 
    #     score = run_vinandro_scoring(receptor)
    #     sys.exit()
    # if program_choice == 'vina_scoring': 
    #     score = run_vina_scoring(receptor)
    #     sys.exit()
