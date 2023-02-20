#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 15:43:49 2022

@author: akshat
"""
import os 
import subprocess
import numpy as np 
from rdkit import Chem
from rdkit.Chem import MolToSmiles as mol2smi
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions


def desalt_smi(smi): 
    '''
    Take the largest molecule (with the most number of atoms) in the presence 
    multiple conjoined smiles (seperated by '.')

    Parameters
    ----------
    smi : str
        Smile string of a valid molecule.

    Returns
    -------
    smi_ : str
        Cleaned up (desalted) smiles string that is the largest component of smi.
    '''

    if '.' in smi: 
        smi_all   = smi.split('.')
        mol_split = [Chem.MolFromSmiles(x) for x in smi_all]
        num_atoms = [m.GetNumAtoms() for m in mol_split]
        idx_large = np.argmax(num_atoms)
        smi_      = smi_all[idx_large]
    else: 
        return smi
        
    return smi_


def enumerate_sterio(smi, asigned=True): 
    '''
    Enumerate all sterioisomers of the provided molecule smi. 
    Note: Only unspecified stereocenters are expanded. 

    Parameters
    ----------
    smi : str
         Valid molecule SMILE string.
    asigned: bool
         if True, isomers will be generated for only the unasigned sterio-locations 
                  for a smile (faster)
         if False, all isomer combinations will be generated, regardless of what is 
                  specified in the input smile (slower)

    Returns
    -------
    sterio_smiles: list of strs.
         A list of valid smile strings, representing sterioisomers.

    '''

    m = Chem.MolFromSmiles(smi)
    
    if asigned == True: # Faster
        opts = StereoEnumerationOptions(unique=True)
    else: 
        opts = StereoEnumerationOptions(unique=True, onlyUnassigned=False)
    
    isomers = tuple(EnumerateStereoisomers(m, options=opts))
    
    sterio_smiles = []
    for smi in sorted(Chem.MolToSmiles(x,isomericSmiles=True) for x in isomers):
        sterio_smiles.append(smi)
        
    return sterio_smiles


def reorderTautomers(m):
    '''
    Generate all posible tautomers of an input molecule. 
    The first tautomer is RdKit's canonical tautomer. 

    Parameters
    ----------
    m : RdKit mol object
        Valid (not-none) mol object generated via RdKit.

    Returns
    -------
    res : list of mol objects
        A collection of all tautomers of input molecule, returned as mol objects.
    '''
    enumerator = rdMolStandardize.TautomerEnumerator()
    canon = enumerator.Canonicalize(m)
    csmi = Chem.MolToSmiles(canon)
    res = [canon]
    tauts = enumerator.Enumerate(m)
    smis = [Chem.MolToSmiles(x) for x in tauts]
    stpl = sorted((x,y) for x,y in zip(smis,tauts) if x!=csmi)
    res += [y for x,y in stpl]
    
    return res


def neutralize_atoms(mol):
    '''
    Neutralize molecule such that the net charge is 0. 
    Used from: https://www.rdkit.org/docs/Cookbook.html
    
    Parameters
    ----------
    mol : RdKit mol object
        Valid (not-none) mol object generated via RdKit.

    Returns
    -------
    mol : RdKit mol object
        Desalted RdKit mol object of orginal molecule.

    '''
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches] 
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol

def check_energy(lig_): 
    """
    Check the quality of a generated structure by computing its total energy using the Open Babel obenergy tool.

    Parameters:
        lig_ (str): the name of the ligand file in PDBQT format.

    Returns:
        total_energy (float): the computed total energy of the ligand in Kcal/mol.
    """
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', lig_]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy


def process_ligand(smi, output_format, asigned_sterio=True, max_sterio=8, max_tautomer=2): 
    '''
    Convert provided smile string into 3-d coordinates, ready for molecular docking. 
    Steps performed include: 
        1. Desalting & neutralization
        2. Tautomer generation 
        3. Generating protonation state (for each tautomer) at pH 7.4
    the 3-d coordinates (one per tautomer) is saved in the 'ligands' directory
    

    Parameters
    ----------
    smi : str
        valid smile string of a ligand.
    output_format: str
        Type of 3D mol file (pdb, pdbqt, sdf, mol2)
    asigned_sterio: bool
         if True, isomers will be generated for only the unasigned sterio-locations 
                  for a smile (faster)
         if False, all isomer combinations will be generated, regardless of what is 
                  specified in the input smile (slower)
    Returns
    -------
    None.

    '''
    output_format = output_format.lower().strip()
    if not(output_format == 'sdf' or output_format == 'mol2' or output_format == 'pdb' or output_format == 'pdbqt' or output_format == 'sd'): 
        print('Unsopported file type requested')
        return 
    
    mol = Chem.MolFromSmiles(smi)
    if mol == None: 
        print('Invalid starting molecule provided! A valid SMILE is required for the program.')
        return 
    
    # Desault the input molecule
    smi = desalt_smi(smi) 
    mol = Chem.MolFromSmiles(smi)
    
    # Neutralize the molecule such that the net charge is 0
    mol = neutralize_atoms(mol) 
    
    # Enumerate sterio-isomers of the molecules: 
    smi = mol2smi(mol)
    
    sterio_smiles = enumerate_sterio(smi, asigned_sterio)
    sterio_smiles = sterio_smiles[0: max_sterio]

    # Enumerate all tautomers of the sterio-isomers: 
    smiles_all = []
    for item in sterio_smiles: 
        tautomers = reorderTautomers(Chem.MolFromSmiles(item))
        tautomers = [mol2smi(x) for x in tautomers]
        smiles_all.extend(tautomers[0:max_tautomer])
    
    for i,smi in enumerate(smiles_all): 
        
        with open('./test.smi', 'w') as f: 
            f.writelines(smi)        

        cmd = ['obabel', 'test.smi', '--gen3d', '-O', './ligands/{}.{}'.format(i, output_format), '-p', '7.4']

        try: 
            command_run = subprocess.run(cmd, capture_output=True, timeout=10)
        except: 
            os.system('rm ./ligands/{}.{}'.format(i, output_format))
            continue 

        if command_run.returncode == 0: 
            energy_val = check_energy('./ligands/{}.{}'.format(i, output_format))
            if energy_val >= 10000: 
                print('Warning: Bad ligand conformer encountered (removing file): {}'.format(smi))
                os.system('rm ./ligands/{}.{}'.format(i, output_format))
        else: 
            os.system('rm ./ligands/{}.{}'.format(i, output_format))

    os.system('rm test.smi')    
    
    
if __name__ == '__main__': 
    A = process_ligand('BrC=CC1OC(C2)(F)C2(Cl)C1.CC.[Cl][Cl]', 'pdbqt', asigned_sterio=True)
    

    
