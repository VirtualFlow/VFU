#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 15:43:49 2022

@author: akshat
"""
import os 
import numpy as np 
from rdkit import Chem
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import MolToSmiles as mol2smi


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

def process_ligand(smi, output_format): 
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

    Returns
    -------
    None.

    '''
    output_format = output_format.lower().strip()
    if not(output_format == 'sdf' or output_format == 'mol2' or output_format == 'pdb' or output_format == 'pdbqt'): 
        print('Unsopported file type requested')
        return 
    
    mol = Chem.MolFromSmiles(smi)
    if mol == None: 
        print('Invalid starting molecule provided! A valid SMILE is required for the program.')
        return 
    
    smi = desalt_smi(smi) # Desault the input molecule
    
    mol = neutralize_atoms(mol) # Neutralize the molecule such that the net charge is 0
    
    tautomers = reorderTautomers(mol) # Generate all posible tautomers of a molecule
    smiles_all = [mol2smi(x) for x in tautomers]
        
    for i,smi in enumerate(smiles_all): 
        
        with open('./test.smi', 'w') as f: 
            f.writelines(smi)        
        os.system('obabel test.smi --gen3d  -O ./ligands/{}.{}'.format(i, output_format)) # Protonation states at pH 7.4 is used!
        
    os.system('rm test.smi')    
    
    
    
if __name__ == '__main__': 
    process_ligand('Oc1c(cccc3)c3nc2ccncc12.CC.[Cl][Cl]', 'sdf')