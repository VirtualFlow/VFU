#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 19 17:23:20 2023

@author: akshat
"""

# Parameters:  
is_selfies     = False 
is_peptide     = False
program_choice = 'qvina' # smina/qvina/qvina-w/vina/vina_carb/vina_xb/gwovina/PLANTS/autodock_gpu/autodock_cpu/EquiBind/rDock/gnina/ledock/idock
                                # /autodock_vina/adfr/AutodockVina_1.2/AutodockZN/flexx/MM-GBSA/MCDock/LigandFit/GalaxyDock3/dock6/FRED/iGemDock/gold
                                 # /glide/rosetta-ligand/M-Dock/SEED/nnscore2/rf-score/molegro/FitDock/PSOVina/smina-scoring/gnina-scoring/AutoDock-Koto
                                 # /LightDock/RLDock/MpSDockZN/ad4_scoring/vinandro_scoring/vina_scoring
                                 
                                 
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
