#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Sep 25 20:01:20 2022

@author: akshat
"""
import os 
import time 
import subprocess
import multiprocessing
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

    print('Note: For use of vina gpu, the receptor needs to be prepared in a specific way. Have a look at the examples provided in https://github.com/ccsb-scripps/AutoDock-GPU & the example dir we provided within executables/vf_gpu_example.zip')
    command = []
    
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')
    if recetor_format[-1] != 'fld' and recetor_format[-2] != 'maps': 
        raise Exception('Receptor needs to be of file type .maps.fld (example: 1stp_protein.maps.fld). Please try again, after incorporating this correction.')
    
    # check for the existence of the executable: 
    try:
        if 'gpu' in program_choice: 
            executable = [x for x in os.listdir('./executables') if 'autodock_gpu' in x][0]
        elif 'cpu' in program_choice: 
            executable = [x for x in os.listdir('./executables') if 'autodock_cpu' in x][0]
        else: 
            raise Exception('Executable must be of format autodock_cpu/gpu')
    except: 
        raise Exception('Executable file autodock_cpu/gpu not found in executables directory')
           
    # Assign the right program for docking:  
    command.append('./executables/{}'.format(program_choice))
    
    # Prepare the ligands: 
    process_ligand(smi, 'pdbqt')
    lig_locations = os.listdir('./ligands/')

    # Get ready for running the files: 
    dlg_file = [x for x in os.listdir('./') if '.dlg' in x][0]
        
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
            results[lig_] = [docking_score, dlg_file]
            
    return results



def run_EquiBind(receptor, smi): 
    files_ls = os.listdir('./')
    if not('data' in files_ls and 'inference.py' in files_ls):
        raise Exception('Please make sure process EquiBind based on the instructions provided in the readme. I could not find the key EquiBind files.')
    
    if receptor.split('.')[2] != 'pdb': 
        raise Exception('For EquiBind, protein file needs to be in pdb file type. Please incorporate this correction')
    
    # Process the ligands
    process_ligand(smi, 'sdf') 

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
    process_ligand(smi, 'sd') 
    lig_locations = os.listdir('./ligands/')
    
    ref_lig = '' # TODO!!!
    raise Exception('Note: a reference ligand file needs to be filled in the line above. Please do so and comment this line!')
    
    # Creation of the prm file: 
    print('Please have a look at the prm parameters. Inside [TODO]; we have assigned some default values.')
    with open('config.prm', 'w') as f: 
        f.writelines('RBT_PARAMETER_FILE_V1.00\n')
        f.writelines('TITLE gart_DUD\n')
        f.writelines('RECEPTOR_FILE {}\n'.format(receptor))
        f.writelines('SECTION MAPPER\n')
        f.writelines('    SITE_MAPPER RbtLigandSiteMapper\n')
        f.writelines('    REF_MOL {}\n'.format(ref_lig))
        f.writelines('    RADIUS 6.0\n')
        f.writelines('    SMALL_SPHERE 1.0\n')
        f.writelines('    MIN_VOLUME 100\n')
        f.writelines('    MAX_CAVITIES 1\n')
        f.writelines('    VOL_INCR 0.0\n')
        f.writelines('   GRIDSTEP 0.5\n')
        f.writelines('END_SECTION\n')
        f.writelines('SECTION CAVITY\n')
        f.writelines('    SCORING_FUNCTION RbtCavityGridSF\n')
        f.writelines('    WEIGHT 1.0\n')
        f.writelines('END_SECTION\n')

    # Cavity generation: 
    os.system('rbcavity -was -d -r config.prm')
    
    results = {}
    
    # Perform docking: 
    os.system('mkdir rDock_outputs')
    for i,lig_ in enumerate(lig_locations): 
        os.system('mkdir rDock_outputs/{}'.format(i))
        lig_path = 'ligands/{}'.format(lig_)
        os.system('rbdock -i {} -o {} -r config.prm -p dock.prm -n 50'.format(lig_path, 'rDock_outputs/{}.sd'.format(i)))
    
        # Read the docking scores: 
        with open('rDock_outputs/{}.sd'.format(i), 'r') as f: 
            lines = f.readlines()
        score = []
        for i,item in enumerate(lines):
            if item.strip() == '>  <SCORE>': 
                score.append(float(lines[i+1]))
        docking_score = min(score)
        results[lig_] = docking_score
    return results



def generate_ledock_file(receptor='pro.pdb',rmsd=1.0,x=[0,0],y=[0,0],z=[0,0], n_poses=10, l_list=[],l_list_outfile='',out='dock.in'):
    rmsd=str(rmsd)
    x=[str(x) for x in x]
    y=[str(y) for y in y]
    z=[str(z) for z in z]
    n_poses=str(n_poses)

    with open(l_list_outfile,'w') as l_out:
        for element in l_list:
            l_out.write(element)
    l_out.close()

    file=[
        'Receptor\n',
        receptor + '\n\n',
        'RMSD\n',
        rmsd +'\n\n',
        'Binding pocket\n',
        x[0],' ',x[1],'\n',
        y[0],' ',y[1],'\n',
        z[0],' ',z[1],'\n\n',
        'Number of binding poses\n',
        n_poses + '\n\n',
        'Ligands list\n',
        l_list_outfile + '\n\n',
        'END']
    
    with open(out,'w') as output:
        for line in file:
            output.write(line)
    output.close()


def run_leDock(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z): 
    # Ensure receptor is in the right format
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')

    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in ledock
    lig_locations = os.listdir('./ligands/')

    results = {}
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        
        generate_ledock_file(receptor=receptor,
                             x=[center_x-size_x, center_x+size_x],
                             y=[center_y-size_y, center_y+size_y],
                             z=[center_z-size_z, center_z+size_z],
                             n_poses=10,
                             rmsd=1.0,
                             l_list= lig_path, 
                             l_list_outfile='ledock_ligand.list',
                             out='dock.in')
        
        ledock_cmd = ['./executables/ledock', 'dock.in']
        ledock_cmd = subprocess.run(ledock_cmd, capture_output=True)
        
        if ledock_cmd.returncode == 0: 
            os.system('cp ./ligands/{}.dok ./outputs/'.format(lig_.split('.')[0]))
            with open('./outputs/{}.dok'.format(lig_.split('.')[0]), 'r') as f: 
                lines = f.readlines()
            lines = [x for x in lines if 'Score' in x]
            scores = []
            for item in lines: 
                A = item.split('Score')[-1].strip().split(': ')[1].split(' ')[0]
                scores.append(float(A))
            results[lig_path] = ['./outputs/{}.dok'.format(lig_.split('.')[0]), min(scores)]
        
        else: 
            results[lig_path] = 'FAIL'
        
        os.system('rm dock.in ledock_ligand.list')
        
        
        
    
def process_idock_output(results): 
    poses_ = os.listdir('./')
    poses_ = [x for x in poses_ if 'pdbqt' in x]
    
    for item in poses_: 
        try: 
            ob_cmd = ['obenergy', './{}'.format(item)]
            command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
            command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
            total_energy         = float(command_obabel_check.split(' ')[-2])
        except: 
            total_energy = 10000 # Calculation has failed. 
        
        if total_energy < 10000: 
            
            # Read the output file: 
            with open('./log.csv', 'r') as f: 
                lines = f.readlines()
            lines = lines[1: ]
            map_ = {}
            for item in lines: 
                A =item.split(',')
                map_[A[0]] = float(A[2])
            
            # Save the result: 
            results[item] = map_[item.split('.')[0]]

            # Copy the file: 
            os.system('cp {} ./outputs/{}'.format(item, item))
        else: 
            os.system('rm {}'.format(item))
            
    os.system('rm log.csv')

    return 

def run_adfr_docking(receptor, smi): 
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')
    
    files_ls = os.listdir('./')
    target_file = [x for x in files_ls if '.trg' in x]
    if len(target_file) == 0: 
        raise Exception('A trg file containing all the parameters is required for running adfr. Please have a look at the tutorial in: https://ccsb.scripps.edu/adfr/documentation/')

    # prepare the ligands:
    process_ligand(smi, 'pdbqt') 
    lig_locations = os.listdir('./ligands/')
    
    print('Using target file: ', target_file[0])
    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        
        cmd = ['adfr', '-t', '{}'.format(target_file[0]), '-l', '{}'.format(lig_path), '--jobName', 'rigid']
        
        # Perform the docking: 
        command_run = subprocess.run(cmd, capture_output=True)
        
        if command_run.returncode != 0: 
            results[lig_path] = 'FAIL'
        
        docking_out = command_run.stdout.decode("utf-8")
        
        docking_scores = []
        for item in docking_out: 
            A = item.split(' ')
            A = [x for x in A if x != '']
            try: 
                a_1, a_2, a_3 = float(A[0]), float(A[1]), float(A[2])
            except: 
                continue
            docking_scores.append(float(a_2))

        results[lig_] = docking_scores

    return results


def run_flexx_docking(receptor, smi): 
    import multiprocessing
    results = {}

    ref_lig = 'ref_lig.mol2' # TODO
    ref_lig = os.listdir('./config')
    if 'ref_lig.mol2' not in ref_lig: 
        raise Exception('A reference ligand by the name of ref_lig.mol2 needs to be copied in the config directory for running flexx')
    
    executable_files = os.listdir('./executables')
    if 'flexx' not in executable_files: 
        raise Exception('The flexx executable was not foung. Please note: the execuatable file (named flexx) needs to be placed inside the directory executables')
        
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Receptor needs to be in pdb format')
    
    # prepare the ligands:
    process_ligand(smi, 'mol2') 
        
    lig_locations = os.listdir('./ligands/')

    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.sdf'.format(lig_.split('.')[0])
        
        os.system('./flexx -i {} -o {} -p {} -r {} --thread-count {}'.format(lig_path, out_path, receptor, './config/ref_lig.mol2', multiprocessing.cpu_count()))
        
        # Check energy of docked pose: 
        total_energy = check_energy(lig_)
        
        if total_energy < 10000: 

            # Read in the docking scores: 
            with open(out_path, 'r') as f: 
                lines = f.readlines()    
            
            for i,item in enumerate(lines): 
                docking_scores = []
                if '>  <docking-score>' in item : 
                    docking_score = float(lines[i+1])
                    docking_scores.append(docking_score)
            
            results[lig_] = [docking_scores, out_path]
        else: 
            results[lig_] = 'Extremely high pose energy encountered/Unsuccessfull execution.'
            
        return results

def run_AutodockZN(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness): 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'pdbqt': 
        raise Exception('Receptor needs to be in pdbqt format. Please try again, after incorporating this correction.')    
        
    if os.path.exists('~/ADFRsuite-1.0/bin/pythonsh') == False: 
        raise Exception('Could not locate ADFRsuite file (ADFRsuite-1.0/bin/pythonsh) in the home directory.')
    if os.path.exists('~/ADFRsuite-1.0/bin/autogrid4') == False: 
        raise Exception('Could not locate ADFRsuite file (ADFRsuite-1.0/bin/autogrid4) in the home directory.')
        
    # prepare the ligands:
    process_ligand(smi, 'pdbqt') 
    lig_locations = os.listdir('./ligands/')
    
    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])


        # Generate affinity maps: 
        os.system('~/ADFRsuite-1.0/bin/pythonsh ./config/prepare_gpf4zn.py -l {} -r {} -o receptor_tz.gpf -p npts={},{},{} -p gridcenter={},{},{} –p parameter_file=./config/AD4Zn.dat'.format(lig_path, receptor, size_x, size_y, size_z, center_x, center_y, center_z))
        os.system('~/ADFRsuite-1.0/bin/autogrid4 -p receptor_tz.gpf -l receptor_tz.glg')
        
        # Run AutoDockVina: 
        cmd = ['./config/AutodockVina_1.2', '--ligand', '{}'.format(lig_path), '--maps', 'receptor_tz', '--scoring', 'ad4', '--exhaustiveness', '{}'.format(exhaustiveness), '--out', '{}'.format(out_path)]
        command_run = subprocess.run(cmd, capture_output=True)

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
            
            results[lig_path] = [docking_out, out_path]
        
        # Delete auxilarry files: 
        os.system('rm receptor_tz.gpf receptor_tz.glg')
        
    return results

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


def run_mcdock(receptor, smi): 
    # Check to ensure receptor is in the right format: 
    if receptor.split('.')[-1] != 'xyz': 
        raise Exception('Please provide the receptor in xyz format for MCDock')
        
    # Check to ensure MCDock executable exists: 
    if os.path.exists('./executables/mcdock'): 
        raise Exception('Executable named mcdock not found in the executables directory')
        
    # Process all ligands: 
    process_ligand(smi, 'xyz') 
    lig_locations = os.listdir('./ligands/')
    
    results = {}

    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.xyz'.format(lig_.split('.')[0])
        
        # Run docking
        os.system('./executables/mcdock --target {} --ligand {}'.format(receptor, lig_path))
        
        # Read in the results: 
        with open('./out.xyz', 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if 'Binding Energy' in x]
        binding_energies = []
        for item in lines: 
            binding_energies.append(float(item.split(' ')[2].split('\t')[0]))    
            
        # Delete/move auxillary files: 
        os.system('rm min.xyz')
        os.system('cp out.xyz {}'.format(out_path))
        os.system('rm out.xyz conformers.xyz')
        
        results[lig_path] = binding_energies
        
    return results

def run_ligand_fit(receptor, smi, center_x, center_y, center_z): 
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Please provide the receptor in pdb format for LigandFit')
    
    if os.path.exists('./executables/ligandfit') == False: 
        raise Exception('Executable named ligandfit not found in the executables directory')
    if os.path.exists('./config/receptor.mtz') == False: 
        raise Exception('Receptor mtz file (titled receptor.mtz) not found in config directory. File is required for running LigandFit')

    # Process all ligands: 
    process_ligand(smi, 'pdb') 
    lig_locations = os.listdir('./ligands/')

    results = {}

    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.xyz'.format(lig_.split('.')[0])

        os.system('./executables/ligandfit data=./config/receptor.mtz model={} ligand={} search_center={},{},{}'.format(receptor, lig_path, center_x, center_y, center_z))
        
        # Read in results: 
        with open('./LigandFit_run_1_/ligand_1_1.log', 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if 'Best score' in x]
        scores = []
        for item in lines: 
            scores.append( float([x for x in item.split(' ') if x != ''][-2]) )
        
        # Remove auxillary file: 
        os.system('rm -rf PDS')
        os.system('cp ./LigandFit_run_1_/ligand_fit_1.pdb {}'.format(out_path))
        os.system('rm -rf LigandFit_run_1_')
        
        results[lig_] = [scores, out_path]
    
    return results

def run_GalaxyDock3(receptor, smi, center_x, center_y, center_z, exhaustiveness): 

    results = {}
    # Check to ensure receptor is in the right format: 
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Please provide the receptor in pdb format for MCDock')    

    # Check to ensure MCDock executable exists: 
    if os.path.exists('./executables/GalaxyDock3'): 
        raise Exception('Executable named GalaxyDock3 not found in the executables directory')
    if os.path.exists('./data'): 
        raise Exception('Data directory (./data) not found. Please have a look at the data directory in https://galaxy.seoklab.org/files/by2hsnvxjf/softwares/galaxydock.html (click the link galaxydock3)')
    
    # Process all ligands: 
    process_ligand(smi, 'mol2') 
    lig_locations = os.listdir('./ligands/')
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        

        # grid_n_elem    : Number of grid points for each directioni. This is should be
        #                  given in odd number. [61 61 61]
        grid_n_elem = [61, 61, 61]
        print('Note: i am setting grid_n_elem to [61, 61, 61]. Please change this default behaviour if need be. ')
        
        # grid_width     : Grid space between points in angstrom. [0.375]
        grid_width  = 0.375
        print('Note: i am setting grid_width to 0.375. Please change this default behaviour if need be. ')

        # Generate the input file: 
        with open('./galaxydock.in', 'w') as f: 
            f.writelines(['!==============================================\n'])
            f.writelines(['! I/O Parameters\n'])
            f.writelines(['!==============================================\n'])
            f.writelines(['data_directory    ./data\n'])
            f.writelines(['infile_pdb        {}\n'.format(receptor)])
            f.writelines(['infile_ligand        {}\n'.format(lig_path)])
            f.writelines(['top_type          polarh\n'])
            f.writelines(['fix_type          all\n'])
            f.writelines(['ligdock_prefix    out\n'])
            f.writelines(['!==============================================\n'])
            f.writelines(['! Grid Options\n'])
            f.writelines(['!==============================================\n'])
            f.writelines(['grid_box_cntr     {} {} {}\n'.format(center_x, center_y, center_z)])
            f.writelines(['grid_n_elem       {} {} {}\n'.format(grid_n_elem[0], grid_n_elem[1], grid_n_elem[2])]) 
            f.writelines(['grid_width        {}\n'.format(grid_width)])   
            f.writelines(['!==============================================\n'])
            f.writelines(['! Energy Parameters\n'])
            f.writelines(['!==============================================\n'])
            f.writelines(['weight_type              GalaxyDock3\n'])
            f.writelines(['!==============================================\n'])
            f.writelines(['! Initial Bank Parameters\n'])
            f.writelines(['!==============================================\n'])    
            f.writelines(['first_bank               rand\n'])
            f.writelines(['max_trial                {}\n'.format(exhaustiveness)])
            f.writelines(['e0max                    1000.0\n'])
            f.writelines(['e1max                    1000000.0\n'])
            f.writelines(['n_proc 1'])
            
        # Run the script: 
        os.system('./executables/GalaxyDock3 galaxydock.in > log')
        
        # Read in results: 
        with open('./out_fb.E.info', 'r') as f: 
            lines = f.readlines()
        lines = lines[3: ]
        docking_scores = []
        for item in lines: 
            try: 
                A = item.split(' ')
                A = [x for x in A if x != '']
                docking_scores.append(float(A[5]))
            except: 
                continue
        
        # Remove auxillary files
        os.system('rm log out_cl.E.info merged_ligand.mol2 out_cl.size.info out_co.info out_fb.E.info out_cl.mol2 out_ib.E.info out_ib.mol2')
        
        # Transfer out_fb.mol2
        os.system('cp {} {}'.format('out_fb.mol2', out_path))
        os.system('rm out_fb.mol2')
        
        results[lig_path] = [out_path, docking_scores]
        
    return results

def run_dock6(receptor, smi): 

    chimera_path  = '/home/akshat/chimera' # Please update the Chimera path 
    dock6_path    = '/home/akshat/dock6'   # Location to the dock6 directory
    
    ref_lig       = '/ref_lig.mol2' # Reference ligand needs to be specified for dock6
    box_padding   = 12.0
    
    results       = {}

    if os.path.exists(chimera_path) == False: 
        raise Exception('Location of Chemira not found (used location from variable chimera_path) when trying to initiate dock6 calculation.')
    if os.path.exists(dock6_path) == False: 
        raise Exception('Location of dock6 not found (used location from variable dock6_path) when trying to initiate dock6 calculation.')


    if os.path.exists(ref_lig) == False: 
        raise Exception('Please specify the location of the reference ligand for dock6 to run')

    # Check to ensure receptor is in the right format: 
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Please provide the receptor in pdb format for dock6')    

    
    # Prepare the receptor using Chimera: 
    os.system('{}/bin/chimera --nogui {} ./config/dockprep.py'.format(chimera_path, receptor))
    
    # doc6 pre-processing
    os.system('{}/bin/sphgen INSPH'.format(dock6_path))
    os.system('{}/bin/sphere_selector rec.sph {} {}'.format(dock6_path, ref_lig, box_padding))
    os.system('{}/bin/showbox < box.in'.format(dock6_path))
    os.system('{}/bin/grid -i grid.in'.format(dock6_path))
    
    # Process all ligands: 
    process_ligand(smi, 'mol2') 
    lig_locations = os.listdir('./ligands/')
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.xyz'.format(lig_.split('.')[0])
        
        # Generate dock6 input file: 
        with open('./dock.in', 'w') as f: 
            f.writelines(['conformer_search_type                                        flex\n'])
            f.writelines(['user_specified_anchor                                        no\n'])
            f.writelines(['limit_max_anchors                                            no\n'])
            f.writelines(['min_anchor_size                                              40\n'])
            f.writelines(['pruning_use_clustering                                       yes\n'])
            f.writelines(['pruning_max_orients                                          100\n'])
            f.writelines(['pruning_clustering_cutoff                                    100\n'])
            f.writelines(['pruning_conformer_score_cutoff                               25.0\n'])
            f.writelines(['pruning_conformer_score_scaling_factor                       1.0\n'])
            f.writelines(['use_clash_overlap                                            no\n'])
            f.writelines(['write_growth_tree                                            no\n'])
            f.writelines(['use_internal_energy                                          yes\n'])
            f.writelines(['internal_energy_cutoff                                       100.0\n'])
            f.writelines(['ligand_atom_file                                             {}\n'.format(lig_path)])
            f.writelines(['limit_max_ligands                                            no\n'])
            f.writelines(['receptor_site_file                                           selected_spheres.sph\n'])
            f.writelines(['max_orientations                                             500\n'])
            f.writelines(['chemical_matching                                            no\n'])
            f.writelines(['use_ligand_spheres                                           no\n'])
            f.writelines(['bump_filter                                                  no\n'])
            f.writelines(['score_molecules                                              yes\n'])
            f.writelines(['contact_score_primary                                        no\n'])
            f.writelines(['contact_score_secondary                                      no\n'])
            f.writelines(['grid_score_primary                                           yes\n'])
            f.writelines(['grid_score_secondary                                         no\n'])
            f.writelines(['grid_score_rep_rad_scale                                     1\n'])
            f.writelines(['grid_score_vdw_scale                                         1\n'])
            f.writelines(['grid_score_grid_prefix                                       grid\n'])
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
            f.writelines(['minimize_anchor                                              yes\n'])
            f.writelines(['minimize_flexible_growth                                     yes\n'])
            f.writelines(['use_advanced_simplex_parameters                              no\n'])
            f.writelines(['simplex_max_cycles                                           1\n'])
            f.writelines(['simplex_score_converge                                       0.1\n'])
            f.writelines(['simplex_cycle_converge                                       1.0\n'])
            f.writelines(['simplex_trans_step                                           1.0\n'])
            f.writelines(['simplex_rot_step                                             0.1\n'])
            f.writelines(['simplex_tors_step                                            10.0\n'])
            f.writelines(['simplex_anchor_max_iterations                                500\n'])
            f.writelines(['simplex_grow_max_iterations                                  500\n'])
            f.writelines(['simplex_grow_tors_premin_iterations                          0\n'])
            f.writelines(['simplex_random_seed                                          0\n'])
            f.writelines(['simplex_restraint_min                                        no\n'])
            f.writelines(['atom_model                                                   all\n'])
            f.writelines(['vdw_defn_file                                                {}/parameters\n'.format(dock6_path)])
            f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(dock6_path)])
            f.writelines(['flex_drive_file                                              {}/parameters/flex_drive.tbl\n'.format(dock6_path)])
            f.writelines(['vdw_defn_file                                                {}/parameters/vdw_AMBER_parm99.defn\n'.format(dock6_path)])
            f.writelines(['flex_defn_file                                               {}/parameters/flex.defn\n'.format(dock6_path)])
            f.writelines(['ligand_outfile_prefix                                        ligand_out\n'])
            f.writelines(['write_orientations                                           no\n'])
            f.writelines(['num_scored_conformers                                        1\n'])
            f.writelines(['rank_ligands                                                 no\n'])
            
        # Run the docking: 
        os.system('{}/bin/dock6 -i dock.in'.format(dock6_path))
        
        # Remove aux files: 
        os.system('rm dock.in')
        dock_file = [x for x in os.listdir('./') if 'ligand_out' in x]
        dock_file = [x for x in dock_file if 'mol2' in x][0]
        os.system('cp {} {}'.format(dock_file, out_path))
        os.system('rm dock_file')
        
        # Save the results: 
        with open('./ligand_out_scored.mol2', 'r') as f: 
            lines = f.readlines()
        docking_score = float(lines[2].split(' ')[-1])
            
        results[lig_path] = [out_path, docking_score]
        os.system('rm ligand_out_scored.mol2')
        
    return results

def run_fred_docking(receptor, smi, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness): 
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Please provide the receptor in pdb format for FRED')    
    
    if os.path.exists('./oe_license.txt') == False: 
        raise Exception('OpenEye licence file (oe_license.txt) not found in working directory')

    process_ligand(smi, 'mol2') 
    lig_locations = os.listdir('./ligands/')
    
    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        
        # Run docking: 
        os.system('python ./config/dock_fred.py --receptor-fn {} --ligand-fn {} --center-x {} --center-y {} --center-z {} --radius {} --num-poses {} --output-fn {}'.format(receptor, lig_path, center_x, center_y, center_z, max([size_x, size_y, size_z]), exhaustiveness, out_path))
        
        results[lig_path] = out_path
    
    return results

def run_iGemDock(receptor, smi, exhaustiveness): 
    if os.path.exists('./executables/mod_ga'): 
        raise Exception('iGemDock executable mod_ga not found in the executables folder')
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Please provide the receptor in pdb format for dock6')    
        
    process_ligand(smi, 'mol2') 
    
    lig_locations = os.listdir('./ligands/')

    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.xyz'.format(lig_.split('.')[0])
        
        # Perform Docking: 
        os.system('./executables/mod_ga {} {} {} -d ./'.format(exhaustiveness, receptor, lig_path))
                
        # Read in the poses: 
        docked_pose = os.listdir('./docked_Pose/')[0]
        os.system('cp {} {}'.format(docked_pose, out_path))
        
        with open(out_path, 'r') as f: 
            lines = f.readlines()
        docking_score = lines[4]
        docking_score = float([x for x in docking_score.split(' ') if x!=''][1])
        
        os.system('rm -rf docked_Pose')
            
        results[lig_path] = [out_path, docking_score]
    
    return results
        
def perform_gold_docking(receptor, smi, size_x, size_y, size_z, center_x, center_y, center_z): 
    if os.path.exists('./executables/gold_auto'): 
        raise Exception('Gold executable gold_auto not found in the executables folder')
    if receptor.split('.')[-1] != 'mol2': 
        raise Exception('Please provide the receptor in mol2 format for gold')        

    process_ligand(smi, 'mol2') 
    
    lig_locations = os.listdir('./ligands/')

    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.xyz'.format(lig_.split('.')[0])

        with open('input.conf', 'a+') as f: 
            f.writelines(['  GOLD CONFIGURATION FILE\n'])        
            f.writelines(['  AUTOMATIC SETTINGS\n'])        
            f.writelines(['autoscale = 1\n'])        
            f.writelines(['  POPULATION\n']) 
            f.writelines(['popsiz = auto\n'])        
            f.writelines(['select_pressure = auto\n'])        
            f.writelines(['n_islands = auto\n'])        
            f.writelines(['maxops = auto\n'])        
            f.writelines(['niche_siz = auto\n'])        
            f.writelines(['  GENETIC OPERATORS\n'])        
            f.writelines(['pt_crosswt = auto\n'])        
            f.writelines(['allele_mutatewt = auto\n'])        
            f.writelines(['migratewt = auto\n'])        
            f.writelines(['  FLOOD FILL\n'])        
            f.writelines(['radius = {}\n'.format(max([size_x, size_y, size_z]))])        
            f.writelines(['origin = {}   {}   {}\n'.format(center_x, center_y, center_z)])
            f.writelines(['do_cavity = 0\n'])        
            f.writelines(['floodfill_center = point\n'])        
            f.writelines(['   DATA FILES\n'])        
            f.writelines(['ligand_data_file {} 10\n'.format(lig_path)])        
            f.writelines(['param_file = DEFAULT\n'])        
            f.writelines(['set_ligand_atom_types = 1\n'])        
            f.writelines(['set_protein_atom_types = 0\n'])        
            f.writelines(['directory = out\n'])        
            f.writelines(['tordist_file = DEFAULT\n'])        
            f.writelines(['make_subdirs = 0\n'])        
            f.writelines(['save_lone_pairs = 1\n'])        
            f.writelines(['fit_points_file = fit_pts.mol2\n'])        
            f.writelines(['read_fitpts = 0\n'])        
            f.writelines(['bestranking_list_filename = bestranking.lst\n'])        
            f.writelines(['   FLAGS\n'])        
            f.writelines(['internal_ligand_h_bonds = 1\n'])        
            f.writelines(['flip_free_corners = 1\n'])        
            f.writelines(['match_ring_templates = 1\n'])        
            f.writelines(['flip_amide_bonds = 0\n'])        
            f.writelines(['flip_planar_n = 1 flip_ring_NRR flip_ring_NHR\n'])        
            f.writelines(['flip_pyramidal_n = 0\n'])        
            f.writelines(['rotate_carboxylic_oh = flip\n'])        
            f.writelines(['use_tordist = 1\n'])        
            f.writelines(['postprocess_bonds = 1\n'])        
            f.writelines(['rotatable_bond_override_file = DEFAULT\n'])        
            f.writelines(['solvate_all = 1\n'])        
            f.writelines(['   TERMINATION\n'])        
            f.writelines(['early_termination = 1\n'])        
            f.writelines(['n_top_solutions = 3\n'])        
            f.writelines(['rms_tolerance = 1.5\n'])        
            f.writelines(['   CONSTRAINTS\n'])        
            f.writelines(['force_constraints = 0\n']) 
            f.writelines(['   COVALENT BONDING\n'])        
            f.writelines(['covalent = 0\n']) 
            f.writelines(['   SAVE OPTIONS\n'])        
            f.writelines(['save_score_in_file = 1\n'])        
            f.writelines(['save_protein_torsions = 1\n'])        
            f.writelines(['  FITNESS FUNCTION SETTINGS\n'])        
            f.writelines(['initial_virtual_pt_match_max = 4\n'])        
            f.writelines(['relative_ligand_energy = 1\n'])        
            f.writelines(['gold_fitfunc_path = goldscore\n'])        
            f.writelines(['score_param_file = DEFAULT\n'])        
            f.writelines(['  PROTEIN DATA\n'])    
            f.writelines(['protein_datafile = {}\n'.format(receptor)])  
        
        # Run docking
        os.system('./executables/gold_auto input.conf')
        
        os.system('cp out/gold_ligand_m1.mol2 {}'.format(out_path))
        
        with open('./out/ligand_m1.rnk', 'r') as f: 
            lines = f.readlines()
            
        docking_score = float([x for x in lines[-1].split(' ') if x!=''][1])
        
        results[lig_path] = [out_path, docking_score]
    
    return results

def run_glide_docking(receptor, center_x, center_y, center_z, size_x, size_y, size_z, smi): 

    print('Note: The path to Schrodinger is specified via the $SCHRODINGER variable, which is assigned during installation.')    
    
    # Receptor prparation with Glide: 
    os.system('$SCHRODINGER/utilities/prepwizard {} receptor.maegz -fillsidechains -addOXT -epik_pH 7 -minimize_adj_h'.format(receptor))
    
    # Create grid file: 
    with open('./glide-grid.in', 'w') as f: 
        f.writelines(['FORCEFIELD   OPLS_2005\n'])
        f.writelines(['GRID_CENTER   {}, {}, {}\n'.format(center_x, center_y, center_z)])
        f.writelines(['GRIDFILE   receptor.zip\n'])
        f.writelines(['INNERBOX   {}, {}, {}\n'.format(size_x, size_y, size_z)])
        f.writelines(['OUTERBOX   {}, {}, {}\n'.format(size_x+15, size_y+15, size_z+15)])
        f.writelines(['RECEP_FILE   receptor.maegz\n'])
        
    # Run LigPrep; ligand processing with Glide: 
    with open('./grid_prep.sh', 'w') as f: 
        f.writelines(['while :\n'])
        f.writelines(['\tdo\n'])
        f.writelines(["\tif [[ `tail -n1 receptor.log|awk '{print $1}'` != 'DONE.' ]]; then\n"])
        f.writelines(['\t\tcontinue\n'])
        f.writelines(['\tfi\n'])
        f.writelines(["\tif [[ `tail -n1 receptor.log|awk '{print $1}'` == 'DONE.' ]]; then\n"])
        f.writelines(['\t\t$SCHRODINGER/glide glide-grid.in\n'])
        f.writelines(['\t\tbreak\n'])
        f.writelines(['\tfi\n'])
        f.writelines(['done\n'])
        
    # Wait 3mins for protein preparation to finish: 
    time.sleep(180) 
    os.system('chmod 777 ./grid_prep.sh')
    os.system('./grid_prep.sh')
    # Wait 3mins for grid preparation to finish: 
    time.sleep(180) 

    process_ligand(smi, 'sd') 
    
    lig_locations = os.listdir('./ligands/')

    results = {}
    
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/complex_{}.maegz'.format(lig_.split('.')[0])
        
        # Run LigPrep; ligand processing with Glide: 
        with open('./lig_prep.sh', 'w') as f: 
            f.writelines(['while :\n'])
            f.writelines(['\tdo\n'])
            f.writelines(["\tif [[ `tail -n1 glide-grid.log|awk '{print $1}'` != 'Total' ]]; then\n"])
            f.writelines(['\t\tcontinue\n'])
            f.writelines(['\tfi\n'])
            f.writelines(["\tif [[ `tail -n1 glide-grid.log|awk '{print $1}'` != 'Total' ]]; then\n"])
            f.writelines(['\t\t$SCHRODINGER/ligprep -isd {} -omae ligand.mae -epik -ph 7\n'.format(lig_path)])
            f.writelines(['\t\tbreak\n'])
            f.writelines(['\tfi\n'])
            f.writelines(['done\n'])
            
        os.system('chmod 777 ./lig_prep.sh; ./lig_prep.sh\n')
        # Wait 30seconds for ligand preparation to finish: 
        time.sleep(30) 
        
        # Perform molecular docking: 
        with open('./glide-dock.in', 'w') as f: 
            f.writelines(['FORCEFIELD   OPLS_2005\n'])
            f.writelines(['GRIDFILE   receptor.zip\n'])
            f.writelines(['LIGANDFILE   ligand.mae\n'])
            f.writelines(['PRECISION   SP\n'])
        
        with open('./docking_run.sh', 'w') as f: 
            f.writelines(['while :\n'])
            f.writelines(['\tdo\n'])
            f.writelines(["\tif [[ `tail -n1 ligand.log|awk '{print $1}'` != 'backend' ]]; then\n"])
            f.writelines(['\t\tcontinue\n'])
            f.writelines(['\tfi\n'])
            f.writelines(["\tif [[ `tail -n1 ligand.log|awk '{print $1}'` != 'backend' ]]; then\n"])
            f.writelines(['\t\t$SCHRODINGER/glide glide-grid.in\n'])
            f.writelines(['\t\tbreak\n'])
            f.writelines(['\tfi\n'])

        os.system('chmod 777 ./docking_run.sh; ./docking_run.sh')
        # Wait 60seconds for docking to finish: 
        time.sleep(60) 

        with open('./glide-dock.csv', 'r') as f: 
            lines = f.readlines()
        docking_scores = [] # Collect the docking score
        for item in lines[1: ]: 
            A = item.split(',')
            docking_scores.append(float(A[5]))

        # Delete aux files: 
        os.system('rm glide-dock.csv glide-dock.log receptor.log glide-grid.log; cp glide-dock_pv.maegz {}; rm glide-dock_pv.maegz'.format(out_path))

        results[smi] = [docking_scores, out_path]
    return results
    
def run_rosetta_docking(receptor, smi, center_x, center_y, center_z, exhaustiveness): 
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
        
    if os.path.exists('$ROSETTA/source/scripts/python/public/molfile_to_params.py') == False: 
        raise Exception('Rosetta file, located in $ROSETTA/source/scripts/python/public/molfile_to_params.py could not be found.')
    if os.path.exists('$ROSETTA/source/bin/rosetta_scripts.default.linuxgccrelease') == False: 
        raise Exception('Rosetta file, located in $ROSETTA/source/bin/rosetta_scripts.default.linuxgccrelease could not be found.')


    # prepare the ligands:
    with open('./test.smi', 'w') as f: 
        f.writelines([smi])
    os.system('obabel ./test.smi --gen3D -O ligand.mol2')
    os.system('rm ./test.smi')
    # Generate conformational library for ligand: 
    os.syetem('obabel ligand.mol2 -O conformers.sdf --conformer --score rmsd --writeconformers --nconf 30')
    os.system('rm ./ligand.mol2')

    os.system('$ROSETTA/source/scripts/python/public/molfile_to_params.py -n LIG -p LIG --conformers-in-one-file conformers.sdf')
    
        
    # run molecular docking: 
    with open('./run_docking.sh', 'w') as f: 
        f.writelines(["$ROSETTA/source/bin/rosetta_scripts.default.linuxgccrelease  \\\n"])
        f.writelines(["	-database $ROSETTA/database \\\n"])
        f.writelines(["\t@ options \\\n"])
        f.writelines(["\t\t-parser:protocol dock.xml \\\n"])
        f.writelines(["\t\t-parser:script_vars X={} Y={} Z={} \\\n".format(center_x, center_y, center_z)])
        f.writelines(["\t\t-in:file:s complex.pdb \\\n"])
        f.writelines(["\t\t-in:file:extra_res_fa LIG.params \\\n"])
        f.writelines(["\t\t-out:nstruct 10 \\\n"])
        f.writelines(["\t\t-out:level {} \\\n".format(exhaustiveness)])
        f.writelines(["\t\t-out:suffix out\n"])
        
    os.system('chmod 777 ./run_docking.sh')
    os.system('./run_docking.sh')
            
    # Collect output files: 
    results = {}
    results[smi] = []
    
    out_files = [x for x in os.listdir('./') if 'complexout' in x]    
    for file in out_files: 
        os.system('cp {} ./outputs/{}'.format(file, file))
        os.system('rm {}'.format(file))
        results[smi].append('./outputs/{}'.format(file))
    
    with open('./scoreout.sc', 'r') as f: 
        lines = f.readlines()
    lines = lines[2: ]
    docking_scores = []
    for item in lines: 
        A = item.split(' ')
        A = [x for x in A if x!='']
        docking_scores.append(float(A[44]))
    
    results[smi].append(docking_scores)

    return results

def run_mdock_docking(receptor, smi): 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'sph': 
        raise Exception('Receptor needs to be in sph format. Please try again, after incorporating this correction.')
    
    ref_lig = './config/ref_lig.pdb'
    if os.path.exists(ref_lig) == False: 
        raise Exception('Reference ligand {} not found. Please try again, after incorporating this correction.'.format(ref_lig))

    mdock_path = '/home/MDock'

    if os.path.exists('{}/bin/MDock_Linux'.format(mdock_path)) == False: 
        raise Exception('MDock path {} not found. Please try again, after incorporating this correction.'.format(mdock_path))

    os.system('{}/bin/get_sph_Linux {} {}'.format(mdock_path, ref_lig, receptor))
    
    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')
      
    results = {}
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        
        with open('./mdock_dock.par', 'w') as f: 
            f.writelines('clash_potential_penalty      |      3.0\n')
            f.writelines('orient_ligand (yes/no)       |      yes\n')
            f.writelines('minimize_ligand (yes/no)     |      yes\n')
            f.writelines('maximum_orientations         |      100\n')
            f.writelines('gridded_score (yes/no)       |      yes\n')
            f.writelines('grid_spacing (0.3~0.5)       |      0.4\n')
            f.writelines('sort_orientations (yes/no)   |      yes\n')
            f.writelines('write_score_total            |      100\n')
            f.writelines('write_orientations (yes/no)  |      yes\n')
            f.writelines('minimization_cycles (1~3)    |      1\n')
            f.writelines('ligand_selectivity (yes/no)  |      no\n')
            f.writelines('box_filename (optional)      |      \n')
            f.writelines('grid_box_size                |      10.0\n')
            f.writelines('sphere_point_filename        |      recn.sph\n')
            
            
        # $MDock/bin/MDock_Linux protein ligand.mol2 -param mdock_dock.par
        os.system('{}/bin/MDock_Linux protein {} -param mdock_dock.par'.format(mdock_path, lig_path))

        os.system('cp ./mdock_dock.mol2 {}'.format(out_path))

        docking_scores = []
        with open('./mdock_dock.out', 'r') as f: 
            lines = f.readlines()
        for item in lines: 
            docking_scores.append( float([x for x in item.split(' ') if x != ''][4]))
            
        results[lig_path] = [out_path, docking_scores]
            
    return results 

def run_seed_docking(receptor, smi): 
    chimera_path = '/home/chimera'
    seed_path    = '/home/SEED'
    if os.path.exists(chimera_path) == False: 
        raise Exception('Chimera path {} not found. Please try again, after updating the Chimera path. '.format(chimera_path))
    if os.path.exists(seed_path) == False: 
        raise Exception('SEED path {} not found. Please try again, after updating the SEED path. '.format(seed_path))
    if os.path.exists(seed_path+'/bin/seed_4') == False: 
        raise Exception('SEED executable {} not found. Please try again, after updating the SEED executable path. '.format(seed_path+'/bin/seed_4'))

    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'mol2': 
        raise Exception('Receptor needs to be in mol2 format. Please try again, after incorporating this correction.')

    # Prepare receptor for SEED        
    os.system('python ./config/mol2seed4_receptor.py {} {} receptor_seed.mol2'.format(receptor, receptor))

    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')
      
    results = {}
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        
        # Preparing Ligand for SEED: 
        os.system('charge=`{}/bin/chimera --nogui --silent {} ./config/charges.py`; antechamber -i {} -fi mol2 -o ligand_gaff.mol2 -fo mol2 -at gaff2 -c gas -rn LIG -nc $charge -pf y'.format(chimera_path, lig_path, lig_path))
        os.system('python ./config/mol2seed4_receptor.py ligand_gaff.mol2 ligand_gaff.mol2 ligand_seed.mol2')
        
        os.system('{}/bin/seed_4 seed.inp > log'.format(seed_path))
        
        os.system('cp ./ligand_seed_best.mol2 {}'.format(out_path))
        
        with open('./seed_best.dat', 'r') as f: 
            lines = f.readlines()
        docking_score = float([x for x in lines[1].split(' ') if x != ''][4])
        results[lig_path] = [docking_score, out_path]
        
        os.system('rm apolar_rec.mol2 apolar_rec_reduc.mol2 length_hb.gen ligand_seed_best.mol2 polar_rec.mol2 polar_rec_reduc.mol2 sas_apolar.pdb seed.out seed_best.dat seed_clus.dat')

    return results


def run_molegro_docking(receptor, smi): 
    
    receptor_format = receptor.split('.')[-1]
    if receptor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
    if os.path.exists(receptor) == False: 
        raise Exception('Recpetion path {} not found.'.format(receptor))
    
    molegro_path = '$HOME/MVD'
    if os.path.exists(molegro_path) == False: 
        raise Exception('Molegro directory path {} not found. Please try again, after updating the location for Molegro. '.format(molegro_path))
        
    ref_ligand = ''
    if os.path.exists(ref_ligand) == False: 
        raise Exception('Molegro requires a reference ligand for docking. Reference ligand pat {} was not found. Please try again, after updating the path.'.format(ref_ligand))
        
    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')
      
    results = {}
    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        
        os.system('cat {} {} > ligands.mol2'.format(ref_ligand, lig_path))
        
        with open('docking.mvdscript', 'a+') as f: 
            f.writelines(['// Molegro Script Job.\n\n'])
            f.writelines(['IMPORT Proteins;Waters;Cofactors FROM {}\n\n'.format(receptor)])
            f.writelines(['PREPARE Bonds=IfMissing;BondOrders=IfMissing;Hydrogens=IfMissing;Charges=Always; TorsionTrees=Always\n\n'])
            f.writelines(['IMPORT All FROM ligands.mol2\n\n'])
            f.writelines(['SEARCHSPACE radius=12;center=Ligand[0]\n\n'])
            f.writelines(['DOCK Ligand[1]\n\n\n'])
            f.writelines(['EXIT'])

        os.system('$Molegro/bin/mvd docking.mvdscript -nogui')
        
        cmd_ = [molegro_path, '.docking.mvdscript', '-nogui']
        cmd_run = subprocess.run(cmd_, capture_output=True)
        cmd_run = cmd_run.stdout.decode("utf-8").split('\n')[-2]
        cmd_run = [x for x in cmd_run if 'Pose:' in x]
        scores = []
        for item in cmd_run: 
            scores.append( float(item.split('Energy')[-1].split(' ')[1][:-2]) )

        out_path_all = []
        mol_files = [x for x in os.listdir('./') if 'mol2' in x]
        for i,item in mol_files: 
            out_path_aug = out_path.split('.mol2')[0] + '_{}_'.format(i) + '.mol2'
            out_path_all.append(out_path_aug)
            os.system('cp {} {}'.format(item, out_path_aug))
        
        del_files = [x for x in os.listdir('./') if '.txt' in x or '.mvdresults' in x]
        for item in del_files: os.system('rm {}'.format(item))
        
        results[lig_path] = [scores, out_path_all]
        
    return results

def run_fitdock_docking(receptor, smi): 

    fitdock_executable = './execuatable/FitDock'
    if os.path.exists(fitdock_executable) == False: 
        raise Exception('FitDock executable path {} not found. Please try again, after adding the executable (named FitDock) in the right directory.'.format(fitdock_executable))
        
    results = {}

    receptor_template = ''
    if os.path.exists(receptor_template) == False: 
        raise Exception('Receptor template path {} not found. Please try again, after specifying the location for the template receptor (right above this line).'.format(fitdock_executable))

    ligand_reference = ''
    if os.path.exists(ligand_reference) == False: 
        raise Exception('Ligand reference path {} not found. Please try again, after specifying the location for the reference ligand (right above this line).'.format(ligand_reference))
        
    if receptor_template.split('.')[-1] != 'pdb': 
        raise Exception('Receptor template needs to be in pdb format. Please try again, after incorporating this correction.')
    if ligand_reference.split('.')[-1] != 'mol2': 
        raise Exception('Reference ligand needs to be in mol2 format. Please try again, after incorporating this correction.')
    if receptor.split('.')[-1] != 'pdb': 
        raise Exception('Processed protein needs to be in pdb format. Please try again, after incorporating this correction.')

    # prepare the ligands:
    process_ligand(smi, 'mol2') # mol2 ligand format is supported in plants
    lig_locations = os.listdir('./ligands/')

    for lig_ in lig_locations: 
        lig_path = 'ligands/{}'.format(lig_)
        out_path = './outputs/pose_{}.mol2'.format(lig_.split('.')[0])
        
        os.system('./execuatable/FitDock -Tprot {} -Tlig {} -Qprot {} -Qlig {} -ot ot.mol2 -os os.mol2 -o o.mol2 > out.log'.format(receptor_template, ligand_reference, receptor, lig_path))
        os.system('cp o.mol2 {}'.format(out_path))
        
        with open('./out.log', 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if 'Binding Score after  EM' in x]
        docking_score = float(lines[0].split(' ')[-2])
            
        results[lig_path] = [out_path, docking_score]
        os.system('rm o.mol2 os.mol2 ot.mol2 out.log')
    
    return results 
        

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
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')
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
    
def run_lightdock_docking(receptor, smi, exhaustiveness): 
    results = {}
    
    # receptor needs to be in mol2 format: 
    recetor_format = receptor.split('.')[-1]
    if recetor_format != 'pdb': 
        raise Exception('Receptor needs to be in pdb format. Please try again, after incorporating this correction.')    

    if os.path.exists('./lightdock') == False: 
        raise Exception('Directory ./lightdock not found. Please try again, after incorporating this correction.')
        
    # prepare the ligands:
    process_ligand(smi, 'pdb') 
    lig_locations = os.listdir('./ligands/')
    
    os.system('cp {} ./'.format(receptor))
    receptor_short = receptor.split('/')[-1]
    
    for i,lig_ in enumerate(lig_locations): 
        lig_path = 'ligands/{}'.format(lig_)
        os.system('cp {} ./'.format(lig_path))
        lig_short = lig_path.split('/')[1]
        
        # LightDock setup:         
        os.system('./lightdock/bin/lightdock3_setup.py {} {} --noxt --noh --now -anm'.format(receptor_short, lig_short))
        # Docking: 
        os.system('./lightdock/bin/lightdock3.py setup.json 100 -c 1 -l 0')
        # Save poses: 
        os.system('./lightdock/bin/lgd_generate_conformations.py {} {} swarm_0/gso_100.out {}'.format(receptor_short, lig_short ,exhaustiveness))

        dir_ls = os.listdir('./swarm_0/')
        complex_all = [x for x in dir_ls if '.pdb' in x]
        
        with open('./swarm_0/gso_100.out', 'r') as f: 
            lines = f.readlines()
        lines = lines[1: ]
        scoring = []
        for item in lines: 
            A = item.split(' ')
            scoring.append(float(A[-1]))

        os.system('cp -a swarm_0 ./outputs/{}'.format(i))
        results[lig_] = [scoring, ['./outputs/{}/'.format(i)+x for x in complex_all]]
        print(results)

        os.system('rm -rf swarm* init')
        os.system('rm lightdock_lig.nm.npy lightdock_peptide.pdb lightdock_peptide_mask.npy lightdock_rec.nm.npy lightdock_receptor.pdb lightdock_receptor_mask.npy setup.json lightdock.info lightdock.info.1 {}'.format(lig_short))
        rem_files = [x for x in os.listdir('./') if 'lightdock' in x and x != 'lightdock']
        for x in rem_files: os.system('rm {}'.format(x))        

    os.system('rm {}'.format(receptor_short))
    return results

def run_RLDock_docking(receptor, smi, exhaustiveness): 
    
    if receptor.split('.')[2] != 'mol2': 
        raise Exception('Receptor file needs to be mol2 file type. Please incorporate this correction')
        
    # Process the ligands
    process_ligand(smi, 'mol2')
    
    results = {}
    lig_locations = os.listdir('./ligands/')
    for i,lig_ in enumerate(lig_locations): 
        lig_path = 'ligands/{}'.format(lig_)
        output_path = './outputs/{}'.format(lig_)
        
        os.system('./executables/RLDOCK -i {} -l {} -c {} -n {} -s config/sphere.dat'.format(receptor, lig_path, exhaustiveness, multiprocessing.cpu_count()))
        
        with open('./output_cluster.mol2', 'r') as f: 
            lines = f.readlines()
        lines = [x for x in lines if '# Total_Energy:' in x]
        docking_scores = []
        for item in lines: 
            docking_scores.append(float(item.split(' ')[-1]))
        
        del_files = [x for x in os.listdir('./') if 'output' in x and x != 'outputs']
        os.system('cp output_cluster.mol2 {}'.format(output_path))
        for file_ in del_files: os.system('rm {}'.format(file_))
        results[lig_path] = [output_path, docking_scores]
        
    return results


def run_MpSDockZN_docking(receptor, smi): 
    results = {}
    
    dock6_path   = ''
    chimera_path = ''    
    box_in_file  = ''
    grid_in_file  = ''
    dock_in_file  = ''

    if receptor.split('.')[2] != 'pdb': 
        raise Exception('Receptor file needs to be pdb file type. Please incorporate this correction')
    
    if os.path.exists(dock6_path) == False: 
        raise Exception('dock6 path not found. Please incorporate this correction')
    if os.path.exists(chimera_path) == False: 
        raise Exception('chimera path not found. Please incorporate this correction')
    if os.path.exists(grid_in_file) == False: 
        raise Exception('grid.in file not found. Please incorporate this correction')
    if os.path.exists(dock_in_file) == False: 
        raise Exception('dock.in file not found. Please incorporate this correction')
        
    # Process the ligands
    process_ligand(smi, 'mol2')
    
    results = {}
    lig_locations = os.listdir('./ligands/')
    for i,lig_ in enumerate(lig_locations): 
        lig_path = 'ligands/{}'.format(lig_)
        output_path = './outputs/{}'.format(lig_)
        
        with open('./run.sh', 'w') as f: 
            f.writelines(['export Chimera={}\n'.format(chimera_path)])
            f.writelines(['export DOCK6={}\n'.format(dock6_path)])
            f.writelines(['charge=`$Chimera/bin/chimera --nogui --silent {} charges.py`\n'.format(lig_path)])
            f.writelines(['antechamber -i {} -fi mol2 -o ligand_input.mol2 -fo mol2 -at sybyl -c gas -rn LIG -nc $charge -pf y\n'.format(lig_path)])
            f.writelines(['$DOCK6/bin/showbox < {} \n'.format(box_in_file)])
            f.writelines(['$DOCK6/bin/grid -i {} \n'.format(grid_in_file)])
            f.writelines(['./executables/MpSDock -i {} \n'.format(dock_in_file)])
            
        os.system('chmod 777 ./run.sh; ./run.sh')
        os.system('cp receptor_input_docked_result.mol2 {}'.format(output_path))
        
        score_all = []
        with open('./receptor_input_docked_result.list', 'r') as f: 
            lines = f.readlines()
        for item in lines: 
            A = item.split(' ')
            A = [x for x in A if x != '']
            try: score_1, score_2, score_3, score_4, score_5 = float(A[0]), float(A[1]), float(A[2]),float(A[3]), float(A[4])
            except: continue 
            final_score = score_1 + score_2 + score_3 + score_4 + score_5
            score_all.append(final_score)
        
        results[lig_path] = [output_path, score_all]
        
    return results 


def check_energy(lig_): 
    # Check the quality of generated structure (some post-processing quality control):
    try: 
        ob_cmd = ['obenergy', './outputs/pose_{}.pdbqt'.format(lig_.split('.')[0])]
        command_obabel_check = subprocess.run(ob_cmd, capture_output=True)
        command_obabel_check = command_obabel_check.stdout.decode("utf-8").split('\n')[-2]
        total_energy         = float(command_obabel_check.split(' ')[-2])
    except: 
        total_energy = 10000 # Calculation has failed. 
        
    return total_energy
