import os 
import subprocess
import multiprocessing 


def get_docking_score():
    with open('./temp_dir/DOCKING_TEST_log.txt', 'r') as f:
        lines = f.readlines()
        
    lines = lines[25].split(' ')
    lines = [x for x in lines if len(x) > 0]
    lines = float(lines[1])
    
    return lines

cpu_count = multiprocessing.cpu_count()



smi = 'CC1=CC2=C(C=C1)C(C)=CC=C2'

# Convert smile -> sdf for docking 
with open('./test.smi', 'w') as f: 
    f.writelines(smi)
    
os.system('obabel ./test.smi --gen3D -O test.sdf')
os.system('rm test.smi')

# Non-flexible smina docking: 
# command = './executables/smina -r ./protein/prot.pdb -l ./test.sdf --exhaustiveness 10 --autobox_ligand ./lig/lig.pdb --autobox_add 3 -o ./temp_dir/DOCKING_TEST.pdb --log ./temp_dir/DOCKING_TEST_log.txt'.format(cpu_count)

# Command for smina (on pi3k failure): 
command = './executables/smina -r ./protein/prot_1.pdb -l ./test.sdf --exhaustiveness 10 --autobox_ligand ./lig/lig_1.pdb --autobox_add 3 -o ./temp_dir/DOCKING_TEST.pdb --log ./temp_dir/DOCKING_TEST_log.txt'


os.system(command)


#ret = os.system(command) 
#ret = subprocess.run(command)

#ret = subprocess.run(['./executables/smina'], capture_output=True, text=True)
# A = subprocess.run(['./executables/smina', '-r', './protein/prot.pdb',
#                     '-l', './test.sdf', '--exhaustiveness', '10', 
#                     '--autobox_ligand', './lig/lig.pdb', '--autobox_add', '3', '-o', './temp_dir/DOCKING_TEST.pdb',
#                     '--log', './temp_dir/DOCKING_TEST_log.txt'], capture_output=True)

# A = subprocess.run(['./executables/smina', '-r', './protein/prot.pdb',
#                     '-l', './test.sdf','--exhaustiveness', '10', 
#                     '--autobox_ligand', './lig/lig.pdb', '--flexdist', '2', 'flexdist_ligand', './test.sdf', '--autobox_add', '3', '-o', './temp_dir/DOCKING_TEST.pdb',
#                     '--log', './temp_dir/DOCKING_TEST_log.txt'], capture_output=True)


#os.system('smina -r ./temp_dir/rec.pdb -l ./test.sdf --cpu {} --exhaustiveness 10 --autobox_ligand ./temp_dir/lig.pdb --autobox_add 3 -o ./temp_dir/DOCKING_TEST.pdb --log ./temp_dir/DOCKING_TEST_log.txt'.format(cpu_count))
