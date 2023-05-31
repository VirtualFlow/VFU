# VF-Unity
Streamlined version of VirtualFlow combining both VFVS and VFLP

## Prerequisites
Please clone the repository using: 
```
git clone git@github.com:VirtualFlow/VF-Unity.git
```
Please ensure that the following packages are installed: 
- [RDKit version 2021.09.5](https://www.rdkit.org/docs/Install.html)
- [Open Babel 3.1.0](https://openbabel.org/docs/dev/Installation/install.html)
- [Python 3.7.13](https://www.python.org/downloads/)
- [SELFIES](https://github.com/aspuru-guzik-group/selfies)


## File Navigator
* `run_vf_unity.py`: Main file that initiates docking calculations. 
* `initiate_calc.py`: File that initiates scoring, docking calculations. 
* `lig_process.py`: Process provided ligand into 3D format compatible with docking program (used by )
* `pose_prediction.py`: File for running pose prediction on processed ligands (used by pose_prediction). 
* `scoring_functions.py`: File for running scoring on already docked ligands (used by initiate_calc.py). 
* `config.txt`: Config file concisting of user definable parameters for running calculation. 
* `/ligands/`: Directory created by VF-Unity which will contain all processed ligands in a ready-to-dock format. 
* `/outputs/`: Directory created by VF-Unity which will contain docked ligand files. 



## Quick Start (Using a config.txt file)
We will be running QuickVina on a processed protein located in the config directory (`5wiu_test.pdbqt`). We provide a config.txt file which contains all parmeters for this simple run. 
Please edit this file based on your preferance: 
```
# The choice of the docking method
# Possible choices:
# AutoDock-Koto, AutodockVina_1.2, AutodockZN, EquiBind, FRED
# FitDock, GalaxyDock3, LigandFit, LightDock, M-Dock
# MCDock, MM-GBSA, PLANTS, PSOVina, RLDock
# SEED, adfr, autodock_cpu, autodock_gpu, AutodockVina_1.1.2
# dock6, flexx, gnina, gnina-scoring
# gold, gwovina, iGemDock, idock, ledock
# molegro, nnscore2, qvina, qvina-w, rDock
# rf-score, rosetta-ligand, smina, smina-scoring, vina
# vina_carb, vina_xb, GlideSP, GlideXP, GlideHTVS
# qvina_gpu, qvina_w_gpu, vina_gpu, vina_gpu_2.0
# The choice of the scoring function
# Possible choices:
# nnscore2, rf-score, smina-scoring, ad4_scoring, vinandro_scoring
# vina_scoring, gnina_scoring, chemplp_scoring, PLP_scoring, PLP95_scoring
# contact_scoring, continuous_scoring, grid_scoring, 
# asp, chemscore, goldscore, plp, mm_gbsa_scoring, Hawkins_gbsa
# Please note: different pose prediction/docking methods can be combined with scoring functions.
# For example: ’qvina+nnscore2’.
# For supported choices/combinations please see the VirtualFlow homepage.

program_choice=qvina+nnscore2

# The x,y&z coordinates of the center of the docking space. The binding space describes the location where a molecule
# is allowed to bind.
center_x=-17.820
center_y=16.140
center_z=-18.643

# The size (in Angstroms) of the docking space in the x,y&z directions.
size_x=20
size_y=20
size_z=20

# How many poses to search for, providing a limit to the maximum number of iterations that a docking program performs in
# search for good poses
# Large exhaustiveness settings lead to increased computational costs.
exhaustiveness=10 


# Molecule (either a string in smiles, selfies or amino-acid sequence) to be used for docking
# If the is_selfies=True, or is_peptide=True, a conversion from selfies->smiles
# and aa-sequence->smiles is performed. 
smi=C1=CC(=CC=C1CSCC2C(C(C(O2)N3C=NC4=C(N=CN=C43)N)O)O)Cl 
is_selfies=False
is_peptide=False


# Location to the prepared receptor file
#  The receptor needs to be in the correct format, supported by the user's selected docking program.
#  Additionally, the file needs to be present in the config directory.
receptor=./config/5wiu_test.pdbqt
```
To execute the program, please run: 
```
python3 run_vf_unity.py
```
We note: 
1. The processed ligands will be located within the newly created ligands directory.
2. The docked output from running QuickVina will be located in the newly created outputs directory.
3. The default behaviour is for the program (QuickVina) is to make use of all available CPUs. 
4. A summary csv file `docking_output.csv` is created (for running QuickVina): 
    ```
    Ligand File,Docking Values,Docking Pose
    7.pdbqt,"-9.9,-9.8,-9.1,-9.0,-9.0,-8.8,-8.8,-8.7,-8.6",./outputs/pose_7.pdbqt
    14.pdbqt,"-9.6,-9.1,-9.0,-8.8,-8.7,-8.6,-8.6,-8.4,-8.4",./outputs/pose_14.pdbqt
    6.pdbqt,"-9.8,-9.7,-9.2,-8.9,-8.8,-8.8,-8.7,-8.7,-8.5",./outputs/pose_6.pdbqt
    9.pdbqt,"-9.3,-9.3,-9.3,-9.1,-9.1,-9.0,-9.0,-8.9,-8.8",./outputs/pose_9.pdbqt
    12.pdbqt,"-8.9,-8.9,-8.8,-8.5,-8.5,-8.5,-8.5,-8.1,-8.0",./outputs/pose_12.pdbqt
    4.pdbqt,"-9.6,-9.3,-9.2,-9.2,-9.0,-9.0,-8.8,-8.7,-8.7",./outputs/pose_4.pdbqt
    3.pdbqt,"-9.9,-9.2,-9.0,-8.9,-8.8,-8.7,-8.6,-8.5,-8.4",./outputs/pose_3.pdbqt
    1.pdbqt,"-9.4,-9.4,-9.1,-9.0,-8.9,-8.8,-8.7,-8.7,-8.7",./outputs/pose_1.pdbqt
    5.pdbqt,"-9.8,-9.5,-9.4,-9.3,-9.3,-9.2,-9.1,-9.1,-9.0",./outputs/pose_5.pdbqt
    10.pdbqt,"-9.6,-9.5,-9.4,-9.3,-9.2,-9.2,-9.1,-8.9,-8.8",./outputs/pose_10.pdbqt
    2.pdbqt,"-9.6,-9.3,-8.9,-8.9,-8.8,-8.8,-8.8,-8.6,-8.5",./outputs/pose_2.pdbqt
    15.pdbqt,"-9.8,-9.3,-9.2,-9.2,-9.1,-9.0,-9.0,-8.8,-8.8",./outputs/pose_15.pdbqt
    13.pdbqt,"-9.1,-9.1,-9.0,-8.8,-8.7,-8.7,-8.7,-8.7,-8.7",./outputs/pose_13.pdbqt
    8.pdbqt,"-9.1,-9.1,-9.1,-8.8,-8.7,-8.7,-8.6,-8.6,-8.5",./outputs/pose_8.pdbqt
    11.pdbqt,"-9.7,-9.5,-9.3,-9.2,-9.2,-9.1,-9.0,-8.9,-8.8",./outputs/pose_11.pdbqt
    0.pdbqt,"-9.6,-9.3,-9.3,-9.0,-8.9,-8.9,-8.7,-8.7,-8.5",./outputs/pose_0.pdbqt
    ```
5. Sepperately, a summary csv file is created (`rescoring_output.csv`) is created when running rescoring (nnscore2): 
    ```
    Docked Ligand,Re-scored Value
    ./outputs/pose_11.pdbqt,Kd = 114.94 fM;Kd = 2.29 pM;Kd = 3.22 fM;Kd = 10.3 fM;Kd = 11.92 fM;Kd = 0.87 fM;Kd = 307.09 fM;Kd = 2.18 pM;Kd = 22.69 pM
    ./outputs/pose_10.pdbqt,Kd = 14.22 fM;Kd = 13.13 fM;Kd = 10.96 pM;Kd = 29.48 fM;Kd = 73.28 fM;Kd = 2.35 fM;Kd = 108.96 fM;Kd = 122.82 fM;Kd = 0.35 fM
    ./outputs/pose_9.pdbqt,Kd = 216.4 fM;Kd = 685.98 fM;Kd = 155.95 fM;Kd = 116.24 fM;Kd = 72.66 fM;Kd = 189.04 pM;Kd = 4.74 pM;Kd = 7.01 pM;Kd = 1.04 nM
    ./outputs/pose_5.pdbqt,Kd = 0.14 fM;Kd = 3.58 pM;Kd = 2.75 fM;Kd = 23.24 fM;Kd = 7.79 pM;Kd = 152.06 nM;Kd = 0.12 fM;Kd = 121.77 fM;Kd = 18.27 pM
    ./outputs/pose_12.pdbqt,Kd = 0.16 fM;Kd = 11.43 fM;Kd = 1.76 pM;Kd = 4.74 pM;Kd = 2.15 pM;Kd = 4.79 pM;Kd = 364.24 fM;Kd = 55.44 fM;Kd = 20.77 fM
    ./outputs/pose_7.pdbqt,Kd = 101.4 fM;Kd = 0.12 fM;Kd = 1.14 fM;Kd = 134.77 fM;Kd = 0.55 fM;Kd = 10.53 nM;Kd = 0.0 fM;Kd = 129.84 fM;Kd = 8.86 nM
    ./outputs/pose_13.pdbqt,Kd = 0.34 fM;Kd = 0.13 fM;Kd = 1.74 pM;Kd = 73.35 fM;Kd = 203.96 fM;Kd = 11.28 pM;Kd = 3.75 pM;Kd = 85.11 fM;Kd = 5.17 pM
    ./outputs/pose_0.pdbqt,Kd = 19.81 fM;Kd = 0.01 fM;Kd = 0.06 fM;Kd = 14.96 fM;Kd = 0.1 fM;Kd = 578.15 fM;Kd = 120.14 fM;Kd = 24.31 fM;Kd = 3.01 pM
    ./outputs/pose_4.pdbqt,Kd = 12.63 fM;Kd = 2.81 pM;Kd = 0.0 fM;Kd = 52.42 fM;Kd = 0.96 fM;Kd = 310.52 fM;Kd = 36.57 pM;Kd = 17.94 fM;Kd = 0.13 fM
    ./outputs/pose_6.pdbqt,Kd = 0.91 fM;Kd = 8.34 pM;Kd = 0.14 fM;Kd = 0.19 fM;Kd = 0.02 fM;Kd = 0.0 fM;Kd = 16.4 fM;Kd = 60.84 fM;Kd = 1.34 pM
    ./outputs/pose_8.pdbqt,Kd = 94.53 fM;Kd = 23.83 pM;Kd = 6.22 pM;Kd = 1.34 fM;Kd = 276.47 fM;Kd = 0.63 fM;Kd = 71.78 fM;Kd = 5.64 fM;Kd = 18.55 pM
    ./outputs/pose_1.pdbqt,Kd = 0.13 fM;Kd = 0.09 fM;Kd = 11.47 pM;Kd = 11.69 fM;Kd = 5.36 pM;Kd = 8.35 pM;Kd = 2.44 pM;Kd = 4.19 pM;Kd = 23.34 pM
    ./outputs/pose_3.pdbqt,Kd = 0.0 fM;Kd = 1.65 fM;Kd = 14.99 fM;Kd = 0.24 fM;Kd = 213.73 pM;Kd = 6.34 pM;Kd = 472.1 fM;Kd = 329.75 pM;Kd = 261.08 fM
    ./outputs/pose_15.pdbqt,Kd = 1.06 pM;Kd = 136.55 fM;Kd = 1.82 pM;Kd = 106.86 fM;Kd = 4.26 pM;Kd = 1.55 pM;Kd = 6.23 pM;Kd = 61.89 fM;Kd = 242.37 pM
    ./outputs/pose_2.pdbqt,Kd = 11.5 fM;Kd = 18.0 fM;Kd = 2.34 pM;Kd = 783.94 fM;Kd = 202.75 fM;Kd = 79.9 fM;Kd = 210.65 fM;Kd = 13.13 fM;Kd = 511.8 fM
    ./outputs/pose_14.pdbqt,Kd = 162.12 fM;Kd = 726.22 fM;Kd = 10.38 fM;Kd = 432.15 fM;Kd = 0.06 fM;Kd = 9.24 pM;Kd = 3.07 pM;Kd = 23.07 nM;Kd = 7.97 pM
    ```
    
## Quick Start (Using a python function call)
```
from run_vf_unity import main 

program_choice   = 'qvina'
scoring_function = 'nnscore2' 
center_x         = -17.820
center_y         = 16.140
center_z         = -18.643
size_x           = 20
size_y           = 20 
size_z           = 20
exhaustiveness   = 10
smi              = 'C1=CC(=CC=C1CSCC2C(C(C(O2)N3C=NC4=C(N=CN=C43)N)O)O)Cl'
is_selfies       = False
is_peptide       = False
receptor         = './config/5wiu_test.pdbqt'

pose_pred_out, re_scored_values = main(program_choice, scoring_function, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor)
```
We note: 
1. The processed ligands will be located within the newly created ligands directory.
2. The docked output from running QuickVina will be located in the newly created outputs directory.
3. The default behaviour is for the program (QuickVina) is to make use of all available CPUs. 
4. The output from running th QuickVina calculation will be stored in the dictionary pose_pred_out.
5. The output from running th QuickVina NNScore2.0 will be stored in the dictionary re_scored_values.
6. No output csv files are created in this case. 

## Running calculations for multiple molecules: 
Please prepare a file concisting of a list of smiles. For example, `molecules.txt`: 
```
Index,Smiles
0,CCCCCCCCCC
1,C1=CC(=CC=C1CSCC2C(C(C(O2)N3C=NC4=C(N=CN=C43)N)O)O)Cl
``` 
The molecules can be run simply using out python function call: 
```
import os 
from run_vf_unity import main 

program_choice   = 'qvina'
scoring_function = '' 
center_x         = -17.820
center_y         = 16.140
center_z         = -18.643
size_x           = 20
size_y           = 20 
size_z           = 20
exhaustiveness   = 10
is_selfies       = False
is_peptide       = False
receptor         = './config/5wiu_test.pdbqt'

with open('./molecules.smi', 'r') as f: 
    lines = f.readlines()
lines = lines[1: ]

for item in lines: 
    idx,smi = item.split(',')
    pose_pred_out, re_scored_values = main(program_choice, scoring_function, center_x, center_y, center_z, size_x, size_y, size_z, exhaustiveness, smi, is_selfies, is_peptide, receptor)
    os.system('rm -rf ligands')
    os.system('cp -a outputs outputs_{}'.format(idx))
```
The corresponding index (column 1) of a molecule in the `molecules.txt` file will be used for storing the results in `outputs_*index*`

## Special Considerations
### Using AutoDock-GPU/CPU
Please compile the code using instructions from: [https://github.com/ccsb-scripps/AutoDock-GPU](https://github.com/ccsb-scripps/AutoDock-GPU). 
After successfull compilation, within the bin directory, an executable will be made (example name: `autodock_gpu_1wi`). Then, the code is ready to run. 
We provide an example inside `./executables/vf_gpu_example.zip`. Inside the directory, a prepared protein-ligand pair is provided and the code can be run 
using: `./autodock_gpu_1wi --ffile 1stp_protein.maps.fld --lfile ./1stp_ligand.pdbqt`


### Using EquiBind
Please download the code using instructions from: [https://github.com/HannesStark/EquiBind](https://github.com/HannesStark/EquiBind). 
Please create a conda enviroment per the instructions of the EquiBind repository. 
Copy paste all files inside the working directory of VF-Unity. 


### Using rDock
Please compile the code using the instruction provided in: [https://rdock.sourceforge.net/installation/](https://rdock.sourceforge.net/installation/). 
We installed rDock using anacond (with `conda install -c bioconda rdock`)

### Using MM-GBSA
Please install AmberTools: [https://ambermd.org/GetAmber.php#ambertools](https://ambermd.org/GetAmber.php#ambertools). We managed performed the download using conda. 
Additionally, please note: the variable `chimera_path` should be updated to location of the Chimera on your system. Chimera can be downloaded: [https://www.cgl.ucsf.edu/chimera/download.html](https://www.cgl.ucsf.edu/chimera/download.html).


### AutoDockZN
Requirments
Please install the ADFR software suite: [https://ccsb.scripps.edu/adfr/downloads/](https://ccsb.scripps.edu/adfr/downloads/). 
Please install meeko with: `pip install meeko`. 
Instructions for running AutoDockZN can be found on their official website: [https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html](https://autodock-vina.readthedocs.io/en/latest/docking_zinc.html).


### AutoDockZN
Requirments
Please install the OpenEye software suite: [https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html](https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html). 

### SEED
Requirments
Please install AmberTools: [https://ambermd.org/GetAmber.php#ambertools](https://ambermd.org/GetAmber.php#ambertools).

### RF-score
Please paste the executable from [https://github.com/oddt/rfscorevs_binary](https://github.com/oddt/rfscorevs_binary) in the executables directory. 

### qvina_gpu
A qvina_gpu executable (of name 'qvina_gpu') should be placed in directory: /executables
Instructions for compilation are provided in https://github.com/DeltaGroupNJUPT/QuickVina2-GPU

### qvina_w_gpu
A qvina_gpu executable (of name 'qvina_w_gpu') should be placed in directory: /executables
Instructions for compilation are provided in https://github.com/DeltaGroupNJUPT/QVina-W-GPU

### vina_gpu
A qvina_gpu executable (of name 'vina_gpu') should be placed in directory: /executables
Instructions for compilation are provided in https://github.com/DeltaGroupNJUPT/Vina-GPU

### vina_gpu_2.0
A qvina_gpu executable (of name 'vina_gpu_2.0') should be placed in directory: /executables
Instructions for compilation are provided in https://github.com/DeltaGroupNJUPT/Vina-GPU-2.0

### ADFR
An adfr executable (of name adfr) should be placed in directory: /executables directory. An executable needs to be compiled based on a user’s system using instructions described in [https://ccsb.scripps.edu/adfr/downloads/](https://ccsb.scripps.edu/adfr/downloads/).

### HDock
An HDock executable (of name hdock) should be placed in directory: /executables directory. 
An createpl executable (of name createpl) should be placed in directory: /executables directory. 


### FRED
Requirments
Please install OpenEye: [https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html](https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html).
Additionally, a valid OpenEye licence is required. Namely, a file named oe_license.txt needs to be placed in the working directory. 

### GalaxyDock3
Please paste the executable from [https://galaxy.seoklab.org/files/by2hsnvxjf/softwares/galaxydock.html](https://galaxy.seoklab.org/files/by2hsnvxjf/softwares/galaxydock.html) in the executables directory. 
Note a data directory (name 'data') is required for successful runs. 

### LightDock
We suggest downloading LightDock from: [https://github.com/lightdock/lightdock](https://github.com/lightdock/lightdock). For this, within the current working directory (VF-Unity), please run: 
```
git clone https://github.com/lightdock/lightdock.git
virtualenv venv
source venv/bin/activate
cd lightdock
pip install -e .
```

### MpSDockZN
Please add the executable (of name `MpSDock`) in executables directory.
Please install AmberTools: [https://ambermd.org/GetAmber.php#ambertools](https://ambermd.org/GetAmber.php#ambertools).
We managed performed the download using conda. 
Please note: the variable `chimera_path` should be updated to location of the Chimera on your system. Chimera can be downloaded: [https://www.cgl.ucsf.edu/chimera/download.html](https://www.cgl.ucsf.edu/chimera/download.html).
A working dock6 download (with a valid licence) is required. The variable `dock6_path` should be updated to location of the Chimera on your system.
A `box.in` input file is required for the program. Please specify the path in variable `box_in_file`.
A `grid.in` input file is required for the program. Please specify the path in variable `grid_in_file`.
A `dock.in` input file is required for the program. Please specify the path in variable `dock_in_file`.

### Running with CovDock
A valid Schrödinger license is required to run CovDock.

### Running with GlideSP/XP/HTVS
A valid Schrödinger license is required to run GlideSP/XP/HTVS.

### Re-scoring with GOLD fitness functions
A valid GOLD license is required to re-score docking results with GOLD fitness functions.  
An executable named `gold_auto` should be placed in the `/executables` directory.

When specifying a GOLD fitness function with which to re-score docking results, a dynamically loadable 
shared object library must be available in the top-level directory of this repository and named correspondingly 
with the parameter defined in the VFU configuration, i.e., one of `plp`, `asp`, `chemscore`, `goldscore`.

### Contributing
If you are interested in contributing to VirtualFlow, whether it is to report a bug or to extend VirtualFlow with your own code, please see the file [CONTRIBUTING.md](CONTRIBUTING.md) and the file [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).



### License
The project ist distributed under the GNU GPL v2.0. Please see the file [LICENSE](LICENSE) for more details. 


### Citation
Gorgulla, Christoph, et al. "VirtualFlow 2.0-The Next Generation Drug Discovery Platform Enabling Adaptive Screens of 69 Billion Molecules." bioRxiv (2023): 2023-04.
