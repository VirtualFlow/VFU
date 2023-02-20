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



## Quick Start (Running a docking calculation)
We will be running QuickVina on a processed protein located in the config directory (`5wiu_test.pdbqt`). We provide a config.txt file which contains all parmeters for this simple run. 
Please edit this file based on your preferance: 
```
# The choice of the docking method
# Possible choices:
# AutoDock-Koto, AutodockVina_1.2, AutodockZN, EquiBind, FRED
# FitDock, GalaxyDock3, LigandFit, LightDock, M-Dock
# MCDock, MM-GBSA, PLANTS, PSOVina, RLDock
# SEED, adfr, autodock_cpu, autodock_gpu, autodock_vina
# dock6, flexx, glide, gnina, gnina-scoring
# gold, gwovina, iGemDock, idock, ledock
# molegro, nnscore2, qvina, qvina-w, rDock
# rf-score, rosetta-ligand, smina, smina-scoring, vina
# vina_carb, vina_xb
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
    ./outputs/pose_11.pdbqt,Kd = 1.0 pM Kd = 672.97 fM Kd = 8.17 pM Kd = 349.7 fM Kd = 4.05 pM Kd = 2.11 pM Kd = 2.0 nM Kd = 66.72 fM Kd = 1.15 fM
    ./outputs/pose_10.pdbqt,Kd = 11.52 fM Kd = 11.49 fM Kd = 266.9 fM Kd = 11.81 fM Kd = 0.0 fM Kd = 7.63 pM Kd = 115.91 fM Kd = 773.49 fM Kd = 1.79 pM
    ./outputs/pose_9.pdbqt,Kd = 30.43 fM Kd = 22.13 pM Kd = 48.66 fM Kd = 1.03 pM Kd = 14.35 pM Kd = 103.81 fM Kd = 10.97 pM Kd = 0.14 fM Kd = 0.18 fM
    ./outputs/pose_5.pdbqt,Kd = 0.19 fM Kd = 0.12 fM Kd = 55.87 fM Kd = 6.36 pM Kd = 2.96 pM Kd = 0.0 fM Kd = 4.77 fM Kd = 7.49 pM Kd = 44.51 pM
    ./outputs/pose_12.pdbqt,Kd = 1.79 fM Kd = 1.5 fM Kd = 2.2 pM Kd = 192.03 pM Kd = 34.51 fM Kd = 15.97 pM Kd = 1.33 pM Kd = 45.29 fM Kd = 9.07 fM
    ./outputs/pose_7.pdbqt,Kd = 279.09 fM Kd = 0.12 fM Kd = 7.22 pM Kd = 9.8 fM Kd = 2.08 fM Kd = 0.01 fM Kd = 0.0 fM Kd = 0.01 fM Kd = 20.22 pM
    ./outputs/pose_13.pdbqt,Kd = 1.0 fM Kd = 0.13 fM Kd = 4.18 fM Kd = 41.14 fM Kd = 560.85 fM Kd = 142.89 fM Kd = 3.61 fM Kd = 376.98 fM Kd = 0.02 fM
    ./outputs/pose_0.pdbqt,Kd = 4.52 fM Kd = 0.27 fM Kd = 3.81 pM Kd = 1.4 fM Kd = 14.76 fM Kd = 20.42 pM Kd = 0.0 fM Kd = 95.3 pM Kd = 370.42 pM
    ./outputs/pose_4.pdbqt,Kd = 22.57 fM Kd = 528.16 fM Kd = 0.01 fM Kd = 0.0 fM Kd = 1.16 fM Kd = 7.17 pM Kd = 14.49 fM Kd = 379.8 pM Kd = 35.22 fM
    ./outputs/pose_6.pdbqt,Kd = 0.5 fM Kd = 20.04 pM Kd = 0.17 fM Kd = 0.05 fM Kd = 10.83 fM Kd = 481.17 fM Kd = 9.3 pM Kd = 67.82 pM Kd = 3.61 pM
    ./outputs/pose_8.pdbqt,Kd = 389.19 fM Kd = 0.68 fM Kd = 5.92 pM Kd = 0.0 fM Kd = 0.13 fM Kd = 62.06 nM Kd = 2.1 pM Kd = 2.79 fM Kd = 96.16 pM
    ./outputs/pose_1.pdbqt,Kd = 0.13 fM Kd = 0.11 fM Kd = 17.96 pM Kd = 73.79 fM Kd = 34.98 fM Kd = 296.29 fM Kd = 4.31 pM Kd = 18.54 fM Kd = 0.25 fM
    ./outputs/pose_3.pdbqt,Kd = 0.16 fM Kd = 50.41 fM Kd = 1.56 fM Kd = 0.18 fM Kd = 398.94 pM Kd = 12.01 fM Kd = 152.27 pM Kd = 1.19 fM Kd = 0.0 fM
    ./outputs/pose_15.pdbqt,Kd = 181.12 fM Kd = 1.11 fM Kd = 0.01 fM Kd = 24.82 pM Kd = 0.8 fM Kd = 54.8 pM Kd = 489.94 pM Kd = 0.0 fM Kd = 13.35 fM
    ./outputs/pose_2.pdbqt,Kd = 1.3 fM Kd = 18.63 fM Kd = 0.0 fM Kd = 19.07 fM Kd = 32.29 pM Kd = 38.18 pM Kd = 77.43 fM Kd = 830.99 fM Kd = 0.02 fM
    ./outputs/pose_14.pdbqt,Kd = 11.13 pM Kd = 5.51 pM Kd = 48.55 fM Kd = 17.8 fM Kd = 20.05 pM Kd = 973.78 fM Kd = 45.21 pM Kd = 856.73 fM Kd = 7.61 pM
     
    ```

## Running in batch 

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

### FRED
Requirments
Please install OpenEye: [https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html](https://docs.eyesopen.com/toolkits/python/quickstart-python/install.html).
Additionally, a valid OpenEye licence is required. Namely, a file named oe_license.txt needs to be placed in the working directory. 

### GalaxyDock3
Please paste the executable from [https://galaxy.seoklab.org/files/by2hsnvxjf/softwares/galaxydock.html](https://galaxy.seoklab.org/files/by2hsnvxjf/softwares/galaxydock.html) in the executables directory. 
Note a data directory (name 'data') is required for successful runs. 

### SMINA Scoring
Please not that the ligand to be scored must be the sole ligand in the file. 

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

### Contributing
If you are interested in contributing to VirtualFlow, whether it is to report a bug or to extend VirtualFlow with your own code, please see the file [CONTRIBUTING.md](CONTRIBUTING.md) and the file [CODE_OF_CONDUCT.md](CODE_OF_CONDUCT.md).



### License
The project ist distributed under the GNU GPL v2.0. Please see the file [LICENSE](LICENSE) for more details. 


### Citation
TODO