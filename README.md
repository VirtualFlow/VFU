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



## Quick Start (Running a single docking calculation)
For running a single example (a docking score of a single SMILES/SELFIES), we will be running `run_vf_unity.py`. Specifically, the following steps can be followed: 
1. Inside the `config` directory, please add a protein file, and a ligand file. We will be performing this demonstration using the protein `prot_1.pdb`, already pasted inside the `config` directory.
2. Please open the file `run_single.py`. At the top of the file are all the parameters for running the script. The parameters are, 
   - is_selfies     = False 
   - program_choice = 'smina' 
   - receptor       = './config/prot_1.pdb'
   - smi            = 'BrC=CC1OC(C2)(F)C2(Cl)C1.CC.[Cl][Cl]'
   - Parameters for executing the docking pose search: 
     - exhaustiveness = 10
     - center_x       = -16  
     - center_y       = 145  
     - center_z       = 27   
     - size_x         = 10  
     - size_y         = 10   
     - size_z         = 10   
3. Please do not worry. We will be describing these paramters one at a time! :)
4. We begin by looking at the paramter `smi`. The variable needs to be set to a valid smile string. To process the SMILE string, the program uses RdKit to (1) desaults, (2) neutralizes and (3) enumerates sterio-isomers. Subsequently, all molecules will be converted to 3D using OpenBabel. 
5. `is_selfies`: If set to True, the program expects the variable `smi` to contain a valied SELFIES string. Please ensure to install SELFIES using `pip install selfies`. 
6. `receptor`: The name of the receptor file. Note: it is assumed that this file name is located within the config directory. 
7. `program_choice`: We currently support the docking programs `smina,qvina,qvina-w,vina,vina_carb,vina_xb,gwovina,PLANTS,autodock_gpu,autodock_cpu,EquiBind,rDock,gnina,ledock,idock,autodock_vina,adfr`. The variable `program_choice` can be set to any of these values. <br />Note: There are special instructions for AutoDock-GPU(CPU),EquiBind & rDock. Please have a look at the Special Considerations section below. 
8. Finally, there are a few parameters for the docking software: `exhaustiveness, center_x/y/z, size_x/y/z`. The exhaustive parameter describes how many poses to search for (we have set it to 10 for this example). center_x/y/z and size_x/y/z describe where the docking is performed in the protein (i.e., the binding spot for the ligand).
9. Alright! I hope the parameters make sense, along with the example we will be running! You are ready to run `python3 run_single.py`
10. The output poses for the docking calculation are saved within the directory `outputs`. The dictionay `results` (from running `run_single.py`) contains the ligand file, the docked pose file name, along with the docking score. 


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