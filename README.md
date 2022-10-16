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



## Running a single example
For running a single example (a docking score of a single SMILES/SELFIES), we will be running `run_single.py`. Specifically, the following steps can be followed: 
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
Additionally, please note: the variable `chimera_path` should be updated to location of the Chimera on your system. Chemera can be downloaded: [https://www.cgl.ucsf.edu/chimera/download.html](https://www.cgl.ucsf.edu/chimera/download.html).




## Questions, problems?
Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (TODO)

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
