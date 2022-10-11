# VF-Unity
Streamlined version of VirtualFlow combining both VFVS and VFLP


## Prerequisites
Please clone the repository using: 
```
git clone git@github.com:VirtualFlow/VF-Unity.git
```



## Special Considerations
### Using AutoDock-GPU
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


## Questions, problems?
Make a github issue 😄. Please be as clear and descriptive as possible. Please feel free to reach
out in person: (TODO)

## License

[Apache License 2.0](https://choosealicense.com/licenses/apache-2.0/)
