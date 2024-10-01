# aBPA: ancient Bacterial Pangenome Analysis

ancient Bacterial Pangenome Analysis is a `nextflow` package.

# 1/ Installation


First you need to have `nextflow` and `conda/mamba` installed. Please visit: https://www.nextflow.io/docs/latest/install.html and follow the instructions. In the case of conda, if you don't have any conda version installed I would reccomend miniforge (https://github.com/conda-forge/miniforge).


After installing `nextflow`, if you wish to be able to call the executable without explicitly define the full PATH while not having to use sudo:


### Open your bashrc file:


>`nano .bashrc`


### Then add this line at the end: 


>`export PATH=$PATH:/home/users/myuser/Softwares/`

### Finally refresh the environment to activate the changes


>`source ~/.bashrc`


### Note that the PATH will be unique in your system. In my case I installed nextflow in the folder /Softwares/


### Then git clone this repository:



>`git clone https://github.com/mudymudy/aBPA/`


After downloading the repository you should see the folders `bin/` `config/` `envs/` and the files `aBPA.nf` and `nextflow.config`.



# 2/ First steps
