# Mathematical modeling of wave processes in a rope system
This is official repository of the following research paper:   


This repository contains the Python code of proposed 2D rope model and some results of experiments.    


## Program usage instructions
You have several options how to use builded program.

1. If you are a **Windows** OS user you can download the dist folder from this repository and run the file *rope_system.exe*  by double clicking or from the command line (which is more convenient and recommended). 
You can enter the following command in your cmd in the folder which contains this file: `rope_system.exe --help` (or `-h`) and you will get a full list of parameters which you can configure in this program.  

Please, don't delete default_model_config.json file. It contains default system setup data. *rope_system.exe* is builded from the *main.py* file from this repository, so is a safe exe.

2. If you are an advanced user who knows how to read a Python code and you have any IDE or text editor for running it, follow the next steps: 
    1. Clone this repository, using `git clone` command. 
    2. Create a Python virtual environment with Python version 3.12+ and libraries from *requirements.txt*.
    3. Run (and modify) code from *main.py* file. 

3. Docker container for running program on Linux or MacOS will be added later.  