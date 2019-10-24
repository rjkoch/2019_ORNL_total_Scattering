# 2019 ORNL Total Scattering School, Diffpy-cmi demonstrations

This contains several demonstrations of Diffpy-CMI.

This project is both a binder project, and a set of pure python files.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/rjkoch/2019_ORNL_total_scattering/master)


For running the pure python files locally on your machine, you will f
first need to install the Diffpy-CMI packages.

## For installation of Diffpy-CMI in Unix based (Linux and Mac OS) systems

1. Install some version of Conda Python 3:
    - Miniconda is recommended, and can be found here:
    
        https://docs.conda.io/en/latest/miniconda.html
2. Install Diffpy-CMI
    - Open a Unix terminal
    - Create a new environment, and install Diffpy-CMI wihin it, using the following command:
            
         `conda create -n diffpy --override-channels --channel defaults --channel diffpy diffpy-cmi`
        - Enter `y` when prompted to install the relevant packages.

    - Activate your new environment with
        
        `conda activate diffpy`

3. Download and unzip these tutorial files to a meaningful directory.
    - Navigate to wherever you have placed the files using `cd /path/to/files`

## For installation of Diffpy-CMI in Windows systems

***Only Windows 10 is currently supported, and Administrator rights are required.***

1. Install the Windows Subsystem for Linux (WSL)
    - Instructions for this can be found here:

        https://docs.microsoft.com/en-us/windows/wsl/install-win10
    - You can choose which version of Linux you'd like to use, but Ubuntu is recomended.

2. Install an Xserver for plotting data.
    - VcXsrv is a good option, and can be found here:

        https://sourceforge.net/projects/vcxsrv/

3. Start your Xserver
    - Left click the Windows start button (Windows flag in lower left corner)
    - Type "xlaunch"
    - Right click on the xlaunch icon and select "Run as Administrator"
    - In the dialog window that opens:
        - Select "Multiple windows," enter -1 in the Display number field, and click "Next"
        - Select "Start no client," and click "Next"
        - Do not change any settings here, and click "Next"
        - If you would like your Xserver to start every time Windows boots up:
            - Click "Save configuration"
            - Navigate to `C:\ProgramData\Microsoft\Windows\Start Menu\Programs\StartUp`
            - Save the file as config.xlaunch
            - If you choose not to save a config.xlaunch file here, you will need
            to start your xserver after each reboot manually.
        - Click "Finish"

4. Tell WSL how to connect to your Xserver
    - Open a WSL command prompt
        - If you installed Ubuntu, click the Windows start button, type Ubuntu, and press enter.
    - Enter the following command:
        
        `echo 'export DISPLAY=localhost:0' >> ~/.bashrc `
    
5. Install some version of Conda Python 3:
    - Miniconda is recommended, and can be found here:
    
        https://docs.conda.io/en/latest/miniconda.html
    - ***Be sure to download the Linux (not Windows) version of 
    whatever conda python you choose to use!*** 

6. Install Diffpy-CMI
    - Open a WSL terminal
    - Create a new environment, and install Diffpy-CMI wihin it, using the following command:
            
         `conda create -n diffpy --override-channels --channel defaults --channel diffpy diffpy-cmi`
        - Enter `y` when prompted to install the relevant packages.

    - Activate your new environment with
        
        `conda activate diffpy`

7. Download and unzip these tutorial files to a meaningful directory.
    - Navigate to wherever you have placed the files using `cd /path/to/files`
    
    
For older versions of Windows, a Linux virtual machine (VM) can be created.
VirtualBox is a good, free, option for doing this, and can be found here:
https://www.virtualbox.org/

Once this is installed, and the VM running, 
follow the instructions above for Unix based installation.

## Running the tutorials
- Tutorials are meant to be run in order
- The relevant files are python files, of the form *.py
- Files contain extensive comments.
- Open a Unix or WSL terminal, and activate the relevant environment with
    
    `conda activate diffpy` 
- Navigate to a tutorial directory, for example, type `cd ./tutorial_1` for tutorial 1
- A tutorial can be run  in one of two ways:
    - Entering `python tutorial_1.py` directly in your terminal
    - In an IPython session:
        - Enter `ipython` in your terminal. This will start an interactive python interpreter.
        - Enter `%run tutorial_1.py`
        
        
***When in doubt, ask!***
