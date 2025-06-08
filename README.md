# 1. Phenalyzer Workflow Overview
Phenalyzer is a pipeline to process and analyze raw data PhenoCycler qptiff files. The workflow consists of two steps with their separate software packages:
- **Step 1: Segment TMA cores in QuPath using CellPose:**
  - Load qptiff files and dearray TMA into cores
  - Run CellPose groovy script 
      - Combine multiple membrane markers for boundary stain
      - Nucleated and anucleated cell profiles are stored
      - Automatically calculate marker-wise mean/median/min/max of each cell's nucleus/cytoplasma/membrane/entire area
  - Optional: run object classification 
  - Export tsv file of segmented cells
- **Step 2: Phenalyzer Data Processing in R (VS Code):**
  - Import tsv files and map metadata
  - Gate cells into subsets
  - Clustering/Heatmap of each subsets
  - Collapse and/or annotate clusters
  - Optional: load spatial proximity masks




# 2. Phenalyzer QuPath CellPose Setup Guide
### 2.1 Install Programs
- **Install python**
  - Download
- **Install CellPose**
  - Download
- **Install QuPath**
  - Download






# 3. Phenalyzer Software R Setup Guide


### 3.1 Windows Preparation

- **Enable Required Features:**
  - Go to Windows Settings → System → Optional Features and ensure **WMIC** is installed.
  - Go to "More Windows Features" and enable:
    - `Virtual Machine Platform`
    - `Windows Subsystem for Linux`

- **Install Visual Studio Code (VS Code):**
  - Download and install from [https://code.visualstudio.com/](https://code.visualstudio.com/).

- **Install WSL and Ubuntu:**
  - Open PowerShell as Administrator and run:
    ```powershell
    wsl --install
    ```
  - It will download Ubuntu by default. Follow prompts and create a user account (no spaces in username).

- **Restart your PC.**

### 3.2 VS Code & WSL Integration

- **Install WSL Extension in VS Code and connect to Ubuntu:**
  - Install [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl).
  - In the bottom left of VS code you can now connect to WSL, which loads Ubuntu

- **Create Target Directory For Project:**
  - The bottom panel of VS Code is minimized by default. You can either pull up the bottom line of VS code with the mouse, or toggle the panel (`Ctrl+J`)
  - By defaulte the terminal will show a bash terminal in your home directory on the Ubuntu side, which is `/home/maba`. 
  - In your user home directory, make a projects folder where we pull the GitHub repo(s). Two commands will create the folder and hop into it:
    ```bash
    mkdir projects
    cd projects
    ```
- **Clone Phenalizer repo:**    
  - You can now clone the Phenalyzer repo from GitHu you can either clone it via the GUI and setting the destination to `/home/maba/projects` and then open it, or use the bash terminal straight away:
    ```bash
    git clone https://github.com/mabar1/Phenalyzer
    cd Phenalyzer
    ```
  - After cloning the repo, you will notice that the terminal is set to 
    ```bash
    maba@DESKTOP-TH5JAGE:~/projects/Phenalyzer$
    ```
  - Since the repo is not cloning any .tsv files, the input folder is not cloned and we create it now manually. All .tsv files go in this folder:
    ```bash
    mkdir QuPath_tsv_inputs
    ```

### 3.3 Install R 4.5 in Ubuntu

- **Make Installer Executable and Install R:**
  - We strongly advice you to use WSL Ubuntu to run R, rather than installing the packages on the Windows side. This however means that we need to install R now also on the virtual machine. This is straight forward using the ```r_installer.sh``` script. After cloning, this script is read-only, and we need to change the permissions of that file in the bash terminal. Make sure you are still in the cloned directory and did not hop accidentally into the input folder. Then change its permission to read, write and execute by typing:
  ```bash
  chmod 777 r_installer.sh
  ```
  - afterwards, installing R is a piece of cake:
  ```bash
  ./r_installer.R
  ```
  - Enter your Ubuntu password as needed. Afterwards you can go and make a tea while occasionally pressing Enter of type yes as required.
  - After installation, restart VS code and verify that the R installation was successfull by calling for it in the bash terminal:
  ```bash
  which R
  ```
  - To which the terminal should return: `/usr/bin/R`.

### 3.4 Prepare VS Code R Environment

- **Install REditorSupport Extension in VS Code.**
  - Install [R Extension for Visual Studio Code](https://marketplace.visualstudio.com/items?itemName=reditorsupport.r) from the Marketplace or by searching for `REditorSupport` in the Extensions tab on the left panel of VS Code.
  - Restart VS Code
- **Install Language Server in R**
  - R needs to be initialized in VS Code. So far we had only one profile launched in the bottom panel, the bash terminal. Before installing any R packages, the R terminal is initialized by either clikcing the R R terminal button at the very bottom right of VS Code, or by clicking the `+` button on the right side, above the cube-icon of bash. Both actions should start R 4.5 in a second terminal. In the R terminal, we install R packages like so:
   ```r
   install.packages("languageserver")
   ```
  - R will attempt to install the packages into the unwritable `usr/local/lib/R/site-library` because we did not specify a library path. It therefore will promt you to use a personal library instead (probably inside your user directory here: `/home/maba/R/x86_64-pc-linux-gnu-library/4.5`). Confirm that.

- **Restart VS Code.**


### 3.5 Install add-ons to improve R Console experience: radian
 While the following three extensions are optional, installing radian is highly recommended. Radian is a python-based frontend that turns your R code more interactive by highlighting syntax or auto-completing your commands. 
- **Install radian (Improved R Console):**
  - radian is installed from the bash console. Since it is a python application wrapping around the R library, pipx will take care of the installation in a minute. We need pipx first, followed by radian:
  ```bash
  sudo apt install pipx
  pipx install radian
  ```
  - **Restart VS Code.**
  - after restarting VS Code, the bash terminal should be able to locate the radian installation. Type in the bash terminal:
    ```bash
    which radian
    ```
  - This should return `/home/maba/.local/bin/radian`.

- **Configure VS Code to Use radian:**
  - Open settings (`Ctrl+Shift+P` → "Preferences: Open Settings (UI)" or `Ctrl+,`).
  - Search for `Rterm`, specifically the field **R > Rterm: Linux** and paste there your path: `/home/maba/.local/bin/radian`.
  - Search for `radian` and enable **R: Bracketed Paste**.
  - Close and restart VS Code.
  - If you initialize now a new R terminal, you will see that the command line changed and starts with a blue `r$>`, indicating radian is active and runs the R terminal frontend.

### 3.6 Install add-ons to improve R Console experience: Plot Viewer
  - httpgd is by now fully integrated into VS Code. As the name suggests, this is an extension that watches your session and sends all plots through a WebSocket to a virtual server that allows you to interact with your plots and keeps all plots stored in one window.
- Open the VS Code Settings and search for `r.plot.useHttpgd` and enable it. All plots will now pop up in the same panel.

### 3.7 Install add-ons to improve R Console experience: Linting (lintr)
  - `lintr` is a R package, so you need to initialize an R console. We install the newest version from GitHub, for that we need a second package `remotes` that is pulling  `lintr`. Both packages will go into the user directory as before: `/home/maba/R/x86_64-pc-linux-gnu-library/4.5`
  ```r
  install.packages("remotes")
  remotes::install_github("r-lib/lintr")
  ```

- Lintr is a good way to write code correctly. However I found lintr to be too strict, so I toned linter down. To customize lintr, we create a starup file that disables a couple of linter settings. This file is created in your homes directory. We use the bash terminal for this:
  - In bash:
    ```bash
    cd ~
    vi ./lintr
    ```
  - This is an empty text file. Paste the following lines. Watch out, linter is *very fussy* about this file: Ensure the closing bracket is on the last line of the last command and there is a blank line at the end. Make sure you hit enter to create that last line:
    ```r
    linters: with_defaults(
      line_length_linter = NULL,
      open_curly_linter = NULL,
      commented_code_linter = NULL,
      trailing_whitespace_linter = NULL)

    ```
  - To store the file, press `Esc`, type `:wq` to save and exit.
  - Restart VS Code. If the terminal has notable delay and does not respond to your commands, the chances are high that the `./lintr` file is wrongly formated.

## 4. Install Phenalyzer Dependencies
- cmake is required to install a multitude of R packages (such as `s2, stars, nloptr, FlowSOM`). We will install cmake from the bash terminal, not R! In bash, type:
- **Install cmake for package dependencies:**
  ```bash
  sudo apt install cmake
  ```
- **Restart VS Code.**
- **Source the `environment_setup.R` script in your R terminal:**
  ```r
  source("environment_setup.R")
  ```
  - This will install all required packages for Phenalyzer. (This may take a while—make a coffee!)

---

**You are now ready to use Phenalyzer!**