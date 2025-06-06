# Phenalyzer Software R Setup Guide

## 1. Windows Preparation

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
  - Follow prompts to download Ubuntu and create a user account (no spaces in username).

- **Restart your PC.**

## 2. VS Code & WSL Integration

- **Install WSL Extension in VS Code:**
  - Install [Remote - WSL](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-wsl).

- **Connect to WSL:**
  - Open VS Code, connect to WSL: Ubuntu.

- **Set Up Project Directory:**
  - In the bash terminal:
    ```bash
    mkdir projects
    cd projects
    ```
  - Clone the Phenalyzer repo from GitHub:
    ```bash
    git clone https://github.com/mabar1/Phenalyzer
    cd Phenalyzer
    ```

## 3. Install R 4.5 in Ubuntu

- **Make Installer Executable:**
  ```bash
  chmod 777 r_installer.sh
  ```
- **Run the Installer:**
  ```bash
  ./r_installer.R
  ```
  - Enter your Ubuntu password as needed.
  - Follow prompts (press Enter or type `y` as required).

- **Verify R Installation:**
  ```bash
  which R
  ```
  - Should return `/usr/bin/R`.

## 4. VS Code R Environment

- **Install REditorSupport Extension in VS Code.**
- **Restart VS Code.**

- **Initialize R Terminal in VS Code:**
  - Click the R terminal button (should start R 4.5).

- **Install Language Server in R:**
  ```r
  install.packages("languageserver")
  ```
  - If prompted about unwritable directories, choose to use your personal library (e.g., `/home/maba/R/x86_64-pc-linux-gnu-library/4.5`).

- **Restart VS Code.**

## 5. Install radian (Improved R Console)

- **Install pipx and radian:**
  ```bash
  sudo apt install pipx
  pipx install radian
  ```

- **Install cmake for package dependencies:**
  ```bash
  sudo apt install cmake
  ```

- **Restart VS Code.**

- **Verify radian Installation:**
  ```bash
  which radian
  ```
  - Should return `/home/maba/.local/bin/radian`.

- **Configure VS Code to Use radian:**
  - Open settings (`Ctrl+Shift+P` → "Preferences: Open Settings (UI)" or `Ctrl+,`).
  - Search for `Rterm` and set **R > Rterm: Linux** to `/home/maba/.local/bin/radian`.
  - Search for `radian` and enable **R: Bracketed Paste**.
  - Close and restart VS Code.

- **Check Radian in Terminal:**
  - Open a new R terminal. The prompt should start with a blue `r$>`, indicating radian is active.

## 6. Enable Plot Viewer

- **In VS Code Settings:**
  - Search for `r.plot.useHttpgd` and enable it.

## 7. Linting Setup (Optional)

- **Install remotes and lintr in R:**
  ```r
  install.packages("remotes")
  remotes::install_github("r-lib/lintr")
  ```

- **Customize lintr (Optional):**
  - In bash:
    ```bash
    cd ~
    vi ./lintr
    ```
  - Paste the following (ensure the closing bracket is on the last line and there is a blank line at the end):
    ```r
    linters: with_defaults(
      line_length_linter = NULL,
      open_curly_linter = NULL,
      commented_code_linter = NULL,
      trailing_whitespace_linter = NULL)

    ```
  - In `vi`, press `Esc`, type `:wq` to save and exit.

## 8. Install Phenalyzer Dependencies

- **Restart VS Code.**
- **Source the `environment_setup.R` script in your R terminal:**
  ```r
  source("environment_setup.R")
  ```
  - This will install all required packages for Phenalyzer. (This may take a while—make a coffee!)

---

**You are now ready to use Phenalyzer!**