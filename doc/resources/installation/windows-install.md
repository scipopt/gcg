# Windows Installation Guide {#windows-install}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Currently, GCG cannot be installed on pure Windows. This is due to unportable commands in the GCG source code. Still, there are **ways to get GCG running on your Windows computer**.
Generally, there are three different options:

1. Install a full Linux on your PC
2. Install a Virtual Machine on your PC
3. Use the Windows Subsystem for Linux (WSL) (see our @ref wsl "guide"))

While the first two options are more comfortable in the long-run, using the **WSL is a quicker variant** and more suitable **for users who are less experienced with Linux**.

## Windows Subsystem for Linux (WSL) {#wsl}
For the mentioned first two cases, there are guides on the internet after which you can simply use the usual @ref install "Installation Guides". 
For the WSL, however, we want to provide a more **extensive guide**.
### 1. Checking Requirements
Using the WSL requires you to have...
- Windows 10 with 2019-Nov Update installed
- Windows Store installed (currently, login can still be skipped)
- [optional] An IDE of your choice (e.g. [Visual Studio Code](https://code.visualstudio.com/))

### 2. Installing the Windows Subsystem
1. **Install** the Windows Subsystem for Linux using the [guide provided by Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install-win10). We tested Ubuntu 20.04 LTS, but it should run on other distributions as well.
2. After the final step (the download from Windows Store), **launch the WSL**, e.g. by typing "Ubuntu" after opening the start menu. It will take a while to install (do not abort the process or close the window).
3. Following the installation, you will be asked you for a **username and a password**. Remember that username and password. 
4. The `home` directory of the Linux subsystem can be reached via the path `\\wsl$\Ubuntu-xx.xx` (can be entered after opening the start menu or in a windows explorer window; replace `xx.xx` by your version). 
We recommend to **create a folder there** (e.g. `myfolder/`) and giving your user permissions using `sudo chown username myfolder/`. This way, you won't need `sudo` for every command.

### 3. Installing GCG and the SCIP Optimization Suite
* Theoretically, you can just follow the usual @ref install "installation guide", as if you were on Linux. 
* Under Windows, we currently recommend using Makefiles (CMake has some compatibility problems in the WSL).
* If you are new to Linux, we recommend to use one of our automated installers: (link located on the installation page), 
  1. Download one of the automated installer scripts
  2. Move it to your `repos` (or similar) folder using the windows explorer, navigating to `\\wsl$\Ubuntu-xx.xx\home\<username>\myfolder\` (replace `xx.xx` by your version and `<username>` by your WSL username)
  3. Make it executable using `sudo chmod +x <installername.sh>` (replace `<installername.sh>` by your installer, which should end with `.sh`).
  4. Execute it using `./<installername.sh>` and follow the installation instructions.

### 4. Using an IDE to develop features for GCG
When using an IDE, you have to remember to **switch your terminal** to the Linux Subsystem instead of the Windows Powershell or possible different installs like the Git Bash for Windows.
Other than that, you can simply **open the folder containing your GCG** within the `\\wsl$\Ubuntu...` paths.

#### Caveats when using WSL
1. The policy management in WSL is not exactly as in usual Linux. There are different things to pay attention to.
    - Always remember to own the folder your repositories reside in. If your folder is called `repos`, execute `chown <username> repos/` from the folder above.
    - When generating the SSH key for git, try to regenerate your SSH key using `sudo ssh-keygen` (with your desired arguments). Note that then the key will be generated as root and can only be opened as root (`sudo`).
2. The package management and DNS resolving is not perfect.
    - Always update your package management system (`sudo apt-get update`) before installing anything.
    - Do not connect to any VPN when downloading packages.
3. Folders containing blank spaces might generate problems. Use underscores instead.
