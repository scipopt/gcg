# Windows Installation Guide {#windows-install}
> This page guides you through **installing GCG on a Windows computer**.

Generally, there are three different options:

1. Use the Microsoft Visual C++ (MSVC) compiler together with CMake (see our @ref msvc "guide").
2. Use the Windows Subsystem for Linux (WSL) (see our @ref wsl "guide").
3. Install a Virtual Machine on your PC.

While the first option will build native Windows executables the other two options rely on Linux. Please note that many scipts shipped with GCG typically require Linux. Hence, if you want to use them, you should choose option 2 or 3.

# Microsoft Visual C++ (MSVC) {#msvc}
In this section, we descibe the steps needed to build GCG using the Microsoft Visual C++ (MSVC) compiler. We do this using the Git installation as an example. The guide assumes that the following tools are installed, configured properly, and accesible through the Micorsoft Windows PowerShell command-line shell:

* Git
* CMake
* MSCV compiler (e.g., Microsoft Visual Studio 2022)

## 1. Install dependencies
We recommend using [vcpkg](https://github.com/microsoft/vcpkg) to install required dependencies.

### Install vcpkg (skip this if you have vcpkg already available on your system)
1. Open a PowerShell window and navigate to a directory that will be used a parent directory for the vcpkg installation.
2. Clone the repository: `git clone git@git.or.rwth-aachen.de:gcg/gcg.git`
3. `cd vcpkg`
4. `./bootstrap-vcpkg.bat -disableMetrics`
5. `./vcpkg integrate install`
6. Set the environment variable `VCPKG_ROOT` to the path of your vcpkg installation (e.g., `C:\vcpkg`). You can set the variable globally (recommended) or for each PowerShell session.

If you wish, you can add the vcpkg directory to your `PATH` environment variable.

### Install the dependencies using vcpkg
1. Open a PowerShell window and navigate into the vcpkg directory (not required if you added it to your `PATH` environment variable)
2. `./vcpkg install readline gmp zlib gsl jansson` (leave out the `./` if you modified your `PATH` environment variable)

## 2. Build GCG
After installing the dependencies we are ready to clone and build GCG.

### Clone GCG
1. Navigate to a proper parent directory.
2. `git clone https://github.com/scipopt/gcg.git`
3. `cd gcg`

### Build GCG
If your installed CMake version is equal to 3.25 or higher, you can configure and build GCG as well as test the build by calling
* `cmake --workflow --preset gcg-windows-release`
By using `--preset gcg-windows-debug` a debug build will be compiled and tested.

Otherwise, you can use the following commands
1. Configure: `cmake -S . -B build -DGCG_DEV_BUILD=ON -DZIMPL=OFF -DIPOPT=OFF -DPAPILO=OFF -DCMAKE_TOOLCHAIN_FILE="$env:VCPKG_ROOT"/scripts/buildsystems/vcpkg.cmake`
2. Build: `cmake --build build --target gcg --config Release`
3. Test: `cmake --build build --target gcg_check --config Release`
You can use `--config Debug` instead of `--config Release` to compile and test debug builds.

# Windows Subsystem for Linux (WSL) {#wsl}
@htmlonly
<iframe width="560" height="315" src="https://www.youtube-nocookie.com/embed/ruItZK9NU6c" style="margin:auto; display:block" frameborder="3" allow="accelerometer; autoplay; clipboard-write; encrypted-media; gyroscope; picture-in-picture" allowfullscreen></iframe>
@endhtmlonly

## 1. Checking Requirements
Using the WSL requires you to have...
- Windows 10 with 2019-Nov Update installed
- Windows Store installed (currently, login can still be skipped)
- [optional] An IDE of your choice (e.g. [Visual Studio Code](https://code.visualstudio.com/))

## 2. Installing the Windows Subsystem
1. **Install** the Windows Subsystem for Linux using the [guide provided by Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install-win10). We tested Ubuntu 20.04 LTS, but it should run on other distributions as well.
2. After the final step (the download from Windows Store), **launch the WSL**, e.g. by typing "Ubuntu" after opening the start menu. It will take a while to install (do not abort the process or close the window).
3. Following the installation, you will be asked you for a **username and a password**. Remember that username and password. 
4. The `home` directory of the Linux subsystem can be reached via the path `\\wsl$\Ubuntu-xx.xx` (can be entered after opening the start menu or in a windows explorer window; replace `xx.xx` by your version). 
We recommend to **create a folder there** (e.g. `myfolder/`) and giving your user permissions using `sudo chown username myfolder/`. This way, you won't need `sudo` for every command.

## 3. Installing GCG and the SCIP Optimization Suite
* Theoretically, you can just follow the usual @ref install "installation guide", as if you were on Linux. 
* If you are new to Linux, we recommend to use one of our automated installers: (link located on the installation page), 
  1. Download one of the automated installer scripts
  2. Move it to your `repos` (or similar) folder using the windows explorer, navigating to `\\wsl$\Ubuntu-xx.xx\home\<username>\myfolder\` (replace `xx.xx` by your version and `<username>` by your WSL username)
  3. Make it executable using `sudo chmod +x <installername.sh>` (replace `<installername.sh>` by your installer, which should end with `.sh`).
  4. Execute it using `./<installername.sh>` and follow the installation instructions.

## 4. Using an IDE to develop features for GCG
When using an IDE, you have to remember to **switch your terminal** to the Linux Subsystem (called "WSL Bash") instead of the Windows Powershell or possible different installs like the Git Bash for Windows.
Other than that, you can simply **open the folder containing your GCG** within the `\\wsl$\Ubuntu...` paths (they are classed as network devices!) using the "Remote" feature of your IDE (lower left corner
in VSCode, see this [guide](https://code.visualstudio.com/docs/remote/wsl)).

### Caveats when using WSL
1. The policy management in WSL is not exactly as in usual Linux. There are different things to pay attention to.
    - Always remember to own the folder your repositories reside in. If your folder is called `repos`, execute `chown <username> repos/` from the folder above.
    - When generating the SSH key for git, try to regenerate your SSH key using `sudo ssh-keygen` (with your desired arguments). Note that then the key will be generated as root and can only be opened as root (`sudo`).
2. The package management and DNS resolving is not perfect.
    - Always update your package management system (`sudo apt-get update`) before installing anything.
    - Do not connect to any VPN when downloading packages.
3. Folders containing blank spaces might generate problems. Use underscores instead.
