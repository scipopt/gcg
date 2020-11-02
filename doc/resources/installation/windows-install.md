# Windows Installation Guide {#windows-install}
> **This page is still in development and may be incomplete. Please excuse any inconveniences.**

Currently, GCG cannot be installed on pure Windows. This is due to unportable commands in the GCG source code. Still, there are ways to get GCG running on your Windows computer.\n
The following guide provides instructions to install GCG and also use an IDE with debugger such as VSCode and Jetbrains.

### Requirements
- Windows 10 with 2019-Nov Update installed
- Windows Store installed (currently, login can still be skipped)
- [optional] An IDE of your choice (e.g. [Visual Studio Code](https://code.visualstudio.com/))

### Windows Subsystem Installation
1. Install the Windows Subsystem for Linux using the [guide provided by Microsoft](https://docs.microsoft.com/en-us/windows/wsl/install-win10). We tested the [Ubuntu 18.04 LTS](https://www.microsoft.com/de-de/p/ubuntu-1804-lts/9n9tngvndl3q), but it should run on other distributions as well.
2. [cmake] As of early 2020, it is required to apply a fix to use `cmake` in the WSL. A guide on how to apply this fix can be found [here](https://www.turek.dev/post/fix-wsl-file-permissions/).
3. Finally, follow the usual @ref install "installation guide", **as if you were on Linux**.
