# Compatibility {#compatibility}
# Supported Operating Systems {#os-support}

|   | 32 bit | 64 bit | Supported Extensions | Comments |
|---|:---:|:---:|---|---|
| Ubuntu 18.04        |   | yes | GMP, hmetis, bliss |  |
| Debian 9            |   | yes | GMP, hmetis, bliss |  |
| macOS Catalina      | - | yes | GMP, hmetis        |  |
| Windows 10 ([WSL](https://docs.microsoft.com/de-de/windows/wsl/install-win10), Ubuntu 18.04) |  | yes | GMP, hmetis, bliss | Apply [required fix](https://www.turek.dev/post/fix-wsl-file-permissions/) for cmake|
| Windows 10 ([Visual Studio](https://visualstudio.microsoft.com/de/)) | no | no |  | tested with Visual Studio 2017 |
| Windows 10 ([minGW](http://www.mingw.org/)) |  |  |  |  |

For Windows, please also refer to our @subpage windows-install .
# Supported Compilers {#compiler-support}
| Compiler | Support |
|:---:     |:---:    |
| gcc      | yes     |


# Supported File Formats {#input-formats}

## Problem file formats
GCG supports all file formats supported by SCIP to read in problems, solution, etc.
The original problem can for example be read in as:

    .lp
    .mps
    .cip

files.

## Decomposition file formats
If GCG is not able to automatically detect a structure suitable to perform a
Dantzig-Wolfe reformulation, you need to specify the structure yourself and make
it available to GCG. There are some file formats for structure information.
Currently, GCG supports

    .dec
    .blk

files. They are documented in reader_dec.h and reader_blk.h, respectively.
You can find examples for `dec` as well as `blk` files in the `check/instances/` directory.

### Automatic decomp-file adding
During detection, GCG will automatically detect similiarly named `.dec`-files
if they are located in the same directory. To prevent that, set the corresponding
mode (for `make test`, see \ref makefiles-args).
