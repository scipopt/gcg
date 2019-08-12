# Compatibility {#compatibility}
# Supported Operating Systems {#os-support}

|   | 32 bit | 64 bit | Supported Extensions | Comments |
|---|:---:|:---:|---|---|
| Windows 10 (Visual Studio) | no | no |  | tested with Visual Studio 2017 |
| Windows 10 (minGW) | ? | ? |  |  |
| Debian 9 | ? | yes | GMP, hmetis, bliss |  |
| Ubuntu 18.04 | ? | yes | GMP, hmetis, bliss |  |
| macOS | ? | ? |   |  |

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


# Known Bugs {#known-bugs}

<table>
  <tr>
    <th>Bug</th>
    <th>Reason</th>
    <th>Version</th>
  </tr>
  <tr>
    <td>Pressing CTRL-C does not always terminate GCG</td>
    <td>GCG creates SCIP instances for subproblems. If CTRL-C is pressed, one of these subinstances is stopped, but not the calling function of GCG.<br></td>
    <td>3.0.0, current master</td>
  </tr>
  <tr>
    <td>As described in the SCIP documentation, one should be able to run a test with
    `SETTINGS="set1,set2"` and SCIP should run the test on both settings. This does not work.</td>
    <td>Unknown</td>
    <td>3.0.0, current master</td>
  </tr>
  <tr>
    <td>If you install different versions of the SCIP Optimization Suite (including GCG), there might occur compatibility issues.</td>
    <td>GCG might link to the most recent version instead of the one that the link points to.<br></td>
    <td>all versions</td>
  </tr>
</table>
