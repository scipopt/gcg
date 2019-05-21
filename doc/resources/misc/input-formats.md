# Supported file formats {#input-formats}

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
