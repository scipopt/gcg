# Makefiles Arguments {#makefiles-args}
You can modify your Makefiles installation by setting flags during compilation. For example, for some visualization
scripts, you'll need `STATISTICS=true` to be set **during compilation**. Important flags are listed here
(a complete list can be found by typing `make help` (GCG-specific arguments) or `make --help`):

### Main Arguments
#### Compilation, debugging and statistics (only for `make`)

    -j [N], --jobs[=N]          Allow N jobs at once; infinite jobs with no arg.
    --debug[=FLAGS]             Print various types of debugging information.
    STATISTICS[=B]              Print additional statistics (esp. for pricing)

#### Additional Packages (only for `make`)

    READLINE[=B]                enables READLINE, required for command line interaction (default: true)
    ZIMPL[=B]                   enables ZIMPL, required to convert .zpl files to .lp/.mps files
    GAMS[=B]                    enables GAMS, required to convert .gms files to .lp/.mps files
    CLIQUER[=B]                 enables the Cliquer
    HMETIS[=B]                  enables hMETIS (Hypergraph & Circuit Partitioning)
    OPENMP[=B]                  enables parallelization using OpenMP
    GSL[=B]                     enables the GNU Scientific Library, needed for one detector

#### Limits and Modes (only for `make test`)

    TEST[=NAME]                 Name of the test set. Default: "short".
    SETTING[=SET]               Choose settings for the test run as defined in settings/SET.set.
    STATISTICS[=B]              Print additional statistics (beware: different to the one above,
                                which is used during compilation), B in {true,false}.
    NODE[=N]                    Limit of Nodes to be opened during branching.
    TIME[=N]                    Time limit for the whole solving
    OPT[=TYPE]                  Choose from {opt, dbg, prf}.
    MODE[=TYPE]                 Use different modes. 0 or none to prevent from using dec files.

### GCG-unspecific arguments
#### Compilation process (errors)

    -i, --ignore-errors         Ignore errors from recipes.
    --warn-undefined-variables  Warn when an undefined variable is referenced.
    -k, --keep-going            Keep going when some targets can't be made.
    -S, --no-keep-going, --stop Turns off -k.

#### Compilation process (speedup)

    -j [N], --jobs[=N]          Allow N jobs at once; infinite jobs with no arg.
    -l [N], --load-average[=N], --max-load[=N]
                                Don't start multiple jobs unless load is below N.
    -O[TYPE], --output-sync[=TYPE]
                                Synchronize output of parallel jobs by TYPE.

#### Compilation process (messages)

    -d                          Print lots of debugging information.
    -s, --silent, --quiet       Don't echo recipes.
    -v, --version               Print the version number of make and exit.
    -w, --print-directory       Print the current directory.
    --no-print-directory        Turn off -w, even if it was turned on implicitly.
    --trace                     Print tracing information.
