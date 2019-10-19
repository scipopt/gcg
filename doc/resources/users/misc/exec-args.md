# Command Line Arguments {#exec-args}

GCG can be started from the command line without the @ref interactive-menu "interactive menu".
For example, this is used by @ref u1 "Use Case 1" of the @ref users "Users Guide". The
following arguments can be used when executing GCG:

    syntax: ./bin/gcg [-l <logfile>] [-q] [-s <settings>] [-f <problem>] [-m <mastersettings>] [-d <decomposition>] [-b <batchfile>] [-c "command"]
      -l <logfile>           : copy output into log file
      -q                     : suppress screen messages
      -s <settings>          : load parameter settings (.set) file
      -m <mastersettings>    : load master parameter settings (.set) file
      -f <problem>           : load and solve problem file
      -d <decomposition>     : load decomposition file
      -o <primref> <dualref> : pass primal and dual objective reference values for validation at the end of the solve
      -b <batchfile>         : load and execute dialog command batch file (can be used multiple times)
      -c "command"           : execute single line of dialog commands (can be used multiple times)
