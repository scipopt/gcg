# do automatic testing # {#testing}
Just like SCIP, GCG supports automatic testing. It uses the same command
structure as SCIP, so if any information is missing,
you can also consult the SCIP documentation.


# Testing Basics #

### How do I test?
#### Simple Testing
If you want to check whether your own modules generate the right results,
you can work with the testset files `*.test` in the folder `/check/testset`
to specify which instances you want to test. The `*.solu` files specify their
respective solution. <i>For further syntax of the `*.solu`-files, please check
the <a href="https://scip.zib.de/doc-6.0.1/html/TEST.php">SCIP documentation</a>.</i><br>
A simple test on `testrun.test` is executed with the following command
in the GCG main folder:

    make TEST=testrun test

Note that the default test is `short.test`, which will be executed if no
`TEST`-argument is given.

#### Testing with statistics
In order to use the test data to create visualizations, the testing has to be
done with statistics enabled. You can do that with:

    make TEST=testrun STATISTICS=true test


#### Testing with defined settings
To test with your own settings, e.g. `my.set` (located in `/settings`)
you can use the following command:

    make TEST=testrun SETTINGS=my test


### What files do I get?
In the folder `/check/results` you will see multiple files per test.
The output files for a test executed by
`make TEST=testrun SETTINGS=my STATISTICS=true test` are as follows:

    check.testrun.gcg-GCGVERSION.linux.x86_64.gnu.opt.spx2.COMPUTERNAME.my.default.{.err,.out,.res}

and the files `{.pav,.tex}` are added after the test was completed.
If you aborted the test before it finished, you can still generate those files
(incompletely) using the command `./evalcheck.sh` in the `/check` folder.

### What do the files contain?
The .out files contain the `stdout`-output and therefore all data that was
printed during the execution of GCG.<br>
The .res files only contain the execution table with
primal bound, dual bound, gap and so on.<br>
The latter file is also given as (incomplete) LaTeX code in the .tex file.<br>

# Using the results #
The last chapter only covered the absolute basics of testing, mainly used
for verification purposes. But heads up, the output files still have a purpose,
even if you verified that your GCG works (and how fast it works).

### Plot performance profiles
Using the tool `plotperprof.py`, you can plot performance profiles like this one:
<img src="">
The usage of this tool is as follows:
