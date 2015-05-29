#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

SOLVER=$1
BINNAME=$2
NAME=$3
TIMELIMIT=$4
SOLFILE=$5
THREADS=$6
MIPGAP=$7

TMPFILE=check.$SOLVER.tmp

echo > $TMPFILE
echo > $SOLFILE

# set threads to given value 
if test $THREADS != 0
then
    echo threads=$THREADS                >> $TMPFILE
fi

# set mipgap to given value
echo miprelstop=$MIPGAP                  >> $TMPFILE

# set timing to wall-clock time and pass time limit
echo cputime=0                           >> $TMPFILE
echo maxtime=-$TIMELIMIT                 >> $TMPFILE

# use deterministic mode (warning if not possible) 
# nothing for Xpress

# set display options
# nothing for Xpress

# read, optimize, display statistics, write solution, and exit
echo readprob $NAME                      >> $TMPFILE
echo                                     >> $TMPFILE
echo mipoptimize                         >> $TMPFILE

# Output some additional information
echo format \"Objective value = %s\" [ if {\$mipstatus==1 \|\| \$mipstatus==2} {\$lpobjval} else {\$mipobjval} ] >> $TMPFILE
echo format \"Best Bound = %s\" \$bestbound                         >> $TMPFILE
echo puts \"Nodes explored = \$nodes\"                              >> $TMPFILE

echo if {\$mipsols \> 0} {writeslxsol -p -m solution.slx}            >> $TMPFILE
# Append objective value to solution file
echo set SLXFH [open solution.slx a]                                >> $TMPFILE
echo puts \$SLXFH \" X =obj= \$mipobjval\"                          >> $TMPFILE
echo if {\(\$mipstatus==1 \|\| \$mipstatus==2\) \&\& \$lpstatus==2} { puts \$SLXFH \"=infeas=\" }  >> $TMPFILE
echo if {\$mipstatus==5} {puts \$SLXFH \"=infeas=\"}                >> $TMPFILE
echo quit                                                           >> $TMPFILE

$BINNAME < $TMPFILE 

if test -f solution.slx
then
    # translate XPRESS solution format into format for solution checker.
    #  The SOLFILE format is a very simple format where in each line 
    #  we have a <variable, value> pair, separated by spaces. 
    #  A variable name of =obj= is used to store the objective value 
    #  of the solution, as computed by the solver. A variable name of 
    #  =infeas= can be used to indicate that an instance is infeasible.
    sed -n "s/\(=infeas=\)\|\( *[A-Za-z]  *\(.*[^ ]  *[0-9.Ee+-][0-9.Ee+-]*\)\)/\1\3/gp" < solution.slx > $SOLFILE
    rm solution.slx
fi

# remove tmp file
rm $TMPFILE

