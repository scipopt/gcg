#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: run_mosek.sh,v 1.5 2011/03/21 13:52:23 bzfgamra Exp $

# Note, this script uses the PyMosek API.
# Therefore, the system environment must be set up correctly:
#  * The environment variable LD_LIBRARY_PATH must contain the path to the mosek library.
#  * The environment variable PYTHONPATH must point to the location of PyMosek module.
# See the MOSEK manual for further help. 

# The PyMosek module is compatible with Python 2.5 and all later 2.X versions.
# See the MOSEK manual for more details.
PYTHON=python2.6

SOLVER=$1
LPS=$2
NAME=$3
TIMELIMIT=$4
SOLFILE=$5
MIPGAP=$6

TMPFILE=check.$SOLVER.tmp

echo > $TMPFILE
echo > $SOLFILE

# read, optimize, display statistics, write solution, and exit
# use deterministic mode (warning if not possible) 

# read, optimize, display statistics, write solution, and exit
echo "import sys"                                                                          >> $TMPFILE
echo "import mosek"                                                                        >> $TMPFILE
echo "def mosekprint(txt):"                                                                >> $TMPFILE
echo "    sys.stdout.write(txt)"                                                           >> $TMPFILE 
echo "env  = mosek.Env()"                                                                  >> $TMPFILE
echo "task = env.Task(0,0)"                                                                >> $TMPFILE
echo "task.set_Stream(mosek.streamtype.log,mosekprint)"                                    >> $TMPFILE
echo "task.readdata(\"$NAME\")"                                                            >> $TMPFILE  
# set mipgap to given value
echo "task.putdouparam(mosek.dparam.mio_tol_rel_gap,$MIPGAP)"                              >> $TMPFILE            
echo "task.putdouparam(mosek.dparam.optimizer_max_time,$TIMELIMIT)"                        >> $TMPFILE            
echo "trmcode = task.optimize()"                                                           >> $TMPFILE  
echo "task.solutionsummary(mosek.streamtype.msg)"                                          >> $TMPFILE
echo "if trmcode not in [mosek.rescode.ok,mosek.rescode.trm_max_time]:"                    >> $TMPFILE 
echo "    aborted = 'TRUE'"                                                                >> $TMPFILE  
echo "else:"                                                                               >> $TMPFILE  
echo "    aborted = 'FALSE'"                                                               >> $TMPFILE  
echo "if trmcode in [mosek.rescode.trm_max_time]:"                                         >> $TMPFILE     
echo "    timeout = 'TRUE'"                                                                >> $TMPFILE  
echo "else:"                                                                               >> $TMPFILE  
echo "    timeout = 'FALSE'"                                                               >> $TMPFILE  
echo "if task.getintinf(mosek.iinfitem.mio_num_relax)>0:"                                  >> $TMPFILE
echo "    dobj = '%24.16e' % task.getdouinf(mosek.dinfitem.mio_obj_bound)"                 >> $TMPFILE 
echo "else:"                                                                               >> $TMPFILE
echo "    dobj = '-infty'"                                                                 >> $TMPFILE
echo "if task.getintinf(mosek.iinfitem.mio_num_int_solutions)>0:"                          >> $TMPFILE 
echo "    pobj = '%-24.16e' % task.getdouinf(mosek.dinfitem.mio_obj_int)"                  >> $TMPFILE
echo "else:"                                                                               >> $TMPFILE
echo "    pobj = '+infty'"                                                                 >> $TMPFILE
echo "s = mosek.soltype.itg"                                                               >> $TMPFILE 
echo "[prosta,solsta] = task.getsolutionstatus(s)"                                         >> $TMPFILE
echo "if prosta in [ mosek.prosta.prim_infeas ]:"                                          >> $TMPFILE 
echo "    pobj = 'infeasible'"                                                             >> $TMPFILE 
echo "    dobj = '-'"                                                                      >> $TMPFILE 
echo "f = open(\"$SOLFILE\",\"wt\")"                                                       >> $TMPFILE
echo "t = mosek.solsta"                                                                    >> $TMPFILE             
echo "if solsta in [t.integer_optimal, t.near_integer_optimal, t.optimal, t.near_optimal, t.prim_feas, t.near_prim_feas ]:" >> $TMPFILE 
echo "     numvar = task.getnumvar()"                                                      >> $TMPFILE
echo "     f.write('=obj=%24.16e\n' % task.getprimalobj(s))"                               >> $TMPFILE 
echo "     for j in range(0,numvar):"                                                      >> $TMPFILE 
echo "          sk,x,sl,su,sn = task.getsolutioni(mosek.accmode.var,j,mosek.soltype.itg)"  >> $TMPFILE
echo "          f.write('%s %24.16e\n' % (task.getvarname(j),x))"                          >> $TMPFILE 
echo "else:"                                                                               >> $TMPFILE
echo "    f.write(\"=infeas=\")"                                                           >> $TMPFILE 
echo "f.close()"                                                                           >> $TMPFILE
echo "print 'MIPLIBsolverversion = %d.%d.%d.%d' % mosek.getversion()"                      >> $TMPFILE
echo "print 'MIPLIBbbnodes = %d' % task.getintinf(mosek.iinfitem.mio_num_relax)"           >> $TMPFILE
echo "print 'MIPLIBdb = %s' % dobj"                                                        >> $TMPFILE
echo "print 'MIPLIBpb = %s' % pobj"                                                        >> $TMPFILE
echo "print 'MIPLIBaborted = %s' % aborted"                                                >> $TMPFILE 
echo "print 'MIPLIBtimeout = %s' %  timeout"                                               >> $TMPFILE 

$PYTHON $TMPFILE
rm $TMPFILE





