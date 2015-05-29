#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*            This file is part of the test engine for MIPLIB2010            *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: parse_gurobi.awk,v 1.3 2011/04/12 12:46:02 bzfgamra Exp $


# set all solver specific data:
#  solver ["?"]
#  solverversion ["?"]
#  solverremark [""]
#  bbnodes [0]
#  db [-infty]
#  pb [+infty]
#  aborted [1]
#  timeout [0]

BEGIN {
  solver = "Gurobi";
}
# solver version
/^Gurobi Optimizer version/ { solverversion = $4; }
# branch and bound nodes
/^Explored/ { bbnodes = $2; }
# infeasible model
/^Model is infeasible/ {
  db = pb;
}
# dual and primal bound
/^Best objective/ {
 if( $3 != "-," )
  pb = $3 + 0.0;
 if( $6 != "-," )
  db = $6 + 0.0;
}
# solving status
/^Explored/ { aborted = 0; }
/^Time limit reached/ { timeout = 1; }
/^Solve interrupted/ { timeout = 1; }

