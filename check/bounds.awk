#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program and library             *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#*  Licensed under the Apache License, Version 2.0 (the "License");          *
#*  you may not use this file except in compliance with the License.         *
#*  You may obtain a copy of the License at                                  *
#*                                                                           *
#*      http://www.apache.org/licenses/LICENSE-2.0                           *
#*                                                                           *
#*  Unless required by applicable law or agreed to in writing, software      *
#*  distributed under the License is distributed on an "AS IS" BASIS,        *
#*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. *
#*  See the License for the specific language governing permissions and      *
#*  limitations under the License.                                           *
#*                                                                           *
#*  You should have received a copy of the Apache-2.0 license                *
#*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#
#@file    bounds.awk
#@brief   evaluate development of dual/primal bounds
#@author  Christian Puchert
#
function abs(x)
{
   return x < 0 ? -x : x;
}
function min(x,y)
{
   return (x) < (y) ? (x) : (y);
}
function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}

BEGIN {
   infty = +1e+20;

   nprobs = 0;

   printf("# Development of bounds and gap during the solving process\n")  >DATFILE;
}

#
# problem name
#
/^@01/ {
   filename = $2;

   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > 12 )
      shortprob = substr(prob, length(prob)-11, 12);
   else
      shortprob = prob;

   # problem type
   bp = match($2, "problems/binpacking");
   cut = match($2, "problems/cuttingstock");
   col = match($2, "problems/coloring");
   cpmp = match($2, "problems/cpmp");
   ga = match($2, "problems/gap");

   if( bp != 0 )
      type = 1;
   else if ( cut != 0 )
      type = 2;
   else if ( col != 0 )
      type = 3;
   else if ( cpmp != 0 )
      type = 4;
   else if ( ga != 0 )
      type = 5;
   else
      type = 0;

   printf("# Instance: %s\n", prob) >DATFILE;

   pb = "      --      ";
   db = "      --      ";
   gap = "    Inf";
}
#
# primal and dual bounds during the solving process
#
/^ *(\*)?[A-Za-z\*] *[0-9]*(\.[0-9]*)?s *\|/ {
   n = split($0, a, "|");
   m = match(a[1],/[0-9]*/);

   timestamp = substr(a[1],3,4);
   currentdb = a[n-2];
   currentpb = a[n-1];
   currentgap = a[n];

   printf("%s %s %s %s %d\n", timestamp, currentdb, currentpb, currentgap, type)  >DATFILE;
}
#
# problem: master or original?
#
/^Original Program statistics:/ { inmasterprob = 0; inoriginalprob = 1; }
/^Master Program statistics:/ { inmasterprob = 1; inoriginalprob = 0; }
#
# solution
#
/^  Primal Bound     :/ {
   if( inoriginalprob && $4 != "infeasible" && $4 != "-" ) {
      pb = sprintf("%13.6e", $4);
   }
}
/^  Dual Bound       :/ {
   if( inoriginalprob && $4 != "-" )
      db = sprintf("%13.6e", $4);
}
/^  Gap              :/ {
   if( inoriginalprob && $4 != "infinite" )
      gap = sprintf("%7.2f%%", $4);
}
#
# time
#
/^Solving Time       :/ { tottime = sprintf("%4.1f", $4) } # for older scip version ( < 2.0.1.3 )
/^  solving          :/ { tottime = sprintf("%4.1f", $3) }

/^=ready=/ {

   # print bounds and gap at the end of the solving process;
   # append two new lines to indicate a new data set
   printf("%s %s %s %s %d\n", tottime, db, pb, gap, type)  >DATFILE;
   printf("\n\n")                                          >DATFILE;
   
   nprobs++;
}

END {

   printf("# Development of bounds and gap during the solving process\n")  >GPFILE;
   printf("# Problem types:\n")                                            >GPFILE;
   printf("# 0: unknown\n")                                                >GPFILE;
   printf("# 1: Bin Packing\n")                                            >GPFILE;
   printf("# 2: Cutting Stock\n")                                          >GPFILE;
   printf("# 3: Coloring\n")                                               >GPFILE;
   printf("# 4: Capacitated p-median\n")                                   >GPFILE;
   printf("# 5: Generalized Assignment\n")                                 >GPFILE;
   printf("\n")

   printf("set terminal pdf\n")             >GPFILE;
   printf("set xlabel \"Time\"\n")          >GPFILE;
   printf("set ylabel \"Value\"\n")         >GPFILE;
   printf("set grid xtics ytics\n")         >GPFILE;
   printf("set datafile missing \"--\"\n")  >GPFILE;
   for( i = 1; i <= nprobs; ++i ) {
      printf("plot for [i=2:3] '%s' index %d using 1:i with lp\n", DATFILE, i-1)  >GPFILE;
   }

}
