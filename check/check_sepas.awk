#!/usr/bin/awk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* Copyright (C) 2010-2012 Operations Research, RWTH Aachen University       *
#*                         Zuse Institute Berlin (ZIB)                       *
#*                                                                           *
#* This program is free software; you can redistribute it and/or             *
#* modify it under the terms of the GNU Lesser General Public License        *
#* as published by the Free Software Foundation; either version 3            *
#* of the License, or (at your option) any later version.                    *
#*                                                                           *
#* This program is distributed in the hope that it will be useful,           *
#* but WITHOUT ANY WARRANTY; without even the implied warranty of            *
#* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
#* GNU Lesser General Public License for more details.                       *
#*                                                                           *
#* You should have received a copy of the GNU Lesser General Public License  *
#* along with this program; if not, write to the Free Software               *
#* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#
#@file    check.awk
#@brief   SCIP Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Alexander Martin
#@author  Timo Berthold
#@author  Robert Waniek
#@author  Gregor Hendel
#@author  Christian Puchert
#@author  Gerald Gamrath
#@author  Martin Bergner
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
   timegeomshift = 10.0;
   nodegeomshift = 100.0;
   pricegeomshift = 5.0;
   lpgeomshift = 2.0;
   pricecallsgeomshift = 100;
   priceprobsgeomshift = 100;
   pricevarsgeomshift = 100;
   sblpgeomshift = 0.0;
   pavshift = 0.0;
   onlyinsolufile = 0;          # should only instances be reported that are included in the .solu file?
   onlyintestfile = 0;          # should only instances be reported that are included in the .test file?  TEMPORARY HACK!
   onlypresolvereductions = 0;  # should only instances with presolve reductions be shown?
   useshortnames = 1;           # should problem name be truncated to fit into column?
   writesolufile = 0;           # should a solution file be created from the results
   printsoltimes = 0;           # should the times until first and best solution be shown
   checksol = 1;                # should the solution check of SCIP be parsed and counted as a fail if best solution is not feasible?
   NEWSOLUFILE = "new_solufile.solu";
   infty = +1e+20;
   headerprinted = 0;

   nprobs = 0;
   ndbprobs = 0;
   sbab = 0;
   slp = 0;
   ssim = 0;
   ssblp = 0;
   stottime = 0.0;
   stimetofirst = 0.0;
   stimetobest = 0.0;
   spricetime = 0.0;
   slptime = 0.0;
   spricecalls = 0.0;
   spriceprobs = 0.0;
   spricevars = 0.0;
   nodegeom = 0.0;
   timegeom = 0.0;
   pricegeom = 0.0;
   pricecallsgeom = 0.0;
   priceprobsgeom = 0.0;
   pricevarsgeom = 0.0;
   lpgeom = 0.0;
   timetofirstgeom = 0.0;
   timetobestgeom = 0.0;
   sblpgeom = 0.0;
   conftimegeom = 0.0;
   basictimegeom = 0.0;
   overheadtimegeom = 0.0;
   shiftednodegeom = nodegeomshift;
   shiftedtimegeom = timegeomshift;
   shiftedpricegeom = pricegeomshift;
   shiftedpricecallsgeom = pricecallsgeomshift;
   shiftedpriceprobsgeom = priceprobsgeomshift;
   shiftedpricevarsgeom = pricevarsgeomshift;
   shiftedlpgeom = pricegeomshift;
   shiftedsblpgeom = sblpgeomshift;
   shiftedconftimegeom = timegeomshift;
   shiftedbasictimegeom = timegeomshift;
   shiftedoverheadtimegeom = timegeomshift;
   shiftedtimetofirstgeom = timegeomshift;
   shiftedtimetobestgeom = timegeomshift;
   timeouttime = 0.0;
   timeouts = 0;
   failtime = 0.0;
   fail = 0;
   pass = 0;
   settings = "default";
   lpsname = "?";
   lpsversion = "?";
   scipversion = "?";
   gcgversion = "?";
   githash = "?";
   gcggithash = "?";
   conftottime = 0.0;
   overheadtottime = 0.0;
   timelimit = 0.0;

   #sepa base stuff 
   basetimegeomshift = 10;
   basecallsgeomshift = 100;
   basefoundgeomshift = 100;
   baseappliedgeomshift = 100;
   baseconvexgeomshift = 1;
   basediffstartgeomshift = 10;
   basediffendgeomshift = 10;

   cgmipbasegeomshift = 1;
   cgmipmastergeomshift = 1;

   cliquebasegeomshift = 1;
   cliquemastergeomshift = 1;
   cmirbasegeomshift = 1;
   cmirmastergeomshift = 1;
   flowcoverbasegeomshift = 1;
   flowcovermastergeomshift = 1;
   gomorybasegeomshift = 1;
   gomorymastergeomshift = 1;
   impliedboundsbasegeomshift = 1;
   impliedboundsmastergeomshift = 1;
   mcfbasegeomshift = 1;
   mcfmastergeomshift = 1;
   oddcyclebasegeomshift = 1;
   oddcyclemastergeomshift = 1;
   strongcgbasegeomshift = 1;
   strongcgmastergeomshift = 1;
   zerohalfbasegeomshift = 1;
   zerohalfmastergeomshift = 1;


   mastertimegeomshift = 10;
   mastercallsgeomshift = 100;
   masterfoundgeomshift = 100;
   masterappliedgeomshift = 100;

   rootdbgeomshift = 1;

   sbasetime = 0.0;
   sbasecalls = 0.0;
   sbasefound = 0.0;
   sbaseapplied = 0.0;
   sbaseconvex = 0.0;
   sbasel1diff = 0.0;
   basetimegeom = 0.0;
   basecallsgeom = 0.0;
   basefoundgeom = 0.0;
   baseappliedgeom = 0.0;
   baseconvexgeom = 0.0;
   basediffstartgeom = 0.0;
   basediffendgeom = 0.0;
   shiftedbasegeom = basetimegeomshift;
   shiftedbasecallsgeom = basecallsgeomshift;
   shiftedbasefoundgeom = basefoundgeomshift;
   shiftedbaseappliedgeom = baseappliedgeomshift;
   shiftedbaseconvexgeom = baseconvexgeomshift;
   shiftedbasediffstartgeom = basediffstartgeomshift;
   shiftedbasediffendgeom = basediffendgeomshift;

   scgmipbase = 0.0;
   cgmipbasegeom = 0.0;
   shiftedcgmipbasegeom = cgmipbasegeomshift;

   scgmipmaster = 0.0;
   cgmipmastergeom = 0.0;
   shiftedcgmipmastergeom = cgmipmastergeomshift;

   scliquebase = 0.0;
   cliquebasegeom = 0.0;
   shiftedcliquebasegeom = cliquebasegeomshift;
   scliquemaster = 0.0;
   cliquemastergeom = 0.0;
   shiftedcliquemastergeom = cliquemastergeomshift;
   scmirbase = 0.0;
   cmirbasegeom = 0.0;
   shiftedcmirbasegeom = cmirbasegeomshift;
   scmirmaster = 0.0;
   cmirmastergeom = 0.0;
   shiftedcmirmastergeom = cmirmastergeomshift;
   sflowcoverbase = 0.0;
   flowcoverbasegeom = 0.0;
   shiftedflowcoverbasegeom = flowcoverbasegeomshift;
   sflowcovermaster = 0.0;
   flowcovermastergeom = 0.0;
   shiftedflowcovermastergeom = flowcovermastergeomshift;
   sgomorybase = 0.0;
   gomorybasegeom = 0.0;
   shiftedgomorybasegeom = gomorybasegeomshift;
   sgomorymaster = 0.0;
   gomorymastergeom = 0.0;
   shiftedgomorymastergeom = gomorymastergeomshift;
   simpliedboundsbase = 0.0;
   impliedboundsbasegeom = 0.0;
   shiftedimpliedboundsbasegeom = impliedboundsbasegeomshift;
   simpliedboundsmaster = 0.0;
   impliedboundsmastergeom = 0.0;
   shiftedimpliedboundsmastergeom = impliedboundsmastergeomshift;
   smcfbase = 0.0;
   mcfbasegeom = 0.0;
   shiftedmcfbasegeom = mcfbasegeomshift;
   smcfmaster = 0.0;
   mcfmastergeom = 0.0;
   shiftedmcfmastergeom = mcfmastergeomshift;
   soddcyclebase = 0.0;
   oddcyclebasegeom = 0.0;
   shiftedoddcyclebasegeom = oddcyclebasegeomshift;
   soddcyclemaster = 0.0;
   oddcyclemastergeom = 0.0;
   shiftedoddcyclemastergeom = oddcyclemastergeomshift;
   sstrongcgbase = 0.0;
   strongcgbasegeom = 0.0;
   shiftedstrongcgbasegeom = strongcgbasegeomshift;
   sstrongcgmaster = 0.0;
   strongcgmastergeom = 0.0;
   shiftedstrongcgmastergeom = strongcgmastergeomshift;
   szerohalfbase = 0.0;
   zerohalfbasegeom = 0.0;
   shiftedzerohalfbasegeom = zerohalfbasegeomshift;
   szerohalfmaster = 0.0;
   zerohalfmastergeom = 0.0;
   shiftedzerohalfmastergeom = zerohalfmastergeomshift;


   smastertime = 0.0;
   smastercalls = 0.0;
   smasterfound = 0.0;
   smasterapplied = 0.0;
   mastertimegeom = 0.0;
   mastercallsgeom = 0.0;
   masterfoundgeom = 0.0;
   masterappliedgeom = 0.0;
   shiftedmastergeom = mastertimegeomshift;
   shiftedmastercallsgeom = mastercallsgeomshift;
   shiftedmasterfoundgeom = masterfoundgeomshift;
   shiftedmasterappliedgeom = masterappliedgeomshift;

   srootdb = 0.0;
   rootdbgeom = 0.0;
   shiftedrootdbgeom = rootdbgeomshift; 

   #initialize paver input file
   if( PAVFILE != "" ) {
      printf("* Trace Record Definition\n") > PAVFILE;
      printf("* InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime\n") > PAVFILE;
      printf("* NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables\n") > PAVFILE;
   }
}
/^IP\// {  # TEMPORARY HACK to parse .test file
   intestfile[$1] = 1;
}
/=opt=/  { solstatus[$2] = "opt"; sol[$2] = $3; }   # get optimum
/=inf=/  { solstatus[$2] = "inf"; }                 # problem infeasible (no feasible solution exists)
/=best=/ { solstatus[$2] = "best"; sol[$2] = $3; }  # get best known solution value
/=unkn=/ { solstatus[$2] = "unkn"; }                # no feasible solution known
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

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];
   vars = 0;
   binvars = 0;
   intvars = 0;
   implvars = 0;
   contvars = 0;
   cons = 0;
   lincons = 0;
   origvars = 0;
   origcons = 0;
   objsense = 0;
   timeout = 0;
   feasible = 0;
   pb = +infty;
   firstpb = +infty;
   db = -infty;
   rootdb = -infty;
   simpiters = 0;
   bbnodes = 0;
   primlps = 0;
   primiter = 0;
   duallps = 0;
   dualiter = 0;
   sblps = 0;
   sbiter = 0;
   tottime = 0.0;
   timetofirst = -1.0;
   timetobest = -1.0;
   inconflict = 0;
   inconstime = 0;
   confclauses = 0;
   confliterals = 0.0;
   conftime = 0.0;
   pricetime = 0.0;
   lptime = 0.0;
   npriceprobs = 0;
   overheadtime = 0.0;
   aborted = 1;
   readerror = 0;
   gapreached = 0;
   sollimitreached = 0;
   memlimitreached = 0;
   nodelimitreached = 0;
   starttime = 0.0;
   endtime = 0.0;
   timelimit = 0.0;
   inoriginalprob = 1;
   inmasterprob = 1;
   incons = 0;
   valgrinderror = 0;
   valgrindleaks = 0;
   bestsolfeas = 1;

   # sepa base stuff
   basetime = 0.0;
   mastertime = 0.0;
   insepas = 1;
   insepabase = 0;
   insepamaster = 0;
}

/@03/ { starttime = $2; }
/@04/ { endtime = $2; }

/^GCG version/ {
   # get GCG version
   gcgversion = $3;
   
      # get git hash 
   if( $(NF-1) == "[GitHash:" ) {
      split($NF, v, "]");
      gcggithash = v[1];
   }
}

/^SCIP version/ {
   # get SCIP version 
   scipversion = $3; 

   # get name of LP solver
   if( $13 == "SoPlex" )
      lpsname = "spx";
   else if( $13 == "CPLEX" )
      lpsname = "cpx";
   else if( $13 == "NONE" )
      lpsname = "none";
   else if( $13 == "Clp" )
      lpsname = "clp";
   else if( $13 == "MOSEK" )
      lpsname = "msk";
   else if( $13 == "Gurobi" )
      lpsname = "grb";
   else if( $13 == "NONE" )
      lpsname = "none";
   else if( $13 == "QSopt" )
      lpsname = "qso";
#   else if( $13 == "???" )
#      lpsname = "xprs";

    # get LP solver version 
   if( NF >= 14 ) {
      split($14, v, "]");
      lpsversion = v[1];
   }

   # get git hash 
   if( $(NF-1) == "[GitHash:" ) {
      split($NF, v, "]");
      githash = v[1];
   }
}
/^SCIP> SCIP> / { $0 = substr($0, 13, length($0)-12); }
/^SCIP> / { $0 = substr($0, 7, length($0)-6); }
/^loaded parameter file/ { settings = $4; sub(/<.*settings\//, "", settings); sub(/\.set>/, "", settings); }
/^GCG> GCG> loaded parameter file/ { settings = $6; sub(/<.*settings\//, "", settings); sub(/\.set>/, "", settings); }
/^parameter <limits\/time> set to/ { timelimit = $5; }
/^limits\/time =/ { timelimit = $3; }
#
# get objective sense
#
/^  Objective sense  :/ {
   if ( $4 == "minimize" )
      objsense = 1;
   if ( $4 == "maximize" )
      objsense = -1;
   # objsense is 0 otherwise
}
#
# problem: master or original?
#
/^Original Program statistics:/ { inmasterprob = 0; inoriginalprob = 1; }
/^Master Program statistics:/ { inmasterprob = 1; inoriginalprob = 0; }
#
# conflict analysis
#
/^Conflict Analysis  :/ { inconflict = 1; }
/^  propagation      :/ {
   if( inconflict == 1 ) {
      conftime += $3; #confclauses += $5 + $7; confliterals += $5 * $6 + $7 * $8;
   }
}
/^  infeasible LP    :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  strong branching :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  pseudo solution  :/ {
   if( inconflict == 1 ) {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  applied globally :/ {
   if( inconflict == 1 ) {
      confclauses += $7; confliterals += $7 * $8;
   }
}
/^  applied locally  :/ {
   if( inconflict == 1 ) {
      confclauses += $7; confliterals += $7 * $8;
   }
}
#
# pricing
#
/^  gcg              :/ {
   if( inmasterprob )
   {
      pricetime = $3;
      pricecall = $5;
      pricevars = $6;
   }
}
#
# sepa base
#
/^Separators         :/ { 
   
   inconflict = 0; 
   insepas = 1;			
}

/^  base             :/ {
   if( inmasterprob && insepas )
   {

      basetime = $3;
      basecalls = $5;
      basefound = $8;
      baseapplied = $9;
   }
}



#
# sepa master
#
/^  master           :/ {
   if( inmasterprob && insepas )
   {
      mastertime = $3;
      mastercalls = $5;
      masterfound = $8;
      masterapplied = $9;
   }
}

/^SepaBase:/ {
 
   insepabase = 1;
   insepamaster = 0;
   baseconvex = $11;
   basediffstart = $12;
   basediffend = $13;   
}

/^SepaMaster/ {

   insepabase = 0;
   insepamaster = 1;
}

/^cgmip / {
   if(insepabase)
   {
      cgmipbase = $2;
   }
   else if(insepamaster)
   {
      cgmipmaster = $2;
   }
}

/^clique / {
   if(insepabase)
   {
      cliquebase = $2;
   }
   else if(insepamaster)
   {
      cliquemaster = $2;
   }
}

/^cmir / {
   if(insepabase)
   {
      cmirbase = $2;
   }
   else if(insepamaster)
   {
      cmirmaster = $2;
   }
}

/^flowcover / {
   if(insepabase)
   {
      flowcoverbase = $2;
   }
   else if(insepamaster)
   {
      flowcovermaster = $2;
   }
}

/^gomory / {
   if(insepabase)
   {
      gomorybase = $2;
   }
   else if(insepamaster)
   {
      gomorymaster = $2;
   }
}

/^impliedbounds / {
   if(insepabase)
   {
      impliedboundsbase = $2;
   }
   else if(insepamaster)
   {
      impliedboundsmaster = $2;
   }
}

/^mcf / {
   if(insepabase)
   {
      mcfbase = $2;
   }
   else if(insepamaster)
   {
      mcfmaster = $2;
   }
}

/^oddcycle / {
   if(insepabase)
   {
      oddcyclebase = $2;
   }
   else if(insepamaster)
   {
      oddcyclemaster = $2;
   }
}
/^strongcg / {
   if(insepabase)
   {
      strongcgbase = $2;
   }
   else if(insepamaster)
   {
      strongcgmaster = $2;
   }
}

/^zerohalf / {
   if(insepabase)
   {
      zerohalfbase = $2;
   }
   else if(insepamaster)
   {
      zerohalfmaster = $2;
   }
}


/^Pricers            :/ { insepas = 0; }


/^Constraint Timings :/ { inconstime = 1; }
#/^  logicor          :/ { if( inconstime == 1 ) { overheadtime += $3; } }
/^Propagators        :/ { inconstime = 0; }
/^  switching time   :/ { overheadtime += $4; }
#
# problem size
#
/^Presolved Problem  :/ { inoriginalprob = 0; }
/^  Variables        :/ {
   if( !inmasterprob )
   {
      if( inoriginalprob )
         origvars = $3;
      else
      {
         vars = $3;
         intvars = $6;
         implvars = $8;
         contvars = $11;
         binvars = vars - intvars - implvars - contvars;
      }
   }
}
/^  Constraints      :/ {
   if( !inmasterprob )
   {
      if( inoriginalprob )
         origcons = $3;
      else
         cons = $3;
   }
   else
   {
#      if( inoriginalprob )
         mcons = $3;
   }
}


/^Chosen decomposition/ {
   blocks = $4;
   rel = $4;
}

/^Matrix has / {
   blocks = $3;
   rel = $6;
}
/^GCG                : Performing Dantzig-Wolfe with [0-9]+ blocks./ {
   blocks = $6;
}
/^  mip              :/ {
   npriceprobs = $4 + $6;
}

#
# count number of linear constraints 
#
/^Constraints        :/ {
   incons = 1;
   lincons = 0;
}
/^  knapsack         :/ {
   if( incons == 1 ) {
      n  = split ($3, a, "+");
      lincons += a[1];
   }
}
/^  setppc           :/ {
   if( incons == 1 ) {
      n  = split ($3, a, "+");
      lincons += a[1];
   }
}
/^  linear           :/ { 
   if( incons  == 1 ) {
      n  = split ($3, a, "+");
      lincons += a[1];
   }
}
/^  logicor          :/ { 
   if( incons == 1 ) {
      n  = split ($3, a, "+");
      lincons += a[1];
   }
}
/^  varbound         :/ { 
   if( incons == 1 ) {
      n  = split ($3, a, "+");
      lincons += a[1];
   }
}
/^Constraint Timings :/ {
   incons = 0;
}

#
# solution
#
/^Original Problem   : no problem exists./ { readerror = 1; }
/^SCIP Status        :/ { aborted = 0; }
/solving was interrupted/ { if( inoriginalprob ) timeout = 1; }
/gap limit reached/ { if( inoriginalprob ) gapreached = 1; }
/solution limit reached/ { if( inoriginalprob ) sollimitreached = 1; }
/memory limit reached/ { if( inoriginalprob ) memlimitreached = 1; }
/node limit reached/ { if( inoriginalprob ) nodelimitreached = 1; }
/problem is solved/ { if( inoriginalprob ) timeout = 0; }
/best solution is not feasible in original problem/  { if( inoriginalprob ) bestsolfeas = 0; }

/^  First Solution   :/ {
   timetofirst = $11;
   firstpb = $4;
}
/^  Primal Bound     :/ {
   if( $4 == "infeasible" ) {
      pb = +infty;
      db = +infty;
      feasible = 0;
   }
   else if( $4 == "-" ) {
      pb = +infty;
      feasible = 0;
   }
   else {
      pb = $4;
      feasible = 1;
      timetobest = $11;
   }
}
/^Dual Bound         :/ { 
   if( $4 != "-" ) 
      db = $4;
}
/^  Final Dual Bound :/ {
   if( $5 != "-" )
      rootdb = $5;
   else
       rootdb = db;  # SCIP most likely finished during root node, perhaps due to a solution limit. the rootdb is NOT printed then, but needed later
}
#
# iterations
#
/^  primal LP        :/ { simpiters += $6; lptime += $4 }
/^  dual LP          :/ { simpiters += $6; lptime += $4 }
/^  barrier LP       :/ { simpiters += $6; }
/^  nodes \(total\)    :/ { bbnodes = $4 }
/^  primal LP        :/ { primlps = $5; primiter = $6; }
/^  dual LP          :/ { duallps = $5; dualiter = $6; }
/^  strong branching :/ { sblps = $5; sbiter = $6; }
#
# time
#
/^Solving Time       :/ { tottime = $4 } # for older scip version ( < 2.0.1.3 )
/^  solving          :/ { tottime = $3 } 
#
# valgrind check
#
/^==[0-9]*== ERROR SUMMARY:/       { valgrinderror = $4 }
/^==[0-9]*==    definitely lost:/  { valgrindleaks += $4 }
/^==[0-9]*==    indirectly lost:/  { valgrindleaks += $4 }
/^==[0-9]*==    possibly lost:/    { valgrindleaks += $4 }
#
# solver status overview (in order of priority): 
# 1) solver broke before returning solution => abort
# 2) solver cut off the optimal solution (solu-file-value is not between primal and dual bound) => fail
#    (especially if problem is claimed to be solved but solution is not the optimal solution)
# 3) solver solved problem with the value in solu-file (if existing) => ok
# 4) solver solved problem which has no (optimal) value in solu-file => solved
#    (since we here don't detect the direction of optimization, it is possible 
#     that a solver claims an optimal solution which contradicts a known feasible solution)
# 5) solver found solution better than known best solution (or no solution was known so far) => better
# 7) solver reached gaplimit or limit of number of solutions => gaplimit, sollimit
# 8) solver reached any other limit (like time or nodes) => timeout
# 9) otherwise => unknown
#
/^=ready=/ {

   #since the header depends on the parameter printsoltimes and settings it is no longer possible to print it in the BEGIN section
   if( !headerprinted ) {
      ntexcolumns = 8 + (2 * printsoltimes);


      if (TEXFILE != "") {
         #print header of tex file table
         printf("\\documentclass[leqno]{article}\n")                      >TEXFILE;
         printf("\\usepackage{a4wide}\n")                                 >TEXFILE;
         printf("\\usepackage{amsmath,amsfonts,amssymb,booktabs}\n")      >TEXFILE;
         printf("\\usepackage{supertabular}\n")                           >TEXFILE;
         printf("\\usepackage{rotating}\n")                               >TEXFILE;
         printf("\\usepackage{lscape}\n")                                 >TEXFILE;
         printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
         printf("\\begin{document}\n\n")                                  >TEXFILE;
         printf("\\begin{center}\n")                                      >TEXFILE;
         printf("\\tiny\n")                                               >TEXFILE;
         printf("\\setlength{\\tabcolsep}{2pt}\n")                        >TEXFILE;
         printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n")    >TEXFILE;
         printf("\\tablehead{\n\\toprule\n")                              >TEXFILE;
         printf("& \\multicolumn{2}{c}{cgmip} & \\multicolumn{2}{c}{clique} & \\multicolumn{2}{c}{cmir} & \\multicolumn{2}{c}{flowcover} & \\multicolumn{2}{c}{gomory} & \\multicolumn{2}{c}{impliedbounds} & \\multicolumn{2}{c}{mcf} & \\multicolumn{2}{c}{oddcycle} & \\multicolumn{2}{c}{strongcg} & \\multicolumn{2}{c}{zerohalf} \\\\")      >TEXFILE;
         printf("\n\\midrule\n")                              >TEXFILE;
         printf("Name             &   b & m & b & m & b & m & b & m & b & m & b & m & b & m & b & m & b & m & b & m ") > TEXFILE;
         printf("\\\\\n") > TEXFILE;
         printf("\\midrule\n}\n")                                         >TEXFILE;
         printf("\\tabletail{\n\\midrule\n")                              >TEXFILE;
         printf("\\multicolumn{%d}{r} \\; continue next page \\\\\n", ntexcolumns) >TEXFILE;
         printf("\\bottomrule\n}\n")                                      >TEXFILE;
         printf("\\tablelasttail{\\bottomrule}\n")                        >TEXFILE;
         printf("\\tablecaption{SCIP with %s settings}\n",settings)       >TEXFILE;
	 printf("\\hspace{-1cm}\n")					  >TEXFILE;
         printf("\\begin{supertabular}{@{\\extracolsep{\\fill}}lrlrlrlrlrlrlrlrlrlrlr") >TEXFILE;
         printf("@{}}\n") > TEXFILE;
      }
      
      #print header of table when this regular expression is matched for the first time
      tablehead1 = "-----------+----+--- Original --+-- Presolved --+------+------+------+--------------+--------------+------+------------ Pricing ----------+-------+--------+-------+-------+";
      tablehead2 = "Name       |     cgmip  |    clique |     cmir  | flowcover |   gomory  | impbounds |    mcf    |  oddcycle |  strongcg |  zerohalf";
      tablehead3 = "-----------+--------+--------+--------+-----------+--------+---------------+--------+----------+----------+---------";
      
      tablehead1 = tablehead1"--------\n";
      tablehead2 = tablehead2"       \n";
      tablehead3 = tablehead3"--------\n";
   
      printf(tablehead1);
      printf(tablehead2);
      printf(tablehead3);

      headerprinted = 1;
   }

   if( (!onlyinsolufile || solstatus[prob] != "") &&
       (!onlyintestfile || intestfile[filename]) ) {

      #avoid problems when comparing floats and integer (make everything float)
      temp = pb;
      pb = 1.0*temp;
      temp = db;
      db = 1.0*temp;
      temp = rootdb;
      rootdb = 1.0*temp;
      
      # if objsense could not be determined so far (output is maybe too old)
      if ( objsense == 0 )
      {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

	 # firstpb and rootdb are used to detect the direction of optimization (min or max)
	 if( timetofirst < 0.0 )
	    temp = pb;
	 else
	    temp = firstpb;
	 firstpb = 1.0*temp;

	 if ( firstpb - rootdb > max(abstol,reltol) )
	    objsense = 1;   # minimize
	 else
	    objsense = -1;  # maximize
      }
      
      # modify primal bound for maximization problems without primal solution
      if ( objsense == -1 && pb >= +infty )
	 pb = -1.0 * pb;

      # modify dual bound for infeasible maximization problems
      if ( objsense == -1 && db >= +infty )
	 db = -1.0 * db;

      # modify root dual bound for infeasible maximization problems
      if ( objsense == -1 && rootdb >= +infty )
	 rootdb = -1.0 * rootdb;

      nprobs++;

      optimal = 0;
      markersym = "\\g";
      if( abs(pb - db) < 1e-06 && pb < infty ) {
         gap = 0.0;
         optimal = 1;
         markersym = "  ";
      }
      else if( abs(db) < 1e-06 )
         gap = -1.0;
      else if( abs(pb) < 1e-06 )
         gap = -1.0;
      else if( pb*db < 0.0 )
         gap = -1.0;
      else if( abs(db) >= +infty )
         gap = -1.0;
      else if( abs(pb) >= +infty )
         gap = -1.0;
      else
         gap = 100.0*abs((pb-db)/min(abs(db),abs(pb)));
      if( gap < 0.0 )
         gapstr = "  --  ";
      else if( gap < 1e+04 )
         gapstr = sprintf("%6.1f", gap);
      else
         gapstr = " Large";
      
      if( vars == 0 )
         probtype = "--";
      else if( lincons < cons )
         probtype = "CIP";
      else if( binvars == 0 && intvars == 0 )
         probtype = "LP";
      else if( contvars == 0 ) {
         if( intvars == 0 && implvars == 0 )
            probtype = "BP";
         else
            probtype = "IP";
      }
      else {
         if( intvars == 0 )
            probtype = "MBP";
         else
            probtype = "MIP";
      }

      if( aborted && endtime - starttime > timelimit && timelimit > 0.0 ) {
         timeout = 1;
         aborted = 0;
         tottime = endtime - starttime;
      }
      else if( gapreached || sollimitreached || memlimitreached || nodelimitreached )
         timeout = 0;

      if( aborted && tottime == 0.0 )
         tottime = timelimit;
      if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

      if( aborted || timetobest < 0.0 ) {
         timetofirst = tottime;
         timetobest = tottime;
      }

      lps = primlps + duallps;
      simplex = primiter + dualiter;
      stottime += tottime;
      stimetofirst += timetofirst;
      stimetobest += timetobest;
      spricetime += pricetime;
      spricecalls += pricecall;
      spriceprobs += npriceprobs;
      spricevars += pricevars;
      slptime += lptime;
      sbab += bbnodes;
      slp += lps;
      ssim += simplex;
      ssblp += sblps;
      conftottime += conftime;
      overheadtottime += overheadtime;
      basictime = tottime - conftime - overheadtime;

      #sepa base stuff 
      sbasetime += basetime;
      sbasecalls += basecalls;
      sbasefound += basefound;
      sbaseapplied += baseapplied;
      sbaseconvex += baseconvex;
      sbasediffstart += basediffstart;
      sbasediffend += basediffend;

      scgmipbase += cgmipbase;
      cgmipbasegeom = cgmipbasegeom^((nprobs-1)/nprobs) * max(cgmipbase, 1.0)^(1.0/nprobs);
      shiftedcgmipbasegeom = shiftedcgmipbasegeom^((nprobs-1)/nprobs) * max(cgmipbase+cgmipbasegeomshift, 1.0)^(1.0/nprobs);

      scgmipmaster += cgmipmaster;
      cgmipmastergeom = cgmipmastergeom^((nprobs-1)/nprobs) * max(cgmipmaster, 1.0)^(1.0/nprobs);
      shiftedcgmipmastergeom = shiftedcgmipmastergeom^((nprobs-1)/nprobs) * max(cgmipmaster+cgmipmastergeomshift, 1.0)^(1.0/nprobs);

      scliquebase += cliquebase;
      cliquebasegeom = cliquebasegeom^((nprobs-1)/nprobs) * max(cliquebase, 1.0)^(1.0/nprobs);
      shiftedcliquebasegeom = shiftedcliquebasegeom^((nprobs-1)/nprobs) * max(cliquebase+cliquebasegeomshift, 1.0)^(1.0/nprobs);
      scliquemaster += cliquemaster;
      cliquemastergeom = cliquemastergeom^((nprobs-1)/nprobs) * max(cliquemaster, 1.0)^(1.0/nprobs);
      shiftedcliquemastergeom = shiftedcliquemastergeom^((nprobs-1)/nprobs) * max(cliquemaster+cliquemastergeomshift, 1.0)^(1.0/nprobs);
      scmirbase += cmirbase;
      cmirbasegeom = cmirbasegeom^((nprobs-1)/nprobs) * max(cmirbase, 1.0)^(1.0/nprobs);
      shiftedcmirbasegeom = shiftedcmirbasegeom^((nprobs-1)/nprobs) * max(cmirbase+cmirbasegeomshift, 1.0)^(1.0/nprobs);
      scmirmaster += cmirmaster;
      cmirmastergeom = cmirmastergeom^((nprobs-1)/nprobs) * max(cmirmaster, 1.0)^(1.0/nprobs);
      shiftedcmirmastergeom = shiftedcmirmastergeom^((nprobs-1)/nprobs) * max(cmirmaster+cmirmastergeomshift, 1.0)^(1.0/nprobs);
      sflowcoverbase += flowcoverbase;
      flowcoverbasegeom = flowcoverbasegeom^((nprobs-1)/nprobs) * max(flowcoverbase, 1.0)^(1.0/nprobs);
      shiftedflowcoverbasegeom = shiftedflowcoverbasegeom^((nprobs-1)/nprobs) * max(flowcoverbase+flowcoverbasegeomshift, 1.0)^(1.0/nprobs);
      sflowcovermaster += flowcovermaster;
      flowcovermastergeom = flowcovermastergeom^((nprobs-1)/nprobs) * max(flowcovermaster, 1.0)^(1.0/nprobs);
      shiftedflowcovermastergeom = shiftedflowcovermastergeom^((nprobs-1)/nprobs) * max(flowcovermaster+flowcovermastergeomshift, 1.0)^(1.0/nprobs);
      sgomorybase += gomorybase;
      gomorybasegeom = gomorybasegeom^((nprobs-1)/nprobs) * max(gomorybase, 1.0)^(1.0/nprobs);
      shiftedgomorybasegeom = shiftedgomorybasegeom^((nprobs-1)/nprobs) * max(gomorybase+gomorybasegeomshift, 1.0)^(1.0/nprobs);
      sgomorymaster += gomorymaster;
      gomorymastergeom = gomorymastergeom^((nprobs-1)/nprobs) * max(gomorymaster, 1.0)^(1.0/nprobs);
      shiftedgomorymastergeom = shiftedgomorymastergeom^((nprobs-1)/nprobs) * max(gomorymaster+gomorymastergeomshift, 1.0)^(1.0/nprobs);
      simpliedboundsbase += impliedboundsbase;
      impliedboundsbasegeom = impliedboundsbasegeom^((nprobs-1)/nprobs) * max(impliedboundsbase, 1.0)^(1.0/nprobs);
      shiftedimpliedboundsbasegeom = shiftedimpliedboundsbasegeom^((nprobs-1)/nprobs) * max(impliedboundsbase+impliedboundsbasegeomshift, 1.0)^(1.0/nprobs);
      simpliedboundsmaster += impliedboundsmaster;
      impliedboundsmastergeom = impliedboundsmastergeom^((nprobs-1)/nprobs) * max(impliedboundsmaster, 1.0)^(1.0/nprobs);
      shiftedimpliedboundsmastergeom = shiftedimpliedboundsmastergeom^((nprobs-1)/nprobs) * max(impliedboundsmaster+impliedboundsmastergeomshift, 1.0)^(1.0/nprobs);
      smcfbase += mcfbase;
      mcfbasegeom = mcfbasegeom^((nprobs-1)/nprobs) * max(mcfbase, 1.0)^(1.0/nprobs);
      shiftedmcfbasegeom = shiftedmcfbasegeom^((nprobs-1)/nprobs) * max(mcfbase+mcfbasegeomshift, 1.0)^(1.0/nprobs);
      smcfmaster += mcfmaster;
      mcfmastergeom = mcfmastergeom^((nprobs-1)/nprobs) * max(mcfmaster, 1.0)^(1.0/nprobs);
      shiftedmcfmastergeom = shiftedmcfmastergeom^((nprobs-1)/nprobs) * max(mcfmaster+mcfmastergeomshift, 1.0)^(1.0/nprobs);
      soddcyclebase += oddcyclebase;
      oddcyclebasegeom = oddcyclebasegeom^((nprobs-1)/nprobs) * max(oddcyclebase, 1.0)^(1.0/nprobs);
      shiftedoddcyclebasegeom = shiftedoddcyclebasegeom^((nprobs-1)/nprobs) * max(oddcyclebase+oddcyclebasegeomshift, 1.0)^(1.0/nprobs);
      soddcyclemaster += oddcyclemaster;
      oddcyclemastergeom = oddcyclemastergeom^((nprobs-1)/nprobs) * max(oddcyclemaster, 1.0)^(1.0/nprobs);
      shiftedoddcyclemastergeom = shiftedoddcyclemastergeom^((nprobs-1)/nprobs) * max(oddcyclemaster+oddcyclemastergeomshift, 1.0)^(1.0/nprobs);
      sstrongcgbase += strongcgbase;
      strongcgbasegeom = strongcgbasegeom^((nprobs-1)/nprobs) * max(strongcgbase, 1.0)^(1.0/nprobs);
      shiftedstrongcgbasegeom = shiftedstrongcgbasegeom^((nprobs-1)/nprobs) * max(strongcgbase+strongcgbasegeomshift, 1.0)^(1.0/nprobs);
      sstrongcgmaster += strongcgmaster;
      strongcgmastergeom = strongcgmastergeom^((nprobs-1)/nprobs) * max(strongcgmaster, 1.0)^(1.0/nprobs);
      shiftedstrongcgmastergeom = shiftedstrongcgmastergeom^((nprobs-1)/nprobs) * max(strongcgmaster+strongcgmastergeomshift, 1.0)^(1.0/nprobs);
      szerohalfbase += zerohalfbase;
      zerohalfbasegeom = zerohalfbasegeom^((nprobs-1)/nprobs) * max(zerohalfbase, 1.0)^(1.0/nprobs);
      shiftedzerohalfbasegeom = shiftedzerohalfbasegeom^((nprobs-1)/nprobs) * max(zerohalfbase+zerohalfbasegeomshift, 1.0)^(1.0/nprobs);
      szerohalfmaster += zerohalfmaster;
      zerohalfmastergeom = zerohalfmastergeom^((nprobs-1)/nprobs) * max(zerohalfmaster, 1.0)^(1.0/nprobs);
      shiftedzerohalfmastergeom = shiftedzerohalfmastergeom^((nprobs-1)/nprobs) * max(zerohalfmaster+zerohalfmastergeomshift, 1.0)^(1.0/nprobs);



      basetimegeom = basetimegeom^((nprobs-1)/nprobs) * max(basetime, 1.0)^(1.0/nprobs);
      basecallsgeom = basecallsgeom^((nprobs-1)/nprobs) * max(basecalls, 1.0)^(1.0/nprobs);
      basefoundgeom = basefoundgeom^((nprobs-1)/nprobs) * max(basefound, 1.0)^(1.0/nprobs);
      baseappliedgeom = baseappliedgeom^((nprobs-1)/nprobs) * max(baseapplied, 1.0)^(1.0/nprobs);
      baseconvexgeom = baseconvexgeom^((nprobs-1)/nprobs) * max(baseconvex, 1.0)^(1.0/nprobs);
      basediffstartgeom = basediffstartgeom^((nprobs-1)/nprobs) * max(basediffstart, 1.0)^(1.0/nprobs);
      basediffendgeom = basediffendgeom^((nprobs-1)/nprobs) * max(basediffend, 1.0)^(1.0/nprobs);
 


      shiftedbasetimegeom = shiftedbasetimegeom^((nprobs-1)/nprobs) * max(basetime+basetimegeomshift, 1.0)^(1.0/nprobs);
      shiftedbasecallsgeom = shiftedbasecallsgeom^((nprobs-1)/nprobs) * max(basecalls+basecallsgeomshift, 1.0)^(1.0/nprobs);
      shiftedbasefoundgeom = shiftedbasefoundgeom^((nprobs-1)/nprobs) * max(basefound+basefoundgeomshift, 1.0)^(1.0/nprobs);
      shiftedbaseappliedgeom = shiftedbaseappliedgeom^((nprobs-1)/nprobs) * max(baseapplied+baseappliedgeomshift, 1.0)^(1.0/nprobs);
      shiftedbaseconvexgeom = shiftedbaseconvexgeom^((nprobs-1)/nprobs) * max(baseconvex+baseconvexgeomshift, 1.0)^(1.0/nprobs);
      shiftedbasediffstartgeom = shiftedbasediffstartgeom^((nprobs-1)/nprobs) * max(basediffstart+basediffstartgeomshift, 1.0)^(1.0/nprobs);
      shiftedbasediffendgeom = shiftedbasediffendgeom^((nprobs-1)/nprobs) * max(basediffend+basediffendgeomshift, 1.0)^(1.0/nprobs);




      smastertime += mastertime;
      smastercalls += mastercalls;
      smasterfound += masterfound;
      smasterapplied += masterapplied;

      mastertimegeom = mastertimegeom^((nprobs-1)/nprobs) * max(mastertime, 1.0)^(1.0/nprobs);
      mastercallsgeom = mastercallsgeom^((nprobs-1)/nprobs) * max(mastercalls, 1.0)^(1.0/nprobs);
      masterfoundgeom = masterfoundgeom^((nprobs-1)/nprobs) * max(masterfound, 1.0)^(1.0/nprobs);
      masterappliedgeom = masterappliedgeom^((nprobs-1)/nprobs) * max(masterapplied, 1.0)^(1.0/nprobs);
 
      shiftedmastertimegeom = shiftedmastertimegeom^((nprobs-1)/nprobs) * max(mastertime+mastertimegeomshift, 1.0)^(1.0/nprobs);
      shiftedmastercallsgeom = shiftedmastercallsgeom^((nprobs-1)/nprobs) * max(mastercalls+mastercallsgeomshift, 1.0)^(1.0/nprobs);
      shiftedmasterfoundgeom = shiftedmasterfoundgeom^((nprobs-1)/nprobs) * max(masterfound+masterfoundgeomshift, 1.0)^(1.0/nprobs);
      shiftedmasterappliedgeom = shiftedmasterappliedgeom^((nprobs-1)/nprobs) * max(masterapplied+masterappliedgeomshift, 1.0)^(1.0/nprobs);

      #if( abs(pb - rootdb) < 1e-06 && pb < infty ) {
      #   rootdb = 0.0;
      #}
      #else if( abs(rootdb) < 1e-06 )
      #   rootdb = -1.0;
      #else if( abs(pb) < 1e-06 )
      #   rootdb = -1.0;
      #else if( pb*rootdb < 0.0 )
      #   rootdb = -1.0;
      #else if( abs(rootdb) >= +infty )
      #   rootdb = -1.0;
      #else if( abs(pb) >= +infty )
      #   rootdb = -1.0;
      #else
      #   rootdb = 100.0*abs((pb-rootdb)/min(abs(rootdb),abs(pb)));

      #if( rootdb < 0.0 )
      #   rootdbstr = "  --  ";
      #else if( rootdb < 1e+04 )
      #   rootdbstr = sprintf("%6.1f", rootdb);
      #else
      #   rootdbstr = " Large";

      #if( rootdb >= 0.0 )
      #{ 
         ndbprobs++;
	 srootdb += rootdb;
         rootdbgeom = rootdbgeom^((ndbprobs-1)/ndbprobs) * max(rootdb, 1.0)^(1.0/ndbprobs);
         shiftedrootdbgeom = shiftedrootdbgeom^((ndbprobs-1)/ndbprobs) * max(rootdb+rootdbgeomshift, 1.0)^(1.0/ndbprobs);
      #}

      nodegeom = nodegeom^((nprobs-1)/nprobs) * max(bbnodes, 1.0)^(1.0/nprobs);
      sblpgeom = sblpgeom^((nprobs-1)/nprobs) * max(sblps, 1.0)^(1.0/nprobs);
      timegeom = timegeom^((nprobs-1)/nprobs) * max(tottime, 1.0)^(1.0/nprobs);
      pricegeom = pricegeom^((nprobs-1)/nprobs) * max(pricetime, 1.0)^(1.0/nprobs);
      pricecallsgeom = pricecallsgeom^((nprobs-1)/nprobs) * max(pricecall, 1.0)^(1.0/nprobs);
      priceprobsgeom = priceprobsgeom^((nprobs-1)/nprobs) * max(npriceprobs, 1.0)^(1.0/nprobs);
      pricevarsgeom = pricevarsgeom^((nprobs-1)/nprobs) * max(pricevars, 1.0)^(1.0/nprobs);
      lpgeom = lpgeom^((nprobs-1)/nprobs) * max(lptime, 1.0)^(1.0/nprobs);

      conftimegeom = conftimegeom^((nprobs-1)/nprobs) * max(conftime, 1.0)^(1.0/nprobs);
      overheadtimegeom = overheadtimegeom^((nprobs-1)/nprobs) * max(overheadtime, 1.0)^(1.0/nprobs);
      basictimegeom = basictimegeom^((nprobs-1)/nprobs) * max(basictime, 1.0)^(1.0/nprobs);

      shiftednodegeom = shiftednodegeom^((nprobs-1)/nprobs) * max(bbnodes+nodegeomshift, 1.0)^(1.0/nprobs);
      shiftedsblpgeom = shiftedsblpgeom^((nprobs-1)/nprobs) * max(sblps+sblpgeomshift, 1.0)^(1.0/nprobs);
      shiftedtimegeom = shiftedtimegeom^((nprobs-1)/nprobs) * max(tottime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedpricegeom = shiftedpricegeom^((nprobs-1)/nprobs) * max(pricetime+pricegeomshift, 1.0)^(1.0/nprobs);
      shiftedpricecallsgeom = shiftedpricecallsgeom^((nprobs-1)/nprobs) * max(pricecall+pricecallsgeomshift, 1.0)^(1.0/nprobs);
      shiftedpriceprobsgeom = shiftedpriceprobsgeom^((nprobs-1)/nprobs) * max(npriceprobs+priceprobsgeomshift, 1.0)^(1.0/nprobs);
      shiftedpricevarsgeom = shiftedpricevarsgeom^((nprobs-1)/nprobs) * max(pricevars+pricevarsgeomshift, 1.0)^(1.0/nprobs);
      shiftedlpgeom = shiftedlpgeom^((nprobs-1)/nprobs) * max(lptime+lpgeomshift, 1.0)^(1.0/nprobs);
      shiftedconftimegeom = shiftedconftimegeom^((nprobs-1)/nprobs) * max(conftime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedoverheadtimegeom = shiftedoverheadtimegeom^((nprobs-1)/nprobs) * max(overheadtime+timegeomshift, 1.0)^(1.0/nprobs);
      shiftedbasictimegeom = shiftedbasictimegeom^((nprobs-1)/nprobs) * max(basictime+timegeomshift, 1.0)^(1.0/nprobs);

      timetobestgeom = timetobestgeom^((nprobs-1)/nprobs) * max(timetobest,1.0)^(1.0/nprobs);
      timetofirstgeom = timetofirstgeom^((nprobs-1)/nprobs) * max(timetofirst,1.0)^(1.0/nprobs);
      shiftedtimetofirstgeom = shiftedtimetofirstgeom^((nprobs-1)/nprobs) * max(timetofirst + timegeomshift, 1.0)^(1.0/nprobs);
      shiftedtimetobestgeom = shiftedtimetobestgeom^((nprobs-1)/nprobs) * max(timetobest + timegeomshift, 1.0)^(1.0/nprobs);

      status = "";
      if( readerror ) {
         status = "readerror";
         failtime += tottime;
         fail++;
      }
      else if( aborted ) {
         status = "abort";
         failtime += tottime;
         fail++;
      }

      else if( checksol && !bestsolfeas ) {
         status = "fail";
         failtime += tottime;
         fail++;
      }
      else if( solstatus[prob] == "opt" ) {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

	 # objsense = 1 -> minimize; objsense = -1 -> maximize
         if( ( objsense == 1 && (db-sol[prob] > reltol || sol[prob]-pb > reltol) ) || ( objsense == -1 && (sol[prob]-db > reltol || pb-sol[prob] > reltol) ) ) {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else {
            if( timeout || gapreached || sollimitreached || memlimitreached || nodelimitreached ) 
	    {
               if( timeout )
                  status = "timeout";
               else if( gapreached )
                  status = "gaplimit";
               else if( sollimitreached )
                  status = "sollimit";
               else if( memlimitreached )
                  status = "memlimit";
               else if( nodelimitreached )
                  status = "nodelimit";

               timeouttime += tottime;
               timeouts++;
            }
            else {
               if( (abs(pb - db) <= max(abstol, reltol)) && abs(pb - sol[prob]) <= reltol ) {
                  status = "ok";
                  pass++;
               }
               else {
                  status = "fail";
                  failtime += tottime;
                  fail++;
               }
            }
         }
      }
      else if( solstatus[prob] == "best" ) {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

	 # objsense = 1 -> minimize; objsense = -1 -> maximize
         if( ( objsense == 1 && db-sol[prob] > reltol) || ( objsense == -1 && sol[prob]-db > reltol) ) {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else {
            if( timeout || gapreached || sollimitreached || memlimitreached || nodelimitreached ) {
               if( (objsense == 1 && sol[prob]-pb > reltol) || (objsense == -1 && pb-sol[prob] > reltol) ) {
                  status = "better";
                  timeouttime += tottime;
                  timeouts++;
               }
               else {
                  if( timeout )
                     status = "timeout";
                  else if( gapreached )
                     status = "gaplimit";
                  else if( sollimitreached )
                     status = "sollimit";
                  else if( memlimitreached )
                     status = "memlimit";
                  else if( nodelimitreached )
                     status = "nodelimit";
                  timeouttime += tottime;
                  timeouts++;
               }
            }
            else {
               if( abs(pb - db) <= max(abstol, reltol) ) {
                  if( abs(firstpb - rootdb) <= max(abstol,reltol) )
                     status = "solved not verified";
                  else
                     status = "solved";
                  pass++;
               }
               else {
                  status = "fail";
                  failtime += tottime;
                  fail++;
               }
            }
         }
      }
      else if( solstatus[prob] == "unkn" ) {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

	 if( timeout || gapreached || sollimitreached || memlimitreached || nodelimitreached ) {
            if( abs(pb) < infty ) {
               status = "better";
               timeouttime += tottime;
               timeouts++;
            }
	    else {
	       if( timeout )
		  status = "timeout";
	       else if( gapreached )
		  status = "gaplimit";
	       else if( sollimitreached )
		  status = "sollimit";
	       else if( memlimitreached )
		  status = "memlimit";
	       else if( nodelimitreached )
		  status = "nodelimit";
	       timeouttime += tottime;
	       timeouts++;
	    }
	 }
         else if( abs(pb - db) <= max(abstol, reltol) ) {
            status = "solved not verified";
            pass++;
         }
         else {
	    status = "unknown";
         }
      }
      else if( solstatus[prob] == "inf" ) {
         if( !feasible ) {
	    if( timeout || memlimitreached || nodelimitreached ) {
	       if( timeout )
		  status = "timeout";
	       else if( memlimitreached )
		  status = "memlimit";
	       else if( nodelimitreached )
		  status = "nodelimit";
	       timeouttime += tottime;
	       timeouts++;
	    }
            else {
               status = "ok";
               pass++;
            }
         }
         else {
            status = "fail";
            failtime += tottime;
            fail++;
         }
      }
      else {
         reltol = 1e-5 * max(abs(pb),1.0);
         abstol = 1e-4;

         if( timeout || gapreached || sollimitreached || memlimit || nodelimit ) {
	    if( timeout )
	       status = "timeout";
	    else if( gapreached )
	       status = "gaplimit";
	    else if( sollimitreached )
	       status = "sollimit";
	    else if( memlimitreached )
	       status = "memlimit";
	    else if( nodelimitreached )
	       status = "nodelimit";
	    timeouttime += tottime;
	    timeouts++;
	 }
         else if( abs(pb - db) < max(abstol,reltol) ) {
            status = "solved not verified";
            pass++;
         }
         else {
               status = "unknown";
         }
      }

      if( valgrinderror > 0 ) {
         status = "fail"
         failtime += tottime;
         fail++;
      } else if( valgrindleaks > 0 ) {
         status = "fail"
         failtime += tottime;
         fail++;
      }

      if( writesolufile ) {
         if( pb == +infty && db == +infty )
            printf("=inf= %-18s\n",prob)>NEWSOLUFILE;
         else if( pb == db )
            printf("=opt= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else if( pb < +infty )
            printf("=best= %-18s %16.9g\n",prob,pb)>NEWSOLUFILE;
         else
            printf("=unkn= %-18s\n",prob)>NEWSOLUFILE;
      }

      #write output to both the tex file and the console depending on whether printsoltimes is activated or not
      if( !onlypresolvereductions || origcons > cons || origvars > vars ) {
         if (TEXFILE != "") {
            printf("%-16s & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d &  %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d",
                   pprob, cgmipbase, cgmipmaster, cliquebase, cliquemaster, cmirbase, cmirmaster, flowcoverbase, flowcovermaster, gomorybase, gomorymaster, impliedboundsbase, impliedboundsmaster, mcfbase, mcfmaster, oddcyclebase, oddcyclemaster, strongcgbase, strongcgmaster, zerohalfbase, zerohalfmaster ) >TEXFILE;
            if( printsoltimes )
               printf(" & %7.1f & %7.1f", timetofirst, timetobest) > TEXFILE;
            printf("\\\\\n") > TEXFILE;
         }
            printf("%-12s  %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d",
                   shortprob, cgmipbase, cgmipmaster, cliquebase, cliquemaster, cmirbase, cmirmaster, flowcoverbase, flowcovermaster, gomorybase, gomorymaster, impliedboundsbase, impliedboundsmaster, mcfbase, mcfmaster, oddcyclebase, oddcyclemaster, strongcgbase, strongcgmaster, zerohalfbase, zerohalfmaster );
         printf("\n");
      }
      if( PAVFILE != "" ) {
         #PAVER output: see http://www.gamsworld.org/performance/paver/pprocess_submit.htm
         if( status == "abort" ) {
            modelstat = 13;
            solverstat = 13;
         } else if( status == "fail" || status == "unknown" ) {
            modelstat = 7;
            solverstat = 1;
         } else if( status == "timeout" ) {
            modelstat = abs(pb) < infty ? 8 : 14;
            solverstat = 3;
         } else if( status == "gaplimit" || status == "better" ) {
            modelstat = 8;
            solverstat = 1;
         } else if( status == "ok" || status == "solved not verified" ) {
            modelstat = 1;
            solverstat = 1;
         } else {
            modelstat = 13;
            solverstat = 13;
         }
         pavprob = prob;
         if( length(pavprob) > 25 )
              pavprob = substr(pavprob, length(pavprob)-24,25);
         if( vars == 0 )
            gamsprobtype = "LP";
         else if( lincons < cons && binvars == 0 && intvars == 0 )
            gamsprobtype = "NLP";
         else if( lincons < cons )
            gamsprobtype = "MINLP";
         else if( binvars == 0 && intvars == 0 )
            gamsprobtype = "LP";
         else
            gamsprobtype = "MIP";
         #InputFileName,ModelType,SolverName,Direction,ModelStatus,SolverStatus,ObjectiveValue,ObjectiveValueEstimate,SolverTime
         #NumberOfNodes,NumberOfIterations,NumberOfEquations,NumberOfVariables
         printf("%s,%s,SCIP_%s,%d,%d,%d,%g,%g,%g,", pavprob, gamsprobtype, settings, objsense == 1 ? 1 : 0, modelstat, solverstat, pb, db, tottime+pavshift) > PAVFILE;
         printf("%d,%d,%d,%d\n", bbnodes, simpiters, cons, vars) > PAVFILE;
      }
   }
}
END {
   shiftednodegeom -= nodegeomshift;
   shiftedsblpgeom -= sblpgeomshift;
   shiftedtimegeom -= timegeomshift;
   shiftedconftimegeom -= timegeomshift;
   shiftedoverheadtimegeom -= timegeomshift;
   shiftedbasictimegeom -= timegeomshift;
   shiftedtimetofirstgeom -= timegeomshift;
   shiftedtimetobestgeom -= timegeomshift;
   shiftedpricecallsgeom -= pricecallsgeomshift;
   shiftedpriceprobsgeom -= priceprobsgeomshift;
   shiftedpricevarsgeom -= pricevarsgeomshift;

   shiftedbasecallsgeom -= basecallsgeomshift;
   shiftedbasefoundgeom -= basefoundgeomshift;
   shiftedbaseappliedgeom -= baseappliedgeomshift;
   shiftedbaseconvexgeom -= baseconvexgeomshift;
   shiftedbasediffstartgeom -= basediffstartgeomshift;
   shiftedbasediffendgeom -= basediffendgeomshift;

   shiftedcgmipbasegeom -= cgmipbasegeomshift;
   shiftedcgmipmastergeom -= cgmipmastergeomshift;

   shiftedcliquebasegeom -= cliquebasegeomshift;
   shiftedcliquemastergeom -= cliquemastergeomshift;
   shiftedcmirbasegeom -= cmirbasegeomshift;
   shiftedcmirmastergeom -= cmirmastergeomshift;
   shiftedflowcoverbasegeom -= flowcoverbasegeomshift;
   shiftedflowcovermastergeom -= flowcovermastergeomshift;
   shiftedgomorybasegeom -= gomorybasegeomshift;
   shiftedgomorymastergeom -= gomorymastergeomshift;
   shiftedimpliedboundsbasegeom -= impliedboundsbasegeomshift;
   shiftedimpliedboundsmastergeom -= impliedboundsmastergeomshift;
   shiftedmcfbasegeom -= mcfbasegeomshift;
   shiftedmcfmastergeom -= mcfmastergeomshift;
   shiftedoddcyclebasegeom -= oddcyclebasegeomshift;
   shiftedoddcyclemastergeom -= oddcyclemastergeomshift;
   shiftedstrongcgbasegeom -= strongcgbasegeomshift;
   shiftedstrongcgmastergeom -= strongcgmastergeomshift;
   shiftedzerohalfbasegeom -= zerohalfbasegeomshift;
   shiftedzerohalfmastergeom -= zerohalfmastergeomshift;



   shiftedmastercallsgeom -= mastercallsgeomshift;
   shiftedmasterfoundgeom -= masterfoundgeomshift;
   shiftedmasterappliedgeom -= masterappliedgeomshift;

   shiftedrootdbgeom -= rootdbgeomshift;

   if (TEXFILE != "" ) {
      printf("\\midrule\n")                                                 >TEXFILE;
      printf("%-11s (%2d) & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d &  %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d & %4d & %-4d",
                   "Total", nprobs, scgmipbase, scgmipmaster, scliquebase, scliquemaster, scmirbase, scmirmaster, sflowcoverbase, sflowcovermaster, sgomorybase, sgomorymaster, simpliedboundsbase, simpliedboundsmaster, smcfbase, smcfmaster, soddcyclebase, soddcyclemaster, sstrongcgbase, sstrongcgmaster, szerohalfbase, szerohalfmaster ) >TEXFILE;

      printf("\\\\\n") > TEXFILE;
      printf("%-11s      & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f %4.1f",
                   "Geom. Mean", cgmipbasegeom, cgmipmastergeom, cliquebasegeom, cliquemastergeom, cmirbasegeom, cmirmastergeom, flowcoverbasegeom, flowcovermastergeom, gomorybasegeom, gomorymastergeom, impliedboundsbasegeom, impliedboundsmastergeom, mcfbasegeom, mcfmastergeom, oddcyclebasegeom, oddcyclemastergeom, strongcgbasegeom, strongcgmastergeom, zerohalfbasegeom, zerohalfmastergeom ) >TEXFILE;

      printf("\\\\\n") > TEXFILE;
      printf("%-13s    & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f & %4.1f %4.1f",
                   "Shifted Geom.", shiftedcgmipbasegeom, shiftedcgmipmastergeom, shiftedcliquebasegeom, shiftedcliquemastergeom, shiftedcmirbasegeom, shiftedcmirmastergeom, shiftedflowcoverbasegeom, shiftedflowcovermastergeom, shiftedgomorybasegeom, shiftedgomorymastergeom, shiftedimpliedboundsbasegeom, shiftedimpliedboundsmastergeom, shiftedmcfbasegeom, shiftedmcfmastergeom, shiftedoddcyclebasegeom, shiftedoddcyclemastergeom, shiftedstrongcgbasegeom, shiftedstrongcgmastergeom, shiftedzerohalfbasegeom, shiftedzerohalfmastergeom ) >TEXFILE;

      printf("\\\\\n") > TEXFILE;
   }
   printf(tablehead3);
   printf("\n");

   tablebottom1 = "-------------------------------------------------------------------------------------------------------------------";
   tablebottom2 = "  ";
   tablebottom3 = "---------------------------------------------------------------------------------------------------------------------------";
   
   
   tablebottom1 = tablebottom1"\n";
   tablebottom2 = tablebottom2"\n";
   tablebottom3 = tablebottom3"\n";
   
   printf(tablebottom1);
   printf(tablebottom2);
   printf(tablebottom3);

      printf("Total (%5d) %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d | %4d %4d",
                   nprobs, scgmipbase, scgmipmaster, scliquebase, scliquemaster, scmirbase, scmirmaster, sflowcoverbase, sflowcovermaster, sgomorybase, sgomorymaster, simpliedboundsbase, simpliedboundsmaster, smcfbase, smcfmaster, soddcyclebase, soddcyclemaster, sstrongcgbase, sstrongcgmaster, szerohalfbase, szerohalfmaster );

   
   printf("\n");
      printf("shifted geom. %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f | %4.1f %4.1f",
                   shiftedcgmipbasegeom, shiftedcgmipmastergeom, shiftedcliquebasegeom, shiftedcliquemastergeom, shiftedcmirbasegeom, shiftedcmirmastergeom, shiftedflowcoverbasegeom, shiftedflowcovermastergeom, shiftedgomorybasegeom, shiftedgomorymastergeom, shiftedimpliedboundsbasegeom, shiftedimpliedboundsmastergeom, shiftedmcfbasegeom, shiftedmcfmastergeom, shiftedoddcyclebasegeom, shiftedoddcyclemastergeom, shiftedstrongcgbasegeom, shiftedstrongcgmastergeom, shiftedzerohalfbasegeom, shiftedzerohalfmastergeom );

   printf("\n");
   printf(tablebottom3);

   if (TEXFILE != "" ) {
      printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
      printf("\\end{supertabular}\n")                                      >TEXFILE;
      printf("\\end{center}\n")                                             >TEXFILE;
      printf("\\end{document}\n")                                           >TEXFILE;
   }

   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 GCG(%s)SCIP(%s)%s(%s):%s\n", gcgversion, scipversion, lpsname, lpsversion, settings);
}
