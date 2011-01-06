#!/bin/gawk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Colum Generation                                 *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id$
#
#@file    check.awk
#@brief   SCIP Check Report Generator
#@author  Thorsten Koch
#@author  Tobias Achterberg
#@author  Alexander Martin
#@author  Timo Berthold
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
   nodegeomshift = 10.0;
   pricegeomshift = 5.0;
   lpgeomshift = 2.0;
   pricecallsgeomshift = 100;
   priceprobsgeomshift = 100;
   pricevarsgeomshift = 100;
   sblpgeomshift = 0.0;
   pavshift = 1.0;
   onlyinsolufile = 0;  # should only instances be reported that are included in the .solu file?
   onlyintestfile = 0;  # should only instances be reported that are included in the .test file?  TEMPORARY HACK!
   onlypresolvereductions = 0;  # should only instances with presolve reductions be shown?
   useshortnames = 1;   # should problem name be truncated to fit into column?

   printf("\\documentclass[leqno]{article}\n")                      >TEXFILE;
   printf("\\usepackage{a4wide}\n")                                 >TEXFILE;
   printf("\\usepackage{amsmath,amsfonts,amssymb,booktabs}\n")      >TEXFILE;
   printf("\\pagestyle{empty}\n\n")                                 >TEXFILE;
   printf("\\begin{document}\n\n")                                  >TEXFILE;
   printf("\\begin{table}[p]\n")                                    >TEXFILE;
   printf("\\begin{center}\n")                                      >TEXFILE;
   printf("\\setlength{\\tabcolsep}{2pt}\n")                        >TEXFILE;
   printf("\\newcommand{\\g}{\\raisebox{0.25ex}{\\tiny $>$}}\n")    >TEXFILE;
   printf("\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lrrrrrrrrrrrrrrr@{}}\n") >TEXFILE;
   printf("\\toprule\n")                                         >TEXFILE;
   printf("Name             &  Conss &   Vars &       Dual Bound &     Primal Bound &  Gap\\%% &   Calls &   Probs &    Vars &  P-Time & LP-Time &      Nodes &     Time \\\\\n") > TEXFILE;
   printf("\\midrule\n")                                         >TEXFILE;
   
   printf("-----------+----+--- Original --+-- Presolved --+------+------+------+--------------+--------------+------+------------ Pricing ----------+-------+--------+-------+-------+-------\n");
   printf("Name       |Type| Conss |  Vars | Conss |  Vars |Blocks| Rel  |MConss|  Dual Bound  | Primal Bound | Gap%% | Calls | Probs |  Vars |  Time |LP-Time|  Iters | Nodes |  Time |       \n");
   printf("-----------+----+-------+-------+-------+-------+------+------+------+--------------+--------------+------+-------+-------+-------+-------+-------+--------+-------+-------+-------\n");

   nprobs = 0;
   sbab = 0;
   slp = 0;
   ssim = 0;
   ssblp = 0;
   stottime = 0.0;
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
   timeouttime = 0.0;
   timeouts = 0;
   failtime = 0.0;
   fail = 0;
   pass = 0;
   settings = "default";
   conftottime = 0.0;
   overheadtottime = 0.0;
   timelimit = 0.0;
}
/^IP\// { # TEMPORARY HACK to parse .test files
   intestfile[$1] = 1;
}
/=opt=/  {  # get optimum
   if( NF >= 3 ) {
      solstatus[$2] = "opt";
      sol[$2] = $3;
   }
   else
      solstatus[$2] = "feas";
}
/=inf=/  { solstatus[$2] = "inf"; sol[$2] = 0.0; } # problem infeasible
/=best=/ { solstatus[$2] = "best"; sol[$2] = $3; } # get best known solution value
/=unkn=/ { solstatus[$2] = "unkn"; }               # no feasible solution known
/^@01/ { 
   filename = $2;

   n  = split ($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > 18 )
      shortprob = substr(prob, length(prob)-17, 18);
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
   origvars = 0;
   origcons = 0;
   timeout = 0;
   feasible = 1;
   pb = +1e20;
   db = -1e20;
   simpiters = 0;
   bbnodes = 0;
   primlps = 0;
   primiter = 0;
   duallps = 0;
   dualiter = 0;
   sblps = 0;
   sbiter = 0;
   tottime = 0.0;
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
   timelimit = 0.0;
   starttime = 0.0;
   endtime = 0.0;
   inoriginalprob = 1;
   inmasterprob = 1;
}
/@03/ { starttime = $2; }
/@04/ { endtime = $2; }
/^SCIP> SCIP> / { $0 = substr($0, 13, length($0)-12); }
/^SCIP> / { $0 = substr($0, 7, length($0)-6); }
/^loaded parameter file/ { settings = $4; sub(/<.*settings\//, "", settings); sub(/\.set>/, "", settings); }
/^parameter <limits\/time> set to/ { timelimit = $5; }
#
# problem: master or original?
#
/^Original Program statistics:/ { inmasterprob = 0; inoriginalprob = 1; }
#
# conflict analysis
#
/^Conflict Analysis  :/ { inconflict = 1; }
/^  propagation      :/ {
   if( inconflict == 1 )
   {
      conftime += $3; #confclauses += $5 + $7; confliterals += $5 * $6 + $7 * $8;
   }
}
/^  infeasible LP    :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  strong branching :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  pseudo solution  :/ {
   if( inconflict == 1 )
   {
      conftime += $4; #confclauses += $6 + $8; confliterals += $6 * $7 + $8 * $9;
   }
}
/^  applied globally :/ {
   if( inconflict == 1 )
   {
      confclauses += $7; confliterals += $7 * $8;
   }
}
/^  applied locally  :/ {
   if( inconflict == 1 )
   {
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
      pricecall = $4;
      pricevars = $5;
   }
}

/^Separators         :/ { inconflict = 0; }
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
      if( inoriginalprob )
         mcons = $3;
   }
}
/^Matrix has / {
   blocks = $3;
   rel = $5;
}
/^  mip              :/ {
   npriceprobs = $4 + $6;
}
#
# solution
#
/^Original Problem   : no problem exists./ { readerror = 1; }
/^SCIP Status        :/ { aborted = 0; }
/solving was interrupted/  { timeout = 1; }
/gap limit reached/ { gapreached = 1; }
/solution limit reached/ { sollimitreached = 1; }
/memory limit reached/ { memlimitreached = 1; }
/problem is solved/    { timeout = 0; }
/^  Primal Bound     :/ {
   if( $4 == "infeasible" )
   {
      pb = 1e+20;
      db = 1e+20;
      feasible = 0;
   }
   else if( $4 == "-" )
   {
      pb = 1e+20;
      feasible = 0;
   }
   else
      pb = $4;
}
/^  Dual Bound       :/ { 
   if( $4 != "-" ) 
      db = $4;
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
/^Solving Time       :/ { tottime = $4 }
#
# Output
#
/^=ready=/ {
   if( (!onlyinsolufile || solstatus[prob] != "") &&
      (!onlyintestfile || intestfile[filename]) )
   {
      nprobs++;

      optimal = 0;
      markersym = "\\g";
      if( abs(pb - db) < 1e-06 )
      {
         gap = 0.0;
         optimal = 1;
         markersym = "  ";
      }
      else if( abs(db) < 1e-06 )
         gap = -1.0;
      else if( pb*db < 0.0 )
         gap = -1.0;
      else if( abs(db) >= 1e+20 )
         gap = -1.0;
      else if( abs(pb) >= 1e+20 )
         gap = -1.0;
      else
         gap = 100.0*abs((pb-db)/db);
      if( gap < 0.0 )
         gapstr = "  --  ";
      else if( gap < 1e+04 )
         gapstr = sprintf("%6.1f", gap);
      else
         gapstr = " Large";

      if( vars == 0 )
         probtype = "--";
      else if( binvars == 0 && intvars == 0 )
         probtype = "LP";
      else if( contvars == 0 )
      {
         if( intvars == 0 && implvars == 0 )
            probtype = "BP";
         else
            probtype = "IP";
      }
      else
      {
         if( intvars == 0 )
            probtype = "MBP";
         else
            probtype = "MIP";
      }

      if( aborted && endtime - starttime > timelimit && timelimit > 0.0 )
      {
         timeout = 1;
         aborted = 0;
         tottime = endtime - starttime;
      }
      else if( gapreached || sollimitreached )
         timeout = 0;
      if( aborted && tottime == 0.0 )
         tottime = timelimit;
      if( timelimit > 0.0 )
         tottime = min(tottime, timelimit);

      lps = primlps + duallps;
      simplex = primiter + dualiter;
      stottime += tottime;
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

      status = "";
      if( readerror )
      {
         status = "readerror";
         failtime += tottime;
         fail++;
      }
      else if( aborted )
      {
         status = "abort";
         failtime += tottime;
         fail++;
      }
      else if( solstatus[prob] == "opt" )
      {
         if( !sollimitreached && !gapreached && !timeout )
            wronganswer = (pb - db > 1e-4 || abs(pb - sol[prob]) > 1e-5*max(abs(pb),1.0));
         else if( pb >= db )
            wronganswer = (db > sol[prob] + 1e-5*max(abs(sol[prob]),1.0) ||
               pb < sol[prob] - 1e-5*max(abs(sol[prob]),1.0));
         else
            wronganswer = (pb > sol[prob] + 1e-5*max(abs(sol[prob]),1.0) ||
               db < sol[prob] - 1e-5*max(abs(sol[prob]),1.0));

         if( wronganswer )
         {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else if( timeout )
         {
            status = "timeout";
            timeouttime += tottime;
            timeouts++;
         }
         else
         {
            status = "ok";
            pass++;
         }
      }
      else if( solstatus[prob] == "feas" || solstatus[prob] == "inf" )
      {
         if( timeout )
            wronganswer = (feasible && solstatus[prob] == "inf");
         else
            wronganswer = (feasible != (solstatus[prob] == "feas"));

         if( wronganswer )
         {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else if( timeout )
         {
            status = "timeout";
            timeouttime += tottime;
            timeouts++;
         }
         else
         {
            status = "ok";
            pass++;
         }
      }
      else if( solstatus[prob] == "best" )
      {
         if( db > sol[prob] + 1e-4 )
         {
            status = "fail";
            failtime += tottime;
            fail++;
         }
         else if( timeout )
         {
            status = "timeout";
            timeouttime += tottime;
            timeouts++;
         }
         else
            status = "unknown";
      }
      else
      {
         if( timeout )
         {
            status = "timeout";
            timeouttime += tottime;
            timeouts++;
         }
         else
            status = "unknown";
      }

      if( !onlypresolvereductions || origcons > cons || origvars > vars )
      {
         printf("%-16s & %6d & %6d & %16.9g & %16.9g & %6s &  %6d &  %6d &  %6d & %7.1f & %7.1f & %s%8d &%s%7.1f \\\\\n",
                pprob, cons, vars, db, pb, gapstr, pricecall, npriceprobs, pricevars, pricetime, lptime,
		markersym, bbnodes, markersym, tottime) >TEXFILE;
         printf("%-12s %-3s %7d %7d %7d %7d %6d %6d %6d %14.9g %14.9g %6s %7d %7d %7d %7.1f %7.1f %8d %7d %7.1f %s\n",
		shortprob, probtype, origcons, origvars, cons, vars, blocks, rel, mcons, db, pb, gapstr, 
		pricecall, npriceprobs, pricevars, pricetime, lptime, simpiters, bbnodes, tottime, status);
      }

      # PAVER output: see http://www.gamsworld.org/performance/paver/pprocess_submit.htm
      if( solstatus[prob] == "opt" || solstatus[prob] == "feas" )
         modelstat = 1;
      else if( solstatus[prob] == "inf" )
         modelstat = 1;
      else if( solstatus[prob] == "best" )
         modelstat = 8;
      else
         modelstat = 1;
      if( status == "ok" || status == "unknown" )
         solverstat = 1;
      else if( status == "timeout" )
         solverstat = 3;
      else
         solverstat = 10;
      pavprob = prob;
      if( length(pavprob) > 25 )
         pavprob = substr(pavprob, length(pavprob)-24,25);
      printf("%s,MIP,SCIP_%s,0,%d,%d,%g,%g\n", pavprob, settings, modelstat, solverstat, pb, tottime+pavshift) > PAVFILE;
   }
}
END {
   shiftednodegeom -= nodegeomshift;
   shiftedsblpgeom -= sblpgeomshift;
   shiftedtimegeom -= timegeomshift;
   shiftedconftimegeom -= timegeomshift;
   shiftedoverheadtimegeom -= timegeomshift;
   shiftedbasictimegeom -= timegeomshift;


   printf("\\midrule\n")                                                 >TEXFILE;
   printf("%-11s (%2d) &        &        &                  &                  &        & %7d & %7d & %7d &%8.1f &%8.1f &  %9d & %8.1f \\\\\n",
          "Total", nprobs, spricecalls, spriceprobs, spricevars, spricetime, slptime, sbab, stottime) >TEXFILE;
   printf("%-11s      &        &        &                  &                  &        & %7d & %7d & %7d &%8.1f &%8.1f &  %9d & %8.1f \\\\\n",
          "Geom. Mean", pricecallsgeom, priceprobsgeom, pricevarsgeom, pricegeom, lpgeom, nodegeom, timegeom) >TEXFILE;
   printf("%-13s    &        &        &                  &                  &        & %7d & %7d & %7d &%8.1f &%8.1f &  %9d & %8.1f \\\\\n",
          "Shifted Geom.", shiftedpricecallsgeom, shiftedpriceprobsgeom, shiftedpricevarsgeom, shiftedpricegeom, 
	  shiftedlpgeom, shiftednodegeom, shiftedtimegeom) >TEXFILE;
   printf("-----------+----+-------+-------+-------+-------+------+------+------+--------------+--------------+------+-------+-------+-------+-------+-------+--------+-------+-------+-------\n");
   printf("\n");
   printf("------------------------------[Nodes]---------------[Time]----------[Pricing-Time]--------[LP-Time]----------[Pricing-Probs]----\n");
   printf("  Cnt  Pass  Time  Fail  total(k)     geom.     total     geom.     total     geom.     total     geom.    total     geom. \n");
   printf("---------------------------------------------------------------------------------------------------------------------------\n");
   printf("%5d %5d %5d %5d %9d %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9.1f %9d %9.1f\n",
	  nprobs, pass, timeouts, fail, sbab / 1000, nodegeom, stottime, timegeom, spricetime, pricegeom, 
	  slptime, lpgeom, spriceprobs, priceprobsgeom);
   printf(" shifted geom. [%4d/%4.1f/%3.1f/%3.1f]%9.1f           %9.1f           %9.1f           %9.1f           %9.1f \n",
	  nodegeomshift, timegeomshift, pricegeomshift, lpgeomshift, shiftednodegeom, shiftedtimegeom, 
	  shiftedpricegeom, shiftedlpgeom, shiftedpriceprobsgeom);
   printf("---------------------------------------------------------------------------------------------------------------------------\n");
   printf("\\bottomrule\n")                                              >TEXFILE;
   printf("\\noalign{\\vspace{6pt}}\n")                                  >TEXFILE;
   printf("\\end{tabular*}\n")                                           >TEXFILE;
   printf("\\caption{%s}\n", settings)                                   >TEXFILE;
   printf("\\end{center}\n")                                             >TEXFILE;
   printf("\\end{table}\n")                                              >TEXFILE;
   printf("\\end{document}\n")                                           >TEXFILE;

   printf("@02 timelimit: %g\n", timelimit);
   printf("@01 SCIP:%s\n", settings);
}
