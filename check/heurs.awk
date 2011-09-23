#!/bin/gawk -f
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
#*                                                                           *
#*                  This file is part of the program                         *
#*          GCG --- Generic Column Generation                                *
#*                  a Dantzig-Wolfe decomposition based extension            *
#*                  of the branch-cut-and-price framework                    *
#*         SCIP --- Solving Constraint Integer Programs                      *
#*                                                                           *
#* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
# $Id: $
#
#@file    heurs.awk
#@brief   Check Report Generator for GCG Heuristics
#@author  Christian Puchert

function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}

BEGIN {
   timegeomshift = 10.0;
   solsgeomshift = 10.0;
   onlyinsolufile = 0;  # should only instances be reported that are included in the .solu file?
   onlyintestfile = 0;  # should only instances be reported that are included in the .test file?  TEMPORARY HACK!
   useshortnames = 1;   # should problem name be truncated to fit into column?
   headerprinted = 0;
   
   nprobs = 0;
   nheurs = 0;
   
   # temporary hack; if the user passes an argument heuristics=... by the commandline, only those will be considered
   considerall = 1;
   nheurs = split(heuristics, heurs, ",");
   if( nheurs > 0 )
      considerall = 0;
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

# ignore heuristics that have been switched off in the settings file
/^heuristics\/[A-Za-z]+\/freq/ {
   n = split($0, a, "=");
   m = split(a[1], b, "/");
   heur = b[2];
   freq = a[2];
   
   if( freq == -1 )
      ignored[heur] = 1;
}

/^@01/ { 
   filename = $2;

   n  = split($2, a, "/");
   m = split(a[n], b, ".");
   prob = b[1];
   if( b[m] == "gz" || b[m] == "z" || b[m] == "GZ" || b[m] == "Z" )
      m--;
   for( i = 2; i < m; ++i )
      prob = prob "." b[i];

   if( useshortnames && length(prob) > 18 )
   {
      shortprob = substr(prob, length(prob)-17, 18);
   }
   else
      shortprob = prob;

   # Escape _ for TeX
   n = split(prob, a, "_");
   pprob = a[1];
   for( i = 2; i <= n; i++ )
      pprob = pprob "\\_" a[i];
      
   aborted = 1;
   readerror = 0;
   inoriginalprob = 1;
}

/^Original Problem   : no problem exists./ { readerror = 1; }
/^SCIP Status        :/ { aborted = 0; }

/^Master Program statistics:/ { inoriginalprob = 0; }
/^Original Program statistics:/ { inoriginalprob = 1; }

# when reached the heuristics section, collect statistics for each heuristic
/^Primal Heuristics  :/ {
   if( inoriginalprob )
   {
      while( getline > 0 && match($0, /^LP                 :/) == 0 )
      {
         heur = $1;
         sub(/:/, "", heur);
      
         if( heur != "LP" && heur != "pseudo" && !ignored[heur] )
         {
#            for( i = 1; i <= nheurs; i++ )
#               if( heurs[i] == heur )
#                  break;
#            if( i == nheurs+1 && considerall == 1 )
            if( considerall == 1 && headerprinted == 0 )
            {
               heurs[nheurs] = heur;
               nheurs++;
            
               stime[heur] = 0.0;
               scalls[heur] = 0;
               sfound[heur] = 0;
            
               timegeom[heur] = 0.0;
               callsgeom[heur] = 0.0;
               foundgeom[heur] = 0.0;
            
               shiftedtimegeom[heur] = timegeomshift;
               shiftedcallsgeom[heur] = solsgeomshift;
               shiftedfoundgeom[heur] = solsgeomshift;
            }
               
            time[heur] = $(NF-2);
            calls[heur] = $(NF-1);
            found[heur] = $NF;
         }
      }
   }
}

/^  First Solution   :/ {
   if( inoriginalprob )
   {
      firstheur = $NF;
      gsub(/[<>)]/, "", firstheur);
   }
}
/^  Primal Bound     :/ {
   if( inoriginalprob )
   {
      if( $4 == "infeasible" || $4 == "-" )
         feasible = 0;
      else
      {
         feasible = 1;
         bestheur = $NF;
         gsub(/[<>)]/, "", bestheur);
      }
   }
}

/^=ready=/ {
   if( (!onlyinsolufile || solstatus[prob] != "") &&
      (!onlyintestfile || intestfile[filename]) &&
      !aborted && !readerror  )
   {
     
      if( headerprinted == 0 )
      {
         tablehead1 = "------------------+";
         tablehead2 = "Name              |";
         tablehead3 = "------------------+";
         for( i = 1; i <= nheurs; i++ )
         {
            tablehead1 = tablehead1"-----------------+";
            tablehead2 = sprintf("%s%-17s|", tablehead2, heurs[i]);
            tablehead3 = tablehead3"-----------------+";
         }
         if( nheurs > 1 )
         {
            tablehead1 = tablehead1"---------------------+";
            tablehead2 = tablehead2"All heuristics       |";
            tablehead3 = tablehead3"---------------------+";
         }
         tablehead1 = tablehead1"-----------------+-----------------\n";
         tablehead2 = tablehead2"First solution   |Best solution    \n";
         tablehead3 = tablehead3"-----------------+-----------------\n";
            
         printf(tablehead1);
         printf(tablehead2);
         printf(tablehead3);
            
         headerprinted = 1;         
      }
   
      nprobs++;
      
      alltime = 0;
      allcalls = 0;
      allfound = 0;
      
      for( i = 1; i <= nheurs; i++ )
      {
         heur = heurs[i];
         
         stime[heur] += time[heur];
         scalls[heur] += calls[heur];
         sfound[heur] += found[heur];
         
         timegeom[heur] = timegeom[heur]^((nprobs-1)/nprobs) * max(time[heur], 1.0)^(1.0/nprobs);
         callsgeom[heur] = callsgeom[heur]^((nprobs-1)/nprobs) * max(calls[heur], 1.0)^(1.0/nprobs);
         foundgeom[heur] = foundgeom[heur]^((nprobs-1)/nprobs) * max(found[heur], 1.0)^(1.0/nprobs);
         
         shiftedtimegeom[heur] = shiftedtimegeom[heur]^((nprobs-1)/nprobs) * max(time[heur]+timegeomshift, 1.0)^(1.0/nprobs);
         shiftedcallsgeom[heur] = shiftedcallsgeom[heur]^((nprobs-1)/nprobs) * max(calls[heur]+solsgeomshift, 1.0)^(1.0/nprobs);
         shiftedfoundgeom[heur] = shiftedfoundgeom[heur]^((nprobs-1)/nprobs) * max(found[heur]+solsgeomshift, 1.0)^(1.0/nprobs);
         
         alltime += time[heur];
         allcalls += calls[heur];
         allfound += found[heur];
      }
      
      if( !onlypresolvereductions || origcons > cons || origvars > vars )
      {
         printf("%-18s ", shortprob);
         for( i = 1; i <= nheurs; i++ )
         {
            heur = heurs[i];
            printf("%7.1f/%4d/%4d ", time[heur], calls[heur], found[heur]);
         }
         if( nheurs > 1 )
            printf("%9.1f/%5d/%5d ", alltime, allcalls, allfound);
         printf("%-17s %-17s\n", firstheur, bestheur);
      }
   }
}

END {
   for( i = 1; i <= nheurs; i++ )
   {
      heur = heurs[i];
   
      shiftedtimegeom[heur] -= timegeomshift;
      shiftedcallsgeom[heur] -= solsgeomshift;
      shiftedfoundgeom[heur] -= solsgeomshift;
   }

   alltime = 0;
   allcalls = 0;
   allfound = 0;

   printf("------------------+");
   for( i = 1; i <= nheurs; i++ )
   {
      printf("-----------------+");
   }
   if( nheurs > 1)
      printf("---------------------+");
   printf("-----------------+-----------------\n");
   
   printf("Total (%4d)       ", nprobs);
   for( i = 1; i <= nheurs; i++ )
   {
      heur = heurs[i];      
      printf("%7.1f/%4d/%4d ", stime[heur], scalls[heur], sfound[heur]);
      
      alltime += stime[heur];
      allcalls += scalls[heur];
      allfound += sfound[heur];

   }
   if( nheurs > 1 )
      printf("%9.1f/%5d/%5d\n", alltime, allcalls, allfound);
   else
      printf("\n");
   
   printf("Geom. Mean         ");
   for( i = 1; i <= nheurs; i++ )
   {
      heur = heurs[i];
      printf("%7.1f/%4d/%4d ", timegeom[heur], callsgeom[heur], foundgeom[heur]);
   }
   printf("                   \n");   
   
   printf("Shifted Mean       ");
   for( i = 1; i <= nheurs; i++ )
   {
      heur = heurs[i];
      printf("%7.1f/%4d/%4d ", shiftedtimegeom[heur], shiftedcallsgeom[heur], shiftedfoundgeom[heur]);
   }
   printf("                   \n");
   printf("\n");   
}
