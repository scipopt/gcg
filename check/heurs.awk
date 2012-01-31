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
#
#@file    heurs.awk
#@brief   Check Report Generator for GCG Heuristics
#@author  Christian Puchert

function max(x,y)
{
   return (x) > (y) ? (x) : (y);
}

function floor(x) {
   y = int(x)
   return (y) > (x) ? (y) - 1 : (y);
}
 
function ceil(x) {
   y = int(x)
   return (y) < (x) ? (y) + 1 : (y);
}

BEGIN {
   timegeomshift = 10.0;
   solsgeomshift = 10.0;
   onlyinsolufile = 0;  # should only instances be reported that are included in the .solu file?
   onlyintestfile = 0;  # should only instances be reported that are included in the .test file?  TEMPORARY HACK!
   useshortnames = 1;   # should problem name be truncated to fit into column?
   headerprinted = 0;

   firstheur = "";
   firstmasterheur = "";
   bestheur = "";
   bestmasterheur = "";
   
   nprobs = 0;
   nmheurs = 0;
   noheurs = 0;
   
   nmheurs = split(masterheurs, mheurs, ",");
   noheurs = split(origheurs, oheurs, ",");
   for( i = 1; i <= nmheurs; ++i )
      heurs[i] = mheurs[i];
   for( i = 1; i <= noheurs; ++i )
      heurs[nmheurs + i] = oheurs[i];
   nheurs = nmheurs + noheurs;
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
   for( i = 2; i <= n; ++i )
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
   while( getline > 0 && match($0, /^LP                 :/) == 0 )
   {
      heur = $1;
      sub(/:/, "", heur);
      
      if( heur != "LP" && heur != "pseudo" && !ignored[heur] )
      {

         # only read statistics if the heuristic was in the input
         if( inoriginalprob )
         {
            for( i = 1; i <= noheurs; ++i )
            {
               if( oheurs[i] == heur )
                  break;
            }
            if( i > noheurs )
               continue;
            i = nmheurs + i;
         }
         else
         {
            for( i = 1; i <= nmheurs; ++i )
               if( mheurs[i] == heur )
                  break;
            if( i > nmheurs )
               continue;
         }

         if( headerprinted == 0 )
         {
            stime[i] = 0.0;
            scalls[i] = 0;
            sfound[i] = 0;
            
            timegeom[i] = 0.0;
            callsgeom[i] = 0.0;
            foundgeom[i] = 0.0;
            
            shiftedtimegeom[i] = timegeomshift;
            shiftedcallsgeom[i] = solsgeomshift;
            shiftedfoundgeom[i] = solsgeomshift;
         }
               
         time[i] = $(NF-2);
         calls[i] = $(NF-1);
         found[i] = $NF;
      }
   }
}

/^  First Solution   :/ {
   if( inoriginalprob )
   {
      firstheur = $NF;
      gsub(/[<>)]/, "", firstheur);
   }
   else
   {
      firstmasterheur = $NF;
      gsub(/[<>)]/, "", firstmasterheur);
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
   else
   {      
      if( $4 != "infeasible" && $4 != "-" )
      {
         bestmasterheur = $NF;
         gsub(/[<>)]/, "", bestmasterheur);
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
         if( nheurs >= 1 )
         {
            tablehead01 = "                  ";
            tablehead02 = "                  ";
            if( nmheurs >= 1 )
            {
               for( i = 1; i <= nmheurs; ++i )
               {
                  tablehead01 = tablehead01"+-----------------";
               }
               tablehead02 = tablehead02"|";
               for( i = 1; i <= floor((18 * nmheurs - 15) / 2); ++i )
                  tablehead02 = tablehead02" ";
               tablehead02 = tablehead02"Master problem";
               for( i = 1; i <= ceil((18 * nmheurs - 15) / 2); ++i )
                  tablehead02 = tablehead02" ";
            }
            if( noheurs >= 1 )
            {
               for( i = 1; i <= noheurs; ++i )
               {
                  tablehead01 = tablehead01"+-----------------";
               }
               tablehead02 = tablehead02"|";
               for( i = 1; i <= floor((18 * noheurs - 17) / 2); ++i )
                  tablehead02 = tablehead02" ";
               tablehead02 = tablehead02"Original problem";
               for( i = 1; i <= ceil((18 * noheurs - 17) / 2); ++i )
                  tablehead02 = tablehead02" ";
            }
            tablehead01 = tablehead01"+\n";
            tablehead02 = tablehead02"|\n";

            printf(tablehead01);
            printf(tablehead02);

         }

         tablehead1 = "------------------+";
         tablehead2 = "Name              |";
         tablehead3 = "------------------+";

         for( i = 1; i <= nheurs; ++i )
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
         tablehead1 = tablehead1"-------------------+-------------------\n";
         tablehead2 = tablehead2"First solution     |Best solution      \n";
         tablehead3 = tablehead3"-------------------+-------------------\n";
            
         printf(tablehead1);
         printf(tablehead2);
         printf(tablehead3);
            
         headerprinted = 1;         
      }
   
      nprobs++;
      
      alltime = 0;
      allcalls = 0;
      allfound = 0;
      
      if( !onlypresolvereductions || origcons > cons || origvars > vars )
      {
         printf("%-18s ", shortprob);

         for( i = 1; i <= nheurs; ++i )
         {
            stime[i] += time[i];
            scalls[i] += calls[i];
            sfound[i] += found[i];
         
            timegeom[i] = timegeom[i]^((nprobs-1)/nprobs) * max(time[i], 1.0)^(1.0/nprobs);
            callsgeom[i] = callsgeom[i]^((nprobs-1)/nprobs) * max(calls[i], 1.0)^(1.0/nprobs);
            foundgeom[i] = foundgeom[i]^((nprobs-1)/nprobs) * max(found[i], 1.0)^(1.0/nprobs);
         
            shiftedtimegeom[i] = shiftedtimegeom[i]^((nprobs-1)/nprobs) * max(time[i]+timegeomshift, 1.0)^(1.0/nprobs);
            shiftedcallsgeom[i] = shiftedcallsgeom[i]^((nprobs-1)/nprobs) * max(calls[i]+solsgeomshift, 1.0)^(1.0/nprobs);
            shiftedfoundgeom[i] = shiftedfoundgeom[i]^((nprobs-1)/nprobs) * max(found[i]+solsgeomshift, 1.0)^(1.0/nprobs);
         
            alltime += time[i];
            allcalls += calls[i];
            allfound += found[i];
      
            printf("%7.1f/%4d/%4d ", time[i], calls[i], found[i]);
         }

         if( nheurs > 1 )
            printf("%9.1f/%5d/%5d ", alltime, allcalls, allfound);

         if( firstheur == "relaxation" && firstmasterheur != "" )
            firstheur = sprintf("m:%s", firstmasterheur);
         if( bestheur == "relaxation" && bestmasterheur != "" )
            bestheur = sprintf("m:%s", bestmasterheur);
         printf("%-19s %-19s\n", firstheur, bestheur);
      }
   }

   firstheur = "";
   firstmasterheur = "";
   bestheur = "";
   bestmasterheur = "";

}

END {
   for( i = 1; i <= nheurs; ++i )
   {
      shiftedtimegeom[i] -= timegeomshift;
      shiftedcallsgeom[i] -= solsgeomshift;
      shiftedfoundgeom[i] -= solsgeomshift;
   }

   alltime = 0;
   allcalls = 0;
   allfound = 0;

   printf("------------------+");
   for( i = 1; i <= nheurs; ++i )
      printf("-----------------+");
   if( nheurs > 1)
      printf("---------------------+");
   printf("-------------------+-------------------\n");
   
   printf("Total (%4d)       ", nprobs);
   for( i = 1; i <= nheurs; ++i )
   {
      printf("%7.1f/%4d/%4d ", stime[i], scalls[i], sfound[i]);
      
      alltime += stime[i];
      allcalls += scalls[i];
      allfound += sfound[i];

   }

   if( nheurs > 1 )
      printf("%9.1f/%5d/%5d\n", alltime, allcalls, allfound);
   else
      printf("\n");
   
   printf("Geom. Mean         ");
   for( i = 1; i <= nheurs; ++i )
   {
      printf("%7.1f/%4d/%4d ", timegeom[i], callsgeom[i], foundgeom[i]);
   }
   printf("                   \n");   
   
   printf("Shifted Mean       ");
   for( i = 1; i <= nheurs; ++i )
   {
      printf("%7.1f/%4d/%4d ", shiftedtimegeom[i], shiftedcallsgeom[i], shiftedfoundgeom[i]);
   }
   printf("                   \n");
   printf("\n");   
}
