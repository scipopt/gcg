/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2026 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   miscvisualization.cpp
 * @brief  miscellaneous methods for visualizations
 * @author Hanna Franzen
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#include "gcg/cons_decomp.h"
#include "gcg/miscvisualization.h"
#include "gcg/class_detprobdata.h"
#include "gcg/class_partialdecomp.h"

#include "scip/scip.h"

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <unistd.h>
#endif

#include <stdlib.h>
#include <sstream>

using namespace gcg;


/* Gives a consistent filename for a (single) partialdec visualization that includes the probname and partialdecID
 *
 * @returns standardized filename
 */
void GCGgetVisualizationFilename(
   GCG* gcg,                        /* scip data structure */
   PARTIALDECOMP* partialdec,       /* partialdec that is to be visualized */
   const char* extension,  /* file extension (to be included in the name) */
   char* filename          /* filename output */
   )
{
   SCIP* scip;
   char* name;
   char detectorchainstring[SCIP_MAXSTRLEN];
   char probname[SCIP_MAXSTRLEN];

   scip = GCGgetOrigprob(gcg);
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);

   /* get detector chain string*/
   partialdec->buildDecChainString(detectorchainstring);

   assert( partialdec != NULL );

   /* print header */
   if(strlen(detectorchainstring) > 0)
   {
      /* if there is a PARTIALDECOMP that was detected in GCG */
      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s-%s-%d-%d%s", name, detectorchainstring, partialdec->getID(),
         partialdec->getNBlocks(), extension);
   }
   else
   {
      /* if there is a PARTIALDECOMP but it was not detected in GCG */
      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s-%d-%d%s", name, partialdec->getID(),
         partialdec->getNBlocks(), extension);
   }

   /* some filenames can still have dots in them (usually from prob name) which can cause confusion.
    * Does not replace characters in the file extension. */
   for(size_t i = 0; i < strlen(filename) - strlen(extension); i++)
   {
      if(filename[i] == '.')
         filename[i] = '-';

      if(filename[i] == '(')
         filename[i] = '-';

      if(filename[i] == ')')
        filename[i] = '-';
   }
}


/* Gives the path of the provided file */
void GCGgetFilePath(
   FILE* file,       /* file */
   char* path        /* buffer containing the path afterward, must be of length PATH_MAX! */
   )
{
   char sympath[SCIP_MAXSTRLEN];
   int filedesc;

   filedesc = fileno(file); /* get link to file descriptor */
   if( filedesc < 0 )
   {
      SCIPerrorMessage("File reading error, no fileno!\n");
      return;
   }
   SCIPsnprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
   realpath(sympath, path);
}
