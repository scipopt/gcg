/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/* This program is free software; you can redistribute it and/or             */
/* modify it under the terms of the GNU Lesser General Public License        */
/* as published by the Free Software Foundation; either version 3            */
/* of the License, or (at your option) any later version.                    */
/*                                                                           */
/* This program is distributed in the hope that it will be useful,           */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of            */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             */
/* GNU Lesser General Public License for more details.                       */
/*                                                                           */
/* You should have received a copy of the GNU Lesser General Public License  */
/* along with this program; if not, write to the Free Software               */
/* Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   class_miscvisualization.cpp
 * @brief  miscellaneous methods for visualizations
 * @author Hanna Franzen
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "cons_decomp.h"
#include "class_miscvisualization.h"
#include "class_seeedpool.h"
#include "class_seeed.h"
#include "wrapper_seeed.h"

#include "scip/scip.h"

#include <unistd.h>
#include <stdlib.h>
#include <sstream>

using namespace gcg;


/** constructor */
MiscVisualization::MiscVisualization(){}

/** destructor */
MiscVisualization::~MiscVisualization(){}

/** gives a consistent filename for a (single) seeed visualization that includes the probname and seeedID
 *
 * @return filename including the extension
 * */
SCIP_RETCODE MiscVisualization::GCGgetVisualizationFilename(
   SCIP* scip,             /**< scip data structure */
   SeeedPtr seeed,         /**< seeed that is to be visualized */
   const char* extension,  /**< file extension */
   char* filename          /**< filename output */
   )
{
   char* name;
   char* detectorchainstring;
   char probname[SCIP_MAXSTRLEN];

   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);

   /* get detector chain string*/
   detectorchainstring = seeed->getDetectorChainString();

   /* print header */
   if( seeed == NULL )
      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s", name);
   else
   {
      (void) SCIPsnprintf(filename, SCIP_MAXSTRLEN, "%s-%s-%d-%d-%s", name, detectorchainstring, seeed->getID(),
         seeed->getNBlocks(), extension);
   }

   return SCIP_OKAY;
}


/** gives the path of the file */
char* MiscVisualization::GCGgetFilePath(
   SCIP* scip,       /**< scip data structure */
   FILE* file        /**< file */
   )
{
   char* pfile = '\0';
   char sympath[SCIP_MAXSTRLEN];
   int filedesc;
   int success;

   filedesc = fileno(file); /* get link to file descriptor */
   if( filedesc < 0 )
   {
      SCIPerrorMessage("File reading error, no fileno!\n");
   }
   snprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
   pfile = realpath(sympath, NULL);
   //success = readlink(sympath, pfile, SCIP_MAXSTRLEN); /* get actual path including extension */
//   if( success < 0 )
//   {
//      SCIPerrorMessage("File reading error!\n");
//   }
   return pfile;
}


/** checks in which seeedpool the seeed with given ID is stored and returns that seeedpool
 *
 * @returns pool: Seeedpool* where the Seeed was found
 */
Seeedpool* MiscVisualization::GCGgetSeeedpoolForSeeed(
   SCIP* scip,       /**< SCIP data structure */
   int seeedid       /**< ID of Seeed */
   )
{
   SEEED_WRAPPER seeedwr;
   Seeedpool* seeedpool;

   GCGgetCurrentSeeedpools(scip, &seeedwr, NULL);
   seeedpool = seeedwr.seeedpool;

   /* find in presolved */
   for( int i = 0; i < seeedpool->getNAncestorSeeeds(); ++i)
   {
      if( seeedpool->getAncestorSeeed(i)!= NULL && seeedpool->getAncestorSeeed(i)->getID() == seeedid )
         return seeedpool;
   }

   for( int i = 0; i < seeedpool->getNIncompleteSeeeds(); ++i)
   {
      if( seeedpool->getIncompleteSeeed(i)->getID() == seeedid )
         return seeedpool;
   }

   for( int i = 0; i < seeedpool->getNFinishedSeeeds(); ++i)
   {
      if( seeedpool->getFinishedSeeed(i)->getID() == seeedid )
         return seeedpool;
   }

   for( int i = 0; i < seeedpool->getNCurrentSeeeds(); ++i)
   {
      if( seeedpool->getCurrentSeeed(i)->getID() == seeedid )
         return seeedpool;
   }

   /* find in unpresolved */
   GCGgetCurrentSeeedpools(scip, NULL, &seeedwr);
   seeedpool = seeedwr.seeedpool;

   if( seeedpool != NULL )
   {
      for( int i = 0; i < seeedpool->getNAncestorSeeeds(); ++i)
      {
         if( seeedpool->getAncestorSeeed(i)!= NULL && seeedpool->getAncestorSeeed(i)->getID() == seeedid )
            return seeedpool;
      }

      for( int i = 0; i < seeedpool->getNIncompleteSeeeds(); ++i)
      {
         if( seeedpool->getIncompleteSeeed(i)->getID() == seeedid )
            return seeedpool;
      }

      for( int i = 0; i < seeedpool->getNFinishedSeeeds(); ++i)
      {
         if( seeedpool->getFinishedSeeed(i)->getID() == seeedid )
            return seeedpool;
      }

      for( int i = 0; i < seeedpool->getNCurrentSeeeds(); ++i)
      {
         if( seeedpool->getCurrentSeeed(i)->getID() == seeedid )
            return seeedpool;
      }
   }

   return NULL;
}
