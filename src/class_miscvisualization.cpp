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

#include "class_miscvisualization.h"
#include "class_seeedpool.h"
#include "class_seeed.h"

#include "scip/scip.h"

#include <unistd.h>

namespace gcg {

/** constructor */
MiscVisualization::MiscVisualization(){

}

/** destructor */
MiscVisualization::~MiscVisualization(){

}

/** gives a consistent filename for a (single) seeed visualization that includes the probname and seeedID
 *
 * @return filename including the extension
 * */
char* MiscVisualization::GCGgetVisualizationFilename(
   SCIP* scip,       /**< scip data structure */
   SeeedPtr seeed,   /**< seeed that is to be visualized */
   char* extension   /**< file extension */
   )
{
   char* name;
   char* detectorchainstring;
   char probname[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char* filename;

   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);

   /*@todo change this for seeeds! */
//   /* get detector chain string*/
//   detectorchainstring = DECdecompGetDetectorChainString(scip, decdecomp);
//
//   /* print header */
//   if( decdecomp == NULL )
//      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", name);
//   else
//   {
//      if(outputPDF)
//         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%s_%d_%d", name, detectorchainstring, DECdecompGetSeeedID(decdecomp),
//            decdecomp->nblocks);
//      else
//         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s-%s-%d-%d", name, detectorchainstring, DECdecompGetSeeedID(decdecomp),
//            decdecomp->nblocks);
//   }

   return outname;
}

/** gives the path of the file */
char* MiscVisualization::GCGgetFilePath(
   SCIP* scip,       /**< scip data structure */
   FILE* file        /**< file */
   )
{
   char* pfile;
   char sympath[SCIP_MAXSTRLEN];
   int filedesc;
   int success;

   filedesc = fileno(file); /* get link to file descriptor */
   if( filedesc < 0 )
   {
      /*@todo error or similar*/
   }
   snprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
   success = readlink(sympath, pfile, SCIP_MAXSTRLEN); /* get actual path including extension */
   if( success < 0 )
   {
      /*@todo error or similar*/
   }
   return pfile;
}

/** gets a pointer to the Seeed with given ID
 *
 * @returns SeeedPtr to Seeed or NULL if there is no Seeed with the given ID
 * @returns pool: Seeedpool* where the Seeed was found
 */
SeeedPtr MiscVisualization::GCGgetSeeed(
   SCIP* scip,       /**< SCIP data structure */
   int seeedid,      /**< ID of Seeed */
   Seeedpool* pool   /**< outputs where the Seeed was found (if not needed input NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SeeedPtr seeed;
   Seeedpool* seeedpool;

   /* get Seeed from seeedid */
   seeed = NULL;
   conshdlr = SCIPfindConshdlr( scip, "decomp" );

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot find Seeed!\n");
      return NULL;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeedpool = conshdlrdata->seeedpool;
   pool = seeedpool;

   if( seeedpool != NULL )
   {
      /* find in presolved */
      for( size_t i = 0; i < seeedpool->ancestorseeeds.size(); ++i)
      {
         if( seeedpool->ancestorseeeds[i]!= NULL && seeedpool->ancestorseeeds[i]->getID() == seeedid )
            return seeedpool->ancestorseeeds[i];
      }

      for( size_t i = 0; i < seeedpool->incompleteSeeeds.size(); ++i)
      {
         if( seeedpool->incompleteSeeeds[i]->getID() == seeedid )
            return seeedpool->incompleteSeeeds[i];
      }

      for( size_t i = 0; i < seeedpool->finishedSeeeds.size(); ++i)
      {
         if( seeedpool->finishedSeeeds[i]->getID() == seeedid )
            return seeedpool->finishedSeeeds[i];
      }
   }

   /* find in unpresolved */
   seeedpool = conshdlrdata->seeedpoolunpresolved;
   pool = seeedpool;

   if( seeedpool != NULL )
   {
      for( size_t i = 0; i < seeedpool->incompleteSeeeds.size(); ++i)
      {
         if( seeedpool->incompleteSeeeds[i]->getID() == seeedid )
            return seeedpool->incompleteSeeeds[i];
      }

      for( size_t i = 0; i < seeedpool->ancestorseeeds.size(); ++i)
      {
         if( seeedpool->ancestorseeeds[i]!= NULL &&  seeedpool->ancestorseeeds[i]->getID() == seeedid )
            return seeedpool->ancestorseeeds[i];
      }

      for( size_t i = 0; i < seeedpool->finishedSeeeds.size(); ++i)
      {
         if( seeedpool->finishedSeeeds[i]->getID() == seeedid )
            return seeedpool->finishedSeeeds[i];
      }
   }

   pool = NULL;
   return seeed;
}

} /* namespace gcg */

