/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   reader_cls.cpp
 * @brief  CLS reader for writing files containing classification data
 * @author Michael Bastubbe
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#include <ctype.h>

#include "gcg/reader_cls.h"
#include "gcg/cons_decomp.h"
#include "gcg/cons_decomp.hpp"
#include "gcg/class_detprobdata.h"
#include "gcg/class_conspartition.h"
#include "gcg/class_varpartition.h"


#define READER_NAME             "clsreader"
#define READER_DESC             "reader for writing classification data"
#define READER_EXTENSION        "cls"
#define DEFAULT_USETRANSFORM    TRUE

struct SCIP_ConshdlrData
{
};



/** data for cls reader */
struct SCIP_ReaderData
{
   GCG* gcg;
   SCIP_Bool usetransform;
};

/*
 * Local methods
 */


/** write classification data */
static
SCIP_RETCODE writeCls(
   GCG*                  gcg,                /**< GCG data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   SCIP* scip;
   SCIP_Bool transformed;
   gcg::DETPROBDATA* detprobdata;

   assert(gcg != NULL);

   scip = GCGgetOrigprob(gcg);

   SCIP_CALL( SCIPgetBoolParam(scip,
         "reading/clsreader/usetransform", &transformed));

   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      transformed = FALSE;

   if( !transformed )
      detprobdata = GCGconshdlrDecompGetDetprobdataOrig(gcg);
   else
      detprobdata = GCGconshdlrDecompGetDetprobdataPresolved(gcg);

   if( detprobdata->conspartitioncollection.empty() )
   {
      GCGconshdlrDecompClassify(gcg, !detprobdata->isAssignedToOrigProb());
      GCGconshdlrDecompCalcCandidatesNBlocks(gcg, !detprobdata->isAssignedToOrigProb());
   }

   SCIPinfoMessage(scip, file, "# a1) <number of partitions>\n" );
   SCIPinfoMessage(scip, file, "# a2) for each partition:\n" );
   SCIPinfoMessage(scip, file, "# b1)    VAR or CONS\n" );
   SCIPinfoMessage(scip, file, "# b2)    <name of partition>\n" );
   SCIPinfoMessage(scip, file, "# b3)    <number of classes>\n" );
   SCIPinfoMessage(scip, file, "# b4)    for each class:\n" );
   SCIPinfoMessage(scip, file, "# c1)       <name of class>: <description of class>\n" );
   SCIPinfoMessage(scip, file, "# c2)       <number of class elements>\n" );
   SCIPinfoMessage(scip, file, "# c3)       for each element of class:\n" );
   SCIPinfoMessage(scip, file, "# d1)          <name of element> (e.g. variable or constraint name, concerning transformed [default] or original problem)\n" );
   SCIPinfoMessage(scip, file, "###########################################\n" );

   /* a */
   SCIPinfoMessage(scip, file, "%d\n", (int) detprobdata->conspartitioncollection.size() + (int) detprobdata->varpartitioncollection.size() );

   for( size_t c = 0; c < detprobdata->conspartitioncollection.size() ; ++c )
   {
      gcg::ConsPartition* partition = detprobdata->conspartitioncollection[c];

      std::vector<std::vector<int> > conssofclasses = std::vector<std::vector<int> >(partition->getNClasses()) ;
      for( int cons = 0; cons < detprobdata->getNConss(); ++cons )
         conssofclasses[partition->getClassOfCons(cons)].push_back(cons);

      /* b1 */
      SCIPinfoMessage(scip, file, "CONS\n" );
      /* b2 */
      SCIPinfoMessage(scip, file, "%s \n", partition->getName());
      /* b3 */
      SCIPinfoMessage(scip, file, "%d\n", partition->getNClasses());
      for( int cl = 0; cl < partition->getNClasses(); ++cl )
      {
         /* c1 */
         SCIPinfoMessage(scip, file, "%s: %s\n", partition->getClassName(cl), partition->getClassDescription(cl));
         /* c2 */
         SCIPinfoMessage(scip, file, "%ld\n",  conssofclasses[cl].size());
         /* c3 */
         for( size_t clm = 0; clm < conssofclasses[cl].size(); ++clm )
         {
            SCIPinfoMessage(scip, file, "%s\n",  SCIPconsGetName(detprobdata->getCons(conssofclasses[cl][clm])));
         }
      }
   }


   for( size_t c = 0; c < detprobdata->varpartitioncollection.size() ; ++c )
   {
      gcg::VarPartition* partition = detprobdata->varpartitioncollection[c];

      std::vector<std::vector<int> > varsofclasses = std::vector<std::vector<int> >(partition->getNClasses()) ;
      for( int var = 0; var < detprobdata->getNVars(); ++var )
         varsofclasses[partition->getClassOfVar(var)].push_back(var);

      /* b1 */
      SCIPinfoMessage(scip, file, "VAR\n" );
      /* b2 */
      SCIPinfoMessage(scip, file, "%s \n", partition->getName());
      /* b3 */
      SCIPinfoMessage(scip, file, "%d\n", partition->getNClasses() );
      for( int cl = 0; cl < partition->getNClasses(); ++cl )
      {
         /* c1 */
         SCIPinfoMessage(scip, file, "%s: %s\n", partition->getClassName(cl), partition->getClassDescription(cl));
         /* c2 */
         SCIPinfoMessage(scip, file, "%d\n", partition->getNVarsOfClasses()[cl] );
         /* c3 */
         for( size_t clm = 0; clm <varsofclasses[cl].size(); ++clm )
         {
            SCIPinfoMessage(scip, file, "%s\n",  SCIPvarGetName(detprobdata->getVar(varsofclasses[cl][clm])));
         }
      }
   }


   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

#define readerCopyCls NULL

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeCls)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData( reader );
   assert( readerdata != NULL );

   SCIPfreeMemory( scip, &readerdata );

   assert( strcmp( SCIPreaderGetName( reader ), READER_NAME ) == 0);
   return SCIP_OKAY;
}

/** problem reading method of reader */
#define readerReadCls NULL

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteCls)
{
   /*lint --e{715}*/
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData( reader );
   assert( readerdata != NULL );
   SCIP_CALL( writeCls( readerdata->gcg, file ) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** includes the cls reader into SCIP */
SCIP_RETCODE GCGincludeReaderCls(
   GCG*                  gcg                 /**< GCG data structure */
   )
{
   SCIP_READERDATA* readerdata = NULL;
   SCIP* origprob = GCGgetOrigprob(gcg);
   assert(origprob != NULL);

   /* create cls reader data */
   SCIP_CALL( SCIPallocMemory(origprob, &readerdata) );
   readerdata->gcg = gcg;

   /* include cls reader */
   SCIP_CALL( SCIPincludeReader(origprob, READER_NAME, READER_DESC, READER_EXTENSION,
      readerCopyCls, readerFreeCls, readerReadCls, readerWriteCls, readerdata) );

   SCIP_CALL( SCIPaddBoolParam(origprob,
      "reading/clsreader/usetransform",
      "should the transformed (and possibly presolved problem) be use or original one",
      &readerdata->usetransform, FALSE, DEFAULT_USETRANSFORM, NULL, NULL) );


   return SCIP_OKAY;
}
