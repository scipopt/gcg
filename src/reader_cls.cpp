/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2018 Operations Research, RWTH Aachen University       */
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

/**@file   reader_cls.cpp
 * @brief  CLS reader for writing files containing classifier data
 * @author Michael Bastubbe
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strcasecmp() */
#endif
#include <ctype.h>

#include "reader_cls.h"
#include "cons_decomp.h"
#include "class_seeedpool.h"
#include "class_consclassifier.h"
#include "class_varclassifier.h"



#define READER_NAME             "clsreader"
#define READER_DESC             "reader for writing classifier data"
#define READER_EXTENSION        "cls"
#define DEFAULT_USETRANSFORM    TRUE

struct SCIP_ConshdlrData
{
   gcg::Seeedpool* seeedpoolunpresolved;
   gcg::Seeedpool* seeedpool;
};



/** data for dec reader */
struct SCIP_ReaderData
{
   SCIP_Bool usetransform;
};

/*
 * Local methods
 */


/** write a DEC file for a given decomposition */
SCIP_RETCODE GCGwriteCls(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   assert(scip != NULL);

   /* TODO
    * format description:
    * a) <number of classifiers>
    * for each classifier:
    *    b1) VAR or CONS
    *    b2) <name of classifier>>
    *    b3) <number of classes>
    *    b4) for each class:
    *       c1) <name of class>: <description of class>
    *       c2) <number of class elements>
    *       c3) for each element of class:
    *          d1) <name of element> (e.g. variable or constraint name, concerning transformed [default] or original problem)
    *
    *
    * 1. Write cons_decomp method(s) to get all relevant data, i.e.
    *    - for all cons and var classifier in both presolved and unpresolved seeedpool:
    *      - Class to index mapping, class names, class descriptions (and maybe class decomp info)
    *      - Cons/var to classes mapping
    *    - index to SCIPcons/SCIPvar mapping from unpresolved seeedpool (and maybe seeedpool)
    *
    * 2. Write an output method to write these information into a file
    *
    */
   SCIP_Bool transformed;
   gcg::Seeedpool* seeedpool;


   assert(scip != NULL);

   SCIP_CALL( SCIPgetBoolParam(scip,
         "reading/clsreader/usetransform", &transformed));


   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      transformed = FALSE;

   if( !transformed && SCIPconshdlrDecompGetSeeedpoolUnpresolvedExtern(scip) == NULL )
      SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);

   if( transformed && SCIPconshdlrDecompGetSeeedpoolExtern(scip) == NULL )
      SCIPconshdlrDecompCreateSeeedpool(scip);

   seeedpool = (gcg::Seeedpool*)(transformed ? SCIPconshdlrDecompGetSeeedpoolExtern(scip) : SCIPconshdlrDecompGetSeeedpoolUnpresolvedExtern(scip));


   SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);

   if( seeedpool->consclassescollection.size() == 0 )
      seeedpool->calcClassifierAndNBlockCandidates(scip);

   SCIPinfoMessage(scip, file, "# a1) <number of classifiers>\n" );
   SCIPinfoMessage(scip, file, "# a2) for each classifier:\n" );
   SCIPinfoMessage(scip, file, "# b1)    VAR or CONS\n" );
   SCIPinfoMessage(scip, file, "# b2)    <name of classifier>\n" );
   SCIPinfoMessage(scip, file, "# b3)    <number of classes>\n" );
   SCIPinfoMessage(scip, file, "# b4)    for each class:\n" );
   SCIPinfoMessage(scip, file, "# c1)       <name of class>: <description of class>\n" );
   SCIPinfoMessage(scip, file, "# c2)       <number of class elements>\n" );
   SCIPinfoMessage(scip, file, "# c3)       for each element of class:\n" );
   SCIPinfoMessage(scip, file, "# d1)          <name of element> (e.g. variable or constraint name, concerning transformed [default] or original problem)\n" );
   SCIPinfoMessage(scip, file, "###########################################\n" );


   /** a */
   SCIPinfoMessage(scip, file, "%d\n", (int) seeedpool->consclassescollection.size() + (int) seeedpool->varclassescollection.size() );

   for( size_t c = 0; c < seeedpool->consclassescollection.size() ; ++c )
   {
      gcg::ConsClassifier* classifier = seeedpool->consclassescollection[c];

      std::vector<std::vector<int> > conssofclasses = std::vector<std::vector<int> >(classifier->getNClasses()) ;
      for( int cons = 0; cons < seeedpool->getNConss(); ++cons )
         conssofclasses[classifier->getClassOfCons(cons)].push_back(cons);

      /** b1 */
      SCIPinfoMessage(scip, file, "CONS\n" );
      /** b2 */
      SCIPinfoMessage(scip, file, "%s \n", classifier->getName());
      /** b3 */
      SCIPinfoMessage(scip, file, "%d\n", classifier->getNClasses() );
      for( int cl = 0; cl < classifier->getNClasses(); ++cl )
      {
         /** c1 */
         SCIPinfoMessage(scip, file, "%s: %s\n", classifier->getClassName(cl), classifier->getClassDescription(cl) );
         /** c2 */
         SCIPinfoMessage(scip, file, "%d\n",  conssofclasses[cl].size() );
         /** c3 */
         for( size_t clm = 0; clm < conssofclasses[cl].size(); ++clm )
         {
            SCIPinfoMessage(scip, file, "%s\n",  SCIPconsGetName( seeedpool->getConsForIndex( conssofclasses[cl][clm])) );
         }
      }
   }


   for( size_t c = 0; c < seeedpool->varclassescollection.size() ; ++c )
   {
      gcg::VarClassifier* classifier = seeedpool->varclassescollection[c];

      std::vector<std::vector<int> > varsofclasses = std::vector<std::vector<int> >(classifier->getNClasses()) ;
      for( int var = 0; var < seeedpool->getNVars(); ++var )
         varsofclasses[classifier->getClassOfVar(var)].push_back(var);

      /** b1 */
      SCIPinfoMessage(scip, file, "VAR\n" );
      /** b2 */
      SCIPinfoMessage(scip, file, "%s \n", classifier->getName());
      /** b3 */
      SCIPinfoMessage(scip, file, "%d\n", classifier->getNClasses() );
      for( int cl = 0; cl < classifier->getNClasses(); ++cl )
      {
         /** c1 */
         SCIPinfoMessage(scip, file, "%s: %s\n", classifier->getClassName(cl), classifier->getClassDescription(cl) );
         /** c2 */
         SCIPinfoMessage(scip, file, "%d\n",  classifier->getNVarsOfClasses()[cl] );
         /** c3 */
         for( size_t clm = 0; clm <varsofclasses[cl].size(); ++clm )
         {
            SCIPinfoMessage(scip, file, "%s\n",  SCIPvarGetName( seeedpool->getVarForIndex( varsofclasses[cl][clm])) );
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
   SCIP_CALL( GCGwriteCls( scip, file ) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** includes the cls reader into SCIP */
SCIP_RETCODE SCIPincludeReaderCls(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create cls reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include cls reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
      readerCopyCls, readerFreeCls, readerReadCls, readerWriteCls, readerdata) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "reading/clsreader/usetransform",
      "should the transformed (and possibly presolved problem) be use or original one",
      &readerdata->usetransform, FALSE, DEFAULT_USETRANSFORM, NULL, NULL) );


   return SCIP_OKAY;
}
