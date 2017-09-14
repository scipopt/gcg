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

/**@file   reader_cls.cpp
 * @brief  CLS reader for writing files containing classifier data
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


#define READER_NAME             "clsreader"
#define READER_DESC             "reader for writing classifier data"
#define READER_EXTENSION        "cls"


/** data for dec reader */
struct SCIP_ReaderData
{
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

   return SCIP_OKAY;
}
