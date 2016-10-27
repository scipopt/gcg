/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       */
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

/**@file   reader_tex.c
 * @brief  tex file reader for writing decomposition details to LaTeX files
 * @author Hanna Franzen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>
#if defined(_WIN32) || defined(_WIN64)
#else
#include <strings.h> /*lint --e{766}*/ /* needed for strcasecmp() */
#endif
#include <ctype.h>
#include <stdio.h>

#include "reader_tex.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"

#include "cons_decomp.h"
#include "pub_decomp.h"

#define READER_NAME             "texreader"
#define READER_DESC             "file reader for writing decomposition details to LaTeX files"
#define READER_EXTENSION        "tex"

#if defined(_WIN32) || defined(_WIN64)
#define LINEBREAK "\r\n"
#else
#define LINEBREAK "\n"
#endif

#define CONSHDLR_DEFAULT        "decomp"

/** data for dec reader */
struct SCIP_ReaderData
{
   const char*       conshdlr; /** constraint handler containing decomposition data */
};

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeTex)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}

/** problem reading method of reader */
static
SCIP_DECL_READERREAD(readerReadTex)
{  /*lint --e{715}*/
   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "Please read in a problem before reading in the corresponding structure file!\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( GCGreadTex(scip, filename, result) );

   return SCIP_OKAY;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteTex)
{  /*lint --e{715}*/
   int ndecomps;

   assert(scip != NULL);
   assert(reader != NULL);

   ndecomps = SCIPconshdlrDecompGetNDecdecomps(scip);

   SCIP_CALL( GCGwriteDecompsToTex(scip, file, SCIPconshdlrDecompGetDecdecomps(scip), &ndecomps, TRUE) );
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}

/* the reader is not supposed to read files,
 * returns a reading error */
SCIP_RETCODE GCGreadTex(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< full path and name of file to read, or NULL if stdin should be used */
   SCIP_RESULT*          result              /**< pointer to store the result of the file reading call */
   )
{
   return SCIP_READERROR;
}

/** write LaTeX code header & begin of document */
static
SCIP_RETCODE writeHeaderCode(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file                /**< File pointer to write to */
   )
{
   SCIPinfoMessage(scip, file, "%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                  This file is part of the program                         * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *          GCG --- Generic Column Generation                                * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                  a Dantzig-Wolfe decomposition based extension            * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                  of the branch-cut-and-price framework                    * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *         SCIP --- Solving Constraint Integer Programs                      * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                         Zuse Institute Berlin (ZIB)                       * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * This program is free software; you can redistribute it and/or             * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * modify it under the terms of the GNU Lesser General Public License        * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * as published by the Free Software Foundation; either version 3            * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * of the License, or (at your option) any later version.                    * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * This program is distributed in the hope that it will be useful,           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * but WITHOUT ANY WARRANTY; without even the implied warranty of            * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * GNU Lesser General Public License for more details.                       * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * You should have received a copy of the GNU Lesser General Public License  * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * along with this program; if not, write to the Free Software               * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.* %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% *                                                                           * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%%                                                                               %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% @author Hanna Franzen                                                         %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\documentclass[a4paper,10pt]{article}                                           %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "%% packages                                                                      %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\usepackage[utf8]{inputenc}                                                     %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{document}                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);

   return SCIP_OKAY;
}

/** write LaTeX code for general decomposition statistics */
static
SCIP_RETCODE writeGeneralStatisticsCode(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File pointer to write to */
   DEC_DECOMP**         decomps,            /**< Decompositions structure */
   int*                 ndecomps            /**< Number of decompositions */
   )
{

/*   SCIPinfoMessage(scip, file, "\\section*{Detection Statistics}                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\addcontentsline{toc}{section}{Detection Statistics}                           %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{tabular}{ll}                                                            %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  \\textbf{Problem}: & %s \\\\                                                  %s", SCIPgetProbName(scip), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Number of found decompositions: & %s \\\\                                     %s", SCIPconshdlrDecompGetNDecdecomps(scip), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\                %s", *ndecomps, LINEBREAK);
   SCIPinfoMessage(scip, file, "\\end{tabular}                                                                  %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\vspace{0.3cm}                                                                 %s", LINEBREAK);
*/
   GCGprintDetectorStatistics(scip, file);
/*
   SCIPinfoMessage(scip, file, "\\vspace{0.3cm}                                                                 %s", LINEBREAK);
*/
   /*@todo get and output more statistics*/



   return SCIP_OKAY;
}

/** write LaTeX code for one decomposition */
static
SCIP_RETCODE writeDecompCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp              /**< Decomposition array pointer */
   )
{
   /*
   SCIPinfoMessage(scip, file, "\\newline                                                                       %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\section*{Decomposition: %s}                                                   %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\addcontentsline{toc}{section}{Detection Statistics}                           %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   @todo get and output statistics
*/
   return SCIP_OKAY;
}

/** write LaTeX code for end of document */
static
SCIP_RETCODE writeEndCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{

   SCIPinfoMessage(scip, file, "\\end{document}                                                                  %s", LINEBREAK);

   return SCIP_OKAY;
}

/** write a tex file for the visualization & statistics of a given set of decomposition */
SCIP_RETCODE GCGwriteDecompsToTex(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP**          decomps,            /**< Decomposition array pointer */
   int*                  ndecomps,           /**< Number of decompositions */
   SCIP_Bool             statistics          /**< if true detection statistics are included in report */
   )
{
   DEC_DECOMP** sorteddecomps;
   int i;

   assert(scip != NULL);
   assert(*ndecomps > 0);

   /*@todo sort decomps into sorteddecomps (just rearrange pointers)*/
   sorteddecomps = decomps;

   SCIP_CALL( writeHeaderCode(scip,file) );

   if(statistics)
   {
      SCIP_CALL( writeGeneralStatisticsCode(scip,file,sorteddecomps,ndecomps) );
   }

   for( i=0; i<*ndecomps; i++ )
   {
      SCIP_CALL( writeDecompCode(scip,file,decomps[i]) );
   }

   SCIP_CALL( writeEndCode(scip,file) );

   return SCIP_OKAY;
}

/** includes the tex file reader in SCIP */
SCIP_RETCODE
SCIPincludeReaderTex(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create dec reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include dec reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeTex, readerReadTex, readerWriteTex, readerdata));

   SCIP_CALL( SCIPaddStringParam(scip,
         "reading/texreader/parameters", "constraint handler containing decomposition data",
         NULL, FALSE, CONSHDLR_DEFAULT, NULL, NULL) );

   return SCIP_OKAY;
}
