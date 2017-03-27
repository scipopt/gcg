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
#include <stdlib.h>
#include <unistd.h>

#include "reader_tex.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "reader_gp.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "struct_decomp.h"


#define READER_NAME             "texreader"
#define READER_DESC             "file reader for writing decomposition details to LaTeX files"
#define READER_EXTENSION        "tex"

#define DEFAULT_USEGP            FALSE
#define DEFAULT_MAXNDECOMPS      50
#define DEFAULT_RETURNTYPE       0
#define DEFAULT_PICTURESONLY     FALSE
#define DEFAULT_DRAFTMODE        FALSE

/** data for dec reader */
struct SCIP_ReaderData
{
   SCIP_Bool       usegp;           /** if true uses gp files as intermediate step */
   int             maxndecomps;     /** maximum number of decompositions to visualize
                                      * (ones with best score first are preferred) */
   int             returntype;      /** output only decompositions of type:47: error
                                      * 0=all types, 1=arrowhead, 2=staircase, 3=diagonal, 4=bordered */
   SCIP_Bool       picturesonly;    /** if true only tex code for the pictures is generated
                                      * (no statistics, no report file) */
   SCIP_Bool       draftmode;       /** if true shows no non-zeroes, recommended if too slow or too memory-intensive */
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
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL,
         "Please read in a problem before reading in the corresponding structure file!\n");
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
   SCIP_READERDATA* readerdata;

   assert(scip != NULL);
   assert(reader != NULL);

   ndecomps = SCIPconshdlrDecompGetNDecdecomps(scip);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL( GCGwriteDecompsToTex(scip, file, SCIPconshdlrDecompGetDecdecomps(scip), &ndecomps, TRUE, TRUE, readerdata) );
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

/** gets number of decompositions of a certain type in a given decomposition structure */
static
SCIP_RETCODE getNDecompsOfType(
   SCIP*                scip,               /**< SCIP data structure */
   DEC_DECOMP**         decomps,            /**< Decompositions structure */
   int*                 ndecomps,           /**< Number of decompositions in the structure */
   DEC_DECTYPE          type,               /**< type that is to be counted */
   int*                 number              /**< number of decomps of given type (resultpointer) */
   )
{
   int i;

   *number = 0;
   for( i = 0; i < *ndecomps; i++ )
   {
      if( DECdecompGetType(decomps[i]) == type )
         *number = *number+1;
   }
   return SCIP_OKAY;
}

/** gets path of file */
static
SCIP_RETCODE getPath(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< Decompositions structure */
   char*                pfile               /**< return path of file */
   )
{
   char sympath[SCIP_MAXSTRLEN];
   int filedesc;
   int success;

   filedesc = fileno(file); /* get link to file descriptor */
   if( filedesc < 0 )
   {
      return SCIP_FILECREATEERROR;
   }
   snprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
   success = readlink(sympath, pfile, SCIP_MAXSTRLEN); /* get actual path including extension */
   if( success < 0 )
   {
      return SCIP_NOFILE;
   }
   return SCIP_OKAY;
}

/** write LaTeX code header, begin of document, general statistics and table of contents */
static
SCIP_RETCODE writeHeaderCode(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File pointer to write to */
   SCIP_Bool            statistics,         /**< if true detection statistics and are included in report */
   DEC_DECOMP**         decomps,            /**< Decompositions structure */
   int*                 ndecomps,           /**< Number of decompositions */
   SCIP_Bool            toc,                /**< if true table of contents is included */
   SCIP_READERDATA*     readerdata          /**< reader specific arguments */
   )
{
   char* pname;
   char ppath[SCIP_MAXSTRLEN];
   int ndecompsoftype;

   strcpy(ppath, (char*) SCIPgetProbName(scip));
   SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);

   SCIPinfoMessage(scip, file, "%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% *                  This file is part of the program                         * \n");
   SCIPinfoMessage(scip, file, "%% *          GCG --- Generic Column Generation                                * \n");
   SCIPinfoMessage(scip, file, "%% *                  a Dantzig-Wolfe decomposition based extension            * \n");
   SCIPinfoMessage(scip, file, "%% *                  of the branch-cut-and-price framework                    * \n");
   SCIPinfoMessage(scip, file, "%% *         SCIP --- Solving Constraint Integer Programs                      * \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% * Copyright (C) 2010-2016 Operations Research, RWTH Aachen University       * \n");
   SCIPinfoMessage(scip, file, "%% *                         Zuse Institute Berlin (ZIB)                       * \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% * This program is free software; you can redistribute it and/or             * \n");
   SCIPinfoMessage(scip, file, "%% * modify it under the terms of the GNU Lesser General Public License        * \n");
   SCIPinfoMessage(scip, file, "%% * as published by the Free Software Foundation; either version 3            * \n");
   SCIPinfoMessage(scip, file, "%% * of the License, or (at your option) any later version.                    * \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% * This program is distributed in the hope that it will be useful,           * \n");
   SCIPinfoMessage(scip, file, "%% * but WITHOUT ANY WARRANTY; without even the implied warranty of            * \n");
   SCIPinfoMessage(scip, file, "%% * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             * \n");
   SCIPinfoMessage(scip, file, "%% * GNU Lesser General Public License for more details.                       * \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% * You should have received a copy of the GNU Lesser General Public License  * \n");
   SCIPinfoMessage(scip, file, "%% * along with this program; if not, write to the Free Software               * \n");
   SCIPinfoMessage(scip, file, "%% * Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301, USA.* \n");
   SCIPinfoMessage(scip, file, "%% *                                                                           * \n");
   SCIPinfoMessage(scip, file, "%% * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * \n");
   SCIPinfoMessage(scip, file, "%%                                                                               \n");
   SCIPinfoMessage(scip, file, "%% @author Hanna Franzen                                                         \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "\\documentclass[a4paper,10pt]{article}                                           \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "%% packages                                                                      \n");
   SCIPinfoMessage(scip, file, "\\usepackage[utf8]{inputenc}                                                     \n");
   SCIPinfoMessage(scip, file, "\\usepackage[hidelinks]{hyperref}                                                \n");
   SCIPinfoMessage(scip, file, "\\usepackage{tikz}                                                               \n");
   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, file, "\\usepackage{gnuplot-lua-tikz}                                                \n");
   }
   SCIPinfoMessage(scip, file, " \\usetikzlibrary{external}                                                      \n");
   SCIPinfoMessage(scip, file, " \\tikzexternalize                                                               \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "\\begin{document}                                                                \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "\\begin{titlepage}                                                               \n");
   SCIPinfoMessage(scip, file, "  \\centering                                                                    \n");
   SCIPinfoMessage(scip, file, "  \\thispagestyle{empty}                                                         \n");
   SCIPinfoMessage(scip, file, "  {\\Huge Report: %s} \\\\ \\today                                               \n",
      pname);

   if( statistics )
   {
      SCIPinfoMessage(scip, file, "                                                                              \n");
      SCIPinfoMessage(scip, file, "\\vspace{2cm}                                                                 \n");
      SCIPinfoMessage(scip, file, "\\begin{tabular}{ll}                                                          \n");
      SCIPinfoMessage(scip, file, "  \\textbf{Problem}: & \\begin{minipage}{0pt}                                 \n");
      SCIPinfoMessage(scip, file, "                         \\begin{verbatim}%s\\end{verbatim}                   \n",
         pname);
      SCIPinfoMessage(scip, file, "                       \\end{minipage} \\\\                                   \n");
      SCIPinfoMessage(scip, file, "  Number of variables in original problem: & %i  \\\\                         \n",
         SCIPgetNOrigVars(scip));
      SCIPinfoMessage(scip, file, "  \\vspace{0.5cm}                                                             \n");
      SCIPinfoMessage(scip, file, "  Number of constraints in original problem: & %i  \\\\                       \n",
         SCIPgetNOrigConss(scip));
      SCIPinfoMessage(scip, file, "  Number of found decompositions: & %i  \\\\                                  \n",
         SCIPconshdlrDecompGetNDecdecomps(scip));
      if( readerdata->returntype != 0 )
      {
          getNDecompsOfType(scip,decomps,ndecomps,readerdata->returntype, &ndecompsoftype);
          SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\          \n",
             ndecompsoftype);
      }
      else
      {
         SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\           \n",
            *ndecomps);
      }
      SCIPinfoMessage(scip, file, "\\end{tabular}                                                                \n");
      SCIPinfoMessage(scip, file, "                                                                              \n");
   }
   SCIPinfoMessage(scip, file, "\\end{titlepage}                                                                 \n");

   if( toc && readerdata->picturesonly == FALSE)
   {
      SCIPinfoMessage(scip, file, "\\thispagestyle{empty}                                                        \n");
      SCIPinfoMessage(scip, file, "\\tableofcontents                                                             \n");
      SCIPinfoMessage(scip, file, "\\newpage                                                                     \n");
   }

   return SCIP_OKAY;
}

/** writes the code for a Tikz visualization of the decomposition into the file
 * works analogously to the SCIPwriteGp function in reader_gp.c */
static
SCIP_RETCODE writeTikz(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp,             /**< Decomposition array pointer */
   SCIP_READERDATA*      readerdata          /**< reader specific arguments */
   )
{
   SCIP_VAR*** subscipvars;
   SCIP_CONS*** subscipconss;
   SCIP_VAR** linkingvars;
   SCIP_CONS** linkingconss;
   SCIP_VAR** vars;
   SCIP_CONS** conss;
   SCIP_HASHMAP* varindexmap;
   SCIP_HASHMAP* consindexmap;
   int* nsubscipvars;
   int* nsubscipconss;
   int* nstairlinkingvars;
   size_t varindex = 1;
   size_t consindex = 1;
   int startx = 0;
   int starty = 0;
   int endx = 0;
   int endy = 0;
   int nlinkingvars;
   int nlinkingconss;
   int i;
   int j;
   int nvars;
   int nconss;
   int maxindvars = 0;
   int maxindcons = 0;
   int maxind = 0;
   double radius = 5;
   float xpoint;
   float ypoint;

   assert(scip != NULL);

   subscipvars = DECdecompGetSubscipvars(decomp);
   nsubscipvars = DECdecompGetNSubscipvars(decomp);
   subscipconss = DECdecompGetSubscipconss(decomp);
   nsubscipconss = DECdecompGetNSubscipconss(decomp);
   linkingvars = DECdecompGetLinkingvars(decomp);
   nlinkingvars = DECdecompGetNLinkingvars(decomp);
   linkingconss = DECdecompGetLinkingconss(decomp);
   nlinkingconss = DECdecompGetNLinkingconss(decomp);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);

   /* --- compute indices for variables & constraints --- */

   varindexmap = NULL;
   consindexmap = NULL;

   if( decomp != NULL )
   {
      /* go through the blocks and create the indices */
      /* @todo add " && DECdecompGetType(decomp) != DEC_DECTYPE_ARROWHEAD" to this in seeed version */
      if( DECdecompGetType(decomp) != DEC_DECTYPE_UNKNOWN && DECdecompGetType(decomp) != DEC_DECTYPE_STAIRCASE )
      {
         SCIP_CALL( SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
            SCIP_CALL( SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );
         for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
         {
            for( j = 0; j < nsubscipvars[i]; ++j )
            {
               assert(subscipvars[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(varindexmap, subscipvars[i][j], (void*)varindex) );
               if( (int)varindex > maxindvars )
                  maxindvars = (int) varindex;
               varindex++;
            }

            for( j = 0; j < nsubscipconss[i]; ++j )
            {
               assert(subscipconss[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(consindexmap, subscipconss[i][j], (void*)consindex) );
               if( (int)consindex > maxindcons )
                  maxindcons = (int) consindex;
               consindex++;
            }
         }

         for( j = 0; j < nlinkingvars; ++j )
         {
            assert(linkingvars[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(varindexmap, linkingvars[j], (void*)varindex) );
            if( (int)varindex > maxindvars )
               maxindvars = (int) varindex;
            varindex++;
         }
         for( j = 0; j < nlinkingconss; ++j )
         {
            assert(linkingconss[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(consindexmap, linkingconss[j], (void*)consindex) );
            if( (int)consindex > maxindcons )
               maxindcons = (int) consindex;
            consindex++;
         }
      }
      /* @todo add " || DECdecompGetType(decomp) == DEC_DECTYPE_ARROWHEAD" to this in seeed version */
      else if( DECdecompGetType(decomp) == DEC_DECTYPE_STAIRCASE )
      {
         /* get the computed index maps instead of computing them here */
         varindexmap = DECdecompGetVarindex(decomp);
         consindexmap = DECdecompGetConsindex(decomp);

         /* determine max indices */
         for(i = 0; i < nvars; i++)
         {
            if(SCIPhashmapExists(varindexmap, vars[i]) &&
               (int)(size_t)(SCIPhashmapGetImage(varindexmap, vars[i]))>maxindvars)
            {
               maxindvars = (int)(size_t)(SCIPhashmapGetImage(varindexmap, vars[i]));
            }
         }
         for(i = 0; i < nconss; i++)
         {
            if(SCIPhashmapExists(consindexmap, conss[i]) &&
               (int)(size_t)(SCIPhashmapGetImage(consindexmap, conss[i]))>maxindcons)
            {
               maxindcons = (int)(size_t)(SCIPhashmapGetImage(consindexmap, conss[i]));
            }
         }

         assert(varindexmap != NULL);
         assert(consindexmap != NULL);
      }
   }

   /* the max indices must be at least one to be compatible with division */
   maxindvars = maxindvars<1?1:maxindvars;
   maxindcons = maxindcons<1?1:maxindcons;
   /* determine the highest index */
   maxind = maxindvars>maxindcons?maxindvars:maxindcons;

   /* --- write header --- */

   SCIPinfoMessage(scip, file, "  \\resizebox{\\textwidth}{!}{                                                   \n");
   SCIPinfoMessage(scip, file, "  \\begin{tikzpicture}                                                           \n");

   /* --- draw grey rectangles with standard outline (black) for the blocks --- */
   /* note: the picture is scaled to the page's textwidth in order to scale down large pictures.
    * Instead of var-/consindex the value of (index/maxindex)*textwidth/height is used
    */

   if( DECdecompGetType(decomp) == DEC_DECTYPE_ARROWHEAD || DECdecompGetType(decomp) == DEC_DECTYPE_BORDERED )
   {
      for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
      {
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file,
            "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
            (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
         startx = endx;
         starty = endy;
      }
      endx += nlinkingvars;
      endy += nlinkingconss;
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
         (0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
         (startx+0.5)/maxindvars, (+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
         (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
   }
   else
   {
      if( DECdecompGetType(decomp) == DEC_DECTYPE_STAIRCASE )
      {
         nstairlinkingvars = DECdecompGetNStairlinkingvars(decomp);
         for( i = 0; i < DECdecompGetNBlocks(decomp)-1; ++i )
         {
            endx += nsubscipvars[i]+nstairlinkingvars[i];
            endy += nsubscipconss[i];
            SCIPinfoMessage(scip, file,
               "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
               (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
            startx = endx-nstairlinkingvars[i];
            starty = endy;
         }
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file,
            "    \\draw [fill=gray] (%f*\\textwidth,%f*\\textheight) rectangle (%f*\\textwidth,%f*\\textheight);\n",
            (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
      }
   }

   /* --- draw black dots for nonzeroes --- */

   /* draw the dots */
   if(!readerdata->draftmode)
   {
      for( i = 0; i < nconss; i++ )
      {
         int ncurvars = GCGconsGetNVars(scip, conss[i]);
         SCIP_VAR** curvars = NULL;

         if( ncurvars > 0 )
         {
            SCIP_CALL( SCIPallocBufferArray( scip, &curvars, ncurvars) );
            SCIP_CALL( GCGconsGetVars(scip, conss[i], curvars, ncurvars) );
         }

         for( j = 0; j < ncurvars; j++ )
         {
            /* if the problem has been created but has not been processed yet, output the whole model */
            if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
            {
               SCIPinfoMessage(scip, file,
                  "                                                                                \n");
               SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth,%f*\\textheight) circle [radius=%f];\n",
                  (SCIPvarGetIndex(curvars[j]))/maxindvars, (i)/maxindcons, radius/maxind);
            }
            else
            {
               /* if there is no decomposition, output the presolved model! */
               if( decomp == NULL || DECdecompGetType(decomp) == DEC_DECTYPE_UNKNOWN )
               {
                  SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth,%f*\\textheight) circle [radius=%f];\n",
                     (SCIPvarGetIndex(curvars[j]))/maxindvars, (i)/maxindcons, radius/maxind);
               }
               /* if there is a decomposition, output the indices derived from the decomposition above*/
               else
               {
                  assert(varindexmap != NULL);
                  assert(consindexmap != NULL);
                  /*@todo make the following if statement into an assertion*/
                  if( SCIPhashmapExists(varindexmap, SCIPvarGetProbvar(curvars[j]))
                     && SCIPhashmapExists(consindexmap, conss[i]))
                  {
                     xpoint =
                        ( (float)(size_t)SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])) )/(float)maxindvars;
                     ypoint = ( (float)(size_t)SCIPhashmapGetImage(consindexmap, conss[i]) )/ (float)maxindcons;
                     SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth,%f*\\textheight) circle [radius=%f];\n",
                        xpoint, ypoint, radius/maxind);
                  }
               }
            }
         }

         SCIPfreeBufferArrayNull(scip, &curvars);
      }
   }

   SCIPinfoMessage(scip, file, "                                                                                \n");

   /* --- write closing --- */

   SCIPinfoMessage(scip, file, "  \\end{tikzpicture}                                                            \n");
   SCIPinfoMessage(scip, file, "  }                                                                             \n");

   /* @todo add " && DECdecompGetType(decomp) != DEC_DECTYPE_ARROWHEAD" to this in seeed version */
   if( DECdecompGetType(decomp) != DEC_DECTYPE_STAIRCASE )
   {
      SCIPhashmapFree(&varindexmap);
      SCIPhashmapFree(&consindexmap);
   }

   return SCIP_OKAY;
}

/** write LaTeX code for one decomposition */
static
SCIP_RETCODE writeDecompCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp,             /**< Decomposition array pointer */
   SCIP_READERDATA*      readerdata          /**< reader specific arguments */
   )
{
   FILE* gpfile;
   DEC_DETECTOR** detectorchain;
   DEC_SCORES scores;
   char* filepath;
   char* pname;
   char ppath[SCIP_MAXSTRLEN];
   char decompname[SCIP_MAXSTRLEN];
   char gpfilename[SCIP_MAXSTRLEN];
   char gpname[SCIP_MAXSTRLEN];
   char pfile[SCIP_MAXSTRLEN];
   char pfilecpy[SCIP_MAXSTRLEN];
   char dectype[SCIP_MAXSTRLEN];
   char detectorchainstring[SCIP_MAXSTRLEN];
   int sizedetectorchain;
   int i;

   /* construct detector chain string*/

   detectorchain = DECdecompGetDetectorChain(decomp);
   sizedetectorchain = DECdecompGetDetectorChainSize(decomp);

   sprintf(detectorchainstring, "%s", DECdetectorGetName(detectorchain[0]));

   for( i=1; i < sizedetectorchain; ++i )
   {
      sprintf(detectorchainstring, "%s-%s",detectorchainstring, DECdetectorGetName(detectorchain[i]) );
   }
   SCIPinfoMessage(scip, NULL, "%s \n", detectorchainstring);

   assert(decomp != NULL);

   (void) SCIPsnprintf(decompname, SCIP_MAXSTRLEN, "%c-%d", detectorchainstring, DECdecompGetNBlocks(decomp));

   if( readerdata->usegp )
   {
      /* --- create a gnuplot file for the decomposition --- */

      /* get path to write to and put it into gpfilename */
      getPath(scip, file, pfile);
      strcpy(pfilecpy, pfile);
      SCIPsplitFilename(pfilecpy, &filepath, NULL, NULL, NULL);
      strcpy(gpfilename, filepath);
      strcat(gpfilename, "/");

      /* get name of file and attach it to gpfilename */
      strcpy(ppath, (char*) SCIPgetProbName(scip));
      SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);
      if( pname != NULL &&  pname[0] != '\0' )
      {
         strcat(gpfilename, pname);
         strcat(gpfilename, "-");
      }

      if( decompname != NULL &&  decompname[0] != '\0' )
      {
         strcat(gpfilename, decompname);
      }
      else
      {
         return SCIP_FILECREATEERROR;
      }
      strcpy(gpname, gpfilename);
      strcat(gpfilename, ".gp");

      /* write gp file for decomp using the gp reader (using the tex output option) */
      gpfile = fopen(gpfilename, "w");
      if( gpfile == NULL )
      {
         return SCIP_FILECREATEERROR;
      }

      SCIPwriteGp(scip, gpfile, decomp, TRUE, FALSE);

      fclose(gpfile);
   }

   /* --- gather information & output them into .tex file --- */


   DECevaluateDecomposition(scip, decomp, &scores);

   if(!readerdata->picturesonly)
   {
      SCIPinfoMessage(scip, file, "\\section*{Decomposition: %s}                                   \n", decompname);
      SCIPinfoMessage(scip, file, "\\addcontentsline{toc}{section}{Decomposition: %s}              \n", decompname);
      SCIPinfoMessage(scip, file, "                                                                \n");
   }
   SCIPinfoMessage(scip, file, "\\begin{figure}[!htb]                                              \n");
   SCIPinfoMessage(scip, file, "  \\begin{center}                                                  \n");
   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, file, "    \\input{%s-%c-%d}                                            \n",
         pname, detectorchainstring, DECdecompGetNBlocks(decomp));
   }
   else
   {
      writeTikz(scip, file, decomp, readerdata);
   }

   SCIPinfoMessage(scip, file, "  \\end{center}                                                    \n");
   SCIPinfoMessage(scip, file, "\\end {figure}                                                     \n");
   if(!readerdata->picturesonly)
   {
      SCIPinfoMessage(scip, file, "                                                                \n");
      SCIPinfoMessage(scip, file, "\\vspace{0.3cm}                                                 \n");
      SCIPinfoMessage(scip, file, "\\begin{tabular}{ll}                                            \n");
      SCIPinfoMessage(scip, file, "  Found by detector: & %s \\\\                                  \n",detectorchainstring);
      switch(DECdecompGetType(decomp))
      {
         case DEC_DECTYPE_ARROWHEAD:
            strcpy(dectype,"arrowhead");
            break;
         case DEC_DECTYPE_STAIRCASE:
            strcpy(dectype,"staircase");
            break;
         case DEC_DECTYPE_DIAGONAL:
            strcpy(dectype,"diagonal");
            break;
         case DEC_DECTYPE_BORDERED:
            strcpy(dectype,"bordered");
            break;
         default:
            strcpy(dectype,"unknown");
            break;
      }
      SCIPinfoMessage(scip, file, "  Type of decomposition: & %s \\\\                                              \n",
         dectype);
      SCIPinfoMessage(scip, file, "  Number of blocks: & %i \\\\                                                   \n",
         DECdecompGetNBlocks(decomp));
      SCIPinfoMessage(scip, file, "  Number of linking variables: & %i \\\\                                        \n",
         DECdecompGetNLinkingvars(decomp));
      SCIPinfoMessage(scip, file, "  Number of linking constraints: & %i \\\\                                      \n",
         DECdecompGetNLinkingconss(decomp));
      SCIPinfoMessage(scip, file, "  Block density score: & %f \\\\                                                \n",
         scores.densityscore);
      SCIPinfoMessage(scip, file, "  Interlinking blocks score: & %f \\\\                                          \n",
         scores.linkingscore);
      SCIPinfoMessage(scip, file, "  Border score: & %f \\\\                                                       \n",
         scores.borderscore);
      SCIPinfoMessage(scip, file, "  \\textbf{Total score:} & \\textbf{%f} \\\\                                    \n",
         scores.totalscore);
      SCIPinfoMessage(scip, file, "\\end{tabular}                                                                  \n");
   }
   SCIPinfoMessage(scip, file, "\\clearpage                                                                     \n");
   SCIPinfoMessage(scip, file, "                                                                                \n");

   return SCIP_OKAY;
}

/** write LaTeX code for end of document */
static
SCIP_RETCODE writeEndCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{

   SCIPinfoMessage(scip, file, "\\end{document}                                                                  \n");

   return SCIP_OKAY;
}

/** makes a new makefile and readme for the given .tex file */
static
SCIP_RETCODE makeMakefileAndReadme(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< Decompositions structure */
   SCIP_READERDATA*     readerdata          /**< reader specific arguments */
   )
{
   FILE* makefile;
   FILE* readme;
   char* filepath;
   char* filename;
   char pfile[SCIP_MAXSTRLEN];
   char pfilecpy[SCIP_MAXSTRLEN];
   char makefilename[SCIP_MAXSTRLEN];
   char readmename[SCIP_MAXSTRLEN];
   char name[SCIP_MAXSTRLEN];
   const char makename[SCIP_MAXSTRLEN] = "makepdf";

   /* --- create a Makefile --- */

   /* get path to write to and put it into makefilename */
   getPath(scip, file, pfile);
   strcpy(pfilecpy, pfile);
   SCIPsplitFilename(pfilecpy, &filepath, &filename, NULL, NULL);
   strcpy(makefilename, filepath);
   strcat(makefilename, "/");
   strcpy(name, makename);
   strcat(name, "_");
   strcat(name, filename);
   strcat(name, ".make");
   strcat(makefilename, name);


   /* open and write makefile */
   makefile = fopen(makefilename, "w");
   if( makefile == NULL )
   {
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, makefile, "                                                                             \n");
   SCIPinfoMessage(scip, makefile, "# latexmk automatically manages the .tex files                               \n");
   SCIPinfoMessage(scip, makefile, "%s.pdf: %s.tex                                                               \n",
      filename, filename);
   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
      SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
      SCIPinfoMessage(scip, makefile, "\t@echo Compiling gp files to tex                                            \n");
      SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
      SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
      SCIPinfoMessage(scip, makefile, "\tgnuplot *.gp                                                               \n");
   }
   SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
   SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
   SCIPinfoMessage(scip, makefile, "\t@echo Compiling tex code. This may take a while.                           \n");
   SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
   SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
   SCIPinfoMessage(scip, makefile,
      "\t@latexmk -pdf -pdflatex=\"pdflatex -interaction=batchmode -shell-escape\" -use-make %s.tex \n", filename);
   SCIPinfoMessage(scip, makefile, "\t@make -f %s clean                                                          \n", name);
   SCIPinfoMessage(scip, makefile, "                                                                             \n");
   SCIPinfoMessage(scip, makefile, "clean:                                                                       \n");
   SCIPinfoMessage(scip, makefile, "\t@latexmk -c                                                                \n");
   SCIPinfoMessage(scip, makefile, "\t@rm -f report_*figure*.*                                                   \n");
   SCIPinfoMessage(scip, makefile, "\t@rm -f *.auxlock                                                           \n");
   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, makefile, "\t@rm -f *.gp                                                             \n");
   }
   SCIPinfoMessage(scip, makefile, "                                                                             \n");
   SCIPinfoMessage(scip, makefile, "cleanall:                                                                    \n");
   SCIPinfoMessage(scip, makefile, "\t@latexmk -C                                                                \n");
   SCIPinfoMessage(scip, makefile, "\t@make -f %s clean                                                          \n", name);

   /* close makefile */
   fclose(makefile);

   /* --- create a readme file --- */

   /* use same file path as the makefile */
   strcpy(readmename, filepath);
   strcat(readmename, "/");
   strcat(readmename, "README_");
   strcat(readmename, makename);

   /* open and write readme */
   readme = fopen(readmename, "w");
   if( readme == NULL )
   {
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, readme, "README: How to create a PDF file from the .tex file(s) using the %s file     \n", name);
   SCIPinfoMessage(scip, readme, "                                                                             \n");
   SCIPinfoMessage(scip, readme, "Instead of using 'make' use 'make -f %s' instead                             \n", name);
   SCIPinfoMessage(scip, readme, "                                                                             \n");
   SCIPinfoMessage(scip, readme, "Clean options:                                                               \n");
   SCIPinfoMessage(scip, readme, "\t'clean' clears all present intermediate files (if any exist)               \n");
   SCIPinfoMessage(scip, readme, "\t'cleanall' clears all generated files INCLUDING .pdf                       \n");

   /* close readme file */
   fclose(readme);
   return SCIP_OKAY;
}

/** writes tex files for the visualization & statistics of a given set of decomposition
 * and writes a Makefile to compile the files with
 */
SCIP_RETCODE GCGwriteDecompsToTex(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP**          decomps,            /**< Decomposition array pointer */
   int*                  ndecomps,           /**< Number of decompositions */
   SCIP_Bool             statistics,         /**< if true detection statistics and are included in report */
   SCIP_Bool             toc,                /**< if true table of contents is included */
   SCIP_READERDATA*      readerdata          /**< reader specific arguments */
   )
{
   FILE* decompfile;
   char* filepath;
   char* filename;
   char pfile[SCIP_MAXSTRLEN];
   char decompname[SCIP_MAXSTRLEN];
   char tempname[SCIP_MAXSTRLEN] = {'\0'};
   char tempstr[SCIP_MAXSTRLEN] = {'\0'};
   char tempc = '\0';
   SCIP_Bool writedecomp;
   int i;
   int maxrounds;
   int ndecompsoftype;

   assert(scip != NULL);
   assert(*ndecomps > 0);

   getPath(scip, file, pfile);
   SCIPsplitFilename(pfile, &filepath, &filename, NULL, NULL);

   /* --- make a makefile and readme file --- */
   makeMakefileAndReadme(scip,file,readerdata);

   /* --- make the tex file(s) --- */

   /* write LaTeX header including title and (optional) statistics & table of contents */
   SCIP_CALL( writeHeaderCode(scip,file,statistics,decomps,ndecomps,toc,readerdata) );

   if( readerdata->returntype != 0 )
   {
      getNDecompsOfType(scip,decomps,ndecomps,readerdata->returntype, &ndecompsoftype);
   }
   else
   {
      ndecompsoftype = *ndecomps;
   }

   /* check if the number of max decomps exceeds the number of available outputs */
   if( readerdata->maxndecomps < ndecompsoftype )
   {
      maxrounds = readerdata->maxndecomps;
   }
   else
   {
      maxrounds = *ndecomps;
   }

   /* write LaTeX code for each decomp starting with the highest score */
   /* note: decomps come sorted from lowest to highest score */
   /* only output such decompositions of the given type */
   for( i = 0; i < *ndecomps && maxrounds > 0; i++ )
   {
      if( decomps[i] != NULL )
      {
         writedecomp = FALSE;
         if( readerdata->returntype == 0 )
            writedecomp = TRUE;
         else if( (unsigned int)readerdata->returntype == DECdecompGetType(decomps[i]) )
            writedecomp = TRUE;

         if( writedecomp == TRUE )
         {
            if(readerdata->picturesonly)
            {
               /* get file path and attach detectorchar + nblocks */
               strcpy(decompname, filepath);
               strcat(decompname, "/");

               strcpy(tempname, filename);
               strcat(tempname, "-");
               tempc = DECdetectorGetChar(DECdecompGetDetector(decomps[i]));
               strcat(tempname, &tempc);
               strcat(tempname, "-");
               sprintf(tempstr,"%d",DECdecompGetNBlocks(decomps[i]));
               strcat(tempname, tempstr);
               tempstr[0] = '\0';

               strcat(decompname, tempname);
               decompfile = fopen(decompname, "w");
               if( decompfile == NULL )
               {
                  return SCIP_FILECREATEERROR;
               }
               /* write decomposition picture into new file */
               SCIP_CALL( writeDecompCode(scip,decompfile,decomps[i],readerdata) );
               fclose(decompfile);

               /* input the decomposition into main file */
               SCIPinfoMessage(scip, file, "    \\input{%s}                                        \n",tempname);
               tempname[0] = '\0';
            }
            else
            {
               /* if picturesonly is false, put decomposition infos into main file */
               SCIP_CALL( writeDecompCode(scip,file,decomps[i],readerdata) );
            }
            maxrounds--;
         }
      }
   }

   /*write an ending for the LaTeX code*/
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

   /* create tex reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include tex reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeTex, readerReadTex, readerWriteTex, readerdata));

   /* include possible parameters */
   SCIP_CALL( SCIPaddBoolParam(scip,
      "reading/texreader/usegp", "if true uses gp files as intermediate step",
      &readerdata->usegp, FALSE, DEFAULT_USEGP, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/texreader/picturesonly",
         "if true only tex code for the pictures is generated (no statistics, no report file)",
         &readerdata->picturesonly, FALSE, DEFAULT_PICTURESONLY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip,
         "reading/texreader/draftmode",
         "if true shows no non-zeroes, recommended if too slow or too memory-intensive",
         &readerdata->draftmode, FALSE, DEFAULT_DRAFTMODE, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "reading/texreader/maxndecomps",
      "maximum number of decompositions to visualize (ones with best score are preferred)",
      &readerdata->maxndecomps, FALSE, DEFAULT_MAXNDECOMPS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip,
      "reading/texreader/returntype",
      "output only decompositions of type 0=all types, 1=arrowhead, 2=staircase, 3=diagonal, 4=bordered",
      &readerdata->returntype, FALSE, DEFAULT_RETURNTYPE, 0, 4, NULL, NULL) );

   return SCIP_OKAY;
}
