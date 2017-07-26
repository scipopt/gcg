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

/**@file   reader_tex.cpp
 * @brief  tex file reader for writing seeeds to LaTeX files
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

#include "reader_tex.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "reader_gp.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "struct_decomp.h"
#include "class_miscvisualization.h"


#define READER_NAME             "texreader"
#define READER_DESC             "LaTeX file writer for seeed visualization"
#define READER_EXTENSION        "tex"

#define DEFAULT_USEGP            FALSE
#define DEFAULT_PICTURESONLY     FALSE

/** data for dec reader */
struct SCIP_ReaderData
{
   SCIP_Bool       usegp;           /** if true uses gp files as intermediate step */
   SCIP_Bool       picturesonly;    /** if true only tex code for the pictures is generated (no statistics included) */
};

/** Getter of parameter usegp */
SCIP_Bool GCGtexGetUseGp(
   SCIP*                scip               /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   readerdata = SCIPreaderGetData(SCIPfindReader(scip, "texreader"));
   return readerdata->usegp;
}

/** Getter of parameter picturesonly */
SCIP_Bool GCGtexGetPicturesonly(
   SCIP*                scip               /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;
   readerdata = SCIPreaderGetData(SCIPfindReader(scip, "texreader"));
   return readerdata->picturesonly;
}

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

/** Problem reading method of reader.
 *  Since the reader is not supposed to read files this returns a reading error. */
static
SCIP_DECL_READERREAD(readerReadTex)
{  /*lint --e{715}*/
   return SCIP_READERROR;
}

/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteTex)
{  /*lint --e{715}*/

   assert(scip != NULL);
   assert(reader != NULL);

   SCIP_CALL( GCGtexWriteHeaderCode(scip,file) );
   SCIP_CALL( GCGtexWriteDecompCode(scip, file, DECgetBestDecomp(scip)) );
   SCIP_CALL( GCGtexWriteEndCode(scip,file) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/* outputs the r, g, b decimal values for the rgb hex input */
static
SCIP_RETCODE getRgbDecFromHex(
   char*    hex,     /**< input hex rgb code of form "#000000" */
   int*     red,     /**< output decimal r */
   int*     green,   /**< output decimal g */
   int*     blue     /**< output decimal b */
   )
{
   char* temp = '\0';
   int check = 0;
   unsigned int r = 0;
   unsigned int g = 0;
   unsigned int b = 0;

   assert( hex[0] == '#' );

   /* remove the # at the beginning */
   strcpy( temp, hex );
   memmove( temp, temp+1, strlen( temp ) );

   /* extract int values from the rest */
   check = sscanf( temp, "%02x%02x%02x", &r, &g, &b );
   assert( check == 3 );

   *red = (int) r;
   *green = (int) g;
   *blue = (int) b;

   return SCIP_OKAY;
}

/** write LaTeX code header & begin of document
 * The proper order in which a tex file is written goes as follows:
 *    -> GCGtexWriteHeaderCode         (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteHeaderCode(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file                /**< File pointer to write to */
   )
{
   SCIP_READERDATA* readerdata;
   char* pname;
   char ppath[SCIP_MAXSTRLEN];

   readerdata = SCIPreaderGetData(SCIPfindReader(scip, "texreader"));

   strcpy(ppath, SCIPgetProbName(scip));
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

   /*@todo add defines of the current default colors here (use getRgbDecFromHex) */

   return SCIP_OKAY;
}

/** write LaTeX code title page that includes general statistics about the problem
 *  * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    -> GCGtexWriteTitlepage          (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteTitlepage(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File pointer to write to */
   int*                 npresenteddecomps   /**< Number of decompositions to be shown in the file or NULL if unknown */
   )
{
   char* pname;
   char ppath[SCIP_MAXSTRLEN];
   int ndecomps;

   ndecomps = SCIPconshdlrDecompGetNFinishedDecomps(scip);
   strcpy(ppath, SCIPgetProbName(scip));
   SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);

   SCIPinfoMessage(scip, file, "\\begin{titlepage}                                                            \n");
   SCIPinfoMessage(scip, file, "  \\centering                                                                 \n");
   SCIPinfoMessage(scip, file, "  \\thispagestyle{empty}                                                      \n");
   SCIPinfoMessage(scip, file, "  {\\Huge Report: %s} \\\\ \\today                                            \n",
      pname);
   SCIPinfoMessage(scip, file, "                                                                              \n");
   SCIPinfoMessage(scip, file, "\\vspace{2cm}                                                                 \n");
   SCIPinfoMessage(scip, file, "\\begin{tabular}{{lp{10cm}}}                                                  \n");
   SCIPinfoMessage(scip, file, "  \\textbf{Problem}: & \\begin{minipage}{10cm}                                \n");
   SCIPinfoMessage(scip, file, "                         \\begin{verbatim}%s\\end{verbatim}                   \n",
      pname);
   SCIPinfoMessage(scip, file, "                       \\end{minipage} \\\\                                   \n");
   SCIPinfoMessage(scip, file, "  Number of variables in original problem: & %i  \\\\                         \n",
      SCIPgetNOrigVars(scip));
   SCIPinfoMessage(scip, file, "  \\vspace{0.5cm}                                                             \n");
   SCIPinfoMessage(scip, file, "  Number of constraints in original problem: & %i  \\\\                       \n",
      SCIPgetNOrigConss(scip));
   SCIPinfoMessage(scip, file, "  Number of found decompositions: & %i  \\\\                                  \n",
      SCIPconshdlrDecompGetNFinishedDecomps(scip));
   if(npresenteddecomps != NULL){
      if( ndecomps > *npresenteddecomps )
      {
         SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\ \n",
            *npresenteddecomps);
      }
      else
      {
         SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\ \n", ndecomps);
      }
   }
   SCIPinfoMessage(scip, file, "\\end{tabular}                                                                \n");
   SCIPinfoMessage(scip, file, "                                                                              \n");
   SCIPinfoMessage(scip, file, "\\end{titlepage}                                                              \n");
   SCIPinfoMessage(scip, file, "\\newpage                                                                     \n");

   return SCIP_OKAY;
}

/** write LaTeX code for table of contents
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    -> GCGtexWriteTableOfContents    (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteTableOfContents(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file                /**< File pointer to write to */
   )
{
   SCIPinfoMessage(scip, file, "\\thispagestyle{empty}                                                        \n");
   SCIPinfoMessage(scip, file, "\\tableofcontents                                                             \n");
   SCIPinfoMessage(scip, file, "\\newpage                                                                     \n");
   SCIPinfoMessage(scip, file, "                                                                              \n");

   return SCIP_OKAY;
}

/** writes the code for a Tikz visualization of the decomposition into the file
 * works analogously to the SCIPwriteGp function in reader_gp.c */
static
SCIP_RETCODE writeTikz(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp,             /**< Decomposition pointer */
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

   SCIPinfoMessage(scip, file, "  \\begin{tikzpicture}                                                           \n");

   /* --- draw grey rectangles with standard outline (black) for the blocks --- */
   /* note: the picture is scaled to the page's textwidth in order to scale down large pictures.
    * Instead of var-/consindex the value of (index/maxindex)*textwidth/height is used
    */

   if( DECdecompGetType(decomp) == DEC_DECTYPE_ARROWHEAD || DECdecompGetType(decomp) == DEC_DECTYPE_BORDERED
       || DECdecompGetType(decomp) == DEC_DECTYPE_DIAGONAL )
   {
      for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
      {
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file,
            "    \\draw [fill=gray] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
            (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
         startx = endx;
         starty = endy;
      }
      endx += nlinkingvars;
      endy += nlinkingconss;
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=orange] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
         (0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=yellow] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
         (startx+0.5)/maxindvars, (+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
      SCIPinfoMessage(scip, file,
         "    \\draw [fill=purple] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
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
               "    \\draw [fill=pink] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
               (startx+0.5)/maxindvars, (starty+0.5)/maxindcons, (endx+0.5)/maxindvars, (endy+0.5)/maxindcons);
            startx = endx-nstairlinkingvars[i];
            starty = endy;
         }
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file,
            "    \\draw [fill=gray] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
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
               SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) circle [radius=%f*0.75];\n",
                  (SCIPvarGetIndex(curvars[j]))/maxindvars, (i)/maxindcons, radius/maxind);
            }
            else
            {
               /* if there is no decomposition, output the presolved model! */
               if( decomp == NULL || DECdecompGetType(decomp) == DEC_DECTYPE_UNKNOWN )
               {
                  SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) circle [radius=%f*0.75];\n",
                     (SCIPvarGetIndex(curvars[j]))/maxindvars, (i)/maxindcons, radius/maxind);
               }
               /* if there is a decomposition, output the indices derived from the decomposition above*/
               else
               {
                  assert(varindexmap != NULL);
                  assert(consindexmap != NULL);
                  if( SCIPhashmapExists(varindexmap, SCIPvarGetProbvar(curvars[j]))
                     && SCIPhashmapExists(consindexmap, conss[i]))
                  {
                     xpoint =
                        ( (float)(size_t)SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])) )/(float)maxindvars;
                     ypoint = ( (float)(size_t)SCIPhashmapGetImage(consindexmap, conss[i]) )/ (float)maxindcons;
                     SCIPinfoMessage(scip, file, "    \\draw [fill] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) circle [radius=%f*0.75];\n",
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

   if( DECdecompGetType(decomp) != DEC_DECTYPE_STAIRCASE )
   {
      SCIPhashmapFree(&varindexmap);
      SCIPhashmapFree(&consindexmap);
   }

   return SCIP_OKAY;
}

/** write LaTeX code for one decomposition
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    -> GCGtexWriteDecompCode         (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteDecompCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp              /**< Decomposition pointer */
   )
{
   FILE* gpfile;
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   DEC_DETECTOR** detectorchain;
   DEC_SCORES scores;
   SCIP_Bool gpdraftwason = FALSE;
   char* filepath;
   char* pname;
   char* detectorchainstring;
   char ppath[SCIP_MAXSTRLEN];
   char decompname[SCIP_MAXSTRLEN];
   char gpfilename[SCIP_MAXSTRLEN];
   char gpname[SCIP_MAXSTRLEN];
   char pfile[SCIP_MAXSTRLEN];
   char pfilecpy[SCIP_MAXSTRLEN];
   char dectype[SCIP_MAXSTRLEN];
   char fulldetectorstring[SCIP_MAXSTRLEN];
   int sizedetectorchain;
   int i;

   assert(decomp != NULL);
   reader = SCIPfindReader(scip, "texreader");
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   /* get detector chain string & full-text string*/
   detectorchainstring = DECdecompGetDetectorChainString(scip, decomp);

   detectorchain = DECdecompGetDetectorChain(decomp);
   sizedetectorchain = DECdecompGetDetectorChainSize(decomp);
   if( detectorchain[0] != NULL)
      sprintf(fulldetectorstring, "%s", DECdetectorGetName(detectorchain[0]));
   else
      sprintf(fulldetectorstring, "%s", "user");
   for( i=1; i < sizedetectorchain; ++i )
   {
      sprintf(fulldetectorstring, "%s, %s",fulldetectorstring, DECdetectorGetName(detectorchain[i]) );
   }

   (void) SCIPsnprintf(decompname, SCIP_MAXSTRLEN, "%s-%d-%d", detectorchainstring, DECdecompGetSeeedID(decomp),
      DECdecompGetNBlocks(decomp));
   /* tex will have problems with the character '_' */
   for(i = 0; i < SCIP_MAXSTRLEN; i++)
   {
      if(decompname[i] == '_'){
         decompname[i] = '-';
      }
   }

   if( readerdata->usegp )
   {
      /* --- create a gnuplot file for the decomposition --- */

      /* get path to write to and put it into gpfilename */
      GCGgetFilePath(file, pfile);
      strcpy(pfilecpy, pfile);
      SCIPsplitFilename(pfilecpy, &filepath, NULL, NULL, NULL);
      strcpy(gpfilename, filepath);
      strcat(gpfilename, "/");

      /* get name of file and attach it to gpfilename */
      strcpy(ppath, SCIPgetProbName(scip));
      SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);
      if( pname != NULL &&  pname[0] != '\0' )
      {
         strcpy(gpname, pname);
         strcat(gpname, "-");
      }

      if( decompname != NULL &&  decompname[0] != '\0' )
      {
         strcat(gpname, decompname);
      }
      else
      {
         return SCIP_FILECREATEERROR;
      }
      strcat(gpfilename, gpname);
      strcat(gpfilename, ".gp");

      /* write gp file for decomp using the gp reader (using the tex output option) */
      gpfile = fopen(gpfilename, "w");
      if( gpfile == NULL )
      {
         return SCIP_FILECREATEERROR;
      }

      /* write gp in the tex draft mode and restore the original parameter afterwards */
      gpdraftwason = GCGgpGetDraftmode(scip);
      GCGgpSetDraftmode(scip, readerdata->draftmode);

      SCIPwriteGp(scip, gpfile, decomp, TRUE, FALSE);

      GCGgpSetDraftmode(scip, gpdraftwason);

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
      SCIPinfoMessage(scip, file, "    \\input{%s-%s-%d-%d}                                          \n",
         pname, detectorchainstring, DECdecompGetSeeedID(decomp), DECdecompGetNBlocks(decomp));
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
      SCIPinfoMessage(scip, file, "\\begin{tabular}{lp{10cm}}                                      \n");
      SCIPinfoMessage(scip, file,
         "  Found by detector(s): & \\begin{minipage}{10cm}\\begin{verbatim}%s\\end{verbatim}\\end{minipage} \\\\ \n",
         fulldetectorstring);
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

/** write LaTeX code for end of document
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    -> GCGtexWriteEndCode            (required)
 *    GCGtexWriteMakefileAndReadme     (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteEndCode(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{

   SCIPinfoMessage(scip, file, "\\end{document}                                                                  \n");

   return SCIP_OKAY;
}

/** makes a new makefile and readme for the given .tex file
 * The proper order in which a tex file is written goes as follows:
 *    GCGtexWriteHeaderCode            (required)
 *    GCGtexWriteTitlepage             (optional)
 *    GCGtexWriteTableOfContents       (optional)
 *    GCGtexWriteDecompCode            (required as often as the number of decompositions you wish to visualize)
 *    GCGtexWriteEndCode               (required)
 *    -> GCGtexWriteMakefileAndReadme  (optional but highly recommended)
 */
SCIP_RETCODE GCGtexWriteMakefileAndReadme(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file                /**< File for which the makefile & readme are generated */
   )
{
   FILE* makefile;
   FILE* readme;
   SCIP_READERDATA* readerdata;
   char* filepath;
   char* filename;
   char pfile[SCIP_MAXSTRLEN];
   char pfilecpy[SCIP_MAXSTRLEN];
   char makefilename[SCIP_MAXSTRLEN];
   char readmename[SCIP_MAXSTRLEN];
   char name[SCIP_MAXSTRLEN];
   const char makename[SCIP_MAXSTRLEN] = "makepdf";

   readerdata = SCIPreaderGetData(SCIPfindReader(scip, "texreader"));

   /* --- create a Makefile --- */

   /* get path to write to and put it into makefilename */
   GCGgetFilePath(file, pfile);
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

   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, makefile, "GPFILES := $(wildcard *.gp)\n");
   }
   SCIPinfoMessage(scip, makefile, "                                                                             \n");
   SCIPinfoMessage(scip, makefile, "# latexmk automatically manages the .tex files                               \n");
   SCIPinfoMessage(scip, makefile, "%s.pdf: %s.tex\n",
      filename, filename);
   if( readerdata->usegp )
   {
      SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
      SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
      SCIPinfoMessage(scip, makefile, "\t@echo Compiling gp files to tex                                            \n");
      SCIPinfoMessage(scip, makefile, "\t@echo                                                                      \n");
      SCIPinfoMessage(scip, makefile, "\t@echo ------------                                                         \n");
      SCIPinfoMessage(scip, makefile, "\t$(SHELL) -ec  'for i in $(GPFILES); \\\n");
      SCIPinfoMessage(scip, makefile, "\t\tdo \\\n");
      SCIPinfoMessage(scip, makefile, "\t\tgnuplot $$i; \\\n");
      SCIPinfoMessage(scip, makefile, "\t\tdone'\n");
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

   return SCIP_OKAY;
}
