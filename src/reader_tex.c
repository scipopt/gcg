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
#include <unistd.h>

#include "reader_tex.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "reader_gp.h"
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

#define PATH_DEFAULT            " "

/** data for dec reader */
struct SCIP_ReaderData
{
   char*       path;     /** path to main file */
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

   SCIP_CALL( GCGwriteDecompsToTex(scip, file, SCIPconshdlrDecompGetDecdecomps(scip), &ndecomps, TRUE, TRUE, TRUE) );
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

/** write LaTeX code header, begin of document, general statistics and table of contents */
static
SCIP_RETCODE writeHeaderCode(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File pointer to write to */
   SCIP_Bool            statistics,         /**< if true detection statistics and are included in report */
   DEC_DECOMP**         decomps,            /**< Decompositions structure */
   int*                 ndecomps,           /**< Number of decompositions */
   SCIP_Bool            toc,                /**< if true table of contents is included */
   SCIP_Bool            useGp               /**< if true Gnuplot will be used for visualization */
   )
{
   char* pname;
   char* ppath;

   ppath = (char*) SCIPgetProbName(scip);
   SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);

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
   SCIPinfoMessage(scip, file, "\\usepackage[hidelinks]{hyperref}                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\usepackage{tikz}                                                               %s", LINEBREAK);
   if(useGp)
   {
      SCIPinfoMessage(scip, file, "\\usepackage{gnuplot-lua-tikz}                                                   %s", LINEBREAK);
   }
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{document}                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{titlepage}                                                               %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  \\centering                                                                    %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  \\thispagestyle{empty}                                                         %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  {\\Huge Report: %s} \\\\ \\today                                               %s", pname, LINEBREAK);

   if(statistics)
   {
      SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "\\vspace{2cm}                                                                    %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "\\begin{tabular}{ll}                                                             %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "  \\textbf{Problem}: & \\begin{minipage}{0pt}                                    %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "                         \\begin{verbatim}%s\\end{verbatim}                      %s", pname, LINEBREAK);
      SCIPinfoMessage(scip, file, "                       \\end{minipage} \\\\                                      %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "  Number of found decompositions: & %i  \\\\                                     %s", SCIPconshdlrDecompGetNDecdecomps(scip), LINEBREAK);
      SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\                 %s", *ndecomps, LINEBREAK);
      SCIPinfoMessage(scip, file, "\\end{tabular}                                                                   %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "                                                                                 %s", LINEBREAK);
   }
   SCIPinfoMessage(scip, file, "\\end{titlepage}                                                                 %s", LINEBREAK);

   if(toc)
   {
      SCIPinfoMessage(scip, file, "\\thispagestyle{empty}                                                           %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "\\tableofcontents                                                                %s", LINEBREAK);
      SCIPinfoMessage(scip, file, "\\newpage                                                                        %s", LINEBREAK);
   }

   return SCIP_OKAY;
}

/** writes the code for a Tikz visualization of the decomposition into the file
 * works analogously to the SCIPwriteGp function in reader_gp.c */
static
SCIP_RETCODE writeTikz(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decomp              /**< Decomposition array pointer */
   )
{
   SCIP_VAR*** subscipvars;
   SCIP_CONS*** subscipconss;
   SCIP_VAR** linkingvars;
   SCIP_CONS** linkingconss;
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
   int nconss;
   double radius = 0.5;

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

   /* --- write header --- */

   SCIPinfoMessage(scip, file, "  \\resizebox{\\textwidth}{!}{                                                  %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  \\begin{tikzpicture}                                                          %s", LINEBREAK);

   /* --- draw grey rectangles with standard outline (black) for the blocks --- */

   if( DECdecompGetType(decomp) == DEC_DECTYPE_ARROWHEAD || DECdecompGetType(decomp) == DEC_DECTYPE_BORDERED )
   {
      for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
      {
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", startx+0.5, starty+0.5, endx+0.5, endy+0.5, LINEBREAK);
         startx = endx;
         starty = endy;
      }
      endx += nlinkingvars;
      endy += nlinkingconss;
      SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", 0.5, starty+0.5, endx+0.5, endy+0.5, LINEBREAK);
      SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", startx+0.5, +0.5, endx+0.5, endy+0.5, LINEBREAK);
      SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", startx+0.5, starty+0.5, endx+0.5, endy+0.5, LINEBREAK);
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
            SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", startx+0.5, starty+0.5, endx+0.5, endy+0.5, LINEBREAK);
            startx = endx-nstairlinkingvars[i];
            starty = endy;
         }
         endx += nsubscipvars[i];
         endy += nsubscipconss[i];
         SCIPinfoMessage(scip, file, "    \\draw [fill=gray] (%f,%f) rectangle (%f,%f);                               %s", startx+0.5, starty+0.5, endx+0.5, endy+0.5, LINEBREAK);
      }
   }

   /* --- draw black dots for the constraints --- */

   /* compute indices */
   varindexmap = NULL;
   consindexmap = NULL;

   if( decomp != NULL )
   {
      /* if we don't have staircase, but something else, go through the blocks and create the indices */
      if( DECdecompGetType(decomp) == DEC_DECTYPE_ARROWHEAD || DECdecompGetType(decomp) == DEC_DECTYPE_BORDERED || DECdecompGetType(decomp) == DEC_DECTYPE_DIAGONAL )
      {
         SCIP_CALL( SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
         SCIP_CALL( SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );

         for( i = 0; i < DECdecompGetNBlocks(decomp); ++i )
         {
            for( j = 0; j < nsubscipvars[i]; ++j )
            {
               assert(subscipvars[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(varindexmap, subscipvars[i][j], (void*)varindex) );
               varindex++;
            }
            for( j = 0; j < nsubscipconss[i]; ++j )
            {
               assert(subscipconss[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(consindexmap, subscipconss[i][j], (void*)consindex) );
               consindex++;
            }
         }

         for( j = 0; j < nlinkingvars; ++j )
         {
            assert(linkingvars[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(varindexmap, linkingvars[j], (void*)varindex) );
            varindex++;
         }
         for( j = 0; j < nlinkingconss; ++j )
         {
            assert(linkingconss[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(consindexmap, linkingconss[j], (void*)consindex) );
            consindex++;
         }
      }
      else if( DECdecompGetType(decomp) == DEC_DECTYPE_STAIRCASE )
      {
         /* in case of staircase get the original indices */
         varindexmap = DECdecompGetVarindex(decomp);
         consindexmap = DECdecompGetVarindex(decomp);

         assert(varindexmap != NULL);
         assert(consindexmap != NULL);
      }
   }

   /* draw the dots */
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
            SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
            SCIPinfoMessage(scip, file, "    \\draw [fill] (%f,%f) circle [radius=%f];                                   %s", SCIPvarGetIndex(curvars[j]), i, radius, LINEBREAK);
         }
         else
         {
            /* if there is no decomposition, output the presolved model! */
            if( decomp == NULL || DECdecompGetType(decomp) == DEC_DECTYPE_UNKNOWN )
            {
               SCIPinfoMessage(scip, file, "    \\draw [fill] (%f,%f) circle [radius=%f];                                   %s", SCIPvarGetIndex(curvars[j]), i, radius, LINEBREAK);
            }
            /* if there is a decomposition, output the indices derived from the decomposition above*/
            else
            {
               assert(varindexmap != NULL);
               assert(consindexmap != NULL);
               /*@todo make the following if statement into an assertion*/
               if(SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])) != NULL && SCIPhashmapGetImage(consindexmap, conss[i]) != NULL)
               {
                  SCIPinfoMessage(scip, file, "    \\draw [fill] (%d,%d) circle [radius=%f];                                   %s", SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])), SCIPhashmapGetImage(consindexmap, conss[i]), radius, LINEBREAK);
               }
            }
         }
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
   }


   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);

   /* --- write closing --- */

   SCIPinfoMessage(scip, file, "  \\end{tikzpicture}                                                            %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  }                                                                             %s", LINEBREAK);

   if( decomp != NULL && DECdecompGetType(decomp) != DEC_DECTYPE_STAIRCASE )
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
   FILE*                 makefile,           /**< File pointer to corresponding makefile */
   DEC_DECOMP*           decomp,             /**< Decomposition array pointer */
   SCIP_Bool             useGp               /**< if true Gnuplot will be used for visualization */
   )
{
   char* filepath;
   char* pname;
   char ppath[SCIP_MAXSTRLEN];
   char decompname[SCIP_MAXSTRLEN];
   char gpfilename[SCIP_MAXSTRLEN];
   char gpname[SCIP_MAXSTRLEN];
   char sympath[SCIP_MAXSTRLEN];
   char pfile[SCIP_MAXSTRLEN];
   FILE* gpfile;
   int filedesc;
   int success;
   DEC_SCORES scores;

   assert(decomp != NULL);
   (void) SCIPsnprintf(decompname, SCIP_MAXSTRLEN, "%c-%d", DECdetectorGetChar(DECdecompGetDetector(decomp)), DECdecompGetNBlocks(decomp));

   if(useGp)
   {
      /* --- create a gnuplot file for the decomposition --- */

      /* get path to write to and put it into gpfilename */
      filedesc = fileno(file); /* get link to file descriptor */
      if(filedesc < 0)
      {
         return SCIP_FILECREATEERROR;
      }
      snprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
      success = readlink(sympath, pfile, SCIP_MAXSTRLEN); /* get actual path including extension */
      if(success < 0)
      {
         return SCIP_NOFILE;
      }
      SCIPsplitFilename(pfile, &filepath, NULL, NULL, NULL);
      strcpy(gpfilename, filepath);
      strcat(gpfilename, "/");

      /* get name of file and attach it to gpfilename */
      strcpy(ppath, (char*) SCIPgetProbName(scip));
      SCIPsplitFilename(ppath, NULL, &pname, NULL, NULL);
      if(pname != NULL &&  pname[0] != '\0')
      {
         strcat(gpfilename, pname);
         strcat(gpfilename, "-");
      }

      if(decompname != NULL &&  decompname[0] != '\0')
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
      if(gpfile == NULL)
      {
         return SCIP_FILECREATEERROR;
      }

      SCIPwriteGp(scip, gpfile, decomp, TRUE, FALSE);

      fclose(gpfile);

      /*@todo add tex of gnuplot file to makefile*/
   }

   /* --- gather further information & output them --- */

   DECevaluateDecomposition(scip, decomp, &scores);

   SCIPinfoMessage(scip, file, "\\section*{Decomposition: %s}                                                   %s", decompname, LINEBREAK);
   SCIPinfoMessage(scip, file, "\\addcontentsline{toc}{section}{Decomposition: %s}                              %s", decompname, LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{figure}[!htb]                                                           %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  \\begin{center}                                                               %s", LINEBREAK);
   if(useGp)
   {
      SCIPinfoMessage(scip, file, "    \\input{%s-%c-%d}                                                           %s", pname, DECdetectorGetChar(DECdecompGetDetector(decomp)), DECdecompGetNBlocks(decomp), LINEBREAK);
   }
   else
   {
      writeTikz(scip, file, decomp);
   }
   SCIPinfoMessage(scip, file, "  \\end{center}                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\end {figure}                                                                  %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\vspace{0.3cm}                                                                 %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\begin{tabular}{lll}                                                           %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "  Found by detector: & %s & \\\\                                                %s", DECdetectorGetName(DECdecompGetDetector(decomp)), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Number of blocks: & %i & \\\\                                                 %s", DECdecompGetNBlocks(decomp), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Number of linking variables: & %i & \\\\                                      %s", DECdecompGetNLinkingvars(decomp), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Number of linking constraints: & %i & \\\\                                    %s", DECdecompGetNLinkingconss(decomp), LINEBREAK);
   SCIPinfoMessage(scip, file, "  Scores: & Total score: & %f \\\\                                              %s", scores.totalscore, LINEBREAK);
   SCIPinfoMessage(scip, file, "  & Block density score: & %f \\\\                                              %s", scores.densityscore, LINEBREAK);
   SCIPinfoMessage(scip, file, "  & Interlinking blocks score: & %f \\\\                                        %s", scores.linkingscore, LINEBREAK);
   SCIPinfoMessage(scip, file, "  & Border score: & %f \\\\                                                     %s", scores.borderscore, LINEBREAK);
   SCIPinfoMessage(scip, file, "\\end{tabular}                                                                  %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "\\newpage                                                                       %s", LINEBREAK);
   SCIPinfoMessage(scip, file, "                                                                                %s", LINEBREAK);

   /*@todo get and output statistics*/

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
   SCIP_Bool             useGp               /**< if true Gnuplot will be used for visualization */
   )
{
   DEC_DECOMP** sorteddecomps;
   FILE* makefile;
   char* filepath;
   char* filename;
   char sympath[SCIP_MAXSTRLEN];
   char pfile[SCIP_MAXSTRLEN];
   char makefilename[SCIP_MAXSTRLEN];
   int filedesc;
   int success;
   int i;

   assert(scip != NULL);
   assert(*ndecomps > 0);

   /* --- create a Makefile --- */

   /* get path to write to and put it into gpfilename */
   filedesc = fileno(file); /* get link to file descriptor */
      if(filedesc < 0)
   {
   return SCIP_FILECREATEERROR;
   }
   snprintf(sympath, SCIP_MAXSTRLEN, "/proc/self/fd/%d", filedesc); /* set symbolic link to file */
   success = readlink(sympath, pfile, SCIP_MAXSTRLEN); /* get actual path including extension */
   if(success < 0)
   {
   return SCIP_NOFILE;
   }
   SCIPsplitFilename(pfile, &filepath, &filename, NULL, NULL);
   strcpy(makefilename, filepath);
   strcat(makefilename, "/");
   strcat(makefilename, "Makefile");

   /* open and write first lines of makefile */
   makefile = fopen(makefilename, "w");
   if(makefile == NULL)
   {
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, makefile, "# LaTeX code might have to be compiled several times                         %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, ".PHONY: %s.pdf all clean                                                     %s", filename, LINEBREAK);
   SCIPinfoMessage(scip, makefile, "                                                                             %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "# dependencies, gp etc here                                                  %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "                                                                             %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "# latexmk automatically manages the .tex files                               %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "%s.pdf: %s.tex                                                               %s", filename, filename, LINEBREAK);
   SCIPinfoMessage(scip, makefile, "\tlatexmk -pdf -pdflatex=\"pdflatex -interaction=nonstopmode\" -use-make %s.tex %s", filename, LINEBREAK);
   SCIPinfoMessage(scip, makefile, "                                                                             %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "clean:                                                                       %s", LINEBREAK);
   SCIPinfoMessage(scip, makefile, "\tlatexmk -CA                                                               %s", LINEBREAK);


   /*@todo write into makefile */


   /* --- make the tex files --- */

   /*@todo sort decomps into sorteddecomps (just rearrange pointers)*/
   sorteddecomps = decomps;

   SCIP_CALL( writeHeaderCode(scip,file,statistics,sorteddecomps,ndecomps,toc,useGp) );

   for( i=0; i<*ndecomps; i++ )
   {
      if(decomps[i] != NULL)
      {
         SCIP_CALL( writeDecompCode(scip,file,makefile,decomps[i], useGp) );
      }
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
         "reading/texreader/path", "path to tex file",
         NULL, FALSE, PATH_DEFAULT, NULL, NULL) );

   return SCIP_OKAY;
}
