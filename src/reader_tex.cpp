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
#include <fstream>

#include "reader_tex.h"
#include "scip_misc.h"
#include "pub_gcgvar.h"
#include "reader_gp.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "struct_decomp.h"
#include "params_visu.h"

#include "class_miscvisualization.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "wrapper_seeed.h"

#define READER_NAME             "texreader"
#define READER_DESC             "LaTeX file writer for seeed visualization"
#define READER_EXTENSION        "tex"

#define DEFAULT_USEGP            FALSE
#define DEFAULT_PICTURESONLY     FALSE

namespace gcg{


/** destructor of reader to free user data (called when SCIP is exiting) */

SCIP_DECL_READERFREE(readerFreeTex)
{
   return SCIP_OKAY;
}

/** Problem reading method of reader.
 *  Since the reader is not supposed to read files this returns a reading error. */

SCIP_DECL_READERREAD(readerReadTex)
{  /*lint --e{715}*/
   return SCIP_READERROR;
}

/** problem writing method of reader */

SCIP_DECL_READERWRITE(readerWriteTex)
{
   int seeedid;

   assert(scip != NULL);
   assert(reader != NULL);

   seeedid = 0;

   /* get seeed to write */
   seeedid = *(DECgetBestSeeed(scip));

   if(seeedid == -1)
   {
      SCIPerrorMessage("Could not find best Seeed!\n");
      *result = SCIP_DIDNOTRUN;
   }
   else
   {

      GCGwriteTexVisualization(scip, file, seeedid, TRUE, FALSE);
      *result = SCIP_SUCCESS;
   }

   return SCIP_OKAY;
}


/* outputs the r, g, b decimal values for the rgb hex input */
static
SCIP_RETCODE getRgbFromHex(
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


/** converts a hex color code into a tex-conform line of code that defines the color as \colorname */
static
char* getTexColorFromHex(
   char* hex,              /* hex code for color */
   const char* colorname   /* name of color */
   )
{
   char* texcode;
   int r;
   int g;
   int b;

   getRgbFromHex( hex, &r, &g, &b );
   texcode = '\0';

   strcat( texcode, "\\definecolor{" );
   strcpy( texcode, colorname );
   strcpy( texcode, "}{RGB}{" );
   snprintf(texcode, SCIP_MAXSTRLEN, "%d", r);
   strcpy( texcode, "," );
   snprintf(texcode, SCIP_MAXSTRLEN, "%d", g);
   strcpy( texcode, "," );
   snprintf(texcode, SCIP_MAXSTRLEN, "%d", b);
   strcpy( texcode, "}" );

   return texcode;
}


/** write LaTeX code header & begin of document to given file */
static
SCIP_RETCODE writeTexHeader(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file                /**< File pointer to write to */
   )
{
   char* pname;
   char ppath[SCIP_MAXSTRLEN];

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
   SCIPinfoMessage(scip, file, "%% * Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       * \n");
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
   SCIPinfoMessage(scip, file, " \\usetikzlibrary{external}                                                      \n");
   SCIPinfoMessage(scip, file, " \\tikzexternalize                                                               \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");

  /* introduce colors of current color scheme */
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorMasterconss(),
      "colormasterconss"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorMastervars(),
      "colormastervars"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorLinking(),
      "colorlinking"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorStairlinking(),
      "colorstairlinking"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorBlock(),
      "colorblock"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorOpen(),
      "coloropen"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorNonzero(),
      "colornonzero"));
   SCIPinfoMessage(scip, file, "%s                            \n", getTexColorFromHex(SCIPvisuGetColorLine(),
      "colorline"));
   SCIPinfoMessage(scip, file, "                                                                                 \n");
   SCIPinfoMessage(scip, file, "\\begin{document}                                                                \n");
   SCIPinfoMessage(scip, file, "                                                                                 \n");

   return SCIP_OKAY;
}


/** write LaTeX code title page that includes general statistics about the problem to given file */
static
SCIP_RETCODE writeTexTitlepage(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File pointer to write to */
   int*                 npresentedseeeds    /**< Number of decompositions to be shown in the file or NULL if unknown */
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
   SCIPinfoMessage(scip, file, "  Number of found finished decompositions: & %i  \\\\                         \n",
      SCIPconshdlrDecompGetNFinishedDecomps(scip));
   SCIPinfoMessage(scip, file, "  Number of found incomplete decompositions: & %i  \\\\                       \n",
      SCIPconshdlrDecompGetNSeeeds(scip) - SCIPconshdlrDecompGetNFinishedDecomps(scip));
   if(npresentedseeeds != NULL){
      if( ndecomps > *npresentedseeeds )
      {
         SCIPinfoMessage(scip, file, "  Number of decompositions presented in this document: & %i \\\\ \n",
            *npresentedseeeds);
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


/** write LaTeX code for table of contents to given file */
static
SCIP_RETCODE writeTexTableOfContents(
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


/** writes line to given file that contains the tikz code for a box with given dimensions
 * and options with a linebreak at the end. */
static
SCIP_RETCODE writeTikzBox(
   SCIP* scip,       /**< SCIP data structure */
   FILE* file,       /**< File pointer to write to */
   int xmax,         /**< maximum x axis value */
   int ymax,         /**< maximum y axis value */
   int x1,           /**< x value of lower left vertex coordinate */
   int y1,           /**< y value of lower left vertex coordinate */
   int x2,           /**< x value of upper right vertex coordinate */
   int y2,           /**< y value of upper right vertex coordinate */
   const char* color /**< color name */
   )
{
   SCIPinfoMessage(scip, file,
      "    \\draw [fill=%s] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) rectangle (%f*\\textwidth*0.75,%f*\\textwidth*0.75);\n",
      color, x1 / xmax, y1 / ymax, x2 / xmax, y2 / ymax);
   return SCIP_OKAY;
}


/*@todo adapt this for tikz, currently this gives gp code */
/** writes line to given file that contains the tikz code for a point with given radius
 * and options with a linebreak at the end. */
static
SCIP_RETCODE writeTikzNonzeros(
   SCIP* scip,             /**< SCIP data structure */
   FILE* file,             /**< filename to write to (including path & extension) */
   Seeed* seeed,           /**< Seeed for which the nonzeros should be visualized */
   Seeedpool* seeedpool,   /**< current Seeedpool */
   float radius            /**< radius of the dots */
   )
{
   std::vector<int> orderToRows(seeed->getNConss(), -1);
   std::vector<int> rowToOrder(seeed->getNConss(), -1);
   std::vector<int> orderToCols(seeed->getNVars(), -1);
   std::vector<int> colsToOrder(seeed->getNVars(), -1);
   int counterrows = 0;
   int countercols = 0;

   /** order of constraints */
   /* master constraints */
   for( int i = 0; i < seeed->getNMasterconss() ; ++i )
   {
      int rowidx = seeed->getMasterconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /* block constraints */
   for( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for(int i = 0; i < seeed->getNConssForBlock(b) ; ++i )
      {
         int rowidx = seeed->getConssForBlock(b)[i];
         orderToRows[counterrows] = rowidx;
         rowToOrder[rowidx] = counterrows;
         ++counterrows;
      }
   }

   /** open constraints */
   for( int i = 0; i < seeed->getNOpenconss() ; ++i )
   {
      int rowidx = seeed->getOpenconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /** order of variables */

   /* linking variables */
   for( int i = 0; i < seeed->getNLinkingvars() ; ++i )
   {
      int colidx = seeed->getLinkingvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* master variables */
   for( int i = 0; i < seeed->getNMastervars() ; ++i )
   {
      int colidx = seeed->getMastervars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* block variables */
   for( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for(int i = 0; i < seeed->getNVarsForBlock(b) ; ++i )
      {
         int colidx = seeed->getVarsForBlock(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
      for(int i = 0; i < seeed->getNStairlinkingvars(b) ; ++i )
      {
         int colidx = seeed->getStairlinkingvars(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
   }

   /* open vars */
   for( int i = 0; i < seeed->getNOpenvars() ; ++i )
   {
      int colidx = seeed->getOpenvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* write scatter plot */
   for( int row = 0; row < seeed->getNConss(); ++row )
   {
      for ( int col = 0; col < seeed->getNVars(); ++col )
      {
         assert( orderToRows[row] != -1 );
         assert( orderToCols[col] != -1 );
         if( seeedpool->getVal( orderToRows[row], orderToCols[col]  ) != 0 )
         {
            SCIPinfoMessage(scip, file,
               "    \\draw [fill] (%f*\\textwidth*0.75,%f*\\textwidth*0.75) circle [radius=%f*0.75];\n",
               col + 0.5, row + 0.5, radius);
         }
      }
   }

   return SCIP_OKAY;
}


static
SCIP_RETCODE writeTexSeeed(
   SCIP* scip,             /**< SCIP data structure */
   FILE* file,             /**< filename (including path) to write to */
   Seeed* seeed,           /**< Seeed for which the nonzeros should be visualized */
   Seeedpool* seeedpool    /**< current Seeedpool */
   )
{
   int rowboxcounter = 0;
   int colboxcounter = 0;
   int nvars;
   int nconss;

   nvars = seeed->getNVars();
   nconss = seeed->getNConss();

   SCIPinfoMessage(scip, file, "\\begin{figure}[!htb]                                              \n");
   SCIPinfoMessage(scip, file, "  \\begin{center}                                                  \n");
   SCIPinfoMessage(scip, file, "  \\begin{tikzpicture}                                             \n");

   /* --- draw boxes ---*/

   /* linking vars */
   writeTikzBox(scip, file, nvars, nconss, 0, 0, seeed->getNLinkingvars(), seeed->getNConss(),
      (const char*) "colorlinking");
   colboxcounter += seeed->getNLinkingvars();

   /* mastervars */
   writeTikzBox(scip, file, nvars, nconss, colboxcounter, 0, seeed->getNMastervars()+colboxcounter, seeed->getNConss(),
      (const char*) "colormastervars");
   colboxcounter += seeed->getNMastervars();

   /* masterconss */
   writeTikzBox(scip, file, nvars, nconss, 0, 0, seeed->getNVars(), seeed->getNMasterconss(),
      (const char*) "colormasterconss");
   rowboxcounter += seeed->getNMasterconss();

   /* blocks */
   for( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      writeTikzBox(scip, file, nvars, nconss, colboxcounter, rowboxcounter,
         colboxcounter + seeed->getNVarsForBlock(b), rowboxcounter + seeed->getNConssForBlock(b),
         (const char*) "colorblock");
      colboxcounter += seeed->getNVarsForBlock(b);

      if( seeed->getNStairlinkingvars(b) != 0 )
      {
         writeTikzBox(scip, file, nvars, nconss, colboxcounter, rowboxcounter,
            colboxcounter + seeed->getNStairlinkingvars(b),
            rowboxcounter + seeed->getNConssForBlock(b) + seeed->getNConssForBlock(b+1),
            (const char*) "colorstairlinking");
      }
      colboxcounter += seeed->getNStairlinkingvars(b);
      rowboxcounter += seeed->getNConssForBlock(b);
   }

   /* open */
   writeTikzBox(scip, file, nvars, nconss, colboxcounter, rowboxcounter, colboxcounter + seeed->getNOpenvars(),
      rowboxcounter+seeed->getNOpenconss(), (const char*) "coloropen" );
   colboxcounter += seeed->getNOpenvars();
   rowboxcounter += seeed->getNOpenconss();

   /* --- draw nonzeros --- */
   writeTikzNonzeros( scip, file, seeed, seeedpool, SCIPvisuGetNonzeroRadius( seeed->getNVars(), seeed->getNConss(), 1 ) );

   SCIPinfoMessage(scip, file, "  \\end{tikzpicture}                                               \n");
   SCIPinfoMessage(scip, file, "  \\end{center}                                                    \n");
   SCIPinfoMessage(scip, file, "\\end {figure}                                                     \n");

   return SCIP_OKAY;
}

static
SCIP_RETCODE writeTexSeeedStatistics(
   SCIP* scip,             /**< SCIP data structure */
   FILE* file,             /**< filename (including path) to write to */
   Seeed* seeed            /**< Seeed for which the nonzeros should be visualized */
   )
{
   DEC_DETECTOR** detectorchain;
   char fulldetectorstring[SCIP_MAXSTRLEN];
   int sizedetectorchain;
   int i;

   /* get detector chain full-text string*/
   detectorchain = seeed->getDetectorchain();
   sizedetectorchain = seeed->getNDetectors();
   if( detectorchain[0] != NULL)
      sprintf(fulldetectorstring, "%s", DECdetectorGetName(detectorchain[0]));
   else
      sprintf(fulldetectorstring, "%s", "user");
   for( i=1; i < sizedetectorchain; ++i )
   {
      sprintf(fulldetectorstring, "%s, %s",fulldetectorstring, DECdetectorGetName(detectorchain[i]) );
   }

   SCIPinfoMessage(scip, file, "                                                                \n");
   SCIPinfoMessage(scip, file, "\\vspace{0.3cm}                                                 \n");
   SCIPinfoMessage(scip, file, "\\begin{tabular}{lp{10cm}}                                      \n");
   SCIPinfoMessage(scip, file,
      "  Found by detector(s): & \\begin{minipage}{10cm}\\begin{verbatim}%s\\end{verbatim}\\end{minipage} \\\\ \n",
      fulldetectorstring);
   SCIPinfoMessage(scip, file, "  Number of blocks: & %i \\\\                                                   \n",
      seeed->getNBlocks());
   SCIPinfoMessage(scip, file, "  Number of master variables: & %i \\\\                                         \n",
      seeed->getNMastervars());
   SCIPinfoMessage(scip, file, "  Number of master constraints: & %i \\\\                                       \n",
      seeed->getNMasterconss());
   SCIPinfoMessage(scip, file, "  Number of linking variables: & %i \\\\                                        \n",
      seeed->getNLinkingvars());
   SCIPinfoMessage(scip, file, "  Number of stairlinking variables: & %i \\\\                                   \n",
      seeed->getNTotalStairlinkingvars());
   SCIPinfoMessage(scip, file, "  Max white score: & %f \\\\                                                    \n",
      seeed->getMaxWhiteScore());
   SCIPinfoMessage(scip, file, "\\end{tabular}                                                                  \n");

   SCIPinfoMessage(scip, file, "\\clearpage                                                                     \n");
   SCIPinfoMessage(scip, file, "                                                                                \n");

   return SCIP_OKAY;
}


/** write LaTeX code for end of document to given file */
static
SCIP_RETCODE writeTexEnding(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{

   SCIPinfoMessage(scip, file, "\\end{document}                                                                  \n");

   return SCIP_OKAY;
}

/** writes a visualization for the given seeed */
SCIP_RETCODE GCGwriteTexVisualization(
   SCIP* scip,             /**< SCIP data structure */
   FILE* file,             /**< filename including path */
   int seeedid,            /**< id of seeed to visualize */
   SCIP_Bool statistics,   /**< additionally to picture show statistics */
   SCIP_Bool usegp         /**< true if the gp reader should be used to visualize the individual seeeds */
   )
{
   MiscVisualization* misc = new MiscVisualization();
   Seeed* seeed;
   Seeedpool* seeedpool = NULL;
   char* gpname;
   char* pdfname;

   /* get seeed */
   seeed = misc->GCGgetSeeedWithPool(scip, seeedid, seeedpool);

   /* write tex code into file */
   writeTexHeader(scip, file);

   if(!usegp)
   {
      writeTexSeeed(scip, file, seeed, seeedpool);
   }
   else
   {
      /* in case a gp file should be generated include it */
      gpname = misc->GCGgetVisualizationFilename(scip, seeed, "gp");
      pdfname = misc->GCGgetVisualizationFilename(scip, seeed, "pdf");

      GCGwriteGpVisualization(scip, gpname, pdfname, seeedid);

      SCIPinfoMessage(scip, file, "\\begin{figure}[!htb]                                              \n");
      SCIPinfoMessage(scip, file, "  \\begin{center}                                                  \n");
      SCIPinfoMessage(scip, file, "    \\input{%s}                                           \n", pdfname);
      SCIPinfoMessage(scip, file, "  \\end{center}                                                    \n");
      SCIPinfoMessage(scip, file, "\\end {figure}                                                     \n");
   }
   if(statistics)
      writeTexSeeedStatistics(scip, file, seeed);
   writeTexEnding(scip, file);

   return SCIP_OKAY;
}


/** writes a report for the given seeeds */
SCIP_RETCODE GCGwriteTexReport(
   SCIP* scip,             /**< SCIP data structure */
   FILE* file,             /**< filename including path */
   int* seeedids,          /**< ids of seeeds to visualize */
   int* nseeeds,            /**< number of seeeds to visualize */
   SCIP_Bool titlepage,    /**< true if a title page should be included in the document */
   SCIP_Bool toc,          /**< true if an interactive table of contents should be included */
   SCIP_Bool statistics,   /**< true if statistics for each seeed should be included */
   SCIP_Bool usegp         /**< true if the gp reader should be used to visualize the individual seeeds */
   )
{
   MiscVisualization* misc = new MiscVisualization();
   Seeed* seeed;
   Seeedpool* seeedpool = NULL;
   char* gpname;
   char* pdfname;

   /* write tex code into file */
   writeTexHeader(scip, file);
   if(titlepage)
      writeTexTitlepage(scip, file, nseeeds);
   if(toc)
      writeTexTableOfContents(scip, file);
   for(int i = 0; i < *nseeeds; i++)
   {
      if(toc)
      {
         char decompname[SCIP_MAXSTRLEN];

         SCIPinfoMessage(scip, file, "\\section*{Decomposition: %s}                                   \n", decompname);
         SCIPinfoMessage(scip, file, "\\addcontentsline{toc}{section}{Decomposition: %s}              \n", decompname);
         SCIPinfoMessage(scip, file, "                                                                \n");
      }
      /* get and write each seeed */
      seeed = misc->GCGgetSeeedWithPool(scip, seeedids[i], seeedpool);
      if(!usegp)
      {
         writeTexSeeed(scip, file, seeed, seeedpool);
      }
      else
      {
         /* in case a gp file should be generated include it */
         gpname = misc->GCGgetVisualizationFilename(scip, seeed, "gp");
         pdfname = misc->GCGgetVisualizationFilename(scip, seeed, "pdf");

         GCGwriteGpVisualization(scip, gpname, pdfname, seeedids[i]);

         SCIPinfoMessage(scip, file, "\\begin{figure}[!htb]                                              \n");
         SCIPinfoMessage(scip, file, "  \\begin{center}                                                  \n");
         SCIPinfoMessage(scip, file, "    \\input{%s}                                           \n", pdfname);
         SCIPinfoMessage(scip, file, "  \\end{center}                                                    \n");
         SCIPinfoMessage(scip, file, "\\end {figure}                                                     \n");
      }
      if(statistics)
         writeTexSeeedStatistics(scip, file, seeed);
   }
   writeTexEnding(scip, file);

   return SCIP_OKAY;
}


/** writes a visualization of the family tree of the current seeedpool */
SCIP_RETCODE GCGwriteTexFamilyTree(
   SCIP* scip,       /**< SCIP data structure */
   FILE* file,       /**< filename including path */
   SCIP_Bool usegp   /**< true if the gp reader should be used to visualize the individual seeeds */
   )
{
   //   char* detectorchainstring;
   //   DEC_DETECTOR** detectorchain;
   //   int sizedetectorchain;
   //
   //   detectorchainstring = seeed->getDetectorChainString();
   //   detectorchain = seeed->getDetectorchain();
   //   sizedetectorchain = seeed->getNDetectors();

   return SCIP_OKAY;
}

/** makes a new makefile and readme for the given .tex file */
SCIP_RETCODE GCGtexWriteMakefileAndReadme(
   SCIP*                scip,               /**< SCIP data structure */
   FILE*                file,               /**< File for which the makefile & readme are generated */
   SCIP_Bool            usegp               /**< true if there are gp files to be included in the makefile */
   )
{
   FILE* makefile;
   FILE* readme;
   char* filepath;
   char* filename;
   char* pfile;
   char pfilecpy[SCIP_MAXSTRLEN];
   char makefilename[SCIP_MAXSTRLEN];
   char readmename[SCIP_MAXSTRLEN];
   char name[SCIP_MAXSTRLEN];
   const char makename[SCIP_MAXSTRLEN] = "makepdf";

   /* --- create a Makefile --- */

   /* get path to write to and put it into makefilename */
   gcg::MiscVisualization* miscvisu = new gcg::MiscVisualization();
   pfile = miscvisu->GCGgetFilePath(scip, file);
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

   if( usegp )
   {
      SCIPinfoMessage(scip, makefile, "GPFILES := $(wildcard *.gp)\n");
   }
   SCIPinfoMessage(scip, makefile, "                                                                             \n");
   SCIPinfoMessage(scip, makefile, "# latexmk automatically manages the .tex files                               \n");
   SCIPinfoMessage(scip, makefile, "%s.pdf: %s.tex\n",
      filename, filename);
   if( usegp )
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
   if( usegp )
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
   /* include tex reader */
   SCIP_CALL(SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION, NULL,
           readerFreeTex, readerReadTex, readerWriteTex, NULL));

   return SCIP_OKAY;
}

} /* namespace gcg */
