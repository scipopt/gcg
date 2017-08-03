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

/**@file   reader_gp.cpp
 * @brief  GP file reader writing seeeds to gnuplot files
 * @author Martin Bergner
 * @author Hanna Franzen
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include <fstream>

#include "reader_gp.h"
#include "scip_misc.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "pub_decomp.h"
#include "params_visu.h"

#include "class_seeed.h"
#include "class_seeedpool.h"

#define READER_NAME             "gpreader"
#define READER_DESC             "gnuplot file writer for seeed visualization"
#define READER_EXTENSION        "gp"

#define READERGP_GNUPLOT_BOXTEMPLATE(i, x1, y1, x2, y2) "set object %d rect from %.1f,%.1f to %.1f,%.1f fc rgb \"grey\"\n", (i), (x1), (y1), (x2), (y2)
#define READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, x1, y1, x2, y2, color) "set object %d rect from %.1f,%.1f to %.1f,%.1f fc rgb \"%s\"\n", (i), (x1), (y1), (x2), (y2), (color)
#define READERGP_GNUPLOT_HEADER(outputname) "set terminal pdf\nset output \"%s.pdf\"\n", (outputname)
#define READERGP_GNUPLOT_RANGES(xmax, ymax) "set xrange [-1:%d]\nset yrange[%d:-1]\n", (xmax), (ymax)
#define READERGP_GNUPLOT_PLOTCMD "plot \"-\" using 1:2:3 notitle with circles fc rgb \"%s\" fill solid\n", DEFAULT_COLOR_NONZERO

#define READERGP_GNUPLOT_HEADER_TEX(outputname) "set terminal tikz\nset output \"%s.tex\"\nunset xtics\nunset ytics\nunset border\nunset key\nset style fill solid 1.0 noborder\nset size ratio -1\n", (outputname)

#define DEFAULT_DRAFTMODE FALSE

namespace gcg{

/** data for dec reader */
struct SCIP_ReaderData
{
   SCIP_Bool       draftmode;       /** if true shows no non-zeroes, recommended if too slow or too memory-intensive */
};

/*
 * Local methods
 */

/** writes gp code to given file that contains all nonzero points */
SCIP_RETCODE writeGpNonzeros(
   Seeed* seeed,           /**< Seeed for which the nonzeros should be visualized */
   Seeedpool* seeedpool,   /**< current Seeedpool */
   const char* filename    /**< filename to write to (including path & extension) */
   )
{
   std::vector<int> orderToRows(seeed->getNConss(), -1);
   std::vector<int> rowToOrder(seeed->getNConss(), -1);
   std::vector<int> orderToCols(seeed->getNVars(), -1);
   std::vector<int> colsToOrder(seeed->getNVars(), -1);
   int counterrows = 0;
   int countercols = 0;
   std::ofstream ofs;

   ofs.open (filename, std::ofstream::out );

   /** order of constraints */
   /* master constraints */
   for ( int i = 0; i < seeed->getNMasterconss() ; ++i )
   {
      int rowidx = seeed->getMasterconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /* block constraints */
   for ( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for (int i = 0; i < seeed->getNConssForBlock(b) ; ++i )
      {
         int rowidx = seeed->getConssForBlock(b)[i];
         orderToRows[counterrows] = rowidx;
         rowToOrder[rowidx] = counterrows;
         ++counterrows;
      }
   }

   /** open constraints */
   for ( int i = 0; i < seeed->getNOpenconss() ; ++i )
   {
      int rowidx = seeed->getOpenconss()[i];
      orderToRows[counterrows] = rowidx;
      rowToOrder[rowidx] = counterrows;
      ++counterrows;
   }

   /** order of variables */

   /* linking variables */
   for ( int i = 0; i < seeed->getNLinkingvars() ; ++i )
   {
      int colidx = seeed->getLinkingvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* master variables */
   for ( int i = 0; i < seeed->getNMastervars() ; ++i )
   {
      int colidx = seeed->getMastervars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* block variables */
   for ( int b = 0; b < seeed->getNBlocks() ; ++b )
   {
      for (int i = 0; i < seeed->getNVarsForBlock(b) ; ++i )
      {
         int colidx = seeed->getVarsForBlock(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
      for (int i = 0; i < seeed->getNStairlinkingvars(b) ; ++i )
      {
         int colidx = seeed->getStairlinkingvars(b)[i];
         orderToCols[countercols] = colidx;
         colsToOrder[colidx] = countercols;
         ++countercols;
      }
   }

   /** open vars */
   for ( int i = 0; i < seeed->getNOpenvars() ; ++i )
   {
      int colidx = seeed->getOpenvars()[i];
      orderToCols[countercols] = colidx;
      colsToOrder[colidx] = countercols;
      ++countercols;
   }

   /* write scatter plot */
   for( int row = 0; row < seeed->getNConss(); ++row )
      for ( int col = 0; col < seeed->getNVars(); ++col )
      {
         assert( orderToRows[row] != -1);
         assert( orderToCols[col] != -1);
         if( seeedpool->getVal( orderToRows[row], orderToCols[col]  ) != 0 )
            ofs << col+0.5 << " " << row+0.5 << std::endl;
      }

   ofs.close();

   return SCIP_OKAY;
}


/** write file header with terminal etc. */
static
SCIP_RETCODE writeFileHeader(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   const char*           outname,            /**< the name of the gnuplot outputname */
   SCIP_Bool             outputPDF           /**< if true give pdf file, if false give tex file instead */
   )
{

   if(outputPDF)
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER(outname));
   else
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER_TEX(outname));

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_RANGES(SCIPgetNVars(scip), SCIPgetNConss(scip)));
   return SCIP_OKAY;
}

/** write decomposition header such as rectangles for blocks etc. */
static
SCIP_RETCODE writeDecompositionHeader(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp           /**< Decomposition pointer */
   )
{
   int i;
   int b;
   int startx;
   int starty;
   int endx;
   int endy;
   int nstairlinkingvars;
   int nmastervars;
   int nvars;
   int nconss;
   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);
   if( decdecomp->type == DEC_DECTYPE_UNKNOWN || decdecomp->nblocks == 0 )
   {
      return SCIP_OKAY;
   }

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   nstairlinkingvars = 0;
   if (decdecomp->nstairlinkingvars != NULL)
   {
      for( b = 0; b < decdecomp->nblocks - 1; ++b )
      {
         nstairlinkingvars += decdecomp->nstairlinkingvars[b];
      }
   }
   nmastervars = 0;
   for( i = 0; i < decdecomp->nlinkingvars; ++i )
   {
      if( (int)(size_t)SCIPhashmapGetImage(decdecomp->vartoblock, decdecomp->linkingvars[i]) == decdecomp->nblocks + 1)
         nmastervars++;
   }
   startx = 0;
   starty = 0;
   i = 1;
   /** write linking var box */
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, startx + 0.5, starty + 0.5, decdecomp->nlinkingvars - nstairlinkingvars - nmastervars + 0.5, nconss + 0.5, DEFAULT_COLOR_LINKING));
   i++;
   startx += decdecomp->nlinkingvars - nstairlinkingvars - nmastervars;
   /** write master var box */
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, startx + 0.5, starty + 0.5, startx + nmastervars + 0.5, nconss + 0.5, DEFAULT_COLOR_MASTERVARS));
   i++;
  startx += nmastervars;
  /** write linking cons box */
  SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, 0 + 0.5, 0 + 0.5, nvars + 0.5, decdecomp->nlinkingconss + 0.5, DEFAULT_COLOR_MASTERCONS));
  i++;
  starty += decdecomp->nlinkingconss;

  endx = startx;
  endy = starty;
  for( b = 0; b < decdecomp->nblocks; ++b )
  {
     endx += decdecomp->nsubscipvars[b];
     endy += decdecomp->nsubscipconss[b];
     SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, startx + 0.5, starty + 0.5, endx + 0.5, endy + 0.5, DEFAULT_COLOR_BLOCK));
     i++;
     if(decdecomp->nstairlinkingvars != NULL )
     {
        if(decdecomp->nstairlinkingvars[b] != 0)
        {
        startx = endx;
        endx += decdecomp->nstairlinkingvars[b];
        SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATECOLORED(i, startx + 0.5, starty + 0.5, endx + 0.5, starty + decdecomp->nsubscipconss[b] + decdecomp->nsubscipconss[b+1] + 0.5, DEFAULT_COLOR_STAIRLINKING));
        i++;
        }
     }
     startx = endx;
     starty = endy;
  }

   return SCIP_OKAY;
}

/** write the plot commands */
static
SCIP_RETCODE writePlotCommands(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   assert(scip != NULL);
   assert(file != NULL);

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_PLOTCMD);
   return SCIP_OKAY;
}

/** write the data optionally using the decomposition data */
static
SCIP_RETCODE writeData(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp           /**< Decomposition pointer */
   )
{
   SCIP_CONS** conss;
   SCIP_HASHMAP* varindexmap;
   SCIP_HASHMAP* consindexmap;
   int* stairlinkingvars;
   int nconss;

   int i;
   int j;

   assert(scip != NULL);
   assert(file != NULL);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   varindexmap = NULL;
   consindexmap = NULL;

   if( decdecomp != NULL )
   {
      size_t varindex = 1;
      size_t consindex = 1;
      SCIP_CALL( SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
      SCIP_CALL( SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );

      /** fillout the array stairlinkingvars */
      SCIP_CALL( SCIPallocBufferArray(scip, &stairlinkingvars, SCIPgetNVars(scip)) );
      for( j = 0; j < SCIPgetNVars(scip); ++j)
         stairlinkingvars[j] = 0;
      for( i = 0; i < decdecomp->nblocks - 1; ++i)
      {
         if(decdecomp->nstairlinkingvars != NULL)
         {
            for( j = 0; j < decdecomp->nstairlinkingvars[i]; ++j)
            {
               assert(SCIPhashmapGetImage(decdecomp->varindex, decdecomp->stairlinkingvars[i][j]) != NULL);
               assert(stairlinkingvars[(int)(size_t)SCIPhashmapGetImage(decdecomp->varindex, decdecomp->stairlinkingvars[i][j])] == 0);
               stairlinkingvars[(int)(size_t)SCIPhashmapGetImage(decdecomp->varindex, decdecomp->stairlinkingvars[i][j])] = 1;
            }
         }
      }

      /** fillout the hashmaps */
      for( j = 0; j < decdecomp->nlinkingconss; ++j )
      {
         assert(decdecomp->linkingconss[j] != NULL);
         SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->linkingconss[j], (void*)consindex) );
         consindex++;
      }
      /** linkingvars */
      for( j = 0; j < decdecomp->nlinkingvars; ++j )
      {
         assert(decdecomp->linkingvars[j] != NULL);
         assert((int)(size_t)SCIPhashmapGetImage(decdecomp->vartoblock, decdecomp->linkingvars[j]) == decdecomp->nblocks + 2 || (int)(size_t)SCIPhashmapGetImage(decdecomp->vartoblock, decdecomp->linkingvars[j]) == decdecomp->nblocks + 1);
         if( (int)(size_t)SCIPhashmapGetImage(decdecomp->vartoblock, decdecomp->linkingvars[j]) == decdecomp->nblocks + 2 && stairlinkingvars[(int)(size_t)SCIPhashmapGetImage(decdecomp->varindex, decdecomp->linkingvars[j])] == 0)
         {
            SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->linkingvars[j], (void*)varindex) );
            varindex++;
         }
      }
      /** mastervars */
      for( j = 0; j < decdecomp->nlinkingvars; ++j )
      {
         assert(decdecomp->linkingvars[j] != NULL);
         if( (int)(size_t)SCIPhashmapGetImage(decdecomp->vartoblock, decdecomp->linkingvars[j]) == decdecomp->nblocks + 1 )
         {
            SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->linkingvars[j], (void*)varindex) );
            varindex++;
         }
      }
      SCIPdebugMessage("Block information:\n");

      for( i = 0; i < decdecomp->nblocks; ++i )
      {
         SCIPdebugPrintf("Block %d:\n", i+1);
         SCIPdebugPrintf("\tVars: %d", decdecomp->nsubscipvars[i]);
         SCIPdebugPrintf("\tConss: %d\n", decdecomp->nsubscipconss[i]);
         for( j = 0; j < decdecomp->nsubscipvars[i]; ++j )
         {
            assert(decdecomp->subscipvars[i][j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->subscipvars[i][j], (void*)varindex) );
            varindex++;
         }
         if( decdecomp->nstairlinkingvars != 0)
         {
            for( j = 0; j < decdecomp->nstairlinkingvars[i]; ++j )
            {
               assert(decdecomp->stairlinkingvars[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->stairlinkingvars[i][j], (void*)varindex) );
               varindex++;
            }
         }
         for( j = 0; j < decdecomp->nsubscipconss[i]; ++j )
         {
            assert(decdecomp->subscipconss[i][j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->subscipconss[i][j], (void*)consindex) );
            consindex++;
         }
      }
      SCIPfreeBufferArray(scip, &stairlinkingvars);
   }

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
         assert(curvars != NULL);

         /* if the problem has been created, output the whole model */
         if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
         {
            SCIPinfoMessage(scip, file, "%d %d 0.25\n", SCIPvarGetIndex(curvars[j]), i);
            continue;
         }

         /* if there is no decomposition, output the presolved model! */
         if( decdecomp == NULL || decdecomp->type == DEC_DECTYPE_UNKNOWN )
         {
            SCIPinfoMessage(scip, file, "%d %d 0.25\n", SCIPvarGetIndex(curvars[j]), i);
         }
         /* if there is a decomposition, output the indices derived from the decomposition above*/
         else
         {
            assert(varindexmap != NULL);
            assert(consindexmap != NULL);
            assert(SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])) != NULL);
            assert(SCIPhashmapGetImage(consindexmap, conss[i]) != NULL);

            SCIPinfoMessage(scip, file, "%d %d 0.25\n",
               SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])),
               SCIPhashmapGetImage(consindexmap, conss[i])
               );
         }
      }

      SCIPfreeBufferArrayNull(scip, &curvars);
   }

   if( decdecomp != NULL && decdecomp->type != DEC_DECTYPE_STAIRCASE )
   {
      SCIPhashmapFree(&varindexmap);
      SCIPhashmapFree(&consindexmap);
   }

   return SCIP_OKAY;
}


/** write trailer of the file */
static
SCIP_RETCODE writeFileTrailer(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< File pointer to write to */
   )
{
   SCIPinfoMessage(scip, file, "e\n");
   return SCIP_OKAY;
}


/*
 * Callback methods of reader
 */

#define readerCopyGp NULL

/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeGp)
{
   SCIP_READERDATA* readerdata;

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIPfreeMemory(scip, &readerdata);

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   return SCIP_OKAY;
}


/** problem reading method of reader */
#define readerReadGp NULL



/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGp)
{
   /*lint --e{715}*/
   assert(scip != NULL);

   SCIP_CALL( SCIPwriteGp(scip, file, DECgetBestDecomp(scip), TRUE, TRUE) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** writes the decomposition to the specific file */
SCIP_RETCODE SCIPwriteGp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition, /**< whether to write decomposed problem */
   SCIP_Bool             outputPDF           /**< if true give pdf file, if false give tex file instead */
   )
{
   char* name;
   char* detectorchainstring;
   char probname[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];

   SCIP_READERDATA* readerdata;
   readerdata = SCIPreaderGetData(SCIPfindReader(scip, "gpreader"));

   assert(scip != NULL);
   assert(file != NULL);

   if( writeDecomposition && decdecomp == NULL )
   {
      SCIPwarningMessage(scip, "Cannot write decomposed problem if decomposition structure empty!");
      writeDecomposition = FALSE;
   }
   /* sanitize filename */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);

   /* get detector chain string*/
   detectorchainstring = DECdecompGetDetectorChainString(scip, decdecomp);

   /* print header */
   if( decdecomp == NULL )
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", name);
   else
   {
      if(outputPDF)
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%s_%d_%d", name, detectorchainstring, DECdecompGetSeeedID(decdecomp),
            decdecomp->nblocks);
      else
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s-%s-%d-%d", name, detectorchainstring, DECdecompGetSeeedID(decdecomp),
            decdecomp->nblocks);
   }

   SCIP_CALL( writeFileHeader(scip, file, outname, outputPDF) );

   /* write decomp information such as rectangles */
   if( writeDecomposition )
      SCIP_CALL( writeDecompositionHeader(scip, file, decdecomp) );

   /* write the plot header*/
   SCIP_CALL( writePlotCommands(scip, file) );

   /* write data (if draftmode is not on) */
   if(!readerdata->draftmode){
      SCIP_CALL( writeData(scip, file, decdecomp) );
   }

   /* write file end */
   SCIP_CALL( writeFileTrailer(scip, file) );
   return SCIP_OKAY;
}


/** includes the gp file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create gp reader data */
   SCIP_CALL( SCIPallocMemory(scip, &readerdata) );

   /* include gp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
      readerCopyGp, readerFreeGp, readerReadGp, readerWriteGp, readerdata) );

   SCIP_CALL( SCIPaddBoolParam(scip,
      "reading/gpreader/draftmode",
      "if true shows no non-zeroes, recommended if too slow or too memory-intensive",
      &readerdata->draftmode, FALSE, DEFAULT_DRAFTMODE, NULL, NULL) );

   return SCIP_OKAY;
}

} /* namespace gcg */
