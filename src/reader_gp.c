/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: reader_gp.c,v 1.21 2010/01/04 20:35:47 bzfheinz Exp $"

/**@file   reader_gp.c
 * @ingroup FILEREADERS 
 * @brief  GP file reader
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>
#include <string.h>

#include "reader_gp.h"
#include "scip_misc.h"

#define READER_NAME             "gpreader"
#define READER_DESC             "gnuplot file writer for matrix visualization"
#define READER_EXTENSION        "gp"

#define READERGP_GNUPLOT_BOXTEMPLATE(i, x1, y1, x2, y2) "set object %d rect from %.1f,%.1f to %.1f,%.1f fc rgb \"grey\"\n", (i), (x1), (y1), (x2), (y2)
#define READERGP_GNUPLOT_HEADER(outputname) "set terminal pdf\nset output \"%s.pdf\"\nunset xtics\nunset ytics\nunset border\n", (outputname)
#define READERGP_GNUPLOT_RANGES(xmax, ymax) "set xrange [0:%d]\nset yrange[%d:0]\n", (xmax), (ymax)
//#define READERGP_GNUPLOT_HEADER(outputname) "set terminal pdf\nset output \"%s.pdf\"\nunset border\n", (outputname)

#define READERGP_GNUPLOT_PLOTCMD "plot \"-\" lt 0 pt 5 ps 0.5 notitle\n"
/*
 * Data structures
 */

/** data for gp reader */
struct SCIP_ReaderData
{
   DECDECOMP* decdecomp;
   SCIP_HASHMAP *vartoindex;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** write file header with terminal etc. */
static
SCIP_RETCODE writeFileHeader(
   SCIP*        scip,                           /**< SCIP data structure */
   FILE*        file,                           /**< File pointer to write to */
   const char*  outname                         /**< the name of the gnuplot outputname */
   )
{

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER(outname));
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_RANGES(SCIPgetNVars(scip), SCIPgetNConss(scip)));
   return SCIP_OKAY;
}

/** write decomposition header such as rectangles for blocks etc. */
static
SCIP_RETCODE writeDecompositionHeader(
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file,                                /**< File pointer to write to */
   DECDECOMP* decdecomp                      /**< Decomposition pointer */
   )
{

   int i;
   int startx;
   int starty;
   int endx;
   int endy;
   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);
   if(decdecomp == NULL || decdecomp->type == DEC_UNKNOWN || decdecomp->nblocks == 0)
   {
      return SCIP_OKAY;
   }

   if(decdecomp->type == DEC_ARROWHEAD || decdecomp->type == DEC_BORDERED)
   {
      startx = 0;
      starty = 0;
      endx = 0;
      endy = 0;

      for( i = 0; i < decdecomp->nblocks; ++i)
      {
         endx += decdecomp->nsubscipvars[i];
         endy += decdecomp->nsubscipconss[i];
         SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+1, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
         startx = endx;
         starty = endy;
      }
      endx += decdecomp->nlinkingvars;
      endy += decdecomp->nlinkingconss;
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+2, 0.5, starty+0.5, endx+0.5, endy+0.5));
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+3, startx+0.5, +0.5, endx+0.5, endy+0.5));
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+4, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
   }
   return SCIP_OKAY;
}

/** write the plot commands */
static
SCIP_RETCODE writePlotCommands(
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file                                 /**< File pointer to write to */
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
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file,                                /**< File pointer to write to */
   DECDECOMP* decdecomp                      /**< Decomposition pointer */

   )
{
   SCIP_VAR** vars;
   SCIP_CONS** conss;

   SCIP_HASHMAP* varindexmap;
   SCIP_HASHMAP* consindexmap;
   int nvars;
   int nconss;

   int i;
   int j;
   size_t varindex;
   size_t consindex;

   assert(scip != NULL);
   assert(file != NULL);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL(SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)));
   SCIP_CALL(SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)));

   if(decdecomp != NULL)
   {
      /* if we don't have staicase, but something else, go through the blocks and create the indices */
      if(decdecomp->type != DEC_STAIRCASE )
      {
         SCIPdebugMessage("Block information:\n");
         varindex = 1;
         consindex = 1;
         for( i = 0; i < decdecomp->nblocks; ++i)
         {
            SCIPdebugPrintf("Block %d:\n", i+1);
            SCIPdebugPrintf("\tVars: %d", decdecomp->nsubscipvars[i]);
            SCIPdebugPrintf("\tConss: %d\n", decdecomp->nsubscipconss[i]);
            for( j = 0; j < decdecomp->nsubscipvars[i]; ++j)
            {
               assert(decdecomp->subscipvars[i][j] != NULL);
               SCIP_CALL(SCIPhashmapInsert(varindexmap, decdecomp->subscipvars[i][j], (void*)varindex));
               varindex++;
            }
            for( j = 0; j < decdecomp->nsubscipconss[i]; ++j)
            {
   //            SCIPinfoMessage(scip, NULL, "%d, %d, %d; ", i, j, decdecomp->nsubscipconss[i] );
               assert(decdecomp->subscipconss[i][j] != NULL);
               SCIP_CALL(SCIPhashmapInsert(consindexmap, decdecomp->subscipconss[i][j], (void*)consindex));
               consindex++;
            }
         }

         SCIPdebugPrintf("Linking:\n");
         SCIPdebugPrintf("\tVars: %d", decdecomp->nlinkingvars);
         SCIPdebugPrintf("\tConss: %d\n\n", decdecomp->nlinkingconss);

         for( j = 0; j < decdecomp->nlinkingvars; ++j)
         {
            assert(decdecomp->linkingvars[j]);
            SCIP_CALL(SCIPhashmapInsert(varindexmap, decdecomp->linkingvars[j], (void*)varindex));
            varindex++;
         }
         for( j = 0; j < decdecomp->nlinkingconss; ++j)
         {
            assert(decdecomp->linkingconss[j]);
            SCIP_CALL(SCIPhashmapInsert(consindexmap, decdecomp->linkingconss[j], (void*)consindex));
            consindex++;
         }

         /* try this fix in order to assign indices to every variable (does not work!)*/
         /*
         for (j = 0; j < SCIPgetNVars(scip); ++j)
         {
            if(SCIPhashmapGetImage(varindexmap, SCIPgetVars(scip)[j]) == NULL)
            {
               SCIP_CALL(SCIPhashmapInsert(varindexmap, SCIPgetVars(scip)[j], (void*)varindex));
               varindex++;
            }
         }
         */
      }
      else if(decdecomp->type == DEC_STAIRCASE)
      {
         varindexmap = decdecomp->varindex;
         consindexmap = decdecomp->consindex;
      }
      for( i = 0; i < nconss; i++)
      {
         vars = SCIPgetVarsXXX(scip, conss[i]);
         nvars = SCIPgetNVarsXXX(scip, conss[i]);

         for( j = 0; j < nvars; j++)
         {
            /* if the problem has been created, output the whole model */
            if( SCIPgetStage(scip) == SCIP_STAGE_PROBLEM )
            {
               SCIPinfoMessage(scip, file, "%d, %d\n", SCIPvarGetIndex(vars[j]), i);
               continue;
            }

            /* somehow this assumes that the problem is presolved or we would have been*/
            if(!SCIPvarIsActive(vars[j]))
            {
              // continue;
            }
            /* if there is no decomposition, output the presolved model! */
            if(decdecomp == NULL)
            {
               SCIPinfoMessage(scip, file, "%d, %d\n", SCIPvarGetIndex(vars[j]), i);
            }
            /* if there is a decomposition, output the indices derived from the decomposition above*/
            else
            {

               assert(SCIPhashmapGetImage(varindexmap, vars[j]) != NULL);
               assert(SCIPhashmapGetImage(consindexmap, conss[i]) != NULL);

               SCIPinfoMessage(scip, file, "%d, %d\n",
                     SCIPhashmapGetImage(varindexmap, vars[j]),
                     SCIPhashmapGetImage(consindexmap, conss[i])
                     );

            }
         }
         SCIPfreeMemoryArrayNull(scip, &vars);
      }
      /* SCIPinfoMessage(scip, NULL, "varindex: %d, consindex: %d", varindex, consindex); */

      if(decdecomp != NULL && decdecomp->type != DEC_STAIRCASE)
      {
         SCIPhashmapFree(&varindexmap);
         SCIPhashmapFree(&consindexmap);
      }


   }

   return SCIP_OKAY;

}


/** write trailer of the file */
static
SCIP_RETCODE writeFileTrailer(
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file                                 /**< File pointer to write to */
   )
{
   SCIPinfoMessage(scip, file, "e\n");
   return SCIP_OKAY;
}

/*
 * Callback methods of reader
 */
/** copy method for reader plugins (called when SCIP copies plugins) */
#if 0
static
SCIP_DECL_READERCOPY(readerCopyGp)
{  /*lint --e{715}*/
   SCIPerrorMessage("method of gp reader not implemented yet\n");
   SCIPABORT(); /*lint --e{527}*/

   return SCIP_OKAY;
}
#else
#define readerCopyGp NULL
#endif
/** destructor of reader to free user data (called when SCIP is exiting) */
static
SCIP_DECL_READERFREE(readerFreeGp)
{
   SCIP_READERDATA* readerdata;

   assert(strcmp(SCIPreaderGetName(reader), READER_NAME) == 0);
   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);
   SCIPfreeMemory(scip, &readerdata);

   return SCIP_OKAY;
}


/** problem reading method of reader */
#define readerReadGp NULL



/** problem writing method of reader */
static
SCIP_DECL_READERWRITE(readerWriteGp)
{
   SCIP_READERDATA* readerdata;
   assert(scip != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   SCIP_CALL(SCIPwriteGp(scip, file, readerdata->decdecomp, TRUE));
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}



/*
 * reader specific interface methods
 */

SCIP_RETCODE SCIPReaderGpSetDecomp(
   SCIP* scip,
   DECDECOMP* decdecomp
   )
{
   SCIP_READER* reader;
   SCIP_READERDATA* readerdata;
   assert(scip != NULL);

   assert(decdecomp != NULL);

   reader = SCIPfindReader(scip, READER_NAME);
   assert(reader != NULL);

   readerdata = SCIPreaderGetData(reader);
   assert(readerdata != NULL);

   readerdata->decdecomp = decdecomp;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPwriteGp(
   SCIP* scip,                                /**< SCIP data structure */
   FILE* file,                                /**< File pointer to write to */
   DECDECOMP* decdecomp,                      /**< Decomposition pointer */
   SCIP_Bool writeDecomposition                  /**< whether to write decomposed problem */
   )
{
   char outname[SCIP_MAXSTRLEN];
   assert(scip != NULL);
   assert(file != NULL);
   if(writeDecomposition)
   {
      if(decdecomp == NULL)
      {
         SCIPwarningMessage("Cannot write decomposed problem if decomposition structure empty!");
         return SCIP_INVALIDDATA;
      }
   }
   /* print header */
   if(decdecomp == NULL)
   {
      SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   }
   else
   {
      SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d", SCIPgetProbName(scip), decdecomp->nblocks);
   }
   SCIP_CALL(writeFileHeader(scip, file, outname));

   /* write decomp information such as rectangles */
   if(writeDecomposition)
   {
      SCIP_CALL(writeDecompositionHeader(scip, file, decdecomp));
   }

   /* write the plot header*/
   SCIP_CALL(writePlotCommands(scip, file));


   /* write data */
   SCIP_CALL(writeData(scip, file, decdecomp));

   /* write file end */
   SCIP_CALL(writeFileTrailer(scip, file));
   return SCIP_OKAY;
}

/** includes the gp file reader in SCIP */
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_READERDATA* readerdata;

   /* create gp reader data */
   SCIP_CALL(SCIPallocMemory(scip, &readerdata));
   
   /* include gp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
         readerCopyGp, readerFreeGp, readerReadGp, readerWriteGp, readerdata) );

   /* add gp reader parameters */
   
   return SCIP_OKAY;
}
