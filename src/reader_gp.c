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

/**@file   reader_gp.c
 * @brief  GP file reader writing gnuplot files
 * @author Martin Bergner
 * @author Michael Feldmeier
 * @todo change output file type based on parameter
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "reader_gp.h"
#include "scip_misc.h"
#include "struct_decomp.h"
#include "cons_decomp.h"
#include "scip/cons_quadratic.h"

#define READER_NAME             "gpreader"
#define READER_DESC             "gnuplot file writer for matrix visualization"
#define READER_EXTENSION        "gp"

#define READERGP_GNUPLOT_BOXTEMPLATE(i, x1, y1, x2, y2) "set object %d rect from %.1f,%.1f to %.1f,%.1f fc rgb \"grey\"\n", (i), (x1), (y1), (x2), (y2)
#define READERGP_GNUPLOT_HEADER(outputname) "set terminal pdf\nset output \"%s.pdf\"\nunset xtics\nunset ytics\nunset border\nunset key\nset style fill solid 1.0 noborder\nset size ratio -1\n", (outputname)
#define READERGP_GNUPLOT_RANGES(xmax, ymax) "set xrange [0:%d]\nset yrange[%d:0]\n", (xmax), (ymax)
#define READERGP_GNUPLOT_PLOTCMD "plot \"-\" using 1:2:3 with circles fc rgb \"black\"\n"

/*
 * Local methods
 */

/** write file header with terminal etc. */
static
SCIP_RETCODE writeFileHeader(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   const char*           outname             /**< the name of the gnuplot outputname */
   )
{

   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER(outname));
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_RANGES(SCIPgetNVars(scip)+1, SCIPgetNConss(scip)+1));
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
   int startx;
   int starty;
   int endx;
   int endy;
   assert(scip != NULL);
   assert(file != NULL);
   assert(decdecomp != NULL);
   if( decdecomp->type == DEC_DECTYPE_UNKNOWN || decdecomp->nblocks == 0 )
   {
      return SCIP_OKAY;
   }

   if( decdecomp->type == DEC_DECTYPE_ARROWHEAD || decdecomp->type == DEC_DECTYPE_BORDERED )
   {
      startx = 0;
      starty = 0;
      endx = 0;
      endy = 0;

      for( i = 0; i < decdecomp->nblocks; ++i )
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

   if( decdecomp->type == DEC_DECTYPE_STAIRCASE )
   {
      startx = 0;
      starty = 0;
      endx = 0;
      endy = 0;

      for( i = 0; i < decdecomp->nblocks-1; ++i )
      {
         endx += decdecomp->nsubscipvars[i]+decdecomp->nstairlinkingvars[i];
         endy += decdecomp->nsubscipconss[i];
         SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+1, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
         startx = endx-decdecomp->nstairlinkingvars[i];
         starty = endy;
      }
      endx += decdecomp->nsubscipvars[i];
      endy += decdecomp->nsubscipconss[i];
      SCIPinfoMessage(scip, file, READERGP_GNUPLOT_BOXTEMPLATE(i+1, startx+0.5, starty+0.5, endx+0.5, endy+0.5));
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
      assert(decdecomp->type == DEC_DECTYPE_ARROWHEAD
         || decdecomp->type == DEC_DECTYPE_BORDERED
         || decdecomp->type == DEC_DECTYPE_DIAGONAL
         || decdecomp->type == DEC_DECTYPE_UNKNOWN
         || decdecomp->type == DEC_DECTYPE_STAIRCASE);

      /* if we don't have staicase, but something else, go through the blocks and create the indices */
      if( decdecomp->type == DEC_DECTYPE_ARROWHEAD || decdecomp->type == DEC_DECTYPE_BORDERED || decdecomp->type == DEC_DECTYPE_DIAGONAL )
      {
         size_t varindex = 1;
         size_t consindex = 1;

         SCIP_CALL( SCIPhashmapCreate(&varindexmap, SCIPblkmem(scip), SCIPgetNVars(scip)) );
         SCIP_CALL( SCIPhashmapCreate(&consindexmap, SCIPblkmem(scip), SCIPgetNConss(scip)) );

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
            for( j = 0; j < decdecomp->nsubscipconss[i]; ++j )
            {
               assert(decdecomp->subscipconss[i][j] != NULL);
               SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->subscipconss[i][j], (void*)consindex) );
               consindex++;
            }
         }

         SCIPdebugPrintf("Linking:\n");
         SCIPdebugPrintf("\tVars: %d", decdecomp->nlinkingvars);
         SCIPdebugPrintf("\tConss: %d\n\n", decdecomp->nlinkingconss);

         for( j = 0; j < decdecomp->nlinkingvars; ++j )
         {
            assert(decdecomp->linkingvars[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(varindexmap, decdecomp->linkingvars[j], (void*)varindex) );
            varindex++;
         }
         for( j = 0; j < decdecomp->nlinkingconss; ++j )
         {
            assert(decdecomp->linkingconss[j] != NULL);
            SCIP_CALL( SCIPhashmapInsert(consindexmap, decdecomp->linkingconss[j], (void*)consindex) );
            consindex++;
         }
      }
      else if( decdecomp->type == DEC_DECTYPE_STAIRCASE )
      {
         varindexmap = decdecomp->varindex;
         consindexmap = decdecomp->consindex;

         assert(varindexmap != NULL);
         assert(consindexmap != NULL);
      }
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
            SCIPinfoMessage(scip, file, "%d %d 0.5\n", SCIPvarGetIndex(curvars[j]), i);
            continue;
         }

         /* if there is no decomposition, output the presolved model! */
         if( decdecomp == NULL || decdecomp->type == DEC_DECTYPE_UNKNOWN )
         {
            SCIPinfoMessage(scip, file, "%d %d 0.5\n", SCIPvarGetIndex(curvars[j]), i);
         }
         /* if there is a decomposition, output the indices derived from the decomposition above*/
         else
         {
            assert(varindexmap != NULL);
            assert(consindexmap != NULL);
            assert(SCIPhashmapGetImage(varindexmap, SCIPvarGetProbvar(curvars[j])) != NULL);
            assert(SCIPhashmapGetImage(consindexmap, conss[i]) != NULL);

            SCIPinfoMessage(scip, file, "%d %d 0.5\n",
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

   SCIP_CALL( SCIPwriteGp(scip, file, DECgetBestDecomp(scip), TRUE) );

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/*
 * Methods to plot Qp Instances
 */


/** writes gnuplot code headerfor QP instance plot */
static
SCIP_RETCODE writeFileHeaderQp(
   SCIP*                 scip,                       /** SCIP data structure */
   FILE*                 file,                       /** File pointer to write to */
   const char*           outname,                    /** Name of the .pdf output file containing the plot  */
   int                   nvariables,                 /** Number of variables */
   int                   nlinearconstraints,         /** Number of linear constraints */
   int                   nquadraticconstraints       /** Number of quadratic constraints */
   )
{
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_HEADER(outname));
   SCIPinfoMessage(scip, file, READERGP_GNUPLOT_RANGES(nvariables+1, nlinearconstraints + nquadraticconstraints * nvariables + 1));
   return SCIP_OKAY;
}


/** writes gnuplot code for a single linear constraint */
static
SCIP_RETCODE writeLinearConstraint(
   SCIP*                 scip,                        /** SCIP data structure */
   FILE*                 file,                        /** File pointer to write to */
   SCIP_CONS*            linearconstraint,            /** Linear constraint */
   int*                  variableposition,            /** Array indicating at which position the variable at [i] is to be plot */
   int                   nvariables,                  /** Number of variables */
   SCIP_VAR**            tempvars,                    /** Allocated array to be used as temporary storage*/
   SCIP_Bool             writereorderedvariables,     /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                      /** y-coordinate in the plot */
   )
{
   int nconsvars;
   int i;

   nconsvars = GCGconsGetNVars( scip, linearconstraint );
   GCGconsGetVars( scip, linearconstraint, tempvars, nconsvars );

   for( i = 0; i < nconsvars; i++ )
   {
      int position;
      if( writereorderedvariables )
      {
         position = variableposition[SCIPvarGetIndex( tempvars[i])];
      }
      else
      {
         position = SCIPvarGetIndex( tempvars[i]);
      }
      SCIPinfoMessage(scip, file, "%d %d 0.5\n",position, yoffset);
   }
   return SCIP_OKAY;
}

/** writes gnuplot code for all linear constraints */
static
SCIP_RETCODE writeLinearConstraints(
   SCIP*                 scip,                       /** SCIP data structure */
   FILE*                 file,                       /** File pointer to write to */
   SCIP_CONS**           linearconstraints,          /** Array of linear constraints */
   int                   nlinearconstraints,         /** Number of linear constraints */
   int*                  variableposition,           /** Array indicating at which position the variable at [i] is to be plot */
   int                   nvariables,                 /** Number of variables */
   SCIP_Bool             writereorderedvariables,    /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                     /** y-coordinate in the plot */
   )
{
   SCIP_VAR** tempvars;
   int i;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tempvars, nvariables) );

   for(i=0; i < nlinearconstraints; i++)
   {
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(linearconstraints[i])), "linear") == 0);
      SCIP_CALL( writeLinearConstraint(scip, file, linearconstraints[i], variableposition, nvariables,  tempvars, writereorderedvariables, yoffset+i) );
   }

   SCIPfreeBlockMemoryArray(scip, &tempvars, nvariables);

   return SCIP_OKAY;
}

/** writes gnuplot code for a quadratic term */
static
SCIP_RETCODE writeQuadraticConstraintQuadraticTerm(
   SCIP*                 scip,                        /** SCIP data structure */
   FILE*                 file,                        /** File pointer to write to */
   SCIP_CONS*            quadraticconstraint,         /** Quadratic constraint */
   SCIP_QUADVARTERM*     quadraticterm,               /** Quadratic term to plot */
   int*                  variableposition,            /** Array indicating at which position the variable at [i] is to be plot */
   SCIP_Bool             writereorderedvariables,     /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                      /** y-coordinate in the plot */
   )
{
   int position;
   SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadraticconstraint, quadraticterm->var, &position) );
   if( writereorderedvariables )
   {
      position = variableposition[position];
   }
   SCIPinfoMessage(scip, file, "%d %d 0.5\n",position, yoffset + position);

   return SCIP_OKAY;
}

/** writes gnuplot code for a bilinear term */
static
SCIP_RETCODE writeQuadraticConstraintBilinearTerm(
   SCIP*                 scip,                         /** SCIP data structure */
   FILE*                 file,                         /** File pointer to write to */
   SCIP_CONS*            quadraticconstraint,          /** Quadratic constraint */
   SCIP_BILINTERM*       bilinearterm,                 /** Bilinear term to plot */
   int*                  variableposition,             /** Array indicating at which position the variable at [i] is to be plot */
   SCIP_Bool             writereorderedvariables,      /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                       /** y-coordinate in the plot */
   )
{
   int positionVar1;
   int positionVar2;
   SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadraticconstraint, bilinearterm->var1, &positionVar1) );
   SCIP_CALL( SCIPfindQuadVarTermQuadratic(scip, quadraticconstraint, bilinearterm->var2, &positionVar2) );
   if( writereorderedvariables )
   {
      positionVar1 = variableposition[positionVar1];
      positionVar2 = variableposition[positionVar2];
   }
   SCIPinfoMessage(scip, file, "%d %d 0.5\n",positionVar1, yoffset + positionVar2);
   SCIPinfoMessage(scip, file, "%d %d 0.5\n",positionVar2, yoffset + positionVar1);

   return SCIP_OKAY;
}

/** writes gnuplot code for a single quadratic constraint */
static
SCIP_RETCODE writeQuadraticConstraint(
   SCIP*                 scip,                       /** SCIP data structure */
   FILE*                 file,                       /** File pointer to write to */
   SCIP_CONS*            quadraticconstraint,        /** Quadratic constraint */
   int*                  variableposition,           /** Array indicating at which position the variable at [i] is to be plot */
   int                   nvariables,                 /** Number of variables */
   SCIP_Bool             writereorderedvariables,    /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                     /** y-coordinate in the plot*/
   )
{
   SCIP_BILINTERM* bilinearterms;
   SCIP_QUADVARTERM* quadraticterms;
   int nbilinearterms;
   int nquadraticterms;
   int i;

   SCIP_CALL( SCIPsortQuadVarTermsQuadratic(scip, quadraticconstraint) );

   bilinearterms = SCIPgetBilinTermsQuadratic(scip, quadraticconstraint);
   quadraticterms = SCIPgetQuadVarTermsQuadratic(scip, quadraticconstraint);
   nbilinearterms = SCIPgetNBilinTermsQuadratic(scip, quadraticconstraint);
   nquadraticterms = SCIPgetNQuadVarTermsQuadratic(scip, quadraticconstraint);

   /* write quadratic Terms */
   for( i = 0; i < nquadraticterms; i++ )
   {
      SCIP_CALL( writeQuadraticConstraintQuadraticTerm(scip, file, quadraticconstraint, &quadraticterms[i], variableposition, writereorderedvariables, yoffset) );
   }

   /* write bilinear Terms */
   for( i = 0; i < nbilinearterms; i++)
   {
      SCIP_CALL( writeQuadraticConstraintBilinearTerm(scip, file, quadraticconstraint, &bilinearterms[i], variableposition, writereorderedvariables, yoffset) );
   }

   return SCIP_OKAY;
}

/** writes gnuplot code for all quadratic constraints */
static
SCIP_RETCODE writeQuadraticConstraints(
   SCIP*                 scip,                       /** SCIP data structure */
   FILE*                 file,                       /** File pointer to write to */
   SCIP_CONS**           quadraticconstraints,       /** Array of quadratic constraints */
   int                   nquadraticconstraints,      /** Number of quadratic constraints */
   int*                  variableposition,           /** Array indicating at which position the variable at [i] is to be plot */
   int                   nvariables,                 /** Number of variables */
   SCIP_Bool             writereorderedvariables,    /** Use standard order of the variables or reposition as indicated in variablePosition? */
   int                   yoffset                     /** y-coordinate in the plot */
   )
{
   int i;
   int offset = yoffset;
   for( i = 0; i < nquadraticconstraints; i++ )
   {
      assert(strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(quadraticconstraints[i])), "quadratic") == 0);
      SCIP_CALL( writeQuadraticConstraint(scip, file, quadraticconstraints[i], variableposition, nvariables, writereorderedvariables, yoffset + nvariables * i) );
   }
   return SCIP_OKAY;
}


/*
 * reader specific interface methods
 */

/** writes gnuplot code to the specific file to plot QP instance*/
SCIP_RETCODE SCIPwriteQpGp(
   SCIP*                 scip,                     /** SCIP data structure */
   FILE*                 file,                     /** File pointer to write to */
   SCIP_CONS**           linearconstraints,        /** Array of linear constraints */
   int                   nlinearconstraints,       /** Number of linear constraints */
   SCIP_CONS**           quadraticconstraints,     /** Array of quadratic constraints */
   int                   nquadraticconstraints,    /** Number of quadratic constraints */
   int*                  variableposition,         /** Array indicating at which position the variable at [i] is to be plot */
   int                   nvariables,               /** Number of variables */
   SCIP_Bool             writereorderedvariables   /** Use standard order of the variables or reposition as indicated in variablePosition? */
   )
{
   char probname[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *name;

   assert(scip != NULL);
   assert(file != NULL);

   if( writereorderedvariables && variableposition == NULL )
   {
      SCIPwarningMessage(scip, "Cannot write reordered structure if position array is empty!");
      writereorderedvariables = FALSE;
   }

   /* sanitize filename */
   (void) SCIPsnprintf(probname, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(probname, NULL, &name, NULL, NULL);
   (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", name);

   /* write header */
   SCIP_CALL( writeFileHeaderQp(scip, file, outname, nvariables, nlinearconstraints, nquadraticconstraints) );

   /* write the plot header*/
   SCIP_CALL( writePlotCommands(scip, file) );

   /* write linear constraints */
   SCIP_CALL( writeLinearConstraints(scip, file, linearconstraints, nlinearconstraints, variableposition, nvariables, writereorderedvariables, 0) );

   /* write quadratic constraints */
   SCIP_CALL( writeQuadraticConstraints(scip, file, quadraticconstraints, nquadraticconstraints, variableposition, nvariables, writereorderedvariables, nlinearconstraints) );

   return SCIP_OKAY;
}


/** writes the decomposition to the specific file */
SCIP_RETCODE SCIPwriteGp(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file,               /**< File pointer to write to */
   DEC_DECOMP*           decdecomp,          /**< Decomposition pointer */
   SCIP_Bool             writeDecomposition  /**< whether to write decomposed problem */
   )
{
   char probname[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *name;

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

   /* print header */
   if( decdecomp == NULL )
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s", name);
   else
      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%c_%d", name, DECdetectorGetChar(decdecomp->detector), decdecomp->nblocks);

   SCIP_CALL( writeFileHeader(scip, file, outname) );

   /* write decomp information such as rectangles */
   if( writeDecomposition )
      SCIP_CALL( writeDecompositionHeader(scip, file, decdecomp) );

   /* write the plot header*/
   SCIP_CALL( writePlotCommands(scip, file) );

   /* write data */
   SCIP_CALL( writeData(scip, file, decdecomp) );

   /* write file end */
   SCIP_CALL( writeFileTrailer(scip, file) );
   return SCIP_OKAY;
}


/** includes the gp file reader into SCIP */
SCIP_RETCODE SCIPincludeReaderGp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   /* include gp reader */
   SCIP_CALL( SCIPincludeReader(scip, READER_NAME, READER_DESC, READER_EXTENSION,
      readerCopyGp, readerFreeGp, readerReadGp, readerWriteGp, NULL) );

   return SCIP_OKAY;
}

/** Compiles a given gnuplot file */
SCIP_RETCODE GCGcompileGpFile(
   SCIP*                  scip,          /**< SCIP data structure */
   char*                  filename       /** path to .gp file */
   )
{
   char command[SCIP_MAXSTRLEN];
   strcpy(command, "gnuplot ");
   strcat(command, filename);
   SCIPinfoMessage(scip, NULL, "%s\n", command);
   system(command);
   return SCIP_OKAY;
}

