/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   dec_arrowheur.c
 * @brief  arrowheur detector
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "dec_arrowheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"

#define DEC_DETECTORNAME      "arrowheur"    /**< name of the detector */
#define DEC_DESC              "enforces arrowhead structures using graph partitioning" /**< description of detector */
#define DEC_PRIORITY          1000           /**< priority of the detector */
#define DEC_DECCHAR           'a'            /**< display character of detector */
#define DEC_ENABLED           TRUE           /**< should detector be called by default */

/* Default parameter settings */
#define DEFAULT_VARWEIGHT         1          /**< weight for variable nodes */
#define DEFAULT_VARWEIGHTBIN      2          /**< weight for binary variable nodes */
#define DEFAULT_VARWEIGHTINT      2          /**< weight for integer variable nodes */
#define DEFAULT_VARWEIGHTCONT     1          /**< weight for continous variable nodes */
#define DEFAULT_VARWEIGHTIMPL     2          /**< weight for implicit integer variable nodes */
#define DEFAULT_CONSWEIGHT        5          /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED          1          /**< random seed for the hmetis call */
#define DEFAULT_TIDY              TRUE       /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES        0.2        /**< percentage of dummy vertices*/
#define DEFAULT_CONSWEIGHT_SETPPC 5          /**< weight for constraint hyperedges that are setpartitioning or covering constraints */
#define DEFAULT_MAXBLOCKS         10         /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS         2          /**< value for the minimum number of blocks to be considered */
#define DEFAULT_ALPHA             0.0        /**< factor for standard deviation of constraint weights */
#define DEFAULT_BETA              0.5        /**< factor of how the weight for equality and inequality constraints is distributed (keep 1/2 for the same on both) */
#define DEFAULT_METIS_UBFACTOR    5.0        /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE     FALSE      /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB  TRUE       /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_REALNAME          FALSE      /**< whether the metis name should be real or temporary */

/*
 * Data structures
 */

/** private detector data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
   SCIP_PTRARRAY* hedges;                    /**< variable array of hyperedges */
   SCIP_INTARRAY* copytooriginal;            /**< array mapping copied to original variables */
   int*           partition;                 /**< array storing vertex partitions */
   int            nvertices;                 /**< number of vertices */
   int*           varpart;                   /**< array storing variable partition */
   char           tempfile[SCIP_MAXSTRLEN];  /**< filename for the metis input file */

   /* weight parameters */
   int       varWeight;             /**< weight of a variable hyperedge */
   int       varWeightBinary;       /**< weight of a binary variable hyperedge */
   int       varWeightContinous;    /**< weight of a continuous variable hyperedge */
   int       varWeightInteger;      /**< weight of an integer variable hyperedge */
   int       varWeightImplint;      /**< weight of an implicit integer variable hyperedge */
   int       consWeight;            /**< weight of a constraint hyperedge */
   int       consWeightSetppc;      /**< weight of a setppc constraint hyperedge */
   SCIP_Real alpha;                 /**< factor for constraint coefficient value standard deviation */
   SCIP_Real beta;                  /**< factor for equality od inequality constraints */

   /* general parameters */
   SCIP_Real dummynodes;      /**< percent of dummy nodes */
   SCIP_Bool tidy;            /**< whether tempory metis files should be cleaned up */
   int       maxblocks;       /**< maximal number of blocks to test */
   int       minblocks;       /**< minimal number of blocks to test */

   /* metis parameters */
   int       randomseed;      /**< metis random seed */
   SCIP_Real metisubfactor;   /**< metis unbalance factor */
   SCIP_Bool metisverbose;    /**< should metis ouput be displayed */
   SCIP_Bool metisuseptyperb; /**< flag to indicate whether metis uses kway or rb partitioning */
   SCIP_Bool realname;        /**< flag to indicate real problem name or temporary filename for metis files */

   /* various data */
   SCIP_CLOCK* metisclock;    /**< clock to measure metis time */
   int         blocks;        /**< indicates the current block */
   SCIP_Bool   found;         /**< indicates whethere a decomposition has been found */
};

enum htype
{
   VARIABLE, CONSTRAINT
};
typedef enum htype hType;


/** hyper edge data structure */
struct hyperedge
{

   hType type;                /**< the type of the hyperegde (is it a split variable or a real constraint) */
   int *variableIds;          /**< the associated variable IDs that appear in the hyperedge */
   int nvariableIds;          /**< number of variable ids */
   int originalId;            /**< the original SCIP ID of this constraint or variable */

   int cost;                  /**< cost of the hyperedge */

};
typedef struct hyperedge HyperEdge;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** detector initialization method */
static
DEC_DECL_INITDETECTOR(initArrowheur)
{

   int i;
   int nvars;
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varpart, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      detectordata->varpart[i] = -1;
   }
   SCIP_CALL( SCIPcreatePtrarray(scip, &detectordata->hedges) );
   SCIP_CALL( SCIPcreateIntarray(scip, &detectordata->copytooriginal) );

   SCIP_CALL( SCIPcreateWallClock(scip, &detectordata->metisclock) );

   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitArrowheur)
{
   DEC_DETECTORDATA* detectordata;
   int i;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);
   /* copy data to decomp structure */
   if( !detectordata->found )
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->varpart);

   /* free dynamic arrays */
   for( i = 0; i <= SCIPgetPtrarrayMaxIdx(scip, detectordata->hedges); ++i )
   {
      HyperEdge* hedge;
      hedge = (HyperEdge*) SCIPgetPtrarrayVal(scip, detectordata->hedges, i);
      assert(hedge != NULL);

      SCIPfreeMemoryArray(scip, &hedge->variableIds);
      SCIPfreeMemory(scip, &hedge);
   }

   SCIP_CALL( SCIPfreeClock(scip, &detectordata->metisclock) );
   SCIP_CALL( SCIPfreePtrarray(scip, &detectordata->hedges) );
   SCIP_CALL( SCIPfreeIntarray(scip, &detectordata->copytooriginal) );
   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** compute weight of a hyperedge */
static
SCIP_RETCODE computeHyperedgeWeight(
   SCIP*             scip,          /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata,  /**< detector data data structure */
   SCIP_CONS*        cons,          /**< constraint belonging to the hyperegde */
   int*              cost           /**< pointer storing the hyperedge cost */
   )
{
   int j;
   int ncurvars;
   SCIP_Bool upgraded;
   SCIP_CONS* upgdcons;
   const char* hdlrname;

   upgraded = FALSE;
   *cost = detectordata->consWeight;
   hdlrname = SCIPconshdlrGetName(SCIPconsGetHdlr(cons));

   if( (strcmp("linear", hdlrname) == 0) )
   {
      SCIP_CALL( SCIPupgradeConsLinear(scip, cons, &upgdcons) );
      if( upgdcons != NULL )
         upgraded = TRUE;
   }
   else
   {
      upgdcons = cons;
   }

   if( upgdcons != NULL )
   {
     hdlrname =  SCIPconshdlrGetName(SCIPconsGetHdlr(upgdcons));
   }
   else
   {
     hdlrname =  SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
   }
   ncurvars = SCIPgetNVarsXXX(scip, cons);

   if( (strcmp("setppc", hdlrname) == 0) )
   {
      switch(SCIPgetTypeSetppc(scip, upgdcons))
      {
      case SCIP_SETPPCTYPE_COVERING:
         *cost = detectordata->consWeightSetppc;
         break;

      case SCIP_SETPPCTYPE_PARTITIONING:
         *cost = detectordata->consWeightSetppc;
         break;
      case SCIP_SETPPCTYPE_PACKING:
         *cost = detectordata->consWeightSetppc;
         break;
      default:
         *cost = detectordata->consWeight;
         break;
      }
   }
   else if( strcmp(hdlrname, "logicor") == 0 )
   {
      *cost = detectordata->consWeightSetppc;
   }
   else
   {
      SCIP_Real mean;
      SCIP_Real variance;
      SCIP_Real stddev;
      SCIP_Real * vals;

      mean = 0.0;
      variance = 0.0;
      vals = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &vals, ncurvars) );
         SCIP_CALL( SCIPgetValsXXX(scip, cons, vals, ncurvars) );
      }

      *cost = detectordata->consWeight;

      /* calculate variety using the normalized variance */
      for( j = 0; j < ncurvars; ++j )
      {
         assert(vals != NULL);
         mean += vals[j] / ncurvars;
      }
      if( ncurvars <= 1 )
      {
         variance = 0.0;
      }
      else
      {
         for( j = 0; j < ncurvars; ++j )
         {
            assert(vals != NULL);
            assert(ncurvars > 1);
            variance += pow((vals[j] - mean), 2.0) / (ncurvars-1);
         }
      }
      assert(variance >= 0);
      stddev = sqrt(variance);
      SCIPfreeBufferArrayNull(scip, &vals);

      /// @todo MAGIC NUMBER 2
      if( SCIPisEQ(scip, SCIPgetRhsXXX(scip, cons), SCIPgetLhsXXX(scip, cons)) )
      {
         /* we are dealing with an equality*/
         /*lint --e{524} */
         *cost = SCIPceil(scip, detectordata->beta*2.0*detectordata->consWeight+detectordata->alpha*stddev);
      }
      else
      {
         /*lint --e{524} */
         *cost = SCIPceil(scip, (1.0-detectordata->beta)*2.0*detectordata->consWeight+detectordata->alpha*stddev);
      }

   }
   if( upgraded == TRUE )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &upgdcons) );
   }
   return SCIP_OKAY;
}

/**
 * builds a graph structure out of the matrix.
 *
 * The function will create an HyperEdge for every constraint and every variable.
 * It will additionally create vertices for every variable and in particular
 * a copy of this variable for every constraint in which the variable has a
 * nonzero coefficient. The copies will be connected by the hyperedge for
 * the particular constraint and all copies of a variable will be connected by
 * the hyperedge belonging to that variable. The weight of these variable
 * hyperedges can be specified.
 *
 * @todo The nonzeroness is not checked, all variables in the variable array are considered
 */
static SCIP_RETCODE buildGraphStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< presolver data data structure */
   )
{
   SCIP_CONS **conss;
   SCIP_PTRARRAY *hedges;
   SCIP_INTARRAY *copytoorig;
   SCIP_INTARRAY** maporigtocopies;
   int *nmaporigtocopies;
   int nconss = SCIPgetNConss( scip );
   int nvars = SCIPgetNVars( scip );
   int id = 0;
   int nhyperedges = 0;
   int nvertices;
   int i;
   int j;
   int varWeight;
   int *copies;

   HyperEdge *hedge;
   assert(scip != NULL);
   assert(detectordata != NULL);

   conss = SCIPgetConss(scip);
   hedges = detectordata->hedges;
   copytoorig = detectordata->copytooriginal;

   SCIP_CALL( SCIPallocMemoryArray(scip, &copies, nconss) );

   /* we need at least nconss + nvars hyperedges */
   SCIP_CALL( SCIPextendPtrarray(scip, hedges,  0, nconss+nvars) );
   /* we have at least nvars may copy vertices */
   SCIP_CALL( SCIPextendIntarray(scip, copytoorig, 0, nvars) );

   /* map the original variable to all of its copies */
   SCIP_CALL( SCIPallocMemoryArray(scip, &maporigtocopies, nvars) );

   /* these are the number of copies for the given variable */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nmaporigtocopies, nvars) );

   /* initialize for every variable the list of copies */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL( SCIPcreateIntarray(scip, &maporigtocopies[i]) );
      nmaporigtocopies[i] = 0;
   }

   /* go through all constraints */
   for( i = 0; i < nconss; ++i )
   {
      int *varids;
      SCIP_VAR **vars;
      /* get the number of nonzeros in this constraint  */

      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );

      /* if there are no variables, skip the constraint */
      if( ncurvars == 0 )
         continue;

      /*
       * may work as is, as we are copying the constraint later regardless
       * if there are variables in it or not
       */
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, ncurvars) );
      SCIP_CALL( SCIPgetVarsXXX(scip, conss[i], vars, ncurvars) );

      /** @todo skip all variables that have a zero coeffient or where all coefficients add to zero */
      /** @todo Do more then one entry per variable actually work? */

      /* allocate a hyperedge for the constraint */
      SCIP_CALL( SCIPallocMemory(scip, &hedge) );
      hedge->type = CONSTRAINT;
      hedge->originalId = i;

      /* compute its weight*/
      SCIP_CALL( computeHyperedgeWeight(scip, detectordata, conss[i], &(hedge->cost)) );

      /* lets collect the variable ids of the variables */
      SCIP_CALL( SCIPallocBufferArray(scip, &varids, ncurvars) );
      hedge->nvariableIds = 0;

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var;
         int varIndex;

         /* if the variable is inactive, skip it */
         if( !SCIPisVarRelevant(vars[j]) )
            continue;

         var = SCIPvarGetProbvar(vars[j]);
         assert(var != NULL);
         varIndex = SCIPvarGetProbindex(var);
         /* assert that the variable is active and not multiaggregated, otherwise, the mapping will be wrong */
         /* the multiaggregation is useless, if we don't presolve, it might be interesting otherwise */
         assert(SCIPvarIsActive(var));
         assert(SCIPvarGetStatus(var) != SCIP_VARSTATUS_MULTAGGR);

         /* add the variable id to the end of the variable id array and update the size */
         varids[hedge->nvariableIds] = id;
         ++(hedge->nvariableIds);

         /* put the copied id (started from 0, should be the index of the nonzero entry) to the end of the map of original to copy ids  */
         SCIP_CALL( SCIPsetIntarrayVal(scip, maporigtocopies[varIndex], nmaporigtocopies[varIndex], id) );
         ++(nmaporigtocopies[varIndex]);
         SCIPdebugMessage("Adding %d at %d to copytoorig.\n", varIndex, id);
         SCIP_CALL( SCIPsetIntarrayVal(scip, copytoorig, id, varIndex) );
         ++id;
         /* Check the mapping here */
#ifdef SCIP_DEBUG
         {
            int k;
            SCIPdebugMessage("Cons %s (%d): ", SCIPconsGetName(conss[i]), i);
            SCIPdebugPrintf("Var %s (%d): ", SCIPvarGetName(var), varIndex);
            for( k = 0; k < nmaporigtocopies[varIndex]; ++k )
            {
               int copy;
               int orig;
               copy = SCIPgetIntarrayVal(scip, maporigtocopies[varIndex], k);
               orig = SCIPgetIntarrayVal(scip, copytoorig, copy);
               SCIPdebugPrintf("%d (%d), ", copy+1, orig);
               assert(varIndex == orig);
            }
            SCIPdebugPrintf("\n");
         }
#endif
      }

      if( hedge->nvariableIds > 1 )
      {
         /* if the hyperedge contains more then 0 variables, add it to the end */
         SCIP_CALL( SCIPduplicateMemoryArray(scip, &hedge->variableIds, varids, hedge->nvariableIds) );
         SCIP_CALL( SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge) );
         ++nhyperedges;
      }
      else
      {
         SCIPfreeMemory(scip, &hedge);
      }
      SCIPfreeBufferArray(scip, &varids);
      SCIPfreeBufferArray(scip, &vars);
   }

   /* build variable hyperedges */
   for( i = 0; i < nvars; i++ )
   {
      /*
       * handle variables only appearing in the objective
       * this is in a sense useless, as it increases the work for metis, but
       * it is consistent this way
       */
      int size;
      size = nmaporigtocopies[i];
      assert( size >= 0);

      /* if the hyperedge contains only 1 variable or less, skip it */
      if( size <= 1 )
         continue;

      switch ( SCIPvarGetType( SCIPgetVars( scip )[i] ) ) {
      case SCIP_VARTYPE_CONTINUOUS:
         varWeight = detectordata->varWeightContinous;
         break;
      case SCIP_VARTYPE_INTEGER:
         varWeight = detectordata->varWeightInteger;
         break;
      case SCIP_VARTYPE_IMPLINT:
         varWeight = detectordata->varWeightImplint;
         break;
      case SCIP_VARTYPE_BINARY:
         varWeight = detectordata->varWeightBinary;
         break;
      default:
         varWeight = detectordata->varWeight;
         break;
      }

      SCIP_CALL( SCIPallocMemory(scip, &hedge) );
      hedge->type = VARIABLE;
      hedge->originalId = i;
      hedge->nvariableIds = 0;
      hedge->variableIds = NULL;
      hedge->cost = varWeight;
      SCIP_CALL( SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge) );
      ++nhyperedges;

      /* link the copies together */
      SCIP_CALL( SCIPallocMemoryArray(scip, &hedge->variableIds, size) );
      hedge->nvariableIds = size;

      SCIPdebugMessage("nvars hedge: ");

      for( j = 0; j < size; ++j )
      {
         hedge->variableIds[j] = SCIPgetIntarrayVal(scip, maporigtocopies[i], j);
         SCIPdebugPrintf("%d, ", hedge->variableIds[j]+1);
      }
      SCIPdebugPrintf("\n");
      SCIP_CALL( SCIPfreeIntarray(scip, &maporigtocopies[i]) );
   }
   SCIPfreeMemoryArray(scip, &maporigtocopies);
   SCIPfreeMemoryArray(scip, &nmaporigtocopies);
   SCIPfreeMemoryArray(scip, &copies);
   nvertices = id;
   detectordata->nvertices = nvertices;

   assert(SCIPgetPtrarrayMaxIdx(scip, hedges)+1 == nhyperedges);
   //   assert(nvertices == SCIPgetIntarrayMaxIdx(scip, copytoorig)+1);


   return SCIP_OKAY;
}

/** will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   SCIP_RESULT*          result              /**< result indicating whether the detection was successful */
   )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];

   int status;
   int i;
   int nvertices;
   int* partition;


   SCIP_FILE *zfile;
   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   remainingtime = DECgetRemainingTime(scip);
   nvertices = detectordata->nvertices;

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      detectordata->varpart[i] = -1;
   }

   /* call metis via syscall as there is no library usable ... */
   if( !SCIPisInfinity(scip, DECgetRemainingTime(scip)) )
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f;hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               remainingtime,
               detectordata->tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      (void) SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               detectordata->tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL( SCIPresetClock(scip, detectordata->metisclock) );
   SCIP_CALL( SCIPstartClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " %d", detectordata->blocks );
   status = system( metiscall );

   SCIP_CALL( SCIPstopClock(scip, detectordata->metisclock) );
   SCIPdebugMessage("time left before metis started: %f, time metis spend %f, remainingtime: %f\n", remainingtime, SCIPgetClockTime(scip, detectordata->metisclock),  remainingtime-SCIPgetClockTime(scip, detectordata->metisclock) );

   /* check error codes */
   if( status == -1 )
   {
      SCIPerrorMessage("System call did not succed: %s\n", strerror( errno ));
      SCIPerrorMessage("Call was %s\n", metiscall);
   }
   else if( status != 0 )
   {

      SCIPerrorMessage("Calling hmetis unsuccessful! See the above error message for more details.\n");
      SCIPerrorMessage("Call was %s\n", metiscall);
   }

   /* exit gracefully in case of errors */
   if( status != 0 )
   {
      return SCIP_ERROR;
   }

   /*
    * parse the output into the vector
    * alloc the memory
    */
   if( detectordata->partition == NULL )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->partition, nvertices) );
   }

   assert(detectordata->partition != NULL);
   partition = detectordata->partition;

   (void) SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", detectordata->tempfile, detectordata->blocks);

   zfile = SCIPfopen(metisout, "r");
   i = 0;
   while( !SCIPfeof(zfile) && i < nvertices )
   {
      int temp;
      if( SCIPfgets(line, SCIP_MAXSTRLEN, zfile) == NULL )
      {
         SCIPerrorMessage("Line could not be read\n");
         return SCIP_READERROR;
      }

      temp = atoi(line);
      assert(temp >= 0 && temp <= detectordata->blocks-1);
      partition[i] = temp;
      i++;
   }

   if( i != nvertices )
   {
      SCIPerrorMessage("Couldn't read partition for all vertices.\n");
      return SCIP_READERROR;
   }

   SCIPfclose(zfile);

   /* if desired delete the temoprary metis file */
   if( detectordata->tidy )
   {
      status = unlink( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s\n", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", detectordata->tempfile);
   }
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** maps the partitions for the disaggregated vertices to the original vertices */
static
SCIP_RETCODE assignBlocksToOriginalVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata        /**< presolver data data structure */
   )
{

   int i;
#ifndef NDEBUG
   int nvars;
#endif
   int *partition;
   int *origpart;
   int nvertices;
   assert(scip != NULL);
   assert(detectordata != NULL);

   nvertices = detectordata->nvertices;
   partition = detectordata->partition;
   origpart = detectordata->varpart;

#ifndef NDEBUG
   nvars = SCIPgetNVars( scip );
#endif

    /* go through the new vertices */
   for( i = 0; i < nvertices ; ++i )
   {
      int originalId;
      /* find out the original id (== index of the var in the vars array) */
      originalId = SCIPgetIntarrayVal(scip, detectordata->copytooriginal, i);

      /* add the id to the set of ids for the original vertex */
      assert(originalId >= 0 && originalId < nvars);
      assert(partition[i] >= 0);
      if( origpart[originalId] == -1 )
      {
         origpart[originalId] = partition[i];
      }
      else if( origpart[originalId] >= 0 )
      {
         if( origpart[originalId] != partition[i] )
         {
            origpart[originalId] = -2;
         }
      }
      assert(origpart[originalId] == -2 || origpart[originalId] >= 0);
      assert(origpart[originalId] <= detectordata->blocks);
   }

   return SCIP_OKAY;
}

/** builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DETECTORDATA*     detectordata,       /**< presolver data data structure */
   DEC_DECOMP*           decdecomp,          /**< decdecomp data structure */
   int                   nblocks,            /**< number of blocks for this decomposition */
   SCIP_RESULT*          result              /**< result pointer */
   )
{
   SCIP_Bool *isVarHandled;
   SCIP_VAR*** subscipvars;
   SCIP_VAR** linkingvars;
   SCIP_CONS*** subscipconss;
   SCIP_CONS** linkingconss;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;

   int *nsubscipvars;
   int *nsubscipconss;
   int nlinkingvars;
   int nlinkingconss;
   int i;
   int j;
   SCIP_CONS **conss;
   int nconss;
   SCIP_VAR **vars;
   int nvars;
   SCIP_Bool emptyblocks = FALSE;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nconss = SCIPgetNConss( scip );
   nvars = SCIPgetNVars( scip );

   conss = SCIPgetConss( scip );
   vars = SCIPgetVars( scip );

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );

   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) );

      nsubscipconss[i] = 0;
      nsubscipvars[i] = 0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   nlinkingconss = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nvars) );
   nlinkingvars = 0;

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nconss) );

   SCIP_CALL( SCIPallocBufferArray(scip, &isVarHandled, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      isVarHandled[i] = FALSE;
   }

   /* go through all of the constraints */
   for( i = 0; i < nconss; i++ )
   {
      int consblock = -1;
      int ncurvars;
      SCIP_VAR **curvars;

      if( strcmp(SCIPconshdlrGetName(SCIPconsGetHdlr(conss[i])), "origbranch") == 0 )
         continue;

      /* sort the variables into corresponding buckets */
      ncurvars = SCIPgetNVarsXXX( scip, conss[i] );
      curvars = NULL;
      if( ncurvars > 0 )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &curvars, ncurvars) );
         SCIP_CALL( SCIPgetVarsXXX( scip, conss[i], curvars, ncurvars) );
      }

      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         int varblock = -1;
         if( !SCIPisVarRelevant(curvars[j]) )
            continue;

         var = SCIPvarGetProbvar(curvars[j]);

         assert(var != NULL);
         assert(SCIPvarIsActive(var));
         assert(!SCIPvarIsDeleted(var));
         /*
          * if the variable has already been handled, we do not need to look
          * at it again and only need to set the constraint
          */
         if( !isVarHandled[SCIPvarGetProbindex(var)] )
         {
            isVarHandled[SCIPvarGetProbindex( var )] = TRUE;
            /*
             * if the following assertion fails, the picture and the mapping is
             * certainly wrong!
             */
            assert(vars[SCIPvarGetProbindex(var)] == var);
            assert(SCIPvarGetProbindex(var) < nvars);
            assert(detectordata->varpart[SCIPvarGetProbindex(var)] < detectordata->blocks);
            /* get the number of blocks the current variable is in */
            assert(detectordata->varpart[SCIPvarGetProbindex(var)] == -2 ||
                   detectordata->varpart[SCIPvarGetProbindex(var)] >= 0);
            /* if the variable is in exactly one block */
            if( detectordata->varpart[SCIPvarGetProbindex(var)] != -2 )
            {
               /* then the partition is given */
               varblock = detectordata->varpart[SCIPvarGetProbindex(var)];
               assert(varblock < detectordata->blocks);
               subscipvars[varblock][nsubscipvars[varblock]] = var;
               //               SCIPdebugMessage("v: %s\n", SCIPvarGetName(var));
               ++(nsubscipvars[varblock]);
            }
            /*
             * if the variable is a linking variable, don't update the constraint
             * block and add the variable to the linking variables
             */
            else
            {
               varblock = detectordata->blocks+1;
               //               SCIPdebugMessage("v: %s\n", SCIPvarGetName(var));
               linkingvars[nlinkingvars] = var;
               ++nlinkingvars;
            }

            /* finally set the hashmap image */
            assert(!SCIPhashmapExists(vartoblock, var));
            SCIP_CALL( SCIPhashmapInsert(vartoblock, var, (void*) (size_t) varblock) );
         }
         else
         {
            varblock = (int)(size_t)SCIPhashmapGetImage(vartoblock, var);
            assert(varblock == detectordata->varpart[SCIPvarGetProbindex(var)] ||  detectordata->varpart[SCIPvarGetProbindex(var)] == -2);
         }


         /*
          * if the variable is not a linking variable, add it to the correct
          * block and update the block of the constraint, if necessary
          */
         if( varblock <= detectordata->blocks )
         {
            /*
             * if the block of the constraint has not been set yet, set it to
             * the block of the current variable
             */
            if( consblock == -1 )
            {
               consblock = varblock;
            }
            /*
             * if the block has been set but the current variable belongs to a
             * different block, we have a linking constraint
             */
            else if( consblock >= 0 && consblock != varblock )
            {
               consblock = -2;
            }

            /*
             * At this point the constraint is either in the border or it
             * belongs to the same block as the current variable
             */
            assert(consblock == -2 ||consblock == varblock);

         }
      }
      SCIPfreeBufferArrayNull(scip, &curvars);

      /*
       *  sort the constraints into the corresponding bucket
       *
       *  if the constraint is linking, put it there
       */
      if( consblock < 0 )
      {
         size_t block;

         block = detectordata->blocks +1;
         linkingconss[nlinkingconss] = conss[i];
         ++nlinkingconss;
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*)(block)) );

      }
      /* otherwise put it in its block */
      else
      {
         subscipconss[consblock][nsubscipconss[consblock]] = conss[i];
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL( SCIPhashmapInsert(constoblock, conss[i], (void*) (size_t) consblock) );
         ++(nsubscipconss[consblock]);
      }
   }

   /*
    * go through all variables and look at the not handled variables and add
    * them to the correct partition
    */

   for( i = 0; i < nvars; i++ )
   {
      int partitionOfVar;
      if( isVarHandled[i] )
         continue;

      partitionOfVar = -1;
      if( detectordata->varpart[i] >= 0 )
      {
         partitionOfVar = detectordata->varpart[i];
      }

      if( partitionOfVar != -1 )
      {
         subscipvars[partitionOfVar][nsubscipvars[partitionOfVar]] = SCIPvarGetProbvar(vars[i]);
         //         SCIPdebugMessage("v: %s\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])));
         ++nsubscipvars[partitionOfVar];
      }
      else
      {
         //         SCIPdebugMessage("v: %s\n", SCIPvarGetName(SCIPvarGetProbvar(vars[i])));
         linkingvars[nlinkingvars] = SCIPvarGetProbvar(SCIPvarGetProbvar(vars[i]));
         ++nlinkingvars;
      }

   }

   SCIPfreeBufferArray(scip, &isVarHandled);

   /* first, make sure that there are constraints in every block, otherwise the hole thing is useless */
   for( i = 0; i < detectordata->blocks; ++i )
   {
      if( nsubscipconss[i] == 0 )
      {
         SCIPdebugMessage("Block %d does not have any constraints!\n", i);
         emptyblocks = TRUE;
      }
   }

   if( !emptyblocks )
   {
      /* copy the local data to the decomp structure */
      DECdecompSetNBlocks(decdecomp, nblocks);
      DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
      SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
      SCIP_CALL( DECdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
      if( nlinkingconss > 0 )
      {
         SCIP_CALL( DECdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
         DECdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
      }
      if( nlinkingvars > 0 )
      {
         DECdecompSetType(decdecomp, DEC_DECTYPE_ARROWHEAD);
         SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars) );
      }
      DECdecompSetVartoblock(decdecomp, vartoblock);
      DECdecompSetConstoblock(decdecomp, constoblock);
   }
   else {
      SCIPhashmapFree(&constoblock);
      SCIPhashmapFree(&vartoblock);
   }
   /* free all local data */
   for( i = 0; i < nblocks; ++i )
   {
      SCIPfreeBufferArray(scip, &subscipconss[i]);
      SCIPfreeBufferArray(scip, &subscipvars[i]);
   }

   SCIPfreeBufferArray(scip, &subscipconss);
   SCIPfreeBufferArray(scip, &subscipvars);
   SCIPfreeBufferArray(scip, &nsubscipconss);
   SCIPfreeBufferArray(scip, &nsubscipvars);

   SCIPfreeBufferArray(scip, &linkingconss);
   SCIPfreeBufferArray(scip, &linkingvars);

   *result = emptyblocks? SCIP_DIDNOTFIND:SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** creates the temporary metis input file */
static
SCIP_RETCODE createMetisFile(
   SCIP*                 scip,               /**< SCIP data struture */
   DEC_DETECTORDATA*     detectordata        /**< detector data structure */
   )
{
   FILE* file;
   int temp_filedes;
   int i;
   int j;
   SCIP_PTRARRAY *hedges;
   int status;
   int nvertices;
   int ndummyvertices;

   hedges = detectordata->hedges;
   nvertices = detectordata->nvertices;
   /*lint --e{524} */
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);

   if( !detectordata->realname )
   {
      (void) SCIPsnprintf(detectordata->tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   }
   else
   {
      (void) SCIPsnprintf(detectordata->tempfile, SCIP_MAXSTRLEN, "gcg-%s-XXXXXX", SCIPgetProbName(scip));
   }

   if( (temp_filedes = mkstemp(detectordata->tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror( errno ));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", detectordata->tempfile);

   file = fdopen(temp_filedes, "w");
   if( file == NULL )
   {
      SCIPerrorMessage("Could not open temporary metis file!\n");
      return SCIP_FILECREATEERROR;
   }

   SCIPinfoMessage(scip, file, "%d %d 1\n", SCIPgetPtrarrayMaxIdx(scip, hedges)+1, nvertices+ndummyvertices);
   for( i = 0; i <= SCIPgetPtrarrayMaxIdx(scip, hedges); i++ )
   {
      HyperEdge* hedge;
      hedge = (HyperEdge*) SCIPgetPtrarrayVal(scip, hedges, i);
      assert(hedge != NULL);
      assert(hedge->nvariableIds != 0);
      SCIPinfoMessage(scip, file, "%d ", hedge->cost);

      for( j = 0; j < hedge->nvariableIds; ++j )
      {
         SCIPinfoMessage(scip, file, "%d ",hedge->variableIds[j] + 1 );
      }
      SCIPinfoMessage(scip, file, "\n");
   }
   status = fclose(file);

   if( status == -1 )
   {
      SCIPerrorMessage("Could not close '%s'\n", detectordata->tempfile);
      return SCIP_WRITEERROR;
   }


   return SCIP_OKAY;
}

/** detection callback method */
static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildArrowhead)
{
   int i;
   int j;
   int ndecs;
   int status;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   ndecs = detectordata->maxblocks-detectordata->minblocks+1;
   *ndecdecomps = 0;

   /* allocate space for output data */
   assert(detectordata->maxblocks >= detectordata->minblocks);
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, ndecs) );

   /* build the hypergraph structure from the original problem */
   SCIP_CALL( buildGraphStructure(scip, detectordata) );

   for( i = 0; i < ndecs; ++i )
   {
      SCIP_CALL( DECdecompCreate(scip, &(*decdecomps)[i]) );
   }

   SCIP_CALL( createMetisFile(scip, detectordata) );

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting Arrowhead structure:");
   for( j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      detectordata->blocks = i;
      /* get the partitions for the new variables from metis */
      SCIP_CALL( callMetis(scip, detectordata, result) );

      if( *result != SCIP_SUCCESS )
      {
         *result = SCIP_DIDNOTFIND;
         return SCIP_OKAY;
      }
      else
      {
         detectordata->found = TRUE;
      }

      /* deduce the partitions for the original variables */
      SCIP_CALL( assignBlocksToOriginalVariables( scip, detectordata) );

      SCIP_CALL( buildTransformedProblem(scip, detectordata, (*decdecomps)[j], i, result) );
      if( *result == SCIP_SUCCESS )
      {
         *ndecdecomps += 1;
         ++j;
      }
   }
   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " done, %d decompositions found.\n",  *ndecdecomps);
   for( i = *ndecdecomps; i < ndecs; ++i )
   {
      DECdecompFree(scip, &((*decdecomps)[i]) );
   }

   SCIP_CALL( SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps) );

   if( detectordata->tidy )
   {
      status = unlink( detectordata->tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** creates the arrowheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionArrowheur(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->blocks = -1;

   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, detectordata, detectAndBuildArrowhead, initArrowheur, exitArrowheur) );


   /* add arrowheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/arrowheur/beta", "factor on how heavy equality (beta) and inequality constraints are measured", &detectordata->beta, FALSE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL ) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/arrowheur/alpha", "factor on how heavy the standard deviation of the coefficients is measured", &detectordata->alpha, FALSE, DEFAULT_ALPHA, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/varWeight", "Weight of a variable hyperedge", &detectordata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/varWeightBinary", "Weight of a binary variable hyperedge", &detectordata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/varWeightContinous", "Weight of a continuos variable hyperedge", &detectordata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/varWeightImplint", "Weight of a implicit integer variable hyperedge", &detectordata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/varWeightInteger", "Weight of a integer variable hyperedge", &detectordata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/arrowheur/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/arrowheur/dummynodes", "percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/arrowheur/consWeightSetppc", "Weight for constraint hyperedges that are setpartitioning or covering constraints", &detectordata->consWeightSetppc, FALSE, DEFAULT_CONSWEIGHT_SETPPC, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddRealParam(scip, "detectors/arrowheur/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/arrowheur/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/arrowheur/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/arrowheur/realname", "Should the problem be used for metis files or a temporary name", &detectordata->realname, FALSE, DEFAULT_REALNAME, NULL, NULL) );

   return SCIP_OKAY;
}
