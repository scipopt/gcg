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
 * @ingroup DETECTORS
 * @brief  arrowheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

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

#define DEC_DETECTORNAME      "arrowheur"   /**< name of the detector */
#define DEC_PRIORITY          1000          /**< priority of the detector */

/* Default parameter settings */
#define DEFAULT_BLOCKS                    2     /**< number of blocks */
#define DEFAULT_VARWEIGHT                 1     /**< weight for variable nodes */
#define DEFAULT_VARWEIGHTBIN              2     /**< weight for binary variable nodes */
#define DEFAULT_VARWEIGHTINT              2     /**< weight for integer variable nodes */
#define DEFAULT_VARWEIGHTCONT             1     /**< weight for continous variable nodes */
#define DEFAULT_VARWEIGHTIMPL             2     /**< weight for implicit integer variable nodes */
#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES                0.2   /**< percentage of dummy vertices*/
#define DEFAULT_CONSWEIGHT_SETPPC         5     /**< weight for constraint hyperedges that are setpartitioning or covering constraints */
#define DEFAULT_MAXBLOCKS                 10    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */
#define DEFAULT_ALPHA                     0.0   /**< factor for standard deviation of constraint weights */
#define DEFAULT_BETA                      0.5   /**< factor of how the weight for equality and inequality constraints is distributed (keep 1/2 for the same on both) */
#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY


#define DWSOLVER_REFNAME(name, blocks, varcont, varint, cons, dummy, alpha, beta, conssetppc)  \
   "%s_%d_%d_%d_%d_%.1f_%.1f_%.1f_%d_ref.txt", \
   (name), (blocks), (varcont), (varint), (cons), (dummy), (alpha), (beta), (conssetppc)

#define GP_NAME(name, blocks, varcont, varint, cons, dummy, alpha, beta, conssetppc)  \
   "%s_%d_%d_%d_%d_%.1f_%.1f_%.1f_%d.gp", \
   (name), (blocks), (varcont), (varint), (cons), (dummy), (alpha), (beta), (conssetppc)

/*
 * Data structures
 */

/** score data structure **/
struct SCIP_ArrowheurScores
{
   SCIP_Real borderscore;
   SCIP_Real minkequicutscore;
   SCIP_Real equicutscorenormalized;
   SCIP_Real densityscore;
   SCIP_Real linkingscore;
};
typedef struct SCIP_ArrowheurScores SCIP_ARROWHEURSCORES;
/** presolver data */
struct DEC_DetectorData
{
   /* Graph stuff for hmetis */
   SCIP_PTRARRAY *hedges;
   SCIP_INTARRAY *copytooriginal;
   int *partition;
   int nvertices;
   int nhyperedges;
   int *varpart;

   /* Stuff to get the dw-solver to work*/
   SCIP_HASHMAP *constolpid;

   SCIP_Bool tidy;
   int blocks;
   int maxblocks;
   int minblocks;
   int varWeight;
   int varWeightBinary;
   int varWeightContinous;
   int varWeightInteger;
   int varWeightImplint;
   int consWeight;
   int randomseed;
   SCIP_Bool found;
   SCIP_Real dummynodes;
   int consWeightSetppc;
   SCIP_Real alpha;
   SCIP_Real beta;

   SCIP_Real metisubfactor;
   SCIP_Bool metisverbose;
   SCIP_Bool metisuseptyperb;
   SCIP_CLOCK *metisclock;
   int priority;
};

enum htype
{
   VARIABLE, CONSTRAINT
};
typedef enum htype hType;

struct hyperedge
{

   hType type;       ///< The type of the hyperegde (is it a split variable or a real constraint)
   int *variableIds; ///< the associated variable IDs that appear in the hyperedge
  // int copyId;       ///< The ids of the associated copy
   int nvariableIds; ///< number of variable ids
   int originalId;   ///< the original SCIP ID of this constraint or variable

   int cost;

};
typedef struct hyperedge HyperEdge;

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

/** Prints the score of the decomposition */
static
SCIP_RETCODE printArrowheurScores(
      SCIP*                 scip,
      DEC_DETECTORDATA*   detectordata,
      SCIP_ARROWHEURSCORES* scores
      )
{
   char name[SCIP_MAXSTRLEN];
   SCIPsnprintf(name, SCIP_MAXSTRLEN, DWSOLVER_REFNAME(SCIPgetProbName(scip),
         detectordata->blocks,
         detectordata->varWeightContinous,
         detectordata->varWeightInteger,
         detectordata->consWeight,
         detectordata->dummynodes,
         detectordata->alpha,
         detectordata->beta,
         detectordata->consWeightSetppc));

   SCIPdebugMessage("SCORES:\t%s\t%s\t%f\t%f\t%f\t%f\t%f\n",
         SCIPgetProbName(scip), name,
         scores->borderscore,
         scores->densityscore,
         scores->linkingscore,
         scores->minkequicutscore,
         scores->equicutscorenormalized);

   return SCIP_OKAY;
}

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
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varpart, nvars));
   for( i = 0; i < nvars; ++i )
   {
      detectordata->varpart[i] = -1;
   }
   SCIP_CALL(SCIPcreatePtrarray(scip, &detectordata->hedges));
   SCIP_CALL(SCIPcreateIntarray(scip, &detectordata->copytooriginal));

   SCIP_CALL(SCIPhashmapCreate(&detectordata->constolpid, SCIPblkmem(scip), nconss));

   /* initialise consttolpid hashmap and hyper edges array */
   for( i = 0; i < nconss; ++i )
   {
      SCIP_CALL(SCIPhashmapInsert(detectordata->constolpid, SCIPgetConss(scip)[i], (void*)(size_t)i));
   }

   SCIP_CALL(SCIPcreateWallClock(scip, &detectordata->metisclock));

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
   if( !detectordata->found)
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   SCIPfreeMemoryArray(scip, &detectordata->partition);
   SCIPfreeMemoryArray(scip, &detectordata->varpart);

   /* free hash map */
   SCIPhashmapFree(&detectordata->constolpid);

   /* free dynamic arrays */
   for( i = 0; i <= SCIPgetPtrarrayMaxIdx(scip, detectordata->hedges); ++i )
   {
      HyperEdge* hedge;
      hedge = (HyperEdge*) SCIPgetPtrarrayVal(scip, detectordata->hedges, i);
      assert(hedge != NULL);

      SCIPfreeMemoryArray(scip, &hedge->variableIds);
      SCIPfreeMemory(scip, &hedge);
   }

   SCIP_CALL(SCIPfreeClock(scip, &detectordata->metisclock));
   SCIPfreePtrarray(scip, &detectordata->hedges);
   SCIPfreeIntarray(scip, &detectordata->copytooriginal);
   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

static
int computeHyperedgeWeight(
   SCIP*             scip,
   DEC_DETECTORDATA* detectordata,
   SCIP_CONS*        cons,
   int*              cost
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

   if((strcmp("linear", hdlrname) == 0))
   {
      SCIP_CALL(SCIPupgradeConsLinear(scip, cons, &upgdcons));
      if(upgdcons != NULL)
         upgraded = TRUE;
   }
   else
   {
      upgdcons = cons;
   }


   if(upgdcons != NULL)
   {
     hdlrname =  SCIPconshdlrGetName(SCIPconsGetHdlr(upgdcons));
   }
   else
   {
     hdlrname =  SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
   }
   ncurvars = SCIPgetNVarsXXX(scip, cons);

   if((strcmp("setppc", hdlrname) == 0))
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
   else if(strcmp(hdlrname, "logicor") == 0)
   {
      *cost = detectordata->consWeightSetppc;
   }
   else
   {
      double mean;
      double variance;
      double stddev;
      SCIP_Real * vals;

      mean = 0.0;
      variance = 0.0;
      vals = SCIPgetValsXXX(scip, cons);

      *cost = detectordata->consWeight;

      /* calculate variety using the normalized variance */
      for( j = 0; j < ncurvars; ++j )
      {
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
            assert(ncurvars > 1);
            variance += pow((vals[j] - mean), 2.0) / (ncurvars-1);
         }
      }
      assert(variance >= 0);
      stddev = sqrt(variance);
      SCIPfreeMemoryArray(scip, &vals);

      // TODO: MAGIC NUMBER 2
      if( SCIPisEQ(scip, SCIPgetRhsXXX(scip, cons), SCIPgetLhsXXX(scip, cons)) )
      {
         /* we are dealing with an equality*/
         *cost = SCIPceil(scip, detectordata->beta*2.0*detectordata->consWeight+detectordata->alpha*stddev);
      }
      else
      {
         *cost = SCIPceil(scip, (1.0-detectordata->beta)*2.0*detectordata->consWeight+detectordata->alpha*stddev);
      }

   }
   if( upgraded == TRUE )
   {
      SCIP_CALL(SCIPreleaseCons(scip, &upgdcons));
   }
   return SCIP_OKAY;
}

/**
 * Builds a graph structure out of the matrix.
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
      SCIP*             scip,          /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata   /**< presolver data data structure */
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

   SCIP_CALL(SCIPallocMemoryArray(scip, &copies, nconss));

   /* we need at least nconss + nvars hyperedges */
   SCIP_CALL(SCIPextendPtrarray(scip, hedges,  0, nconss+nvars));
   /* we have at least nvars may copy vertices */
   SCIP_CALL(SCIPextendIntarray(scip, copytoorig, 0, nvars));

   /* map the original variable to all of its copies */
   SCIP_CALL(SCIPallocMemoryArray(scip, &maporigtocopies, nvars));

   /* these are the number of copies for the given variable */
   SCIP_CALL(SCIPallocMemoryArray(scip, &nmaporigtocopies, nvars));

   /* initialize for every variable the list of copies */
   for( i = 0; i < nvars; ++i )
   {
      SCIP_CALL(SCIPcreateIntarray(scip, &maporigtocopies[i]));
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
      vars = SCIPgetVarsXXX( scip, conss[i] );

      /* TODO: skip all variables that have a zero coeffient or where all coefficients add to zero */
      /* TODO: Do more then one entry per variable actually work? */

      /* allocate a hyperedge for the constraint */
      SCIP_CALL(SCIPallocMemory(scip, &hedge));
      hedge->type = CONSTRAINT;
      hedge->originalId = i;

      /* compute its weight*/
      SCIP_CALL(computeHyperedgeWeight(scip, detectordata, conss[i], &(hedge->cost)));

      /* lets collect the variable ids of the variables */
      SCIP_CALL(SCIPallocMemoryArray(scip, &varids, ncurvars));
      hedge->nvariableIds = 0;

      for( j = 0; j < ncurvars; ++j )
      {
         SCIP_VAR* var;
         int varIndex;

         /* if the variable is inactive, skip it */
         if( !isVarRelevant(vars[j]) )
         {
            continue;
         }
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
         SCIP_CALL(SCIPsetIntarrayVal(scip, maporigtocopies[varIndex], nmaporigtocopies[varIndex], id));
         ++(nmaporigtocopies[varIndex]);
         SCIP_CALL(SCIPsetIntarrayVal(scip, copytoorig, id, varIndex));
         ++id;
         /* Check the mapping here */
#ifdef SCIP_DEBUG
         {
            int k;
            SCIPdebugPrintf("Cons: %d: ", i);
            SCIPdebugPrintf("Var: %d: ", varIndex);
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
         SCIP_CALL(SCIPduplicateMemoryArray(scip, &hedge->variableIds, varids, hedge->nvariableIds));
         SCIP_CALL(SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge));
         ++nhyperedges;
      }
      else
      {
         SCIPfreeMemory(scip, &hedge);
      }
      SCIPfreeMemoryArray(scip, &varids);
      SCIPfreeMemoryArray(scip, &vars);
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

      SCIP_CALL(SCIPallocMemory(scip, &hedge));
      hedge->type = VARIABLE;
      hedge->originalId = i;
      hedge->nvariableIds = 0;
      hedge->variableIds = NULL;
      hedge->cost = varWeight;
      SCIP_CALL(SCIPsetPtrarrayVal(scip, hedges, nhyperedges, hedge));
      ++nhyperedges;

      /* link the copies together */
      SCIP_CALL(SCIPallocMemoryArray(scip, &hedge->variableIds, size));
      hedge->nvariableIds = size;

      for( j = 0; j < size; ++j )
      {
         hedge->variableIds[j] = SCIPgetIntarrayVal(scip, maporigtocopies[i], j);
         SCIPdebugPrintf("%d, ", hedge->variableIds[j]+1);
      }
      SCIPdebugPrintf("\n");
      SCIP_CALL(SCIPfreeIntarray(scip, &maporigtocopies[i]));
   }
   SCIPfreeMemoryArray(scip, &maporigtocopies);
   SCIPfreeMemoryArray(scip, &nmaporigtocopies);
   SCIPfreeMemoryArray(scip, &copies);
   nvertices = id;
   detectordata->nvertices = nvertices;

   assert(SCIPgetPtrarrayMaxIdx(scip, hedges)+1 == nhyperedges);
   assert(nvertices == SCIPgetIntarrayMaxIdx(scip, copytoorig)+1);


   return SCIP_OKAY;
}

/** Will call hmetis via a system call */
static
SCIP_RETCODE callMetis(
      SCIP*              scip,          /**< SCIP data struture */
      DEC_DETECTORDATA*  detectordata,  /**< presolver data data structure */
      SCIP_RESULT*       result         /**< result indicating whether the detection was successful */
      )
{
   char metiscall[SCIP_MAXSTRLEN];
   char metisout[SCIP_MAXSTRLEN];
   char line[SCIP_MAXSTRLEN];
   char tempfile[SCIP_MAXSTRLEN];

   int status;
   int i;
   int j;
   int nvertices;
   int ndummyvertices;
   int* partition;

   SCIP_PTRARRAY *hedges;
   SCIP_FILE *zfile;
   FILE* file;
   int temp_filedes = -1;
   SCIP_Real remainingtime;

   assert(scip != NULL);
   assert(detectordata != NULL);

   *result = SCIP_DIDNOTRUN;

   remainingtime = DECgetRemainingTime(scip);

   if( remainingtime <= 0 )
   {
      return SCIP_OKAY;
   }

   for( i = 0; i < SCIPgetNVars(scip); ++i )
   {
      detectordata->varpart[i] = -1;
   }

   hedges = detectordata->hedges;
   nvertices = detectordata->nvertices;
   ndummyvertices = SCIPceil(scip, detectordata->dummynodes*nvertices);

   SCIPsnprintf(tempfile, SCIP_MAXSTRLEN, "gcg-metis-XXXXXX");
   if( (temp_filedes = mkstemp(tempfile)) < 0 )
   {
      SCIPerrorMessage("Error creating temporary file: %s\n", strerror( errno ));
      return SCIP_FILECREATEERROR;
   }

   SCIPdebugMessage("Temporary filename: %s\n", tempfile);

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
      SCIPerrorMessage("Could not close '%s'\n", tempfile);
      return SCIP_WRITEERROR;
   }

   /* call metis via syscall as there is no library usable ... */
   if( !SCIPisInfinity(scip, DECgetRemainingTime(scip)) )
   {
      SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"ulimit -t %.0f;hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               remainingtime,
               tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }
   else
   {
      SCIPsnprintf(metiscall, SCIP_MAXSTRLEN, "zsh -c \"hmetis %s %d -seed %d -ptype %s -ufactor %f %s\"",
               tempfile,
               detectordata->blocks,
               detectordata->randomseed,
               detectordata->metisuseptyperb ? "rb" : "kway",
               detectordata->metisubfactor,
               detectordata->metisverbose ? "" : "> /dev/null" );
   }

   SCIP_CALL(SCIPresetClock(scip, detectordata->metisclock));
   SCIP_CALL(SCIPstartClock(scip, detectordata->metisclock));
   SCIPdebugMessage("Calling metis with: %s\n", metiscall);
   status = system( metiscall );

   SCIP_CALL(SCIPstopClock(scip, detectordata->metisclock));
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
      if( detectordata->tidy )
      {
         status = unlink( tempfile );
         if( status == -1 )
         {
            SCIPerrorMessage("Could not remove metis input file: ", strerror( errno ));
            return SCIP_WRITEERROR;
         }
      }
      return SCIP_ERROR;
   }

   /*
    * parse the output into the vector
    * alloc the memory
    */
   if( detectordata->partition == NULL )
   {
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->partition, nvertices));
   }

   assert(detectordata->partition != NULL);
   partition = detectordata->partition;

   SCIPsnprintf(metisout, SCIP_MAXSTRLEN, "%s.part.%d", tempfile, detectordata->blocks);

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
   SCIPfclose(zfile);

   /* if desired delete the temoprary metis file */
   if( detectordata->tidy )
   {
      status = unlink( tempfile );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis input file: %s\n", strerror( errno ));
         return SCIP_WRITEERROR;
      }
      status = unlink( metisout );
      if( status == -1 )
      {
         SCIPerrorMessage("Could not remove metis output file: %s\n", strerror( errno ));
         return SCIP_WRITEERROR;
      }
   }
   else
   {
      SCIPinfoMessage(scip, NULL, "Temporary file is in: %s\n", tempfile);
   }
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** Maps the partitions for the disaggregated vertices to the original vertices */
static
SCIP_RETCODE assignBlocksToOriginalVariables(
      SCIP*             scip,          /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata   /**< presolver data data structure */
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
      else if(origpart[originalId] >= 0)
      {
         if(origpart[originalId] != partition[i])
         {
            origpart[originalId] = -2;
         }
      }
      assert(origpart[originalId] == -2 || origpart[originalId] >= 0);
      assert(origpart[originalId] <= detectordata->blocks);
   }

   return SCIP_OKAY;
}

/** Builds the transformed problem in the new scip instance */
static SCIP_RETCODE buildTransformedProblem(
   SCIP*                    scip,           /**< SCIP data structure */
   DEC_DETECTORDATA*        detectordata,   /**< presolver data data structure */
   DECDECOMP*               decdecomp,      /**< decdecomp data structure */
   int                      nblocks,        /**< number of blocks for this decomposition */
   SCIP_ARROWHEURSCORES*    score,           /**< scores */
   SCIP_RESULT*             result
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
   assert(score != NULL);

   nconss = SCIPgetNConss( scip );
   nvars = SCIPgetNVars( scip );

   conss = SCIPgetConss( scip );
   vars = SCIPgetVars( scip );

   SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars, nblocks) );

   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipconss, nblocks) );
   SCIP_CALL( SCIPallocBufferArray(scip, &nsubscipvars, nblocks) );

   for (i = 0; i < nblocks; ++i )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipconss[i], nconss) );
      SCIP_CALL( SCIPallocBufferArray(scip, &subscipvars[i], nvars) );

      nsubscipconss[i] = 0;
      nsubscipvars[i] = 0;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &linkingconss, nconss) );
   nlinkingconss = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &linkingvars, nconss) );
   nlinkingvars = 0;

   SCIP_CALL( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), nconss) );

   score->minkequicutscore = 0;
   score->equicutscorenormalized = 0;

   SCIP_CALL(SCIPallocMemoryArray(scip, &isVarHandled, nvars));
   for( i = 0; i < nvars; ++i )
   {
      isVarHandled[i] = FALSE;
   }

   /* go through all of the constraints */
   for( i = 0; i < nconss; i++ )
   {
      long int consblock = -1;

      /* sort the variables into corresponding buckets */
      int ncurvars = SCIPgetNVarsXXX( scip, conss[i] );
      SCIP_VAR **curvars = SCIPgetVarsXXX( scip, conss[i] );
      for( j = 0; j < ncurvars; j++ )
      {
         SCIP_VAR* var;
         long int varblock = -1;
         if( !isVarRelevant(curvars[j]) )
         {
//            SCIPprintVar(scip, curvars[j], NULL);
            continue;
         }
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
            SCIP_CALL(SCIPhashmapInsert(vartoblock, var, (void*)varblock));
         }
         else
         {
            varblock = (long int)SCIPhashmapGetImage(vartoblock, var);
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
      SCIPfreeMemoryArrayNull(scip, &curvars);

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
         SCIP_CALL(SCIPhashmapInsert(constoblock, conss[i], (void*)(block)));

      }
      /* otherwise put it in its block */
      else
      {
         subscipconss[consblock][nsubscipconss[consblock]] = conss[i];
         assert(!SCIPhashmapExists(constoblock, conss[i]));
         SCIP_CALL(SCIPhashmapInsert(constoblock, conss[i], (void*)consblock));
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
      if( detectordata->varpart[i] < 0 )
      {
         switch(SCIPvarGetType(SCIPvarGetProbvar(vars[i])))
         {
         case SCIP_VARTYPE_BINARY:
            score->minkequicutscore += detectordata->varWeightBinary;
            break;
         case SCIP_VARTYPE_CONTINUOUS:
            score->minkequicutscore += detectordata->varWeightContinous;
            break;
         case SCIP_VARTYPE_IMPLINT:
            score->minkequicutscore += detectordata->varWeightImplint;
            break;
         case SCIP_VARTYPE_INTEGER:
            score->minkequicutscore += detectordata->varWeightInteger;
            break;
         default:
            break;
         }
      }
      if( isVarHandled[i] )
      {
         continue;

      }

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

   SCIPfreeMemoryArray(scip, &isVarHandled);

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
      DECdecdecompSetNBlocks(decdecomp, nblocks);
      DECdecdecompSetType(decdecomp, DEC_DECTYPE_BORDERED);
      SCIP_CALL( DECdecdecompSetSubscipvars(scip, decdecomp, subscipvars, nsubscipvars) );
      SCIP_CALL( DECdecdecompSetSubscipconss(scip, decdecomp, subscipconss, nsubscipconss) );
      SCIP_CALL( DECdecdecompSetLinkingconss(scip, decdecomp, linkingconss, nlinkingconss) );
      if( nlinkingvars > 0 )
         SCIP_CALL( DECdecdecompSetLinkingvars(scip, decdecomp, linkingvars, nlinkingvars) );
      DECdecdecompSetVartoblock(decdecomp, vartoblock);
      DECdecdecompSetConstoblock(decdecomp, constoblock);
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

   *result = emptyblocks? SCIP_DIDNOTFIND:SCIP_SUCCESS;
   return SCIP_OKAY;
}

static
SCIP_RETCODE evaluateDecomposition(
      SCIP*                 scip,           /**< SCIP data structure */
      DEC_DETECTORDATA*     detectordata,   /**< presolver data data structure */
      DECDECOMP*            decdecomp,      /**< decomposition data structure */
      SCIP_ARROWHEURSCORES* score           /**< returns the score of the decomposition */
      )
{
   char name[SCIP_MAXSTRLEN];
   long int matrixarea;
   long int borderarea;
   int nvars;
   int nconss;
   int i;
   int j;
   int k;
   /*   int blockarea; */
   SCIP_Real varratio;
   int* nzblocks;
   int nblocks;
   int* nlinkvarsblocks;
   int* nvarsblocks;
   SCIP_Real* blockdensities;
   int* blocksizes;
   SCIP_Real density;

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(score != NULL);

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);

   nblocks = DECdecdecompGetNBlocks(decdecomp);

   /* get the right name */
   SCIPsnprintf(name, SCIP_MAXSTRLEN,
         DWSOLVER_REFNAME(SCIPgetProbName(scip),
                          nblocks,
                          detectordata->varWeightContinous,
                          detectordata->varWeightInteger,
            detectordata->consWeight,
                          detectordata->dummynodes,
                          detectordata->alpha,
                          detectordata->beta,
                          detectordata->consWeightSetppc)
      );

   SCIP_CALL(SCIPallocMemoryArray(scip, &nzblocks, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nlinkvarsblocks, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blockdensities, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &blocksizes, nblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &nvarsblocks, nblocks));
   /*
    * 3 Scores
    *
    * - Area percentage (min)
    * - block density (max)
    * - \pi_b {v_b|v_b is linking}/#vb (min)
    */

   /* calculate matrix area */
   matrixarea = nvars*nconss;

   /* calculate slave sizes, nonzeros and linkingvars */
   for( i = 0; i < nblocks; ++i )
   {
      SCIP_CONS** curconss;
      int ncurconss;
      int nvarsblock;
      SCIP_Bool *ishandled;

      SCIP_CALL(SCIPallocMemoryArray(scip, &ishandled, nvars));
      nvarsblock = 0;
      nzblocks[i] = 0;
      nlinkvarsblocks[i] = 0;
      for( j = 0; j < nvars; ++j )
      {
         ishandled[j] = FALSE;
      }
      curconss = DECdecdecompGetSubscipconss(decdecomp)[i];
      ncurconss = DECdecdecompGetNSubscipconss(decdecomp)[i];

      for( j = 0; j < ncurconss; ++j )
      {
         SCIP_VAR** curvars;
         SCIP_VAR* var;
         int ncurvars;

         curvars = SCIPgetVarsXXX(scip, curconss[j]);
         ncurvars = SCIPgetNVarsXXX(scip, curconss[j]);
         for( k = 0; k < ncurvars; ++k )
         {
            long int block;
            if( !isVarRelevant(curvars[k]) )
               continue;

            var = SCIPvarGetProbvar(curvars[k]);
            assert(var != NULL);
            assert(SCIPvarIsActive(var));
            assert(!SCIPvarIsDeleted(var));
            ++(nzblocks[i]);
            assert(SCIPhashmapExists(DECdecdecompGetVartoblock(decdecomp), var));
            block = (long int) SCIPhashmapGetImage(DECdecdecompGetVartoblock(decdecomp), var);

            if(block == detectordata->blocks+1 && ishandled[SCIPvarGetProbindex(var)] == FALSE)
            {
               ++(nlinkvarsblocks[i]);
            }
            ishandled[SCIPvarGetProbindex(var)] = TRUE;
         }

         SCIPfreeMemoryArray(scip, &curvars);
      }

      for( j = 0; j < nvars; ++j )
      {
         if( ishandled[j] )
         {
            ++nvarsblock;
         }
      }

      blocksizes[i] = nvarsblock*ncurconss;
      nvarsblocks[i] = nvarsblock;
      if(blocksizes[i] > 0)
      {
         blockdensities[i] = 1.0*nzblocks[i]/blocksizes[i];
      }
      else
      {
         blockdensities[i] = 0.0;
      }

      assert(blockdensities[i] >= 0 && blockdensities[i] <= 1.0);
      SCIPfreeMemoryArray(scip, &ishandled);
   }

   /* calculate border area */
   borderarea = DECdecdecompGetNLinkingconss(decdecomp)*nvars+DECdecdecompGetNLinkingvars(decdecomp)*(nconss-DECdecdecompGetNLinkingconss(decdecomp));

   /*   blockarea = 0; */
   density = 1E20;
   varratio = 1.0;
   for( i = 0; i < nblocks; ++i )
   {
      /* calculate block area */
      /* blockarea += blocksizes[i]; */


      /* calculate density */
      density = MIN(density, blockdensities[i]);

      /* calculate linking var ratio */
      if( DECdecdecompGetNLinkingvars(decdecomp) > 0 )
      {
         varratio *= 1.0*nlinkvarsblocks[i]/DECdecdecompGetNLinkingvars(decdecomp);
      }
      else
      {
         varratio = 0;
      }
   }

   score->linkingscore = (0.5+0.5*varratio);
   score->borderscore = (1.0*(borderarea)/matrixarea);
   score->densityscore = (1-density);

   SCIPfreeMemoryArray(scip, &nzblocks);
   SCIPfreeMemoryArray(scip, &nlinkvarsblocks);
   SCIPfreeMemoryArray(scip, &blockdensities);
   SCIPfreeMemoryArray(scip, &blocksizes);
   SCIPfreeMemoryArray(scip, &nvarsblocks);

   return SCIP_OKAY;

}


static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildArrowhead)
{

   SCIP_ARROWHEURSCORES* scores;
   SCIP_Real* cumscores;
   int i;
   int j;
   int ndecs;

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
   SCIP_CALL( SCIPallocMemoryArray(scip, &scores, ndecs) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &cumscores, ndecs) );

   /* build the hypergraph structure from the original problem */
   SCIP_CALL(buildGraphStructure(scip, detectordata));

   for(i = 0; i < ndecs; ++i)
   {
      SCIP_CALL_ABORT( DECdecdecompCreate(scip, &(*decdecomps)[i]) );
   }


   for( j = 0, i = detectordata->minblocks; i <= detectordata->maxblocks; ++i )
   {
      detectordata->blocks = i;
      /* get the partitions for the new variables from metis */
      SCIP_CALL(callMetis(scip, detectordata, result));

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

      SCIP_CALL( buildTransformedProblem(scip, detectordata, (*decdecomps)[j], i, &scores[j], result) );
      if( *result == SCIP_SUCCESS )
      {
         SCIP_CALL( evaluateDecomposition(scip, detectordata, (*decdecomps)[j], &scores[j]) );
         SCIP_CALL( printArrowheurScores(scip, detectordata, &scores[j]) );

         cumscores[j] = scores[j].borderscore*scores[j].linkingscore*scores[j].densityscore;
         *ndecdecomps += 1;
         ++j;
      }
   }

   SCIPsortRealPtr(cumscores, (void**) *decdecomps, *ndecdecomps);

   for( i = *ndecdecomps; i < ndecs; ++i )
   {
      DECdecdecompFree(scip, &((*decdecomps)[i]) );
   }

   SCIP_CALL(SCIPreallocMemoryArray(scip, decdecomps, *ndecdecomps));
   SCIPfreeMemoryArray(scip, &cumscores);
   SCIPfreeMemoryArray(scip, &scores);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** returns the priority of the detector */
static
DEC_DECL_GETPRIORITY(getPriority)
{
   DEC_DETECTOR* arrowheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   arrowheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(arrowheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(arrowheur), DEC_DETECTORNAME) == 0);
   return detectordata->priority;
}

/** creates the arrowheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionArrowheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL(SCIPallocMemory(scip, &detectordata));

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->partition = NULL;
   detectordata->blocks = -1;

   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildArrowhead, initArrowheur, exitArrowheur, getPriority));


   /* add arrowheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "arrowheur/beta", "factor on how heavy equality (beta) and inequality constraints are measured", &detectordata->beta, FALSE, DEFAULT_BETA, 0.0, 1.0, NULL, NULL ));
   SCIP_CALL(SCIPaddRealParam(scip, "arrowheur/alpha", "factor on how heavy the standard deviation of the coefficients is measured", &detectordata->alpha, FALSE, DEFAULT_ALPHA, 0.0, 1E20, NULL, NULL ));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/varWeight", "Weight of a variable hyperedge", &detectordata->varWeight, FALSE, DEFAULT_VARWEIGHT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/varWeightBinary", "Weight of a binary variable hyperedge", &detectordata->varWeightBinary, FALSE, DEFAULT_VARWEIGHTBIN, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/varWeightContinous", "Weight of a continuos variable hyperedge", &detectordata->varWeightContinous, FALSE, DEFAULT_VARWEIGHTCONT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/varWeightImplint", "Weight of a implicit integer variable hyperedge", &detectordata->varWeightImplint, FALSE, DEFAULT_VARWEIGHTIMPL, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/varWeightInteger", "Weight of a integer variable hyperedge", &detectordata->varWeightInteger, FALSE, DEFAULT_VARWEIGHTINT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/consWeight", "Weight of a constraint hyperedge", &detectordata->consWeight, FALSE, DEFAULT_CONSWEIGHT, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddBoolParam(scip, "arrowheur/tidy", "Whether to clean up temporary files", &detectordata->tidy, FALSE, DEFAULT_TIDY, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/randomseed", "random seed for hmetis", &detectordata->randomseed, FALSE, DEFAULT_RANDSEED, -1, INT_MAX, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "arrowheur/dummynodes", "percentage of dummy nodes for metis", &detectordata->dummynodes, FALSE, DEFAULT_DUMMYNODES, 0.0, 1.0, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/consWeightSetppc", "Weight for constraint hyperedges that are setpartitioning or covering constraints", &detectordata->consWeightSetppc, FALSE, DEFAULT_CONSWEIGHT_SETPPC, 0, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddRealParam(scip, "arrowheur/ubfactor", "Unbalance factor for metis", &detectordata->metisubfactor, FALSE, DEFAULT_METIS_UBFACTOR, 0.0, 1E20, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "arrowheur/metisverbose", "Should the metis output be displayed", &detectordata->metisverbose, FALSE, DEFAULT_METIS_VERBOSE, NULL, NULL ));
   SCIP_CALL(SCIPaddBoolParam(scip, "arrowheur/metisuseptyperb", "Should the rb or kway method be used for partitioning by metis", &detectordata->metisuseptyperb, FALSE, DEFAULT_METISUSEPTYPE_RB, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "arrowheur/priority", "random seed for hmetis", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));

   return SCIP_OKAY;
}
