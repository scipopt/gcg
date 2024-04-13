/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2023 Operations Research, RWTH Aachen University       */
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

/**@file    sepa_xyz.c
 *
 * @brief   xyz separator for master problem (put your description here)
 * @author  Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
//#define SCIP_DEBUG
#include <assert.h>

#include "mastercutdata.h"
#include "gcg.h"
#include "pricer_gcg.h"
#include "relax_gcg.h"
#include "sepa_subsetrow.h"
#include "struct_mastercutdata.h"
#include "struct_sepagcg.h"
#include "type_mastercutdata.h"
#include "type_sepagcg.h"




#define SEPA_NAME              "subsetrow"
#define SEPA_DESC              "subsetrow separator"
#define SEPA_PRIORITY                 0
#define SEPA_FREQ                    10
#define SEPA_MAXBOUNDDIST           1.0
#define SEPA_USESSUBSCIP          FALSE /**< does the separator use a secondary SCIP instance? */
#define SEPA_DELAY                FALSE /**< should separation method be delayed, if other separators found cuts? */
#define STARTMAXCUTS                50

/*
 * Data structures
 */

/* TODO: fill in the necessary separator data */

/** separator data */
struct SCIP_SepaData
{
   GCG_MASTERCUTDATA**     mastercuts;
   GCG_MASTERCUTDATA**     newcuts;
   SCIP_HASHMAP*           rowidxmap;
   int                     nmastercuts;
   int                     nnewcuts;
   int                     newcutssize;
   int                     mastercutssize;
};


/*
 * Local methods
 */
/** allocates enough memory to hold more cuts */
static
SCIP_RETCODE ensureSizeMastercuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   int                   size                /**< new size of cut arrays */
)
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->mastercuts != NULL);
   assert(sepadata->nmastercuts <= sepadata->mastercutssize);
   assert(sepadata->nmastercuts >= 0);

   if( sepadata->mastercutssize < size )
   {
      int newmaxcuts = SCIPcalcMemGrowSize(scip, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->mastercuts), sepadata->mastercutssize, newmaxcuts) );
      sepadata->mastercutssize = newmaxcuts;
   }
   assert(sepadata->mastercutssize >= size);

   return SCIP_OKAY;
}

/** allocates enough memory to hold more cuts */
static
SCIP_RETCODE ensureSizeNewcuts(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SEPADATA*        sepadata,           /**< separator data structure */
   int                   size                /**< new size of cut arrays */
)
{
   assert(scip != NULL);
   assert(sepadata != NULL);
   assert(sepadata->newcuts != NULL);
   assert(sepadata->nnewcuts <= sepadata->newcutssize);
   assert(sepadata->nnewcuts >= 0);

   if( sepadata->newcutssize < size )
   {
      int newmaxcuts = SCIPcalcMemGrowSize(scip, size);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &(sepadata->newcuts), sepadata->newcutssize, newmaxcuts) );
      sepadata->newcutssize = newmaxcuts;
   }
   assert(sepadata->newcutssize >= size);

   return SCIP_OKAY;
}


/*
 * Callback methods of separator
 */

/* TODO: Implement all necessary separator methods. The methods with an #if 0 ... #else #define ... are optional */


#define sepaCopySubsetrow NULL
#define sepaExitSubsetrow NULL
#define sepaExecsolSubsetrow NULL

static
SCIP_DECL_SEPAEXITSOL(sepaExitsolSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* row;
   int i;

   assert(scip != NULL);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa);
   assert(sepadata != NULL);

   SCIPdebugMessage("clear gcg sepa subsetrow data\n");
   /* release all the master cuts and empty the hashmap tracking the indices of the cuts */
   for( i = 0; i < sepadata->nmastercuts; i++ )
   {
      row = GCGmastercutGetRow(sepadata->mastercuts[i]);
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   sepadata->nmastercuts = 0;
   SCIPhashmapRemoveAll(sepadata->rowidxmap);

   for( i = 0; i < sepadata->nnewcuts; i++ )
   {
      row = GCGmastercutGetRow(sepadata->newcuts[i]);
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }
   sepadata->nnewcuts = 0;

   return SCIP_OKAY;
}


/** destructor of separator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_SEPAFREE(sepaFreeSubsetrow)
{  /*lint --e{715}*/

   SCIP_SEPADATA* sepadata;

   assert(scip != NULL);
   assert(sepa != NULL);

   SCIPdebugMessage("free gcg sepa subsetrow sepa data\n");
   sepadata = SCIPsepaGetData(sepa);
   SCIPfreeBlockMemoryArray(scip, &(sepadata->mastercuts), sepadata->mastercutssize);
   SCIPfreeBlockMemoryArray(scip, &(sepadata->newcuts), sepadata->newcutssize);
   SCIPhashmapFree(&sepadata->rowidxmap);
   SCIPfreeBlockMemory(scip, &sepadata);

   return SCIP_OKAY;
}


/** LP solution separation method of separator */
static
SCIP_DECL_SEPAEXECLP(sepaExeclpSubsetrow)
{  /*lint --e{715}*/
   SCIPdebugMessage("sepaExeclpSubsetrow\n");

   return SCIP_OKAY;
}


/*
 * Callback methods of GCG separator
 */


/** method for determining if cut was generated by this separator */
static
GCG_DECL_SEPAGETMASTERCUTDATA(gcgsepaGetMastercutdataSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   int rowindex;
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(row != NULL);
   assert(GCGisMaster(scip));

   sepadata = SCIPsepaGetData(sepa->separator);
   assert(sepadata != NULL);
   if( SCIPhashmapExists(sepadata->rowidxmap, row) )
   {
      rowindex = SCIPhashmapGetImageInt(sepadata->rowidxmap, row);
      *result = sepadata->mastercuts[rowindex];
      assert( GCGmastercutGetRow(*result) == row );
      SCIPdebugMessage("Found mastercutdata for row %s\n", SCIProwGetName(row));
      return SCIP_OKAY;
   }

   /* if for some reason the corresponding mastercutdata was not found: return NULL */
   SCIPdebugMessage("Did not find mastercutdata for row %s\n", SCIProwGetName(row));
   *result = NULL;
   return SCIP_OKAY;
}


/** method for making newly generated cuts available */
static
GCG_DECL_SEPAUPDATENEWCUTS(gcgsepaUpdateNewCutsSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   SCIP_ROW* row;
   int i;
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));

   sepadata = SCIPsepaGetData(sepa->separator);
   assert(sepadata != NULL);
   SCIPdebugMessage("%i cuts currently in newcuts\n", sepadata->nnewcuts);
   /* transfer the cuts which were actually added to the LP, to the master cuts array */
   for( i = 0; i < sepadata->nnewcuts; i++ )
   {
      row = GCGmastercutGetRow(sepadata->newcuts[i]);
      if ( GCGmastercutIsActive(sepadata->newcuts[i]) )
      {
         ensureSizeMastercuts(scip, sepadata, sepadata->nmastercuts + 1);
         sepadata->mastercuts[sepadata->nmastercuts] = sepadata->newcuts[i];
         SCIP_CALL( SCIPhashmapInsertInt(sepadata->rowidxmap, row, sepadata->nmastercuts) );
         ++(sepadata->nmastercuts);
         SCIPcaptureRow(scip, row);
         SCIPdebugMessage("cut %s was transferred to mastercuts\n", SCIProwGetName(row));
      }
      SCIPreleaseRow(scip, &row);
   }

   /* after all relevant cuts were transferred, we can clear the newcuts array */
   sepadata->nnewcuts = 0;
   SCIPdebugMessage("cleared newcuts after transferring active cuts\n");

   return SCIP_OKAY;
}


/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPAGETCUTS(gcgsepaGetCutsSubsetrow)
{
   SCIP_SEPADATA* sepadata;
   assert(scip != NULL);
   assert(sepa != NULL);

   sepadata = SCIPsepaGetData(sepa->separator);
   assert(sepadata != NULL);
   *cuts = sepadata->mastercuts;
   *ncuts = sepadata->nmastercuts;

   return SCIP_OKAY;
}


/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPAGETCOLCOEFFICIENT(gcgsepaGetColCoefficientSubsetrow)
{
   SCIPdebugMessage("gcgsepaGetColCoefficientSubsetrow\n");
   assert(gcgcol != NULL);
   return SCIP_OKAY;
}

/** method for adding new master variable to cut */
static
GCG_DECL_SEPAGETVARCOEFFICIENT(gcgsepaGetVarCoefficientSubsetrow)
{
   SCIPdebugMessage("gcgsepaGetColCoefficientSubsetrow\n");
   return SCIP_OKAY;
}

/** get all the cuts generated by this master separator */
static
GCG_DECL_SEPASETOBJECTIVE(gcgsepaSetObjectiveSubsetrow)
{
   SCIPdebugMessage("gcgsepaSetObjectiveSubsetrow\n");

   return SCIP_OKAY;
}


/*
 * MASTER separator specific interface methods
 */


/** initialization method of separator (called after problem was transformed) */
static
SCIP_DECL_SEPAINIT(sepaInitSubsetrow)
{
   SCIP* origscip;
   assert(scip != NULL);
   assert(sepa != NULL);
   assert(GCGisMaster(scip));
   origscip = GCGmasterGetOrigprob(scip);
   assert(origscip != NULL);

   SCIPdebugMessage("initialize gcg sepa subsetrow\n");
   /* creates the subsetrow gcg separator and includes in the relaxator data of the original problem */
   SCIP_CALL(GCGrelaxIncludeSeparator(origscip, sepa, gcgsepaGetMastercutdataSubsetrow, gcgsepaUpdateNewCutsSubsetrow,
                                      gcgsepaGetCutsSubsetrow, gcgsepaGetColCoefficientSubsetrow, gcgsepaSetObjectiveSubsetrow));

   return SCIP_OKAY;
}


/*
 * separator specific interface methods
 */

/** creates the scip sepa of the subsetrow separator and includes it in SCIP */
SCIP_RETCODE SCIPincludeSepaSubsetrow(
   SCIP*                 scip                /**< SCIP data structure */
)
{
   SCIP_SEPADATA* sepadata;
   SCIP_SEPA* sepa;

   /* create subsetrow separator data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &sepadata) );
   sepadata->nmastercuts = 0;
   sepadata->nnewcuts = 0;
   sepadata->mastercutssize = SCIPcalcMemGrowSize(scip, STARTMAXCUTS);
   sepadata->newcutssize = SCIPcalcMemGrowSize(scip, STARTMAXCUTS);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->mastercuts), sepadata->mastercutssize) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sepadata->newcuts), sepadata->newcutssize) );
   SCIP_CALL( SCIPhashmapCreate(&sepadata->rowidxmap, SCIPblkmem(scip), sepadata->mastercutssize) );
   sepa = NULL;
   SCIP_CALL( SCIPincludeSepaBasic(scip, &sepa, SEPA_NAME, SEPA_DESC, SEPA_PRIORITY, SEPA_FREQ, SEPA_MAXBOUNDDIST,
                                   SEPA_USESSUBSCIP, SEPA_DELAY,
                                   sepaExeclpSubsetrow, sepaExecsolSubsetrow,
                                   sepadata) );
   assert(sepa != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetSepaFree(scip, sepa, sepaFreeSubsetrow) );
   SCIP_CALL( SCIPsetSepaInit(scip, sepa, sepaInitSubsetrow) );
   //SCIP_CALL( SCIPsetSepaExit(scip, sepa, sepaExitSubsetrow) );
   SCIP_CALL( SCIPsetSepaExitsol(scip, sepa, sepaExitsolSubsetrow) );

   return SCIP_OKAY;
}