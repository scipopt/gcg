/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2024 Operations Research, RWTH Aachen University       */
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

/**@file   mastersepacut.c
 * @brief
 * @author Chantal Reinartz Groba
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include <scip/type_scip.h>

#include "gcg.h"
#include "mastersepacut.h"
#include "struct_mastersepacutdata.h"


/** frees data of subset row cut */
static
SCIP_RETCODE freeSubsetRowCutData(
   SCIP*                      masterscip,    /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUTDATA**    data           /**< pointer to data of subset row cut */
   )
{
   if( *data == NULL )
      return SCIP_OKAY;

   if( (*data)->data.subsetrowcutdata.n > 0 )
   {
      SCIPfreeBlockMemoryArrayNull(masterscip, &((*data)->data.subsetrowcutdata.weights), (*data)->data.subsetrowcutdata.n);
      SCIPfreeBlockMemoryArrayNull(masterscip, &((*data)->data.subsetrowcutdata.conssindices), (*data)->data.subsetrowcutdata.n);
      (*data)->data.subsetrowcutdata.n = 0;
   }

   assert((*data)->data.subsetrowcutdata.weights == NULL);
   assert((*data)->data.subsetrowcutdata.conssindices == NULL);
   assert((*data)->data.subsetrowcutdata.n == 0);

   SCIPfreeBlockMemory(masterscip, data);
   *data = NULL;

   return SCIP_OKAY;
}

/** frees master separator cut */
static
SCIP_RETCODE freeMasterSepaCut(
   SCIP*                masterscip,       /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUT**  mastersepacut     /**< pointer to master separator cut */
)
{
   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(mastersepacut != NULL);
   assert(*mastersepacut != NULL);

#ifdef SCIP_DEBUG
   SCIP_ROW* row;
   SCIP_CALL( GCGmastercutGetRow((*mastersepacut)->mastercutdata, &row) );
   SCIPdebugMessage("free master separator cut: free cut for row %s\n", SCIProwGetName(row));
#endif

   if( (*mastersepacut)->knownvarhistory != NULL )
   {
      SCIPdebugMessage("free mastersepacut: var history is freed\n");
      SCIP_CALL( GCGvarhistoryFreeReference(masterscip, &((*mastersepacut)->knownvarhistory)) );
   }
   assert((*mastersepacut)->knownvarhistory == NULL);

   if( (*mastersepacut)->mastercutdata != NULL )
   {
      SCIP_CALL( GCGmastercutFree(masterscip, &((*mastersepacut)->mastercutdata)) );
   }
   assert((*mastersepacut)->mastercutdata == NULL);

   if( (*mastersepacut)->cuttype == GCG_MASTERSEPACUTTYPE_SUBSETROW )
   {
      SCIP_CALL( freeSubsetRowCutData(masterscip, &((*mastersepacut)->data)) );
   }
   assert((*mastersepacut)->data == NULL);

   SCIPfreeBlockMemory(masterscip, mastersepacut);
   *mastersepacut = NULL;

   return SCIP_OKAY;
}

/** increases usage counter of master separator cut */
SCIP_RETCODE GCGcaptureMasterSepaCut(
   GCG_MASTERSEPACUT*      mastersepacut      /**< master separator cut */
)
{
   assert(mastersepacut != NULL);
   assert(mastersepacut->nuses >= 0);

#ifdef SCIP_DEBUG
   SCIP_ROW* row;
   SCIP_CALL( GCGmastercutGetRow(mastersepacut->mastercutdata, &row) );
   SCIPdebugMessage("capture master separator cut: row %s now has %i nuses\n", SCIProwGetName(row), mastersepacut->nuses);
#endif
   (mastersepacut->nuses)++;

   return SCIP_OKAY;
}

/** decreases usage counter of master separator cut, and frees memory if necessary */
SCIP_RETCODE GCGreleaseMasterSepaCut(
   SCIP*                   masterscip,      /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUT**     mastersepacut    /**< pointer to master separator cut */
)
{
   assert(masterscip != NULL);
   assert(mastersepacut != NULL);
   assert(*mastersepacut != NULL);
   assert((*mastersepacut)->nuses >= 0);

#ifdef SCIP_DEBUG
   SCIP_ROW* row;
   SCIP_CALL( GCGmastercutGetRow((*mastersepacut)->mastercutdata, &row) );
   SCIPdebugMessage("release master separator cut: row %s now has %i nuses\n", SCIProwGetName(row), (*mastersepacut)->nuses);
#endif

   ((*mastersepacut)->nuses)--;
   if( (*mastersepacut)->nuses == 0 )
   {
      SCIP_CALL( freeMasterSepaCut(masterscip, mastersepacut) );
   }

   *mastersepacut = NULL;

   return SCIP_OKAY;
}

/**< creates master separator cut */
SCIP_RETCODE GCGcreateMasterSepaCut(
   SCIP*                   masterscip,          /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUT**     mastersepacut,       /**< pointer to store master separator cut */
   GCG_MASTERSEPACUTTYPE   mastersepacuttype,   /**< type of master separator cut */
   GCG_MASTERCUTDATA*      mastercutdata,       /**< master cut data */
   GCG_VARHISTORY*         varhistory,          /**< variable history */
   GCG_MASTERSEPACUTDATA*  mastersepacutdata    /**< master separator cut data */
   )
{
   assert(masterscip != NULL);
   assert(mastercutdata != NULL);
   assert(mastersepacut != NULL);
   assert(GCGisMaster(masterscip));

   SCIP_CALL( SCIPallocBlockMemory(masterscip, mastersepacut) );
   (*mastersepacut)->mastercutdata = mastercutdata;
   (*mastersepacut)->nuses = 0;
   (*mastersepacut)->knownvarhistory = varhistory;
   (*mastersepacut)->cuttype = mastersepacuttype;
   (*mastersepacut)->data = mastersepacutdata;

   SCIP_CALL( GCGcaptureMasterSepaCut(*mastersepacut) );

   return SCIP_OKAY;
}

/**< returns the master cut data of the master separator cut */
GCG_MASTERCUTDATA* GCGmastersepacutGetMasterCutData(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
   )
{
   assert(mastersepacut != NULL);

   return mastersepacut->mastercutdata;
}

/**< returns the variable history of the master separator cut */
GCG_VARHISTORY* GCGmastersepacutGetVarHistory(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
   )
{
   assert(mastersepacut != NULL);

   return mastersepacut->knownvarhistory;
}

/**< returns the cut type of the master separator cut */
GCG_MASTERSEPACUTTYPE GCGmastersepacutGetCutType(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
)
{
   assert(mastersepacut != NULL);

   return mastersepacut->cuttype;
}

/**< returns the data of the master separator cut */
GCG_MASTERSEPACUTDATA* GCGmastersepacutGetData(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
)
{
   assert(mastersepacut != NULL);

   return mastersepacut->data;
}

/**< set the variable history of master separator cut */
SCIP_RETCODE GCGmastersepacutSetVarHistory(
   SCIP*                   masterscip,       /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUT**     mastersepacut     /**< pointer to master separator cut */
   )
{
   assert(masterscip != NULL);
   assert(mastersepacut != NULL);
   assert(*mastersepacut != NULL);

#ifdef SCIP_DEBUG
   SCIP_ROW* row;
   SCIP_CALL( GCGmastercutGetRow((*mastersepacut)->mastercutdata, &row) );
   SCIPdebugMessage("set var history: set history for row %s\n", SCIProwGetName(row));
#endif

   SCIP_CALL( GCGvarhistoryCopyReference(masterscip, &((*mastersepacut)->knownvarhistory), GCGgetCurrentVarhistoryReference(masterscip)) );
   return SCIP_OKAY;
}


// SUBSETROW CUT SPECIFIC METHODS

/**< creates a subset row cut */
SCIP_RETCODE GCGcreateSubsetRowCut(
   SCIP*                   masterscip,            /**< SCIP data structure (master problem) */
   GCG_MASTERSEPACUT**     mastersepacut,         /**< pointer to store master separator cut */
   GCG_MASTERCUTDATA*      mastercutdata,         /**< mastercutdata associated with the cut */
   GCG_VARHISTORY*         varhistory,            /**< variables history of subset row cut*/
   SCIP_Real*              weights,               /**< weights which were used to create the cut */
   int*                    indices,               /**< indices of constraints used to create the cut */
   int                     n                      /**< number of constraints used to create the cut */
   )
{
   GCG_MASTERSEPACUTDATA* data;

   assert(masterscip != NULL);
   assert(GCGisMaster(masterscip));
   assert(mastersepacut != NULL);
   assert(n >= 0);

   SCIP_CALL( SCIPallocBlockMemory(masterscip, &data) );
   data->data.subsetrowcutdata.n = n;
   data->data.subsetrowcutdata.weights = NULL;
   data->data.subsetrowcutdata.conssindices = NULL;

   if( n > 0 )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &(data->data.subsetrowcutdata.weights), weights, n) );
      SCIP_CALL( SCIPduplicateBlockMemoryArray(masterscip, &(data->data.subsetrowcutdata.conssindices), indices, n) );
   }

   SCIP_CALL( GCGcreateMasterSepaCut(masterscip, mastersepacut, GCG_MASTERSEPACUTTYPE_SUBSETROW,
                                     mastercutdata, varhistory, data) );
   return SCIP_OKAY;
}

/**< returns TRUE or FALSE whether cut is a subset row cut */
SCIP_Bool GCGmastersepacutIsSubsetRow(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
   )
{
   assert(mastersepacut != NULL);

   return mastersepacut->cuttype == GCG_MASTERSEPACUTTYPE_SUBSETROW;
}

/**< returns the number of weights of subset row cut */
int GCGsubsetrowCutGetNWeights(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
   )
{
   GCG_MASTERSEPACUTDATA* data;

   assert(mastersepacut != NULL);
   assert(GCGmastersepacutIsSubsetRow(mastersepacut));

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->data.subsetrowcutdata.n;
}

/**< returns the weights of subset row cut */
SCIP_Real* GCGsubsetrowCutGetWeights(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
)
{
   GCG_MASTERSEPACUTDATA* data;

   assert(mastersepacut != NULL);
   assert(GCGmastersepacutIsSubsetRow(mastersepacut));

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->data.subsetrowcutdata.weights;
}

/**< returns the constraint indices of subset row cut */
int* GCGsubsetrowCutGetConssIndices(
   GCG_MASTERSEPACUT*      mastersepacut     /**< master separator cut */
)
{
   GCG_MASTERSEPACUTDATA* data;

   assert(mastersepacut != NULL);
   assert(GCGmastersepacutIsSubsetRow(mastersepacut));

   data = GCGmastersepacutGetData(mastersepacut);
   assert(data != NULL);

   return data->data.subsetrowcutdata.conssindices;
}
