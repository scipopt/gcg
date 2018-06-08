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

/**@file   dec_colors.cpp
 * @ingroup DETECTORS
 * @brief  detector assigning color classes to constraints and try combinations of colors in the master
 * @author Martin Bergner
 * @todo   allow to set range of subsets
 * @todo   add parameters for min/max subsets
 * @todo   allow for a fine grained control (ignore rhs, lhs and only consider constraint handler?)
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "dec_colors.h"
#include "cons_decomp.h"
#include "scip_misc.h"
#include "pub_decomp.h"

#include <set>
#include <vector>
#include <bitset>
#include <iostream>
#include <algorithm>

/* constraint handler properties */
#define DEC_DETECTORNAME          "colors"    /**< name of detector */
#define DEC_DESC                  "Detector according to color classes" /**< description of detector*/
#define DEC_FREQCALLROUND         1           /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
#define DEC_MAXCALLROUND          INT_MAX     /** last round the detector gets called                            */
#define DEC_MINCALLROUND          0           /** first round the detector gets called                             */
#define DEC_FREQCALLROUNDORIGINAL 1           /** frequency the detector gets called in detection loop while detecting the original problem   */
#define DEC_MAXCALLROUNDORIGINAL  INT_MAX     /** last round the detector gets called while detecting the original problem                            */
#define DEC_MINCALLROUNDORIGINAL  0           /** first round the detector gets called while detecting the original problem    */
#define DEC_PRIORITY              0              /**< priority of the detector */

#define DEC_DECCHAR              'k'            /**< display character of detector */

#define DEC_ENABLED              FALSE          /**< should the detection be enabled */
#define DEC_ENABLED_ORIGINAL     FALSE          /**< should the detection of the original problem be enabled */
#define DEC_ENABLEDFINISHING     FALSE          /**< should the finishing be enabled */
#define DEC_ENABLEDPOSTPROCESSING FALSE          /**< should the finishing be enabled */
#define DEC_SKIP                 FALSE          /**< should detector be skipped if others found detections */
#define DEC_USEFULRECALL         FALSE       /**< is it useful to call this detector on a descendant of the propagated seeed */
#define DEC_LEGACYMODE           FALSE       /**< should (old) DETECTSTRUCTURE method also be used for detection */

/*
 * Data structures
 */

/** constraint handler data */
struct DEC_DetectorData
{
};


/*
 * Local methods
 */

struct ConsData {
   SCIP* scip;
   SCIP_Real lhs;
   SCIP_Real rhs;
   const char* conshdlrname;

   ConsData(SCIP* scip_, SCIP_CONS* cons) {
      scip = scip_;
      lhs = GCGconsGetLhs(scip, cons);
      rhs = GCGconsGetRhs(scip, cons);
      conshdlrname = SCIPconshdlrGetName(SCIPconsGetHdlr(cons));
   }

   void print(void)
   {
      SCIPdebugMessage("Data: %s, lhs %.3f, rhs %.3f\n", conshdlrname, lhs, rhs);
   }
};

typedef struct ConsData CONSDATA;

/* put your local methods here, and declare them static */
#ifdef SCIP_DEBUG
static
void printSubset(
   std::vector<bool> bit_mask
   )
{

   static int cnt = 0;
   std::cout << ++cnt << ". [ ";
   for( std::size_t i = 0; i < bit_mask.size(); ++i )
      if( bit_mask[i] )
         std::cout << i << ' ';
   std::cout << "]\n";
}
#endif


static
SCIP_DECL_SORTPTRCOMP(sortCons)
{
   CONSDATA* dat1 = (CONSDATA*)elem1;
   CONSDATA* dat2 = (CONSDATA*)elem2;

   int result = strncmp(dat1->conshdlrname, dat2->conshdlrname, (size_t) SCIP_MAXSTRLEN);
   if( result == 0)
   {
         if( SCIPisLT(dat1->scip, dat1->lhs, dat2->lhs) )
            result = -1;
         else if (SCIPisGT(dat1->scip, dat1->lhs, dat2->lhs) )
            result = 1;
   }
   if( result == 0)
   {
         if( SCIPisLT(dat1->scip, dat1->rhs, dat2->rhs) )
            result = -1;
         else if (SCIPisGT(dat1->scip, dat1->rhs, dat2->rhs) )
            result = 1;
   }

   //SCIPdebugMessage("1 %s 2\n", result == 0 ? "=": (result < 0 ? "<" : ">"));
   return result;
}

static
SCIP_RETCODE assignConsColors(
      SCIP*                 scip,               /**< SCIP data structure */
      SCIP_CONS**           conss,              /**< constraints to assign */
      int                   nconss,             /**< number of constraints to assign */
      int*                  colors,             /**< assignment of constrainst to colors */
      int*                  ncolors             /**< number of colors */
)
{ /*lint -esym(593, data)*/

   CONSDATA** colordata = NULL;
   int pos;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(nconss > 0);
   assert(colors != NULL);
   assert(ncolors != NULL);

   SCIP_CALL( SCIPallocMemoryArray(scip, &colordata, nconss) );
   *ncolors = 0;
   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons = conss[i];
      CONSDATA* data = new CONSDATA(scip, cons);

      if( !SCIPsortedvecFindPtr((void**)colordata, sortCons, data, *ncolors, &pos) )
      {
         SCIPsortedvecInsertPtr( (void**) colordata, sortCons, data, ncolors, &pos );
      }
      else
         delete data;
   }

   for( int i = 0; i < *ncolors; ++i)
   {
      colordata[i]->print();
   }

   for( int i = 0; i < nconss; ++i )
   {
      SCIP_CONS* cons = conss[i];
      CONSDATA* data = new CONSDATA(scip, cons);

      (void) SCIPsortedvecFindPtr( (void**) colordata, sortCons, data, *ncolors, &pos);
      colors[i] = pos;
      SCIPdebugMessage("Conss <%s> has color %d\n", SCIPconsGetName(conss[i]), pos);
      delete data;
   }

   SCIPfreeMemoryArray(scip, &colordata);
   SCIPdebugMessage("%d colors found\n", *ncolors);
   return SCIP_OKAY;
}


/** creates array for constraints in the master */
static
SCIP_RETCODE createMasterconssArray(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS***          masterconss,        /**< pointer to return masterconss array */
   int*                  nmasterconss,       /**< pointer to return number of masterconss */
   int*                  colors,             /**< array of colors for each constraint */
   std::set<int>         colorset,              /**< color to put into master */
   SCIP_Bool*            masterisempty,      /**< indicates if the master is empty (all constraints in the subproblem) */
   SCIP_Bool*            pricingisempty      /**< indicates if the pricing is empty (all constraints in the master) */
   )
{
   int i;

   SCIP_CONS** conss;
   int nconss;

   assert(scip != NULL);
   assert(masterconss != NULL);
   assert(nmasterconss != NULL);
   assert(masterisempty != NULL);
   assert(pricingisempty != NULL);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPallocMemoryArray(scip, masterconss, nconss) );
   *nmasterconss = 0;

   for (i = 0; i < nconss; ++i)
   {
      if( colorset.find(colors[i]) != colorset.end() ) /*lint !e1702 std::pair*/
      {
         SCIPdebugMessage("Constraint <%s> to be placed in master.\n", SCIPconsGetName(conss[i]));
         (*masterconss)[*nmasterconss] = conss[i];
         ++(*nmasterconss);
      }
   }

   if( *nmasterconss > 0 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, masterconss, *nmasterconss) );
   }
   else
   {
      SCIPfreeMemoryArray(scip, masterconss);
      *nmasterconss = 0;
      *masterconss = NULL;
   }

   *pricingisempty = (*nmasterconss == nconss);
   *masterisempty = (*nmasterconss == 0);

   return SCIP_OKAY;
}

static
bool nextBitmask(
   std::vector<bool>& bit_mask
   )
{ /*lint -esym(1793,std::_Bit_reference::operator=)*/
   std::size_t i = 0;
   for( ; (i < bit_mask.size()) && bit_mask[i]; ++i )
      bit_mask[i] = false;

   if( i < bit_mask.size() )
   {
      bit_mask[i] = true;
      return true;
   } else
      return false;
}

static
std::set<int> getSetFromBits(
   std::vector<bool> bits
)
{
   std::set<int> set;
   for( size_t i = 0; i < bits.size(); ++i )
   {
      if( bits[i] )
         (void) set.insert(int(i));
   }
   return set;
}

static
int nChooseK(
   int n,
   int k
   )
{
   int result = 1;
   for( int i = 1; i <= k; i++ )
   {
      result *= n - (k - i);
      result /= i;
   }
   return result;
}

/** looks for colors components in the constraints in detectordata */
static
SCIP_RETCODE findColorsComponents(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP***         decomps,            /**< vector of decomposition structure */
   int*                  ndecomps,           /**< number of decompositions */
   SCIP_RESULT*          result              /**< result pointer to indicate success oder failure */
   )
{

   SCIP_Bool masterisempty;
   SCIP_Bool pricingisempty;
   int* colors = NULL;
   int ncolors = 0;
   SCIP_CONS** conss = NULL;
   int nconss = 0;
   assert(scip != NULL);
   assert(decomps != NULL);
   assert(ndecomps != NULL);
   assert(result != NULL);

   *ndecomps = 0;
   *decomps = 0;
   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( SCIPallocMemoryArray(scip, &colors, SCIPgetNConss(scip)) ); /*lint !e666 */

   SCIP_CALL( assignConsColors(scip, SCIPgetConss(scip), SCIPgetNConss(scip), colors, &ncolors) );

   std::vector<bool> bit_mask((size_t)ncolors);

   SCIP_CALL( SCIPallocMemoryArray(scip, decomps, 1) ); /*lint !e506*/

   int nbits = 2;

   for( int subsetsize = 2; subsetsize <= nbits; ++subsetsize )
   {
      int size = nChooseK(ncolors, subsetsize);

      SCIP_CALL( SCIPreallocMemoryArray(scip, decomps, (size_t)*ndecomps + size) );

      do
      {
         if( std::count(bit_mask.begin(), bit_mask.end(), true) != subsetsize ) /*lint !e864*/
            continue;

         std::set<int> colorset = getSetFromBits(bit_mask);
#ifdef SCIP_DEBUG
         SCIPdebugMessage("Colors:");
         for( std::set<int>::iterator it = colorset.begin(); it != colorset.end(); ++it)
         {
            SCIPdebugPrintf(" %d", *it);
         }
         SCIPdebugPrintf("\n");
         printSubset(bit_mask);
#endif
         SCIP_CALL( createMasterconssArray(scip, &conss, &nconss, colors, colorset, &masterisempty, &pricingisempty) );

         SCIP_CALL( DECcreateDecompFromMasterconss(scip, &((*decomps)[*ndecomps]), conss, nconss) );
         ++(*ndecomps);

      } while( nextBitmask(bit_mask) );
   }
   if( *ndecomps > 0 )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, decomps, *ndecomps) );
      *result = SCIP_SUCCESS;
   }

   SCIPfreeMemoryArrayNull(scip, &colors);

   return SCIP_OKAY;
}



/** destructor of detector to free user data (called when GCG is exiting) */
static
DEC_DECL_FREEDETECTOR(detectorFreeColors)
{  /*lint --e{715}*/
   DEC_DETECTORDATA *detectordata;

   assert(scip != NULL);
   assert(detector != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   SCIPfreeMemory(scip, &detectordata);

   return SCIP_OKAY;
}

/** detector structure detection method, tries to detect a structure in the problem */
static
DEC_DECL_DETECTSTRUCTURE(detectorDetectColors)
{ /*lint -e715*/
   *result = SCIP_DIDNOTFIND;

   *ndecdecomps = 0;
   *decdecomps = NULL;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Detecting colored structure:");

   SCIP_CALL( findColorsComponents(scip, decdecomps, ndecdecomps, result) );

   if( *result == SCIP_SUCCESS )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " found %d decompositions.\n", *ndecdecomps);
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, " not found.\n");
   }

   if( *ndecdecomps == 0 )
   {
      SCIPfreeMemoryArrayNull(scip, decdecomps);
   }

   return SCIP_OKAY;
}

#define detectorPropagateSeeedColors NULL
#define detectorFinishSeeedColors NULL
#define detectorPostprocessSeeedColors NULL
#define detectorExitColors NULL
#define detectorInitColors NULL

#define setParamAggressiveColors NULL
#define setParamDefaultColors NULL
#define setParamFastColor NULL



/*
 * detector specific interface methods
 */

/** creates the handler for colors constraints and includes it in SCIP */
extern "C"
SCIP_RETCODE SCIPincludeDetectorColors(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   DEC_DETECTORDATA* detectordata;

   /* create colors constraint handler data */
   detectordata = NULL;

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );
   assert(detectordata != NULL);


   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_FREQCALLROUND, DEC_MAXCALLROUND, DEC_MINCALLROUND, DEC_FREQCALLROUNDORIGINAL, DEC_MAXCALLROUNDORIGINAL, DEC_MINCALLROUNDORIGINAL, DEC_PRIORITY, DEC_ENABLED, DEC_ENABLED_ORIGINAL, DEC_ENABLEDFINISHING, DEC_ENABLEDPOSTPROCESSING, DEC_SKIP, DEC_USEFULRECALL, DEC_LEGACYMODE,
         detectordata, detectorDetectColors, detectorFreeColors, detectorInitColors, detectorExitColors, detectorPropagateSeeedColors, NULL, NULL, detectorFinishSeeedColors, detectorPostprocessSeeedColors, setParamAggressiveColors, setParamDefaultColors, setParamFastColor) );


   /* add colors constraint handler parameters */

   return SCIP_OKAY;
}
