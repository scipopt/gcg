/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2019 Operations Research, RWTH Aachen University       */
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

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Michael Bastubbe
 * @author Hanna Franzen
 *
 * This constraint handler will run all registered structure detectors in a loop. They will find partial decompositions in a loop iteration until the decompositions are full
 * or the maximum number of detection rounds is reached.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include "cons_decomp.h"

#include "reader_gp.h"

#include <assert.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <utility>
#include <regex>
#include <vector>

#include <scip/clock.h>
#include <scip/def.h>
#include <scip/pub_cons.h>
#include <scip/pub_dialog.h>
#include <scip/pub_message.h>
#include <scip/pub_misc.h>
#include <scip/type_clock.h>
#include <scip/type_cons.h>
#include <scip/type_dialog.h>
#include <scip/type_message.h>
#include <scip/type_misc.h>
#include <scip/type_paramset.h>
#include <scip/type_result.h>
#include <scip/type_retcode.h>
#include <scip/type_scip.h>
#include <scip/type_set.h>
#include "class_consclassifier.h"
#include "class_seeed.h"
#include "class_seeedpool.h"
#include "class_varclassifier.h"
#include "class_miscvisualization.h"
#include "pub_decomp.h"
#include "type_decomp.h"
#include "wrapper_seeed.h"
#include "reader_tex.h"
#include "scip_misc.h"
#include "relax_gcg.h"



typedef gcg::Seeed* SeeedPtr;


/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"   /**< name of constraint handler */
#define CONSHDLR_DESC          "constraint handler for structure detection"   /**< description of constraint handler */
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                          *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define MAXNDECOMPS 5000                /**< indicates whether to create a decomposition with all constraints in the master if no other specified */

#define DEFAULT_CREATEBASICDECOMP FALSE /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
#define DEFAULT_DUALVALRANDOMMETHOD 1   /**< default value for method to dual initialization of dual values for strong decomposition: 1) naive, 2) expected equal, 3) expected overestimation */
#define DEFAULT_COEFFACTORORIGVSRANDOM 0.5 /**< default value for convex coefficient for orig dual val (1-this coef is factor for random dual value)  */

#define DEFAULT_BLOCKNUMBERCANDSMEDIANVARSPERCONS FALSE     /**< should for block number candidates calculation the medianvarspercons calculation be considered */

#define DEFAULT_ALLOWCLASSIFIERDUPLICATES FALSE       /** if false each new (conss- and vars-) classifier is checked for being a duplicate of an existing one, if so it is not added and NBOT statistically recognized*/
#define DEFAULT_MAXDETECTIONROUNDS 1    /**< maximal number of detection rounds */
#define DEFAULT_MAXNCLASSESLARGEPROBS 5   /** maximum number of classes allowed for large (nvars+nconss > 50000) MIPs for detectors, classifier with more classes are reduced to the maximum number of classes */
#define DEFAULT_MAXNCLASSES 9            /** maximum number of classes allowed for detectors, classifier with more classes are reduced to the maximum number of classes */
#define DEFAULT_MAXNCLASSESFORNBLOCKCANDIDATES 18                 /** maximum number of classes a classifier can have to be used for voting nblockcandidates */
#define DEFAULT_ENABLEORIGDETECTION FALSE                         /**< indicates whether to start detection for the original problem */
#define DEFAULT_CONSSADJCALCULATED                    TRUE        /**< indicates whether conss adjacency datastructures should be calculated, this might slow down initialization, but accelerating refinement methods*/
#define DEFAULT_ENABLEORIGCLASSIFICATION              FALSE       /**< indicates whether to start detection for the original problem */
#define DEFAULT_CONSSCLASSNNONZENABLED                TRUE        /**<  indicates whether constraint classifier for nonzero entries is enabled */
#define DEFAULT_CONSSCLASSNNONZENABLEDORIG            TRUE       /**<  indicates whether constraint classifier for nonzero entries is enabled for the original problem */

#define DEFAULT_CONSSCLASSSCIPCONSTYPEENABLED         TRUE        /**< indicates whether constraint classifier for scipconstype is enabled */
#define DEFAULT_CONSSCLASSSCIPCONSTYPEENABLEDORIG     TRUE       /**< indicates whether constraint classifier for scipconsstype is enabled for the original problem */

#define DEFAULT_AGGREGATIONLIMITNCONSSPERBLOCK        300        /**< if this limit on the number of constraints of a block is exceeded the aggregation information for this block is not calculated */
#define DEFAULT_AGGREGATIONLIMITNVARSPERBLOCK         300        /**< if this limit on the number of variables of a block is exceeded the aggregation information for this block is not calculated */


#define DEFAULT_CONSSCLASSMIPLIBCONSTYPEENABLED         TRUE      /**< indicates whether constraint classifier for miplib consstype is enabled */
#define DEFAULT_CONSSCLASSMIPLIBCONSTYPEENABLEDORIG     TRUE     /**< indicates whether constraint classifier for miplib consstype is enabled for the original problem */

#define DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLED     FALSE       /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled */
#define DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLEDORIG FALSE       /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem */

#define DEFAULT_CONSSCLASSLEVENSHTEINENABLED          FALSE       /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled */
#define DEFAULT_CONSSCLASSLEVENSHTEINENABLEDORIG      FALSE       /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem */

#define DEFAULT_VARCLASSSCIPVARTYPESENABLED           TRUE        /**< indicates whether variable classifier for scipvartypes is enabled */
#define DEFAULT_VARCLASSSCIPVARTYPESENABLEDORIG       TRUE       /**< indicates whether variable classifier for scipvartypes is enabled for the original problem */
#define DEFAULT_BENDERSONLYCONTSUBPR     FALSE      /**< indicates whether only decomposition with only continuous variables in the subproblems should be searched*/
#define DEFAULT_BENDERSONLYBINMASTER     FALSE      /**< indicates whether only decomposition with only binary variables in the master should be searched */

#define DEFAULT_VARCLASSOBJVALSENABLED                TRUE        /**< indicates whether variable classifier for objective function values is enabled */
#define DEFAULT_VARCLASSOBJVALSENABLEDORIG            TRUE       /**< indicates whether variable classifier for objective function values is enabled for the original problem */

#define DEFAULT_VARCLASSOBJVALSIGNSENABLED            TRUE        /**< indicates whether variable classifier for objective function value signs is enabled */
#define DEFAULT_VARCLASSOBJVALSIGNSENABLEDORIG        TRUE       /**< indicates whether variable classifier for objective function value signs is enabled for the original problem */

#define DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER    10000       /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis default */
#define AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER  80000      /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis aggressive */
#define FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER       2000        /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis fast */

#define DEFAULT_ONLYLEGACYMODE                        FALSE    /**< indicates whether detection should only consist of legacy mode detection */
#define DEFAULT_LEGACYMODE                            FALSE    /**< indicates whether detection should consist of legacy mode detection */
#define DEFAULT_STAIRLINKINGHEUR                      FALSE    /**< indicates whether heuristic to reassign linking vars to stairlinking in legacy mode should be activated */

#define DEFAULT_DETECTBENDERS                        FALSE    /**< indicates whether benders detection mode is enabled */


/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP*           useddecomp;                              /**< decomposition structures that was/will be used*/
   DEC_DECOMP**          decdecomps;                              /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;                               /**< array of structure detectors */
   int*                  priorities;                              /**< priorities of the detectors */
   int                   ndetectors;                              /**< number of detectors */
   SCIP_CLOCK*           detectorclock;                           /**< clock to measure detection time */
   SCIP_CLOCK*           completedetectionclock;                  /**< clock to measure detection time */
   SCIP_Bool             hasrun;                                  /**< flag to indicate whether we have already detected */
   int                   ndecomps;                                /**< number of decomposition structures  */
   int                   sizedecomps;                             /**< size of the decomp and complete seeeds array */
   int                   sizeincompleteseeeds;                    /**< size of the incomplete seeeds array */
   int                   maxndetectionrounds;                     /**< maximum number of detection loop rounds  */
   int                   strongdetectiondualvalrandommethod;      /**< method to dual initialization of dual values for strong decomposition: 1) naive, 2) expected equal, 3) expected overestimation */
   SCIP_Real             coeffactororigvsrandom;                  /**< convex coefficient for orig dual val (1-this coef is factor for random dual value)  */
   SCIP_Bool             blocknumbercandsmedianvarspercons;       /**< should for block number candidates calculation the medianvarspercons calculation be considered */
   int                   maxnclassesfornblockcandidates;          /**< maximum number of classes a classifier can have to be used for voting nblockcandidates */
   int                   maxnclassesperclassifier;                /**< maximum number of classes allowed for detectors, classifier with more classes are reduced to the maximum number of classes */
   int                   maxnclassesperclassifierforlargeprobs;   /**< maximum number of classes allowed for large (nvars+nconss > 50000) MIPs for detectors, classifier with more classes are reduced to the maximum number of classes */
   int                   weightinggpresolvedoriginaldecomps;      /**< weighing method for comparing presolved and original decompositions (see corresponding enum)   */
   int                   aggregationlimitnconssperblock;          /**< if this limit on the number of constraints of a block is exceeded the aggregation information for this block is not calculated */
   int                   aggregationlimitnvarsperblock;           /**< if this limit on the number of variables of a block is exceeded the aggregation information for this block is not calculated */
   SCIP_Bool             createbasicdecomp;                       /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
   SCIP_Bool             allowclassifierduplicates;               /**< indicates whether classifier duplicates are allowed (for statistical reasons) */
   SCIP_Bool             conssadjcalculated;                      /**< indicates whether conss adjacency datastructures should be calculated, this might slow down initialization, but accelerating refinement methods*/
   SCIP_Bool             enableorigdetection;                     /**< indicates whether to start detection for the original problem */
   SCIP_Bool             enableorigclassification;                /**< indicates whether to start constraint classification for the original problem */
   SCIP_Bool             conssclassnnonzenabled;                  /**< indicates whether constraint classifier for nonzero entries is enabled */
   SCIP_Bool             conssclassnnonzenabledorig;              /**< indicates whether constraint classifier for nonzero entries is enabled for the original problem */
   SCIP_Bool             conssclassnconstypeenabled;              /**< indicates whether constraint classifier for scipconstype is enabled */
   SCIP_Bool             conssclassnconstypeenabledorig;          /**< indicates whether constraint classifier for scipconstype is enabled for the original problem */
   SCIP_Bool             conssclassnmiplibconstypeenabled;        /**< indicates whether constraint classifier for miplib constype is enabled */
   SCIP_Bool             conssclassnmiplibconstypeenabledorig;    /**< indicates whether constraint classifier for miplib constype is enabled for the original problem */
   SCIP_Bool             consnamenonumbersenabled;                /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled */
   SCIP_Bool             consnamenonumbersenabledorig;            /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem */
   SCIP_Bool             conssclasslevenshteinabled;              /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled */
   SCIP_Bool             conssclasslevenshteinenabledorig;        /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem */
   SCIP_Bool             varclassvartypesenabled;                 /**< indicates whether variable classifier for scipvartypes is enabled */
   SCIP_Bool             varclassvartypesenabledorig;             /**< indicates whether variable classifier for scipvartypes is enabled for the original problem */
   SCIP_Bool             bendersonlycontsubpr;                    /**< indicates whether only decomposition with only continuous variables in the subproblems should be searched*/
   SCIP_Bool             bendersonlybinmaster;                    /**< indicates whether only decomposition with only binary variables in the master should be searched */
   SCIP_Bool             detectbenders;                           /**< indicates wethher or not benders detection mode is enabled */
   SCIP_Bool             varclassobjvalsenabled;                  /**< indicates whether variable classifier for objective function values is enabled */
   SCIP_Bool             varclassobjvalsenabledorig;              /**< indicates whether variable classifier for objective function values is enabled for the original problem */
   SCIP_Bool             varclassobjvalsignsenabled;              /**< indicates whether variable classifier for objective function value signs is enabled */
   SCIP_Bool             varclassobjvalsignsenabledorig;          /**< indicates whether variable classifier for objective function value signs is enabled for the original problem */
   SCIP_Bool             onlylegacymode;                          /**< indicates whether detection should only consist of legacy mode detection, this is sufficient to enable it */
   SCIP_Bool             legacymodeenabled;                       /**< indicates whether the legacy detection mode (detection before v3.0) additionally*/
   SCIP_Bool             stairlinkingheur;                        /**< indicates whether heuristic to reassign linking vars to stairlinking in legacy mode should be activated */
   int**                 candidatesNBlocks;                       /**< pointer to store candidates for number of blocks calculated by the seeedpool(s) */
   int*                  nCandidates;                             /**< number of candidates for number of blocks calculated by the seeedpool */
   int                   ncallscreatedecomp;                      /**< debugging method for counting the number of calls of created decompositions */
   gcg::Seeedpool*		 seeedpool;                               /**< seeedpool that manages the detection process for the presolved transformed problem */
   gcg::Seeedpool*       seeedpoolunpresolved;                    /**< seeedpool that manages the detection process of the unpresolved problem */
   SeeedPtr*             allrelevantfinishedseeeds;               /**< collection  of all relevant seeeds ( i.e. all seeeds w.r.t. copies ) */
   SeeedPtr*             incompleteseeeds;                        /**< collection of incomplete seeeds originating from incomplete decompositions given by the users */
   int                   nallrelevantseeeds;                      /**< number  of all relevant seeeds ( i.e. all seeeds w.r.t. copies ) */
   int                   nincompleteseeeds;                       /**< number  of incomplete seeeds originating from incomplete decompositions given by the users */
   SCIP_Bool             consnamesalreadyrepaired;                /**< stores whether or not    */
   SCIP_Bool             unpresolveduserseeedadded;               /**< stores whether or not an unpresolved user seeed was added */
   int                   seeedcounter;                            /**< counts the number of seeeds, used for seeed ids */
   std::vector<std::pair<SeeedPtr, SCIP_Real> >* candidates;      /**< vector containing the pairs of candidate list of decomps (to visualize, write, consider for family tree, consider for solving etc.) sorted according to  */
   int                   currscoretype;                           /**< indicates which score should be used for comparing (partial) decompositions
                                                                          0:max white,
                                                                          1: border area,
                                                                          2:classic,
                                                                          3:max foreseeing white,
                                                                          4: ppc-max-white,
                                                                          5:max foreseeing white with aggregation info,
                                                                          6: ppc-max-white with aggregation info,
                                                                          7: experimental benders score */
   std::vector<int>*       userblocknrcandidates;                 /**< vector to store block number candidates that were given by user */
   SCIP_Bool               nonfinalfreetransform;                 /**< help bool to notify a nonfinal free transform (needed if presolving is revoked, e.g. if unpresolved decomposition is used, and transformation is not successful) */
   SeeedPtr                seeedtowrite;                          /**< help pointer as interface for writing partial decompositions */
};


/** parameter how to modify scores when comparing decompositions for original and presolved problem
 * (which might differ in size) */
enum weightinggpresolvedoriginaldecomps{
   NO_MODIF = 0,           /**< no modification */
   FRACTION_OF_NNONZEROS,  /**< scores are weighted according to ratio of number nonzeros, the more the worse */
   FRACTION_OF_NROWS,      /**< scores are weighted according to ratio of number nonzeros, the more the worse */
   FAVOUR_PRESOLVED        /**< decompositions for presolved problems are always favoured over decompositions of original problem */
};

/*
 * Local methods
 */


/** local method to handle store a complete seeed in the unpresolved seeedpool
 *
 * @returns SCIP status */
static
SCIP_RETCODE  SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(
     SCIP* scip,     /**< SCIP data structure */
     SeeedPtr  seeed /**< pointer to seeed */
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool success;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert( seeed->isComplete() );
   assert( seeed->isFromUnpresolved() );

   conshdlrdata->seeedpoolunpresolved->addSeeedToFinished(seeed, &success);

   if( !success )
     SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Decomposition to add is already known to gcg!\n");

   return SCIP_OKAY;
   }

/** local method to handle store a complete seeed in the presolved seeedpool
 *
 * @returns SCIP status */
static
SCIP_RETCODE  SCIPconshdlrDecompAddCompleteSeeedForPresolved(
     SCIP* scip,     /**< SCIP data structure */
     SeeedPtr  seeed /**< pointer to seeed */
   )
{
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Bool success;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

      if( conshdlr == NULL )
      {
         SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
         return SCIP_ERROR;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

     assert( seeed->isComplete() );
     assert( !seeed->isFromUnpresolved() );

     conshdlrdata->seeedpool->addSeeedToFinished(seeed, &success);

     if( !success )
        SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Decomposition to add is already known to gcg!\n");

      return SCIP_OKAY;
}

/** local method to handle store a seeed as partial seeed to unpresolved seeedpool
 *
 * @returns SCIP status */
static
SCIP_RETCODE  SCIPconshdlrDecompAddPartialSeeedForUnpresolved(
     SCIP* scip,     /**< SCIP data structure */
     SeeedPtr  seeed /**< pointer to seeed */
   )
{
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;
      SCIP_Bool success;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

      if( conshdlr == NULL )
      {
         SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
         return SCIP_ERROR;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

     assert( !seeed->isComplete() );
     assert( seeed->isFromUnpresolved() );

     conshdlrdata->seeedpoolunpresolved->addSeeedToIncomplete(seeed, &success);

     if( !success )
        SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Decomposition to add is already known to gcg!\n");

     return SCIP_OKAY;
  }

/** local method to handle store a seeed as partial seeed to presolved seeedpool
 *
 * @returns SCIP status*/
static
SCIP_RETCODE  SCIPconshdlrDecompAddPartialSeeedForPresolved(
     SCIP* scip,     /**< SCIP data structure */
     SeeedPtr  seeed /**< pointer to partial seeed */
   )
{
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      SCIP_Bool success;
      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

      if( conshdlr == NULL )
      {
         SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
         return SCIP_ERROR;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

     assert( !seeed->isComplete() );
     assert( !seeed->isFromUnpresolved() );

     conshdlrdata->seeedpool->addSeeedToIncomplete(seeed, &success);

     if( !success )
        SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Decomposition to add is already known to gcg!\n");

     return SCIP_OKAY;
}


/** local method to handle store a seeed in the correct seeedpool
 *
 * @returns SCIP status */
static
SCIP_RETCODE  SCIPconshdlrDecompAddSeeed(
     SCIP* scip,        /**< SCIP data structure */
     SeeedPtr  seeed    /**< seeed pointer */
   )
{
      if( seeed->isComplete() )
      {
         if( seeed->isFromUnpresolved() )
            SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(scip, seeed);
         else
            SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, seeed);
      }
      else
      {
         if( seeed->isFromUnpresolved() )
            SCIPconshdlrDecompAddPartialSeeedForUnpresolved(scip, seeed);
         else
            SCIPconshdlrDecompAddPartialSeeedForPresolved(scip, seeed);
      }

      return SCIP_OKAY;
}

/** local method to find a seeed for a given id in presolved seeedpool or NULL if no seeed with such id is found
 *
 * @returns seeed pointer to seeed with given id */
static
SeeedPtr  SCIPconshdlrDecompGetSeeedFromPresolved(
     SCIP* scip,     /**< SCIP data structure */
     int  seeedid    /**< seeed id */
   )
{
      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

      if( conshdlr == NULL )
      {
         SCIPerrorMessage("Decomp constraint handler is not included, cannot find Seeed!\n");
         return NULL;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      if( conshdlrdata->seeedpool == NULL )
         return NULL;

      for( int i = 0; i < conshdlrdata->seeedpool->getNAncestorSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getAncestorSeeed( i ) != NULL && conshdlrdata->seeedpool->getAncestorSeeed( i )->getID() == seeedid )
            return conshdlrdata->seeedpool->getAncestorSeeed( i );
      }

      for( int i = 0; i < conshdlrdata->seeedpool->getNIncompleteSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getIncompleteSeeed( i )->getID() == seeedid )
            return conshdlrdata->seeedpool->getIncompleteSeeed( i );
      }

      for( int i = 0; i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getFinishedSeeed( i )->getID() == seeedid )
            return conshdlrdata->seeedpool->getFinishedSeeed( i );
      }

      return NULL;

}

/** local method to find a seeed for a given id in unpresolved seeedpool or NULL if no seeed with such id is found
 * @returns seeed pointer to seeed with given id */
static
SeeedPtr  SCIPconshdlrDecompGetSeeedFromUnpresolved(
     SCIP* scip,     /**< SCIP data structure */
     int  seeedid    /**< seeed id */
   ){

      SCIP_CONSHDLR* conshdlr;
      SCIP_CONSHDLRDATA* conshdlrdata;

      conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

      if( conshdlr == NULL )
      {
         SCIPerrorMessage("Decomp constraint handler is not included, cannot find Seeed!\n");
         return NULL;
      }

      conshdlrdata = SCIPconshdlrGetData(conshdlr);
      assert(conshdlrdata != NULL);

      if(conshdlrdata->seeedpoolunpresolved == NULL)
         return NULL;

      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNIncompleteSeeeds(); ++i)
      {
         if(  conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i )->getID() == seeedid )
            return conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i );
      }

      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNAncestorSeeeds(); ++i)
          {
             if( conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )!= NULL &&  conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )->getID() == seeedid )
                return conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i );
          }

      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i )->getID() == seeedid )
            return conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i );
      }

      return NULL;
}


/** local method to find a seeed for a given id or NULL if no seeed with such id is found
 * @returns seeed pointer of seeed with given id */
static
SeeedPtr SCIPconshdlrDecompGetSeeed(
   SCIP* scip,    /**< SCIP data structure */
   int  seeedid   /**< seeed id */
   )
{
   SeeedPtr seeed = NULL;

   seeed =  SCIPconshdlrDecompGetSeeedFromPresolved(scip, seeedid);

   if( seeed == NULL)
      return SCIPconshdlrDecompGetSeeedFromUnpresolved(scip, seeedid);
   else
      return seeed;
}


/** local method to get all seeeds
 * @returns vector of all seeeds*/
static
std::vector<SeeedPtr> getSeeeds(
   SCIP* scip    /**< SCIP data structure */
   )
{
   std::vector<SeeedPtr> allseeeds;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpool != NULL )
   {
      for( int i = 0; i < conshdlrdata->seeedpool->getNIncompleteSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getIncompleteSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpool->getIncompleteSeeed( i ));
      }

      for( int i = 0; i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getFinishedSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpool->getFinishedSeeed( i ));
      }
      for( int i = 0; i < conshdlrdata->seeedpool->getNAncestorSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getAncestorSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpool->getAncestorSeeed( i ));
      }
   }

   if( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNIncompleteSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i ));
      }

      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ));
      }
      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNAncestorSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i ));
      }
   }

   return allseeeds;
}


/** local method to get all seeeds but ancestor seeeds (all finished ones and the ones that are still incomplete)
 * @returns vector of all seeeds that might be chosen for optimization later */
static
std::vector<SeeedPtr> getLeafSeeeds(
   SCIP* scip    /**< SCIP data structure */
   )
{
   std::vector<SeeedPtr> allseeeds;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpool != NULL )
   {
      for( int i = 0; i < conshdlrdata->seeedpool->getNIncompleteSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getIncompleteSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpool->getIncompleteSeeed( i ));
      }

      for( int i = 0; i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpool->getFinishedSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpool->getFinishedSeeed( i ));
      }
   }

   if( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNIncompleteSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i ));
      }

      for( int i = 0; i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i)
      {
         if( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ) != NULL )
            allseeeds.push_back(conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ));
      }
   }

   return allseeeds;
}


/** gets all selected seeeds
 * @returns vector of all selected seeeds */
static
std::vector<SeeedPtr> getSelectedSeeeds(
   SCIP*          scip  /**< SCIP data structure */
   )
{
   std::vector<SeeedPtr> seeeds = getSeeeds(scip);

   std::vector<SeeedPtr> selectedseeeds;
   for(auto seeed : seeeds)
   {
      if(seeed->isSelected())
      {
         selectedseeeds.push_back(seeed);
      }
   }
   return selectedseeeds;
}


SCIP_RETCODE SCIPconshdlrdataDecompUnselectAll(
   SCIP*          scip
   )
{
   std::vector<SeeedPtr> seeeds = getSeeeds(scip);

   for(auto seeed : seeeds)
   {
      seeed->setSelected(false);
   }

   return SCIP_OKAY;
}


/** locally used macro to help with sorting, comparator */
struct sort_pred {
    bool operator()(const std::pair<SeeedPtr, SCIP_Real> &left, const std::pair<SeeedPtr, SCIP_Real> &right) {
        return left.second > right.second;
    }
};


#ifdef ADDONEBLOCKDECOMP
/**
 * create a 'decomposition' consisting of only one single block; used if no other decomposition was found
 *
 * @returns SCIP return code
 */
static
SCIP_RETCODE createOneBlockDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_HASHMAP* newconstoblock;
   DEC_DECOMP* newdecomp;
   SCIP_CONS** conss;
   int nconss;
   int i;

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   SCIP_CALL( SCIPhashmapCreate(&newconstoblock, SCIPblkmem(scip), nconss ) );

   /* assign each constraint to (the only) block 1 */
   for( i = 0; i < nconss; i++ )
   {
      assert(!SCIPhashmapExists(newconstoblock, conss[i]));
      SCIP_CALL( SCIPhashmapInsert(newconstoblock, conss[i], (void*) (size_t) 1) );
   }

   /* create the decomposition data structure and add it to SCIP */
   SCIP_CALL( DECdecompCreate(scip, &newdecomp) );
   assert(newdecomp != NULL);
   SCIP_CALL( DECfilloutDecompFromConstoblock(scip, newdecomp, newconstoblock, 1, FALSE) );

   SCIP_CALL( SCIPconshdlrDecompAddDecdecomp(scip, newdecomp) );

   SCIP_CALL( SCIPhashmapFree(&newconstoblock ) );

   SCIP_CALL( DECdecompFree(scip, &newdecomp ));

   return SCIP_OKAY;
}

#endif

/*
 * Callback methods of constraint handler
 */

/** initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->hasrun = FALSE;
   conshdlrdata->seeedpool = NULL;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      detector->dectime = 0.;
      if( detector->initDetector != NULL )
      {
         SCIPdebugMessage("Calling initDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->initDetector)(scip, detector) );
      }
   }

   return SCIP_OKAY;
}



/** deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);
   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->useddecomp != NULL )
      SCIP_CALL( DECdecompFree(scip, &conshdlrdata->useddecomp) );


   if( conshdlrdata->ndecomps > 0 && conshdlrdata->decdecomps != NULL )
   {
      for( int dec = 0; dec < conshdlrdata->ndecomps; ++dec )
      {

         DECdecompFree(scip, &conshdlrdata->decdecomps[conshdlrdata->ndecomps - dec - 1]);
      }

      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->decdecomps, conshdlrdata->ndecomps);
      conshdlrdata->ndecomps = 0;
      conshdlrdata->decdecomps = NULL;
   }


   conshdlrdata->hasrun = FALSE;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      SCIPfreeMemoryArrayNull(scip, &detector->decomps);
      if( detector->exitDetector != NULL )
      {
         SCIPdebugMessage("Calling exitDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->exitDetector)(scip, detector) );
      }
   }

   delete conshdlrdata->seeedpool;
   conshdlrdata->seeedpool = NULL;

   if( !conshdlrdata->nonfinalfreetransform )
   {
      if( conshdlrdata->seeedpoolunpresolved != NULL )
         delete conshdlrdata->seeedpoolunpresolved;
      conshdlrdata->seeedpoolunpresolved = NULL;
   }

   SCIPconshdlrdataDecompUnselectAll(scip);


   return SCIP_OKAY;
}

/** destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->completedetectionclock) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      if( detector->freeDetector != NULL )
      {
         SCIPdebugMessage("Calling freeDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->freeDetector)(scip, detector) );
      }
      SCIPfreeBlockMemory(scip, &detector);
   }

   /* @todo: This is also done in consExitDecomp() and therefore probably makes no sense here. */
   if ( conshdlrdata->useddecomp != NULL )
      SCIP_CALL( DECdecompFree(scip, &conshdlrdata->useddecomp) );

   if( conshdlrdata->candidates != NULL )
         delete conshdlrdata->candidates;

   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);

   delete conshdlrdata->userblocknrcandidates;

   SCIPfreeMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFORELAX(consEnforeDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsDecomp)
{  /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckDecomp)
{
   /*lint --e{715}*/
   *result = SCIP_FEASIBLE;
   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockDecomp)
{  /*lint --e{715}*/
   return SCIP_OKAY;
}

/*
 * constraint specific interface methods
 */

/*
 * @brief creates the constraint handler for decomp and includes it in SCIP
 * @returns scip return code
 */
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create decomp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);

   conshdlrdata->useddecomp = NULL;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;
   conshdlrdata->hasrun = FALSE;
   conshdlrdata->decdecomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->maxndetectionrounds = 0;
   conshdlrdata->maxnclassesperclassifier = 0;
   conshdlrdata->maxnclassesperclassifierforlargeprobs = 0;
   conshdlrdata->aggregationlimitnconssperblock = 0;
   conshdlrdata->aggregationlimitnvarsperblock = 0;
   conshdlrdata->enableorigdetection = FALSE;
   conshdlrdata->seeedpoolunpresolved = NULL;
   conshdlrdata->seeedpool = NULL;
   conshdlrdata->ncallscreatedecomp = 0;

   conshdlrdata->consnamesalreadyrepaired = FALSE;

   conshdlrdata->unpresolveduserseeedadded = FALSE;
   conshdlrdata->candidates = new std::vector<std::pair<SeeedPtr, SCIP_Real > >(0);
   conshdlrdata->sizedecomps = 10;
   conshdlrdata->seeedcounter = 0;
   conshdlrdata->currscoretype = scoretype::MAX_WHITE;
   conshdlrdata->nonfinalfreetransform = FALSE;
   conshdlrdata->userblocknrcandidates = new std::vector<int>(0);
   conshdlrdata->seeedtowrite = NULL;

   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->completedetectionclock) );
   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp, consLockDecomp,
         conshdlrdata) );
   assert(conshdlr != FALSE);

   SCIP_CALL( SCIPsetConshdlrEnforelax(scip, conshdlr, consEnforeDecomp) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDecomp) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitDecomp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitDecomp) );


   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/decomp/createbasicdecomp", "indicates whether to create a decomposition with all constraints in the master if no other specified", &conshdlrdata->createbasicdecomp, FALSE, DEFAULT_CREATEBASICDECOMP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/allowclassifierduplicates/enabled", "indicates whether classifier duplicates are allowed (for statistical reasons)", &conshdlrdata->allowclassifierduplicates, FALSE, DEFAULT_ALLOWCLASSIFIERDUPLICATES, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/conssadjcalculated", "conss adjecency datastructures should be calculated", &conshdlrdata->conssadjcalculated, FALSE, DEFAULT_CONSSADJCALCULATED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/origprob/enabled", "indicates whether to start detection for the original problem", &conshdlrdata->enableorigdetection, FALSE, DEFAULT_ENABLEORIGDETECTION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/origprob/classificationenabled", "indicates whether to classify constraints and variables for the original problem", &conshdlrdata->enableorigclassification, FALSE, DEFAULT_ENABLEORIGCLASSIFICATION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", "indicates whether constraint classifier for nonzero entries is enabled", &conshdlrdata->conssclassnnonzenabled, FALSE, DEFAULT_CONSSCLASSNNONZENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/nnonzeros/origenabled", "indicates whether constraint classifier for nonzero entries is enabled for the original problem", &conshdlrdata->conssclassnnonzenabledorig, FALSE, DEFAULT_CONSSCLASSNNONZENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/scipconstype/enabled", "indicates whether constraint classifier for scipconstype is enabled", &conshdlrdata->conssclassnconstypeenabled, FALSE, DEFAULT_CONSSCLASSSCIPCONSTYPEENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/scipconstype/origenabled", "indicates whether constraint classifier for scipconsstype is enabled for the original problem", &conshdlrdata->conssclassnconstypeenabledorig, FALSE, DEFAULT_CONSSCLASSSCIPCONSTYPEENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/miplibconstype/enabled", "indicates whether constraint classifier for miplib constypes is enabled", &conshdlrdata->conssclassnmiplibconstypeenabled, FALSE, DEFAULT_CONSSCLASSMIPLIBCONSTYPEENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/miplibconstype/origenabled", "indicates whether constraint classifier for miplib consstype is enabled for the original problem", &conshdlrdata->conssclassnmiplibconstypeenabledorig, FALSE, DEFAULT_CONSSCLASSMIPLIBCONSTYPEENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", "indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled", &conshdlrdata->consnamenonumbersenabled, FALSE, DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamenonumbers/origenabled", "indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem", &conshdlrdata->consnamenonumbersenabledorig, FALSE, DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", "indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled", &conshdlrdata->conssclasslevenshteinabled, FALSE, DEFAULT_CONSSCLASSLEVENSHTEINENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamelevenshtein/origenabled", "indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem", &conshdlrdata->conssclasslevenshteinenabledorig, FALSE, DEFAULT_CONSSCLASSLEVENSHTEINENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/scipvartype/enabled", "indicates whether variable classifier for scipvartypes is enabled", &conshdlrdata->varclassvartypesenabled, FALSE, DEFAULT_VARCLASSSCIPVARTYPESENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/scipvartype/origenabled", "indicates whether variable classifier for scipvartypes is enabled for the original problem", &conshdlrdata->varclassvartypesenabledorig, FALSE, DEFAULT_VARCLASSSCIPVARTYPESENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/objectivevalues/enabled", "indicates whether variable classifier for objective function values is enabled", &conshdlrdata->varclassobjvalsenabled,
      FALSE, DEFAULT_VARCLASSOBJVALSENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/blocknumbercandidates/medianvarspercons", "should for block number candidates calcluation the medianvarspercons calculation be considered", &conshdlrdata->blocknumbercandsmedianvarspercons,
      FALSE, DEFAULT_BLOCKNUMBERCANDSMEDIANVARSPERCONS, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/onlycontsubpr", "indicates whether only decomposition with only continiuous variables in the subproblems should be searched", &conshdlrdata->bendersonlycontsubpr, FALSE, DEFAULT_BENDERSONLYCONTSUBPR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/onlybinmaster", "indicates whether only decomposition with only binary variables in the master should be searched", &conshdlrdata->bendersonlybinmaster, FALSE, DEFAULT_BENDERSONLYBINMASTER, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/objectivevalues/origenabled", "indicates whether variable classifier for objective function values is enabled for the original problem", &conshdlrdata->varclassobjvalsenabledorig, FALSE, DEFAULT_VARCLASSOBJVALSENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/objectivevaluesigns/enabled", "indicates whether variable classifier for objective function value signs is enabled", &conshdlrdata->varclassobjvalsignsenabled, FALSE, DEFAULT_VARCLASSOBJVALSIGNSENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/objectivevaluesigns/origenabled", "indicates whether variable classifier for objective function value signs is enabled for the original problem", &conshdlrdata->varclassobjvalsignsenabledorig, FALSE, DEFAULT_VARCLASSOBJVALSIGNSENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/legacymode/onlylegacymode", "indicates whether detection should only consist of legacy mode detection", &conshdlrdata->onlylegacymode, FALSE, DEFAULT_ONLYLEGACYMODE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/legacymode/enabled", "indicates whether detection consist of legacy mode detection", &conshdlrdata->legacymodeenabled, FALSE, DEFAULT_LEGACYMODE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/legacymode/stairlinkingheur", "indicates whether heuristic to reassign linking vars to stairlinking in legacy mode should be activated", &conshdlrdata->stairlinkingheur, FALSE, DEFAULT_STAIRLINKINGHEUR, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/enabled", "indicates whether benders detection is enabled", &conshdlrdata->detectbenders, FALSE, DEFAULT_DETECTBENDERS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/strong_detection/dualvalrandommethod",
      "Method for random dual values use for strong decomposition: 1) naive, 2) expected equality exponential distributed, 3) expected overestimation exponential distributed ", &conshdlrdata->strongdetectiondualvalrandommethod, FALSE,
      DEFAULT_DUALVALRANDOMMETHOD, 1, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddRealParam(scip, "detection/strong_detection/coeffactororigvsrandom",
      " convex coefficient for orig dual val (1-this coef is factor for random dual value) ", &conshdlrdata->coeffactororigvsrandom, FALSE,
      DEFAULT_COEFFACTORORIGVSRANDOM, 0., 1., NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxrounds",
      "Maximum number of detection loop rounds", &conshdlrdata->maxndetectionrounds, FALSE,
      DEFAULT_MAXDETECTIONROUNDS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxnclassesfornblockcandidates",
      "Maximum number of classes a classifier can have to be used for voting nblockcandidates", &conshdlrdata->maxnclassesfornblockcandidates, FALSE,
      DEFAULT_MAXNCLASSESFORNBLOCKCANDIDATES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxnclassesperclassifier",
      "Maximum number of classes per classifier", &conshdlrdata->maxnclassesperclassifier, FALSE,
      DEFAULT_MAXNCLASSES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/aggregation/limitnconssperblock",
      "if this limit on the number of constraints of a block is exceeded the aggregation information for this block is not calculated ", &conshdlrdata->aggregationlimitnconssperblock, FALSE,
      DEFAULT_AGGREGATIONLIMITNCONSSPERBLOCK, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/aggregation/limitnvarsperblock",
      "if this limit on the number of variables of a block is exceeded the aggregation information for this block is not calculated ", &conshdlrdata->aggregationlimitnvarsperblock, FALSE,
      DEFAULT_AGGREGATIONLIMITNVARSPERBLOCK, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxnclassesperclassifierforlargeprobs",
      "Maximum number of classes per classifier for large problems (nconss + nvars >= 50000)", &conshdlrdata->maxnclassesperclassifierforlargeprobs, FALSE,
      DEFAULT_MAXNCLASSESLARGEPROBS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/origprob/weightinggpresolvedoriginaldecomps",
      "Weighting method when comparing decompositions for presolved and unpresolved problem", &conshdlrdata->weightinggpresolvedoriginaldecomps, TRUE,
      NO_MODIF, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/scoretype",
         "indicates which score should be used for comparing (partial) decompositions (0: max white, 1: border area, 2:classic, 3:max foreseeing white, 4: ppc-max-white, 5:max foreseeing white with aggregation info, 6: ppc-max-white with aggregation info, 7: experimental benders score, 8:strong decomposition score): ", &conshdlrdata->currscoretype, FALSE,
         scoretype::SETPART_FWHITE, 0, 8, NULL, NULL) );

   assert(conshdlrdata->candidates != NULL);

   return SCIP_OKAY;
}


/**
 * finds a non duplicate constraint name of the form c_{a} with minimal natural number {a}
 * @return non duplicate constraint name of the form c_{a} with minimal natural number {a}
 */
static
int findGenericConsname(
   SCIP*              scip,                  /**< SCIP data structure */
   int                startcount,            /**< natural number, lowest candidate number to test */
   char*              consname,              /**< char pointer to store the new non-duplicate name */
   int                namelength             /**< max length of the name */
   )
{
   int candidatenumber;

   candidatenumber = startcount;

   /* terminates since there are only finitely many constraints and i (for c_i) increases every iteration */
   while( TRUE )
   {
      char candidatename[SCIP_MAXSTRLEN] = "c_";
      char number[20];
      sprintf(number, "%d", candidatenumber );
      strcat(candidatename, number );

      if ( SCIPfindCons( scip, candidatename ) == NULL )
      {
         strncpy(consname, candidatename, namelength - 1);
         return candidatenumber;
      }
      else
         ++candidatenumber;
   }
   return -1;
}


/*
 * method to eliminate duplicate constraint names and name unnamed constraints
 * @return SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(
   SCIP*                 scip                   /* SCIP data structure */
   )
{
   long int startcount;
   SCIP_CONS** conss;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   startcount = 1;

   if( conshdlrdata->consnamesalreadyrepaired )
      return SCIP_OKAY;

   unordered_map<std::string, bool> consnamemap;

   SCIPdebugMessage("start repair conss \n ");

   conss = SCIPgetConss(scip);

   for( int i = 0; i < (int) SCIPgetNConss(scip); ++i )
   {
      SCIP_CONS* cons = conss[i];

      SCIPdebugMessage( "cons name: %s\n ", SCIPconsGetName(cons));

      if( SCIPconsGetName(cons) == NULL || strcmp(SCIPconsGetName(cons), "") == 0 || consnamemap[SCIPconsGetName(cons)] )
      {
         if( SCIPgetStage(scip) <= SCIP_STAGE_PROBLEM )
         {
            char newconsname[SCIP_MAXSTRLEN];
            startcount = findGenericConsname(scip, startcount, newconsname, SCIP_MAXSTRLEN ) + 1;
            SCIPdebugMessage( "Change consname to %s\n", newconsname );
            SCIPchgConsName(scip, cons, newconsname );
            consnamemap[newconsname] = true;
         }
         else
         {
            if ( SCIPconsGetName(cons) == NULL )
               SCIPwarningMessage(scip, "Name of constraint is NULL \n");
            else if ( strcmp(SCIPconsGetName(cons), "") == 0 )
               SCIPwarningMessage(scip, "Name of constraint is not set \n");
            else
               SCIPwarningMessage(scip, "Constraint name duplicate: %s \n", SCIPconsGetName(cons) );
         }
      }
      else
      {
         consnamemap[SCIPconsGetName(cons)] = true;
      }

      SCIPdebugMessage( " number of elements: %d \n " , (int) consnamemap.size() );
   }

   conshdlrdata->consnamesalreadyrepaired = TRUE;

   return SCIP_OKAY;
}


/*
 * @brief sets (and adds) the decomposition structure
 * @note this method should only be called if there is no seeed for this decomposition
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /* SCIP data structure */
   DEC_DECOMP*           decdecomp           /* DEC_DECOMP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SeeedPtr seeed;

   if( conshdlrdata->seeedpool == NULL )
      SCIPconshdlrDecompCreateSeeedpool(scip);

   DECdecompSetPresolved(decdecomp, TRUE);

   SCIP_CALL( conshdlrdata->seeedpool->createSeeedFromDecomp(decdecomp, &seeed) );

   SCIP_CALL( SCIPconshdlrDecompAddSeeed(scip, seeed) );

   DECdecompFree(scip, &decdecomp);

   return SCIP_OKAY;
}


/** Checks whether benders detection is enabled
 *
 * @returns true if benders is enabled, false otherwise */
static
SCIP_Bool SCIPconshdlrDecompDetectBenders(
   SCIP*                   scip  /**< SCIP data structure */
   )
{
   SCIP_Bool benders;

   SCIPgetBoolParam(scip, "detection/benders/enabled", &benders);

   return benders;
}

/* Checks whether the currently best candidate is from the unpresolved seeedpool
 *
 * @returns true if best candidate is unpresolved, false otherwise */
SCIP_Bool SCIPconshdlrDecompIsBestCandidateUnpresolved(
   SCIP*                   scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->candidates->size() == 0 )
      return FALSE;

   return conshdlrdata->candidates->at(0).first->isFromUnpresolved();
}


/* \brief returns an array containing all decompositions
 *
 *  Updates the decdecomp decomposition structure by converting all finished seeeds into decompositions and replacing the
 *  old list in the conshdlr.
 *
 *  @returns decomposition array
 *   */
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int decompcounter;

   decompcounter = 0;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   for( int i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      DECdecompFree(scip, &conshdlrdata->decdecomps[conshdlrdata->ndecomps - i - 1]);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->decdecomps, conshdlrdata->ndecomps);

   SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &conshdlrdata->decdecomps, SCIPconshdlrDecompGetNDecdecomps(scip) ) );

   conshdlrdata->ndecomps = SCIPconshdlrDecompGetNDecdecomps(scip) ;

   if( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      for( int i = 0; i <  conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
      {
         gcg::Seeed* seeed = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed(i);
         conshdlrdata->seeedpoolunpresolved->createDecompFromSeeed(seeed, &conshdlrdata->decdecomps[decompcounter] );
         ++decompcounter;
      }
   }
   if( conshdlrdata->seeedpool != NULL )
   {
      for( int i = 0; i <  conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
      {
         gcg::Seeed* seeed = conshdlrdata->seeedpool->getFinishedSeeed(i);
         conshdlrdata->seeedpool->createDecompFromSeeed(seeed, &conshdlrdata->decdecomps[decompcounter] );
         ++decompcounter;
      }
   }
   return conshdlrdata->decdecomps;
}

/* gets the number of decompositions (= amount of finished seeeds)
 *
 * @returns number of decompositions */
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int ndecomps;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ndecomps = 0;

   if( conshdlrdata->seeedpoolunpresolved != NULL )
      ndecomps += conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds();

   if( conshdlrdata->seeedpool != NULL )
      ndecomps += conshdlrdata->seeedpool->getNFinishedSeeeds();

   return ndecomps;
}

/*
 * @brief returns the data of the provided detector
 * @returns data of the provided detector
 */
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /* detector data structure */
   )
{
   assert(detector != NULL);
   return detector->decdata;
}


/*
 * Gets the number of constraints that were active while detecting the decomposition originating from the seeed with the
 * given id, this method is used to decide if the problem has changed since detection, if so the aggregation information
 * needs to be recalculated

 * @returns number of constraints that were active while detecting the decomposition
 */
int SCIPconshdlrDecompGetNFormerDetectionConssForID(
   SCIP*                 scip,               /* SCIP data structure */
   int                   id                  /* id of the seeed */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   gcg::Seeed* seeed;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeed = conshdlrdata->seeedpool->findFinishedSeeedByID(id);
   currseeedpool = conshdlrdata->seeedpool;

   if ( seeed == NULL )
   {
      seeed = conshdlrdata->seeedpoolunpresolved->findFinishedSeeedByID(id);
      currseeedpool = conshdlrdata->seeedpoolunpresolved;
   }

   /* seeed is not found hence we should not trust the isomorph information from detection */
   if (seeed == NULL)
      return -1;

   return currseeedpool->getNConss();
}


/*
 * @brief creates the seeedpool for the presolved problem
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpool(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpool == NULL )
      conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));

   return SCIP_OKAY;
}

/*
 * @brief creates the seeedpool for the unpresolved problem
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpoolUnpresolved(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpoolunpresolved == NULL )
      conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPconshdlrDecompGetSeeedpoolUnpresolved(
   SCIP*                 scip,
   SEEED_WRAPPER*        sw
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sw->seeedpool = conshdlrdata->seeedpoolunpresolved;

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPconshdlrDecompGetSeeedpool(
   SCIP*                 scip,
   SEEED_WRAPPER*        sw
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   sw->seeedpool = conshdlrdata->seeedpool;

   return SCIP_OKAY;

}


/*
 * @brief counts up the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently increased
 */
int SCIPconshdlrDecompIncreaseAndGetNCallsCreateDecomp(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ++conshdlrdata->ncallscreatedecomp;

   return conshdlrdata->ncallscreatedecomp;
}

/*
 * @brief decreases the counter for created decompositions and returns it
 * @returns number of created decompositions that was recently decreased
 */
int SCIPconshdlrDecompDecreaseAndGetNCallsCreateDecomp(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   --conshdlrdata->ncallscreatedecomp;

   return conshdlrdata->ncallscreatedecomp;
}


/*
 * @brief returns the name of the provided detector
 * @returns name of the given detector
 */
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /* detector data structure */
   )
{
   assert(detector != NULL);
   return detector->name;
}

/*
 * @brief searches for the detector with the given name and returns it or NULL if detector is not found
 * @returns detector pointer or NULL if detector with given name is not found
 */
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /* SCIP data structure */
   const char*           name                /* name of the detector */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( strcmp(detector->name, name) == 0 )
      {
         return detector;
      }
   }

   return NULL;
}

/*
 * @brief includes one detector
 * @param scip scip data structure
 * @param name name of the detector
 * @param decchar char that is used in detector chain history for this detector
 * @param description describing main idea of this detector
 * @param freqCallRound frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0
 * @param maxCallRound last detection round the detector gets called
 * @param minCallRound first round the detector gets called (offset in detection loop)
 * @param freqCallRoundOriginal frequency the detector gets called in detection loop while detecting of the original problem
 * @param maxCallRoundOriginal last round the detector gets called while detecting of the original problem
 * @param minCallRoundOriginal first round the detector gets called (offset in detection loop) while detecting of the original problem
 * @param priority  priority of the detector
 * @param enabled whether the detector should be enabled by default
 * @param enabledOriginal whether the detector should be enabled by default for detecting the original problem
 * @param enabledFinishing whether the finishing should be enabled
 * @param enabledPostprocessing whether the postprocessing should be enabled
 * @param skip whether the detector should be skipped if others found structure
 * @param usefulRecall is it useful to call this detector on a descendant of the propagated seeed
 * @param legacymode whether (old) DETECTSTRUCTURE method should also be used for detection
 * @param detectordata the associated detector data (or NULL)
 * @param DEC_DECL_DETECTSTRUCTURE((*detectStructure))   the method that will detect the structure (may be NULL), only used in legecy detection mode
 * @param DEC_DECL_FREEDETECTOR((*freeDetector)) destructor of detector (or NULL)
 * @param DEC_DECL_INITDETECTOR((*initDetector)) initialization method of detector (or NULL)
 * @param DEC_DECL_EXITDETECTOR((*exitDetector)) deinitialization method of detector (or NULL)
 * @param DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)) method to refine a partial decomposition inside detection loop (or NULL)
 * @param DEC_DECL_FINISHSEEED((*finishSeeedDetector)) method to complete a partial decomposition when called in detection loop (or NULL)
 * @param DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)) method to postprocess a complete decomposition, called after detection loop (or NULL)
 * @param DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)) method that is called if the detection emphasis setting aggressive is chosen
 * @param DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector))  method that is called if the detection emphasis setting default is chosen
 * @param DEC_DECL_SETPARAMFAST((*setParamFastDetector))  method that is called if the detection emphasis setting fast is chosen
 * @returns scip return code
 */
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,                   /* SCIP data structure */
   const char*           name,                   /* name of the detector */
   const char            decchar,                /* display character of the detector */
   const char*           description,            /* description of the detector */
   int                   freqCallRound,          /* frequency the detector gets called in detection loop ,ie it is called
                                                      in round r if and only if minCallRound <= r <= maxCallRound AND
                                                      (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound,           /* last round the detector gets called                              */
   int                   minCallRound,           /* first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal,  /* frequency the detector gets called in detection loop while detecting
                                                      of the original problem */
   int                   maxCallRoundOriginal,   /* last round the detector gets called while detecting of the original
                                                      problem */
   int                   minCallRoundOriginal,   /* first round the detector gets called (offset in detection loop) while
                                                      detecting of the original problem */
   int                   priority,               /* priority of the detector                                           */
   SCIP_Bool             enabled,                /* whether the detector should be enabled by default                  */
   SCIP_Bool             enabledOriginal,        /* whether the detector should be enabled by default for detecting the
                                                      original problem */
   SCIP_Bool             enabledFinishing,       /* whether the finishing should be enabled */
   SCIP_Bool             enabledPostprocessing,  /* whether the postprocessing should be enabled */
   SCIP_Bool             skip,                   /* whether the detector should be skipped if others found structure   */
   SCIP_Bool             usefulRecall,           /* is it useful to call this detector on a descendant of the propagated
                                                      seeed */
   SCIP_Bool             legacymode,             /* whether (old) DETECTSTRUCTURE method should also be used for
                                                      detection */
   DEC_DETECTORDATA*     detectordata,           /* the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)), /* the method that will detect the structure (must not be NULL)*/
   DEC_DECL_FREEDETECTOR((*freeDetector)),       /* destructor of detector (or NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),       /* initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector)),       /* deinitialization method of detector (or NULL) */
   DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)),   /* propagation method of detector (or NULL) */
   DEC_DECL_FINISHSEEED((*finishSeeedDetector)),               /* finish method of detector (or NULL) */
   DEC_DECL_POSTPROCESSSEEED((*postprocessSeeedDetector)),     /* postprocess method of detector (or NULL) */
   DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)), /* set method for aggressive parameters of detector
                                                                    (or NULL) */
   DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),       /* set method for default parameters of detector
                                                                    (or NULL) */
   DEC_DECL_SETPARAMFAST((*setParamFastDetector))              /* set method for fast parameters of detector (or NULL) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPallocBlockMemory(scip, &detector) );
   assert(detector != NULL);

   SCIPdebugMessage("Adding detector %i: %s\n", conshdlrdata->ndetectors+1, name);

#ifndef NDEBUG
   assert(DECfindDetector(scip, name) == NULL);
#endif

   /* set meta data of detector */
   detector->decdata = detectordata;
   detector->name = name;
   detector->description = description;
   detector->decchar = decchar;

   /* set memory handling and detection functions */
   detector->freeDetector = freeDetector;
   detector->initDetector = initDetector;
   detector->exitDetector = exitDetector;
   detector->detectStructure = detectStructure;

   /* set functions for editing seeeds */
   detector->propagateSeeed = propagateSeeedDetector;
   detector->finishSeeed = finishSeeedDetector;
   detector->postprocessSeeed = postprocessSeeedDetector;

   /* initialize parameters */
   detector->setParamAggressive =  setParamAggressiveDetector;
   detector->setParamDefault =  setParamDefaultDetector;
   detector->setParamFast =  setParamFastDetector;
   detector->freqCallRound = freqCallRound;
   detector->maxCallRound = maxCallRound;
   detector->minCallRound = minCallRound;
   detector->freqCallRoundOriginal = freqCallRoundOriginal;
   detector->maxCallRoundOriginal = maxCallRoundOriginal;
   detector->minCallRoundOriginal= minCallRoundOriginal;
   detector->priority = priority;
   detector->enabled = enabled;
   detector->enabledOrig = enabledOriginal;
   detector->enabledFinishing = enabledFinishing;
   detector->enabledPostprocessing = enabledPostprocessing;
   detector->skip = skip;
   detector->usefulRecall = usefulRecall;
   detector->legacymode = legacymode;
   detector->overruleemphasis = FALSE;
   detector->ndecomps = 0;
   detector->decomps = NULL;
   detector->dectime = 0.;

   /* add and initialize all parameters accessable from menu */
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabled), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origenabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled for detecting in the original problem", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabledOrig), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled for finishing of incomplete decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabledFinishing), FALSE, enabledFinishing, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled for postprocessing of finished decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabledPostprocessing), FALSE, enabledPostprocessing, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/skip", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be skipped if others found decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->skip), FALSE, skip, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/usefullrecall", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be called on descendants of the current seeed", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->usefulRecall), FALSE, usefulRecall, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/legacymode", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether (old) DETECTSTRUCTURE method of detector <%s> should also be used for detection", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->legacymode), FALSE, legacymode, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/overruleemphasis", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether emphasis settings for detector <%s> should be overruled by normal settings", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->overruleemphasis), FALSE, FALSE, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/freqcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->freqCallRound), FALSE, freqCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/maxcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "maximum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->maxCallRound), FALSE, maxCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/mincallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "minimum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->minCallRound), FALSE, minCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origfreqcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->freqCallRoundOriginal), FALSE, freqCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origmaxcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "maximum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->maxCallRoundOriginal), FALSE, maxCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/origmincallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "minimum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->minCallRoundOriginal), FALSE, minCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/priority", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "priority of detector <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->priority), FALSE, priority, INT_MIN, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, (size_t)conshdlrdata->ndetectors+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->priorities,(size_t) conshdlrdata->ndetectors+1) );

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   return SCIP_OKAY;
}

/*
 * @brief returns the remaining time of scip that the decomposition may use
 * @returns remaining  time that the decompositon may use
 */
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}

/*
 * checks if two pricing problems are identical based on information from detection
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompArePricingprobsIdenticalForSeeedid(
   SCIP*                scip,       /* SCIP data structure */
   int                  seeedid,    /* id of seeed */
   int                  probnr1,    /* block id of first problem */
   int                  probnr2,    /* block id of second problem */
   SCIP_Bool*           identical   /* output Bool, true if identical */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeed* seeed;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeed = conshdlrdata->seeedpool->findFinishedSeeedByID(seeedid);

   if ( seeed == NULL )
   {
      seeed = conshdlrdata->seeedpoolunpresolved->findFinishedSeeedByID(seeedid);
   }


   if( seeed->getNReps() == 0 )
   {
      SCIPdebugMessage("calc aggregation information for seeed!\n");
      seeed->calcAggregationInformation();
   }

   assert(seeed != NULL);

   if( seeed->getRepForBlock(probnr1) == seeed->getRepForBlock(probnr2) )
      *identical = TRUE;
   else
      *identical = FALSE;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " block %d and block %d are represented by %d and %d hence they are identical=%d.\n", probnr1, probnr2, seeed->getRepForBlock(probnr1), seeed->getRepForBlock(probnr2), *identical );

   return SCIP_OKAY;
}

/*
 * @brief for two identical pricing problems a corresponding varmap is created
 * @param scip scip data structure
 * @param hashorig2pricingvar  mapping from orig to pricingvar
 * @param seeedid id of the partial decompostion for which the pricing problems are checked for identity
 * @param probnr1 index of first block
 * @param probnr2 index of second block
 * @param scip1 subscip of first block
 * @param scip2 subscip of second block
 * @param varmap mapping from orig to pricingvar
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateVarmapForSeeedId(
   SCIP*                scip,                /* SCIP data structure */
   SCIP_HASHMAP**       hashorig2pricingvar, /* mapping from orig to pricingvar  */
   int                  seeedid,             /* id of seeed for which to create the varmap */
   int                  probnr1,             /* block id of first problem */
   int                  probnr2,             /* block id of second problem */
   SCIP*                scip1,               /* SCIP data structure for first problem */
   SCIP*                scip2,               /* SCIP data structure for second problem */
   SCIP_HASHMAP*        varmap               /* output varmap */
   )
{
   gcg::Seeed* seeed;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;

   int blockid1;
   int blockid2;
   int representative;
   int repid1;
   int repid2;
   int nblocksforrep;
   std::vector<int> pidtopid;

   repid1 = -1;
   repid2 = -1;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeed = conshdlrdata->seeedpool->findFinishedSeeedByID(seeedid);
   currseeedpool = conshdlrdata->seeedpool;

   if ( seeed == NULL )
   {
      seeed = conshdlrdata->seeedpoolunpresolved->findFinishedSeeedByID(seeedid);
      currseeedpool = conshdlrdata->seeedpoolunpresolved;
   }

   assert(seeed != NULL);

   if( probnr1 > probnr2 )
   {
      blockid1 = probnr2;
      blockid2 = probnr1;
   }
   else
   {
      blockid1 = probnr1;
      blockid2 = probnr2;
   }

   representative = seeed->getRepForBlock(blockid1);
   assert( representative == seeed->getRepForBlock(blockid2) );
   nblocksforrep = (int) seeed->getBlocksForRep(representative).size();

   /* find index in representatives */
   for( int i = 0; i < nblocksforrep; ++i )
   {
      if( seeed->getBlocksForRep(representative)[i] == blockid1 )
         repid1 = i;
      if( seeed->getBlocksForRep(representative)[i] == blockid2 )
      {
         repid2 = i;
         break;
      }
   }

   /* blockid1 should be the representative */
   if( repid1 != 0 )
   {
      SCIPhashmapFree(&varmap);
      varmap = NULL;
      SCIPwarningMessage(scip, NULL, "blockid1 should be the representative (hence has id=0 in reptoblocksarray but in fact has %d) \n", repid1);
      return SCIP_OKAY;
   }

   pidtopid = seeed->getRepVarmap(representative, repid2);

   for( int v = 0; v < SCIPgetNVars(scip2); ++v )
   {
      SCIP_VAR* var1;
      SCIP_VAR* var2;
      SCIP_VAR* var1orig;
      SCIP_VAR* var2orig;
      int var1origid;
      int var2origid;
      int var1originblockid;
      int var2originblockid;

      var2 = SCIPgetVars(scip2)[v];
      assert(var2 != NULL);
      var2orig = GCGpricingVarGetOriginalVar(var2);
      assert(var2orig!=NULL);
      var2origid = currseeedpool->getIndexForVar(var2orig) ;
      assert(var2origid>=0);
      var2originblockid = seeed->getVarProbindexForBlock(var2origid, blockid2) ;
      assert(var2originblockid >= 0);
      var1originblockid = pidtopid[var2originblockid];
      assert(var1originblockid>=0);
      var1origid = seeed->getVarsForBlock(blockid1)[var1originblockid];
      assert(var1origid>=0);
      var1orig = currseeedpool->getVarForIndex(var1origid) ;
      assert(var1orig != NULL);
      var1 = (SCIP_VAR*) SCIPhashmapGetImage(hashorig2pricingvar[blockid1], (void*) var1orig ) ;
      assert(var1 != NULL);

      SCIPhashmapInsert(varmap, (void*) var2, (void*) var1);
   }

   return SCIP_OKAY;
}


/*
 * @brief returns whether or not an unpresolved (untransformed) decompositions exists in the data structures
 * @returns SCIP return code
 */
SCIP_Bool SCIPconshdlrDecompUnpresolvedSeeedExists(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if ( conshdlrdata->seeedpoolunpresolved == NULL )
      return FALSE;


   return ( conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() > 0 );

}


/*
 * @brief add block number user candidate (user candidates are prioritized over found ones)
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompAddBlockNumberCandidate(
   SCIP*                 scip,                /* SCIP data structure */
   int                   blockNumberCandidate /* new block number candidate */
   ){

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->userblocknrcandidates->push_back(blockNumberCandidate);

   if( conshdlrdata->seeedpool != NULL ){
      conshdlrdata->seeedpool->addUserCandidatesNBlocks(blockNumberCandidate);
   }

   if( conshdlrdata->seeedpoolunpresolved != NULL )
         conshdlrdata->seeedpoolunpresolved->addUserCandidatesNBlocks(blockNumberCandidate);

   return SCIP_OKAY;
}

/*
 * @brief returns the number of block candidates given by the user
 * @returns number of block candidates given by the user
 */
int SCIPconshdlrDecompGetNBlockNumberCandidates(
  SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   return (int) conshdlrdata->userblocknrcandidates->size();
}


/*
 * @brief returns block number user candidate with given index
 * @param scip SCIP data structure
 * @param index index of block number user candidate that should be returned
 * @returns block number user candidate with given index
 */
int SCIPconshdlrDecompGetBlockNumberCandidate(
   SCIP*                 scip,                /* SCIP data structure */
   int                   index                /* index in conshdlrdata->userblocknrcandidates vector */
    ){
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   return conshdlrdata->userblocknrcandidates->at(index);
}

/*
 * @brief returns the total detection time, including classification, score computation, etc.
 * @param scip SCIP data structure
 * @returns total detection time
 */
SCIP_Real SCIPconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip   /* SCIP data structure */
    )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Real totaltime;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   totaltime = SCIPgetClockTime(scip,  conshdlrdata->completedetectionclock );

   /*
   if( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      totaltime += conshdlrdata->seeedpoolunpresolved->classificationtime;
      totaltime += conshdlrdata->seeedpoolunpresolved->nblockscandidatescalctime;
      totaltime += conshdlrdata->seeedpoolunpresolved->postprocessingtime;
      totaltime += conshdlrdata->seeedpoolunpresolved->scorecalculatingtime;
      totaltime += conshdlrdata->seeedpoolunpresolved->translatingtime;
   }
   if( conshdlrdata->seeedpool != NULL )
   {
      totaltime += conshdlrdata->seeedpool->classificationtime;
      totaltime += conshdlrdata->seeedpool->nblockscandidatescalctime;
      totaltime += conshdlrdata->seeedpool->postprocessingtime;
      totaltime += conshdlrdata->seeedpool->scorecalculatingtime;
      totaltime += conshdlrdata->seeedpool->translatingtime;
   }
   */

   return totaltime;
}


/*
 * @brief translates unpresolved seeed to a complete presolved one
 * @param scip SCIP data structure
 * @param success  at least one unpresolved seeed could not be translated in a complete presolved one
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompTranslateAndAddCompleteUnpresolvedSeeeds(
   SCIP*                 scip,        /* SCIP data structure */
   SCIP_Bool*            success      /* at least one unpresolved seeed could be translated in a complete presolved one */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* seeedpool;
   gcg::Seeedpool* seeedpoolunpresolved;
   std::vector<SeeedPtr> seeedstotranslate(0);
   std::vector<SeeedPtr> seeedstranslated(0);
   std::vector<SeeedPtr>::iterator seeediter;
   std::vector<SeeedPtr>::iterator seeediterend;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   *success = FALSE;

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);



   seeedpool = conshdlrdata->seeedpool;
   seeedpoolunpresolved = conshdlrdata->seeedpoolunpresolved;

   if( seeedpool == NULL )
      SCIPconshdlrDecompCreateSeeedpool(scip);

   seeedpool = conshdlrdata->seeedpool;

   assert(seeedpool != NULL);
   assert(seeedpoolunpresolved != NULL);

   for( int i = 0; i < seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
   {
      SeeedPtr finseeed = seeedpoolunpresolved->getFinishedSeeed(i);
      if( finseeed->isComplete() )
      {
         assert( finseeed->checkConsistency( ) );
         seeedstotranslate.push_back(finseeed);
      }
   }


   seeedpool->translateSeeeds(seeedpoolunpresolved, seeedstotranslate, seeedstranslated);

   seeediter = seeedstranslated.begin();
   seeediterend = seeedstranslated.end();


   for(; seeediter != seeediterend; ++seeediter )
   {
      seeedpool->prepareSeeed( *seeediter);
      if( (*seeediter)->isComplete() )
      {
         SCIP_CALL(SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, *seeediter ) );
         *success = TRUE;
      }
      else
      {
         (*seeediter)->completeByConnected();
         if ( (*seeediter)->isComplete() )
         {
            SCIP_CALL(SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, *seeediter ) );
            *success = TRUE;
         }
         else
            SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForPresolved(scip, *seeediter ) );
      }
   }

   return SCIP_OKAY;
}

/** method to adapt score for unpresolved decomps
 * @TODO: change score for some parameter settings
 *
 * @returns new score */
static
SCIP_Real SCIPconshdlrDecompAdaptScore(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_Real         oldscore /**< current score (to be updated) */
   )
{
   SCIP_Real score = oldscore;
   int method;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   gcg::Seeedpool* seeedpool;
   gcg::Seeedpool* seeedpoolunpresolved;

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeedpool = conshdlrdata->seeedpool;
   seeedpoolunpresolved = conshdlrdata->seeedpoolunpresolved;

   SCIP_CALL(SCIPgetIntParam(scip, "detection/origprob/weightinggpresolvedoriginaldecomps", &method) );

   if( method == FRACTION_OF_NNONZEROS )
   {
      if ( seeedpool == NULL || seeedpoolunpresolved == NULL )
         return score;

      score *= (SCIP_Real) seeedpoolunpresolved->getNNonzeros() / seeedpool->getNNonzeros();
   }

   if( method == FRACTION_OF_NROWS )
   {
      if ( seeedpool == NULL || seeedpoolunpresolved == NULL )
         return score;

      score *= (SCIP_Real) seeedpoolunpresolved->getNConss() / seeedpool->getNConss();

   }

   if( method == FAVOUR_PRESOLVED )
   {
      score += 1.;
   }

   return score;
}

/*
 * @brief returns whether or not there exists at least one (complete or incomplete) decomposition
 * @param scip SCIP data structure
 * @returns TRUE if there exists at least one (complete or incomplete) decomposition
 */
SCIP_Bool SCIPconshdlrDecompHasDecomp(
   SCIP*    scip     /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
         }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return ( (conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() > 0 )  ||
      ( conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNIncompleteSeeeds() > 0 ) ||
      ( conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->getNIncompleteSeeeds() > 0 )  ||
      ( conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() > 0 )      ) ;
}


/*
 * @brief initilizes the candidates data structures with selected seeeds
 * (or all if there are no selected seeeds) and sort them according to the current scoretype
 * @param scip SCIP data structure
 * @param updatelist whether or not the seeed list should be updated
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip             /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   std::vector<SeeedPtr>::iterator seeediter;
   std::vector<SeeedPtr>::iterator seeediterend;

   std::vector<SeeedPtr> tofinishpresolved(0);
   std::vector<SeeedPtr> tofinishunpresolved(0);
   std::vector<SeeedPtr> selectedseeeds(0);
   std::vector<SeeedPtr> finished(0);
   std::vector<SeeedPtr> finishedunpresolved(0);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot manage decompositions!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("Starting decomposition candidate choosing \n");

   assert(conshdlrdata->candidates != NULL);

   conshdlrdata->candidates->clear();

   assert( SCIPconshdlrDecompCheckConsistency(scip) );

   std::vector<SeeedPtr> selectedlist = getSelectedSeeeds(scip);
   for( size_t selid = 0; selid < selectedlist.size(); ++selid )
   {
      selectedseeeds.push_back(selectedlist.at(selid) );
   }

   if ( selectedseeeds.size() == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL, "currently no decomposition is selected, hence every known decomposition is considered: \n");
      selectedseeeds = getLeafSeeeds(scip);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL,  "number of considered decompositions: %d \n", selectedseeeds.size() );
   }

   /* if there are selected decomps, check if some of them needs to be finished and do so */
   seeediter = selectedseeeds.begin();
   seeediterend = selectedseeeds.end();

   for( ; seeediter != seeediterend; ++seeediter)
   {
      if( !(*seeediter)->isComplete() && (*seeediter)->isFromUnpresolved() )
      {
         tofinishunpresolved.push_back(*seeediter);
      }

      if( !(*seeediter)->isComplete() && !(*seeediter)->isFromUnpresolved() )
      {
         tofinishpresolved.push_back(*seeediter);
      }
   }

   finished = conshdlrdata->seeedpool->finishIncompleteSeeeds(tofinishpresolved);
   if( conshdlrdata->seeedpoolunpresolved != NULL )
      finishedunpresolved = conshdlrdata->seeedpoolunpresolved->finishIncompleteSeeeds(tofinishunpresolved);

   seeediter = selectedseeeds.begin();
   seeediterend = selectedseeeds.end();


   /* get decomp candidates and calculate corresponding score (possibly weighted for unpresolved) */
   for( ; seeediter != seeediterend; ++seeediter )
   {
      SeeedPtr seeed = *seeediter ;
      if( seeed->isComplete() && !seeed->isFromUnpresolved() )
      {
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, seeed->getScore(SCIPconshdlrDecompGetScoretype(scip)) ) );
      }
      if( seeed->isComplete() && seeed->isFromUnpresolved() )
      {
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, SCIPconshdlrDecompAdaptScore(scip, seeed->getScore(SCIPconshdlrDecompGetScoretype(scip)) ) ) );
      }
   }

   seeediter = finished.begin();
   seeediterend = finished.end();

   for( ; seeediter != seeediterend; ++seeediter )
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, (*seeediter)->getScore(SCIPconshdlrDecompGetScoretype(scip)) )  );
   }

   seeediter = finishedunpresolved.begin();
   seeediterend = finishedunpresolved.end();

   for( ; seeediter != seeediterend; ++seeediter )
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, SCIPconshdlrDecompAdaptScore(scip, (*seeediter)->getScore(SCIPconshdlrDecompGetScoretype(scip)) ) ) );
   }

   /* sort decomp candidates according score */
   std::sort( conshdlrdata->candidates->begin(), conshdlrdata->candidates->end(), sort_pred() );

   return SCIP_OKAY;
}

/*@todo cleanup: is SCIPconshdlrDecompAddLegacymodeDecompositions still used?*/
/**
 * @brief calls old detectStructure methods of chosen detectors, translates the resulting decompositions
 *  into seeeds and adds these seeeds to (presolved) seeedpool
 * @param scip SCIP data structure
 * @param result was the legacy call successful
 * @returns SCIP return code
 */
static
SCIP_RETCODE SCIPconshdlrDecompAddLegacymodeDecompositions(
   SCIP* scip,
   SCIP_RESULT* result
   )
{
   /* access to relevant data structures of the conshdlr */
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* seeedpool;

   /* detector data and their detection results */
   char detectorchaininfo[SCIP_MAXSTRLEN];
   int d;
   DEC_DETECTOR* detector;
   DEC_DECOMP** decdecomps;
   int ndecdecomps;
   SCIP_CLOCK* detectorclock;
   SCIP_RESULT decResult;

   /* decompositions and seeeds */
   gcg::SeeedPtr dummyAncestor;
   int dec;
   gcg::SeeedPtr seeed;
   int dupcount;
   SCIP_Bool legacyenabled;
   SCIP_Bool onlylegacy;


   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* check whether legacymode of at least one detector is enabled */
   SCIPgetBoolParam(scip, "detection/legacymode/enabled", &legacyenabled);
   SCIPgetBoolParam(scip, "detection/legacymode/onlylegacymode", &onlylegacy);

   if( !legacyenabled && !onlylegacy )
      return SCIP_OKAY;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Start legacy mode detection.\n");

   /* do transformations and initializations if necessary */
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb( scip ) );

   if ( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
      SCIP_CALL( SCIPpresolve( scip ) );

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   if ( conshdlrdata->seeedpool == NULL )
      conshdlrdata->seeedpool = new gcg::Seeedpool( scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip) );

   seeedpool = conshdlrdata->seeedpool;

   dummyAncestor = new gcg::Seeed( scip, seeedpool->getNewIdForSeeed(), seeedpool );
   seeedpool->addSeeedToAncestor( dummyAncestor );

   SCIPdebugMessagePrint(scip, "Checking %d detectors for legacy mode.\n", conshdlrdata->ndetectors);

   /* for each detector: check whether legacymode is enabled */
   for( d = 0; d < conshdlrdata->ndetectors; ++d )
   {
      decdecomps = NULL;
      ndecdecomps = -1;
      detector = conshdlrdata->detectors[d];
      assert(detector != NULL);

      if( detector->legacymode )
      {
         if( detector->detectStructure == NULL )
         {
            SCIPverbMessage( scip, SCIP_VERBLEVEL_NORMAL , NULL,
               "Legacy mode is not supported by detector <%s>.\n", detector->name );
         }
         else
         {
            SCIPverbMessage( scip, SCIP_VERBLEVEL_NORMAL , NULL,
               "Start legacy mode detection for detector <%s>.\n", detector->name );

            /* measure time detector needs for detecting decompositions */
            SCIPcreateClock( scip, & detectorclock );
            SCIP_CALL_ABORT( SCIPstartClock( scip, detectorclock ) );

            /* call old detectStructure callback method */
            SCIP_CALL( (*detector->detectStructure)( scip, detector->decdata, &decdecomps, &ndecdecomps, &decResult ) );

            SCIP_CALL_ABORT( SCIPstopClock( scip, detectorclock ) );

            if( decResult == SCIP_SUCCESS )
            {
               /* check for duplicates and redundant information */
               for( dec = 0; dec < ndecdecomps; ++dec )
               {
                  assert( decdecomps[dec] != NULL );
               }
               if( ndecdecomps > 2 )
               {
                  int nunique = DECfilterSimilarDecompositions(scip, decdecomps, ndecdecomps);

                  for( dec = nunique; dec < ndecdecomps; ++dec )
                  {
                     SCIP_CALL( DECdecompFree(scip, &(decdecomps[dec])) );
                     decdecomps[dec] = NULL;
                  }

                  ndecdecomps = nunique;
               }

               SCIPdebugMessagePrint( scip, "Translate %d non-redundant decompositions into seeeds.\n", ndecdecomps );

               /* set up detectorchaininfo */
               SCIPsnprintf( detectorchaininfo, SCIP_MAXSTRLEN, "%c(lgc)", detector->decchar );

               dupcount = 0;

               /* translate found decompositions to seeeds and add them to (presolved) seeedpool */
               for( dec = 0; dec < ndecdecomps; ++dec )
               {
                  seeedpool->createSeeedFromDecomp( decdecomps[dec], &seeed );

                  seeed->setDetectorChainString( detectorchaininfo );

                  /* set statistical data */
                  seeed->setDetectorPropagated(detector);
                  /* @todo this is actually the whole detector time! in propagate seeeds methods, time of each seeed
                   * is set after detecting this seeed */
                  seeed->addClockTime( SCIPgetClockTime( scip, detectorclock ) );
                  seeed->addDecChangesFromAncestor( dummyAncestor );
                  seeed->setLegacymode( true );

                  SCIP_Bool success = TRUE;
                  seeedpool->addSeeedToFinished( seeed, &success );

                  if( success == FALSE )
                  {
                     ++dupcount;
                  }
               }

               if ( dupcount > 0 )
               {
                  SCIPdebugMessagePrint( scip, "%d of the resulting seeeds are already contained in the seeedpool.\n", dupcount );
               }

               SCIPfreeClock( scip, & detectorclock );
            }
            else
            {
               SCIPdebugPrintf( "Failure!\n" );
            }
            SCIPfreeMemoryArrayNull( scip, &decdecomps ); // @todo necessary/correct?
         }
      }
   }

   seeedpool->sortFinishedForScore();

   SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Finished legacy mode detection.\n");

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}



/* Checks whether
 *  1) the predecessors of all finished seeeds in both seeedpools can be found
 *
 *  @returns true if seeed information is consistent */
SCIP_Bool SCIPconshdlrDecompCheckConsistency(
   SCIP* scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* 1) the predecessors of all finished seeeds in both seeedpools can be found */
   if( conshdlrdata->seeedpool != NULL)
   {
      for( i = 0; i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
      {
         SeeedPtr seeed = conshdlrdata->seeedpool->getFinishedSeeed( i );

         for( int j = 0; j < seeed->getNAncestors(); ++j )
         {
            int id = seeed->getAncestorID( j );
            if( SCIPconshdlrDecompGetSeeed(scip, id) == NULL )
            {
               SCIPwarningMessage(scip, "Warning: presolved seeed %d has an ancestor (id: %d) that is not found! \n", seeed->getID(), id );
               return FALSE;
            }
         }
      }
   }


   if( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      for( i = 0; i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
      {
         SeeedPtr seeed = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i );

         for( int j = 0; j < seeed->getNAncestors(); ++j )
         {
            int id = seeed->getAncestorID( j );
            if( SCIPconshdlrDecompGetSeeed(scip, id) == NULL )
            {
               SCIPwarningMessage(scip, "Warning: unpresolved seeed %d has an ancestor (id: %d) that is not found! \n", seeed->getID(), id );
               return FALSE;
            }
         }
      }
   }

   return TRUE;
}

/* Gets the next seeed id managed by cons_decomp
 * @returns the next seeed id managed by cons_decomp */
int SCIPconshdlrDecompGetNextSeeedID(
   SCIP* scip   /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   Seeed_Wrapper sw;
   int testnextid;
   bool found = false;
   while(found == false)
   {
      testnextid = conshdlrdata->seeedcounter + 1;
      GCGgetSeeedFromID(scip, &testnextid, &sw);
      conshdlrdata->seeedcounter = testnextid;
      if(sw.seeed == NULL)
      {
         found = true;
      }
   }
   return conshdlrdata->seeedcounter;
}


/*
 * sort the finished decompositions according to the currently chosen score in the according datastructures for the
 * presolved and original problem
 * @returns scip return code
 */
SCIP_RETCODE DECconshdlrDecompSortDecompositionsByScore(
   SCIP*       scip     /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpool != NULL )
      conshdlrdata->seeedpool->sortFinishedForScore();


   if( conshdlrdata->seeedpoolunpresolved != NULL )
      conshdlrdata->seeedpoolunpresolved->sortFinishedForScore();

   return SCIP_OKAY;
}

/* interface method to detect the structure including presolving
 * @returns SCIP return code */
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_RESULT*          result              /* Result pointer to indicate whether some structure was found */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool onlylegacymode;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedpool != NULL )
   {
      delete conshdlrdata->seeedpool;
      conshdlrdata->seeedpool = NULL;
   }

   *result = SCIP_DIDNOTRUN;

   if( SCIPgetNOrigVars(scip) == 0 && SCIPgetNOrigConss(scip) == 0 )
      return SCIP_OKAY;

   /* if the original problem should be solved, then no decomposition will be performed */
   if( GCGgetDecompositionMode(scip) == DEC_DECMODE_ORIGINAL )
      return SCIP_OKAY;

   SCIP_CALL(SCIPresetClock(scip, conshdlrdata->completedetectionclock));
   SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock));

   /* check whether only legacy mode should be executed */
   SCIPgetBoolParam(scip, "detection/legacymode/onlylegacymode", &onlylegacymode);

   SCIPdebugMessage("start only legacy mode? %s \n", (onlylegacymode ? "yes": "no") );
   if( !onlylegacymode )
   {
      std::vector<std::pair<int, int>> candidatesNBlocks(0); /* collection of different variable class distributions */
      std::vector<gcg::ConsClassifier*> consClassDistributions; /* collection of different constraint class distributions */
      std::vector<gcg::VarClassifier*> varClassDistributions; /* collection of different variable class distributions */

      std::vector<SCIP_CONS*> indexToCons; /* stores the corresponding scip constraints pointer */
      std::vector<gcg::SeeedPtr> seeedsunpresolved(0); /* seeeds that were found for the unpresolved problem */
      int i;
      SCIP_Bool presolveOrigProblem;
      SCIP_Bool calculateOrigDecomps;
      SCIP_Bool classifyOrig;

      /* indicate wether only for the unpresolved problem the detection shpuld take place */
      SCIP_Bool detectonlyorig;

      assert(scip != NULL);

      presolveOrigProblem = TRUE;
      detectonlyorig = FALSE;

      SCIPgetBoolParam(scip, "detection/origprob/enabled", &calculateOrigDecomps);
      SCIPgetBoolParam(scip, "detection/origprob/classificationenabled", &classifyOrig);

      /* get data of the seeedpool with original vars and conss */
      SCIPdebugMessage("is seeedpoolunpresolved not initialized yet but needed ? %s -> %s create it \n", (conshdlrdata->seeedpoolunpresolved == NULL ? "yes" : "no"), (conshdlrdata->seeedpoolunpresolved == NULL ? "" : "Do not")  );

      /* scip is not presolved yet => only detect for original problem */
      if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         detectonlyorig = TRUE;

      if ( conshdlrdata->seeedpoolunpresolved == NULL && ( classifyOrig || calculateOrigDecomps || detectonlyorig) )
         conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));         /*< seeedpool with original variables and constraints */


      SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock));

      SCIPdebugMessage("is stage < transformed ? %s -> do %s transformProb() ", (SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED ? "yes" : "no"), (SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED ? "" : "not")  );
      if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
         SCIP_CALL(SCIPtransformProb(scip));

      SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock));

      /* get block number candidates and conslcassifier for original problem*/
      if( classifyOrig || detectonlyorig )
      {
         SCIPdebugMessage("classification for orig problem enabled: calc classifier and nblock candidates \n" );
         conshdlrdata->seeedpoolunpresolved->calcClassifierAndNBlockCandidates(scip);
         candidatesNBlocks = conshdlrdata->seeedpoolunpresolved->getSortedCandidatesNBlocksFull();
         if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
                  conshdlrdata->seeedpoolunpresolved->printBlockcandidateInformation(scip, NULL);
      }
      else
         SCIPdebugMessage("classification for orig problem disabled \n" );

      /* detection for original problem */
      if( calculateOrigDecomps || detectonlyorig )
      {
         SCIPdebugMessage("start finding decompositions for original problem!\n" );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "start finding decompositions for original problem!\n");
         seeedsunpresolved = conshdlrdata->seeedpoolunpresolved->findSeeeds();
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "finished finding decompositions for original problem!\n");
         SCIPdebugMessage("finished finding decompositions for original problem!\n" );
      } else
         SCIPdebugMessage("finding decompositions for original problem is NOT enabled!\n" );

      /* get the cons and var classifier for translating them later*/
      if( classifyOrig )
      {
         for( i = 0; i < conshdlrdata->seeedpoolunpresolved->getNConsClassifiers(); ++i )
         {
            gcg::ConsClassifier* classifier = new gcg::ConsClassifier(
               conshdlrdata->seeedpoolunpresolved->getConsClassifier(i));
            consClassDistributions.push_back(classifier);
         }
         for( i = 0; i < conshdlrdata->seeedpoolunpresolved->getNVarClassifiers(); ++i )
         {
            gcg::VarClassifier* classifier = new gcg::VarClassifier(conshdlrdata->seeedpoolunpresolved->getVarClassifier(i));
            varClassDistributions.push_back(classifier);
         }
      }

      SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock));

      if( !detectonlyorig )
      {

         /*Presolving*/
         if( presolveOrigProblem )
            SCIP_CALL(SCIPpresolve(scip));

         /* detection for presolved problem */

         if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

            /* presolving removed all constraints or variables */
            if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
               conshdlrdata->hasrun = TRUE;

            *result = SCIP_DIDNOTRUN;
            return SCIP_OKAY;
         }

         /* start detection clocks */
         SCIP_CALL(SCIPresetClock(scip, conshdlrdata->detectorclock));
         SCIP_CALL(SCIPstartClock(scip, conshdlrdata->detectorclock));
         SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock));
         if( conshdlrdata->seeedpool == NULL )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "start creating seeedpool for current problem \n");
            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "created seeedpool for current problem, n detectors: %d \n", conshdlrdata->ndetectors);
         }
         else
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "seeedpool is not NULL \n");

         conshdlrdata->seeedpool->calcClassifierAndNBlockCandidates(scip);

         /* get block number candidates and translate orig classification and found seeeds (if any) to presolved problem */
         if( calculateOrigDecomps ||  classifyOrig )
         {
            std::vector<gcg::Seeed*> translatedSeeeds(0);
            std::vector<gcg::ConsClassifier*> translatedConsDistributions(0);
            std::vector<gcg::VarClassifier*> translatedVarDistributions(0);

            conshdlrdata->seeedpool->translateSeeedData(conshdlrdata->seeedpoolunpresolved, seeedsunpresolved,
               translatedSeeeds, consClassDistributions, translatedConsDistributions, varClassDistributions,
               translatedVarDistributions);

            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL , NULL, "number of translated original seeeds: %d \n " , translatedSeeeds.size() );

            conshdlrdata->seeedpool->populate(translatedSeeeds);

            for( size_t d = 0; d < translatedConsDistributions.size(); ++d )
               conshdlrdata->seeedpool->addConsClassifier(translatedConsDistributions[d]);

            for( size_t d = 0; d < translatedVarDistributions.size(); ++d )
               conshdlrdata->seeedpool->addVarClassifier(translatedVarDistributions[d]);

            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL , NULL, "finished translate seeed method!\n");

            for( size_t c = 0; c < candidatesNBlocks.size(); ++c )
               conshdlrdata->seeedpool->addCandidatesNBlocksNVotes(candidatesNBlocks[c].first, candidatesNBlocks[c].second );
         }
      }

     for( int j = 0; j < (int) consClassDistributions.size(); ++j )
        delete consClassDistributions[j];

     for( int j = 0; j < (int) varClassDistributions.size(); ++j )
        delete varClassDistributions[j];

     if( !detectonlyorig )
     {
        conshdlrdata->seeedpool->findDecompositions();
        conshdlrdata->seeedpool->sortFinishedForScore();
        SCIP_CALL(SCIPstopClock(scip, conshdlrdata->detectorclock));
     }

      if( conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() > 0 )
         *result = SCIP_SUCCESS;

      if( conshdlrdata->seeedpoolunpresolved != NULL &&  conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() > 0 )
         *result = SCIP_SUCCESS;

      SCIPdebugMessage("Detection took %fs\n", SCIPgetClockTime( scip, conshdlrdata->detectorclock));

   } /* end of if( !onlylegacy ) */


   if( conshdlrdata->seeedpool != NULL && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
      conshdlrdata->seeedpool->printBlockcandidateInformation(scip, NULL);

   SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock) );
   SCIPconshdlrDecompAddLegacymodeDecompositions( scip, result );
   SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock) );

   /* display timing statistics */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Detection Time: %.2f\n", SCIPconshdlrDecompGetCompleteDetectionTime(scip));
   /** @todo put this output to the statistics output */

   if( *result == SCIP_DIDNOTRUN )
   {
      return SCIP_OKAY;
   }

   /* show that we done our duty */
   conshdlrdata->hasrun = TRUE;
   *result = SCIP_SUCCESS;
   SCIPconshdlrDecompChooseCandidatesFromSelected(scip);

   return SCIP_OKAY;
}


/* writes all selected decompositions */
SCIP_RETCODE DECwriteSelectedDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension          /**< extension for decompositions */
   )
{
   MiscVisualization* misc = new MiscVisualization();
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char outname[SCIP_MAXSTRLEN];
   char tempstring[SCIP_MAXSTRLEN];
   SCIP_Bool nodecomps;

   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nodecomps = ( conshdlrdata->seeedpool == NULL && conshdlrdata->seeedpoolunpresolved == NULL );

   if( nodecomps )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   std::vector<SeeedPtr> selectedseeeds = getSelectedSeeeds(scip);
   if( selectedseeeds.size() == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition selected.\n");
      return SCIP_OKAY;
   }

   for( size_t selid = 0; selid < selectedseeeds.size(); ++selid )
   {
      SeeedPtr seeed;
      seeed = selectedseeeds.at(selid);

      misc->GCGgetVisualizationFilename(scip, seeed, extension, tempstring);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s.%s", directory, tempstring, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s.%s", tempstring, extension);
      }

      conshdlrdata->seeedtowrite = seeed;

      if ( seeed->isFromUnpresolved() )
         SCIP_CALL_QUIET( SCIPwriteOrigProblem(scip, outname, extension, FALSE) );
      else
         SCIP_CALL_QUIET( SCIPwriteTransProblem(scip, outname, extension, FALSE) );

      conshdlrdata->seeedtowrite = NULL;
   }

   return SCIP_OKAY;
}


/* write out all known decompositions
 * @returns SCIP return code  */
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /* SCIP data structure */
   char*                 directory,          /* directory for decompositions */
   char*                 extension,          /* extension for decompositions */
   SCIP_Bool             original,           /* should decomps for original problem be written */
   SCIP_Bool             presolved           /* should decomps for preoslved problem be written */

   )
{
   MiscVisualization* misc = new MiscVisualization();
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char outname[SCIP_MAXSTRLEN];
   char tempstring[SCIP_MAXSTRLEN];
   int i;
   SCIP_Bool nodecomps;

   int maxtowrite;
   int nwritten;

   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   maxtowrite = -1;
   nwritten = 0;

   nodecomps = ( conshdlrdata->seeedpool == NULL && conshdlrdata->seeedpoolunpresolved == NULL );

   nodecomps = nodecomps || ( !presolved && !original );

   nodecomps = nodecomps || (
      ( presolved && conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() == 0) &&
      ( original && conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() == 0)
       )
      ;


   if( presolved && conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   SCIPgetIntParam(scip, "visual/nmaxdecompstowrite", &maxtowrite );

   /* write presolved decomps */
   for( i = 0; presolved && conshdlrdata->seeedpool!= NULL && i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
   {
      SeeedPtr seeed;

      seeed = conshdlrdata->seeedpool->getFinishedSeeed( i );

      misc->GCGgetVisualizationFilename(scip, seeed, extension, tempstring);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s.%s", directory, tempstring, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s.%s", tempstring, extension);
      }

      conshdlrdata->seeedtowrite = seeed;

      SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );

      ++nwritten;

      conshdlrdata->seeedtowrite = NULL;

      if( maxtowrite != -1 && nwritten >= maxtowrite )
         break;


   }

   /* write orig decomps */
   for( i = 0; original && conshdlrdata->seeedpoolunpresolved != NULL &&  i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() ; ++i )
   {
      SeeedPtr seeed;

      seeed = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i );

      misc->GCGgetVisualizationFilename(scip, seeed, extension, tempstring);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s.%s", directory, tempstring, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s.%s", tempstring, extension);
      }

      conshdlrdata->seeedtowrite = seeed;

      SCIP_CALL_QUIET( SCIPwriteOrigProblem(scip, outname, extension, FALSE) );
      ++nwritten;
      conshdlrdata->seeedtowrite = NULL;

      if( maxtowrite != -1 && nwritten >= maxtowrite )
         break;
   }


     return SCIP_OKAY;
}


/* Gets whether the detection already took place
 * @returns true if detection took place, false otherwise */
SCIP_Bool GCGdetectionTookPlace(
   SCIP*  scip /* SCIP data structure */
     ){
   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;

    assert(scip != NULL);

    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert(conshdlr != NULL);

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    return (conshdlrdata->seeedpool != NULL ) || (conshdlrdata->seeedpoolunpresolved != NULL );
}


/* Gets the number of all detectors
 * @returns number of detectors */
int SCIPconshdlrDecompGetNDetectors(
   SCIP* scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;

    assert(scip != NULL);

    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert(conshdlr != NULL);

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    return conshdlrdata->ndetectors;
}


/* Gets an array of all detectors
 *
 * @returns array of detectors */
DEC_DETECTOR** SCIPconshdlrDecompGetDetectors(
   SCIP* scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;

    assert(scip != NULL);

    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert(conshdlr != NULL);

    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    return conshdlrdata->detectors;
}


/* @brief switch nonfinalfreetransform to TRUE
 *
 * used before calling SCIPfreeTransform(), if called to revoke presolving (e.g. if unpresolved decomposition is used, and
 * transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection
 * information should be freed
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompNotifyNonFinalFreeTransform(
   SCIP*                scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nonfinalfreetransform = TRUE;

   return SCIP_OKAY;
}


/* @brief switch nonfinalfreetransform to FALSE
 * used after calling SCIPfreeTransform() if called to revoke presolving (e.g. if unpresolved decomposition is used,
 * and transformation is not successful), this seems mandatory to decide during consExitDecomp if the original detection
 * information should be freed
 * @returns scip return code
 */
SCIP_RETCODE SCIPconshdlrDecompNotifyFinishedNonFinalFreeTransform(
   SCIP*                scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->nonfinalfreetransform = FALSE;

   return SCIP_OKAY;
}



/* gets an array of all seeeds that are currently considered relevant
 * @params seeedswr  output of the relevant seeeds (don't forget to free the individual wrappers after use)
 * @params nseeeds   amount of seeeds that are put in the array
 * @retruns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompGetAllRelevantSeeeds(
   SCIP* scip,                /* SCIP data structure */
   SEEED_WRAPPER** seeedswr,  /* seeed wrapper array for output of the relevant seeeds (don't forget to free the individual wrappers after use) */
   int* nseeeds               /* amount of seeeds that are put in the output array */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   assert(seeedswr != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* get the current max id */
   int maxid  = 0;

   for( int i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNAncestorSeeeds(); ++i )
   {
      if( conshdlrdata->seeedpool->getAncestorSeeed( i ) != NULL &&
         conshdlrdata->seeedpool->getAncestorSeeed( i )->getID() > maxid )
         maxid = conshdlrdata->seeedpool->getAncestorSeeed( i )->getID();
   }

   for( int i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNAncestorSeeeds(); ++i )
   {
      if( conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i ) != NULL &&
         conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )->getID() > maxid )
         maxid = conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )->getID();
   }

   for( int i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpool->getFinishedSeeed( i ) != NULL &&
            conshdlrdata->seeedpool->getFinishedSeeed( i )->getID() > maxid )
            maxid = conshdlrdata->seeedpool->getFinishedSeeed( i )->getID();
      }

      for( int i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ) != NULL &&
            conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i )->getID() > maxid )
            maxid = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i )->getID();
      }

   /* initialize the output array with NULL */
   *nseeeds = maxid+1;

   for( int i = 0; i < *nseeeds; i++ )
   {
      SCIP_CALL( SCIPallocBlockMemory( scip, &(seeedswr[i]) ) );
      seeedswr[i]->seeed = NULL;
      seeedswr[i]->seeedpool = NULL; /* not needed, initialization just for safety reasons */
   }

   /* fill the output array with relevant seeeds */
   for( int i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNAncestorSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i ) == NULL ||
            conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )->getID() < 0  )
            continue;
         seeedswr[conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i )->getID()]->seeed =
            conshdlrdata->seeedpoolunpresolved->getAncestorSeeed( i );
      }

   for( int i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNAncestorSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpool->getAncestorSeeed( i ) == NULL ||
            conshdlrdata->seeedpool->getAncestorSeeed( i )->getID() < 0  )
            continue;
         seeedswr[conshdlrdata->seeedpool->getAncestorSeeed( i )->getID()]->seeed =
            conshdlrdata->seeedpool->getAncestorSeeed( i );
      }

   for( int i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ) == NULL ||
            conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i )->getID() < 0  )
            continue;
         seeedswr[conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i )->getID()]->seeed =
            conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i );
      }

   for( int i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
      {
         if( conshdlrdata->seeedpool->getFinishedSeeed( i ) == NULL ||
            conshdlrdata->seeedpool->getFinishedSeeed( i )->getID() < 0  )
            continue;
         seeedswr[conshdlrdata->seeedpool->getFinishedSeeed( i )->getID()]->seeed =
            conshdlrdata->seeedpool->getFinishedSeeed( i );
      }

   return SCIP_OKAY;
}


/* write the seeed in conshdlrdata->seeedtowrite into the file as dec
 *
 * @returns SCIP return code */
SCIP_RETCODE SCIPconshdlrDecompWriteDec(
   SCIP*     scip,         /* SCIP data structure */
   FILE*     file,         /* file for output */
   SCIP_Bool transformed,  /* is the problem transformed yet */
   SCIP_RESULT* result     /* result of writing dec */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   Seeedpool* seeedpool;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   seeedpool = NULL;
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( transformed )
   {
      if (conshdlrdata->seeedpool == NULL )
         conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
      seeedpool = conshdlrdata->seeedpool;
   }
   else
   {
      if (conshdlrdata->seeedpoolunpresolved == NULL )
         conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
      seeedpool = conshdlrdata->seeedpoolunpresolved;
   }

   if( conshdlrdata->seeedtowrite != NULL )
   {
      conshdlrdata->seeedtowrite->writeAsDec(file, seeedpool, result);
      return SCIP_OKAY;
   }

   if( conshdlrdata->candidates->size() == 0 )
   {
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip);
   }


   if( conshdlrdata->candidates->size() == 0 )
   {
      SCIPwarningMessage(scip, "There are no candidate decompositions!\n");
      return SCIP_OKAY;
   }

   conshdlrdata->candidates->at(0).first->writeAsDec(file, seeedpool, result);


   return SCIP_OKAY;
}


/* Runs the bender detector to create a block matrix and outputs its visualization as .png file
 * @returns SCIP return code*/
SCIP_RETCODE SCIPconshdlrDecompWriteMatrix(
   SCIP*                 scip,               /* scip data structure */
   const char*           filename,           /* filename the output should be written to (including directory) */
   const char*           workfolder,         /* directory in which should be worked */
   SCIP_Bool             originalmatrix      /* should the original (or transformed) matrix be written */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   Seeedpool* seeedpool;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   seeedpool = NULL;
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "start creating seeedpool \n");
   if( !originalmatrix )
   {
      if (conshdlrdata->seeedpool == NULL )
         conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
      seeedpool = conshdlrdata->seeedpool;
   }
   else
   {
      if (conshdlrdata->seeedpoolunpresolved == NULL )
         conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
      seeedpool = conshdlrdata->seeedpoolunpresolved;
   }

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "finished creating seeedpool \n");
   SCIP_CALL( seeedpool->writeMatrix(filename, workfolder ) );

   return SCIP_OKAY;

}

/* Gets the best known decomposition
 * @note caller has to free returned DEC_DECOMP
 * @returns the decomposition if available and NULL otherwise */
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);

   DEC_DECOMP* decomp;
   gcg::Seeedpool* seeedpool;
   gcg::Seeedpool* seeedpoolunpresolved;
   SeeedPtr seeed;


   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if ( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      return NULL;

   if( conshdlrdata->seeedpool == NULL )
      conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));

   seeedpool = conshdlrdata->seeedpool;
   seeedpoolunpresolved = conshdlrdata->seeedpoolunpresolved;

   if( conshdlrdata->candidates->size() == 0 && conshdlrdata->useddecomp == NULL)
   {
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip);
      if (conshdlrdata->candidates->size() == 0)
         return NULL;
   }


   if( conshdlrdata->useddecomp != NULL )
      return conshdlrdata->useddecomp;

   seeed = conshdlrdata->candidates->at( 0 ).first;

   SCIPdebugMessage("In get bestdecomp\n");

   if( SCIPconshdlrDecompIsBestCandidateUnpresolved(scip) )
   {
      std::vector<SeeedPtr> seeedtotranslate(0);
      std::vector<SeeedPtr> translatedseeeds(0);

      seeedtotranslate.push_back(seeed);
      seeedpool->translateSeeeds(seeedpoolunpresolved, seeedtotranslate, translatedseeeds);
      seeed = translatedseeeds[0];
   }

   seeedpool->createDecompFromSeeed(seeed, &decomp);


   return decomp;
}

/* Gets the currently considered best seeed
 * @returns the Seeed of the best Seeed if available and seeedwrapper->seeed = NULL otherwise */
SCIP_RETCODE DECgetSeeedToWrite(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_Bool             transformed,        /* is the problem transformed yet */
   SEEED_WRAPPER*        seeedwrapper        /* seeed wrapper to output */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   int dec;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->seeedtowrite != NULL )
   {
      seeedwrapper->seeed = conshdlrdata->seeedtowrite;

      return SCIP_OKAY;
   }

   if( conshdlrdata->candidates->size() == 0 )
   {
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip);
   }

   if( conshdlrdata->candidates->size() == 0 )
   {
      SCIPwarningMessage(scip, "There are no candidate decompositions!\n");
      seeedwrapper->seeed = NULL;
      return SCIP_OKAY;
   }

   for( dec = 0; dec  < (int) conshdlrdata->candidates->size(); ++dec )
   {
      if ( conshdlrdata->candidates->at(dec).first->isFromUnpresolved() == !transformed )
         break;
   }
   if( dec !=  (int) conshdlrdata->candidates->size() )
      seeedwrapper->seeed = conshdlrdata->candidates->at(dec).first;
   else
      SCIPwarningMessage(scip, "There is no candidate decomposition for the %s problem we can write information for!\n",
         transformed ? "transformed" : "untransformed");

   return SCIP_OKAY;
}

/* writes out a list of all detectors
 * @returns nothing */
void DECprintListOfDetectors(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int ndetectors;
   int i;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ndetectors = conshdlrdata->ndetectors;

   SCIPdialogMessage(scip, NULL, " detector             char priority enabled  description\n");
   SCIPdialogMessage(scip, NULL, " --------------       ---- -------- -------  -----------\n");

   for( i = 0; i < ndetectors; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s", conshdlrdata->detectors[i]->name);
      SCIPdialogMessage(scip, NULL,  "    %c", conshdlrdata->detectors[i]->decchar);
      SCIPdialogMessage(scip, NULL,  " %8d", conshdlrdata->detectors[i]->priority);
      SCIPdialogMessage(scip, NULL,  " %7s", conshdlrdata->detectors[i]->enabled ? "TRUE" : "FALSE");
      SCIPdialogMessage(scip, NULL,  "  %s\n", conshdlrdata->detectors[i]->description);
   }
}


/* Gets whether the detection has been performed
 * @returns whether the detection has been performed */
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->hasrun;
}

/* Gets the character of the detector
 * @returns detector character */
char DECdetectorGetChar(
   DEC_DETECTOR*         detector            /* pointer to detector */
)
{
   if( detector == NULL )
     return '0';
   else
      return detector->decchar;
}

/* gets all currently finished decomps
 * @note: the array is allocated and needs to be freed after use!
 * @returns an array containg all finished decompositions
 * */
DEC_DECOMP** SCIPconshdlrDecompGetFinishedDecomps(
   SCIP*     scip /* SCIP data structure */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   DEC_DECOMP** decomps;

   int ndecomps;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   ndecomps = SCIPconshdlrDecompGetNFinishedDecomps(scip);

   SCIP_CALL_ABORT( SCIPallocMemoryArray(scip, &decomps, ndecomps) );

   if ( conshdlrdata->seeedpool != NULL )
   {
      for( int i = 0; i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
      {
         DEC_DECOMP* decomp;
         SCIP_CALL_ABORT(conshdlrdata->seeedpool->createDecompFromSeeed(conshdlrdata->seeedpool->getFinishedSeeed( i ),
            &decomp ) );

         decomps[i] = decomp;
      }
   }
   if ( conshdlrdata->seeedpoolunpresolved != NULL )
   {
      /*  atm presolved decomposition are allowed (for visualization, by family tree, write all decomps etc)*/
      for( int i = 0; i <  conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds(); ++i )
      {
         DEC_DECOMP* decomp;
         int offset = 0;
         conshdlrdata->seeedpoolunpresolved->createDecompFromSeeed(conshdlrdata->seeedpoolunpresolved->getFinishedSeeed( i ),
            &decomp );

         if( conshdlrdata->seeedpool != NULL )
            offset = conshdlrdata->seeedpool->getNFinishedSeeeds();
         decomps[i + offset] = decomp;
      }
   }
   return decomps;
}

/* Gets the number of all finished Seeeds
 * @returns number of finished Seeeds */
int SCIPconshdlrDecompGetNFinishedDecomps(
   SCIP*       scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* case if there is no seeedpool in data yet */
   if( conshdlrdata->seeedpool == NULL && conshdlrdata->seeedpoolunpresolved == NULL )
      return 0;

   if( conshdlrdata->seeedpool == NULL )
      return conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds();

   if( conshdlrdata->seeedpoolunpresolved == NULL )
         return conshdlrdata->seeedpool->getNFinishedSeeeds();

   /* all other cases */
   return (int) conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() + conshdlrdata->seeedpool->getNFinishedSeeeds();
}

/* Gets the number of all seeeds
 * @returns number of Seeeds */
int SCIPconshdlrDecompGetNSeeeds(
   SCIP*       scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int nseeeds;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   nseeeds = 0;
   if( conshdlrdata->seeedpoolunpresolved != NULL )
      nseeeds += (int)
         conshdlrdata->seeedpoolunpresolved->getNAncestorSeeeds() +
         conshdlrdata->seeedpoolunpresolved->getNCurrentSeeeds() +
         conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds();



   if( conshdlrdata->seeedpool != NULL )
      nseeeds += (int)
      conshdlrdata->seeedpool->getNAncestorSeeeds() +
      conshdlrdata->seeedpool->getNCurrentSeeeds() +
      conshdlrdata->seeedpool->getNFinishedSeeeds();

   return nseeeds;
}

/* display statistics about detectors
 * @returns SCIP return code */
SCIP_RETCODE GCGprintDetectorStatistics(
   SCIP*                 scip,               /* SCIP data structure */
   FILE*                 file                /* output file or NULL for standard output */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int j;
   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);
   assert(scip != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Detector statistics:       time     number     blocks\n");
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file,
         "  %-10.10s       :   %8.2f %10d    ", conshdlrdata->detectors[i]->name, conshdlrdata->detectors[i]->dectime,
         conshdlrdata->detectors[i]->ndecomps );
      for( j = 0; j < conshdlrdata->detectors[i]->ndecomps; ++j )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %d",
            DECdecompGetNBlocks(conshdlrdata->detectors[i]->decomps[j]));
      }
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
   }
   return SCIP_OKAY;
}

/** resets the parameters to their default value
 * @returns SCIP return code */
static
SCIP_RETCODE setDetectionDefault(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   int i;
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL (SCIPsetIntParam(scip, "detection/maxrounds", 2) );
   SCIP_CALL (SCIPsetBoolParam(scip, "detection/origprob/enabled", FALSE) );

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER )
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", FALSE) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIP_Result result;

      char paramname[SCIP_MAXSTRLEN];
      SCIP_Bool paramval;
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
         "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPresetParam(scip, paramname) );

      result = SCIP_DIDNOTRUN;
      if( conshdlrdata->detectors[i]->setParamDefault != NULL )
         conshdlrdata->detectors[i]->setParamDefault(scip, conshdlrdata->detectors[i], &result);
      if( !quiet )
      {
         SCIP_Bool written = FALSE;

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         if( written )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

      }
   }

   return SCIP_OKAY;
}

/** sets the parameters to aggressive values
 *
 * @returns SCIP return code */
static
SCIP_RETCODE setDetectionAggressive(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);


   SCIP_CALL (SCIPsetIntParam(scip, "detection/maxrounds", 3) );

   SCIP_CALL (SCIPsetBoolParam(scip, "detection/origprob/enabled", TRUE) );

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER)
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", FALSE) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
      {
         SCIP_Result result;

         result = SCIP_DIDNOTRUN;
         if( conshdlrdata->detectors[i]->setParamAggressive != NULL )
            conshdlrdata->detectors[i]->setParamAggressive(scip, conshdlrdata->detectors[i], &result);

         if( !quiet )
         {
            char paramname[SCIP_MAXSTRLEN];
            SCIP_Bool paramval;
            SCIP_Bool written = FALSE;

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
               "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
               "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
               "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            if( written )
               SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

         }
      }

   return SCIP_OKAY;
}

/** disables detectors
 *
 * @returns SCIP return code */
static
SCIP_RETCODE setDetectionOff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715}*/
   int i;
   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPsetBoolParam(scip, paramname, FALSE) );
      if( !quiet )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = FALSE\n", paramname);
      }
   }

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void)SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
         "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL(SCIPsetBoolParam(scip, paramname, FALSE));
      if( !quiet )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = FALSE\n", paramname);
      }
   }

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void)SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detection/detectors/%s/legacymode", conshdlrdata->detectors[i]->name);

      SCIP_CALL(SCIPsetBoolParam(scip, paramname, FALSE));
      if( !quiet )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = FALSE\n", paramname);
      }
   }

   return SCIP_OKAY;
}

/** sets the parameters to fast values
 *
 * @returns SCIP return code  */
static
SCIP_RETCODE setDetectionFast(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLRDATA*    conshdlrdata,       /**< constraint handler data structure */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
   )
{ /*lint --e{715} */
   int i;

   assert(scip != NULL);
   assert(conshdlrdata != NULL);

   SCIP_CALL (SCIPsetIntParam(scip, "detection/maxrounds", 1) );


   SCIP_CALL (SCIPsetBoolParam(scip, "detection/origprob/enabled", FALSE) );
   SCIP_CALL (SCIPsetBoolParam(scip, "detection/origprob/classificationenabled", FALSE) );


   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER )
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", FALSE) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIP_Result result;

      result = SCIP_DIDNOTRUN;
      if( conshdlrdata->detectors[i]->overruleemphasis )
         continue;

      if( conshdlrdata->detectors[i]->setParamFast != NULL )
         conshdlrdata->detectors[i]->setParamFast(scip, conshdlrdata->detectors[i], &result);
      if( !quiet )
      {
         char paramname[SCIP_MAXSTRLEN];
         SCIP_Bool paramval;
         SCIP_Bool written = FALSE;

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         if( written )
            SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
      }
   }

   return SCIP_OKAY;
}

/* sets detector parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all detector parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for detection is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the detectors produce more decompositions
 *  - SCIP_PARAMSETTING_OFF which turns off all detection
 *
 *   @returns SCIP return code
 */
SCIP_RETCODE GCGsetDetection(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /* parameter settings */
   SCIP_Bool             quiet               /* should the parameter be set quiet (no output) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   switch( paramsetting )
   {
   case SCIP_PARAMSETTING_AGGRESSIVE:
      SCIP_CALL( setDetectionAggressive(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_OFF:
      SCIP_CALL( setDetectionOff(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_FAST:
      SCIP_CALL( setDetectionFast(scip, conshdlrdata, quiet) );
      break;
   case SCIP_PARAMSETTING_DEFAULT:
      SCIP_CALL( setDetectionDefault(scip, conshdlrdata, quiet) );
      break;
   default:
      SCIPerrorMessage("The given paramsetting is invalid!\n");
      break;
   }

   return SCIP_OKAY;
}


/* returns wrapped Seeed with given id
 *
 * @returns SCIP status */
SCIP_RETCODE GCGgetSeeedFromID(
   SCIP*          scip,       /* SCIP data structure */
   int*           seeedid,    /* id of Seeed */
   SEEED_WRAPPER* seeedwr     /* wrapper for output Seeed */
   )
{
   SeeedPtr s;

   s = SCIPconshdlrDecompGetSeeed(scip, *seeedid);
   seeedwr->seeed = s;

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPconshdlrDecompRefineAndAddSeeed(
  SCIP* scip,
  Seeed_Wrapper* sw
  )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   char const* usergiveninfo;
   char const* presolvedinfo;

   SeeedPtr seeed = sw->seeed;
   assert( seeed != NULL );

   /* get seeedpool (create only does something if seeedpool does not exist) */
   Seeedpool* currseeedpool;
   if(seeed->isFromUnpresolved())
   {
      SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);
      currseeedpool = conshdlrdata->seeedpoolunpresolved;
   }
   else
   {
      SCIPconshdlrDecompCreateSeeedpool(scip);
      currseeedpool = conshdlrdata->seeedpool;
   }

   if(seeed->getSeeedpool() == NULL)
   {
      seeed->setSeeedpool(currseeedpool);
   }

   seeed->flushBooked();

   if( seeed->shouldCompletedByConsToMaster() )
   {
      for(int opencons = 0; opencons < seeed->getNOpenconss(); ++opencons)
         seeed->bookAsMasterCons( seeed->getOpenconss()[opencons] );
      seeed->flushBooked();
   }

   currseeedpool->prepareSeeed(seeed);

   if( !seeed->checkConsistency() )
   {
      delete seeed;
      SCIPwarningMessage(scip, "Your seeed was rejected because of inconsistencies! \n");
      return SCIP_OKAY;
   }
   seeed->buildDecChainString();
   if( seeed->isComplete() )
   {
      if( !seeed->shouldCompletedByConsToMaster() )
         seeed->setUsergiven( USERGIVEN::COMPLETE );

      if( !seeed->isFromUnpresolved() )
      {
         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, seeed));
      }
      /* stems from unpresolved problem */
      else
      {

         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(scip, seeed) );

         /* if seeedpool for presolved problem already exist try to translate seeed */
         if ( conshdlrdata->seeedpool != NULL )          {
            std::vector<Seeed*> seeedtotranslate(0);
            std::vector<Seeed*> newseeeds(0);
            seeedtotranslate.push_back(seeed);
            conshdlrdata->seeedpool->translateSeeeds(conshdlrdata->seeedpoolunpresolved, seeedtotranslate, newseeeds);
            if( newseeeds.size() != 0 )
            {
               SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, newseeeds[0]) );
            }
         }
      }
   }
   else
   {
      assert( !seeed->shouldCompletedByConsToMaster() );
      seeed->setUsergiven( USERGIVEN::PARTIAL );

      if ( !seeed->isFromUnpresolved() )
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForPresolved(scip, seeed) );
      else
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForUnpresolved(scip, seeed) );
   }

   /* set statistics */
   {
      int nvarstoblock = 0;
      int nconsstoblock = 0;

      for ( int b = 0; b < seeed->getNBlocks(); ++b )
      {
         nvarstoblock += seeed->getNVarsForBlock(b);
         nconsstoblock += seeed->getNConssForBlock(b);
      }
      seeed->setDetectorPropagated(NULL);

      seeed->addClockTime(0.);
      seeed->addPctVarsFromFree( (nvarstoblock + seeed->getNMastervars() +seeed->getNLinkingvars())/(SCIP_Real) seeed->getNVars()  );
      seeed->addPctVarsToBlock((nvarstoblock )/(SCIP_Real) seeed->getNVars() );
      seeed->addPctVarsToBorder( (seeed->getNMastervars() +seeed->getNLinkingvars())/(SCIP_Real) seeed->getNVars() ) ;
      seeed->addPctConssToBorder( (seeed->getNMasterconss() ) / (SCIP_Real) seeed->getNConss() ) ;
      seeed->addPctConssFromFree( (seeed->getNMasterconss() + nconsstoblock ) / (SCIP_Real) seeed->getNConss() ) ;
      seeed->addPctConssToBlock( (nconsstoblock ) / (SCIP_Real) seeed->getNConss() );
      seeed->addNNewBlocks(seeed->getNBlocks());
   }

   seeed->findVarsLinkingToMaster();
   seeed->findVarsLinkingToStairlinking();


   if( seeed->getUsergiven() == USERGIVEN::PARTIAL )
      usergiveninfo = "partial";
   if( seeed->getUsergiven() == USERGIVEN::COMPLETE )
      usergiveninfo = "complete";
   if( seeed->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER )
         usergiveninfo = "complete";
   if( seeed->isFromUnpresolved() )
         presolvedinfo = "unpresolved";
   else presolvedinfo = "presolved";

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " added %s decomp for %s problem with %d blocks and %d masterconss, %d linkingvars, "
      "%d mastervars, and max white score of %s %f \n", usergiveninfo, presolvedinfo,
      seeed->getNBlocks(), seeed->getNMasterconss(),
      seeed->getNLinkingvars(), seeed->getNMastervars(), (seeed->isComplete() ? " " : " at best "),
      seeed->getScore(SCORETYPE::MAX_WHITE) );

   return SCIP_OKAY;
}


/* prints blockcandiateinformation in following format:
 * NCANDIDATES
 * CANDIDATE : NVOTES for each candidate
 *
 * @returns SCIP status
 */
SCIP_RETCODE GCGprintBlockcandidateInformation(
   SCIP*                 scip,               /* SCIP data structure */
   FILE*                 file                /* output file or NULL for standard output */
)
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* seeedpool;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeedpool = (conshdlrdata->seeedpool == NULL ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );

   if( seeedpool == NULL )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), NULL, "No block number candidates are calculated yet, consider detecting first..  \n" );
   else
      seeedpool->printBlockcandidateInformation(scip, file);

   return SCIP_OKAY;
}


/* Prints complete detection time to given file
 *
 * @returns SCIP status */
SCIP_RETCODE GCGprintCompleteDetectionTime(
 SCIP*                 givenscip,          /* SCIP data structure */
 FILE*                 file                /* output file or NULL for standard output */
)
{
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "DETECTIONTIME   \n" );
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%f \n", (SCIP_Real) SCIPconshdlrDecompGetCompleteDetectionTime(givenscip) );

   return SCIP_OKAY;
}



/* prints classifier information in following format:
 * NCLASSIFIER
 * CLASSIFIERNAME  for each classifier
 * NCLASSES
 * CLASSNAME  for each class
 * NMEMBERS
 *
 * @returns SCIP status
 */
SCIP_RETCODE GCGprintClassifierInformation(
   SCIP*                 scip,               /* SCIP data structure */
   FILE*                 file                /* output file or NULL for standard output */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* seeedpool;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   seeedpool = (conshdlrdata->seeedpool == NULL ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );

   seeedpool->printClassifierInformation(scip, file);

   return SCIP_OKAY;
}


/* prints block candiate information in following format:
 * NDECOMPS
 * NBLOCKS  for each decomp
 * NCONSS for each block
 * NVARS
 * end for each block
 * NMASTERCONSS
 * NLINKINGVARS
 * NMASTERVARS
 * NSTAIRLINKINGVARS
 * MAXWHITESCORE
 * CLASSICALSCORE
 * HASSETPARTITIONINGMASTER
 * NDETECTORS
 * DETECTORNAME for each detector
 * NCONSCLASSIFIERS
 * CONSCLASSIFIERNAME for each consclassifier
 * nCLASSESMASTER
 * CLASSNAME for each class
 * NVARCLASSIFIERS
 * VARCLASSIFIERNAME for each varclassifier
 * nCLASSESMASTER
 * CLASSNAME for each class
 *
 * @returns SCIP status
 */
SCIP_RETCODE GCGprintDecompInformation(
   SCIP*                 scip,   /* SCIP data structure */
   FILE*                 file    /* output file or NULL for standard output */
   )
{
   std::vector<gcg::Seeed*>::const_iterator seeediter;
   std::vector<gcg::Seeed*>::const_iterator seeediterend;

   assert( SCIPconshdlrDecompCheckConsistency(scip) );

   std::vector<SeeedPtr> seeedlist = getLeafSeeeds(scip);
   seeediter = seeedlist.begin();
   seeediterend = seeedlist.end();

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "DECOMPINFO  \n" );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (int) seeedlist.size() );

   for( ; seeediter != seeediterend; ++seeediter)
   {
      gcg::Seeed* seeed;
      int nblocks = (*seeediter)->getNBlocks();

      seeed = *seeediter;

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "NEWDECOMP  \n" );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (*seeediter)->getNBlocks() );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (*seeediter)->getID() );
      for( int block = 0; block < nblocks; ++block )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNConssForBlock(block) );
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNVarsForBlock(block) );
      }
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNMasterconss() );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNLinkingvars() );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNMastervars() );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", seeed->getNTotalStairlinkingvars() );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  seeed->getMaxWhiteScore() );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  seeed->getScore(scoretype::CLASSIC) );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  seeed->getScore(scoretype::MAX_FORESSEEING_WHITE) );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n",  seeed->hasSetppccardMaster() );

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (int) seeed->getDetectorchainVector( ).size() );

      for( int detector = 0; detector <(int) seeed->getDetectorchainVector( ).size(); ++ detector )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s\n",
            DECdetectorGetName(seeed->getDetectorchainVector( )[detector]) );
      }
      seeed->printClassifierInformation(scip, file);
   }

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPconshdlrDecompSetScoretype(
   SCIP*  scip,
   SCORETYPE sctype
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->currscoretype = sctype;

   return SCIP_OKAY;
}


SCORETYPE SCIPconshdlrDecompGetScoretype(
   SCIP* scip
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return static_cast<scoretype>(conshdlrdata->currscoretype);
}


char* SCIPconshdlrDecompGetScoretypeShortName(
   SCIP*       scip,
   SCORETYPE   sctype
   )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;

   switch(sctype)
   {
   case scoretype::MAX_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maxwhi");
      break;
   case scoretype::CLASSIC:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classi");
      break;
   case scoretype::BORDER_AREA:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "border");
      break;
   case scoretype::MAX_FORESSEEING_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "forswh");
      break;
   case scoretype::MAX_FORESEEING_AGG_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "fawh");
      break;
   case scoretype::SETPART_FWHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfwh");
      break;
   case scoretype::SETPART_AGG_FWHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfawh");
      break;
   case scoretype::BENDERS:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "bender");
      break;
   case scoretype::STRONG_DECOMP:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "strode");
      break;
   default:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "");
   }

   SCIP_CALL_ABORT( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN) );
   return copy;
}


char* SCIPconshdlrDecompGetScoretypeDescription(
   SCIP*       scip,
   SCORETYPE   sctype
   )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;

   switch(sctype)
   {
   case scoretype::MAX_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum white area score (white area is nonblock and nonborder area)");
      break;
   case scoretype::CLASSIC:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classical score");
      break;
   case scoretype::BORDER_AREA:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "minimum border score (i.e. minimizes fraction of border area score)");
      break;
   case scoretype::MAX_FORESSEEING_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing white area score (considering copied linking vars and their master conss; white area is nonblock and nonborder area)");
      break;
   case scoretype::MAX_FORESEEING_AGG_WHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing white area score with aggregation infos (considering copied linking vars and their master conss; white area is nonblock and nonborder area)");
      break;
   case scoretype::SETPART_FWHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing white area score (convex combination of maximum foreseeing white area score and rewarding if master contains only setppc and cardinality constraints)");
      break;
   case scoretype::SETPART_AGG_FWHITE:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing white area score with aggregation information (convex combination of maximum foreseeing white area score and rewarding if a master contains only setppc and cardinality constraints)");
      break;
   case scoretype::BENDERS:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "experimental score to evaluate benders decompositions");
      break;
   case scoretype::STRONG_DECOMP:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "strong decomposition score");
      break;
   default:
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "");
      break;
   }

   SCIP_CALL_ABORT( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN ) );
   return copy;
}


SCIP_RETCODE SCIPconshdlrDecompGetSeeedLeafList(
   SCIP*          scip,
   int**          idlist,
   int*           listlength
   )
{
   std::vector<SeeedPtr> seeeds = getLeafSeeeds(scip);

   *listlength = (int) seeeds.size();

   for(int i = 0; i < (int) seeeds.size(); i++)
   {
      (*idlist)[i] = seeeds[i]->getID();
   }

   return SCIP_OKAY;
}

int SCIPconshdlrDecompGetNSeeedLeafs(
   SCIP*          scip
   )
{
   std::vector<SeeedPtr> seeeds = getLeafSeeeds(scip);
   return (int) seeeds.size();
}


SCIP_RETCODE SCIPconshdlrDecompGetSelectedSeeeds(
   SCIP*          scip,  /* SCIP data structure */
   int**          idlist,
   int*           listlength
   )
{
   /* get list of selected seeeds */
   std::vector<SeeedPtr> selectedseeeds = getSelectedSeeeds(scip);
   /* set the length of the pointer array to the list size */
   *listlength = (int) selectedseeeds.size();

   /* build a pointer list of ids from the seeed list */
   for(int i = 0; i < (int) selectedseeeds.size(); i++)
   {
      (*idlist)[i] = selectedseeeds[i]->getID();
   }

   return SCIP_OKAY;
}


SCIP_Bool SCIPconshdlrDecompGetSelectExists(
   SCIP*          scip  /* SCIP data structure */
   )
{
   /* get list of selected seeeds */
   std::vector<SeeedPtr> selectedseeeds = getSelectedSeeeds(scip);
   
   /* determine whether its size is 0*/
   return (selectedseeeds.size() == 0) ? false : true;
}


/* gets number of seeeds for public interface (pub_decomp.h)*/
int DECgetNDecomps(
   SCIP* scip
   )
{
   return SCIPconshdlrDecompGetNDecdecomps(scip);
}


int GCGgetNBlocksBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its blocks */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNBlocks();
}


int GCGgetNMasterConssBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its master conss */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNMasterconss();
}


int GCGgetNMasterVarsBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its master vars */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNMastervars();
}


int GCGgetNLinkingVarsBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its linking vars */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNLinkingvars();
}


int GCGgetNStairlinkingVarsBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its stairlinking vars */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNTotalStairlinkingvars();
}


float GCGgetScoreBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns its score in respect to the current score type */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getScore(SCIPconshdlrDecompGetScoretype(scip));
}


char* GCGgetDetectorHistoryBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns its detector history */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   assert(seeed != NULL);
   return seeed->getDetectorChainString();
}

SCIP_Bool GCGisPresolvedBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns whether it is presolved */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return !(seeed->isFromUnpresolved());
}


int GCGgetNOpenConssBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its open conss */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNOpenconss();
}


int GCGgetNOpenVarsBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns the number of its open vars */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->getNOpenvars();
}


SCIP_Bool GCGisSelectedBySeeedId(
   SCIP* scip,
   int id
   )
{
   /* get seeed and returns whether it is currently selected */
   Seeed* seeed = SCIPconshdlrDecompGetSeeed(scip, id);
   return seeed->isSelected();
}