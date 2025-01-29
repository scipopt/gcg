/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
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

/**@file   cons_decomp.cpp
 * 
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Michael Bastubbe
 * @author Hanna Franzen
 * @author William Ma
 *
 * This constraint handler will run all registered structure detectors in a loop. They will find partial decompositions in a loop iteration until the decompositions are full
 * or the maximum number of detection rounds is reached.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <deque>
#include <iomanip>
#include <sstream>
#include <utility>
#include <regex>
#include <vector>
#include <list>

#include <scip/clock.h>
#include <scip/pub_cons.h>
#include <scip/pub_misc.h>
#include <scip/type_paramset.h>
#include <scip/scipdefplugins.h>

#include "class_partialdecomp.h"
#include "class_detprobdata.h"
#include "miscvisualization.h"
#include "wrapper_partialdecomp.h"
#include "scip_misc.h"
#include "relax_gcg.h"
#include "decomp.h"
#include "cons_decomp.hpp"
#include "struct_consclassifier.h"
#include "struct_varclassifier.h"
#include "struct_decomp.h"
#include "struct_score.h"

/* constraint handler properties */
#define CONSHDLR_NAME                                 "decomp"    /**< name of constraint handler */
#define CONSHDLR_DESC                                 "constraint handler for structure detection"   /**< description of constraint handler */
#define CONSHDLR_ENFOPRIORITY                         0           /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY                        0           /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ                            -1          /**< frequency for using all instead of only the useful constraints in separation,
                                                                   *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS                            FALSE       /**< should the constraint handler be skipped, if no constraints are available? */

#define DEFAULT_ENABLED                               TRUE        /**< indicates whether detection is enabled */

#define DEFAULT_BLOCKNUMBERCANDSMEDIANVARSPERCONS     FALSE       /**< should for block number candidates calculation the medianvarspercons calculation be considered */

#define DEFAULT_MAXDETECTIONROUNDS                    1           /**< maximal number of detection rounds */
#define DEFAULT_MAXDETECTIONTIME                      600         /**< maximum detection time in seconds */
#define DEFAULT_POSTPROCESS                           TRUE        /**< indicates whether to postprocess full decompositions */
#define DEFAULT_MAXNCLASSESLARGEPROBS                 5           /**< maximum number of classes allowed for large (nvars+nconss > 50000) MIPs for detectors, partitions with more classes are reduced to the maximum number of classes */
#define DEFAULT_MAXNCLASSES                           9           /**< maximum number of classes allowed for detectors, partitions with more classes are reduced to the maximum number of classes */
#define DEFAULT_MAXNCLASSESFORNBLOCKCANDIDATES        18          /**< maximum number of classes a partition can have to be used for voting nblockcandidates */
#define DEFAULT_ENABLEORIGDETECTION                   TRUE        /**< indicates whether to start detection for the original problem */

/* classifier */

#define DEFAULT_ALLOWPARTITIONDUPLICATES             FALSE       /**< if false each new (conss- and vars-) partitions is checked for being a duplicate of an existing one, if so it is not added and NBOT statistically recognized*/
#define DEFAULT_CLASSIFY                              TRUE        /**< indicates whether classification is enabled */
#define DEFAULT_ENABLEORIGCLASSIFICATION              TRUE        /**< indicates whether to start detection for the original problem */

#define DEFAULT_BENDERSONLYCONTSUBPR                  FALSE       /**< indicates whether only decomposition with only continuous variables in the subproblems should be searched*/
#define DEFAULT_BENDERSONLYBINMASTER                  FALSE       /**< indicates whether only decomposition with only binary variables in the master should be searched */

#define DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER    10000       /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis default */
#define AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER 80000       /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis aggressive */
#define FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER       2000        /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis fast */

#define DEFAULT_DETECTBENDERS                         FALSE       /**< indicates whether benders detection mode is enabled */

#define DEFAULT_SCORE                                 "spfwh"      /**< score shortname that is activated */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             enabled;                                 /**< indicates whether detection is enabled */

   std::vector<PARTIALDECOMP*>*   partialdecs;                    /**< list of all existing partialdecs */
   std::unordered_map<int, PARTIALDECOMP*>*   partialdecsbyid;    /**< list of all existing partialdecs */
   int                   partialdeccounter;                       /**< counts the number of created partialdecs, used to determine next partialdec id
                                                                       (NOT amount of currently existing partialdecs as partialdecs might be deleted) */
   GCG_DECOMP**          decomps;                                 /**< array of decomposition structures */
   int                   ndecomps;                                /**< number of decomposition structures (size of decomps)*/

   GCG_CONSCLASSIFIER**  consclassifiers;                         /**< array of cons classifiers */
   int                   nconsclassifiers;                        /**< number of cons classifiers */
   int*                  consclassifierpriorities;                /**< priorities of the cons classifiers */
   GCG_VARCLASSIFIER**   varclassifiers;                          /**< array of var classifiers */
   int                   nvarclassifiers;                         /**< number of var classifiers */
   int*                  varclassifierpriorities;                 /**< priorities of the var classifiers */

   GCG_DETECTOR**        detectors;                               /**< array of structure detectors */
   int*                  priorities;                              /**< priorities of the detectors */
   int                   ndetectors;                              /**< number of detectors */
   GCG_DETECTOR**        propagatingdetectors;                    /**< array of detectors able to propagate partial decompositions */
   int                   npropagatingdetectors;                   /**< number of detectors able to propagate partial decompositions (size of propagatingdetectors) */
   GCG_DETECTOR**        finishingdetectors;                      /**< array of detectors able to finish partial decompositions */
   int                   nfinishingdetectors;                     /**< number of detectors able to finish partial decompositions (size of finishingdetectors) */
   GCG_DETECTOR**        postprocessingdetectors;                 /**< array of detectors able to postprocess decompositions */
   int                   npostprocessingdetectors;                /**< number of detectors able to postprocess decompositions (size of postprocessingdetectors) */

   SCIP_CLOCK*           detectorclock;                           /**< clock to measure detection time */
   SCIP_CLOCK*           completedetectionclock;                  /**< clock to measure detection time */
   SCIP_Bool             hasrunoriginal;                          /**< flag to indicate whether we have already detected (original problem) */
   SCIP_Bool             hasrun;                                  /**< flag to indicate whether we have already detected */
   int                   maxndetectionrounds;                     /**< maximum number of detection loop rounds */
   int                   maxdetectiontime;                        /**< maximum detection time [sec] */
   SCIP_Bool             postprocess;                             /**< indicates whether to postprocess full decompositions */
   SCIP_Bool             blocknumbercandsmedianvarspercons;       /**< should for block number candidates calculation the medianvarspercons calculation be considered */
   int                   maxnclassesfornblockcandidates;          /**< maximum number of classes a partition can have to be used for voting nblockcandidates */
   int                   maxnclassesperpartition;                /**< maximum number of classes allowed for detectors, partition with more classes are reduced to the maximum number of classes */
   int                   maxnclassesperpartitionforlargeprobs;   /**< maximum number of classes allowed for large (nvars+nconss > 50000) MIPs for detectors, partition with more classes are reduced to the maximum number of classes */
   int                   weightinggpresolvedoriginaldecomps;      /**< weighing method for comparing presolved and original decompositions (see corresponding enum)   */

   SCIP_Bool             classify;                                /**< indicates whether classification should take place */
   SCIP_Bool             allowpartitionduplicates;               /**< indicates whether partition duplicates are allowed (for statistical reasons) */
   SCIP_Bool             enableorigdetection;                     /**< indicates whether to start detection for the original problem */
   SCIP_Bool             enableorigclassification;                /**< indicates whether to start constraint classification for the original problem */

   SCIP_Bool             bendersonlycontsubpr;                    /**< indicates whether only decomposition with only continuous variables in the subproblems should be searched*/
   SCIP_Bool             bendersonlybinmaster;                    /**< indicates whether only decomposition with only binary variables in the master should be searched */
   SCIP_Bool             detectbenders;                           /**< indicates wethher or not benders detection mode is enabled */

   int                   ncallscreatedecomp;                      /**< debugging method for counting the number of calls of created decompositions */

   gcg::DETPROBDATA*		 detprobdatapres;                         /**< detprobdata containing data for the presolved transformed problem */
   gcg::DETPROBDATA*     detprobdataorig;                         /**< detprobdata containing data for the original problem */

   /* score data */
   GCG_SCORE**           scores;                                  /**< array of scores */
   int                   nscores;                                 /**< number of scores */
   char*                 currscore;                               /**< currently chosen score shortname */
   SCIP_CLOCK*           scoreclock;                              /**< clock tracking the total score calculation time */

   PARTIALDECOMP*        partialdectowrite;                       /**< pointer enabling the use of SCIPs writeProb/writeTransProb function for writing partial decompositions*/

   SCIP_Bool             consnamesalreadyrepaired;                /**< stores whether or not */
   std::vector<int>*     userblocknrcandidates;                   /**< vector to store block number candidates that were given by user */
   SCIP_Bool             freeorig;                   /**< help bool to notify a nonfinal free transform (needed if presolving is revoked, e.g. if orig decomposition is used, and transformation is not successful) */
};


/** parameter how to modify scores when comparing decompositions for original and presolved problem
 * (which might differ in size) */
enum weightinggpresolvedoriginaldecomps{
   NO_MODIF = 0,           /**< no modification */
   FRACTION_OF_NNONZEROS,  /**< scores are weighted according to ratio of number nonzeros, the more the worse */
   FRACTION_OF_NROWS,      /**< scores are weighted according to ratio of number rows, the more the worse */
   FAVOUR_PRESOLVED        /**< decompositions for presolved problems are always favoured over decompositions of original problem */
};


/** locally used macro to help with sorting, comparator */
struct sort_pred {
   bool operator()(const std::pair<PARTIALDECOMP*, SCIP_Real> &left, const std::pair<PARTIALDECOMP*, SCIP_Real> &right) {
      return left.second > right.second;
   }
};


//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

/*
 * Local methods
 */


/**
 * @brief local function to get the conshdlr data of the current conshdlr
 *
 * @returns the conshdlrdata iff it exists, else NULL
 * @note returns NULL iff conshdlr or its data does not exist
 */
static
SCIP_CONSHDLRDATA* getConshdlrdata(
   SCIP* scip     /**< SCIP data structure */
   )
{
   /* get current conshdlr, identified by the name define */
   assert(scip != NULL);
   SCIP_CONSHDLR* conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   /* check if conshdlr was found */
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot get its data!\n");
      return NULL;
   }

   /* this will result in NULL if none is stored */
   return SCIPconshdlrGetData(conshdlr);
}


/** @brief local method to handle store a partialdec in the correct detprobdata
 *
 * @returns SCIP status */
static
SCIP_RETCODE addPartialdec(
   SCIP* scip,                   /**< SCIP data structure */
   PARTIALDECOMP*  partialdec    /**< partialdec pointer */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   bool success;

   if( partialdec->isComplete() )
   {
      if( partialdec->isAssignedToOrigProb() )
         success = conshdlrdata->detprobdataorig->addPartialdecToFinished(partialdec);
      else
         success = conshdlrdata->detprobdatapres->addPartialdecToFinished(partialdec);
   }
   else
   {
      if( partialdec->isAssignedToOrigProb() )
         success = conshdlrdata->detprobdataorig->addPartialdecToOpen(partialdec);
      else
         success = conshdlrdata->detprobdatapres->addPartialdecToOpen(partialdec);
   }

   if( !success )
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "Decomposition to add is already known to gcg!\n");

   return SCIP_OKAY;
}


PARTIALDECOMP* GCGconshdlrDecompGetPartialdecFromID(
   SCIP* scip,          /**< SCIP data structure */
   int partialdecid     /**< partialdec id */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   auto itr = conshdlrdata->partialdecsbyid->find(partialdecid);
   if( itr != conshdlrdata->partialdecsbyid->end() )
   {
      return itr->second;
   }

   return NULL;
}


/** @brief translates a vector of PARTIALDECOMP pointers into an array of their ids
 * @returns SCIP return code
 */
static
SCIP_RETCODE partialdecVecToIdArray(
   std::vector<PARTIALDECOMP*>& partialdecs, /**< vector of partialdecs (input) */
   int**          idlist,                    /**< array of ids (output) */
   int*           listlength                 /**< length of id array (output) */
   )
{
   *listlength = (int) partialdecs.size();
   int i;
   for(i = 0; i < (int) partialdecs.size(); i++)
   {
      (*idlist)[i] = partialdecs[i]->getID();
   }

   return SCIP_OKAY;
}


/** @brief gets all selected partialdecs
 * @returns vector of all selected partialdecs */
static
std::vector<PARTIALDECOMP*> getSelectedPartialdecs(
   SCIP* scip,  /**< SCIP data structure */
   std::vector<PARTIALDECOMP*>& selectedpartialdecs  /**< vector of partialdecs (input) */
   )
{
   /* get all partialdecs */
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   /* look for selected ones and add them */
   for(auto partialdec : *conshdlrdata->partialdecs)
   {
      if(partialdec->isSelected())
      {
         selectedpartialdecs.push_back(partialdec);
      }
   }
   return selectedpartialdecs;
}


/** @brief gets vector of all finished partialdecs
 */
static
void getFinishedPartialdecs(
   SCIP*          scip,  /**< SCIP data structure */
   std::vector<PARTIALDECOMP*>& finishedpartialdecs /**< will contain finished partialdecs */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   for(auto partialdec : *conshdlrdata->partialdecs)
   {
      if(partialdec->isComplete())
         finishedpartialdecs.push_back(partialdec);
   }
}


/**
 * method to unselect all decompositions, called in consexit, and when the partialdeclist is updated
 * (especially if new (partial ones) are added )
 *
 *@returns SCIP return code
 */
static
SCIP_RETCODE unselectAllPartialdecs(
   SCIP*          scip     /**< SCIP data structure */
   )
{
   /* get all partialdecs */
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   /* set each partialdec to not selected */
   for(auto partialdec : *conshdlrdata->partialdecs)
   {
      partialdec->setSelected(false);
   }

   return SCIP_OKAY;
}


/** @brief initializes a new detection data structure
 * 
 * @returns new detection data structure
 */
static
PARTIALDEC_DETECTION_DATA* createPartialdecDetectionData(
   gcg::DETPROBDATA *detprobdata,   /**< detprobdata to point to */
   PARTIALDECOMP* partialdec        /**< partialdec to pass to detector */
   )
{
   PARTIALDEC_DETECTION_DATA *partialdecdetdata;
   partialdecdetdata = new PARTIALDEC_DETECTION_DATA();
   partialdecdetdata->detprobdata = detprobdata;
   partialdecdetdata->nnewpartialdecs = 0;
   partialdecdetdata->workonpartialdec = new gcg::PARTIALDECOMP( partialdec );
   return partialdecdetdata;
}


/**
 * @brief resets/creates the detprobdata for the given problem
 * @returns scip return code
 */
static
SCIP_RETCODE resetDetprobdata(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             original           /**< whether to do this for the original problem */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb(scip) );

   // for the orig detprobdata, reset only the current partialdecs as the rest will stay the same in any case
   if(original)
   {
      if( conshdlrdata->detprobdataorig == NULL )
         conshdlrdata->detprobdataorig = new gcg::DETPROBDATA(scip, original);
      conshdlrdata->detprobdataorig->clearCurrentPartialdecs();
   }
   // for the presolved problem, the detprobdata is deleted and a new one is created
   else
   {
      assert(SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED);
      if( conshdlrdata->detprobdatapres != NULL )
         delete conshdlrdata->detprobdatapres;
      conshdlrdata->detprobdatapres = new gcg::DETPROBDATA(scip, original);
   }

   return SCIP_OKAY;
}


/** @brief delets the detection data structure
 * 
 * frees all relevant data within the construct before deleting the main structure
 * @returns SCIP return code
 */
static
SCIP_RETCODE deletePartialdecDetectionData(
   SCIP* scip,                      /**< SCIP data structure */
   PARTIALDEC_DETECTION_DATA* data  /**< data to delete */
   )
{
   SCIPfreeMemoryArrayNull( scip, &(data->newpartialdecs) );
   delete data->workonpartialdec;
   data->newpartialdecs = NULL;
   data->nnewpartialdecs = 0;
   delete data;

   return SCIP_OKAY;
}


/** @brief constructs partialdecs using the registered detectors
 * 
 * Takes the current partialdecs in the detprobdata as root,
 * propagates, finishes and postprocesses in rounds.
 * 
 *  @return user has to free partialdecs */
static
SCIP_Retcode detect(
   SCIP* scip,                      /**< SCIP data structure */
   gcg::DETPROBDATA *detprobdata    /**< detprobdata for problem the detection should be performed on */
   )
{
   SCIP_RESULT result = SCIP_DIDNOTFIND;
   int maxndetectionrounds;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   SCIP_CALL_ABORT( SCIPgetIntParam( scip, "detection/maxrounds", &maxndetectionrounds ) );

   // Fill partialdecs vector into deque
   std::deque<PARTIALDECOMP*> partialdecqueue;
   for( auto partialdecomp: detprobdata->getOpenPartialdecs() )
   {
      partialdecqueue.push_back(partialdecomp);
   }

   // TODO while not strg-c
   while (!partialdecqueue.empty() && (
      conshdlrdata->maxdetectiontime == 0 ||
      SCIPgetClockTime(scip, conshdlrdata->detectorclock) < conshdlrdata->maxdetectiontime
   ))
   {
      gcg::PARTIALDECOMP *partialdec = partialdecqueue.front();
      partialdecqueue.pop_front();

      // Check if max round reached for this partialdec
      if((int) partialdec->getDetectorchain().size() >= maxndetectionrounds)
      {
         continue;
      }

      // TODO (priority for "promising pairs")
      for (int j = 0; j < conshdlrdata->npropagatingdetectors; j++)
      {
         GCG_DETECTOR* detector;
         detector = conshdlrdata->propagatingdetectors[j];

         if( !detector->enabled )
         {
            continue;
         }

         /* skip detector if it should not be recalled */
         if( !detector->usefulRecall && partialdec->isPropagatedBy( detector ) )
            continue;

         PARTIALDEC_DETECTION_DATA* partialdecdetdata = createPartialdecDetectionData(detprobdata, partialdec);

         // PROPAGATE
         SCIP_CALL( detector->propagatePartialdec(scip, detector, partialdecdetdata, &result) );
         detector->dectime += partialdecdetdata->detectiontime;

         // Handle found Partialdecs
         for( int k = 0; k < partialdecdetdata->nnewpartialdecs; ++ k )
         {
            PARTIALDECOMP* newpartialdec = partialdecdetdata->newpartialdecs[k];
            newpartialdec->setDetectorPropagated(detector);
            newpartialdec->prepare();
            newpartialdec->addDecChangesFromAncestor(partialdec);

            // If is already complete => store for POSTPROCESSING
            if( newpartialdec->isComplete() )
            {
               if( !detprobdata->addPartialdecToFinished(newpartialdec) )
               {
                  delete newpartialdec;
                  newpartialdec = NULL;
               }
            }
            else
            {
               // Store for further PROPAGATION
               if( detprobdata->addPartialdecToOpen(newpartialdec) )
                  partialdecqueue.push_back(newpartialdec);
               else
               {
                  delete newpartialdec;
                  newpartialdec = NULL;
               }
            }
            if( newpartialdec != NULL )
               detprobdata->addPartialdecToAncestor(newpartialdec);
         }
         deletePartialdecDetectionData(scip, partialdecdetdata);
      }
      detprobdata->addPartialdecToAncestor(partialdec);
   }

   // FINISH partialdecs
   for( auto partialdecomp: detprobdata->getOpenPartialdecs() )
   {
      for(int l = 0; l < conshdlrdata->nfinishingdetectors; l++)
      {
         GCG_DETECTOR *finishingdetector = conshdlrdata->finishingdetectors[l];
         if( !finishingdetector->enabledFinishing )
         {
            continue;
         }

         PARTIALDEC_DETECTION_DATA *finishingdata = createPartialdecDetectionData(detprobdata, partialdecomp);
         SCIP_CALL( finishingdetector->finishPartialdec(scip, finishingdetector, finishingdata, &result) );
         finishingdetector->dectime += finishingdata->detectiontime;

         for(int finished = 0; finished < finishingdata->nnewpartialdecs; ++finished)
         {
            PARTIALDECOMP* newpartialdec = finishingdata->newpartialdecs[finished];
            newpartialdec->deleteEmptyBlocks(false);
            newpartialdec->setDetectorFinished(finishingdetector);

            newpartialdec->prepare();
            newpartialdec->addDecChangesFromAncestor(partialdecomp);
            if( !detprobdata->addPartialdecToFinished(newpartialdec) )
               delete newpartialdec;
         }
         deletePartialdecDetectionData(scip, finishingdata);
      }
   }

   // POSTPROCESSING of finished partialdecs
   if(conshdlrdata->postprocess)
   {
      SCIP_CLOCK* postprocessingclock;
      SCIPcreateClock( scip, & postprocessingclock );
      SCIP_CALL_ABORT( SCIPstartClock( scip, postprocessingclock ) );
      auto& finishedpartialdecs = detprobdata->getFinishedPartialdecs();
      int numpostprocessed = 0;
      int nfinished = (int) finishedpartialdecs.size();
      for( int i = 0; i < nfinished; ++i )
      {
          auto postpartialdec = finishedpartialdecs[i];

         // Check if postprocessing is enabled globally
         for( int d = 0; d < conshdlrdata->npostprocessingdetectors; ++d )
         {
            GCG_DETECTOR* postdetector = conshdlrdata->postprocessingdetectors[d];
            /* if the postprocessing of the detector is not enabled go on with the next detector */
            if( !postdetector->enabledPostprocessing )
            {
               continue;
            }

            PARTIALDEC_DETECTION_DATA* partialdecdetdata = createPartialdecDetectionData(detprobdata, postpartialdec);

            SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "call finisher for detector %s \n", GCGdetectorGetName( postdetector ) );

            // POSTPROCESS
            SCIP_CALL( postdetector->postprocessPartialdec( scip, postdetector, partialdecdetdata, &result ) );
            postdetector->dectime += partialdecdetdata->detectiontime;

            for( int finished = 0; finished < partialdecdetdata->nnewpartialdecs; ++finished )
            {
               PARTIALDECOMP* newpartialdec = partialdecdetdata->newpartialdecs[finished];

               newpartialdec->setDetectorPropagated(postdetector);
               newpartialdec->setFinishedByFinisher(true );
               newpartialdec->prepare();
               newpartialdec->addDecChangesFromAncestor(postpartialdec);

               if( !detprobdata->addPartialdecToFinished(newpartialdec) )
                  delete newpartialdec;
               else
                  numpostprocessed += 1;

            }
            deletePartialdecDetectionData(scip, partialdecdetdata);
         }
      }

      SCIP_CALL_ABORT( SCIPstopClock( scip, postprocessingclock ) );

      detprobdata->postprocessingtime += SCIPgetClockTime( scip, postprocessingclock );
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "POSTPROCESSING of decompositions. Added %d new decomps. \n", numpostprocessed);
      SCIP_CALL_ABORT( SCIPfreeClock( scip, & postprocessingclock ) );
   }
   else
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "POSTPROCESSING disabled\n");
   } // end for postprocessing finished partialdecs

   // STATS
   // count the successful refinement calls for each detector

   // finished partialdecs
   for( PARTIALDECOMP* partialdec: detprobdata->getFinishedPartialdecs() )
   {
      assert( partialdec->checkConsistency() );
      assert( partialdec->getNOpenconss() == 0 );
      assert( partialdec->getNOpenvars() == 0 );

      for( int d = 0; d < conshdlrdata->ndetectors; ++ d )
      {
         if( partialdec->isPropagatedBy(conshdlrdata->detectors[d]) )
         {
            conshdlrdata->detectors[d]->ndecomps += 1;
            conshdlrdata->detectors[d]->ncompletedecomps += 1;
         }
      }
   }

   // open partialdecs
   for( PARTIALDECOMP* partialdec: detprobdata->getOpenPartialdecs() )
   {
      assert( partialdec->checkConsistency() );

      for( int d = 0; d < conshdlrdata->ndetectors; ++ d )
      {
         if( partialdec->isPropagatedBy(conshdlrdata->detectors[d]) )
            conshdlrdata->detectors[d]->ndecomps += 1;
      }
   }

   /* preliminary output detector stats */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Found %d finished decompositions.\n",
      detprobdata->getNFinishedPartialdecs());
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Measured running time per detector:\n");

   for( int i = 0; i < conshdlrdata->ndetectors; ++ i )
   {
      if( conshdlrdata->detectors[i]->ncompletedecomps > 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL,
            "Detector %-25.25s worked on %8d finished decompositions and took a total time of %10.3f\n",
            GCGdetectorGetName( conshdlrdata->detectors[i]), conshdlrdata->detectors[i]->ncompletedecomps,
            conshdlrdata->detectors[i]->dectime);
      }
   }

   return SCIP_OKAY;
}


/** @brief initialization method of constraint handler (called after problem was transformed) */
static
SCIP_DECL_CONSINIT(consInitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* the detection hast not run yet */
   conshdlrdata->hasrun = FALSE;
   conshdlrdata->hasrunoriginal = FALSE;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      GCG_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      /* add the detector time of each detector to 0 */
      detector->dectime = 0.;
      /* init the detectors */
      if( detector->initDetector != NULL )
      {
         SCIPdebugMessage("Calling initDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->initDetector)(scip, detector) );
      }
   }

   return SCIP_OKAY;
}


/** @brief deinitialization method of constraint handler (called before transformed problem is freed) */
static
SCIP_DECL_CONSEXIT(consExitDecomp)
{ /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(conshdlr != NULL);
   assert(scip != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* remove all decomps */
   if( conshdlrdata->ndecomps > 0 && conshdlrdata->decomps != NULL )
   {
      for( int dec = 0; dec < conshdlrdata->ndecomps; ++dec )
      {
         GCGdecompFree(scip, &conshdlrdata->decomps[conshdlrdata->ndecomps - dec - 1]);
      }

      /* remove decomp array structure */
      SCIPfreeBlockMemoryArray(scip, &conshdlrdata->decomps, conshdlrdata->ndecomps);
      conshdlrdata->ndecomps = 0;
      conshdlrdata->decomps = NULL;
   }

   /* reset the run */
   conshdlrdata->hasrun = FALSE;

   /* release the detectors' data sets */
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      GCG_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      if( detector->exitDetector != NULL )
      {
         SCIPdebugMessage("Calling exitDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->exitDetector)(scip, detector) );
      }
   }

   GCGconshdlrDecompFreeDetprobdata(scip);

   /* remove selection of partialdecs */
   unselectAllPartialdecs(scip);

   return SCIP_OKAY;
}


/** @brief destructor of constraint handler to free constraint handler data (called when SCIP is exiting) */
static
SCIP_DECL_CONSFREE(consFreeDecomp)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /* free detection time clocks */
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPfreeClock(scip, &conshdlrdata->completedetectionclock) );

   /* free all detectors */
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      GCG_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      if( detector->freeDetector != NULL )
      {
         SCIPdebugMessage("Calling freeDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->freeDetector)(scip, detector) );
      }

      BMSfreeMemoryArray(&detector->name);
      BMSfreeMemoryArray(&detector->description);

      SCIPfreeBlockMemory(scip, &detector);
   }

   /* free all consclassifiers */
   for( i = 0; i < conshdlrdata->nconsclassifiers; ++i )
   {
      GCG_CONSCLASSIFIER *consclassifier = conshdlrdata->consclassifiers[i];
      assert(consclassifier != NULL);

      if( consclassifier->freeClassifier != NULL )
      {
         SCIPdebugMessage("Calling freeClassifier of consclassifier %s\n", consclassifier->name);
         SCIP_CALL( (*consclassifier->freeClassifier)(scip, consclassifier) );
      }

      BMSfreeMemoryArray(&consclassifier->name);
      BMSfreeMemoryArray(&consclassifier->description);

      SCIPfreeBlockMemory(scip, &conshdlrdata->consclassifiers[i]);
   }

   /* free all varclassifiers */
   for( i = 0; i < conshdlrdata->nvarclassifiers; ++i )
   {
      GCG_VARCLASSIFIER *varclassifier = conshdlrdata->varclassifiers[i];
      assert(varclassifier != NULL);

      if( varclassifier->freeClassifier != NULL )
      {
         SCIPdebugMessage("Calling freeClassifier of varclassifier %s\n", varclassifier->name);
         SCIP_CALL( (*varclassifier->freeClassifier)(scip, varclassifier) );
      }

      BMSfreeMemoryArray(&varclassifier->name);
      BMSfreeMemoryArray(&varclassifier->description);

      SCIPfreeBlockMemory(scip, &conshdlrdata->varclassifiers[i]);
   }

   /* free all scores */
   for( i = 0; i < conshdlrdata->nscores; ++i )
   {
      GCG_SCORE* score = conshdlrdata->scores[i];
      assert(score != NULL);

      if( score->scorefree != NULL )
      {
         SCIPdebugMessage("Calling freeScore of score %s\n", score->name);
         SCIP_CALL( (*score->scorefree)(scip, score) );
      }

      BMSfreeMemoryArray(&score->name);
      BMSfreeMemoryArray(&score->shortname);
      BMSfreeMemoryArray(&score->description);

      SCIPfreeBlockMemory(scip, &conshdlrdata->scores[i]);
   }
   SCIPfreeClock(scip, &conshdlrdata->scoreclock);

   /* remove all remaining data */
   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeMemoryArray(scip, &conshdlrdata->propagatingdetectors);
   SCIPfreeMemoryArray(scip, &conshdlrdata->finishingdetectors);
   SCIPfreeMemoryArray(scip, &conshdlrdata->postprocessingdetectors);
   SCIPfreeMemoryArray(scip, &conshdlrdata->consclassifiers);
   SCIPfreeMemoryArray(scip, &conshdlrdata->consclassifierpriorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->varclassifiers);
   SCIPfreeMemoryArray(scip, &conshdlrdata->varclassifierpriorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->scores);

   delete conshdlrdata->userblocknrcandidates;
   delete conshdlrdata->partialdecs;
   delete conshdlrdata->partialdecsbyid;

   /* remove the data structure of the conshdlr */
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


/**
 * @brief finds a non duplicate constraint name of the form c_{a} with minimal natural number {a}
 * @return id of constraint
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
      /* create new name candidate */
      char candidatename[SCIP_MAXSTRLEN] = "c_";
      char number[20];
      snprintf(number, sizeof(number), "%d", candidatenumber );
      strcat(candidatename, number );

      /* check candidate, if it is not free increase counter for candidate number */
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


/** @brief creates a partialdec for a given decomposition
 *
 * @returns SCIP return code
 * */
static
SCIP_RETCODE createPartialdecFromDecomp(
   SCIP* scip,                   /**< SCIP data structure */
   GCG_DECOMP* decomp,           /**< decomposition the partialdec is created for */
   PARTIALDECOMP** newpartialdec /**< the new partialdec created from the decomp */
   )
{
   assert( decomp != NULL );
   assert( GCGdecompCheckConsistency( scip, decomp ) );

   /* get detprobdata for var & cons indices */
   DETPROBDATA* detprobdata = decomp->presolved ? GCGconshdlrDecompGetDetprobdataPresolved(scip) : GCGconshdlrDecompGetDetprobdataOrig(scip);

   /* create new partialdec and initialize its data */
   PARTIALDECOMP* partialdec = new PARTIALDECOMP(scip, !decomp->presolved);
   partialdec->setNBlocks( GCGdecompGetNBlocks(decomp));

   SCIP_CONS** linkingconss = GCGdecompGetLinkingconss(decomp);
   int nlinkingconss = GCGdecompGetNLinkingconss(decomp);
   SCIP_HASHMAP* constoblock = GCGdecompGetConstoblock(decomp);
   int nblock;

   /* set linking conss */
   for( int c = 0; c < nlinkingconss; ++c )
   {
      partialdec->fixConsToMaster(detprobdata->getIndexForCons(linkingconss[c]));
   }

   /* set block conss */
   for( int c = 0; c < detprobdata->getNConss(); ++c )
   {
      nblock = (int) (size_t) SCIPhashmapGetImage(constoblock, (void*) (size_t) detprobdata->getCons(c));
      if( nblock >= 1 && nblock <= partialdec->getNBlocks() )
      {
         partialdec->fixConsToBlock(c, nblock - 1);
      }
   }

   SCIP_VAR*** stairlinkingvars = GCGdecompGetStairlinkingvars(decomp);

   if( stairlinkingvars != NULL )
   {
      int* nstairlinkingvars = GCGdecompGetNStairlinkingvars(decomp);
      int varindex;

      /* set stairlinkingvars */
      for( int b = 0; b < partialdec->getNBlocks(); ++b )
      {
         for( int v = 0; v < nstairlinkingvars[b]; ++v )
         {
            if( stairlinkingvars[b][v] != NULL )
            {
               varindex = detprobdata->getIndexForVar(stairlinkingvars[b][v]);
               partialdec->fixVarToStairlinking(varindex, b);
            }
         }
      }
   }

   /* set other vars */
   SCIP_HASHMAP* vartoblock = GCGdecompGetVartoblock(decomp);
   if( vartoblock != NULL )
   {
      for( int v = 0; v < detprobdata->getNVars(); ++v )
      {
         nblock = (int) (size_t) SCIPhashmapGetImage(vartoblock,
            (void*) (size_t) SCIPvarGetProbvar(detprobdata->getVar(v)));
         if( nblock == partialdec->getNBlocks() + 2 && !partialdec->isVarStairlinkingvar(v))
         {
            partialdec->fixVarToLinking(v);
         }
         else if( nblock == partialdec->getNBlocks() + 1 )
         {
            partialdec->fixVarToMaster(v);
         }
         else if( nblock >= 1 && nblock <= partialdec->getNBlocks())
         {
            partialdec->fixVarToBlock(v, nblock - 1);
         }
      }
   }

   partialdec->sort();

   /* now all conss and vars should be assigned */
   assert( partialdec->isComplete() );

   /*set all detector-related information*/
   for( int i = 0; i < GCGdecompGetDetectorChainSize(decomp); ++i )
   {
      partialdec->setDetectorPropagated(GCGdecompGetDetectorChain(decomp)[i]);
      partialdec->addClockTime(GCGdecompGetDetectorClockTimes(decomp)[i]);
      partialdec->addPctConssFromFree(1 - *(GCGdecompGetDetectorPctConssFromOpen(decomp)));
      partialdec->addPctConssToBlock(*(GCGdecompGetDetectorPctConssToBlock(decomp)));
      partialdec->addPctConssToBorder(*(GCGdecompGetDetectorPctConssToBorder(decomp)));
      partialdec->addPctVarsFromFree(1 - *(GCGdecompGetDetectorPctVarsFromOpen(decomp)));
      partialdec->addPctVarsToBlock(*(GCGdecompGetDetectorPctVarsToBlock(decomp)));
      partialdec->addPctVarsToBorder(*(GCGdecompGetDetectorPctVarsToBorder(decomp)));
      partialdec->addNNewBlocks(*(GCGdecompGetNNewBlocks(decomp)));
   }

   /* calc maxwhitescore and hashvalue */
   partialdec->prepare();
   partialdec->calcStairlinkingVars();

   *newpartialdec = partialdec;
   return SCIP_OKAY;
}


/**
 * @brief creates a decomposition GCG_DECOMP structure for a given partialdec
 * @returns scip return code
 */
static
SCIP_RETCODE createDecompFromPartialdec(
   SCIP* scip,                /**< SCIP data structure */
   PARTIALDECOMP* partialdec, /**< partialdec the decomposition is created for */
   GCG_DECOMP** newdecomp     /**< the new decomp created from the partialdec */
   )
{
   SCIP_HASHMAP* vartoblock = NULL;
   SCIP_HASHMAP* constoblock = NULL;
   SCIP_HASHMAP* varindex = NULL;
   SCIP_HASHMAP* consindex = NULL;
   SCIP_VAR*** stairlinkingvars = NULL;
   SCIP_VAR*** subscipvars = NULL;
   SCIP_VAR** linkingvars = NULL;
   SCIP_CONS** linkingconss = NULL;
   SCIP_CONS*** subscipconss = NULL;

   int* nsubscipconss = NULL;
   int* nsubscipvars = NULL;
   int* nstairlinkingvars = NULL;
   int nlinkingvars;

   int varcounter = 1; /* in varindex counting starts with 1 */
   int conscounter = 1; /* in consindex counting starts with 1 */
   int counterstairlinkingvars = 0;
   int modifier;
   int nlinkingconss;
   int ndeletedblocks = 0;
   int nmastervarsfromdeleted = 0;
   int v;
   int c;

   assert( partialdec->checkConsistency() );

   DETPROBDATA* detprobdata = partialdec->getDetprobdata();
   assert(detprobdata != NULL);

   std::vector<SCIP_Bool> isblockdeleted = std::vector<SCIP_Bool>(partialdec->getNBlocks(), FALSE);
   std::vector<int> ndeletedblocksbefore = std::vector<int>(partialdec->getNBlocks(), 0);
   std::vector<int> mastervaridsfromdeleted = std::vector<int>(0);
   std::vector<SCIP_Var*> mastervarsfromdeleted = std::vector<SCIP_VAR*>(0);
   std::vector<SCIP_CONS*> relevantconss = detprobdata->getRelevantConss();
   std::vector<SCIP_VAR*> relevantvars = detprobdata->getRelevantVars();
   std::vector<SCIP_VAR*> origfixedtozerovars = detprobdata->getOrigVarsFixedZero();

   /* create decomp data structure */
   SCIP_CALL_ABORT( GCGdecompCreate( scip, newdecomp ) );

   GCGdecompSetPresolved(*newdecomp, !partialdec->isAssignedToOrigProb());

   /* find out if for some blocks all conss have been deleted */
   for( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      SCIP_Bool iscurblockdeleted = TRUE;
      for( c = 0; c < partialdec->getNConssForBlock( b ); ++c )
      {
         int consid = partialdec->getConssForBlock( b )[c];
         SCIP_CONS* scipcons = relevantconss[consid];

         if( scipcons != NULL && !SCIPconsIsDeleted( scipcons)  )
         {
            iscurblockdeleted = FALSE;
            break;
         }
      }
      if ( iscurblockdeleted )
      {
         ++ndeletedblocks;
         isblockdeleted[b] = TRUE;
         for( int b2 = b+1; b2 < partialdec->getNBlocks(); ++b2)
         {
            ++ndeletedblocksbefore[b2];
         }
         /* store deletion information of included vars */
         for( v = 0; v < partialdec->getNVarsForBlock( b ); ++v )
         {
            int varid = partialdec->getVarsForBlock( b )[v];
            SCIP_VAR* scipvar = relevantvars[varid];
            mastervaridsfromdeleted.push_back(varid);
            mastervarsfromdeleted.push_back(scipvar);
            ++nmastervarsfromdeleted;
         }
      }
   }

   /* set nblocks */
   GCGdecompSetNBlocks( *newdecomp, partialdec->getNBlocks() - ndeletedblocks );

   if( partialdec->getNBlocks() - ndeletedblocks == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "All blocks have been deleted since only deleted constraints are contained, no reformulation is done.\n");
   }

   /* prepare constraints data structures */
   if( partialdec->getNMasterconss() != 0 )
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &linkingconss, partialdec->getNMasterconss()) );
   else
      linkingconss = NULL;

   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nsubscipconss, partialdec->getNBlocks() - ndeletedblocks) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipconss, partialdec->getNBlocks() - ndeletedblocks) );

   SCIP_CALL_ABORT( SCIPhashmapCreate(&constoblock, SCIPblkmem(scip), partialdec->getNConss()) );
   SCIP_CALL_ABORT( SCIPhashmapCreate(&consindex, SCIPblkmem(scip), partialdec->getNConss()) );

   /* set linking constraints */
   modifier = 0;
   nlinkingconss = partialdec->getNMasterconss();
   for( c = 0; c < partialdec->getNMasterconss(); ++c )
   {
      int consid = partialdec->getMasterconss()[c];
      SCIP_CONS* scipcons = relevantconss[consid];
      if( partialdec->isAssignedToOrigProb() )
         SCIPgetTransformedCons(scip, scipcons, &scipcons);

      if( scipcons == NULL || SCIPconsIsDeleted(scipcons) || SCIPconsIsObsolete(scipcons) )
      {
         --nlinkingconss;
         ++modifier;
      }
      else
      {
         linkingconss[c-modifier] = scipcons;
         SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, scipcons, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
         SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipcons, (void*) (size_t) conscounter ) );
         conscounter ++;
      }
   }

   if( nlinkingconss != 0 )
      GCGdecompSetLinkingconss(scip, *newdecomp, linkingconss, nlinkingconss);
   else
      linkingconss = NULL;

   /* set block constraints */
   for( int b = 0; b < partialdec->getNBlocks(); ++ b )
   {
      if( isblockdeleted[b] )
         continue;
      modifier = 0;
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipconss[b-ndeletedblocksbefore[b]], partialdec->getNConssForBlock(b)) );
      nsubscipconss[b-ndeletedblocksbefore[b]] = partialdec->getNConssForBlock(b);
      for( c = 0; c < partialdec->getNConssForBlock(b); ++c )
      {
         int consid = partialdec->getConssForBlock(b)[c];
         SCIP_CONS* scipcons = relevantconss[consid];
         if( partialdec->isAssignedToOrigProb() )
            SCIPgetTransformedCons(scip, scipcons, &scipcons);

         if( scipcons == NULL || SCIPconsIsDeleted(scipcons) )
         {
            --nsubscipconss[b-ndeletedblocksbefore[b]];
            ++modifier;
         }
         else
         {
            assert( scipcons != NULL );
            subscipconss[b-ndeletedblocksbefore[b]][c-modifier] = scipcons;
            SCIPdebugMessage("Set cons %s to block %d + 1 - %d in cons to block\n", SCIPconsGetName(scipcons), b, ndeletedblocksbefore[b] );
            SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, scipcons, (void*) (size_t)(b + 1 - ndeletedblocksbefore[b])) );
            SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipcons, (void*) (size_t) conscounter) );
            conscounter ++;
         }
      }
   }

   /* assign all open conss that might be left */
   for( c = 0; c < SCIPgetNConss(scip); ++c )
   {
      if( !GCGisConsGCGCons(SCIPgetConss(scip)[c]) )
      {
         if(!SCIPhashmapExists(constoblock, SCIPgetConss(scip)[c]))
         {
            SCIP_CONS* scipcons = SCIPgetConss(scip)[c];
            SCIP_CALL_ABORT( SCIPhashmapInsert(constoblock, scipcons, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
            SCIP_CALL_ABORT( SCIPhashmapInsert(consindex, scipcons, (void*) (size_t) conscounter) );
            conscounter ++;
         }
      }
   }

   GCGdecompSetSubscipconss(scip, *newdecomp, subscipconss, nsubscipconss);
   GCGdecompSetConstoblock(*newdecomp, constoblock);
   GCGdecompSetConsindex(*newdecomp, consindex);
   /* finished setting constraint data structures */

   /* prepare constraints data structures */
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nsubscipvars, partialdec->getNBlocks() - ndeletedblocks) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &subscipvars, partialdec->getNBlocks() - ndeletedblocks) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &nstairlinkingvars, partialdec->getNBlocks() - ndeletedblocks) );
   SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &stairlinkingvars, partialdec->getNBlocks() -ndeletedblocks) );

   SCIP_CALL_ABORT( SCIPhashmapCreate(&vartoblock, SCIPblkmem(scip), partialdec->getNVars() + (int) origfixedtozerovars.size()) );
   SCIP_CALL_ABORT( SCIPhashmapCreate(&varindex, SCIPblkmem(scip), partialdec->getNVars() + (int) origfixedtozerovars.size()) );

   /* set linkingvars */
   nlinkingvars = partialdec->getNLinkingvars() + partialdec->getNMastervars() + partialdec->getNTotalStairlinkingvars() + nmastervarsfromdeleted + (int) origfixedtozerovars.size();

   if( nlinkingvars != 0 )
      SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &linkingvars, nlinkingvars) );
   else
      linkingvars = NULL;

   for( v = 0; v < partialdec->getNLinkingvars(); ++ v )
   {
      int var = partialdec->getLinkingvars()[v];
      SCIP_VAR* scipvar = SCIPvarGetProbvar(relevantvars[var]);
      assert( scipvar != NULL );

      linkingvars[v] = scipvar;
      SCIPdebugMessage( "Set var %s to block %d + 2 - %d in var to block\n", SCIPvarGetName(scipvar), partialdec->getNBlocks(), ndeletedblocks );
      SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t)(partialdec->getNBlocks() + 2 - ndeletedblocks)) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
      varcounter ++;
   }

   for( v = 0; v < partialdec->getNMastervars(); ++ v )
   {
      int var = partialdec->getMastervars()[v];
      SCIP_VAR* scipvar = SCIPvarGetProbvar( relevantvars[var] );
      linkingvars[v + partialdec->getNLinkingvars()] = scipvar;
      SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
      varcounter ++;
   }

   for( v = 0; v < nmastervarsfromdeleted; ++v)
   {
      SCIP_VAR* var;
      var = SCIPvarGetProbvar(mastervarsfromdeleted[v]);

      linkingvars[partialdec->getNMastervars() + partialdec->getNLinkingvars() + v] = var;
      SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, var, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, var, (void*) (size_t) varcounter) );
      varcounter ++;
   }

   for( v = 0; v < (int) origfixedtozerovars.size(); ++v)
   {
      SCIP_VAR* var;
      var = origfixedtozerovars[v];

      linkingvars[partialdec->getNMastervars() + partialdec->getNLinkingvars() + nmastervarsfromdeleted + v] = var;
      SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, var, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
      SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, var, (void*) (size_t) varcounter) );
      varcounter ++;
   }

   /* set block variables */
   for( int b = 0; b < partialdec->getNBlocks(); ++ b )
   {
      if( isblockdeleted[b] )
         continue;

      if( partialdec->getNVarsForBlock( b ) > 0 )
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, & subscipvars[b -ndeletedblocksbefore[b]], partialdec->getNVarsForBlock(b)) );
      else
         subscipvars[b-ndeletedblocksbefore[b]] = NULL;

      if( partialdec->getNStairlinkingvars( b ) > 0 )
         SCIP_CALL_ABORT( SCIPallocBufferArray(scip, & stairlinkingvars[b-ndeletedblocksbefore[b]], partialdec->getNStairlinkingvars(b)) );
      else
         stairlinkingvars[b-ndeletedblocksbefore[b]] = NULL;

      nsubscipvars[b-ndeletedblocksbefore[b]] = partialdec->getNVarsForBlock(b);
      nstairlinkingvars[b-ndeletedblocksbefore[b]] = partialdec->getNStairlinkingvars(b);

      for( v = 0; v < partialdec->getNVarsForBlock(b); ++ v )
      {
         int var = partialdec->getVarsForBlock(b)[v];
         SCIP_VAR* scipvar = SCIPvarGetProbvar(relevantvars[var]);
         assert( scipvar != NULL );

         subscipvars[b-ndeletedblocksbefore[b]][v] = scipvar;
         SCIPdebugMessage("Set var %s to block %d + 1 - %d in var to block\n", SCIPvarGetName(scipvar), b, ndeletedblocksbefore[b] );
         assert( !SCIPhashmapExists(vartoblock, scipvar) || SCIPhashmapGetImage(vartoblock, scipvar) == (void*) (size_t)(b + 1 - ndeletedblocksbefore[b]) );
         SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t)(b + 1 - ndeletedblocksbefore[b])) );
         SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
         varcounter ++;
      }

      for( v = 0; v < partialdec->getNStairlinkingvars(b); ++ v )
      {
         int var = partialdec->getStairlinkingvars(b)[v];
         SCIP_VAR* scipvar = SCIPvarGetProbvar(relevantvars[var]);
         assert( scipvar != NULL );

         stairlinkingvars[b-ndeletedblocksbefore[b]][v] = scipvar;
         linkingvars[partialdec->getNLinkingvars() + partialdec->getNMastervars() + nmastervarsfromdeleted + counterstairlinkingvars] = scipvar;
         SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t)(partialdec->getNBlocks() + 2 - ndeletedblocks)) );
         SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
         varcounter ++;
         counterstairlinkingvars ++;
      }
   }

   /* if any var is open yet put it in master (in case they were not in relevantvars) */
   for( v = 0; v < SCIPgetNVars(scip); ++v )
   {
      if(!SCIPhashmapExists(vartoblock, SCIPgetVars(scip)[v]))
      {
         SCIP_VAR* scipvar = SCIPvarGetProbvar( SCIPgetVars(scip)[v] );
         SCIP_CALL_ABORT( SCIPhashmapInsert(vartoblock, scipvar, (void*) (size_t)(partialdec->getNBlocks() + 1 - ndeletedblocks)) );
         SCIP_CALL_ABORT( SCIPhashmapInsert(varindex, scipvar, (void*) (size_t) varcounter) );
         varcounter ++;
      }
   }
   
   GCGdecompSetSubscipvars(scip, *newdecomp, subscipvars, nsubscipvars);
   GCGdecompSetStairlinkingvars(scip, *newdecomp, stairlinkingvars, nstairlinkingvars);
   GCGdecompSetLinkingvars(scip, *newdecomp, linkingvars, nlinkingvars, (int) origfixedtozerovars.size(), partialdec->getNMastervars() + nmastervarsfromdeleted);
   GCGdecompSetVarindex(*newdecomp, varindex);
   GCGdecompSetVartoblock(*newdecomp, vartoblock);

   /* //////////////////// free stuff ////////////////////////////// */

   /* free vars stuff */
   SCIPfreeBufferArrayNull(scip, &linkingvars);
   for( int b = partialdec->getNBlocks() - 1 - ndeletedblocks; b >= 0; --b )
   {
      if( nstairlinkingvars[b] != 0 )
      {
         SCIPfreeBufferArrayNull(scip, &stairlinkingvars[b]);
      }
   }

   SCIPfreeBufferArrayNull( scip, &stairlinkingvars);
   SCIPfreeBufferArrayNull( scip, &nstairlinkingvars);

   for( int b = partialdec->getNBlocks() - 1 - ndeletedblocks; b >= 0; --b )
   {
      if( nsubscipvars[b] != 0 )
      {
         SCIPfreeBufferArrayNull(scip, &subscipvars[b]);
      }
   }

   SCIPfreeBufferArrayNull( scip, &subscipvars);
   SCIPfreeBufferArrayNull( scip, &nsubscipvars);

   /* free constraints */
   for( int b = partialdec->getNBlocks() - 1 - ndeletedblocks; b >= 0; --b )
   {
      SCIPfreeBufferArrayNull( scip, &subscipconss[b]);
   }
   SCIPfreeBufferArrayNull( scip, &subscipconss);

   SCIPfreeBufferArrayNull( scip, &nsubscipconss);
   SCIPfreeBufferArrayNull( scip, &linkingconss);

   /* set detectorchain */
   GCGdecompSetDetectorChain(scip, (*newdecomp), partialdec->getDetectorchain().data(), (int) partialdec->getDetectorchain().size());

   /* set last detector in chain as detector that "found" this decomposition */
   if(partialdec->getNDetectors() > 0)
      GCGdecompSetDetector(*newdecomp, partialdec->getDetectorchain().back());

   /* set statistical detector chain data */
   GCGdecompSetPartialdecID(*newdecomp, partialdec->getID());
   if( partialdec->getNDetectors() > 0 )
   {
      GCGdecompSetDetectorClockTimes(scip, *newdecomp, partialdec->getDetectorClockTimes().data());
      GCGdecompSetDetectorPctVarsToBorder(scip, *newdecomp, partialdec->getPctVarsToBorderVector().data());
      GCGdecompSetDetectorPctVarsToBlock(scip, *newdecomp, partialdec->getPctVarsToBlockVector().data());
      GCGdecompSetDetectorPctVarsFromOpen(scip, *newdecomp, partialdec->getPctVarsFromFreeVector().data());
      GCGdecompSetDetectorPctConssToBorder(scip, *newdecomp, partialdec->getPctConssToBorderVector().data());
      GCGdecompSetDetectorPctConssToBlock(scip, *newdecomp, partialdec->getPctConssToBlockVector().data());
      GCGdecompSetDetectorPctConssFromOpen(scip, *newdecomp, partialdec->getPctConssFromFreeVector().data());
      GCGdecompSetNNewBlocks(scip, *newdecomp, partialdec->getNNewBlocksVector().data());
   }

   /* set dectype */
   int newnlinkingvars = GCGdecompGetNLinkingvars((*newdecomp));
   int newnlinkingconss = GCGdecompGetNLinkingconss((*newdecomp));

   if( newnlinkingvars == partialdec->getNTotalStairlinkingvars() && newnlinkingconss == 0 && newnlinkingvars > 0 )
   {
      GCGdecompSetType((*newdecomp), GCG_DECTYPE_STAIRCASE);
   }
   else if( newnlinkingvars > 0 || partialdec->getNTotalStairlinkingvars() > 0 )
   {
      GCGdecompSetType((*newdecomp), GCG_DECTYPE_ARROWHEAD);
   }
   else if( newnlinkingconss > 0 )
   {
      GCGdecompSetType((*newdecomp), GCG_DECTYPE_BORDERED);
   }
   else if( newnlinkingconss == 0 && partialdec->getNTotalStairlinkingvars() == 0 )
   {
      GCGdecompSetType((*newdecomp), GCG_DECTYPE_DIAGONAL);
   }
   else
   {
      GCGdecompSetType((*newdecomp), GCG_DECTYPE_UNKNOWN);
   }

   /* set max white score */
   SCIPdebugMessage(" partialdec maxwhitescore: %f\n", partialdec->getMaxWhiteScore());

   GCGdecompSetMaxWhiteScore(scip, *newdecomp, partialdec->getMaxWhiteScore() );

   /* set detector string */
   char buffer[SCIP_MAXSTRLEN];
   partialdec->buildDecChainString(buffer);
   GCGdecompSetDetectorChainString(scip, *newdecomp, buffer);

   if( !partialdec->isAssignedToOrigProb() )
      SCIP_CALL(GCGdecompAddRemainingConss(scip, *newdecomp) );

   assert(GCGdecompCheckConsistency(scip, *newdecomp) );

   return SCIP_OKAY;
}


/** @brief sorts all registered partialdecs according to score, descending */
static
void sortPartialdecs(
   SCIP* scip  /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   std::sort(conshdlrdata->partialdecs->begin(), conshdlrdata->partialdecs->end(), [&](PARTIALDECOMP* a, PARTIALDECOMP* b) {return (a->getScore(GCGgetCurrentScore(scip)) > b->getScore(GCGgetCurrentScore(scip))); });
}


/** @brief method to adapt score for orig decomps
 * @todo: change score for some parameter settings
 *
 * @returns new score */
static
SCIP_Real GCGconshdlrDecompAdaptScore(
   SCIP*             scip,    /**< SCIP data structure */
   SCIP_Real         oldscore /**< current score (to be updated) */
   )
{
   SCIP_Real score = oldscore;
   int method;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   SCIP_CALL(SCIPgetIntParam(scip, "detection/origprob/advanced/weightinggpresolvedoriginaldecomps", &method) );

   if( method == FRACTION_OF_NNONZEROS )
   {
      if ( conshdlrdata->detprobdatapres == NULL || conshdlrdata->detprobdataorig == NULL )
         return score;

      score *= (SCIP_Real) conshdlrdata->detprobdataorig->getNNonzeros() / conshdlrdata->detprobdatapres->getNNonzeros();
   }

   if( method == FRACTION_OF_NROWS )
   {
      if ( conshdlrdata->detprobdatapres == NULL || conshdlrdata->detprobdataorig == NULL )
         return score;

      score *= (SCIP_Real) conshdlrdata->detprobdataorig->getNConss() / conshdlrdata->detprobdatapres->getNConss();
   }

   if( method == FAVOUR_PRESOLVED )
   {
      score += 1.;
   }

   return score;
}


/** @brief adds constraint partitions with a reduced number of classes */
static
void reduceConsclasses(
   SCIP*                scip,               /**< SCIP data structure */
   gcg::DETPROBDATA*    detprobdata         /**< classification is for problem to which these data correspond */
   )
{
   /* set the number of classes the partitions should be reduced to */
   int maxnclasses = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   if( detprobdata->getNConss() + detprobdata->getNVars() >= 50000 )
      maxnclasses = conshdlrdata->maxnclassesperpartitionforlargeprobs;
   else
      maxnclasses = conshdlrdata->maxnclassesperpartition;

   for( int partitionid = 0; partitionid < detprobdata->getNConsPartitions(); ++ partitionid )
   {
      ConsPartition* newpartition = detprobdata->getConsPartition(partitionid)->reduceClasses(maxnclasses);

      if( newpartition != NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Added reduced version of conspartition %s with %d  different constraint classes \n",
                         detprobdata->getConsPartition(partitionid)->getName(), maxnclasses);
         detprobdata->addConsPartition(newpartition);
      }
   }
}


/** @brief adds variable partitions with a reduced number of classes */
static
void reduceVarclasses(
   SCIP*                scip,                /**< SCIP data structure */
   gcg::DETPROBDATA *detprobdata             /**< classification is for problem to which these data correspond */
   )
{
   /* set the number of classes the partitions should be reduced to */
   int maxnclasses = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   if( detprobdata->getNConss() + detprobdata->getNVars() >= 50000 )
      maxnclasses = conshdlrdata->maxnclassesperpartitionforlargeprobs;
   else
      maxnclasses = conshdlrdata->maxnclassesperpartition;

   for( int partitionid = 0; partitionid < detprobdata->getNVarPartitions(); ++ partitionid )
   {
      VarPartition* newpartition = detprobdata->getVarPartition(partitionid)->reduceClasses(maxnclasses);

      if( newpartition != NULL )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " Added reduced version of varpartition %s with %d different variable classes\n",
                         detprobdata->getVarPartition(partitionid)->getName(), maxnclasses );
         detprobdata->addVarPartition(newpartition);
      }
   }
}


/** @brief sets detection/enabled setting
 * @returns SCIP return code */
static
SCIP_RETCODE setDetectionEnabled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             quiet,              /**< should the parameter be set quietly (no output) */
   SCIP_Bool             enabled             /**< should the detection be enabled */
   )
{
   SCIP_CALL( SCIPsetBoolParam(scip, "detection/enabled", enabled) );
   if( !quiet )
   {
      SCIPinfoMessage(scip, NULL, "detection/enabled = %s\n", enabled ? "TRUE" : "FALSE");
   }
   return SCIP_OKAY;
}


/** @brief resets the parameters to their default value
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

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER )
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", FALSE) );

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
         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
      }
   }

   setDetectionEnabled(scip, quiet, TRUE);

   return SCIP_OKAY;
}


/** @brief sets the parameters to aggressive values
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

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER)
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", FALSE) );

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

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
      }
   }

   setDetectionEnabled(scip, quiet, TRUE);

   return SCIP_OKAY;
}


/** @brief disables detectors
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
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPsetBoolParam(scip, paramname, FALSE) );
      if( !quiet )
      {
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detection/detectors/%s/postprocessingenabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPsetBoolParam(scip, paramname, FALSE) );
      if( !quiet )
      {
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   setDetectionEnabled(scip, quiet, FALSE);

   return SCIP_OKAY;
}


/** @brief sets the parameters to fast values
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

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/miplibconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetStage(scip) >= SCIP_STAGE_PROBLEM && SCIPgetNVars(scip) + SCIPgetNConss(scip) < FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER )
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/classification/consclassifier/consnamelevenshtein/enabled", FALSE) );

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

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
                             "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
      }
   }

   setDetectionEnabled(scip, quiet, TRUE);

   return SCIP_OKAY;
}


/** @brief method to calculate the greatest common divisor
 * @returns greatest common divisor
 */
static
int gcd(
   int a,   /**< first value */
   int b    /**< second value */
   )
{
   return b == 0 ? a : gcd( b, a % b );
}


 /**  */
static
SCIP_DECL_PARAMCHGD(paramChgdScore)
{
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   char* shortname;

   shortname = SCIPparamGetString(param);

   for( int i = 0; i < conshdlrdata->nscores; ++i )
   {
      GCG_SCORE* score;
      score = conshdlrdata->scores[i];

      assert(score != NULL);

      if( strcmp(score->shortname, shortname) == 0 )
      {
         return SCIP_OKAY;
      }
   }
 
   SCIPdialogMessage(scip, NULL, "The score <%s> does not exist.\n", shortname);
   SCIPdialogMessage(scip, NULL, "The default score <%s> is selected.\n", DEFAULT_SCORE);

   SCIPsetStringParam(scip, "detection/scores/selected", DEFAULT_SCORE);
 
   return SCIP_OKAY;
}


/* --------------- public functions --------------- */

/* prints block candidate information (gcg.h) */
SCIP_RETCODE GCGprintBlockcandidateInformation(
   SCIP*                 scip,
   FILE*                 file
   )
{
   gcg::DETPROBDATA* detprobdata;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   detprobdata = (conshdlrdata->detprobdatapres == NULL ? conshdlrdata->detprobdataorig : conshdlrdata->detprobdatapres );

   if( detprobdata == NULL )
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), NULL, "No block number candidates are calculated yet, consider detecting first..  \n" );
   else
      detprobdata->printBlockcandidateInformation(scip, file);

   return SCIP_OKAY;
}


/* prints detection time (gcg.h) */
SCIP_RETCODE GCGprintCompleteDetectionTime(
   SCIP*                 givenscip,
   FILE*                 file
   )
{
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "DETECTIONTIME   \n" );
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(givenscip), file, "%f \n", (SCIP_Real) GCGconshdlrDecompGetCompleteDetectionTime(givenscip) );

   return SCIP_OKAY;
}


/* prints partition information (gcg.h) */
SCIP_RETCODE GCGprintPartitionInformation(
   SCIP*                 scip,
   FILE*                 file
   )
{
   gcg::DETPROBDATA* detprobdata;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   detprobdata = (conshdlrdata->detprobdatapres == NULL ? conshdlrdata->detprobdataorig : conshdlrdata->detprobdatapres );

   detprobdata->printPartitionInformation(file);

   return SCIP_OKAY;
}


/* prints decomp information (gcg.h) */
SCIP_RETCODE GCGprintDecompInformation(
   SCIP*                 scip,
   FILE*                 file
   )
{
   std::vector<gcg::PARTIALDECOMP*>::const_iterator partialdeciter;
   std::vector<gcg::PARTIALDECOMP*>::const_iterator partialdeciterend;

   assert( GCGconshdlrDecompCheckConsistency(scip) );

   std::vector<PARTIALDECOMP*> partialdeclist;
   getFinishedPartialdecs(scip, partialdeclist);
   partialdeciter = partialdeclist.begin();
   partialdeciterend = partialdeclist.end();

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "DECOMPINFO  \n" );
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (int) partialdeclist.size() );

   for( ; partialdeciter != partialdeciterend; ++partialdeciter )
   {
      gcg::PARTIALDECOMP* partialdec;
      int nblocks = (*partialdeciter)->getNBlocks();

      partialdec = *partialdeciter;

      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "NEWDECOMP  \n");
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (*partialdeciter)->getNBlocks());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (*partialdeciter)->getID());
      for( int block = 0; block < nblocks; ++block )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNConssForBlock(block));
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNVarsForBlock(block));
      }
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNMasterconss());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNLinkingvars());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNMastervars());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", partialdec->getNTotalStairlinkingvars());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  partialdec->getMaxWhiteScore());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  partialdec->getScore(GCGconshdlrDecompFindScore(scip, "classic")));
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%f\n",  partialdec->getScore(GCGconshdlrDecompFindScore(scip, "max foreseeing white")));
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n",  partialdec->hasSetppccardMaster());
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (int) partialdec->getDetectorchain( ).size());
      for( int detector = 0; detector <(int) partialdec->getDetectorchain().size(); ++ detector )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%s\n",
                               GCGdetectorGetName(partialdec->getDetectorchain()[detector]));
      }
      partialdec->printPartitionInformation(scip, file);
   }

   return SCIP_OKAY;
}


/* gets number of partialdecs for public interface (pub_decomp.h)*/
int GCGgetNDecomps(
   SCIP* scip
   )
{
   return GCGconshdlrDecompGetNDecomps(scip);
}


char GCGdetectorGetChar(
   GCG_DETECTOR*         detector
   )
{
   if( detector == NULL )
      return '0';
   else
      return detector->decchar;
}


GCG_DETECTORDATA* GCGdetectorGetData(
   GCG_DETECTOR*         detector
   )
{
   assert(detector != NULL);
   return detector->decdata;
}


const char* GCGdetectorGetName(
   GCG_DETECTOR*         detector
   )
{
   assert(detector != NULL);
   return detector->name;
}


SCIP_RETCODE GCGdetectStructure(
   SCIP*                 scip,
   SCIP_RESULT*          result
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   SCIP_CALL(SCIPresetClock(scip, conshdlrdata->completedetectionclock));
   SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock));

   *result = SCIP_DIDNOTRUN;

   if( conshdlrdata->detprobdataorig == NULL )
   {
       resetDetprobdata(scip, TRUE);
   }

   /* if there are no entries there cannot be a structure */
   if( SCIPgetNOrigVars(scip) == 0 && SCIPgetNOrigConss(scip) == 0 )
      return SCIP_OKAY;

   /* if the original problem should be solved, then no decomposition will be performed */
   if( GCGgetDecompositionMode(scip) == GCG_DECMODE_ORIGINAL )
      return SCIP_OKAY;

   /* if there should not be any detection, stop at this point */
   if(!conshdlrdata->enabled)
   {
      return SCIP_OKAY;
   }

   ////////////////////
   ///// ORIGINAL /////
   ////////////////////

   SCIP_Bool calculateOrigDecomps;
   SCIP_Bool classifyOrig;

   SCIPgetBoolParam(scip, "detection/origprob/enabled", &calculateOrigDecomps);
   SCIPgetBoolParam(scip, "detection/origprob/classificationenabled", &classifyOrig);
   if (SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED) {
      if (calculateOrigDecomps) {
         // if there is no root partialdec yet, add root partialdec
         if( conshdlrdata->detprobdataorig->getOpenPartialdecs().empty() )
         {
            PARTIALDECOMP* rootpartialdec = new PARTIALDECOMP(scip, true);
            bool success = conshdlrdata->detprobdataorig->addPartialdecToOpen(rootpartialdec);
            if( !success )
            {
               SCIPerrorMessage("Could not add root partialdecomp to the pool of open decompositions.");
               return SCIP_ERROR;
            }
         }

         SCIPdebugMessage("is stage < transformed ? %s -> do %s transformProb() ", (SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED ? "yes" : "no"), (SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED ? "" : "not")  );

         if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
            SCIP_CALL(SCIPtransformProb(scip));

         // TODO ref classify (introduce flag for each detector and check whether this is needed)

         // CLASSIFICATION
         if( classifyOrig )
         {
            GCGconshdlrDecompClassify(scip, FALSE);
            if( SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
               conshdlrdata->detprobdataorig->printBlockcandidateInformation(scip, NULL);
         }
         else
            SCIPdebugMessage("classification for orig problem disabled \n" );

         // BLOCK CANDIDATES
         GCGconshdlrDecompCalcCandidatesNBlocks(scip, FALSE);

         // FIND DECOMPOSITIONS
         SCIPdebugMessage("start finding decompositions for original problem!\n" );
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "start finding decompositions for original problem!\n");
         SCIP_CALL(SCIPresetClock(scip, conshdlrdata->detectorclock));
         SCIP_CALL(SCIPstartClock(scip, conshdlrdata->detectorclock));
         detect(scip, conshdlrdata->detprobdataorig);
         SCIP_CALL(SCIPstopClock(scip, conshdlrdata->detectorclock));
         SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "finished finding decompositions for original problem!\n");
         SCIPdebugMessage("finished finding decompositions for original problem!\n" );
      }
      else
         SCIPdebugMessage("finding decompositions for original problem is NOT enabled!\n");

      std::vector<SCIP_CONS*> indexToCons; /* stores the corresponding scip constraints pointer */
      std::vector<gcg::PARTIALDECOMP*> partialdecsorig(0); /* partialdecs that were found for the orig problem */

      SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock));
      conshdlrdata->hasrunoriginal = TRUE;
      conshdlrdata->detprobdataorig->freeTemporaryData();
   }
   else 
   { // SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED
      //////////////////////////////////////
      ///// TRANSFORMED / IS PRESOLVED /////
      //////////////////////////////////////

      /* detection for presolved problem */

      if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

         /* presolving removed all constraints or variables */
         if (SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0)
            conshdlrdata->hasrun = TRUE;

         *result = SCIP_DIDNOTRUN;
         return SCIP_OKAY;
      }

      /* start detection clocks */
      SCIP_CALL(SCIPresetClock(scip, conshdlrdata->completedetectionclock));
      SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock));

      /* Classification */
      GCGconshdlrDecompClassify(scip, TRUE);
      GCGconshdlrDecompCalcCandidatesNBlocks(scip, TRUE);

      /* add block number candidates of the original problem */
      if( conshdlrdata->detprobdataorig )
      {
         for( auto& candidatesNBlock : conshdlrdata->detprobdataorig->candidatesNBlocks )
            conshdlrdata->detprobdatapres->addCandidatesNBlocksNVotes(candidatesNBlock.first, candidatesNBlock.second);
      }

      /* add root partialdec */
      if( conshdlrdata->detprobdatapres->getOpenPartialdecs().empty() )
      {
         PARTIALDECOMP* rootpartialdec = new PARTIALDECOMP(scip, false);

         bool success = conshdlrdata->detprobdatapres->addPartialdecToOpen(rootpartialdec);
         assert(success);
         if( !success )
         {
            SCIPerrorMessage("Could not add root decomposition.");
            *result = SCIP_DIDNOTRUN;
            return SCIP_ERROR;
         }
      }

      SCIP_CALL(SCIPresetClock(scip, conshdlrdata->detectorclock));
      SCIP_CALL(SCIPstartClock(scip, conshdlrdata->detectorclock));
      detect(scip, conshdlrdata->detprobdatapres);
      SCIP_CALL(SCIPstopClock(scip, conshdlrdata->detectorclock));
      conshdlrdata->detprobdatapres->sortFinishedForScore();
      SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock));
      conshdlrdata->hasrun = TRUE;
      conshdlrdata->detprobdatapres->freeTemporaryData();
   }

   /////////////////////////
   ///// EVAL SUCCESS //////
   /////////////////////////

   if( conshdlrdata->detprobdatapres != NULL && conshdlrdata->detprobdatapres->getNFinishedPartialdecs() > 0 )
      *result = SCIP_SUCCESS;

   if( conshdlrdata->detprobdataorig != NULL &&  conshdlrdata->detprobdataorig->getNFinishedPartialdecs() > 0 )
      *result = SCIP_SUCCESS;

   SCIPdebugMessage("Detection took %fs\n", SCIPgetClockTime( scip, conshdlrdata->detectorclock));

   if( conshdlrdata->detprobdatapres != NULL && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
      conshdlrdata->detprobdatapres->printBlockcandidateInformation(scip, NULL);

   /* display timing statistics */
   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, "Detection Time: %.2f\n", GCGconshdlrDecompGetCompleteDetectionTime(scip));
   /** @todo put this output to the statistics output */

   if( *result == SCIP_DIDNOTRUN )
   {
      return SCIP_OKAY;
   }

   /* show that we have done our duty */
   *result = SCIP_SUCCESS;

   return SCIP_OKAY;
}


GCG_CONSCLASSIFIER* GCGfindConsClassifier(
   SCIP*                 scip,
   const char*           name
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nconsclassifiers; ++i )
   {
      GCG_CONSCLASSIFIER *classifier;
      classifier = conshdlrdata->consclassifiers[i];
      assert(classifier != NULL);
      if( strcmp(classifier->name, name) == 0 )
      {
         return classifier;
      }
   }

   return NULL;
}


GCG_VARCLASSIFIER* GCGfindVarClassifier(
   SCIP*                 scip,
   const char*           name
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nvarclassifiers; ++i )
   {
      GCG_VARCLASSIFIER *classifier;
      classifier = conshdlrdata->varclassifiers[i];
      assert(classifier != NULL);
      if( strcmp(classifier->name, name) == 0 )
      {
         return classifier;
      }
   }

   return NULL;
}


GCG_DETECTOR* GCGfindDetector(
   SCIP*                 scip,
   const char*           name
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      GCG_DETECTOR* detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);
      if( strcmp(detector->name, name) == 0 )
      {
         return detector;
      }
   }

   return NULL;
}


GCG_SCORE* GCGconshdlrDecompFindScore(
   SCIP*                 scip,
   const char*           name
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(name != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nscores; ++i )
   {
      GCG_SCORE* score;
      score = conshdlrdata->scores[i];
      assert(score != NULL);
      if( strcmp(score->name, name) == 0 )
      {
         return score;
      }
   }

   return NULL;
}


GCG_SCORE* GCGconshdlrDecompFindScoreByShortname(
   SCIP*                 scip,
   const char*           shortname
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   assert(scip != NULL);
   assert(shortname != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
      return NULL;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   for( i = 0; i < conshdlrdata->nscores; ++i )
   {
      GCG_SCORE* score;
      score = conshdlrdata->scores[i];
      assert(score != NULL);
      if( strcmp(score->shortname, shortname) == 0 )
      {
         return score;
      }
   }

   return NULL;
}


GCG_DECOMP* GCGgetBestDecomp(
   SCIP*                 scip,
   SCIP_Bool             printwarnings
   )
{
   GCG_DECOMP* decomp;
   PARTIALDECOMP* partialdec;
   std::vector<std::pair<PARTIALDECOMP*, SCIP_Real> > candidates;

   if ( SCIPgetStage(scip) < SCIP_STAGE_PROBLEM )
      return NULL;

   GCGconshdlrDecompChooseCandidatesFromSelected(scip, candidates, FALSE, printwarnings);
   if ( candidates.empty() )
      return NULL;

   partialdec = candidates[0].first;

   assert( !partialdec->isAssignedToOrigProb() );
   assert( partialdec->isComplete() );

   createDecompFromPartialdec(scip, partialdec, &decomp);

   return decomp;
}


PARTIALDECOMP* GCGgetPartialdecToWrite(
   SCIP*                         scip,
   SCIP_Bool                     transformed
   )
{
   int dec;
   std::vector<std::pair<PARTIALDECOMP*, SCIP_Real> > candidates;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   // see if this is a call from functions like e.g. /see GCGwriteAllDecomps
   if( conshdlrdata->partialdectowrite != NULL )
   {
      return conshdlrdata->partialdectowrite;
   }

   GCGconshdlrDecompChooseCandidatesFromSelected(scip, candidates, !transformed, TRUE);

   // if none was found, output "pure" problem
   if( candidates.empty() )
   {
      int id = GCGconshdlrDecompAddMatrixPartialdec(scip, transformed);
      return GCGconshdlrDecompGetPartialdecFromID(scip, id);
   }

   // get the index of the next fitting candidate
   for( dec = 0; dec < (int) candidates.size(); ++dec )
   {
      if( candidates[dec].first->isAssignedToOrigProb() != transformed )
         break;
   }
   // return the candidate iff one was found
   if( dec != (int) candidates.size() )
      return candidates[dec].first;

   return NULL;
}


SCIP_RETCODE GCGgetPartialdecToWrite(
   SCIP*                         scip,
   SCIP_Bool                     transformed,
   PARTIALDECOMP_WRAPPER*        partialdecwrapper
   )
{
   partialdecwrapper->partialdec = GCGgetPartialdecToWrite(scip, transformed);

   return SCIP_OKAY;
}


SCIP_Real GCGgetRemainingTime(
   SCIP*                 scip 
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}


SCIP_RETCODE GCGincludeConsClassifier(
   SCIP*                 scip,
   const char*           name,
   const char*           description,
   int                   priority,
   SCIP_Bool             enabled,
   GCG_CLASSIFIERDATA*   classifierdata,
   GCG_DECL_FREECONSCLASSIFIER((*freeClassifier)),
   GCG_DECL_CONSCLASSIFY((*classify))
   )
{
   GCG_CONSCLASSIFIER* classifier = NULL;

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   SCIP_CALL( SCIPallocBlockMemory(scip, &classifier) );
   assert(classifier != NULL);

   SCIPdebugMessage("Adding classifier %i: %s\n", conshdlrdata->nconsclassifiers+1, name);

   SCIP_ALLOC( BMSduplicateMemoryArray(&classifier->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&classifier->description, description, strlen(description)+1) );
   
   classifier->priority = priority;
   classifier->enabled = enabled;
   classifier->clsdata = classifierdata;

   classifier->freeClassifier = freeClassifier;
   classifier->classify = classify;

   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/classification/consclassifier/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether constraint classifier for <%s> is enabled", description);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(classifier->enabled), FALSE, enabled, NULL, NULL) );

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->consclassifiers, (size_t)conshdlrdata->nconsclassifiers+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->consclassifierpriorities,(size_t) conshdlrdata->nconsclassifiers+1) );

   conshdlrdata->consclassifiers[conshdlrdata->nconsclassifiers] = classifier;
   conshdlrdata->nconsclassifiers = conshdlrdata->nconsclassifiers+1;

   return SCIP_OKAY;
}


SCIP_RETCODE GCGincludeDetector(
   SCIP*                 scip,
   const char*           name,
   const char            decchar,
   const char*           description,
   int                   freqCallRound,
   int                   maxCallRound,
   int                   minCallRound,
   int                   freqCallRoundOriginal,
   int                   maxCallRoundOriginal,
   int                   minCallRoundOriginal,
   int                   priority,
   SCIP_Bool             enabled,
   SCIP_Bool             enabledFinishing,
   SCIP_Bool             enabledPostprocessing,
   SCIP_Bool             skip,
   SCIP_Bool             usefulRecall,
   GCG_DETECTORDATA*     detectordata,
   GCG_DECL_FREEDETECTOR((*freeDetector)),
   GCG_DECL_INITDETECTOR((*initDetector)),
   GCG_DECL_EXITDETECTOR((*exitDetector)),
   GCG_DECL_PROPAGATEPARTIALDEC((*propagatePartialdecDetector)),
   GCG_DECL_FINISHPARTIALDEC((*finishPartialdecDetector)),
   GCG_DECL_POSTPROCESSPARTIALDEC((*postprocessPartialdecDetector)),
   GCG_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),
   GCG_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),
   GCG_DECL_SETPARAMFAST((*setParamFastDetector))
   )
{
   GCG_DETECTOR* detector = NULL;
   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   if( conshdlrdata == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, &detector) );
   assert(detector != NULL);

   SCIPdebugMessage("Adding detector %i: %s\n", conshdlrdata->ndetectors+1, name);

#ifndef NDEBUG
   assert(GCGfindDetector(scip, name) == NULL);
#endif

   /* set meta data of detector */
   detector->decdata = detectordata;
   SCIP_ALLOC( BMSduplicateMemoryArray(&detector->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&detector->description, description, strlen(description)+1) );
   detector->decchar = decchar;

   /* set memory handling and detection functions */
   detector->freeDetector = freeDetector;
   detector->initDetector = initDetector;
   detector->exitDetector = exitDetector;

   /* set functions for editing partialdecs */
   detector->propagatePartialdec = propagatePartialdecDetector;
   detector->finishPartialdec = finishPartialdecDetector;
   detector->postprocessPartialdec = postprocessPartialdecDetector;

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
   detector->enabledFinishing = enabledFinishing;
   detector->enabledPostprocessing = enabledPostprocessing;
   detector->skip = skip;
   detector->usefulRecall = usefulRecall;
   detector->overruleemphasis = FALSE;
   detector->ndecomps = 0;
   detector->ncompletedecomps = 0;
   detector->dectime = 0.;

   /* add and initialize all parameters accessable from menu */
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/detectors/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabled), FALSE, enabled, NULL, NULL) );

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
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be called on descendants of the current partialdec", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->usefulRecall), FALSE, usefulRecall, NULL, NULL) );

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

   /* add to detector array */
   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   // TODO ref replace propagatePartialdec by more appropriate name once partialdecs are introduced
   /* add to propagating detector array if appropriate */
   if( detector->propagatePartialdec != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->propagatingdetectors, (size_t)conshdlrdata->npropagatingdetectors+1) );
      conshdlrdata->propagatingdetectors[conshdlrdata->npropagatingdetectors] = detector;
      conshdlrdata->npropagatingdetectors = conshdlrdata->npropagatingdetectors+1;
   }

   // TODO ref replace finishPartialdec by more appropriate name once partialdecs are introduced
   /* add to finishing detector array if appropriate */
   if( detector->finishPartialdec != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->finishingdetectors, (size_t)conshdlrdata->nfinishingdetectors+1) );
      conshdlrdata->finishingdetectors[conshdlrdata->nfinishingdetectors] = detector;
      conshdlrdata->nfinishingdetectors = conshdlrdata->nfinishingdetectors+1;
   }

   // TODO ref replace postprocessPartialdec by more appropriate name once partialdecs are introduced
   /* add to postprocessing detector array if appropriate */
   if( detector->postprocessPartialdec != NULL )
   {
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->postprocessingdetectors, (size_t)conshdlrdata->npostprocessingdetectors+1) );
      conshdlrdata->postprocessingdetectors[conshdlrdata->npostprocessingdetectors] = detector;
      conshdlrdata->npostprocessingdetectors = conshdlrdata->npostprocessingdetectors+1;
   }

   return SCIP_OKAY;
}


SCIP_RETCODE GCGincludeVarClassifier(
   SCIP*                 scip,
   const char*           name,
   const char*           description,
   int                   priority,
   SCIP_Bool             enabled,
   GCG_CLASSIFIERDATA*   classifierdata,
   GCG_DECL_FREEVARCLASSIFIER((*freeClassifier)),
   GCG_DECL_VARCLASSIFY((*classify))
   )
{
   GCG_VARCLASSIFIER* classifier = NULL;

   assert(scip != NULL);
   assert(name != NULL);
   assert(description != NULL);

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   SCIP_CALL( SCIPallocBlockMemory(scip, &classifier) );
   assert(classifier != NULL);

   SCIPdebugMessage("Adding classifier %i: %s\n", conshdlrdata->nvarclassifiers+1, name);

   SCIP_ALLOC( BMSduplicateMemoryArray(&classifier->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&classifier->description, description, strlen(description)+1) );
   
   classifier->priority = priority;
   classifier->enabled = enabled;
   classifier->clsdata = classifierdata;

   classifier->freeClassifier = freeClassifier;
   classifier->classify = classify;

   char setstr[SCIP_MAXSTRLEN];
   char descstr[SCIP_MAXSTRLEN];
   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detection/classification/varclassifier/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether variable classifier for <%s> is enabled", description);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(classifier->enabled), FALSE, enabled, NULL, NULL) );

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->varclassifiers, (size_t)conshdlrdata->nvarclassifiers+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->varclassifierpriorities,(size_t) conshdlrdata->nvarclassifiers+1) );

   conshdlrdata->varclassifiers[conshdlrdata->nvarclassifiers] = classifier;
   conshdlrdata->nvarclassifiers = conshdlrdata->nvarclassifiers+1;

   return SCIP_OKAY;
}

char* GCGgetCurrentScoreShortname(
   SCIP*                 scip
   )
{
   assert(scip != NULL);

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   assert(conshdlrdata != NULL);

   return conshdlrdata->currscore;
}

GCG_SCORE* GCGgetCurrentScore(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;

   conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   char* shortname = GCGgetCurrentScoreShortname(scip);

   for( i = 0; i < conshdlrdata->nscores; ++i )
   {
      GCG_SCORE* score;
      score = conshdlrdata->scores[i];
      assert(score != NULL);
      if( strcmp(score->shortname, shortname) == 0 )
      {
         return score;
      }
   }

   return NULL;
}

SCIP_RETCODE GCGconshdlrDecompIncludeScore(
   SCIP*                 scip,
   const char*           name,
   const char*           shortname,
   const char*           description,
   GCG_SCOREDATA*        scoredata,
   GCG_DECL_SCOREFREE    ((*scorefree)),
   GCG_DECL_SCORECALC    ((*scorecalc))
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   GCG_SCORE* score = NULL;

   assert(scip != NULL);
   assert(name != NULL);
   assert(shortname != NULL);
   assert(description != NULL);

   conshdlrdata = getConshdlrdata(scip);

   SCIP_CALL( SCIPallocBlockMemory(scip, &score) );
   assert(score != NULL);

   SCIPdebugMessage("Adding score %i: %s\n", conshdlrdata->nscores+1, name);

   SCIP_ALLOC( BMSduplicateMemoryArray(&score->name, name, strlen(name)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&score->shortname, shortname, strlen(shortname)+1) );
   SCIP_ALLOC( BMSduplicateMemoryArray(&score->description, description, strlen(description)+1) );
   score->scoredata = scoredata;

   score->scorefree = scorefree;
   score->scorecalc = scorecalc;

   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->scores, (size_t)conshdlrdata->nscores+1) );

   conshdlrdata->scores[conshdlrdata->nscores] = score;
   conshdlrdata->nscores = conshdlrdata->nscores+1;

   return SCIP_OKAY;
}

void GCGprintListOfDetectors(
   SCIP*                 scip
   )
{
   int ndetectors;
   int i;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   ndetectors = conshdlrdata->ndetectors;

   SCIPdialogMessage(scip, NULL, " detector             char priority enabled  description\n");
   SCIPdialogMessage(scip, NULL, " --------------       ---- -------- -------  -----------\n");

   for( i = 0; i < ndetectors; ++i )
   {
      SCIPdialogMessage(scip, NULL,  " %-20s    ", conshdlrdata->detectors[i]->name);
      if( strlen(conshdlrdata->detectors[i]->name) > 20 )
         SCIPdialogMessage(scip, NULL, "\n %20s    ", "-->");
      SCIPdialogMessage(scip, NULL,  "%c ", conshdlrdata->detectors[i]->decchar);
      SCIPdialogMessage(scip, NULL,  "%8d ", conshdlrdata->detectors[i]->priority);
      SCIPdialogMessage(scip, NULL,  "%7s  ", conshdlrdata->detectors[i]->enabled ? "TRUE" : "FALSE");
      SCIPdialogMessage(scip, NULL,  "%s", conshdlrdata->detectors[i]->description);
      SCIPdialogMessage(scip, NULL, "\n");
   }
}


SCIP_RETCODE GCGwriteAllDecomps(
   SCIP*                 scip,
   char*                 directory,
   char*                 extension,
   SCIP_Bool             original,
   SCIP_Bool             presolved
   )
{
   char outname[SCIP_MAXSTRLEN];
   char tempstring[SCIP_MAXSTRLEN];
   int i;

   int maxtowrite;
   int nwritten;

   assert(extension != NULL);

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   maxtowrite = -1;
   nwritten = 0;

   if( presolved && conshdlrdata->detprobdatapres != NULL && conshdlrdata->detprobdatapres->getNFinishedPartialdecs() == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   SCIPgetIntParam(scip, "visual/nmaxdecompstowrite", &maxtowrite );

   /* write all finished partialdecs */
   std::vector<PARTIALDECOMP*> partialdecs;
   getFinishedPartialdecs(scip, partialdecs);
   for( i = 0; i < (int) partialdecs.size(); i++)
   {
      PARTIALDECOMP* partialdec = partialdecs.at(i);

      /* get filename */
      GCGgetVisualizationFilename(scip, partialdec, extension, tempstring);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s.%s", directory, tempstring, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s.%s", tempstring, extension);
      }

      /* mark partialdec */
      conshdlrdata->partialdectowrite = partialdec;

      /* these functions write the marked partialdec */
      if( partialdec->isAssignedToOrigProb() )
         SCIP_CALL_QUIET( SCIPwriteOrigProblem(scip, outname, extension, FALSE) );
      else
         SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );

      ++nwritten;
      conshdlrdata->partialdectowrite = NULL;

      /* stop if max is reached */
      if( maxtowrite != -1 && nwritten >= maxtowrite )
         break;
   }

   return SCIP_OKAY;
}


SCIP_RETCODE GCGwriteSelectedDecomps(
   SCIP*                 scip,
   char*                 directory,
   char*                 extension
   )
{
   char outname[SCIP_MAXSTRLEN];
   char tempstring[SCIP_MAXSTRLEN];

   assert(extension != NULL);

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->partialdecs->empty() )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   std::vector<PARTIALDECOMP*> selectedpartialdecs;
   getSelectedPartialdecs(scip, selectedpartialdecs);
   if( selectedpartialdecs.empty() )
   {
      SCIPwarningMessage(scip, "No decomposition selected.\n");
      return SCIP_OKAY;
   }

   for(auto & partialdec : selectedpartialdecs)
   {
      GCGgetVisualizationFilename(scip, partialdec, extension, tempstring);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s.%s", directory, tempstring, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s.%s", tempstring, extension);
      }

      conshdlrdata->partialdectowrite = partialdec;

      if ( partialdec->isAssignedToOrigProb() )
         SCIP_CALL_QUIET( SCIPwriteOrigProblem(scip, outname, extension, FALSE) );
      else
         SCIP_CALL_QUIET( SCIPwriteTransProblem(scip, outname, extension, FALSE) );

      conshdlrdata->partialdectowrite = NULL;
   }

   return SCIP_OKAY;
}


int GCGconshdlrDecompAddBasicPartialdec(
   SCIP* scip,
   SCIP_Bool presolved
   )
{
   PARTIALDECOMP* partialdec = new PARTIALDECOMP(scip, !presolved);
   partialdec->setNBlocks(0);
   partialdec->assignOpenConssToMaster();
   partialdec->prepare();
   addPartialdec(scip, partialdec);
   return partialdec->getID();
}


void GCGconshdlrDecompAddCandidatesNBlocks(
   SCIP* scip,
   SCIP_Bool origprob,
   int candidate
   )
{
   gcg::DETPROBDATA* detprobdata;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   detprobdata = origprob ? conshdlrdata->detprobdataorig : conshdlrdata->detprobdatapres;

   if( candidate > 1 )
   {
      bool alreadyIn = false;
      for( size_t i = 0; i < detprobdata->candidatesNBlocks.size(); ++ i )
      {
         if( detprobdata->candidatesNBlocks[i].first == candidate )
         {
            alreadyIn = true;
            ++detprobdata->candidatesNBlocks[i].second;
            break;
         }
      }
      if( !alreadyIn )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, "added block number candidate: %d \n", candidate );
         detprobdata->candidatesNBlocks.emplace_back(candidate, 1);
      }
   }
}


SCIP_RETCODE GCGconshdlrDecompAddDecomp(
   SCIP*                 scip,
   GCG_DECOMP*           decomp,
   SCIP_Bool             select
   )
{
   PARTIALDECOMP* partialdec;

   if( decomp->presolved && SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
   {
      SCIPerrorMessage("Problem is not presolved yet.");
      return SCIP_ERROR;
   }

   SCIP_CALL( createPartialdecFromDecomp(scip, decomp, &partialdec) );
   SCIP_CALL( addPartialdec(scip, partialdec) );
   partialdec->setSelected(select);

   return SCIP_OKAY;
}


int GCGconshdlrDecompAddMatrixPartialdec(
   SCIP* scip,
   SCIP_Bool presolved
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   DETPROBDATA* detprobdata;
   // check if detection took place (this might be used in e.g. write problem)
   if(presolved)
   {
      if(!GCGconshdlrDecompPresolvedDetprobdataExists(scip))
         resetDetprobdata(scip, false);
      detprobdata = conshdlrdata->detprobdatapres;
   }
   else
   {
      if(!GCGconshdlrDecompOrigDetprobdataExists(scip))
         resetDetprobdata(scip, true);
       detprobdata = conshdlrdata->detprobdataorig;
   }
   
   assert(detprobdata != NULL);

   PARTIALDECOMP* matrixpartialdec;

   matrixpartialdec = new PARTIALDECOMP(scip, !presolved);
   matrixpartialdec->setNBlocks(1);

   for( int i = 0; i < detprobdata->getNConss(); ++i )
      matrixpartialdec->fixConsToBlock(i,0);

   for( int i = 0; i < detprobdata->getNVars(); ++i )
      matrixpartialdec->fixVarToBlock(i,0);

   matrixpartialdec->sort();

   detprobdata->addPartialdecToFinishedUnchecked(matrixpartialdec);

   return matrixpartialdec->getID();
}


SCIP_RETCODE GCGconshdlrDecompAddPreexistingDecomp(
   SCIP*                 scip,
   GCG_DECOMP*           decomp
   )
{
   PARTIALDECOMP* partialdec;

   if( decomp->presolved && SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
   {
      SCIPerrorMessage("Problem is not presolved yet.");
      return SCIP_ERROR;
   }

   SCIP_CALL( createPartialdecFromDecomp(scip, decomp, &partialdec) );
   GCGconshdlrDecompAddPreexisitingPartialDec(scip, partialdec);

   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompAddPreexisitingPartialDec(
   SCIP* scip,
   PARTIALDECOMP* partialdec
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   bool assignedconss = false;
   assert(conshdlrdata != NULL);
   assert( partialdec != NULL );

   if( partialdec->shouldCompletedByConsToMaster() )
   {
      auto& openconss = partialdec->getOpenconssVec();
      for( auto itr = openconss.cbegin(); itr != openconss.cend(); )
      {
         itr = partialdec->fixConsToMaster(itr);
         assignedconss = true;
      }
      partialdec->sort();
   }

   partialdec->prepare();
#ifndef NDEBUG
   if( partialdec->getUsergiven() == USERGIVEN::COMPLETE || partialdec->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER )
      assert( partialdec->isComplete() );
#endif

   if( partialdec->isComplete() )
   {
      if( !assignedconss )
         partialdec->setUsergiven( USERGIVEN::COMPLETE );
      addPartialdec(scip, partialdec);

      /* if detprobdata for presolved problem already exist try to translate partialdec */
      if ( conshdlrdata->detprobdatapres != NULL && partialdec->isAssignedToOrigProb())
      {
         std::vector<PARTIALDECOMP*> partialdectotranslate;
         partialdectotranslate.push_back(partialdec);
         std::vector<PARTIALDECOMP*> newpartialdecs = conshdlrdata->detprobdatapres->translatePartialdecs(conshdlrdata->detprobdataorig, partialdectotranslate);
         if( !newpartialdecs.empty() )
         {
            addPartialdec(scip, newpartialdecs[0]);
         }
      }
   }
   else
   {
      partialdec->setUsergiven( USERGIVEN::PARTIAL );
      addPartialdec(scip, partialdec);
   }

   /* set statistics */
   int nvarstoblock = 0;
   int nconsstoblock = 0;

   for ( int b = 0; b < partialdec->getNBlocks(); ++b )
   {
      nvarstoblock += partialdec->getNVarsForBlock(b);
      nconsstoblock += partialdec->getNConssForBlock(b);
   }

   partialdec->findVarsLinkingToMaster();
   partialdec->findVarsLinkingToStairlinking();

   char const* usergiveninfo;
   char const* presolvedinfo;

   if( partialdec->getUsergiven() == USERGIVEN::PARTIAL )
      usergiveninfo = "partial";
   if( partialdec->getUsergiven() == USERGIVEN::COMPLETE )
      usergiveninfo = "complete";
   if( partialdec->getUsergiven() == USERGIVEN::COMPLETED_CONSTOMASTER )
      usergiveninfo = "complete";
   if( partialdec->isAssignedToOrigProb() )
      presolvedinfo = "original";
   else presolvedinfo = "presolved";

   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " added %s decomp for %s problem with %d blocks and %d masterconss, %d linkingvars, "
                                                    "%d mastervars, and max white score of %s %f \n", usergiveninfo, presolvedinfo,
                   partialdec->getNBlocks(), partialdec->getNMasterconss(),
                   partialdec->getNLinkingvars(), partialdec->getNMastervars(), (partialdec->isComplete() ? " " : " at best "),
                   partialdec->getScore(GCGconshdlrDecompFindScore(scip, "max white")) );

   return SCIP_OKAY;
}


SCIP_CLOCK* GCGconshdlrDecompGetScoreClock(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   return conshdlrdata->scoreclock;
}


SCIP_Real GCGconshdlrDecompGetScoreTotalTime(
   SCIP* scip
)
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   return SCIPgetClockTime(scip, conshdlrdata->scoreclock);
}


void GCGconshdlrDecompAddUserCandidatesNBlocks(
   SCIP* scip,
   int candidate
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   // make sure there is a problem
   if(SCIPgetStage(scip) < SCIP_STAGE_PROBLEM)
   {
      SCIPdialogMessage(scip, NULL,
         "Please add a problem before adding block candidates.\n");
      return;
   }

   conshdlrdata->userblocknrcandidates->push_back(candidate);

   SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "added user block number candidate: %d \n", candidate );
}


SCIP_RETCODE GCGconshdlrDecompArePricingprobsIdenticalForPartialdecid(
   SCIP*                scip,
   int                  partialdecid,
   int                  probnr1,
   int                  probnr2,
   SCIP_Bool*           identical
   )
{
   gcg::PARTIALDECOMP* partialdec;

   partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);
   assert(partialdec != NULL);
   assert(partialdec->isComplete());

   if( !partialdec->aggInfoCalculated() )
   {
      SCIPdebugMessage("calc aggregation information for partialdec!\n");
      // ignore limits of detection since we do this only for this specific partialdec
      partialdec->calcAggregationInformation(true);
   }

   if( partialdec->getRepForBlock(probnr1) == partialdec->getRepForBlock(probnr2) )
      *identical = TRUE;
   else
      *identical = FALSE;

   SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL, NULL, " block %d and block %d are represented by %d and %d hence they are identical=%d.\n", probnr1, probnr2, partialdec->getRepForBlock(probnr1), partialdec->getRepForBlock(probnr2), *identical );

   return SCIP_OKAY;
}


void GCGconshdlrDecompCalcCandidatesNBlocks(
   SCIP* scip,
   SCIP_Bool transformed
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   gcg::DETPROBDATA* detprobdata = (transformed ? conshdlrdata->detprobdatapres : conshdlrdata->detprobdataorig );
   /* strategy: for every subset of constraint classes and variable classes calculate gcd (greatest common divisors)
    * of the corresponding number of constraints/variables assigned to this class */

   SCIP_CLOCK* nblockcandnclock;
   SCIPcreateClock( scip, & nblockcandnclock );
   SCIPstartClock( scip, nblockcandnclock );

   /* if  distribution of classes exceeds this number it is skipped */
   int maximumnclasses;
   SCIP_Bool medianvarspercons;

   /* used for nvars / medianofnvars per conss */
   std::vector<int> nvarspercons(0);
   std::list<int>::iterator iter;
   int candidate = -1;

   medianvarspercons = conshdlrdata->blocknumbercandsmedianvarspercons;
   maximumnclasses = conshdlrdata->maxnclassesfornblockcandidates;

   /* firstly, iterate over all conspartitions */
   for( auto partition : detprobdata->conspartitioncollection )
   {
      /* check if there are too many classes in this distribution and skip it if so */
      if( partition->getNClasses() > maximumnclasses )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " the current consclass distribution includes %d classes but only %d are allowed for GCGconshdlrDecompCalcCandidatesNBlocks()\n", partition->getNClasses(), maximumnclasses );
         continue;
      }

      /* get necessary data of current partition */
      std::vector<std::vector<int>> subsetsOfConstypes = partition->getAllSubsets(true, true, true );
      std::vector<int> nConssOfClasses = partition->getNConssOfClasses();

      /* start with the cardinalities of the consclasses as candidates */
      for( size_t i = 0; i < nConssOfClasses.size(); ++ i )
      {
         GCGconshdlrDecompAddCandidatesNBlocks( scip, detprobdata->isAssignedToOrigProb(), nConssOfClasses[i] );
      }

      /* continue with gcd of all cardinalities in this subset */
      for( size_t subset = 0; subset < subsetsOfConstypes.size(); ++ subset )
      {
         int greatestCD = 1;

         if( subsetsOfConstypes[subset].empty() || subsetsOfConstypes[subset].size() == 1 )
            continue;

         greatestCD = gcd( nConssOfClasses[subsetsOfConstypes[subset][0]], nConssOfClasses[subsetsOfConstypes[subset][1]] );

         for( size_t i = 2; i < subsetsOfConstypes[subset].size(); ++ i )
         {
            greatestCD = gcd( greatestCD, nConssOfClasses[subsetsOfConstypes[subset][i]] );
         }

         GCGconshdlrDecompAddCandidatesNBlocks( scip, detprobdata->isAssignedToOrigProb(), greatestCD );
      }
   }

   /* secondly, iterate over all varpartitions */
   for( auto partition : detprobdata->varpartitioncollection )
   {
      /* check if there are too many classes in this distribution and skip it if so */
      if( partition->getNClasses() > maximumnclasses )
      {
         SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " the current varclass distribution includes %d classes but only %d are allowed for GCGconshdlrDecompCalcCandidatesNBlocks()\n", partition->getNClasses(), maximumnclasses );
         continue;
      }

      /* get necessary data of current partition */
      std::vector<std::vector<int>> subsetsOfVartypes = partition->getAllSubsets(true, true, true, true );
      std::vector<int> nVarsOfClasses = partition->getNVarsOfClasses();

      /* start with the cardinalities of the varclasses as candidates */
      for( int nVarsOfClasse : nVarsOfClasses )
      {
         GCGconshdlrDecompAddCandidatesNBlocks( scip, detprobdata->isAssignedToOrigProb(), nVarsOfClasse );
      }

      /* continue with gcd of all cardinalities in this subset */
      for(auto& subsetsOfVartype : subsetsOfVartypes)
      {
         int greatestCD = 1;

         if( subsetsOfVartype.empty() || subsetsOfVartype.size() == 1 )
            continue;

         greatestCD = gcd( nVarsOfClasses[subsetsOfVartype[0]], nVarsOfClasses[subsetsOfVartype[1]] );

         for( size_t i = 2; i < subsetsOfVartype.size(); ++ i )
         {
            greatestCD = gcd( greatestCD, nVarsOfClasses[subsetsOfVartype[i]] );
         }

         GCGconshdlrDecompAddCandidatesNBlocks( scip, detprobdata->isAssignedToOrigProb(), greatestCD );
      }
   }

   /* block number candidate could be nvars / median of nvarsinconss  only calculate if desired*/
   if( medianvarspercons )
   {
      for( int c = 0; c < detprobdata->getNConss(); ++c )
      {
         nvarspercons.push_back( detprobdata->getNVarsForCons(c) );
      }
      std::sort(nvarspercons.begin(), nvarspercons.end() );
      candidate = (int) detprobdata->getNVars() / nvarspercons[(int)detprobdata->getNConss()/2];

      GCGconshdlrDecompAddCandidatesNBlocks( scip, detprobdata->isAssignedToOrigProb(), candidate);
   }

   SCIPstopClock( scip, nblockcandnclock );
   detprobdata->nblockscandidatescalctime = SCIPgetClockTime( scip, nblockcandnclock );
   SCIPfreeClock( scip, & nblockcandnclock );
}


SCIP_Bool GCGconshdlrDecompCheckConsistency(
   SCIP* scip  /* SCIP data structure */
   )
{
   int i;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   /* check whether the predecessors of all finished partialdecs in both detprobdatas can be found */
   if( conshdlrdata->detprobdatapres != NULL)
   {
      for( i = 0; i < conshdlrdata->detprobdatapres->getNFinishedPartialdecs(); ++i )
      {
         PARTIALDECOMP* partialdec = conshdlrdata->detprobdatapres->getFinishedPartialdec( i );

         for( int j = 0; j < partialdec->getNAncestors(); ++j )
         {
            int id = partialdec->getAncestorID( j );
            if(GCGconshdlrDecompGetPartialdecFromID(scip, id) == NULL )
            {
               SCIPwarningMessage(scip, "Warning: presolved partialdec %d has an ancestor (id: %d) that is not found! \n", partialdec->getID(), id );
               return FALSE;
            }
         }
      }
   }

   if( conshdlrdata->detprobdataorig != NULL )
   {
      for( i = 0; i < conshdlrdata->detprobdataorig->getNFinishedPartialdecs(); ++i )
      {
         PARTIALDECOMP* partialdec = conshdlrdata->detprobdataorig->getFinishedPartialdec( i );

         for( int j = 0; j < partialdec->getNAncestors(); ++j )
         {
            int id = partialdec->getAncestorID( j );
            if(GCGconshdlrDecompGetPartialdecFromID(scip, id) == NULL )
            {
               SCIPwarningMessage(scip, "Warning: orig partialdec %d has an ancestor (id: %d) that is not found! \n", partialdec->getID(), id );
               return FALSE;
            }
         }
      }
   }

   return TRUE;
}


SCIP_RETCODE GCGconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,
   std::vector<std::pair<gcg::PARTIALDECOMP*, SCIP_Real> >& candidates,
   SCIP_Bool original,
   SCIP_Bool printwarnings
   )
{
   std::vector<PARTIALDECOMP*>::iterator partialdeciter;
   std::vector<PARTIALDECOMP*>::iterator partialdeciterend;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   if( conshdlrdata == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot manage decompositions!\n");
      return SCIP_ERROR;
   }

   if( (!original && conshdlrdata->detprobdatapres == NULL) || (original && conshdlrdata->detprobdataorig == NULL) )
   {
       return SCIP_OKAY;
   }

   SCIPdebugMessage("Starting decomposition candidate choosing \n");

   assert( GCGconshdlrDecompCheckConsistency(scip) );

   std::vector<PARTIALDECOMP*> selectedpartialdecs;
   getSelectedPartialdecs(scip, selectedpartialdecs);

   if( selectedpartialdecs.empty() )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL, "currently no decomposition is selected, hence every known decomposition is considered: \n");
      if( original )
      {
         selectedpartialdecs = conshdlrdata->detprobdataorig->getFinishedPartialdecs();
      }
      else
      {
         selectedpartialdecs = conshdlrdata->detprobdatapres->getFinishedPartialdecs();
      }
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL,  "number of considered decompositions: %ld \n", selectedpartialdecs.size() );
   }

   partialdeciter = selectedpartialdecs.begin();
   partialdeciterend = selectedpartialdecs.end();

   /* get decomp candidates and calculate corresponding score */
   for( ; partialdeciter != partialdeciterend; ++partialdeciter )
   {
      PARTIALDECOMP* partialdec = *partialdeciter;
      if( !original && partialdec->isAssignedToOrigProb() )
      {
         partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdec->getTranslatedpartialdecid());
         assert(partialdec != NULL);
         if( !partialdec->isComplete() )
         {
             if( printwarnings )
             {
                SCIPwarningMessage(scip, "A selected decomposition (id=%d) of the orig. problem is ignored since its translation is incomplete.\n", partialdec->getID());
             }
             continue;
         }
      }

      assert(partialdec != NULL);
      if( partialdec->isComplete() )
      {
         candidates.emplace_back(partialdec, partialdec->getScore(GCGgetCurrentScore(scip)));
      }
      else if( printwarnings )
      {
          SCIPwarningMessage(scip, "A selected decomposition (id=%d) is ignored since it is incomplete.\n", partialdec->getID());
      }
   }

   /* sort decomp candidates according to score */
   std::sort( candidates.begin(), candidates.end(), sort_pred() );

   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompClassify(
   SCIP*                scip,
   SCIP_Bool            transformed
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   gcg::DETPROBDATA* detprobdata;
   if (transformed) {
      detprobdata = conshdlrdata->detprobdatapres;
   } else{
      detprobdata = conshdlrdata->detprobdataorig;
   }

   // Start clock
   SCIP_CLOCK* classificationclock;
   SCIPcreateClock( scip, & classificationclock );
   SCIP_CALL( SCIPstartClock( scip, classificationclock ) );

   // Cons classifiers
   for(int i = 0; i < conshdlrdata->nconsclassifiers; i++) {
      GCG_CONSCLASSIFIER *classifier;
      SCIP_Bool enabled;
      classifier = conshdlrdata->consclassifiers[i];

      char setting[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(setting, SCIP_MAXSTRLEN, "detection/classification/consclassifier/%s/enabled", classifier->name);

      SCIPgetBoolParam( scip, setting, &enabled );
      if ( enabled ) {
         classifier->classify(scip, classifier, transformed);
      }
   }

   // Var classifiers
   for(int i = 0; i < conshdlrdata->nvarclassifiers; i++) {
      GCG_VARCLASSIFIER *classifier;
      SCIP_Bool enabled;
      classifier = conshdlrdata->varclassifiers[i];

      char setting[SCIP_MAXSTRLEN];
      (void) SCIPsnprintf(setting, SCIP_MAXSTRLEN, "detection/classification/varclassifier/%s/enabled", classifier->name);

      SCIPgetBoolParam( scip, setting, &enabled );
      if ( enabled ) {
         classifier->classify(scip, classifier, transformed);
      }
   }

   // Reduce number of classes
   reduceConsclasses(scip, detprobdata);
   reduceVarclasses(scip, detprobdata);

   // Stop clock
   SCIP_CALL( SCIPstopClock( scip, classificationclock ) );
   detprobdata->classificationtime += SCIPgetClockTime( scip, classificationclock );
   SCIPfreeClock( scip, & classificationclock );

   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompCreateVarmapForPartialdecId(
   SCIP*                scip,
   SCIP_HASHMAP**       hashorig2pricingvar,
   int                  partialdecid,
   int                  probnr1,
   int                  probnr2,
   SCIP*                scip1,
   SCIP*                scip2,
   SCIP_HASHMAP*        varmap
   )
{
   gcg::PARTIALDECOMP* partialdec;
   gcg::DETPROBDATA* currdetprobdata;

   int blockid1;
   int blockid2;
   int representative;
   int repid1;
   int repid2;
   int nblocksforrep;
   std::vector<int> pidtopid;

   repid1 = -1;
   repid2 = -1;

   partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);
   assert(partialdec != NULL);
   assert(partialdec->isComplete());
   currdetprobdata = partialdec->getDetprobdata();

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

   representative = partialdec->getRepForBlock(blockid1);
   assert( representative == partialdec->getRepForBlock(blockid2) );
   nblocksforrep = (int) partialdec->getBlocksForRep(representative).size();

   /* find index in representatives */
   for( int i = 0; i < nblocksforrep; ++i )
   {
      if( partialdec->getBlocksForRep(representative)[i] == blockid1 )
         repid1 = i;
      if( partialdec->getBlocksForRep(representative)[i] == blockid2 )
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
      SCIPwarningMessage(scip, "blockid1 should be the representative (hence has id=0 in reptoblocksarray but in fact has %d) \n", repid1);
      return SCIP_OKAY;
   }

   pidtopid = partialdec->getRepVarmap(representative, repid2);

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
      var2origid = currdetprobdata->getIndexForVar(var2orig) ;
      assert(var2origid>=0);
      var2originblockid = partialdec->getVarProbindexForBlock(var2origid, blockid2) ;
      assert(var2originblockid >= 0);
      var1originblockid = pidtopid[var2originblockid];
      assert(var1originblockid>=0);
      var1origid = partialdec->getVarsForBlock(blockid1)[var1originblockid];
      assert(var1origid>=0);
      var1orig = currdetprobdata->getVar(var1origid) ;
      assert(var1orig != NULL);
      var1 = (SCIP_VAR*) SCIPhashmapGetImage(hashorig2pricingvar[blockid1], (void*) var1orig ) ;
      assert(var1 != NULL);

      SCIPhashmapInsert(varmap, (void*) var2, (void*) var1);
   }

   return SCIP_OKAY;
}


int GCGconshdlrDecompDecreaseNCallsCreateDecomp(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   --conshdlrdata->ncallscreatedecomp;

   return conshdlrdata->ncallscreatedecomp;
}


void GCGconshdlrDecompDeregisterPartialdecs(
   SCIP* scip,
   SCIP_Bool original
   )
{
   SCIP_CONSHDLRDATA *conshdlrdata = getConshdlrdata(scip);

   for( int i = (int) conshdlrdata->partialdecs->size() - 1; i >= 0; --i)
   {
      PARTIALDECOMP* partialdec = (*conshdlrdata->partialdecs)[i];
      if( partialdec->isAssignedToOrigProb() == (bool) original )
      {
         // ~PARTIALDECOMP will clean up references
         delete partialdec;
      }
   }
}


void GCGconshdlrDecompDeregisterPartialdec(
   SCIP* scip,
   PARTIALDECOMP* partialdec
   )
{
   int i = 0;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   assert(conshdlrdata != NULL);
   assert(partialdec != NULL);

   int id = partialdec->getID();

   // remove partialdec from list
   for(i = (int) conshdlrdata->partialdecs->size() - 1; i >= 0 ; i--)
   {
      // as registering checks for dublicates,
      // assumption that the partialdecs are unique in the list
      if(partialdec == conshdlrdata->partialdecs->at(i))
      {
         // delete item at index i
         conshdlrdata->partialdecs->erase(conshdlrdata->partialdecs->begin() + i);
         break;
      }
   }

   conshdlrdata->partialdecsbyid->erase(id);

   // remove partialdec from ancestors of all other partialdecs to avoid NULL pointers
   for(i = 0; i < (int) conshdlrdata->partialdecs->size(); i++)
   {
      conshdlrdata->partialdecs->at(i)->removeAncestorID(id);
   }
}


void GCGconshdlrDecompFreeDetprobdata(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   /* remove the presolved detprobdata */
   delete conshdlrdata->detprobdatapres;
   conshdlrdata->detprobdatapres = NULL;

   /* if parameter is set, free orig detprobdata */
   if( conshdlrdata->freeorig )
   {
      delete conshdlrdata->detprobdataorig;
      conshdlrdata->detprobdataorig = NULL;
      conshdlrdata->hasrunoriginal = FALSE;
   }
}


void GCGconshdlrDecompFreeOrigOnExit(
   SCIP* scip,
   SCIP_Bool free
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   conshdlrdata->freeorig = free;
}


int GCGconshdlrDecompGetBlockNumberCandidate(
   SCIP*                 scip,
   int                   index
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   assert(index < (int) conshdlrdata->userblocknrcandidates->size());
   return conshdlrdata->userblocknrcandidates->at(index);
}


SCIP_Real GCGconshdlrDecompGetCompleteDetectionTime(
   SCIP*                 scip
   )
{
   SCIP_Real totaltime;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   totaltime = SCIPgetClockTime(scip, conshdlrdata->completedetectionclock );

   return totaltime;
}


GCG_DECOMP** GCGconshdlrDecompGetDecomps(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);
   int decompcounter = 0;
   int i;

   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      GCGdecompFree(scip, &conshdlrdata->decomps[conshdlrdata->ndecomps - i - 1]);
   }

   SCIPfreeBlockMemoryArray(scip, &conshdlrdata->decomps, conshdlrdata->ndecomps);

   SCIP_CALL_ABORT( SCIPallocBlockMemoryArray(scip, &conshdlrdata->decomps, GCGconshdlrDecompGetNDecomps(scip) ) );

   conshdlrdata->ndecomps = GCGconshdlrDecompGetNDecomps(scip);

   sortPartialdecs(scip);
   for(i = 0; i < (int) conshdlrdata->partialdecs->size(); i++)
   {
      gcg::PARTIALDECOMP* partialdec = conshdlrdata->partialdecs->at(i);
      createDecompFromPartialdec(scip, partialdec, &conshdlrdata->decomps[decompcounter] );
      ++decompcounter;
   }

   return conshdlrdata->decomps;
}


std::string GCGconshdlrDecompGetDetectorHistoryByPartialdecId(
   SCIP* scip,
   int id
   )
{
   char buffer[SCIP_MAXSTRLEN];
   /* get partialdec and returns its detector history */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   assert(partialdec != NULL);
   partialdec->buildDecChainString(buffer);
   return std::string(buffer);
}


GCG_DETECTOR** GCGconshdlrDecompGetDetectors(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->detectors;
}


GCG_SCORE** GCGconshdlrDecompGetScores(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->scores;

}


GCG_CONSCLASSIFIER** GCGconshdlrDecompGetConsClassifiers(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->consclassifiers;
}


GCG_VARCLASSIFIER** GCGconshdlrDecompGetVarClassifiers(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->varclassifiers;
}


DETPROBDATA* GCGconshdlrDecompGetDetprobdataOrig(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( !GCGconshdlrDecompOrigDetprobdataExists(scip) )
      resetDetprobdata(scip, true);

   return conshdlrdata->detprobdataorig;
}


DETPROBDATA* GCGconshdlrDecompGetDetprobdataPresolved(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( !GCGconshdlrDecompPresolvedDetprobdataExists(scip) )
      resetDetprobdata(scip, false);

   return conshdlrdata->detprobdatapres;
}


SCIP_RETCODE GCGconshdlrDecompGetFinishedPartialdecsList(
   SCIP*          scip,
   int**          idlist,
   int*           listlength
   )
{
   std::vector<PARTIALDECOMP*> partialdecs;
   getFinishedPartialdecs(scip, partialdecs);
   partialdecVecToIdArray(partialdecs, idlist, listlength);

   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompGetPartialdecsList(
        SCIP*          scip,
        int**          idlist,
        int*           listlength
)
{
    SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
    assert( conshdlrdata != NULL );
    partialdecVecToIdArray(*conshdlrdata->partialdecs, idlist, listlength);

    return SCIP_OKAY;
}


int GCGconshdlrDecompGetNBlockNumberCandidates(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return (int) conshdlrdata->userblocknrcandidates->size();
}


int GCGconshdlrDecompGetNBlocksByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its blocks */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNBlocks();
}


int GCGconshdlrDecompGetNDecomps(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   int ndecomps;

   assert(conshdlrdata != NULL);

   ndecomps = 0;

   for(auto partialdec : *conshdlrdata->partialdecs)
   {
      if(partialdec->isComplete())
         ndecomps++;
   }

   return ndecomps;
}


int GCGconshdlrDecompGetNDetectors(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->ndetectors;
}

int GCGconshdlrDecompGetNConsClassifiers(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nconsclassifiers;
}

int GCGconshdlrDecompGetNVarClassifiers(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nvarclassifiers;
}

int GCGconshdlrDecompGetNScores(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->nscores;
}


int GCGconshdlrDecompGetNextPartialdecID(
   SCIP* scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   conshdlrdata->partialdeccounter = conshdlrdata->partialdeccounter + 1;

   assert(GCGconshdlrDecompGetPartialdecFromID(scip, conshdlrdata->partialdeccounter) == NULL);

   return conshdlrdata->partialdeccounter;
}


int GCGconshdlrDecompGetNFormerDetectionConssForID(
   SCIP*                 scip,
   int                   id
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   // look for the partialdec and check whether it is complete
   auto itr = conshdlrdata->partialdecsbyid->find(id);
   if( itr != conshdlrdata->partialdecsbyid->end() &&  itr->second->isComplete() )
   {
      return itr->second->getDetprobdata()->getNConss();
   }

   /* partialdec is not found hence we should not trust the isomorph information from detection */
   return -1;
}


int GCGconshdlrDecompGetNLinkingVarsByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its linking vars */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNLinkingvars();
}


int GCGconshdlrDecompGetNMasterConssByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its master conss */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNMasterconss();
}


int GCGconshdlrDecompGetNMasterVarsByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its master vars */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNMastervars();
}


int GCGconshdlrDecompGetNOpenConssByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its open conss */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNOpenconss();
}


int GCGconshdlrDecompGetNOpenVarsByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its open vars */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNOpenvars();
}


unsigned int GCGconshdlrDecompGetNFinishedPartialdecsOrig(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
      return 0;

   return (int) conshdlrdata->detprobdataorig->getNFinishedPartialdecs();
}


unsigned int GCGconshdlrDecompGetNFinishedPartialdecsTransformed(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdatapres == NULL )
      return 0;

   return (int) conshdlrdata->detprobdatapres->getNFinishedPartialdecs();
}


unsigned int GCGconshdlrDecompGetNOpenPartialdecsOrig(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
      return 0;

   return (int) conshdlrdata->detprobdataorig->getNOpenPartialdecs();
}


unsigned int GCGconshdlrDecompGetNOpenPartialdecsTransformed(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdatapres == NULL )
      return 0;

   return (int) conshdlrdata->detprobdatapres->getNOpenPartialdecs();
}


unsigned int GCGconshdlrDecompGetNPartialdecs(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return (int) conshdlrdata->partialdecs->size();
}


unsigned int GCGconshdlrDecompGetNPartialdecsOrig(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
      return 0;

   return (int) conshdlrdata->detprobdataorig->getNPartialdecs();
}


unsigned int GCGconshdlrDecompGetNPartialdecsTransformed(
   SCIP*       scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdatapres == NULL )
      return 0;

   return (int) conshdlrdata->detprobdatapres->getNPartialdecs();
}


int GCGconshdlrDecompGetNStairlinkingVarsByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns the number of its stairlinking vars */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->getNTotalStairlinkingvars();
}


std::vector<PARTIALDECOMP*>* GCGconshdlrDecompGetPartialdecs(
   SCIP*          scip  /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return conshdlrdata->partialdecs;
}


SCIP_RETCODE GCGconshdlrDecompGetPartialdecFromID(
   SCIP*          scip,
   int            partialdecid,
   PARTIALDECOMP_WRAPPER* pwr
   )
{
   PARTIALDECOMP* s;

   s = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);
   pwr->partialdec = s;

   return SCIP_OKAY;
}


SCIP_Real GCGconshdlrDecompGetScoreByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns its score in respect to the current score type */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);

   return partialdec->getScore(GCGgetCurrentScore(scip));
}


SCIP_RETCODE GCGconshdlrDecompGetSelectedPartialdecs(
   SCIP*          scip,
   int**          idlist,
   int*           listlength
   )
{
   /* get list of selected partialdecs */
   std::vector<PARTIALDECOMP*> selectedpartialdecs;
   getSelectedPartialdecs(scip, selectedpartialdecs);
   /* set the length of the pointer array to the list size */
   *listlength = (int) selectedpartialdecs.size();

   /* build a pointer list of ids from the partialdec list */
   for(int i = 0; i < (int) selectedpartialdecs.size(); i++)
   {
      (*idlist)[i] = selectedpartialdecs[i]->getID();
   }

   return SCIP_OKAY;
}


int GCGconshdlrDecompIncreaseNCallsCreateDecomp(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   ++conshdlrdata->ncallscreatedecomp;

   return conshdlrdata->ncallscreatedecomp;
}


SCIP_Bool GCGconshdlrDecompIsPresolvedByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns whether it is presolved */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return !partialdec->isAssignedToOrigProb();
}


SCIP_Bool GCGconshdlrDecompIsSelectedByPartialdecId(
   SCIP* scip,
   int id
   )
{
   /* get partialdec and returns whether it is currently selected */
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, id);
   return partialdec->isSelected();
}


SCIP_Bool GCGconshdlrDecompOrigDetprobdataExists(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);
   
   if(conshdlrdata->detprobdataorig == NULL)
      return FALSE;
   
   return TRUE;
}


SCIP_Bool GCGconshdlrDecompOrigPartialdecExists(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
      return FALSE;

   return conshdlrdata->detprobdataorig->getNPartialdecs() > 0;
}


SCIP_Bool GCGconshdlrDecompPresolvedDetprobdataExists(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);
   
   if(conshdlrdata->detprobdatapres == NULL)
      return FALSE;
   
   return TRUE;
}


SCIP_RETCODE GCGconshdlrDecompPrintDetectorStatistics(
   SCIP*                 scip,
   FILE*                 file
   )
{
   int i;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Detector statistics:       time     #decompositions   #complete decompositions\n");
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file,
                            "  %-17.16s:   %8.2f          %10d                 %10d\n",
                            conshdlrdata->detectors[i]->name, conshdlrdata->detectors[i]->dectime,
                            conshdlrdata->detectors[i]->ndecomps, conshdlrdata->detectors[i]->ncompletedecomps );
   }
   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompPrintScoreStatistics(
   SCIP*                 scip,
   FILE*                 file
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Score statistics:\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  Total Time       :   %8.2f\n",
                         GCGconshdlrDecompGetScoreTotalTime(scip));
   return SCIP_OKAY;
}


void GCGconshdlrDecompRegisterPartialdec(
   SCIP* scip,
   PARTIALDECOMP* partialdec
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);

   // do not register a partialdecomp multiple times
   if( conshdlrdata->partialdecsbyid->find(partialdec->getID()) == conshdlrdata->partialdecsbyid->end())
   {
      conshdlrdata->partialdecs->push_back(partialdec);
      conshdlrdata->partialdecsbyid->emplace(partialdec->getID(), partialdec);
   }
}


SCIP_RETCODE GCGconshdlrDecompSelectPartialdec(
   SCIP* scip,          /**< SCIP data structure */
   int partialdecid,    /**< id of partialdecomp */
   SCIP_Bool select     /**< select/unselect */
   )
{
   PARTIALDECOMP* partialdec = GCGconshdlrDecompGetPartialdecFromID(scip, partialdecid);
   if( partialdec )
      partialdec->setSelected(select);
   else
      return SCIP_INVALIDDATA;
   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompSetDetection(
   SCIP*                 scip,
   SCIP_PARAMSETTING     paramsetting,
   SCIP_Bool             quiet 
)
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
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

SCIP_RETCODE GCGconshdlrDecompTranslateNBestOrigPartialdecs(
   SCIP*                 scip,
   int                   n,
   SCIP_Bool             completeGreedily
)
{
   std::vector<std::pair<PARTIALDECOMP*, SCIP_Real> > candidates;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
   {
      resetDetprobdata(scip, TRUE);
      resetDetprobdata(scip, FALSE);
      return SCIP_OKAY;
   }

   if( conshdlrdata->detprobdatapres == NULL )
      resetDetprobdata(scip, FALSE);

   if( conshdlrdata->detprobdataorig->getNOpenPartialdecs() == 0 && conshdlrdata->detprobdataorig->getNFinishedPartialdecs() == 0 )
   {
      return SCIP_OKAY;
   }

   GCGconshdlrDecompChooseCandidatesFromSelected(scip, candidates, TRUE, TRUE);
   if ( !candidates.empty() )
   {
      n = MIN(n, (int)candidates.size());

      std::vector<PARTIALDECOMP *> origpartialdecs(n);
      for( int i = 0; i < n; ++i )
         origpartialdecs[i] = candidates[i].first;

      std::vector<PARTIALDECOMP *> partialdecstranslated = conshdlrdata->detprobdatapres->translatePartialdecs(
         conshdlrdata->detprobdataorig, origpartialdecs);

      if( !partialdecstranslated.empty())
      {
         PARTIALDECOMP *newpartialdec = partialdecstranslated[0];
         if( completeGreedily && !newpartialdec->isComplete() )
            newpartialdec->completeGreedily();
         SCIP_CALL(addPartialdec(scip, newpartialdec));
      }
   }

   return SCIP_OKAY;
}


SCIP_RETCODE GCGconshdlrDecompTranslateOrigPartialdecs(
   SCIP*                 scip
   )
{
   std::vector<PARTIALDECOMP*>::iterator partialdeciter;
   std::vector<PARTIALDECOMP*>::iterator partialdeciterend;

   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->detprobdataorig == NULL )
   {
      resetDetprobdata(scip, TRUE);
      resetDetprobdata(scip, FALSE);
      return SCIP_OKAY;
   }

   if( conshdlrdata->detprobdatapres == NULL )
      resetDetprobdata(scip, FALSE);

   if( conshdlrdata->detprobdataorig->getNOpenPartialdecs() == 0 && conshdlrdata->detprobdataorig->getNFinishedPartialdecs() == 0 )
   {
       return SCIP_OKAY;
   }

   std::vector<PARTIALDECOMP*> partialdecstranslated = conshdlrdata->detprobdatapres->translatePartialdecs(conshdlrdata->detprobdataorig);

   partialdeciter = partialdecstranslated.begin();
   partialdeciterend = partialdecstranslated.end();

   for(; partialdeciter != partialdeciterend; ++partialdeciter )
   {
      SCIP_CALL(addPartialdec(scip, *partialdeciter) );
   }

   return SCIP_OKAY;
}


SCIP_Bool GCGdetectionTookPlace(
   SCIP*  scip,
   SCIP_Bool original
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
   assert(conshdlrdata != NULL);

   return original == TRUE ? conshdlrdata->hasrunoriginal : conshdlrdata->hasrun;
}


SCIP_RETCODE SCIPconshdlrDecompRepairConsNames(
   SCIP*                 scip
   )
{
   long int startcount;
   SCIP_CONS** conss;
   SCIP_CONSHDLRDATA* conshdlrdata = getConshdlrdata(scip);
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


SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip
   )
{
   SCIP_CONSHDLR* conshdlr = NULL;
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;

   /* create decomp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);

   /* initialize conshdlrdata */
   conshdlrdata->partialdecs = new std::vector<PARTIALDECOMP*>();
   conshdlrdata->partialdecsbyid = new std::unordered_map<int, PARTIALDECOMP*>();
   conshdlrdata->partialdeccounter = 0;
   conshdlrdata->decomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->detectors = NULL;
   conshdlrdata->priorities = NULL;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->propagatingdetectors = NULL;
   conshdlrdata->npropagatingdetectors = 0;
   conshdlrdata->finishingdetectors = NULL;
   conshdlrdata->nfinishingdetectors = 0;
   conshdlrdata->postprocessingdetectors = NULL;
   conshdlrdata->npostprocessingdetectors = 0;

   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPcreateClock(scip, &conshdlrdata->completedetectionclock) );

   conshdlrdata->hasrun = FALSE;
   conshdlrdata->hasrunoriginal = FALSE;
   conshdlrdata->ncallscreatedecomp = 0;
   conshdlrdata->detprobdatapres = NULL;
   conshdlrdata->detprobdataorig = NULL;
   conshdlrdata->partialdectowrite = NULL;
   conshdlrdata->consnamesalreadyrepaired = FALSE;
   conshdlrdata->userblocknrcandidates = new std::vector<int>();
   conshdlrdata->freeorig = TRUE;

   /* classifiers */
   conshdlrdata->consclassifiers = NULL;
   conshdlrdata->nconsclassifiers = 0;
   conshdlrdata->consclassifierpriorities = NULL;
   conshdlrdata->varclassifiers = NULL;
   conshdlrdata->nvarclassifiers = 0;
   conshdlrdata->varclassifierpriorities = NULL;

   /* score data */
   conshdlrdata->scores = NULL;
   conshdlrdata->currscore = NULL;
   conshdlrdata->nscores = 0;
   conshdlrdata->scoreclock = NULL;
   SCIPcreateClock(scip, &conshdlrdata->scoreclock);

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

   /* add menu parameters for detection */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/enabled", "Enables detection", &conshdlrdata->enabled, FALSE, DEFAULT_ENABLED, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/postprocess", "Enables postprocessing of complete decompositions", 
                              &conshdlrdata->postprocess, FALSE, DEFAULT_POSTPROCESS, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxrounds",
                              "Maximum number of detection loop rounds", &conshdlrdata->maxndetectionrounds, FALSE,
                              DEFAULT_MAXDETECTIONROUNDS, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxtime",
                              "Maximum detection time in seconds", &conshdlrdata->maxdetectiontime, FALSE,
                              DEFAULT_MAXDETECTIONTIME, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/origprob/enabled", "Enables detection for the original problem", 
                              &conshdlrdata->enableorigdetection, FALSE, DEFAULT_ENABLEORIGDETECTION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/origprob/weightinggpresolvedoriginaldecomps",
                              "Weighting method when comparing decompositions for presolved and orig problem", &conshdlrdata->weightinggpresolvedoriginaldecomps, TRUE,
                              NO_MODIF, 0, 3, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/onlycontsubpr", "If enabled only decomposition with only continiuous variables in the subproblems are searched",
                              &conshdlrdata->bendersonlycontsubpr, FALSE, DEFAULT_BENDERSONLYCONTSUBPR, NULL, NULL) );
   
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/onlybinmaster", "If enabled only decomposition with only binary variables in the master are searched",
                              &conshdlrdata->bendersonlybinmaster, FALSE, DEFAULT_BENDERSONLYBINMASTER, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/benders/enabled", "Enables benders detection", 
                              &conshdlrdata->detectbenders, FALSE, DEFAULT_DETECTBENDERS, NULL, NULL) );

   /* classification */
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/classification/enabled", "Enables classification", &conshdlrdata->classify, FALSE, DEFAULT_CLASSIFY, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/classification/allowduplicates", "If enabled partition duplicates are allowed (for statistical reasons)",
                               &conshdlrdata->allowpartitionduplicates, FALSE, DEFAULT_ALLOWPARTITIONDUPLICATES, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/origprob/classificationenabled", "Enables classification for the original problem",
                              &conshdlrdata->enableorigclassification, FALSE, DEFAULT_ENABLEORIGCLASSIFICATION, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/classification/maxnclassesperpartition",
                              "Maximum number of classes per partition", &conshdlrdata->maxnclassesperpartition, FALSE,
                              DEFAULT_MAXNCLASSES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddIntParam(scip, "detection/classification/maxnclassesperpartitionforlargeprobs",
                              "Maximum number of classes per partition for large problems (nconss + nvars >= 50000)", &conshdlrdata->maxnclassesperpartitionforlargeprobs, FALSE,
                              DEFAULT_MAXNCLASSESLARGEPROBS, 0, INT_MAX, NULL, NULL) );

   /* block numbers */
   SCIP_CALL( SCIPaddIntParam(scip, "detection/blocknrcandidates/maxnclasses",
                              "Maximum number of classes a partition can use for voting nblockcandidates", &conshdlrdata->maxnclassesfornblockcandidates, FALSE,
                              DEFAULT_MAXNCLASSESFORNBLOCKCANDIDATES, 0, INT_MAX, NULL, NULL) );

   SCIP_CALL( SCIPaddBoolParam(scip, "detection/blocknrcandidates/medianvarspercons", 
                              "Enables the use of medianvarspercons calculation for block number candidates calculation ",
                              &conshdlrdata->blocknumbercandsmedianvarspercons,
                              FALSE, DEFAULT_BLOCKNUMBERCANDSMEDIANVARSPERCONS, NULL, NULL) );

   /* scores */
   SCIP_CALL( SCIPaddStringParam(scip, "detection/scores/selected",
                                 "Sets the score calculation for comparing (partial) decompositions (use score shortname)",
                                 &conshdlrdata->currscore,
                                 FALSE, DEFAULT_SCORE, paramChgdScore, NULL) );

   return SCIP_OKAY;
}
