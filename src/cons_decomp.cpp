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

/**@file   cons_decomp.c
 * @ingroup CONSHDLRS
 * @brief  constraint handler for structure detection
 * @author Martin Bergner
 * @author Christian Puchert
 * @author Michael Bastubbe
 *
 * This constraint handler will run all registered structure detectors in
 * in an iterative scheme increasing priority until the first detector finds a suitable structure.
 *
 * It also manages the explore menu.
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


   SeeedPtr              curruserseeed;                           /**< help pointer for reader and toolbox to iteratively build (partial) decomposition */
   SeeedPtr              lastuserseeed;                           /**< help pointer for toolbox to revoke last changes to curruserseeed */

   SCIP_Bool             consnamesalreadyrepaired;                /**< stores whether or not    */

   SCIP_Bool             unpresolveduserseeedadded;               /**< stores whether or not an unpresolved user seeed was added */

   /* new data fields for selection/exploration management */
   int                    startidvisu;                            /**< when displaying the list of decomps, this is the starting index */
   int                    selectvisulength;                       /**< number of decompositions to be displayed at once */
   std::vector<SeeedPtr>* listall;                                /**< vector containing the current list of decomps (to visualize, write, consider for family tree, consider for solving etc. )*/
   std::vector<int>*      selected;                               /**< vector containing the indices of selected decompositions */
   SCIP_Bool              selectedexists;                         /**< are there some selected decompositions */
   int                    seeedcounter;                           /**< counts the number of seeeds, used for seeed ids */
   std::vector<std::pair<SeeedPtr, SCIP_Real> >* candidates;      /**< vector containing the pairs of candidate list of decomps (to visualize, write, consider for family tree, consider for solving etc.) sorted according to  */
   int                    currscoretype;                          /**< indicates which score should be used for comparing (partial) decompositions (
                                                                          0:max white,
                                                                          1: border area,
                                                                          2:classic,
                                                                          3:max foreseeing white,
                                                                          4: ppc-max-white,
                                                                          5:max foreseeing white with aggregation info,
                                                                          6: ppc-max-white with aggregation info,
                                                                          7: experimental benders score */

   SCIP_Bool               nonfinalfreetransform;                 /**< help bool to notify a nonfinal free transform (needed if presolving is revoked, e.g. if unpresolved decomposition is used, and transformation is not successful) */
   std::vector<int>*       userblocknrcandidates;                 /**< vector to store block number candidates that were given by user */
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


/**
 * method to unselect all decompositions, called in consexit, and when the seeedlist is updated
 * (especially if new (partial) are added )
 *
 *@returns SCIP status
 */
static
SCIP_RETCODE SCIPconshdlrdataDecompUnselectAll(
   SCIP*          scip  /**< SCIP data structure */
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

   std::vector<int>::const_iterator selectediter = conshdlrdata->selected->begin();
   std::vector<int>::const_iterator selectediterend = conshdlrdata->selected->end();

   for( ; selectediter != selectediterend; ++selectediter )
   {
      conshdlrdata->listall->at(*selectediter)->setSelected(false);
   }

   conshdlrdata->selected->clear();

   conshdlrdata->selectedexists = FALSE;

   return SCIP_OKAY;
}

/**
 * Gets the currently selected scoretype
 * @returns the currently selected scoretype
 */
static
SCORETYPE SCIPconshdlrdataGetScoretype(
   SCIP_CONSHDLRDATA* conshdlrdata  /**< conshdlr data */
   )
{
   return  static_cast<scoretype>(conshdlrdata->currscoretype);
}


/**
 * Gets the shortname of the given scoretype
 *
 * @returns the shortname of the given Scoretype
 */
static
char*  SCIPconshdlrDecompGetScoretypeShortName(
   SCIP*       scip,    /**< SCIP data structure */
   SCORETYPE   sctype   /**< scoretype */
   )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;
   /* set detector chain info string */
   SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "") ;


   if( sctype == scoretype::MAX_WHITE )
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maxwhi") ;

   if( sctype == scoretype::CLASSIC )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classi") ;

   if( sctype == scoretype::BORDER_AREA )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "border") ;

   if( sctype == scoretype::MAX_FORESSEEING_WHITE )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "forswh") ;

   if( sctype == scoretype::MAX_FORESEEING_AGG_WHITE )
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "fawh") ;


   if( sctype == scoretype::SETPART_FWHITE )
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfwh ") ;

   if( sctype == scoretype::SETPART_AGG_FWHITE )
           SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "spfawh") ;


   if( sctype == scoretype::BENDERS )
              SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "bender") ;


   SCIP_CALL_ABORT ( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN ) );

   return copy;

}

/*!
 * returns the description of the given scoretype
 *
 * @returns description of the scoretype
 */
static
char*  SCIPconshdlrDecompGetScoretypeDescription(
   SCIP*       scip,    /**< SCIP data structure */
   SCORETYPE   sctype   /**< scoretype */
      )
{
   char scoretypename[SCIP_MAXSTRLEN];
   char* copy;

   /* set detector chain info string */
   SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "") ;

   if( sctype == scoretype::MAX_WHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum white area score (i.e. maximize fraction of white area score; white area is nonblock and nonborder area, stairlinking variables count as linking)") ;

   if( sctype == scoretype::CLASSIC)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "classical score") ;

   if( sctype == scoretype::BORDER_AREA)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "minimum border score (i.e. minimizes fraction of border area score; )")  ;

   if( sctype == scoretype::MAX_FORESSEEING_WHITE)
         SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing  white area score (i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking)")  ;

   if( sctype == scoretype::MAX_FORESEEING_AGG_WHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "maximum foreseeing  white area score with aggregation information(i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking)")  ;


   if( sctype == scoretype::SETPART_FWHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing  white area score (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )")  ;

   if( sctype == scoretype::SETPART_AGG_FWHITE)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "setpartitioning maximum foreseeing white area score with aggregation information (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )")  ;

   if( sctype == scoretype::BENDERS)
      SCIPsnprintf( scoretypename, SCIP_MAXSTRLEN, "experimental score to evaluate benders decompositions")  ;


   SCIP_CALL_ABORT ( SCIPduplicateBlockMemoryArray(scip, &copy, scoretypename, SCIP_MAXSTRLEN ) );

   return copy;

}


/** local method to handle store a complete seeed in the unpresolved seeedpool
 *
 * @returns SCIP status */
static
SCIP_RETCODE  SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(
     SCIP* scip,     /**< SCIP data structure */
     SeeedPtr  seeed /**< pointer to seeed */
   ){

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
        SCIPinfoMessage(scip, NULL, "Decomposition to add is already known to gcg!\n");

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
        SCIPinfoMessage(scip, NULL, "Decomposition to add is already known to gcg!\n");

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
        SCIPinfoMessage(scip, NULL, "Decomposition to add is already known to gcg!\n");

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
        SCIPinfoMessage(scip, NULL, "Decomposition to add is already known to gcg!\n");

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
   conshdlrdata->listall->clear();


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

   delete conshdlrdata->selected;
   delete conshdlrdata->listall;
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

   conshdlrdata->curruserseeed = NULL;
   conshdlrdata->lastuserseeed = NULL;
   conshdlrdata->unpresolveduserseeedadded = FALSE;
   conshdlrdata->startidvisu = 0;
   conshdlrdata->selectvisulength = 10;
   conshdlrdata->listall = new std::vector<SeeedPtr>(0, NULL);
   conshdlrdata->selected = new std::vector<int>(0, -1);
   conshdlrdata->candidates = new std::vector<std::pair<SeeedPtr, SCIP_Real > >(0);
   conshdlrdata->selectedexists = FALSE;
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
         "indicates which score should be used for comparing (partial) decompositions (0:max white, 1: border area, 2:classic, 3:max foreseeing white, 4: ppc-max-white, 5:max foreseeing white with aggregation info, 6: ppc-max-white with aggregation info, 7: experimental benders score): ", &conshdlrdata->currscoretype, FALSE,
         scoretype::SETPART_FWHITE, 0, 7, NULL, NULL) );


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


/** Shows header for seeed information in explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompShowListExtractHeader(
   SCIP*                   scip  /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   int ndetectedpresolved;
   int ndetectedunpresolved;
   int nuserpresolvedfull;
   int nuserpresolvedpartial;
   int nuserunpresolvedfull;
   int nuserunpresolvedpartial;

   char* scorename;

   size_t i;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );

   ndetectedpresolved = 0;
   ndetectedunpresolved = 0;
   nuserpresolvedfull = 0;
   nuserpresolvedpartial = 0;
   nuserunpresolvedfull = 0;
   nuserunpresolvedpartial = 0;

   /* count corresponding seeeds */
   for ( i = 0; i < conshdlrdata->listall->size(); ++i )
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->listall->at(i);
      if( seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::NOT && !seeed->isFromUnpresolved() )
         ++ndetectedpresolved;
      if( seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::NOT && seeed->isFromUnpresolved() )
         ++ndetectedunpresolved;
      if( seeed->isComplete() && ( seeed->getUsergiven() == gcg::USERGIVEN::COMPLETE || seeed->getUsergiven() == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && !seeed->isFromUnpresolved() )
         ++nuserpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::PARTIAL && !seeed->isFromUnpresolved() )
         ++nuserpresolvedpartial;
      if( seeed->isComplete() && ( seeed->getUsergiven() == gcg::USERGIVEN::COMPLETE || seeed->getUsergiven() == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && seeed->isFromUnpresolved() )
         ++nuserunpresolvedfull;
      if( !seeed->isComplete() && seeed->getUsergiven() == gcg::USERGIVEN::PARTIAL && seeed->isFromUnpresolved() )
         ++nuserunpresolvedpartial;

   }

   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, "============================================================================================= ");
   SCIPdialogMessage(scip, NULL, "\n");
   SCIPdialogMessage(scip, NULL, "Summary              presolved       original \n");
   SCIPdialogMessage(scip, NULL, "                     ---------       -------- \n");
   SCIPdialogMessage(scip, NULL, "detected             ");
   SCIPdialogMessage(scip, NULL, "%9d       ", ndetectedpresolved );
   SCIPdialogMessage(scip, NULL, "%8d\n", ndetectedunpresolved );
   SCIPdialogMessage(scip, NULL, "user given (partial) ");
   SCIPdialogMessage(scip, NULL, "%9d       ", nuserpresolvedpartial );
   SCIPdialogMessage(scip, NULL, "%8d\n", nuserunpresolvedpartial );
   SCIPdialogMessage(scip, NULL, "user given (full)    ");
   SCIPdialogMessage(scip, NULL, "%9d       ", nuserpresolvedfull );
   SCIPdialogMessage(scip, NULL, "%8d\n", nuserunpresolvedfull );

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");
   SCIPdialogMessage(scip, NULL, "   id   nbloc  nmacon  nlivar  nmavar  nstlva  %.6s  history  pre  nopcon  nopvar  usr  sel \n", scorename );
   SCIPdialogMessage(scip, NULL, " ----   -----  ------  ------  ------  ------  ------  -------  ---  ------  ------  ---  --- \n");

   SCIPfreeBlockMemoryArrayNull(scip, &scorename, SCIP_MAXSTRLEN);


   return SCIP_OKAY;
}

/** Shows information about the current user seeed in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompShowCurrUserSeeedInfo(
   SCIP*                   scip     /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if ( conshdlrdata->curruserseeed->isFromUnpresolved() )
      conshdlrdata->curruserseeed->displaySeeed();
   else
      conshdlrdata->curruserseeed->displaySeeed();


   return SCIP_OKAY;
}

/** Shows detailed information about seeeds in explore menu
 *
 *@returns SCIP status
 * */
static
SCIP_RETCODE SCIPconshdlrDecompShowListExtract(
   SCIP*                 scip               /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   size_t i;

   for( i = conshdlrdata->startidvisu; i < (size_t) conshdlrdata->startidvisu + (size_t) conshdlrdata->selectvisulength && i < conshdlrdata->listall->size(); ++i)
   {
      SeeedPtr seeed;

      seeed = conshdlrdata->listall->at(i);

      assert( seeed->checkConsistency( ) );

      SCIPdialogMessage(scip, NULL, " %4d   ", i );
      SCIPdialogMessage(scip, NULL, "%5d  ", seeed->getNBlocks() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMasterconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNLinkingvars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMastervars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNTotalStairlinkingvars() );
      if( seeed->isComplete() )
         SCIPdialogMessage(scip, NULL, "%.4f  ",  seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) );
      else
         SCIPdialogMessage(scip, NULL, "<=%.2f  ", seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) );
      SCIPdialogMessage(scip, NULL, "%7s  ", seeed->getDetectorChainString() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->isFromUnpresolved() ? "no" : "yes")  );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenvars() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->getUsergiven() == gcg::USERGIVEN::NOT ? "no" : "yes")   );
      SCIPdialogMessage(scip, NULL, "%3s  \n", (seeed->isSelected() ? "yes" : "no")  );
   }

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

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


/** Shows information about the explore screen and its abbreviations
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompShowLegend(
   SCIP* scip  /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   char * scorename;
   char * scoredescr;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   scorename = SCIPconshdlrDecompGetScoretypeShortName(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );
   scoredescr = SCIPconshdlrDecompGetScoretypeDescription(scip, SCIPconshdlrdataGetScoretype(conshdlrdata) );


   SCIPdialogMessage(scip, NULL, "List of included detectors for decompositions histories: \n" );

   SCIPdialogMessage(scip, NULL, "\n%30s    %4s\n", "detector" , "char"  );
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "--------" , "----"  );

   for( int det = 0; det < conshdlrdata->ndetectors; ++det )
   {
      DEC_DETECTOR* detector;

      detector = conshdlrdata->detectors[det];

      SCIPdialogMessage(scip, NULL, "%30s    %4c\n", DECdetectorGetName(detector), DECdetectorGetChar(detector)  );
   }
   SCIPdialogMessage(scip, NULL, "%30s    %4s\n", "given by user" , "U"  );

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "List of abbreviations of decomposition table \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "abbreviation", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "id", "id of the decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nbloc", "number of blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmacon", "number of master constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nlivar", "number of linking variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nmavar", "number of master variables (do not occur in blocks)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nstlva", "number of stairlinking variables (disjoint from linking variables)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", scorename, scoredescr);
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "history", "list of detector chars worked on this decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "pre", "is this decomposition for the presolved problem");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopcon", "number of open constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopvar", "number of open variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "usr", "was this decomposition given by the user");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sel", "is this decomposition selected at the moment");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

   SCIPfreeBlockMemoryArrayNull(scip, &scorename, SCIP_MAXSTRLEN);
   SCIPfreeBlockMemoryArrayNull(scip, &scoredescr, SCIP_MAXSTRLEN);

   return SCIP_OKAY;
}


/** Shows help for the user toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompShowToolboxInfo(
   SCIP* scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "Options to proceed: \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "option", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "conss", "assign unassigned constraints to master/blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "vars", "assign unassigned variables to master(only)/linking/blocks");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "refine", "refine implicit constraint and variables assignments");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "finish", "choose a finishing detector that completes the decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "quit the modification process and returns to main menu");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "undo", "last modification is undone (atm only the last modification can be undone)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "shows a visualization of the current decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "propagate", "list all detectors that can propagate the current seeed and apply one to propagate it");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "finish", "list all detectors that can finish the current seeed and apply one to finish it");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "postprocess", "apply postprocessing to a finished seeed by selecting a suitable postprocessor");
   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

   return SCIP_OKAY;
}


/** Modifies the number of presented seeeds in the explore menu via dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompModifyNVisualized(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int newval;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the maximum number of decompositions displayed at once in the table [%d]:\n",
      conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   newval = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      newval = atoi(ntovisualize);

   if (newval != 0)
      conshdlrdata->selectvisulength = newval;

   return SCIP_OKAY;
}


/** Shows a visualization of the seeed specified by the user via the dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompSelectVisualize(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be visualized:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = -1;
   if( commandlen != 0 )
      idtovisu = atoi(ntovisualize);

   /* check whether ID is in valid range */
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( commandlen == 0 || idtovisu < 0 || idtovisu >= (int)conshdlrdata->listall->size() )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_OKAY;
   }

   conshdlrdata->listall->at(idtovisu)->showVisualisation();

   return SCIP_OKAY;
}


/**
 * Calculates and displays the strong decomposition score for this decomposition in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPconshdlrDecompSelectCalcStrongDecompositionScore(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntocalcstrong;
   SCIP_Bool endoffile;
   int idtocalcstrong;
   int commandlen;

   assert( scip != NULL );
   conshdlr = SCIPfindConshdlr( scip, CONSHDLR_NAME );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData( conshdlr );
   assert( conshdlrdata != NULL );

   /* read the id of the decomposition to be calculate strong decomp score */
   SCIPdialogMessage( scip, NULL,
      "Please specify the id of the decomposition that should be evaluated by strong decomposition score:\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntocalcstrong, &endoffile ) );
   commandlen = strlen( ntocalcstrong );

   idtocalcstrong = -1;
   if( commandlen != 0 )
   {
      std::stringstream convert( ntocalcstrong );
      convert >> idtocalcstrong;

      if ( idtocalcstrong == 0 && ntocalcstrong[0] != '0' )
      {
         idtocalcstrong = -1;
      }
   }

   /* call calculation strong decomp score method according to chosen parameters */
   if( 0 <= idtocalcstrong && idtocalcstrong < (int)conshdlrdata->listall->size() )
   {
      SCIP_Real score;
      gcg::Seeedpool* seeedpool = ( conshdlrdata->listall->at( idtocalcstrong )->isFromUnpresolved() ?
         conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );
      seeedpool->calcStrongDecompositionScore(conshdlrdata->listall->at( idtocalcstrong ), &score);
      SCIPdialogMessage( scip, NULL, "Strong decomposition score of this decomposition is %f.", score) ;
   }
   else
   {
      SCIPdialogMessage( scip, NULL, "This is not an existing id." );
   }

   return SCIP_OKAY;
}


/**
 * Displays information about a seeed that is chosen by the user in a dialog.
 *
 * @returns SCIP status
 */
static
SCIP_RETCODE SCIPconshdlrDecompSelectInspect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntoinspect;
   char* ndetaillevel;
   SCIP_Bool endoffile;
   int idtoinspect;
   int detaillevel;

   int commandlen;

   assert( scip != NULL );
   conshdlr = SCIPfindConshdlr( scip, CONSHDLR_NAME );
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData( conshdlr );
   assert( conshdlrdata != NULL );

   /* read the id of the decomposition to be inspected */
   SCIPdialogMessage( scip, NULL, "Please specify the id of the decomposition to be inspected:\n");
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ntoinspect, &endoffile ) );
   commandlen = strlen( ntoinspect );

   idtoinspect = -1;
   if( commandlen != 0 )
      idtoinspect = atoi( ntoinspect );

   /* check whether ID is in valid range */
   if( idtoinspect < 0 || idtoinspect >= (int)conshdlrdata->listall->size() )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

   /* read the desired detail level; for wrong input, it is set to 1 by default */
   SCIPdialogMessage( scip, NULL,
      "Please specify the detail level:\n  0 - brief overview\n  1 - block and detector info (default)\n  2 - cons and var assignments\n" );
   SCIP_CALL( SCIPdialoghdlrGetWord( dialoghdlr, dialog, " ", &ndetaillevel, &endoffile ) );
   commandlen = strlen( ndetaillevel );

   detaillevel = 1;
   if( commandlen != 0 )
   {
      std::stringstream convert( ndetaillevel );
      convert >> detaillevel;

      if ( detaillevel < 0 || ( detaillevel == 0 && ndetaillevel[0] != '0' ) )
      {
         detaillevel = 1;
      }
   }

   /* call displayInfo method according to chosen parameters */
   assert( 0 <= idtoinspect && idtoinspect < (int)conshdlrdata->listall->size() );
   conshdlrdata->listall->at( idtoinspect )->displayInfo( detaillevel );

   return SCIP_OKAY;
}


/** shows a visualization of current user seeed
 *
 * @returns SCip status*/
static
SCIP_RETCODE SCIPconshdlrDecompSelectVisualizeCurrentUserSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->curruserseeed->showVisualisation();

   return SCIP_OKAY;
}


/** Lets the user select a seeed to modify in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompToolboxChoose(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntochoose;
   SCIP_Bool endoffile;
   int idtochoose;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the (partial) decomposition to be chosen for modification:\n",
      conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntochoose, &endoffile) );
   commandlen = strlen(ntochoose);

   idtochoose = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtochoose = atoi(ntochoose);

   if ( commandlen == 0 || idtochoose < 0 || idtochoose >= (int)conshdlrdata->listall->size() )
   {
      SCIPdialogMessage( scip, NULL, "This id is out of range." );
      return SCIP_PARAMETERWRONGVAL;
   }

   if( conshdlrdata->curruserseeed != NULL )
      delete conshdlrdata->curruserseeed;

   conshdlrdata->curruserseeed = new gcg::Seeed( conshdlrdata->listall->at(idtochoose) );

   return SCIP_OKAY;
}

/** Lets the user select decompositions from the explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompExploreSelect(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;
   SeeedPtr toselect;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be selected:\n", conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   toselect = conshdlrdata->listall->at(idtovisu);

   toselect->setSelected(!toselect->isSelected() );

   if( !toselect->isSelected() )
   {
      conshdlrdata->selected->erase(  find( conshdlrdata->selected->begin(), conshdlrdata->selected->end(), idtovisu) );
   }
   else
   {
      std::cout << "is selected!" << toselect->isSelected() <<std::endl;
      conshdlrdata->selected->push_back(idtovisu);
      assert(toselect->isSelected());
   }

   conshdlrdata->selectedexists = (conshdlrdata->selected->size() > 0);

   return SCIP_OKAY;
}

/** Shows help section of explore menu
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompShowHelp(
   SCIP* scip  /**< SCIP data structure */
   )
{
   assert(scip != NULL);

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "List of selection commands \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "command", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "-------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "select", "selects/unselects decomposition with given id");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "modify", "modify an existing decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "create", "create a new decomposition");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "back", "displays the preceding decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "next", "displays the subsequent decompositions (if there are any)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "top", "displays the first decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "end", "displays the last decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "legend", "displays the legend for table header and history abbreviations");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "help", "displays this help");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "dispNEntries", "modifies the number of displayed decompositions ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "finishes decomposition explorer and goes back to main menu");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "experimental feature: visualizes the specified decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "inspect", "displays detailed information for the specified decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "calc_strong", "calculates and displays the strong decomposition score for this decomposition");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");

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

/*
 * @brief method too handle user input for "explore" command
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecSelect(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
   /* while user has not aborted: show current list extract */

   while ( !finished )
   {
      int commandlen;

      SCIP_CALL( SCIPconshdlrDecompShowListExtractHeader(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowListExtract(scip) );


      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Please enter command or decomposition id to select (or \"h\" for help) : \nGCG/explore> ", &command, &endoffile) );

      commandlen = strlen(command);

      if( strncmp( command, "back", commandlen) == 0 )
      {
         conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
         if(conshdlrdata->startidvisu < 0 )
            conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
         if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
            conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }

      if( strncmp( command, "quit", commandlen) == 0 )
      {
         finished = TRUE;
         SCIP_CALL(SCIPconshdlrDecompChooseCandidatesFromSelected(scip, FALSE) );
         continue;
      }

      if( strncmp( command, "legend", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompShowLegend(scip) );
         continue;
      }

      if( strncmp( command, "dispNEntries", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompModifyNVisualized(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompShowHelp(scip) );
         continue;
      }

      if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectVisualize(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "inspect", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompSelectInspect( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "calc_strong", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompSelectCalcStrongDecompositionScore( scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompExploreSelect(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "toolbox", commandlen) == 0 )
      {
	 //@todo deprecated, use create/modify instead
         SCIP_CALL( SCIPconshdlrDecompExecToolbox(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
      if( strncmp( command, "modify", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompExecToolboxModify(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
      if( strncmp( command, "create", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompExecToolboxCreate(scip, dialoghdlr, dialog) );
         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );
         continue;
      }
   }

   return SCIP_OKAY;
}


/** Lets user modify conss during modification of seeed in toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompToolboxModifyConss(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;
    SCIP_Bool         matching;
    char* consregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert( conshdlr != NULL );
    matching = FALSE;


    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    SeeedPtr seeed  = conshdlrdata->curruserseeed;
    gcg::Seeedpool* seeedpool;
    std::vector<int> matchingconss  = std::vector<int>(0);

    seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
    /* Does user want to modify existing or create a new partial decomposition ?*/
    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
       "Please specify a regular expression (modified ECMAScript regular expression grammar) matching the names of unassigned constraints you want to assign : \nGCG/toolbox> ",
       &consregex, &endoffile) );

    /* case distinction: */

    std::regex expr;
    try  {
       expr = std::regex(consregex);
    }
    catch (const std::regex_error& e) {
       std::cout << "regex_error caught: " << e.what() << '\n';
       if (e.code() == std::regex_constants::error_brack) {
          std::cout << "The code was error_brack\n";
       }
    }

    for( int oc = 0; oc < seeed->getNOpenconss(); ++oc )
    {
       const char* consname;

       consname = SCIPconsGetName(  seeedpool->getConsForIndex(seeed->getOpenconss()[oc] ) );


       if( std::regex_match(consname, expr) )
       {
          matching = TRUE;
          matchingconss.push_back(seeed->getOpenconss()[oc]);
          SCIPdebugMessage(" consname %s matches regex %s \n", consname, consregex );
       } else
          SCIPdebugMessage(" consname %s does not match regex %s \n", consname, consregex);
    }

    if( !matching )
    {
       SCIPdialogMessage(scip, NULL, " There are no unassigned constraints with names matching given regular expression. Return to toolbox main menu.\n");
       return SCIP_OKAY;
    }

    if( conshdlrdata->lastuserseeed != NULL)
       delete conshdlrdata->lastuserseeed;
    conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;


    if( matchingconss.size() > 10 )
       SCIPdebugMessage(" There are %d unassigned constraints with names matching given regular expression. Showing the first 10:\n", (int) matchingconss.size());
    else
       SCIPdebugMessage(" There are %d unassigned constraints with names matching given regular expression: \n", (int) matchingconss.size());

    for( size_t mc = 0 ; mc < 10 && mc < matchingconss.size(); ++mc )
       SCIPdialogMessage(scip, NULL, " %s \n", SCIPconsGetName( seeedpool->getConsForIndex( matchingconss[mc] ) ));

    SCIPdialogMessage(scip, NULL, "\n Should these constraints be added to: \n");
    SCIPdialogMessage(scip, NULL, " master \n");
    SCIPdialogMessage(scip, NULL, " block (to be specified) \n");
    SCIPdialogMessage(scip, NULL, " nothing (return to toolbox main menu)? \n");


    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify how to proceed: \nGCG/toolbox> ", &command, &endoffile) );

    commandlen = strlen(command);

    /* case distinction: */
    if( strncmp( command, "master", commandlen) == 0 )
    {
       for( size_t mc = 0 ;  mc < matchingconss.size(); ++mc )
       {
          seeed->bookAsMasterCons( matchingconss[mc] );
       }
    }
    else if( strncmp( command, "block", commandlen) == 0 )
    {
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify the block number these constraints should be assigned to: \nGCG/toolbox> ", &command2, &endoffile) );
       char* tail;
       int blockid = strtol(command2, &tail, 10);
       for( size_t mc = 0 ;  mc < matchingconss.size(); ++mc )
       {
          seeed->bookAsBlockCons( matchingconss[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    seeed->flushBooked();

   return SCIP_OKAY;
}

/** Lets user specify how to finish the modified seeed while using the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPconshdlrDecompToolboxModifyFinish(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         choosenfinisher;

   char* command;
   SCIP_Bool endoffile;
   char* tail;
   int finisherid;
   SEEED_PROPAGATION_DATA* seeedPropData;
   DEC_DETECTOR* finisher;
   SCIP_Result result;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SeeedPtr seeed  = conshdlrdata->curruserseeed;
   gcg::Seeedpool* seeedpool;
   std::vector<int> matchingvars  = std::vector<int>(0);

   seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   choosenfinisher = FALSE;
   while ( !choosenfinisher )
   {
       SCIPdialogMessage(scip, NULL, " Available finisher: \n");
       /* 1) print out available finisher */
       SCIPdialogMessage(scip, NULL, "%d :  %s \n", -1, "abort" );
       for( int fi = 0; fi < seeedpool->getNFinishingDetectors(); ++fi )
       {
          SCIPdialogMessage(scip, NULL, "%d :  %s \n", fi, DECdetectorGetName(seeedpool->getFinishingDetectorForIndex(fi) ) );
       }

       /* Does user want to modify existing or create a new partial decomposition ?*/
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
          "Please specify the index of the finisher to use : \nGCG/toolbox> ", &command, &endoffile) );

       finisherid = strtol(command, &tail, 10);

       if( finisherid >= seeedpool->getNFinishingDetectors() || finisherid < -1 )
       {
            SCIPdialogMessage(scip, NULL, "The specified id is invalid \n"  );
            continue;
       }
       choosenfinisher = TRUE;
   }

   seeedPropData = new SEEED_PROPAGATION_DATA();
   seeedPropData->seeedpool = seeedpool;
   seeedPropData->nNewSeeeds = 0;
   seeedPropData->seeedToPropagate = new gcg::Seeed(conshdlrdata->curruserseeed);

   if( conshdlrdata->lastuserseeed != NULL)
      delete conshdlrdata->lastuserseeed;
   conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;

   finisher = seeedpool->getFinishingDetectorForIndex(finisherid);
   finisher->finishSeeed(scip, finisher, seeedPropData, &result);

   delete conshdlrdata->curruserseeed;

   for( int i = 0; i <  seeedPropData->nNewSeeeds; ++i)
   {
      delete seeedPropData->newSeeeds[i];
   }

   delete seeedPropData->seeedToPropagate;
   delete seeedPropData;

   return SCIP_OKAY;
}

/** Lets user modify vars during use of the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompToolboxModifyVars(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{

   SCIP_CONSHDLR* conshdlr;
    SCIP_CONSHDLRDATA* conshdlrdata;
    SCIP_Bool         matching;
    char* varregex;
    char* command;
    char* command2;
    SCIP_Bool endoffile;
    int commandlen;

    assert(scip != NULL);
    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    assert( conshdlr != NULL );
    matching = FALSE;


    conshdlrdata = SCIPconshdlrGetData(conshdlr);
    assert(conshdlrdata != NULL);

    SeeedPtr seeed  = conshdlrdata->curruserseeed;
    gcg::Seeedpool* seeedpool;
    std::vector<int> matchingvars  = std::vector<int>(0);

    seeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
    /* Does user want to modify existing or create a new partial decomposition ?*/
    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
       "Please specify a regular expression (modified ECMAScript regular expression grammar) matching the names of unassigned variables you want to assign : \nGCG/toolbox> ",
       &varregex, &endoffile) );

    /* case distinction: */

    std::regex expr;
    try  {
       expr = std::regex(varregex);
    }
    catch (const std::regex_error& e) {
       std::cout << "regex_error caught: " << e.what() << '\n';
       if (e.code() == std::regex_constants::error_brack) {
          SCIPdebugMessage("The code was error_brack\n");
       }
    }

    for( int oc = 0; oc < seeed->getNOpenvars(); ++oc )
    {
       const char* varname;

       varname = SCIPvarGetName(  seeedpool->getVarForIndex(seeed->getOpenvars()[oc] ) );

       SCIPdebugMessage("check var %s for regex %s \n", varname, varregex);

       if( std::regex_match(varname, expr) )
       {
          matching = TRUE;
          matchingvars.push_back(seeed->getOpenvars()[oc]);
          SCIPdebugMessage( " varname %s matches regex %s \n", varname, varregex );
       } else
          SCIPdebugMessage(" varname %s does not match regex %s \n", varname, varregex);
    }

    if( !matching )
    {
       SCIPdialogMessage(scip, NULL,
          " There are no unassigned variables with names matching given regular expression. Return to toolbox main menu.\n");
       return SCIP_OKAY;
    }

    if( conshdlrdata->lastuserseeed != NULL)
       delete conshdlrdata->lastuserseeed;
    conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;


    if( matchingvars.size() > 10 )
       SCIPdialogMessage(scip, NULL,
          " There are %d unassigned variables with names matching given regular expression. Showing the first 10:\n",
          matchingvars.size());
    else
       SCIPdialogMessage(scip, NULL, " There are %d unassigned variables with names matching given regular expression: \n",
          matchingvars.size());

    for( size_t mc = 0 ; mc < 10 && mc < matchingvars.size(); ++mc )
       SCIPdialogMessage(scip, NULL, " %s \n", SCIPvarGetName( seeedpool->getVarForIndex( matchingvars[mc] ) ));

    SCIPdialogMessage(scip, NULL, "\n Should these variables be added to: \n");
    SCIPdialogMessage(scip, NULL, " master-only (static) \n");
    SCIPdialogMessage(scip, NULL, " linking \n");
    SCIPdialogMessage(scip, NULL, " block (to be specified) \n");
    SCIPdialogMessage(scip, NULL, " nothing (return to toolbox main menu)? \n");

    SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please specify how to proceed: \nGCG/toolbox> ", &command,
       &endoffile) );

    commandlen = strlen(command);

    /* case distinction: */
    if( strncmp( command, "master", commandlen) == 0 )
    {
       for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
       {
          seeed->bookAsMasterVar( matchingvars[mc] );
       }
    } else
       if( strncmp( command, "linking", commandlen) == 0 )
           {
              for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
              {
                 seeed->bookAsLinkingVar( matchingvars[mc] );
              }
           }
    else if( strncmp( command, "block", commandlen) == 0 )
    {
       SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
          "Please specify the block number these variables should be assigned to: \nGCG/toolbox> ", &command2, &endoffile) );
       char* tail;
       int blockid = strtol(command2, &tail, 10);
       for( size_t mc = 0 ;  mc < matchingvars.size(); ++mc )
       {
          seeed->bookAsBlockVar( matchingvars[mc], blockid );
       }
    }
    else
       return SCIP_OKAY;

    seeed->flushBooked();
    seeed->deleteEmptyBlocks(true);

   return SCIP_OKAY;
}

/** Apply propagation, finishing or postprocessing to the current user seeed via dialog
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompToolboxActOnSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog,     /**< dialog for user input management */
   toolboxtype             action      /**< what to do: can be set to PROPAGATE, FINISH or POSTPROCESS */
   )
{
   char* command;
   int commandlen;
   SCIP_Bool endoffile;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Result result;
   DEC_Detector** detectors;
   int ndetectors;
   int i, j;
   SEEED_PROPAGATION_DATA* seeedPropData;
   gcg::Seeedpool* seeedpool;
   SCIP_Bool finished, displayinfo;
   char stri[SCIP_MAXSTRLEN];
   const char* actiontype;

   /* set string for dialog */
   if( action == PROPAGATE )
     actiontype = "propagated";
   else if( action == FINISH )
      actiontype = "finished";
   else if( action == POSTPROCESS )
      actiontype = "postprocessed";
   else
      actiontype = "UNDEFINED_ACTION";

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   if( action == POSTPROCESS && conshdlrdata->curruserseeed->isComplete() == FALSE ) 
   {
      SCIPinfoMessage(scip, NULL, "The currently selected seeed is not finished, postprocessing not possible.\n");
      return SCIP_OKAY;
   }

   if( conshdlrdata->ndetectors == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector available!\n\n");
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPallocBufferArray(scip, &detectors, conshdlrdata->ndetectors) );

   /* determine the detectors that implement the specified callback */
   ndetectors = 0;
   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      if( (action == PROPAGATE && conshdlrdata->detectors[i]->propagateFromToolbox)
       || (action == FINISH && conshdlrdata->detectors[i]->finishFromToolbox)
       || (action == POSTPROCESS && conshdlrdata->detectors[i]->postprocessSeeed) )
      {
         detectors[ndetectors] = conshdlrdata->detectors[i];
         ++ndetectors;
      }
   }

   if( ndetectors == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No detector implements this callback, returning!\n\n");
      return SCIP_OKAY;
   }

   /* build seeed propagation data needed in callbacks */
   seeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;

   seeedPropData = new SEEED_PROPAGATION_DATA();
   seeedPropData->seeedpool = seeedpool;
   seeedPropData->nNewSeeeds = 0;
   seeedPropData->seeedToPropagate = new gcg::Seeed(conshdlrdata->curruserseeed);
   seeedPropData->seeedToPropagate->setSeeedpool(seeedpool);
   if( action != POSTPROCESS )
   {
      SCIP_CALL( SCIPallocMemoryArray(scip, &(seeedPropData->newSeeeds), 1) );
      seeedPropData->newSeeeds[0] = NULL;
   }

   /* user dialog to select wanted detector, apply it and handle the returned seeeds, if any */
   finished = FALSE;
   while( !finished )
   {
      result = SCIP_DIDNOTFIND;
      /* list the detectors implementing the specified callback by name with a leading number */
      j = 1;
      SCIPinfoMessage(scip, NULL, "Available detectors:\n");
      for( i = 0; i < ndetectors; ++i )
      {
         SCIPinfoMessage(scip, NULL, "%d)", j);
         SCIPinfoMessage(scip, NULL, "%s\n", detectors[i]->name);
         ++j;
      }
      commandlen = 0;
      while( commandlen == 0 )
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
            "Type in the name or number of the detector that you want to use (or \"none\"): \nGCG/toolbox> ", &command,
            &endoffile) );
         commandlen = strlen(command);
      }

      if( !strncmp( command, "none", commandlen) == 0 && !strncmp( command, "quit", commandlen) == 0 )
      {
         for( i = 0; i < ndetectors; ++i )
         {
            sprintf(stri, "%d", i+1); //used for matching numberings in the list, off-by-one since detectors start with 0
            if( strncmp( command, detectors[i]->name, commandlen) == 0 || strncmp( command, stri, commandlen ) == 0 )
            {
               if( action == PROPAGATE )
                  SCIP_CALL( detectors[i]->propagateFromToolbox(scip, detectors[i], seeedPropData, &result, dialoghdlr, dialog) );
               else if( action == FINISH )
                  SCIP_CALL( detectors[i]->finishFromToolbox(scip, detectors[i], seeedPropData, &result, dialoghdlr, dialog) );
               else if( action == POSTPROCESS )
                  SCIP_CALL( detectors[i]->postprocessSeeed(scip, detectors[i], seeedPropData, &result) );
               break;
            }
         }
      }
      else
      {
         finished = TRUE;
         continue;
      }
      if( result == SCIP_SUCCESS )
      {
         if( action != POSTPROCESS )
         {
            SCIPinfoMessage(scip, NULL, "Considering implicits of newly found seeed(s)...\n");
            for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
            {
               assert(seeedPropData->newSeeeds[i] != NULL);
               seeedPropData->newSeeeds[i]->considerImplicits( ); //There may be open vars/cons left that were not matched
            }
            
            SCIPinfoMessage(scip, NULL, "\nSeeed was successfully %s, %d potentially new seeed(s) found.\n", actiontype,
               seeedPropData->nNewSeeeds);
            
            displayinfo = TRUE;
            if( seeedPropData->nNewSeeeds > 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, 
                     "More than one seeed found. Do you want to display information about all found seeeds anyway? (\"yes\"/\"no\")?\nGCG/toolbox> ",
                     &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "no", commandlen) == 0 )
               {
                  displayinfo = FALSE;
               }
               else if( strncmp( command, "quit", commandlen) == 0 )
               {
                  finished = TRUE;
                  continue;
               }
            }

            if( displayinfo )
            {
               for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
               {
                  seeedPropData->newSeeeds[i]->displayInfo( 0 );
               }
            }

            if( seeedPropData->nNewSeeeds == 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, 
                     "Do you want to visualize the new seeed (\"yes\"/\"no\")?\nGCG/toolbox> ", &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "yes", commandlen) == 0 )
               {
                  SCIP_CALL( SCIPconshdlrDecompSelectVisualize(scip, dialoghdlr, dialog ) );
               }
               else if( strncmp( command, "quit", commandlen) == 0 )
               {
                  finished = TRUE;
                  continue;
               }
            }

            SCIPinfoMessage(scip, NULL, "\nSaving newly found seeeds...\n\n");
            for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
            {
               conshdlrdata->curruserseeed = new gcg::Seeed( seeedPropData->newSeeeds[i] );
               SCIP_CALL( SCIPconshdlrDecompUserSeeedFlush(scip) );
               assert(conshdlrdata->curruserseeed == NULL);
            }

            if( seeedPropData->nNewSeeeds == 1 )
            {
               commandlen = 0;
               while( commandlen == 0 )
               {
                  SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, 
                     "\nDo you want to continue the decomposition with the new Seeed (\"continue\"), or continue with the previous Seeed (\"previous\")?\nGCG/toolbox> ",
                     &command, &endoffile) );
                  commandlen = strlen(command);
               }
               if( strncmp( command, "continue", commandlen) == 0 )
               {
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->newSeeeds[0]);
               }
               else
               {
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
               }
            }
            else
            {
               conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
            }
            finished = TRUE;
            continue;
         }
         else if( action == POSTPROCESS )
         {
            SCIPinfoMessage(scip, NULL, "\nSeeed successfully %s. %d seeed(s) found in the process.\n", actiontype,
               seeedPropData->nNewSeeeds);

            commandlen = 0;
            while( commandlen == 0 )
            {
               SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
                  "Do you want to save all found seeeds (\"all\") or none (\"none\")?\nGCG/toolbox> ", &command, &endoffile) );
               commandlen = strlen(command);
            }
            if( strncmp(command, "all", commandlen) == 0 )
            {
               SCIPinfoMessage(scip, NULL, "Storing seeeds...\n");
               for( i = 0; i < seeedPropData->nNewSeeeds; ++i )
               {
                  conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->newSeeeds[i]);
                  SCIP_CALL( SCIPconshdlrDecompUserSeeedFlush(scip) );
               }
               conshdlrdata->curruserseeed = new gcg::Seeed(seeedPropData->seeedToPropagate);
               SCIPinfoMessage(scip, NULL, "\nAll seeeds stored successfully!\n");
            }
            finished = TRUE;
            continue;
         }
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "Seeed could not be %s.\n", actiontype);

         commandlen = 0;
         while( commandlen == 0 )
         {
            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
               "Do you want to select another detector (\"detector\") or return to the previous menu (\"previous\")?\nGCG/toolbox> ",
               &command, &endoffile) );
            commandlen = strlen(command);
         }
         if( strncmp( command, "detector", commandlen) == 0 )
         {
            continue;
         }
         else
         {
            finished = TRUE;
            continue;
         }
      }
   }

   SCIPfreeMemoryArrayNull( scip, &(seeedPropData->newSeeeds) );
   delete seeedPropData->seeedToPropagate;
   seeedPropData->newSeeeds = NULL;
   seeedPropData->nNewSeeeds = 0;
   delete seeedPropData;

   SCIPfreeBufferArray(scip, &detectors);
   return SCIP_OKAY;
}


/** Finishes a seeed created/modified in the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPconshdlrDecompToolboxFinishSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, FINISH);
}


/** Propagates a seeed created/modified in the toolbox
 *
 * @returns SCIP status */
static
SCIP_RETCODE SCIPconshdlrDecompToolboxPropagateSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, PROPAGATE);
}


/** Postprocesses a seeed created/modified in the toolbox
 *
 * @returns SCIP status*/
static
SCIP_RETCODE SCIPconshdlrDecompToolboxPostprocessSeeed(
   SCIP*                   scip,       /**< SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /**< dialog handler for user input management */
   SCIP_DIALOG*            dialog      /**< dialog for user input management */
   )
{
   return SCIPconshdlrDecompToolboxActOnSeeed(scip, dialoghdlr, dialog, POSTPROCESS);
}


/*
 * @brief method to handle and moderate user input for modifying decompositions
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolboxModify(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;
   int commandlen;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool selectedsomeseeed;

   selectedsomeseeed = TRUE;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }
   /* 1) update list of interesting seeeds */
   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );

   /* 2) while user has not aborted: show current list extract */
   while ( !finished )
   {

      SCIP_CALL( SCIPconshdlrDecompShowListExtractHeader(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowListExtract(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please choose an existing partial decomposition for modification (type \"choose <id>\" or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen = strlen(command);

      /* case distinction: */
      if( strncmp( command, "back", commandlen) == 0 )
      {
         conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
         if(conshdlrdata->startidvisu < 0 )
            conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "next", commandlen) == 0 )
      {
         conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
         if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
            conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }
      if( strncmp( command, "top", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = 0;
         continue;
      }
      if( strncmp( command, "end", commandlen) == 0 )
      {
         conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
         continue;
      }

      if( strncmp( command, "quit", commandlen) == 0 )
      {
         finished = TRUE;
         selectedsomeseeed = FALSE;
         continue;
      }

      if( strncmp( command, "choose", commandlen) == 0 )
      {
         SCIP_RETCODE retcode = SCIPconshdlrDecompToolboxChoose(scip, dialoghdlr, dialog );
	 if (retcode != SCIP_OKAY) 
	 {
	    selectedsomeseeed = FALSE;
	    continue;
	 }
	 else
	 {
	    selectedsomeseeed = TRUE;
	    finished = TRUE;
	    break;
	 }
      }

      if( strncmp( command, "abort", commandlen) == 0 )
      {
         finished = TRUE;
         selectedsomeseeed = FALSE;
         continue;
      }

      if( strncmp( command, "change number displayed", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompModifyNVisualized(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "help", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompShowHelp(scip) );
         continue;
      }

      if( strncmp( command, "visualize", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectVisualize(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "finishseeed", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }
   finished = FALSE;
   while ( !finished && selectedsomeseeed )
   {
      int commandlen2;
      SCIP_Bool success;

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyFinish(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         gcg::Seeedpool* seeedpool;
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         if( seeedpool == NULL )

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
      if( strncmp( command, "finishseeed", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }
   return SCIP_OKAY;
}


/*
 * @brief method to handle and moderate user input for creating new decompositions by the user
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return data structure
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolboxCreate(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   char* command;
   SCIP_Bool endoffile;
   SCIP_Bool         finished;
   int commandlen;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }

   /* create new decomposition */
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ", &command, &endoffile) );
   commandlen = strlen(command);

   if( conshdlrdata->curruserseeed != NULL )
      delete conshdlrdata->curruserseeed;

   gcg::Seeedpool* seeedpool;
   SCIP_Bool isfromunpresolved;

   while( (strncmp( command, "presolved", commandlen) != 0 && strncmp( command, "unpresolved", commandlen) != 0) || commandlen == 0)
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Invalid input. Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ", &command, &endoffile) );
      commandlen = strlen(command);
   }

   /* case distinction: */
   if( strncmp( command, "presolved", commandlen) == 0 )
   {
      isfromunpresolved = FALSE;
      if (conshdlrdata->seeedpool != NULL )
         seeedpool = conshdlrdata->seeedpool;
      else
      {
         if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
         {
            SCIPinfoMessage(scip, NULL, "Problem is not presolved yet. Please presolve it first!\n");
            return SCIP_OKAY;
         }

         conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpool;
      }
   }
   else
   {
      isfromunpresolved = TRUE;
      if ( conshdlrdata->seeedpoolunpresolved == NULL )
         conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
       seeedpool = conshdlrdata->seeedpoolunpresolved;

   }
   if( seeedpool == NULL )
   {
      if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED )

      {
         if (conshdlrdata->seeedpool == NULL )
            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpool;
      }
      else
      {
         if ( conshdlrdata->seeedpoolunpresolved == NULL)
            conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE, SCIPconshdlrDecompDetectBenders(scip));
         seeedpool = conshdlrdata->seeedpoolunpresolved;
      }

   }
   conshdlrdata->curruserseeed = new gcg::Seeed( scip, SCIPconshdlrDecompGetNextSeeedID(scip), seeedpool );
   conshdlrdata->curruserseeed->setIsFromUnpresolved(isfromunpresolved);
   finished = FALSE;
   while ( !finished )
   {
      int commandlen2;
      SCIP_Bool success;

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ", &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->curruserseeed->isFromUnpresolved() )
            seeedpool = conshdlrdata->seeedpoolunpresolved;
         else
            seeedpool = conshdlrdata->seeedpool;
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         if( seeedpool == NULL )

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }
   return SCIP_OKAY;
}

/*
 * method to handle and moderate user input for creating new decompositions
 * and modifying existing decompositions by the user
 *
 * @param scip SCIP data structure
 * @param dialoghdlr dialog handler to handle user input
 * @param dialog dialog to handle user input
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompExecToolbox(
   SCIP*                   scip,       /* SCIP data structure */
   SCIP_DIALOGHDLR*        dialoghdlr, /* dialog handler for user input management */
   SCIP_DIALOG*            dialog      /* dialog for user input management */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_Bool         finished;
   char* command;
   SCIP_Bool endoffile;
   int commandlen;
   SCIP_Bool selectedsomeseeed;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );
   finished = FALSE;

   selectedsomeseeed = TRUE;
   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT )
   {
      SCIPinfoMessage(scip, NULL, "No problem is loaded. Please read in a model first.\n");
      return SCIP_OKAY;
   }
   if( (int)conshdlrdata->listall->size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "No decompositions available. Please detect first.\n");
      return SCIP_OKAY;
   }
   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
   {
      SCIP_CALL( SCIPtransformProb(scip) );
      SCIPinfoMessage(scip, NULL, "Applied tranformation to problem.\n");
   }

   commandlen = 0;

   /* Does user want to modify existing or create a new partial decomposition ?*/
   while( (strncmp( command, "modify", commandlen) != 0 && strncmp( command, "create", commandlen) != 0) || commandlen == 0)
   {
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Do you want to modify an existing (\"modify\") or create a new partial decomposition (\"create\")? : \nGCG/toolbox> ",
         &command, &endoffile) );
      commandlen = strlen(command);
   }


   /* case distinction: */
   if( strncmp( command, "modify", commandlen) == 0 )
   {
      /* 1) update list of interesting seeeds */

         SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );


         /* 2) while user has not aborted: show current list extract */

         while ( !finished )
         {
            int commandlen2;

            SCIP_CALL( SCIPconshdlrDecompShowListExtractHeader(scip) );

            SCIP_CALL( SCIPconshdlrDecompShowListExtract(scip) );

            SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
               "Please choose an existing partial decomposition for modification (type \"choose <id>\" or \"h\" for help) : \nGCG/toolbox> ",
               &command, &endoffile) );

            commandlen2 = strlen(command);

            /* case distinction: */
            if( strncmp( command, "back", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu -= conshdlrdata->selectvisulength;
               if(conshdlrdata->startidvisu < 0 )
                  conshdlrdata->startidvisu = 0;
               continue;
            }
            if( strncmp( command, "next", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu += conshdlrdata->selectvisulength;
               if( conshdlrdata->startidvisu > (int) conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
                  conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
               continue;
            }
            if( strncmp( command, "top", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu = 0;
               continue;
            }
            if( strncmp( command, "end", commandlen2) == 0 )
            {
               conshdlrdata->startidvisu = conshdlrdata->listall->size() - conshdlrdata->selectvisulength ;
               continue;
            }

            if( strncmp( command, "quit", commandlen2) == 0 )
            {
               finished = TRUE;
               selectedsomeseeed = FALSE;
               continue;
            }


            if( strncmp( command, "choose", commandlen2) == 0 )
            {
               SCIP_RETCODE retcode = SCIPconshdlrDecompToolboxChoose(scip, dialoghdlr, dialog );
	       if (retcode != SCIP_OKAY) 
	       {
		  selectedsomeseeed = FALSE;
		  continue;
	       }
	       else
	       {
		  finished = TRUE;
		  break;
	       }
            }


            if( strncmp( command, "abort", commandlen2) == 0 )
            {
               finished = TRUE;
               selectedsomeseeed = FALSE;
               continue;
            }

            if( strncmp( command, "change number displayed", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPconshdlrDecompModifyNVisualized(scip, dialoghdlr, dialog) );
               continue;
            }

            if( strncmp( command, "help", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPconshdlrDecompShowHelp(scip) );
               continue;
            }

            if( strncmp( command, "visualize", commandlen2) == 0 )
            {
               SCIP_CALL(SCIPconshdlrDecompSelectVisualize(scip, dialoghdlr, dialog ) );
               continue;
            }

            if( strncmp( command, "propagate", commandlen2) == 0 )
            {
               SCIP_CALL( SCIPconshdlrDecompToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
               continue;
            }

            if( strncmp( command, "finishseeed", commandlen2) == 0 )
            {
               SCIP_CALL( SCIPconshdlrDecompToolboxFinishSeeed(scip, dialoghdlr, dialog) );
               continue;
            }

            if( strncmp( command, "postprocess", commandlen2) == 0 )
            {
               SCIP_CALL( SCIPconshdlrDecompToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
               continue;
            }
         }
   } /* finished yes == modify */
   else
   {
      /* create new decomposition */
      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ",
         &command, &endoffile) );
      commandlen = strlen(command);

      if( conshdlrdata->curruserseeed != NULL )
         delete conshdlrdata->curruserseeed;

      gcg::Seeedpool* seeedpool;
      SCIP_Bool isfromunpresolved;

      while( (strncmp( command, "presolved", commandlen) != 0 && strncmp( command, "unpresolved", commandlen) != 0) || commandlen == 0)
      {
         SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
            "Invalid input. Should the new partial decomposition be for the presolved or the unpresolved problem? (type \"presolved\" or \"unpresolved\") : \nGCG/toolbox> ",
            &command, &endoffile) );
         commandlen = strlen(command);
      }

      /* case distinction: */
      if( strncmp( command, "presolved", commandlen) == 0 )
      {
         isfromunpresolved = FALSE;
         if (conshdlrdata->seeedpool != NULL )
            seeedpool = conshdlrdata->seeedpool;
         else
         {
            if( SCIPgetStage(scip) < SCIP_STAGE_PRESOLVED )
            {
               SCIPinfoMessage(scip, NULL, "Problem is not presolved yet. Please presolve it first!\n");
               return SCIP_OKAY;
            }

            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpool;
         }
      }
      else
      {
         isfromunpresolved = TRUE;
         if ( conshdlrdata->seeedpoolunpresolved == NULL )
            conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE,
               SCIPconshdlrDecompDetectBenders(scip));
          seeedpool = conshdlrdata->seeedpoolunpresolved;

      }
      if( seeedpool == NULL )
      {
         if( SCIPgetStage(scip) >= SCIP_STAGE_PRESOLVED )

         {
            if (conshdlrdata->seeedpool == NULL )
               conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpool;
         }
         else
         {
            if ( conshdlrdata->seeedpoolunpresolved == NULL)
               conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE,
                  SCIPconshdlrDecompDetectBenders(scip));
            seeedpool = conshdlrdata->seeedpoolunpresolved;
         }

      }
      conshdlrdata->curruserseeed = new gcg::Seeed( scip, SCIPconshdlrDecompGetNextSeeedID(scip), seeedpool );
      conshdlrdata->curruserseeed->setIsFromUnpresolved(isfromunpresolved);
   }

   /* curruserseeed is ready to modify */

   finished = FALSE;
   while ( !finished && selectedsomeseeed )
   {
      int commandlen2;
      SCIP_Bool success;

      SCIP_CALL( SCIPconshdlrDecompShowCurrUserSeeedInfo(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowToolboxInfo(scip) );

      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog,
         "How do you want to proceed the with the current decomposition? (or \"h\" for help) : \nGCG/toolbox> ",
         &command, &endoffile) );

      commandlen2 = strlen(command);

      /* case distinction: */
      if( strncmp( command, "conss", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyConss(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "vars", commandlen2) == 0 )
      {
         SCIPconshdlrDecompToolboxModifyVars(scip, dialoghdlr, dialog);
         continue;
      }
      if( strncmp( command, "refine", commandlen2) == 0 )
      {
         if( conshdlrdata->lastuserseeed != NULL)
            delete conshdlrdata->lastuserseeed;
         conshdlrdata->lastuserseeed = new gcg::Seeed( conshdlrdata->curruserseeed) ;
         conshdlrdata->curruserseeed->considerImplicits();
         continue;
      }

      if( strncmp( command, "quit", commandlen2) == 0 )
      {
         gcg::Seeedpool* seeedpool;
         if( !conshdlrdata->curruserseeed->isFromUnpresolved() && conshdlrdata->seeedpool == NULL )
            SCIPconshdlrDecompCreateSeeedpool(scip);

         seeedpool = ( conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool);
         assert( seeedpool != NULL );

         conshdlrdata->curruserseeed->sort();
         conshdlrdata->curruserseeed->considerImplicits();
         conshdlrdata->curruserseeed->calcHashvalue();
         assert( conshdlrdata->curruserseeed->checkConsistency() );



         if( conshdlrdata->curruserseeed->isComplete() )
         {
            seeedpool->addSeeedToFinished(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         } else
         {
            seeedpool->addSeeedToIncomplete(conshdlrdata->curruserseeed, &success);
            if( !success )
            {
               delete conshdlrdata->curruserseeed;
            }
         }
         conshdlrdata->curruserseeed = NULL;
         finished = TRUE;


         continue;
      }

      if( strncmp( command, "undo", commandlen2) == 0 )
      {
         if ( conshdlrdata->lastuserseeed == NULL )
            SCIPdialogMessage(scip, NULL, " nothing to be undone \n");
         else
         {
            delete conshdlrdata->curruserseeed;
            conshdlrdata->curruserseeed = conshdlrdata->lastuserseeed;
            conshdlrdata->lastuserseeed = NULL;
         }
         continue;
      }


      if( strncmp( command, "visualize", commandlen2) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectVisualizeCurrentUserSeeed(scip, dialoghdlr, dialog ) );
         continue;
      }

      if( strncmp( command, "propagate", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPropagateSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
      if( strncmp( command, "finish", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxFinishSeeed(scip, dialoghdlr, dialog) );
         continue;
      }

      if( strncmp( command, "postprocess", commandlen2) == 0 )
      {
         SCIP_CALL( SCIPconshdlrDecompToolboxPostprocessSeeed(scip, dialoghdlr, dialog) );
         continue;
      }
   }

   return SCIP_OKAY;
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

/*
 * @brief help method to access seeedpool for unpresolved problem
 * @TODO: consider deleting this method will be deleted if the corresponding wrapper classes are introduced
 * @returns pointer to seeedpool wrapper data structure
 */
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolUnpresolvedExtern(
   SCIP*                 scip                /* SCIP data structure */
   ){
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SEEEDPOOL_WRAPPER* help;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   help = (SEEEDPOOL_WRAPPER*) conshdlrdata->seeedpoolunpresolved;

   return (SEEEDPOOL_WRAPPER*) help;
}

/*
 * @brief help method to access seeedpool for transformed problem
 * @TODO: consider deleting this method will be deleted if the corresponidng wrapper classes are introduced
 * @returns pointer to seeedpool wrapper data structure
 */
SEEEDPOOL_WRAPPER* SCIPconshdlrDecompGetSeeedpoolExtern(
   SCIP*                 scip                /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SEEEDPOOL_WRAPPER* help;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   help = (SEEEDPOOL_WRAPPER*) conshdlrdata->seeedpool;

   return (SEEEDPOOL_WRAPPER*) help;

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
 * @param DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)) method to refine a partial decomposition when called by user from console (or NULL)
 * @param DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)) method to complete a partial decomposition when called by user from console (or NULL)
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
   DEC_DECL_PROPAGATEFROMTOOLBOX((*propagateFromToolboxDetector)),   /* propagation from toolbox method of detector
                                                                          (or NULL) */
   DEC_DECL_FINISHFROMTOOLBOX((*finishFromToolboxDetector)),   /* finish from toolbox method of detector (or NULL) */
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

   detector->decdata = detectordata;
   detector->name = name;
   detector->description = description;
   detector->decchar = decchar;

   detector->freeDetector = freeDetector;
   detector->initDetector = initDetector;
   detector->exitDetector = exitDetector;
   detector->detectStructure = detectStructure;

   detector->propagateSeeed = propagateSeeedDetector;
   detector->propagateFromToolbox = propagateFromToolboxDetector;
   detector->finishFromToolbox = finishFromToolboxDetector;
   detector->finishSeeed = finishSeeedDetector;
   detector->postprocessSeeed = postprocessSeeedDetector;
   detector->setParamAggressive =  setParamAggressiveDetector;
   detector->setParamDefault =  setParamDefaultDetector;
   detector->setParamFast =  setParamFastDetector;

   detector->decchar = decchar;

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
 * @brief creates a user seeed for the problem
 *  @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompCreateUserSeeed(
   SCIP*                 scip,               /* SCIP data structure */
   SCIP_Bool             presolved,          /* should the user seeed be created for the presolved problem */
   SCIP_Bool             markedincomplete    /* should the user seeed be a partial one */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed != NULL )
   {
      SCIPwarningMessage(scip, "there is a current user seeed, it is going to be flushed..!\n");
      SCIP_CALL( SCIPconshdlrDecompUserSeeedFlush(scip) );
   }

   currseeedpool = presolved ? conshdlrdata->seeedpool : conshdlrdata->seeedpoolunpresolved;

   assert( currseeedpool != NULL );
   assert( conshdlrdata->curruserseeed == NULL );

   conshdlrdata->curruserseeed = new gcg::Seeed(scip, currseeedpool->getNewIdForSeeed(), currseeedpool );

      conshdlrdata->curruserseeed->setIsFromUnpresolved( !presolved );

   if( markedincomplete )
      conshdlrdata->curruserseeed->setUsergiven(USERGIVEN::PARTIAL);
   else
      conshdlrdata->curruserseeed->setUsergiven(USERGIVEN::COMPLETED_CONSTOMASTER);

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
 * @brief method to update the list of incomplete decompositions in "explore" submenu ( this list changes due to new decompositions,  modified, decompositions or changes of the score
 *
 *  unselects all seeeds and updates list:
 * 1) add presolved finished
 * 2) add presolved unfinished
 * 3) add unpresolved finished
 * 4) add unpresolved partial
 *
 * @param scip SCIP data structure
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUpdateSeeedlist(
   SCIP*          scip     /* SCIP data structure */
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

   assert( SCIPconshdlrDecompCheckConsistency(scip) );

   conshdlrdata->startidvisu = 0;
   SCIPconshdlrdataDecompUnselectAll(scip);
   conshdlrdata->listall->clear();


   if( conshdlrdata->hasrun && conshdlrdata->seeedpool == NULL && conshdlrdata->seeedpoolunpresolved == NULL)
      return SCIP_OKAY;

   /* sort decomposition and finished seeeds according to max white score */
   SCIP_CALL( DECconshdlrDecompSortDecompositionsByScore(scip) );

   /* add seeeds to list */
   /* 1) add presolved finished */
   for( i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNFinishedSeeeds(); ++i )
    {
       SeeedPtr seeed;
       seeed = conshdlrdata->seeedpool->getFinishedSeeed( i );

       conshdlrdata->listall->push_back(seeed);
    }

   /* 2) add presolved unfinished */
   for( i = 0; conshdlrdata->seeedpool != NULL && i < conshdlrdata->seeedpool->getNIncompleteSeeeds(); ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpool->getIncompleteSeeed( i );

      conshdlrdata->listall->push_back(seeed);
   }


   /* 3) add unpresolved finished */
   for( i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() ; ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed(i);
      seeed->setIsFromUnpresolved( TRUE );
      conshdlrdata->listall->push_back(seeed);
   }

   /* 4) add unpresolved partial */
   for( i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNIncompleteSeeeds(); ++i )
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpoolunpresolved->getIncompleteSeeed( i );
      seeed->setIsFromUnpresolved( TRUE );

      conshdlrdata->listall->push_back(seeed);
   }


   return SCIP_OKAY;
}


/*
 * @brief sets the number of blocks
 *
 * set the number of blocks in the current user seeed (which is used for user input (read or modify) )
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /* SCIP data structure */
   int                   nblocks              /* number of blocks */
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

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   conshdlrdata->curruserseeed->setNBlocks(nblocks);

   return SCIP_OKAY;
}


/*
 * @brief sets a constraint by name to a block in the current user seeed
 * @param scip SCIP data structure
 * @param consname name of the constraint that should be set to a block
 * @param blockid index of the block the constraint should be assigned to
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToBlock(
   SCIP*                 scip,                /* SCIP data structure */
   const char*           consname,            /* name of the constraint */
   int                   blockid              /* block index ( counting from 0) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONS* cons;
   gcg::Seeedpool* currseeedpool;
   int consindex;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);

   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   cons = conshdlrdata->curruserseeed->isFromUnpresolved() ? (SCIPfindOrigCons(scip, consname ) == NULL ? SCIPfindCons(scip, consname ): SCIPfindOrigCons(scip, consname )) : SCIPfindCons(scip, consname );
   consindex = currseeedpool->getIndexForCons( cons ) ;

   if( blockid >= conshdlrdata->curruserseeed->getNBlocks() )
         conshdlrdata->curruserseeed->setNBlocks(blockid+1);
   conshdlrdata->curruserseeed->bookAsBlockCons(consindex, blockid);

   return SCIP_OKAY;
}


/*
 * @brief sets a constraint by name to master in the current user seeed
 * @param consname of the constraint that should be set to master
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToMaster(
   SCIP*                 scip,                /* SCIP data structure */
   const char*           consname             /* name of constraint */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   int consindex;
   SCIP_CONS* cons;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;

   cons = conshdlrdata->curruserseeed->isFromUnpresolved() ? SCIPfindOrigCons(scip, consname ) : SCIPfindCons(scip, consname );
   consindex = currseeedpool->getIndexForCons( cons );

   conshdlrdata->curruserseeed->bookAsMasterCons(consindex);

   return SCIP_OKAY;

}


/*
 * @brief sets a variable by name to a block in the current user seeed
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToBlock(
   SCIP*                 scip,                /* SCIP data structure */
   const char*           varname,             /* name of the variable */
   int                   blockid              /* block index ( counting from 0) */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   int varindex;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

   if( blockid >= conshdlrdata->curruserseeed->getNBlocks() )
      conshdlrdata->curruserseeed->setNBlocks(blockid+1);
   conshdlrdata->curruserseeed->bookAsBlockVar(varindex, blockid);

   return SCIP_OKAY;
}


/*
 * @brief sets a variable by name to the master in the current user seeed
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToMaster(
   SCIP*                 scip,                /* SCIP data structure */
   const char*           varname              /* name of the variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   int varindex;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   currseeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

   conshdlrdata->curruserseeed->bookAsMasterVar(varindex);

   return SCIP_OKAY;

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
 * @brief returns the total detection time
 * @param scip SCIP data structure
 * @returns total detection time
 */
SCIP_Real SCIPconshdlrDecompGetCompleteDetectionTime(
    SCIP*                 scip   /* SCIP data structure */
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

   return SCIPclockGetTime( conshdlrdata->completedetectionclock );
}


/*
 * @brief sets a variable by name to the linking variables in the current user seeed
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToLinking(
   SCIP*                 scip,                /* SCIP data structure */
   const char*           varname              /* name of the variable */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   int varindex;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

      currseeedpool = conshdlrdata->curruserseeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
      varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

      conshdlrdata->curruserseeed->bookAsLinkingVar(varindex);

      return SCIP_OKAY;

}


/*
 * finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if
 * complete etc
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedFlush(
   SCIP*                 scip                 /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* currseeedpool;
   SeeedPtr        seeed;

   char const *            usergiveninfo;
   char const *            presolvedinfo;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one..!\n");
      return SCIP_OKAY;
   }

   seeed = conshdlrdata->curruserseeed;
   currseeedpool = seeed->isFromUnpresolved() ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   seeed->setSeeedpool(currseeedpool);
   seeed->flushBooked();

   if( seeed->shouldCompletedByConsToMaster() )
   {
      for( int opencons = 0; opencons < seeed->getNOpenconss(); ++opencons)
         seeed->bookAsMasterCons( seeed->getOpenconss()[opencons] );
      seeed->flushBooked();
   }

   seeed->considerImplicits();
   currseeedpool->prepareSeeed(conshdlrdata->curruserseeed);

   if( !seeed->checkConsistency() )
   {
      SCIPconshdlrDecompUserSeeedReject(scip);
      SCIPwarningMessage(scip, "seeed that was given by the user was rejected because of inconsistencies! \n");
      return SCIP_OKAY;
   }
   conshdlrdata->curruserseeed->buildDecChainString();
   if( conshdlrdata->curruserseeed->isComplete() )
   {
      if( !seeed->shouldCompletedByConsToMaster() )
         conshdlrdata->curruserseeed->setUsergiven( gcg::USERGIVEN::COMPLETE );
      /* stems from presolved problem? */
      if( !conshdlrdata->curruserseeed->isFromUnpresolved() )
      {
         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForPresolved(scip, conshdlrdata->curruserseeed));

      }
      /* stems from unpresolved problem */
      else
      {

         SCIP_CALL( SCIPconshdlrDecompAddCompleteSeeedForUnpresolved(scip, seeed) );

         if ( conshdlrdata->seeedpool != NULL ) /* seeedpool for presolved problem already exist try to translate seeed */
         {
            std::vector<SeeedPtr> seeedtotranslate(0);
            std::vector<SeeedPtr> newseeeds(0);
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
      conshdlrdata->curruserseeed->setUsergiven( gcg::USERGIVEN::PARTIAL );

      if ( !conshdlrdata->curruserseeed->isFromUnpresolved() )
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForPresolved(scip, conshdlrdata->curruserseeed) );
      else
         SCIP_CALL(SCIPconshdlrDecompAddPartialSeeedForUnpresolved(scip, conshdlrdata->curruserseeed) );
   }

   /* set statistics */
   {
      int nvarstoblock = 0;
      int nconsstoblock = 0;

      for ( int b = 0; b < conshdlrdata->curruserseeed->getNBlocks(); ++b )
      {
         nvarstoblock += conshdlrdata->curruserseeed->getNVarsForBlock(b);
         nconsstoblock += conshdlrdata->curruserseeed->getNConssForBlock(b);
      }
      conshdlrdata->curruserseeed->setDetectorPropagated(NULL);

      conshdlrdata->curruserseeed->addClockTime(0.);
      conshdlrdata->curruserseeed->addPctVarsFromFree( (nvarstoblock + conshdlrdata->curruserseeed->getNMastervars() +conshdlrdata->curruserseeed->getNLinkingvars())/(SCIP_Real) conshdlrdata->curruserseeed->getNVars()  );
      conshdlrdata->curruserseeed->addPctVarsToBlock((nvarstoblock )/(SCIP_Real) conshdlrdata->curruserseeed->getNVars() );
      conshdlrdata->curruserseeed->addPctVarsToBorder( (conshdlrdata->curruserseeed->getNMastervars() +conshdlrdata->curruserseeed->getNLinkingvars())/(SCIP_Real) conshdlrdata->curruserseeed->getNVars() ) ;
      conshdlrdata->curruserseeed->addPctConssToBorder( (conshdlrdata->curruserseeed->getNMasterconss() ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() ) ;
      conshdlrdata->curruserseeed->addPctConssFromFree( (conshdlrdata->curruserseeed->getNMasterconss() + nconsstoblock ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() ) ;
      conshdlrdata->curruserseeed->addPctConssToBlock( (nconsstoblock ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() );
      conshdlrdata->curruserseeed->addNNewBlocks(conshdlrdata->curruserseeed->getNBlocks());
   }

   conshdlrdata->curruserseeed->findVarsLinkingToMaster();
   conshdlrdata->curruserseeed->findVarsLinkingToStairlinking();


   if( conshdlrdata->curruserseeed->getUsergiven() == gcg::USERGIVEN::PARTIAL )
      usergiveninfo = "partial";
   if( conshdlrdata->curruserseeed->getUsergiven() == gcg::USERGIVEN::COMPLETE )
      usergiveninfo = "complete";
   if( conshdlrdata->curruserseeed->getUsergiven() == gcg::USERGIVEN::COMPLETED_CONSTOMASTER )
         usergiveninfo = "complete";
   if( conshdlrdata->curruserseeed->isFromUnpresolved() )
         presolvedinfo = "unpresolved";
   else presolvedinfo = "presolved";


   SCIPverbMessage(scip, SCIP_VERBLEVEL_HIGH, NULL, " added %s decomp for %s problem with %d blocks and %d masterconss, %d linkingvars, "
      "%d mastervars, and max white score of %s %f \n", usergiveninfo, presolvedinfo,
      conshdlrdata->curruserseeed->getNBlocks(), conshdlrdata->curruserseeed->getNMasterconss(),
      conshdlrdata->curruserseeed->getNLinkingvars(), conshdlrdata->curruserseeed->getNMastervars(), (conshdlrdata->curruserseeed->isComplete() ? " " : " at best "),
      conshdlrdata->curruserseeed->getScore(SCORETYPE::MAX_WHITE) );

   conshdlrdata->curruserseeed = NULL;

   return SCIP_OKAY;
}

/*
 * @brief rejects and deletes the current user seeed
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedReject(
   SCIP*                 scip                 /* SCIP data structure */
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

   if( conshdlrdata->curruserseeed == NULL )
   {
      SCIPwarningMessage(scip, "there is no current user seeed, you have to create one  before you can reject it\n");
      return SCIP_OKAY;
   }

   delete conshdlrdata->curruserseeed;

   conshdlrdata->curruserseeed = NULL;

   return SCIP_OKAY;
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
 * @brief returns if there is a decomposition that is currently selected by the user (done in explore menu)
 * @param scip SCIP data structure
 * @returns TRUE if there is a decomposition that is currently selected by the user (done in explore menu)
 */
SCIP_Bool SCIPconshdlrDecompExistsSelected(
   SCIP* scip     /* SCIP data structure */
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

    return conshdlrdata->selectedexists;
}


/*
 * @brief initilizes the candidates data structures with selected seeeds (or all if there are no selected seeeds) and sort them according to the current scoretype
 * @param scip SCIP data structure
 * @param updatelist whether or not the seeed list should be updated
 * @returns SCIP return code
 */
SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,             /* SCIP data structure */
   SCIP_Bool updatelist    /* whether to update seeed list beforehand */
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
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdebugMessage("Starting decomposition candidate choosing \n");

   assert(conshdlrdata->candidates != NULL);

   conshdlrdata->candidates->clear();

   if( updatelist )
      SCIP_CALL(SCIPconshdlrDecompUpdateSeeedlist(scip) );

   for( size_t selid = 0; selid < conshdlrdata->selected->size(); ++selid )
   {
      selectedseeeds.push_back(conshdlrdata->listall->at(conshdlrdata->selected->at(selid) ) );
   }

   if ( selectedseeeds.size() == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL, "currently no decomposition is selected, hence every known decomposition is considered: \n");
      selectedseeeds = *conshdlrdata->listall;
      SCIPverbMessage(scip, SCIP_VERBLEVEL_FULL,  NULL,  "number that is examined: %d \n", selectedseeeds.size() );
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
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) ) );
      }
      if( seeed->isComplete() && seeed->isFromUnpresolved() )
      {
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, SCIPconshdlrDecompAdaptScore(scip, seeed->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) ) ) );
      }
   }

   seeediter = finished.begin();
   seeediterend = finished.end();

   for( ; seeediter != seeediterend; ++seeediter )
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, (*seeediter)->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) )  );
   }

   seeediter = finishedunpresolved.begin();
   seeediterend = finishedunpresolved.end();

   for( ; seeediter != seeediterend; ++seeediter )
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, SCIPconshdlrDecompAdaptScore(scip, (*seeediter)->getScore(SCIPconshdlrdataGetScoretype(conshdlrdata)) ) ) );
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
 *  2) selected list is syncron with selected information in seeeds
 *  3) selected exists is synchronized with seleced list
 *
 *  @returns true if seeed information is consistent */
SCIP_Bool SCIPconshdlrDecompCheckConsistency(
   SCIP* scip  /* SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   int i;
   int selectedcounter;

   std::vector<int> livingnoncompleteseeedids(0); /* this is a vector of seeed ids that should be living (living: there is no complete seeed having ) */
   std::vector<int>::const_iterator selectediter;
   std::vector<int>::const_iterator selectediterend;

   std::vector<SeeedPtr>::const_iterator seeediter;
   std::vector<SeeedPtr>::const_iterator seeediterend;


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

   /* 2) selected list is syncronized with selected information in seeeds */
   selectediter = conshdlrdata->selected->begin();
   selectediterend = conshdlrdata->selected->end();

   selectedcounter = 0;

   for( ; selectediter != selectediterend; ++selectediter )
   {
      SeeedPtr seeed = conshdlrdata->listall->at(*selectediter);

      if( !seeed->isSelected() )
      {
         SCIPwarningMessage(scip, "Warning: seeed %d is not selected but in slected list  \n", seeed->getID() );
         return FALSE;
      }
   }

   seeediter = conshdlrdata->listall->begin();
   seeediterend = conshdlrdata->listall->end();

   for( ; seeediter != seeediterend; ++seeediter )
   {
      if( (*seeediter)->isSelected() )
         ++selectedcounter;
   }

   if( (size_t) selectedcounter != conshdlrdata->selected->size() )
   {
      SCIPwarningMessage(scip, "Warning: there are selected seeeds not part of the list (selectedcounter: %d, nselected list> %d) \n", selectedcounter, (int) conshdlrdata->selected->size() );
      return FALSE;
   }

   /* 3) selected exists is syncronized with seleced list */

   if( conshdlrdata->selectedexists != (conshdlrdata->selected->size() > 0) )
   {
      SCIPwarningMessage(scip, "Warning: selectedexists is %d but number of selected is %d   \n", conshdlrdata->selectedexists, conshdlrdata->selected->size() );
      return FALSE;
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

   return ++conshdlrdata->seeedcounter;
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
      SCIPdebugMessage("is seeedpoolunpresolved not initilized yet but needed ? %s -> %s create it \n", (conshdlrdata->seeedpoolunpresolved == NULL ? "yes" : "no"), (conshdlrdata->seeedpoolunpresolved == NULL ? "" : "Do not")  );

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
         if( conshdlrdata->seeedpoolunpresolved != NULL && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
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
            SCIPdebugMessagePrint(scip, "start creating seeedpool for current problem \n");
            conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE, SCIPconshdlrDecompDetectBenders(scip));
            SCIPdebugMessagePrint(scip, "created seeedpool for current problem, n detectors: %d \n", conshdlrdata->ndetectors);
         }
         else
            SCIPdebugMessagePrint(scip, "seeedpool is not NULL \n");

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
        SCIP_CALL(SCIPstopClock(scip, conshdlrdata->detectorclock));
     }

      if( conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() > 0 )
         *result = SCIP_SUCCESS;

      if( conshdlrdata->seeedpoolunpresolved != NULL &&  conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() > 0 )
         *result = SCIP_SUCCESS;

      SCIPdebugMessage("Detection took %fs\n", SCIPclockGetTime(conshdlrdata->detectorclock));

   } /* end of if( !onlylegacy ) */


   if( conshdlrdata->seeedpool != NULL && SCIPgetVerbLevel(scip) >= SCIP_VERBLEVEL_FULL )
      conshdlrdata->seeedpool->printBlockcandidateInformation(scip, NULL);

   SCIP_CALL(SCIPstartClock(scip, conshdlrdata->completedetectionclock) );
   SCIPconshdlrDecompAddLegacymodeDecompositions( scip, result );
   SCIP_CALL(SCIPstopClock(scip, conshdlrdata->completedetectionclock) );

   if( *result == SCIP_DIDNOTRUN )
   {
      return SCIP_OKAY;
   }

   /* show that we done our duty */
   conshdlrdata->hasrun = TRUE;
   *result = SCIP_SUCCESS;
   SCIPconshdlrDecompChooseCandidatesFromSelected(scip, TRUE);

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


   if(  conshdlrdata->selected->size() == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition selected.\n");
      return SCIP_OKAY;
   }


   for( size_t selid = 0; selid < conshdlrdata->selected->size(); ++selid )
   {
      SeeedPtr seeed;

      seeed = conshdlrdata->listall->at(conshdlrdata->selected->at(selid) );

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
 * used before calling SCIPfreeTransform(),, if called to revoke presolving (e.g. if unpresolved decomposition is used, and
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
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip, TRUE);
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

   SCIPinfoMessage(scip, NULL, "start creating seeedpool \n");
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

   SCIPinfoMessage(scip, NULL, "finished creating seeedpool \n");
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
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip, TRUE);
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
      SCIPconshdlrDecompChooseCandidatesFromSelected(scip, TRUE);
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


/*
 * Gets the current scoretype
 * @returns the current scoretype */
scoretype SCIPconshdlrDecompGetCurrScoretype(
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

   return  static_cast<scoretype>(conshdlrdata->currscoretype);

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
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         if( written )
            SCIPinfoMessage(scip, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

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
               SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
               "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
               "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            if( written )
               SCIPinfoMessage(scip, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");

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
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
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
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      char paramname[SCIP_MAXSTRLEN];
      (void)SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detection/detectors/%s/legacymode", conshdlrdata->detectors[i]->name);

      SCIP_CALL(SCIPsetBoolParam(scip, paramname, FALSE));
      if( !quiet )
      {
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
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
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN,
            "detection/detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         if( written )
            SCIPinfoMessage(scip, NULL, "\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
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
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   std::vector<gcg::Seeed*>::const_iterator seeediter;
   std::vector<gcg::Seeed*>::const_iterator seeediterend;

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );

   seeediter = conshdlrdata->listall->begin();
   seeediterend = conshdlrdata->listall->end();

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "DECOMPINFO  \n" );
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "%d\n", (int) conshdlrdata->listall->size() );

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
