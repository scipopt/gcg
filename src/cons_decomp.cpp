/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2017 Operations Research, RWTH Aachen University       */
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
 * increasing priority until the first detector finds a suitable structure.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

//#define SCIP_DEBUG

#include <assert.h>
#include <iostream>
#include <stdio.h>
#include <sstream>
#include "cons_decomp.h"
#include "dec_connected.h"
#include "gcg.h"
#include "struct_detector.h"
#include "string.h"
#include "scip_misc.h"
#include "scip/clock.h"
#include "class_seeed.h"
#include "class_seeedpool.h"


#include <vector>
#include <iomanip>
#include <queue>
#include <fstream>
#include <algorithm>

#include <iterator>

typedef gcg::Seeed* SeeedPtr;



/* constraint handler properties */
#define CONSHDLR_NAME          "decomp"
#define CONSHDLR_DESC          "constraint handler for structure detection"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY        0 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ           -1 /**< frequency for using all instead of only the useful constraints in separation,
                                          *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS        FALSE /**< should the constraint handler be skipped, if no constraints are available? */

#define MAXNDECOMPS 5000                /**< indicates whether to create a decomposition with all constraints in the master if no other specified */

#define DEFAULT_CREATEBASICDECOMP FALSE /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
#define DEFAULT_MAXDETECTIONROUNDS 2    /**< maximal number of detection rounds */
#define DEFAULT_ENABLEORIGDETECTION TRUE /**< indicates whether to start detection for the original problem */

#define DEFAULT_CONSSCLASSNNONZENABLED                TRUE    /**<  indicates whether constraint classifier for nonzero entries is enabled */
#define DEFAULT_CONSSCLASSNNONZENABLEDORIG            TRUE     /**<  indicates whether constraint classifier for nonzero entries is enabled for the original problem */

#define DEFAULT_CONSSCLASSSCIPCONSTYPEENABLED         TRUE     /**< indicates whether constraint classifier for scipconstype is enabled */
#define DEFAULT_CONSSCLASSSCIPCONSTYPEENABLEDORIG     TRUE    /**< indicates whether constraint classifier for scipconsstype is enabled for the original problem */

#define DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLED     FALSE    /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled */
#define DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLEDORIG TRUE     /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem */

#define DEFAULT_CONSSCLASSLEVENSHTEINENABLED          FALSE    /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled */
#define DEFAULT_CONSSCLASSLEVENSHTEINENABLEDORIG      TRUE     /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem */

#define DEFAULT_VARCLASSSCIPVARTYPESENABLED           TRUE     /**< indicates whether variable classifier for scipvartypes is enabled */
#define DEFAULT_VARCLASSSCIPVARTYPESENABLEDORIG       TRUE     /**< indicates whether variable classifier for scipvartypes is enabled for the original problem */

#define DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER    10000    /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis default */
#define AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER  80000    /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis aggressive */
#define FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER       2000     /**< deactivate levenshtein constraint classifier if nrows + ncols exceeds this value for emphasis fast */

/*
 * Data structures
 */

/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP**          decdecomps;                        /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;                         /**< array of structure detectors */
   int*                  priorities;                        /**< priorities of the detectors */
 //  std::vector<SCIP_HASHMAP*> initalpartialdecomps;                      /**< possible incomplete decompositions given by user */
   int                   ndetectors;                        /**< number of detectors */
   SCIP_CLOCK*           detectorclock;                     /**< clock to measure detection time */
   SCIP_Bool             hasrun;                            /**< flag to indicate whether we have already detected */
   int                   ndecomps;                          /**< number of decomposition structures  */
   int                   sizedecomps;                       /**< size of the decomp and complete seeeds array */
   int                   sizeincompleteseeeds;              /**< size of the incomplete seeeds array */
   int                   maxndetectionrounds;               /**< maximum number of detection loop rounds  */
   int                   weightinggpresolvedoriginaldecomps; /**< weighing method for comparing presovled and original decompositions (see corresponding enum)   */
   SCIP_Bool             createbasicdecomp;                 /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
   SCIP_Bool             enableorigdetection;               /**< indicates whether to start detection for the original problem */
   SCIP_Bool             conssclassnnonzenabled;            /**< indicates whether constraint classifier for nonzero entries is enabled */
   SCIP_Bool             conssclassnnonzenabledorig;        /**< indicates whether constraint classifier for nonzero entries is enabled for the original problem */
   SCIP_Bool             conssclassnconstypeenabled;        /**< indicates whether constraint classifier for scipconstype is enabled */
   SCIP_Bool             conssclassnconstypeenabledorig;    /**< indicates whether constraint classifier for scipconstype is enabled for the original problem */
   SCIP_Bool             consnamenonumbersenabled;          /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled */
   SCIP_Bool             consnamenonumbersenabledorig;      /**< indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem */
   SCIP_Bool             conssclasslevenshteinabled;        /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled */
   SCIP_Bool             conssclasslevenshteinenabledorig;  /**< indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem */
   SCIP_Bool             varclassvartypesenabled;           /**< indicates whether variable classifier for scipvartypes is enabled */
   SCIP_Bool             varclassvartypesenabledorig;       /**< indicates whether variable classifier for scipvartypes is enabled for the original problem */

   int**                 candidatesNBlocks;                 /**< pointer to store candidates for number of blocks calculated by the seeedpool */
   int*                  nCandidates;
   SCIP_HASHMAP*         consToIndex;                       /**< hashmap from constraints to indices, to be filled */
   int*                  nConss;

   gcg::Seeedpool*		 seeedpool;                         /** seeedpool that manages the detection  process for the presolved transformed problem */
   gcg::Seeedpool*       seeedpoolunpresolved;              /** seeedpool that manages the deetction of the unpresolved problem */
   SeeedPtr*             allrelevantfinishedseeeds;         /** collection  of all relevant seeeds ( i.e. all seeeds w.r.t. copies ) */
   SeeedPtr*             incompleteseeeds;                  /** collection of incomplete seeeds originatging from incomplete decompostions given by the users */
   int                   nallrelevantseeeds;                /** number  of all relevant seeeds ( i.e. all seeeds w.r.t. copies ) */
   int                   nincompleteseeeds;                 /** number  of incomplete seeeds originatging from incomplete decompostions given by the users */
   SCIP_HASHMAP*         seeedtodecdecomp;                  /**< hashmap from seeeds to the corresponding decdecomp (or NULL if the seeed is incomplete)  */
   SCIP_HASHMAP*         decdecomptoseeed;                  /**< hashmap from decompositions to the corresponding seeed */

   SeeedPtr              curruserseeed;
   SCIP_Bool             unpresolveduserseeedadded;         /**< stores whether or not an unpresolved user seeed was added */

   /** new data fields for selection management */
   int                    startidvisu;                       /** when displaying the list of decomps, this is the starting index */
   int                    selectvisulength;                  /** number of decompositions to be displayed at once */
   std::vector<SeeedPtr>* listall;                           /** vector containing the current list of decomps to visualize*/
   std::vector<int>*      selected;                          /** vector containing the indices of selected decompositions */
   SCIP_Bool              selectedexists;                    /** are there some selected decompositions */

   int                    seeedcounter;                      /** counts the number of seeeds, used for seeed ids */
   int                    helpvisucounter;                   /** help counter for family tree visualization to iterate the heights */

   std::vector<std::pair<SeeedPtr, SCIP_Real> >* candidates;

};

enum weightinggpresolvedoriginaldecomps{
   NO_MODIF = 0,
   FRACTION_OF_NNONZEROS,
   FRACTION_OF_NROWS,
   FAVOUR_PRESOLVED
};


/*
 * Local methods
 */

std::string getSeeedFolderLatex( SeeedPtr seeed )
{
   std::stringstream decompfilename;
   decompfilename << "dec" << seeed->getID() << ".pdf";

   return decompfilename.str();
}


SCIP_Bool unfinishedchildexists(std::vector<SCIP_Bool> const& childsfinished)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return true;
   }
   return false;
}

int getfirstunfinishedchild(std::vector<SCIP_Bool> const& childsfinished, std::vector<int> const& childs)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return childs[s];
   }
   return -1;
}

int getfirstunfinishedchildid(std::vector<SCIP_Bool> const& childsfinished, std::vector<int> const& childs)
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
         return (int)s;
   }
   return -1;
}


/**
 * @return is nextchild the last unfinished child
 */
SCIP_Bool finishnextchild( std::vector<int>& childs, std::vector<SCIP_Bool>& childsfinished, int child )
{
   for( size_t s = 0; s < childsfinished.size(); ++s )
   {
      if( !childsfinished[s] )
      {
         assert(childs[s] == child);
         childsfinished[s] = TRUE;
         return s == childsfinished.size() - 1;
      }
   }
   return FALSE;
}

std::string writeSeeedDetectorChainInfoLatex( SeeedPtr seeed, int currheight, int visucounter )
{
   std::stringstream line;
   std::string relposition;
   int position = visucounter % 3;
   if( position == 0 )
      relposition = "above";
   else if ( position == 1)
      relposition = "";
   else if ( position == 2)
      relposition = "below";
   else
      relposition = "below left";

   if ( currheight != 1)
      relposition = "";

   if ( (size_t) currheight >  seeed->detectorchaininfo.size() )
      line << "edge from parent node [" << relposition << "] {no info" << seeed->getID() << "-" << currheight -1 << " } " ;
   else
   {
      std::string oldinfo = seeed->detectorchaininfo[ currheight - 1];
      /** take latexified detctorchaininfo */
      size_t index = 0;
      while (true) {
         /* Locate the substring to replace. */
         index = oldinfo.find("_", index);
         if (index == std::string::npos)
            break;
         if ( index > 0 &&   oldinfo.at(index-1) == '\\' )
         {
            ++index;
            continue;
         }

         /* Make the replacement. */
         oldinfo.replace(index, 1, "\\_");

         /* Advance index forward so the next iteration doesn't pick it up as well. */
         index += 2;
      }
      std::cout << "oldinfo: " << oldinfo << std::endl;

      line << "edge from parent node [" << relposition << "] {" << oldinfo <<"} " ;
   }


   return line.str();
}


std::string writeSeeedInfoLatex( SeeedPtr seeed )
{
   std::stringstream line;
   line << "\\node[below = \\belowcaptionskip of s" << seeed->getID() << "] (caps" << seeed->getID() << ") {\\scriptsize " << seeed->getShortCaption() << "}; " << std::endl;

   return line.str();
}


std::string writeSeeedIncludeLatex( SeeedPtr seeed, std::string workfolder )
{
   std::stringstream line;
   line << " (s" << seeed->getID() << ") { \\includegraphics[width=0.15\\textwidth]{" << getSeeedFolderLatex(seeed) << "} }" << std::endl;

   return line.str();
}


/** local method to handle storage of finished decompositions and corresponding seeeds */
static
SCIP_RETCODE SCIPstoreSeeedAndDecomp(
  SCIP*                 scip,
  SeeedPtr              seeed,
  DEC_DECOMP*           dec
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR*     conshdlr;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /** check if reallocation is needed */
   if ( conshdlrdata->ndecomps == conshdlrdata->sizedecomps )
   {
      conshdlrdata->sizedecomps  = SCIPcalcMemGrowSize(scip, conshdlrdata->sizedecomps + 1);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->decdecomps, (size_t) conshdlrdata->sizedecomps ) );
      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->allrelevantfinishedseeeds, (size_t) conshdlrdata->sizedecomps ) );
   }

   conshdlrdata->allrelevantfinishedseeeds[conshdlrdata->ndecomps] = seeed;
   conshdlrdata->decdecomps[conshdlrdata->ndecomps] = dec;
   ++conshdlrdata->ndecomps;

   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->seeedtodecdecomp, seeed, (void*) dec ) );
   SCIP_CALL( SCIPhashmapInsert(conshdlrdata->decdecomptoseeed,  dec, (void*) seeed ) );


   assert(SCIPhashmapExists( conshdlrdata->decdecomptoseeed,  dec) );

   return SCIP_OKAY;
}

/** local method to handle storage of finished decompositions and corresponding seeeds */
static
SCIP_RETCODE SCIPstoreIncompleteSeeed(
  SCIP*                 scip,
  SeeedPtr              seeed
)
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR*     conshdlr;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /** check if reallocation is needed */
   if ( conshdlrdata->nincompleteseeeds == conshdlrdata->sizeincompleteseeeds )
   {
      conshdlrdata->sizeincompleteseeeds  = SCIPcalcMemGrowSize(scip, conshdlrdata->sizeincompleteseeeds + 1);

      SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->incompleteseeeds, (size_t) conshdlrdata->sizeincompleteseeeds ) );
   }

   conshdlrdata->incompleteseeeds[conshdlrdata->nincompleteseeeds] = seeed;

   ++conshdlrdata->nincompleteseeeds;

   return SCIP_OKAY;
}


struct sort_pred {
    bool operator()(const std::pair<SeeedPtr, SCIP_Real> &left, const std::pair<SeeedPtr, SCIP_Real> &right) {
        return left.second < right.second;
    }
};


/**
 * create a 'decomposition' consisting of only one single block; used if no other decomposition was found
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

   return SCIP_OKAY;
}

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

   if( conshdlrdata->ndecomps > 0 )
   {
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
      conshdlrdata->decdecomps = NULL;
      conshdlrdata->ndecomps = 0;
   }

   if( conshdlrdata->nallrelevantseeeds > 0 )
   {
      for( i = 0; i < conshdlrdata->nallrelevantseeeds; ++i )
      {
         delete conshdlrdata->allrelevantfinishedseeeds[i];
      }

   }

   if( conshdlrdata->allrelevantfinishedseeeds != NULL )
      SCIPfreeMemoryArray(scip, &conshdlrdata->allrelevantfinishedseeeds);
   conshdlrdata->allrelevantfinishedseeeds = NULL;
   conshdlrdata->nallrelevantseeeds = 0;

   conshdlrdata->hasrun = FALSE;

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      DEC_DETECTOR *detector;
      detector = conshdlrdata->detectors[i];
      assert(detector != NULL);

      detector->ndecomps = 0;
      SCIPfreeMemoryArrayNull(scip, &detector->decomps);
      if( detector->exitDetector != NULL )
      {
         SCIPdebugMessage("Calling exitDetector of %s\n", detector->name);
         SCIP_CALL( (*detector->exitDetector)(scip, detector) );
      }
   }


   delete conshdlrdata->seeedpool;

   conshdlrdata->seeedpool = NULL;
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
   if( conshdlrdata->ndecomps > 0 )
   {
      for( i = 0; i < conshdlrdata->ndecomps; ++i )
      {
         SCIP_CALL( DECdecompFree(scip, &conshdlrdata->decdecomps[i]) );
      }
      SCIPfreeMemoryArray(scip, &conshdlrdata->decdecomps);
   }

   if( conshdlrdata->seeedpool != NULL )
      delete conshdlrdata->seeedpool;

   if( conshdlrdata->seeedpoolunpresolved != NULL )
      delete conshdlrdata->seeedpoolunpresolved;

   if( conshdlrdata->candidates != NULL )
         delete &conshdlrdata->candidates;

   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->incompleteseeeds );
   SCIPfreeMemoryArrayNull( scip, &conshdlrdata->allrelevantfinishedseeeds) ;
   SCIPfreeMemoryArray(scip, &conshdlrdata->priorities);
   SCIPfreeMemoryArray(scip, &conshdlrdata->detectors);
   SCIPfreeMemoryArrayNull(scip, &conshdlrdata->decdecomps);

   delete &conshdlrdata->selected;
   delete &conshdlrdata->listall;

   SCIPfreeMemory(scip, &conshdlrdata);

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

/** creates the handler for decomp constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrDecomp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   /* create decomp constraint handler data */
   SCIP_CALL( SCIPallocMemory(scip, &conshdlrdata) );
   assert(conshdlrdata != NULL);

   conshdlrdata->decdecomps = NULL;
   conshdlrdata->ndecomps = 0;
   conshdlrdata->ndetectors = 0;
   conshdlrdata->priorities = NULL;
   conshdlrdata->detectors = NULL;
   conshdlrdata->hasrun = FALSE;
   conshdlrdata->maxndetectionrounds = 0;
   conshdlrdata->enableorigdetection = FALSE;
   conshdlrdata->seeedpoolunpresolved = NULL;
   conshdlrdata->seeedpool = NULL;
   conshdlrdata->allrelevantfinishedseeeds = NULL;
   conshdlrdata->incompleteseeeds = NULL;
   conshdlrdata->nallrelevantseeeds = 0;
   conshdlrdata->nincompleteseeeds = 0;
   conshdlrdata->decdecomptoseeed = NULL;
   conshdlrdata->seeedtodecdecomp = NULL;
   conshdlrdata->curruserseeed = NULL;
   conshdlrdata->unpresolveduserseeedadded = FALSE;
   conshdlrdata->startidvisu = 0;
   conshdlrdata->selectvisulength = 10;
   conshdlrdata->listall = new std::vector<SeeedPtr>(0, NULL);
   conshdlrdata->selected = new std::vector<int>(0, -1);
   conshdlrdata->candidates = new std::vector<std::pair<SeeedPtr, SCIP_Real > >(0);
   conshdlrdata->selectedexists = FALSE;
   conshdlrdata->sizedecomps = 10;
   conshdlrdata->sizeincompleteseeeds = 10;
   conshdlrdata->seeedcounter = 0;

   SCIP_CALL( SCIPallocMemoryArray( scip, &conshdlrdata->decdecomps, conshdlrdata->sizedecomps) );
   SCIP_CALL( SCIPallocMemoryArray( scip, &conshdlrdata->allrelevantfinishedseeeds, conshdlrdata->sizedecomps) );
   SCIP_CALL( SCIPallocMemoryArray( scip, &conshdlrdata->incompleteseeeds, conshdlrdata->sizeincompleteseeeds) );

   SCIP_CALL( SCIPcreateWallClock(scip, &conshdlrdata->detectorclock) );

   SCIP_CALL( SCIPhashmapCreate( &conshdlrdata->decdecomptoseeed, SCIPblkmem(scip), MAXNDECOMPS ) );
   SCIP_CALL( SCIPhashmapCreate( &conshdlrdata->seeedtodecdecomp, SCIPblkmem(scip), MAXNDECOMPS ) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpDecomp, consEnfopsDecomp, consCheckDecomp, consLockDecomp,
         conshdlrdata) );
   assert(conshdlr != FALSE);

   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeDecomp) );
   SCIP_CALL( SCIPsetConshdlrInit(scip, conshdlr, consInitDecomp) );
   SCIP_CALL( SCIPsetConshdlrExit(scip, conshdlr, consExitDecomp) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/decomp/createbasicdecomp", "indicates whether to create a decomposition with all constraints in the master if no other specified", &conshdlrdata->createbasicdecomp, FALSE, DEFAULT_CREATEBASICDECOMP, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/origprob/enabled", "indicates whether to start detection for the original problem", &conshdlrdata->enableorigdetection, FALSE, DEFAULT_ENABLEORIGDETECTION, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", "indicates whether constraint classifier for nonzero entries is enabled", &conshdlrdata->conssclassnnonzenabled, FALSE, DEFAULT_CONSSCLASSNNONZENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/nnonzeros/origenabled", "indicates whether constraint classifier for nonzero entries is enabled for the original problem", &conshdlrdata->conssclassnnonzenabledorig, FALSE, DEFAULT_CONSSCLASSNNONZENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/scipconstype/enabled", "indicates whether constraint classifier for scipconstype is enabled", &conshdlrdata->conssclassnnonzenabled, FALSE, DEFAULT_CONSSCLASSSCIPCONSTYPEENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/scipconstype/origenabled", "indicates whether constraint classifier for scipconsstype is enabled for the original problem", &conshdlrdata->conssclassnnonzenabledorig, FALSE, DEFAULT_CONSSCLASSSCIPCONSTYPEENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", "indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled", &conshdlrdata->consnamenonumbersenabled, FALSE, DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamenonumbers/origenabled", "indicates whether constraint classifier for constraint names (remove digits; check for identity) is enabled for the original problem", &conshdlrdata->consnamenonumbersenabledorig, FALSE, DEFAULT_CONSSCLASSCONSNAMENONUMBERENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", "indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled", &conshdlrdata->conssclasslevenshteinabled, FALSE, DEFAULT_CONSSCLASSLEVENSHTEINENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/consclassifier/consnamelevenshtein/origenabled", "indicates whether constraint classifier for constraint names (according to levenshtein distance graph) is enabled for the original problem", &conshdlrdata->conssclasslevenshteinenabledorig, FALSE, DEFAULT_CONSSCLASSLEVENSHTEINENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/scipvartype/enabled", "indicates whether variable classifier for scipvartypes is enabled", &conshdlrdata->varclassvartypesenabled, FALSE, DEFAULT_VARCLASSSCIPVARTYPESENABLED, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detection/varclassifier/scipvartype/origenabled", "indicates whether variable classifier for scipvartypes is enabled for the original problem", &conshdlrdata->varclassvartypesenabledorig, FALSE, DEFAULT_VARCLASSSCIPVARTYPESENABLEDORIG, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/maxrounds",
      "Maximum number of detection loop rounds", &conshdlrdata->maxndetectionrounds, FALSE,
      DEFAULT_MAXDETECTIONROUNDS, 0, INT_MAX, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detection/origprob/weightinggpresolvedoriginaldecomps",
      "Weighting method when comparing decompositions for presolved and unpresolved problem", &conshdlrdata->weightinggpresolvedoriginaldecomps, TRUE,
      NO_MODIF, 0, 3, NULL, NULL) );


   assert(conshdlrdata->candidates != NULL);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDecompShowListExtractHeader(
   SCIP*                   scip
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

   size_t i;

   ndetectedpresolved = 0;
   ndetectedunpresolved = 0;
   nuserpresolvedfull = 0;
   nuserpresolvedpartial = 0;
   nuserunpresolvedfull = 0;
   nuserunpresolvedpartial = 0;


   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   /** count corresponding seeeds */
   for ( i = 0; i < conshdlrdata->listall->size(); ++i )
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->listall->at(i);
      if( seeed->isComplete() && seeed->usergiven == gcg::USERGIVEN::NOT && !seeed->isfromunpresolved )
         ++ndetectedpresolved;
      if( seeed->isComplete() && seeed->usergiven == gcg::USERGIVEN::NOT && seeed->isfromunpresolved )
         ++ndetectedunpresolved;
      if( seeed->isComplete() && ( seeed->usergiven == gcg::USERGIVEN::COMPLETE || seeed->usergiven == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && !seeed->isfromunpresolved )
         ++nuserpresolvedfull;
      if( !seeed->isComplete() && seeed->usergiven == gcg::USERGIVEN::PARTIAL && !seeed->isfromunpresolved )
         ++nuserpresolvedpartial;
      if( seeed->isComplete() && ( seeed->usergiven == gcg::USERGIVEN::COMPLETE || seeed->usergiven == gcg::USERGIVEN::COMPLETED_CONSTOMASTER) && seeed->isfromunpresolved )
         ++nuserunpresolvedfull;
      if( !seeed->isComplete() && seeed->usergiven == gcg::USERGIVEN::PARTIAL && seeed->isfromunpresolved )
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
   SCIPdialogMessage(scip, NULL, "   id   nbloc  nmacon  nlivar  nmavar  nstlva  maxwhi  history  pre  nopcon  nopvar  usr  sel \n");
   SCIPdialogMessage(scip, NULL, " ----   -----  ------  ------  ------  ------  ------  -------  ---  ------  ------  ---  --- \n");



   return SCIP_OKAY;
}


/** sets (and adds) the decomposition structure;
 * this method should only be called if there is no seeed for this decomposition
 *
 * **/
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

      assert( seeed->checkConsistency(conshdlrdata->seeedpool) ); // TODO correct?

      SCIPdialogMessage(scip, NULL, " %4d   ", i );
      SCIPdialogMessage(scip, NULL, "%5d  ", seeed->getNBlocks() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMasterconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNLinkingvars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNMastervars() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNTotalStairlinkingvars() );
      if( seeed->isComplete() )
         SCIPdialogMessage(scip, NULL, "%.4f  ", 1. - seeed->getMaxWhiteScore() );
      else
         SCIPdialogMessage(scip, NULL, "<=%.2f  ", 1. - seeed->getMaxWhiteScore() );
      SCIPdialogMessage(scip, NULL, "%7s  ", seeed->detectorchainstring );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->isfromunpresolved ? "no" : "yes")  );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenconss() );
      SCIPdialogMessage(scip, NULL, "%6d  ", seeed->getNOpenvars() );
      SCIPdialogMessage(scip, NULL, "%3s  ", (seeed->usergiven == gcg::USERGIVEN::NOT ? "no" : "yes")   );
      SCIPdialogMessage(scip, NULL, "%3s  \n", (seeed->isSelected() ? "yes" : "no")  );
   }

   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

   return SCIP_OKAY;
}


/** sets (and adds) the decomposition structure;
 * this method should only be called if there is no seeed for this decomposition
 *
 * **/
SCIP_RETCODE SCIPconshdlrDecompAddDecdecomp(
   SCIP*                 scip,               /**< SCIP data structure */
   DEC_DECOMP*           decdecomp           /**< DEC_DECOMP data structure */
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

   SCIP_CALL( conshdlrdata->seeedpool->createSeeedFromDecomp(decdecomp, &seeed) );

   SCIP_CALL( SCIPstoreSeeedAndDecomp(scip, seeed, decdecomp) );

   return SCIP_OKAY;
}

///** sets (and adds) the decomposition structure **/
//SCIP_RETCODE SCIPconshdlrDecompAddConsToBlock(
//   SCIP*                 scip,               /**< SCIP data structure */
//   SCIP_HASHMAP*         consToBlock,        /**< possible incomplete detection info */
//   SCIP_HASHMAP*         varsToBlock        /**< possible incomplete detection info stored as two hashmaps*/
//   )
//{
//   SCIP_CONSHDLR* conshdlr;
//   SCIP_CONSHDLRDATA* conshdlrdata;
//   assert(scip != NULL);
//   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
//   assert( conshdlr != NULL );
//
//   conshdlrdata = SCIPconshdlrGetData(conshdlr);
//   assert(conshdlrdata != NULL);
//
//   return SCIP_ERROR;
//
////   conshdlrdata->initalpartialdecomps.push_back(consToBlock);
//
//      }

SCIP_RETCODE SCIPconshdlrDecompShowLegend(
   SCIP* scip
   )
{
   int det;

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

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

//   SCIPdialogMessage(scip, NULL, "   id   nbloc  nmacon  nlivar  nmavar  nstlva  maxwhi  history  pre  nopcon  nopvar"
//      "  usr"
//      "  sel \n");
//   SCIPdialogMessage(scip, NULL, " ----   -----  ------  ------  ------  ------  ------  -------  ---  ------  ------"
//      "  ---"
//      "  --- \n");

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
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "maxwhi", "maximum white area score (i.e. maximize fraction of white area score; white area is nonblock and nonborder area, stairlinking variables count as linking)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "history", "list of detector chars worked on this decomposition ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "pre", "is this decomposition for the presolved problem");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopcon", "number of open constraints");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "nopvar", "number of open variables");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "usr", "was this decomposition given by the user");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "sel", "is this decomposition selected at the moment");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");
   return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDecompModifyNVisualized(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
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

   SCIPdialogMessage(scip, NULL, "Please specify the maximum number of decompositions displayed at once in the table [%d]:\n", conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   newval = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      newval = atoi(ntovisualize);

   if (newval != 0)
      conshdlrdata->selectvisulength = newval;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDecompSelectVisualize(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
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

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be visualized:\n", conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   gcg::Seeedpool* seeedpool = (conshdlrdata->listall->at(idtovisu)->isfromunpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );
   conshdlrdata->listall->at(idtovisu)->showScatterPlot(seeedpool);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDecompSelectSelect(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   char* ntovisualize;
   SCIP_Bool endoffile;
   int idtovisu;
   SeeedPtr toselect;
   gcg::Seeedpool* seeedpool;

   int commandlen;

   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIPdialogMessage(scip, NULL, "Please specify the id of the decomposition to be visualized:\n", conshdlrdata->selectvisulength );
   SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, " ", &ntovisualize, &endoffile) );
   commandlen = strlen(ntovisualize);

   idtovisu = conshdlrdata->selectvisulength;
   if( commandlen != 0)
      idtovisu = atoi(ntovisualize);

   seeedpool = (conshdlrdata->listall->at(idtovisu)->isfromunpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool );
   toselect = conshdlrdata->listall->at(idtovisu);

   toselect->setSelected(!toselect->isSelected() );

   if( !toselect->isSelected() )
   {
      conshdlrdata->selected->erase(  find( conshdlrdata->selected->begin(), conshdlrdata->selected->end(), idtovisu) );
   }
   else
   {
      conshdlrdata->selected->push_back(idtovisu);
   }

   conshdlrdata->selectedexists = (conshdlrdata->selected->size() > 0);

   return SCIP_OKAY;
}


SCIP_RETCODE SCIPconshdlrDecompShowHelp(
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


   SCIPdialogMessage(scip, NULL, "============================================================================================= \n");

//   SCIPdialogMessage(scip, NULL, "   id   nbloc  nmacon  nlivar  nmavar  nstlva  maxwhi  history  pre  nopcon  nopvar"
//      "  usr"
//      "  sel \n");
//   SCIPdialogMessage(scip, NULL, " ----   -----  ------  ------  ------  ------  ------  -------  ---  ------  ------"
//      "  ---"
//      "  --- \n");

   SCIPdialogMessage(scip, NULL, "\n" );

   SCIPdialogMessage(scip, NULL, "List of selection commands \n" );
   SCIPdialogMessage(scip, NULL, "\n" );
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "command", "description");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "-------", "-----------");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "select", "selects/unselects decomposition with given id");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "back", "displays the preceding decompositions (if there are some)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "next", "displays the subsequent decompositions (if there are some)");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "top", "displays the first decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "end", "displays the last decompositions");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "legend", "displays the legend for table header and history abbreviations");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "help", "displays this help");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "modify", "modifies the number of displayed decompositions ");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "quit", "finishes selection and goes back to main menu");
   SCIPdialogMessage(scip, NULL, "%30s     %s\n", "visualize", "experimental feature: visualizes the specified decomposition ");

   SCIPdialogMessage(scip, NULL, "\n============================================================================================= \n");


   return SCIP_OKAY;
}

SCIP_Bool SCIPconshdlrDecompIsBestCandidateUnpresolved(
   SCIP*                   scip
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;


   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   return conshdlrdata->candidates->at(0).first->isfromunpresolved;
}


SCIP_RETCODE SCIPconshdlrDecompExecSelect(
   SCIP*                   scip,
   SCIP_DIALOGHDLR*        dialoghdlr,
   SCIP_DIALOG*            dialog )
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

   /** 1) update list of interesting seeeds */

   SCIP_CALL( SCIPconshdlrDecompUpdateSeeedlist(scip) );


   /** 2) while user has not aborted: show current list extract */

   while ( !finished )
   {
      int commandlen;

      SCIP_CALL( SCIPconshdlrDecompShowListExtractHeader(scip) );

      SCIP_CALL( SCIPconshdlrDecompShowListExtract(scip) );



      SCIP_CALL( SCIPdialoghdlrGetWord(dialoghdlr, dialog, "Please enter selection command or decomposition id to select (or \"h\" for help) : \nGCG/select> ", &command, &endoffile) );

      commandlen = strlen(command);

      /** case distinction: */
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
         if(conshdlrdata->startidvisu > conshdlrdata->listall->size() - conshdlrdata->selectvisulength )
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

      if( strncmp( command, "modify", commandlen) == 0 )
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

      if( strncmp( command, "select", commandlen) == 0 )
      {
         SCIP_CALL(SCIPconshdlrDecompSelectSelect(scip, dialoghdlr, dialog ) );
         continue;
      }
   }

   return SCIP_OKAY;
}



/** returns the decomposition structure **/
DEC_DECOMP** SCIPconshdlrDecompGetDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->decdecomps;
}

/** returns the decomposition structure **/
int SCIPconshdlrDecompGetNDecdecomps(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->ndecomps;
}

/** returns the data of the provided detector */
DEC_DETECTORDATA* DECdetectorGetData(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->decdata;

}


/** returns the seeedpool **/
gcg::Seeedpool* SCIPconshdlrDecompGetSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->seeedpool;
}

/** returns the seeedpool for the unpresolved problem **/
gcg::Seeedpool* SCIPconshdlrDecompGetSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   assert(scip != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   return conshdlrdata->seeedpoolunpresolved;
}

/** creates the seeedpool **/
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpool(
   SCIP*                 scip                /**< SCIP data structure */
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
      conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE);


   return SCIP_OKAY;
}

/** creates the seeedpool **/
SCIP_RETCODE SCIPconshdlrDecompCreateSeeedpoolUnpresolved(
   SCIP*                 scip                /**< SCIP data structure */
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
      conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE);

   return SCIP_OKAY;
}



/** returns the name of the provided detector */
const char* DECdetectorGetName(
   DEC_DETECTOR*         detector            /**< detector data structure */
   )
{
   assert(detector != NULL);
   return detector->name;
}

/** searches for the detector and returns it or returns NULL if detector is not found*/
DEC_DETECTOR* DECfindDetector(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           name                /**< name of the detector */
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

/** includes the detector */
SCIP_RETCODE DECincludeDetector(
   SCIP*                 scip,                   /**< SCIP data structure */
   const char*           name,                   /**< name of the detector */
   const char            decchar,                /**< display character of the detector */
   const char*           description,            /**< description of the detector */
   int                   freqCallRound,          /** frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 */
   int                   maxCallRound,           /** last round the detector gets called                              */
   int                   minCallRound,           /** first round the detector gets called (offset in detection loop) */
   int                   freqCallRoundOriginal,  /** frequency the detector gets called in detection loop while detecting of the original problem */
   int                   maxCallRoundOriginal,   /** last round the detector gets called while detecting of the original problem */
   int                   minCallRoundOriginal,   /** first round the detector gets called (offset in detection loop) while detecting of the original problem */
   int                   priority,               /**< priority of the detector                                           */
   SCIP_Bool             enabled,                /**< whether the detector should be enabled by default                  */
   SCIP_Bool             enabledOriginal,        /**< whether the detector should be enabled by default for detecting the original problem */
   SCIP_Bool             enabledFinishing,       /**< whether the finishing should be enabled */
   SCIP_Bool             skip,                   /**< whether the detector should be skipped if others found structure   */
   SCIP_Bool             usefulRecall,           /** is it useful to call this detector on a descendant of the propagated seeed */
   DEC_DETECTORDATA*     detectordata,           /**< the associated detector data (or NULL) */
   DEC_DECL_DETECTSTRUCTURE((*detectStructure)), /**< the method that will detect the structure (must not be NULL)*/
   DEC_DECL_FREEDETECTOR((*freeDetector)),       /**< destructor of detector (or NULL) */
   DEC_DECL_INITDETECTOR((*initDetector)),       /**< initialization method of detector (or NULL) */
   DEC_DECL_EXITDETECTOR((*exitDetector)),       /**< deinitialization method of detector (or NULL) */
   DEC_DECL_PROPAGATESEEED((*propagateSeeedDetector)),
   DEC_DECL_FINISHSEEED((*finishSeeedDetector)),
   DEC_DECL_SETPARAMAGGRESSIVE((*setParamAggressiveDetector)),
   DEC_DECL_SETPARAMDEFAULT((*setParamDefaultDetector)),
   DEC_DECL_SETPARAMFAST((*setParamFastDetector))
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
   assert(detectStructure != NULL);
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   assert(detectStructure != NULL);

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
   detector->finishSeeed = finishSeeedDetector;
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
   detector->skip = skip;
   detector->usefulRecall = usefulRecall;
   detector->ndecomps = 0;
   detector->decomps = NULL;
   detector->dectime = 0.;

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/enabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabled), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origenabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled for detecting in the original problem", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabledOrig), FALSE, enabled, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> is enabled for finishing of incomplete decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->enabledFinishing), FALSE, enabledFinishing, NULL, NULL) );


   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/skip", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be skipped if others found decompositions", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->skip), FALSE, skip, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/usefullrecall", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "flag to indicate whether detector <%s> should be called on descendants of the current seeed", name);
   SCIP_CALL( SCIPaddBoolParam(scip, setstr, descstr, &(detector->usefulRecall), FALSE, usefulRecall, NULL, NULL) );


   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/freqcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "frequency the detector gets called in detection loop ,ie it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->freqCallRound), FALSE, freqCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/maxcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "maximum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->maxCallRound), FALSE, maxCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/mincallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "minimum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->minCallRound), FALSE, minCallRound, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origfreqcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "frequency the detector gets called in detection loop,i.e., it is called in round r if and only if minCallRound <= r <= maxCallRound AND  (r - minCallRound) mod freqCallRound == 0 <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->freqCallRoundOriginal), FALSE, freqCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origmaxcallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "maximum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->maxCallRoundOriginal), FALSE, maxCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/origmincallround", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "minimum round the detector gets called in detection loop <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->minCallRoundOriginal), FALSE, minCallRoundOriginal, 0, INT_MAX, NULL, NULL) );

   (void) SCIPsnprintf(setstr, SCIP_MAXSTRLEN, "detectors/%s/priority", name);
   (void) SCIPsnprintf(descstr, SCIP_MAXSTRLEN, "priority of detector <%s>", name);
   SCIP_CALL( SCIPaddIntParam(scip, setstr, descstr, &(detector->priority), FALSE, priority, INT_MIN, INT_MAX, NULL, NULL) );


   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->detectors, (size_t)conshdlrdata->ndetectors+1) );
   SCIP_CALL( SCIPreallocMemoryArray(scip, &conshdlrdata->priorities,(size_t) conshdlrdata->ndetectors+1) );

   conshdlrdata->detectors[conshdlrdata->ndetectors] = detector;
   conshdlrdata->ndetectors = conshdlrdata->ndetectors+1;

   return SCIP_OKAY;

}

/** returns the remaining time of scip that the decomposition may use */
SCIP_Real DECgetRemainingTime(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_Real timelimit;
   assert(scip != NULL);
   SCIP_CALL_ABORT(SCIPgetRealParam(scip, "limits/time", &timelimit));
   if( !SCIPisInfinity(scip, timelimit) )
      timelimit -= SCIPgetSolvingTime(scip);
   return timelimit;
}

/** creates a user seeed for the presolved problem **/
SCIP_RETCODE SCIPconshdlrDecompCreateUserSeeed(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_Bool             presolved           /**< should the user seeed be created for the presolved problem */
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

   conshdlrdata->curruserseeed = new gcg::Seeed(scip, currseeedpool->getNewIdForSeeed(), currseeedpool->getNDetectors(), currseeedpool->getNConss(), currseeedpool->getNVars() );


   conshdlrdata->curruserseeed->stemsFromUnpresolved = !presolved;

   return SCIP_OKAY;
}

SCIP_Bool SCIPconshdlrDecompUnpresolvedUserSeeedAdded(
   SCIP*                 scip                /**< SCIP data structure */
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

   return conshdlrdata->unpresolveduserseeedadded;
}

SCIP_RETCODE SCIPconshdlrdataDecompUnselectAll(
   SCIP*          scip
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   size_t i;

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

SCIP_RETCODE   SCIPconshdlrDecompPopulateSelected(
   SCIP*       scip
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   std::vector<SeeedPtr>  unfinishedunpresolved(0);
   std::vector<SeeedPtr>  unfinishedpresolved(0);
   SCIP_Bool selectedexists;
   size_t i;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   assert( SCIPconshdlrDecompCheckConsistency(scip) );

   selectedexists = SCIPconshdlrDecompExistsSelected(scip);

   /** check existence of seeedpools */
   if( conshdlrdata->seeedpoolunpresolved == NULL)
      SCIPconshdlrDecompCreateSeeedpoolUnpresolved(scip);

   if( conshdlrdata->seeedpool == NULL)
      SCIPconshdlrDecompCreateSeeedpool(scip);


   /** get unfinished unpresolved */
   for( i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->incompleteSeeeds.size() ; ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpoolunpresolved->incompleteSeeeds[i];
      seeed->isfromunpresolved = TRUE;

      if( seeed->isSelected() || (!selectedexists && seeed->usergiven != gcg::USERGIVEN::NOT && !seeed->isComplete() ) )
         unfinishedunpresolved.push_back(seeed);
   }

   /** enable orig detection if there are unpresolved decompostions */
   if( unfinishedunpresolved.size() > 0)
         SCIPsetBoolParam(scip, "detection/origprob/enabled", TRUE );


   /** get unfinished presolved */
   for( i = 0; i < (size_t) conshdlrdata->nincompleteseeeds ; ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->incompleteseeeds[i];

      if( seeed->isSelected() || (!selectedexists && seeed->usergiven != gcg::USERGIVEN::NOT && !seeed->isComplete() ) )
         unfinishedpresolved.push_back(seeed);
   }

   /** clear old seeds   */
   conshdlrdata->seeedpoolunpresolved->currSeeeds.clear();
   conshdlrdata->seeedpool->currSeeeds.clear();
   conshdlrdata->seeedpoolunpresolved->incompleteSeeeds.clear();
   conshdlrdata->seeedpool->incompleteSeeeds.clear();

   /** populate unpresolved */
   for( i = 0 ; i < unfinishedunpresolved.size() ; ++i)
   {
      conshdlrdata->seeedpoolunpresolved->populate(unfinishedunpresolved);
   }

   /** populate presolved */
   for( i = 0 ; i < unfinishedpresolved.size() ; ++i)
   {
      conshdlrdata->seeedpool->populate(unfinishedpresolved);
   }

   return SCIP_OKAY;
}




SCIP_RETCODE SCIPconshdlrDecompUpdateSeeedlist(
   SCIP*          scip
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



   /** sort decomposition and finished seeeds according to max white score */
   SCIP_CALL( DECconshdlrDecompSortDecompositionsByScore(scip) );

   /** add seeeds to list */
   /** 1) add presolved finished */
   for( i = 0; i < conshdlrdata->ndecomps ; ++i)
    {
       SeeedPtr seeed;
       seeed = conshdlrdata->allrelevantfinishedseeeds[i];

       conshdlrdata->listall->push_back(seeed);
    }

   /** 2) add presolved unfinished */
   for( i = 0; i < conshdlrdata->nincompleteseeeds ; ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->incompleteseeeds[i];

      conshdlrdata->listall->push_back(seeed);
   }

   /** 3) add unpresolved finished */
   for( i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() ; ++i)
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpoolunpresolved->getFinishedSeeed(i);
      seeed->isfromunpresolved = TRUE;

      conshdlrdata->listall->push_back(seeed);
   }

   /** 4) add unpresolved partial */
   for( i = 0; conshdlrdata->seeedpoolunpresolved != NULL && i < conshdlrdata->seeedpoolunpresolved->incompleteSeeeds.size() ; ++i )
   {
      SeeedPtr seeed;
      seeed = conshdlrdata->seeedpoolunpresolved->incompleteSeeeds[i];
      seeed->isfromunpresolved = TRUE;

      conshdlrdata->listall->push_back(seeed);
   }


   return SCIP_OKAY;
}


/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetnumberOfBlocks(
   SCIP*                 scip,                /**< SCIP data structure */
   int                   nblocks               /**< number of blocks */
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

/** returns whether there is an user seeed  */
SCIP_Bool SCIPconshdlrDecompUserSeeedIsActive(
   SCIP*                 scip                /**< SCIP data structure */
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

   return conshdlrdata->curruserseeed != NULL;
}


/** sets the number of blocks */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsDefaultMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   SCIP_Bool             consdefaulttomaster  /**< are not specified constraints set to master for default */
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


   conshdlrdata->curruserseeed->usergiven = gcg::USERGIVEN::COMPLETED_CONSTOMASTER;


   return SCIP_OKAY;
}



/** sets a constraint by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname,            /**< name of the constraint */
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

   currseeedpool = conshdlrdata->curruserseeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   cons = conshdlrdata->curruserseeed->stemsFromUnpresolved ? SCIPfindOrigCons(scip, consname ) : SCIPfindCons(scip, consname );
   consindex = currseeedpool->getIndexForCons( cons ) ;

   conshdlrdata->curruserseeed->bookAsBlockCons(consindex, blockid);

   return SCIP_OKAY;
}


/** sets a constraint by name to the master in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetConsToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           consname
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

   currseeedpool = conshdlrdata->curruserseeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;

   cons = conshdlrdata->curruserseeed->stemsFromUnpresolved ? SCIPfindOrigCons(scip, consname ) : SCIPfindCons(scip, consname );
   consindex = currseeedpool->getIndexForCons( cons );

   conshdlrdata->curruserseeed->bookAsMasterCons(consindex);

   return SCIP_OKAY;

}


/** sets a variable by name to a block in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToBlock(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname,             /**< name of the variable */
   int                   blockid              /**< block index ( counting from 0) */
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

   currseeedpool = conshdlrdata->curruserseeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

   conshdlrdata->curruserseeed->bookAsBlockVar(varindex, blockid);

   return SCIP_OKAY;
}


/** sets a variable by name to the master in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToMaster(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
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

   currseeedpool = conshdlrdata->curruserseeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
   varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

   conshdlrdata->curruserseeed->bookAsMasterVar(varindex);

   return SCIP_OKAY;

}


/** sets a variable by name to the linking variables in the current user seeed */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedSetVarToLinking(
   SCIP*                 scip,                /**< SCIP data structure */
   const char*           varname              /**< name of the variable */
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

      currseeedpool = conshdlrdata->curruserseeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;
      varindex = currseeedpool->getIndexForVar( SCIPfindVar(scip, varname ) );

      conshdlrdata->curruserseeed->bookAsLinkingVar(varindex);

      return SCIP_OKAY;

}


/** finalizes and flushes the current user seeed, i.e. consider implicits, calc hashvalue, construct decdecomp if complete etc */
SCIP_RETCODE SCIPconshdlrDecompUserSeeedFlush(
   SCIP*                 scip                 /**< SCIP data structure */
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

   currseeedpool = seeed->stemsFromUnpresolved ? conshdlrdata->seeedpoolunpresolved : conshdlrdata->seeedpool;

   seeed->flushBooked();

   if( seeed->shouldCompletedByConsToMaster() )
   {
      for( int opencons = 0; opencons < seeed->getNOpenconss(); ++opencons)
         seeed->bookAsMasterCons( seeed->getOpenconss()[opencons] );
      seeed->flushBooked();
   }

   currseeedpool->prepareSeeed(conshdlrdata->curruserseeed);


   if( conshdlrdata->curruserseeed->isComplete() )
   {
      if( !seeed->shouldCompletedByConsToMaster() )
         conshdlrdata->curruserseeed->usergiven = gcg::USERGIVEN::COMPLETE;
      /** stems from presolved problem? */
      if( !conshdlrdata->curruserseeed->stemsFromUnpresolved )
      {
         DEC_DECOMP* newdecomp;

         SCIP_CALL( conshdlrdata->seeedpool->createDecompFromSeeed(seeed, &newdecomp) );

         SCIP_CALL( SCIPstoreSeeedAndDecomp(scip, conshdlrdata->curruserseeed, newdecomp) );

      }
      /** stems from unpresolved problem */
      else
      {
         conshdlrdata->seeedpoolunpresolved->addSeeedToFinished(seeed);
         conshdlrdata->unpresolveduserseeedadded = TRUE;
      }

   }
   else
   {
      assert( !seeed->shouldCompletedByConsToMaster() );
      conshdlrdata->curruserseeed->usergiven = gcg::USERGIVEN::PARTIAL;

      if ( !conshdlrdata->curruserseeed->stemsFromUnpresolved )
         SCIP_CALL(SCIPstoreIncompleteSeeed(scip, conshdlrdata->curruserseeed) );
      else
         conshdlrdata->seeedpoolunpresolved->incompleteSeeeds.push_back(conshdlrdata->curruserseeed);

   }

   /** set statistics */


   {
      int nvarstoblock = 0;
      int nconsstoblock = 0;

      for ( int b = 0; b < conshdlrdata->curruserseeed->getNBlocks(); ++b )
      {
         nvarstoblock += conshdlrdata->curruserseeed->getNVarsForBlock(b);
         nconsstoblock += conshdlrdata->curruserseeed->getNConssForBlock(b);
      }
      conshdlrdata->curruserseeed->setDetectorPropagated(NULL);

      conshdlrdata->curruserseeed->detectorClockTimes.push_back(0.);
      conshdlrdata->curruserseeed->pctVarsFromFree.push_back( (nvarstoblock + conshdlrdata->curruserseeed->getNMastervars() +conshdlrdata->curruserseeed->getNLinkingvars())/(SCIP_Real) conshdlrdata->curruserseeed->getNVars()  );
      conshdlrdata->curruserseeed->pctVarsToBlock.push_back((nvarstoblock )/(SCIP_Real) conshdlrdata->curruserseeed->getNVars() );
      conshdlrdata->curruserseeed->pctVarsToBorder.push_back( (conshdlrdata->curruserseeed->getNMastervars() +conshdlrdata->curruserseeed->getNLinkingvars())/(SCIP_Real) conshdlrdata->curruserseeed->getNVars() ) ;
      conshdlrdata->curruserseeed->pctConssToBorder.push_back( (conshdlrdata->curruserseeed->getNMasterconss() ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() ) ;
      conshdlrdata->curruserseeed->pctConssFromFree.push_back( (conshdlrdata->curruserseeed->getNMasterconss() + nconsstoblock ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() ) ;
      conshdlrdata->curruserseeed->pctConssToBlock.push_back( (nconsstoblock ) / (SCIP_Real) conshdlrdata->curruserseeed->getNConss() );
      conshdlrdata->curruserseeed->nNewBlocks.push_back(conshdlrdata->curruserseeed->getNBlocks());
   }

   conshdlrdata->curruserseeed->findVarsLinkingToMaster(currseeedpool);
   conshdlrdata->curruserseeed->findVarsLinkingToStairlinking(currseeedpool);


   if( conshdlrdata->curruserseeed->usergiven == gcg::USERGIVEN::PARTIAL )
      usergiveninfo = "partial";
   if( conshdlrdata->curruserseeed->usergiven == gcg::USERGIVEN::COMPLETE )
      usergiveninfo = "complete";
   if( conshdlrdata->curruserseeed->usergiven == gcg::USERGIVEN::COMPLETED_CONSTOMASTER )
         usergiveninfo = "complete";
   if( conshdlrdata->curruserseeed->stemsFromUnpresolved )
         presolvedinfo = "unpresolved";
   else presolvedinfo = "presolved";

   conshdlrdata->curruserseeed->buildDecChainString();



   SCIPinfoMessage(scip, NULL, " added %s decomp for %s problem with %d blocks and %d masterconss, %d linkingvars, "
      "%d mastervars, and max white score of %s %f \n", usergiveninfo, presolvedinfo,
      conshdlrdata->curruserseeed->getNBlocks(), conshdlrdata->curruserseeed->getNMasterconss(),
      conshdlrdata->curruserseeed->getNLinkingvars(), conshdlrdata->curruserseeed->getNMastervars(), (conshdlrdata->curruserseeed->isComplete() ? " " : " at best "),
      conshdlrdata->curruserseeed->getMaxWhiteScore() );

   conshdlrdata->curruserseeed = NULL;

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPconshdlrDecompTranslateAndAddCompleteUnpresolvedSeeeds(
   SCIP*                 scip,                 /**< SCIP data structure */
   SCIP_Bool*            success               /** at least one unpresolved seeed coud be tranlsate in a complete presolved one */
   ){

   DEC_DECOMP* newdecomp;

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
         assert( finseeed->checkConsistency(seeedpoolunpresolved) );
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
         SCIP_CALL( conshdlrdata->seeedpool->createDecompFromSeeed( *seeediter, &newdecomp) );

         SCIP_CALL(SCIPstoreSeeedAndDecomp(scip, *seeediter, newdecomp) );

         *success = TRUE;
         SCIPdebugMessagePrint(scip, " SUCCESS: unpresolved complete seeed did translate to complete presolved one \n");
      }
      else {
         SCIPdebugMessagePrint(scip, " unpresolved complete seeed did not translate to complete presolved one \n");
      }
   }

   return SCIP_OKAY;
}

/** method to adapt score for unpresolved decomps @TODO: change score for some paramter settings*/
SCIP_Real SCIPconshdlrDecompAdaptScore(
   SCIP*             scip,
   SCIP_Real         oldscore
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

SCIP_Bool SCIPconshdlrDecompHasDecomp(
   SCIP*    scip
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

   return (conshdlrdata->ndecomps > 0 || conshdlrdata->nincompleteseeeds > 0 ||
      ( conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() > 0 ) ||
      ( conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->incompleteSeeeds.size() > 0 ) ) ;
}

/** returns TRUE iff there is at least one full decompositions */
SCIP_Bool SCIPconshdlrDecompHasCompleteDecomp(
   SCIP*    scip
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

   return (conshdlrdata->ndecomps > 0 ||  (conshdlrdata->seeedpoolunpresolved != NULL && conshdlrdata->seeedpoolunpresolved->finishedSeeeds.size() > 0 ) ) ;
}


SCIP_Bool SCIPconshdlrDecompExistsSelected(
   SCIP* scip
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

    return conshdlrdata->selectedexists;
}

SCIP_RETCODE SCIPconshdlrDecompChooseCandidatesFromSelected(
   SCIP* scip,
   SCIP_Bool updatelist
   ){

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   gcg::Seeedpool* seeedpool;
   gcg::Seeedpool* seeedpoolunpresolved;
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

   SCIPinfoMessage(scip, NULL, "Starting decomposition candidate choosing \n");

   assert(conshdlrdata->candidates != NULL);

  // std::vector<std::pair<SeeedPtr, SCIP_Real> > candidates(0);
   conshdlrdata->candidates->clear();



   seeedpool = conshdlrdata->seeedpool;
   seeedpoolunpresolved = conshdlrdata->seeedpoolunpresolved;

   if( updatelist )
      SCIP_CALL(SCIPconshdlrDecompUpdateSeeedlist(scip) );

   for( size_t selid = 0; selid < conshdlrdata->selected->size(); ++selid )
   {
      selectedseeeds.push_back(conshdlrdata->listall->at(conshdlrdata->selected->at(selid) ) );
   }

   if ( selectedseeeds.size() == 0 )
   {
      SCIPinfoMessage(scip, NULL, "nothing selected, choose everything as selected: \n");
      selectedseeeds = *conshdlrdata->listall;
      SCIPinfoMessage(scip, NULL, "number that is examined: %d \n", selectedseeeds.size() );
   }

   /** if there are selected decomps, check if some of them needs to be finished and do so */
   seeediter = selectedseeeds.begin();
   seeediterend = selectedseeeds.end();

   for( ; seeediter != seeediterend; ++seeediter)
   {
      if( !(*seeediter)->isComplete() && (*seeediter)->isfromunpresolved )
      {
         tofinishunpresolved.push_back(*seeediter);
      }

      if( !(*seeediter)->isComplete() && !(*seeediter)->isfromunpresolved )
      {
         tofinishpresolved.push_back(*seeediter);
      }
   }

   finished = conshdlrdata->seeedpool->finishIncompleteSeeeds(tofinishpresolved);
   finishedunpresolved = conshdlrdata->seeedpoolunpresolved->finishIncompleteSeeeds(tofinishunpresolved);


   seeediter = selectedseeeds.begin();
   seeediterend = selectedseeeds.end();

   /** get decomp candidates and calculate corresponding score (possibly weighted for unpresolved) */
   for( ; seeediter != seeediterend; ++seeediter )
   {
      SeeedPtr seeed = *seeediter ;
      if( seeed->isComplete() && !seeed->isfromunpresolved )
      {
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, seeed->getMaxWhiteScore() ) );
      }
      if( seeed->isComplete() && seeed->isfromunpresolved )
      {
         conshdlrdata->candidates->push_back( std::pair<SeeedPtr, SCIP_Real>(seeed, SCIPconshdlrDecompAdaptScore(scip, seeed->getMaxWhiteScore() ) ) );
      }
   }



   seeediter = finished.begin();
   seeediterend = finished.end();

   for( ; seeediter != seeediterend; ++seeediter)
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, (*seeediter)->getMaxWhiteScore() )  );
   }

   seeediter = finishedunpresolved.begin();
   seeediterend = finishedunpresolved.end();

   for( ; seeediter != seeediterend; ++seeediter)
   {
      conshdlrdata->candidates->push_back(std::pair<SeeedPtr, SCIP_Real>(*seeediter, SCIPconshdlrDecompAdaptScore(scip, (*seeediter)->getMaxWhiteScore() ) ) );
   }

   /* sort decomp candidates according score */
   std::sort( conshdlrdata->candidates->begin(), conshdlrdata->candidates->end(), sort_pred() );




   return SCIP_OKAY;
}



/** 1) all finished seeeds are part of the hash maps and */
/** 2) hash maps are syncronized  */
/** 3) all incomplete seeeds in cons_decomp are not preds of finished decomps (incomplete seeeds should be deleted in that case) */
/** 4) all incomplete seeeds in seeedpoolunpresolved are not preds of finished seeeds in complete seeeds in seeedpoolunpresolved and any seeeds in cons_decomp or in seeedpool*/
/** 5) there are no finished seeeds in seeedpool */
/** 6) selected list is syncron with selected information in seeeds */
/** 7) selected exists is syncronized with seleced list */


SCIP_Bool SCIPconshdlrDecompCheckConsistency(
   SCIP* scip
   ){

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   int ncompleteseeeds;
   int nincompleteseeeds;
   int ncompleteseeedsunpresolved;
 //  int nincompleteseeedsunpresolved;

   int i;
   int selectedcounter;

   std::vector<int> livingnoncompleteseeedids(0); /** this is a vector of seeed ids that should be living (living: there is no complete seeed having ) */
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

   ncompleteseeeds = conshdlrdata->ndecomps;
   nincompleteseeeds = conshdlrdata->nincompleteseeeds;
   ncompleteseeedsunpresolved = (conshdlrdata->seeedpoolunpresolved == NULL ? 0 :  conshdlrdata->seeedpoolunpresolved->getNFinishedSeeeds() );
//   nincompleteseeedsunpresolved = (conshdlrdata->seeedpoolunpresolved == NULL ? 0 :  conshdlrdata->seeedpoolunpresolved->currSeeeds.size() ); /** @todo: check if is this wanted*/



   /** 1) all finished seeeds are part of the hash maps and */
   /** 2) hash maps are syncronized  */

   for( i = 0; i < conshdlrdata->ndecomps ; ++i)
   {
      DEC_DECOMP* dec;
      SeeedPtr seeed;

      if( !SCIPhashmapExists( conshdlrdata->decdecomptoseeed, conshdlrdata->decdecomps[i] ) )
      {
         SCIPwarningMessage(scip, "Warning: Inconsistency in cons_decomp data structure: %d-th decomp is missing in decompToSeeed \n", i);
         return FALSE;
      }

      if( !SCIPhashmapExists( conshdlrdata->seeedtodecdecomp, (void*) conshdlrdata->allrelevantfinishedseeeds[i] ) )
      {
         SCIPwarningMessage(scip, "Warning: Inconsistency in cons_decomp data structure: %d-th seeed is missing in seeedtodecomp \n", i);
         return FALSE;
      }

      seeed = (SeeedPtr) SCIPhashmapGetImage(conshdlrdata->decdecomptoseeed, (void*) conshdlrdata->decdecomps[i]);
      dec = (DEC_DECOMP*) SCIPhashmapGetImage(conshdlrdata->seeedtodecdecomp, (void*) conshdlrdata->allrelevantfinishedseeeds[i]);

//      if ( !seeed->checkConsistency(conshdlrdata->seeedpool) )
//         return FALSE;
//
//      if ( !DECdecompCheckConsistency(scip, dec) )
//         return FALSE;

      if( dec != conshdlrdata->decdecomps[i] )
      {
         SCIPwarningMessage(scip, "Warning: Inconsistency in cons_decomp data structure: %d-th decomp and complete seeed are out of sync (decs[i] != seeedtodec(seeeds[i])) \n", i);
         return FALSE;
      }

      if( seeed != conshdlrdata->allrelevantfinishedseeeds[i] )
      {
         SCIPwarningMessage(scip, "Warning: Inconsistency in cons_decomp data structure: %d-th decomp and complete seeed are out of sync (seeeds[i] != dectoseeed(decs[i])) \n", i);
         return FALSE;
      }
   }

   /** 3) all incomplete seeeds in cons_decomp and complete seeeds in seeedpoolunpresolved  are not preds of finished decomps (incomplete seeeds should be deleted in that case) */

   for( i = 0; i < nincompleteseeeds; ++i )
   {
      livingnoncompleteseeedids.push_back( conshdlrdata->incompleteseeeds[i]->getID() );
   }

   for( i = 0; i < ncompleteseeedsunpresolved; ++i )
   {
      livingnoncompleteseeedids.push_back( conshdlrdata->seeedpoolunpresolved->getFinishedSeeed(i)->getID() );
   }

   std::sort(livingnoncompleteseeedids.begin(), livingnoncompleteseeedids.end());

   for ( i = 0; i < ncompleteseeeds && livingnoncompleteseeedids.size() > 0; ++i )
   {
      SeeedPtr seeed = conshdlrdata->allrelevantfinishedseeeds[i];

      std::vector<int>::const_iterator iter = seeed->listofancestorids.begin();
      std::vector<int>::const_iterator iterend = seeed->listofancestorids.end();

      for( ; iter != iterend; ++iter)
      {

         std::vector<int>::const_iterator loweriter = std::lower_bound( livingnoncompleteseeedids.begin(), livingnoncompleteseeedids.end(), *iter);
         if ( loweriter!= livingnoncompleteseeedids.end() &&  *loweriter == *iter )
         {
            SCIPwarningMessage(scip, "Warning: Inconsistency in cons_decomp data structure: finished seeed with id %d and has living incomplete pred with id %d  \n", seeed->getID(), *iter );
            return FALSE;
         }
      }
   }

   /** 4) all incomplete seeeds in seeedpoolunpresolved are not preds of finished seeeds in complete seeeds in seeedpoolunpresolved and any seeeds in cons_decomp or in seeedpool*/



   /** 5) there are no finished seeeds in seeedpool */

   if ( conshdlrdata->seeedpool != NULL && conshdlrdata->seeedpool->getNFinishedSeeeds() > 0 )
   {
      SCIPwarningMessage(scip, "Warning: There are finished seeeds in seeedpool that should have been stores in cons_decomp  \n" );
      return FALSE;

   }

   /** 6) selected list is syncronized with selected information in seeeds */

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
      SCIPwarningMessage(scip, "Warning: there are selected seeeds not part of the list  \n" );
      return FALSE;
   }


   /** 7) selected exists is syncronized with seleced list */

   if( conshdlrdata->selectedexists != (conshdlrdata->selected->size() > 0) )
   {
      SCIPwarningMessage(scip, "Warning: selectedexists is %d but number of selected is %d   \n", conshdlrdata->selectedexists, conshdlrdata->selected->size() );
      return FALSE;
   }



   return TRUE;
}

/** returns the next seeed id managed by cons_decomp */
   int SCIPconshdlrDecompGetNextSeeedID(
     SCIP* scip
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

      return ++conshdlrdata->seeedcounter;
   }


SCIP_RETCODE DECconshdlrDecompSortDecompositionsByScore(
   SCIP*       scip
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   SCIP_Real* scores;
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);

   if( conshdlr == NULL )
   {
      SCIPerrorMessage("Decomp constraint handler is not included, cannot add detector!\n");
      return SCIP_ERROR;
   }

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   SCIP_CALL_ABORT(SCIPallocBufferArray(scip, &scores, conshdlrdata->ndecomps) );

   for (int i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      assert(DECdecompCheckConsistency(scip, conshdlrdata->decdecomps[i] ));
      scores[i] = DECgetMaxWhiteScore(scip, conshdlrdata->decdecomps[i]);
   }

   SCIPsortRealPtr(scores, (void**)conshdlrdata->decdecomps, conshdlrdata->ndecomps);
   SCIPsortRealPtr(scores, (void**)conshdlrdata->allrelevantfinishedseeeds, conshdlrdata->ndecomps);

   SCIPfreeBufferArray(scip, &scores);

   return SCIP_OKAY;
}

/** interface method to detect the structure including presolving */
SCIP_RETCODE DECdetectStructure(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_RESULT*          result              /**< Result pointer to indicate whether some structure was found */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   conshdlrdata->seeedpool = NULL;

   /** get data of the seeedpool with original vars and conss */
   if ( conshdlrdata->seeedpoolunpresolved == NULL )
      conshdlrdata->seeedpoolunpresolved = new gcg::Seeedpool(scip, CONSHDLR_NAME, FALSE);         /**< seeedpool with original variables and constraints */

   std::vector<int> candidatesNBlocks;                            /**< candidates for number of blocks */
   std::vector<gcg::ConsClassifier*> consClassDistributions;         /**< collection of different constraint class distributions */
   std::vector<gcg::VarClassifier*> varClassDistributions;           /**< collection of different variable class distributions */
   std::vector<SCIP_CONS*> indexToCons;                           /**< stores the corresponding scip constraints pointer */

   std::vector<gcg::SeeedPtr> seeedsunpresolved(0);                    /**< seeeds that were found for the unpresolved problem */

   SCIP_Real* scores;
   int i;

   SCIP_Bool presolveOrigProblem;
   SCIP_Bool calculateOrigDecomps;

   assert(scip != NULL);

   presolveOrigProblem = TRUE;

   SCIPgetBoolParam(scip, "detection/origprob/enabled", &calculateOrigDecomps);



   if( SCIPgetStage(scip) < SCIP_STAGE_TRANSFORMED )
      SCIP_CALL( SCIPtransformProb(scip) );

//<<<<<<< HEAD
//   if( conshdlrdata->ndecomps == 0)
//   {
//      candidatesNBlocks = seeedpoolunpresolved.getSortedCandidatesNBlocks();
//      seeedpoolunpresolved.calcClassifierAndNBlockCandidates(scip);
//   }
//=======

   /** get block number candidates and conslcassifier for original problem*/
   conshdlrdata->seeedpoolunpresolved->calcClassifierAndNBlockCandidates(scip);

   candidatesNBlocks = conshdlrdata->seeedpoolunpresolved->getSortedCandidatesNBlocks();

   /** detection for original problem */
   if( conshdlrdata->ndecomps == 0 && calculateOrigDecomps )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL , NULL, "start finding decompositions for original problem!\n");
      seeedsunpresolved = conshdlrdata->seeedpoolunpresolved->findSeeeds();
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL , NULL, "finished finding decompositions for original problem!\n");

   }

   /** get the cons and var classifier for translating them later*/
   for( i = 0; i < conshdlrdata->seeedpoolunpresolved->getNConsClassifiers(); ++i )
   {
      gcg::ConsClassifier* classifier = new gcg::ConsClassifier( conshdlrdata->seeedpoolunpresolved->getConsClassifier(i) );
      consClassDistributions.push_back( classifier );
   }
   for( i = 0; i < conshdlrdata->seeedpoolunpresolved->getNVarClassifiers(); ++i )
   {
      gcg::VarClassifier* classifier = new gcg::VarClassifier( conshdlrdata->seeedpoolunpresolved->getVarClassifier(i) );
      varClassDistributions.push_back( classifier );
   }

   //Presolving
   if(presolveOrigProblem)
      SCIP_CALL( SCIPpresolve(scip) );


   /** detection for presolved problem */

   if( SCIPgetStage(scip) == SCIP_STAGE_INIT || SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_DIALOG, NULL, "No problem exists, cannot detect structure!\n");

      /** presolving removed all constraints or variables */
      if( SCIPgetNVars(scip) == 0 || SCIPgetNConss(scip) == 0 )
         conshdlrdata->hasrun = TRUE;

      *result = SCIP_DIDNOTRUN;
      return SCIP_OKAY;
   }

   /** start detection clocks */
   SCIP_CALL( SCIPresetClock(scip, conshdlrdata->detectorclock) );
   SCIP_CALL( SCIPstartClock(scip, conshdlrdata->detectorclock) );

   /** @TODO: deprecated? */
   if( conshdlrdata->ndecomps > 0 )
   {
      for( int i = 0; i < conshdlrdata->ndecomps; ++i )
         SCIP_CALL( DECdecompTransform(scip, conshdlrdata->decdecomps[i]) );
   }

   if( conshdlrdata->seeedpool == NULL )
   {
      SCIPdebugMessagePrint(scip, "create seeedpool for current problem, n detectors: %d \n", conshdlrdata->ndetectors);

      conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE);
      SCIPdebugMessagePrint(scip, "created seeedpool for current problem, n detectors: %d \n", conshdlrdata->ndetectors);
   }
   else
      SCIPdebugMessagePrint(scip, "seeedpool is not NULL \n");

     conshdlrdata->seeedpool->calcClassifierAndNBlockCandidates(scip);

   /** get block number candidates and translate orig classification and found seeeds (if any) to presolved problem */
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL , NULL, "started translate seeed method!\n");
      std::vector<gcg::Seeed*> translatedSeeeds(0);
      std::vector<gcg::ConsClassifier*> translatedConsDistributions(0);
      std::vector<gcg::VarClassifier*> translatedVarDistributions(0);

      conshdlrdata->seeedpool->translateSeeedData( conshdlrdata->seeedpoolunpresolved, seeedsunpresolved, translatedSeeeds,
         consClassDistributions, translatedConsDistributions, varClassDistributions, translatedVarDistributions );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL , NULL, "number of translated original seeeds: %d \n " , translatedSeeeds.size() );

      conshdlrdata->seeedpool->populate(translatedSeeeds);

      for ( size_t d = 0; d < translatedConsDistributions.size(); ++d )
         conshdlrdata->seeedpool->addConsClassifier( translatedConsDistributions[d] );

      for ( size_t d = 0; d < translatedVarDistributions.size(); ++d )
         conshdlrdata->seeedpool->addVarClassifier( translatedVarDistributions[d] );

      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL , NULL, "finished translate seeed method!\n");

      for( size_t c = 0; c < candidatesNBlocks.size(); ++c )
         conshdlrdata->seeedpool->addCandidatesNBlocks(candidatesNBlocks[c]);
   }

   conshdlrdata->seeedpool->findDecompositions();

   /** these are the first finished decompositions to add */

   for( i = 0; i < conshdlrdata->seeedpool->getNDecompositions(); ++i )
   {
	     SCIP_CALL( SCIPstoreSeeedAndDecomp(scip, conshdlrdata->seeedpool->getFinishedSeeed(i), conshdlrdata->seeedpool->getDecompositions()[i] ) );
   }

	  conshdlrdata->seeedpool->clearFinishedSeeeds();

   SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
   SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);
   //	  seeedpool.freeCurrSeeeds();

   /* evaluate all decompositions */
   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      DEC_SCORES score;
      score.totalscore = 0.0;

      SCIP_CALL( DECevaluateDecomposition(scip, conshdlrdata->decdecomps[i], &score) );
   }


   SCIP_CALL( SCIPstopClock(scip, conshdlrdata->detectorclock) );

   if( conshdlrdata->ndecomps > 0 )
   {
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "Chosen decomposition with %d blocks of type %s.\n",
         DECdecompGetNBlocks(conshdlrdata->decdecomps[0]), DECgetStrType(DECdecompGetType(conshdlrdata->decdecomps[0])));
  //    GCGsetStructDecdecomp(scip, conshdlrdata->decdecomps[0]);
      *result = SCIP_SUCCESS;
   }
   else
   {
      assert(conshdlrdata->ndecomps == 0);
      SCIPverbMessage(scip, SCIP_VERBLEVEL_NORMAL, NULL, "No decomposition found -- solving with one single block.\n");
      SCIP_CALL( createOneBlockDecomp(scip) );
      *result = SCIP_SUCCESS;
   }
   SCIPdebugMessage("Detection took %fs\n", SCIPclockGetTime(conshdlrdata->detectorclock));

   /* show that we done our duty */
   conshdlrdata->hasrun = TRUE;

   SCIPconshdlrDecompChooseCandidatesFromSelected(scip, TRUE);

//   SCIPhashmapFree( &consToIndex );

   return SCIP_OKAY;
}

///** interface method to detect the structure */
//SCIP_RETCODE DECgetSeeedpoolData(
//   SCIP*                 scip,                     /**< SCIP data structure */
//   int**                 candidatesNBlocks,        /**< pointer to store candidates for number of blocks calculated by the seeedpool */
//   int*                  nCandidates,              /**< pointer to store number of candidates for number of blocks calculated by the seeedpool */
//   int***                consClasses,             /**< pointer to store the  collection of different constraint class distributions */
//   int*                  nConsClassDistributions, /**< pointer to store number of constraint class distributions */
//   int**                 nClassesOfDistribution,   /**< pointer to store numbers of classes of the distributions */
//   SCIP_HASHMAP*         consToIndex,              /**< hashmap from constraints to indices, to be filled */
//   int*                  nConss                    /**< pointer to store number of constraints */
//   )
//{
//   std::cout << "SEEEDPOOL 1" << std::endl;
//   gcg::Seeedpool seeedpool(scip, CONSHDLR_NAME, TRUE);
//
//   std::vector<int> candidatesNBlocksVector = seeedpool.getCandidatesNBlocks();
//   *nCandidates = (int)candidatesNBlocksVector.size();
//   SCIP_CALL( SCIPallocMemoryArray(scip, &(candidatesNBlocks), *nCandidates) );
//   for( int i = 0; i < *nCandidates; ++i )
//      (*candidatesNBlocks)[i] = candidatesNBlocksVector[i];
//
//   std::vector<std::vector<int>> consClassesVector;
//   for( int i = 0; i < seeedpool.getNConssClassDistributions(); ++i )
//      consClassesVector.push_back(seeedpool.getConssClassDistributionVector(i));
//   *nConsClassDistributions = seeedpool.getNConssClassDistributions();
//   SCIP_CALL( SCIPallocMemoryArray(scip, &(consClasses), *nConsClassDistributions) );
//   SCIP_CALL( SCIPallocMemoryArray(scip, &(nClassesOfDistribution), *nConsClassDistributions) );
//   for( int i = 0; i < *nConsClassDistributions; ++i )
//   {
//      (*nClassesOfDistribution)[i] = seeedpool.getNClassesOfDistribution(i);
//      SCIP_CALL( SCIPallocMemoryArray(scip, &(consClasses[i]), (*nClassesOfDistribution)[i]) );
//   }
//   for( int i = 0; i < (int)consClassesVector.size(); ++i )
//   {
//      for( int j = 0; j < (int)consClassesVector[i].size(); ++j )
//      {
//         (*consClasses)[i][j] = consClassesVector[i][j];
//      }
//   }
//
//   *nConss = seeedpool.getNConss();
//   SCIP_CALL_ABORT( SCIPhashmapCreate( &consToIndex, SCIPblkmem(scip), *nConss ) );
//   for( int i = 0; i < seeedpool.getNConss(); ++i )
//      SCIP_CALL_ABORT( SCIPhashmapInsert(consToIndex, seeedpool.getConsForIndex(i), (void*) (size_t) i ) );
//
//
//   return SCIP_OKAY;
//}

/** write
 *  out all detected or provided decompositions */
SCIP_RETCODE DECwriteAllDecomps(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 directory,          /**< directory for decompositions */
   char*                 extension           /**< extension for decompositions */
   )
{
   char name[SCIP_MAXSTRLEN];
   char outname[SCIP_MAXSTRLEN];
   char *pname;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   DEC_DETECTOR *detector;
   DEC_DECOMP *decomp;
   DEC_DECOMP *tmp;
   int i;
//   int j;

   assert(scip != NULL);
   assert(extension != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

   if( conshdlrdata->ndecomps == 0 )
   {
      SCIPwarningMessage(scip, "No decomposition available.\n");
      return SCIP_OKAY;
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s", SCIPgetProbName(scip));
   SCIPsplitFilename(name, NULL, &pname, NULL, NULL);

   tmp = conshdlrdata->decdecomps[0];

   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      decomp = conshdlrdata->decdecomps[i];

      assert(decomp != NULL);
      if( directory != NULL )
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s_%d.%s", directory, pname, i, extension);
      }
      else
      {
         (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d.%s", pname, i, extension);
      }
      conshdlrdata->decdecomps[0] = decomp;
      SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
   }

//   for( i = 0; i < conshdlrdata->ndetectors; ++i )
//   {
//      detector =  conshdlrdata->detectors[i];
//      assert(detector != NULL);
//
//      for( j = 0; j < detector->ndecomps; ++j )
//      {
//         decomp = detector->decomps[j];
//         assert(decomp != NULL);
//         if( directory != NULL )
//         {
//            (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s/%s_%c_%d_%d.%s", directory, pname, detector->decchar, DECdecompGetNBlocks(decomp), j, extension);
//         }
//         else
//         {
//            (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%c_%d_%d.%s", pname, detector->decchar, DECdecompGetNBlocks(decomp), j, extension);
//
//         }
//         conshdlrdata->decdecomps[0] = decomp;
//         SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
//      }
//   }

   /** further, get all read in decompositions */
   for( i = 0; i < conshdlrdata->ndecomps; ++i )
   {
      decomp = conshdlrdata->decdecomps[i];
      detector =  DECdecompGetDetector(decomp);

      if( detector != NULL )
         continue;

      (void) SCIPsnprintf(outname, SCIP_MAXSTRLEN, "%s_%d.%s", pname, DECdecompGetNBlocks(decomp), extension);

      conshdlrdata->decdecomps[0] = decomp;
      SCIP_CALL( SCIPwriteTransProblem(scip, outname, extension, FALSE) );
   }

   conshdlrdata->decdecomps[0] = tmp;

   return SCIP_OKAY;
}

std::vector<SeeedPtr> SCIPconshdlrDecompGetAllRelevantSeeeds(
   SCIP* scip
   )
{

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);

 //  std::vector<SeeedPtr> relevantseeeds(0);

   int maxid  = 0;
   std::vector<SeeedPtr> tmpAllRelevantSeeeds(0);

   for ( size_t i = 0; i < conshdlrdata->seeedpool->allrelevantseeeds.size(); ++i )
   {
      if( conshdlrdata->seeedpool->allrelevantseeeds[i] != NULL && conshdlrdata->seeedpool->allrelevantseeeds[i]->getID() > maxid )
         maxid = conshdlrdata->seeedpool->allrelevantseeeds[i]->getID();
   }

   for ( size_t i = 0; i < conshdlrdata->seeedpoolunpresolved->allrelevantseeeds.size(); ++i )
   {
      if( conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i] != NULL &&  conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i]->getID() > maxid )
         maxid = conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i]->getID();
   }

   tmpAllRelevantSeeeds = std::vector<SeeedPtr>(maxid+1, NULL );

   for ( size_t i = 0; i < conshdlrdata->seeedpoolunpresolved->allrelevantseeeds.size(); ++i )
      {
         if ( conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i] == NULL || conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i]->getID() < 0  )
            continue;
         tmpAllRelevantSeeeds[conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i]->getID()] = conshdlrdata->seeedpoolunpresolved->allrelevantseeeds[i];
      }

   for ( size_t i = 0; i < conshdlrdata->seeedpool->allrelevantseeeds.size(); ++i )
      {
         if ( conshdlrdata->seeedpool->allrelevantseeeds[i] == NULL || conshdlrdata->seeedpool->allrelevantseeeds[i]->getID() < 0  )
            continue;
         tmpAllRelevantSeeeds[conshdlrdata->seeedpool->allrelevantseeeds[i]->getID()] = conshdlrdata->seeedpool->allrelevantseeeds[i];
      }

   return tmpAllRelevantSeeeds;
}




SCIP_RETCODE SCIPconshdlrDecompWriteFamilyTreeLatexFile(
   SCIP*       scip,
   const char* filename,                                 /* filename the output should be written to (including directory) */
   const char* workfolder,                               /* directory in which should be worked */
   std::vector<SeeedPtr> seeeds,                         /* vector of seeed pointers the  family tree should be constructed for */
   SCIP_Bool draft
){

   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert(scip != NULL);

   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   assert(conshdlr != NULL);

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert(conshdlrdata != NULL);


   std::ofstream ofs;
   int curr = -1;
   int currheight = 0;
   SCIP_Real firstsibldist = -1.;

   std::stringstream preambel;
   std::string closing = "\\end{tikzpicture}\n\\end{document}";

   /* collectiotn of treeseeds */
   std::vector<SeeedPtr> treeseeeds(0);
   std::vector<int> treeseeedids(0);
   std::vector<SeeedPtr> allrelevantseeeds = SCIPconshdlrDecompGetAllRelevantSeeeds(scip);

   std::vector<SCIP_Bool> isseeedintree(allrelevantseeeds.size(), FALSE );

   int root = -1;
   std::vector<int> parents(allrelevantseeeds.size(), -1);
   std::vector<std::vector<int> > childs (allrelevantseeeds.size(), std::vector<int>(0));
   std::vector<std::vector<SCIP_Bool> > childsfinished(allrelevantseeeds.size(), std::vector<SCIP_Bool>(0));
   std::vector<SCIP_Bool> visited(allrelevantseeeds.size(), FALSE);


   conshdlrdata->helpvisucounter = 0;

   /** check allrelevant seeeds **/
   for( size_t s = 0; s < allrelevantseeeds.size(); ++s )
   {
      assert(allrelevantseeeds[s] == NULL || (int) s == allrelevantseeeds[s]->getID() );
   }

   /** 1) find relevant seeeds in tree and build tree */
   for( size_t s = 0; s < seeeds.size(); ++s )
   {
      int currid;
      if ( seeeds[s] == NULL )
         continue;
      currid = seeeds[s]->getID();
      if( !isseeedintree[seeeds[s]->getID()] )
      {
         isseeedintree[seeeds[s]->getID()] = TRUE;
         treeseeeds.push_back( seeeds[s]);
         treeseeedids.push_back(seeeds[s]->getID());
      }
      else
         break;

      for( size_t i = 0; i < seeeds[s]->listofancestorids.size(); ++i )
      {
         int ancestorid;
         ancestorid = seeeds[s]->listofancestorids[seeeds[s]->listofancestorids.size() - i -1];
         parents[currid] = ancestorid;
         childs[ancestorid].push_back(currid);
         childsfinished[ancestorid].push_back(FALSE);

         if( !isseeedintree[ancestorid] )
         {
            isseeedintree[ancestorid] = TRUE;
            assert(allrelevantseeeds[ancestorid] != NULL);
            treeseeeds.push_back( allrelevantseeeds[ancestorid] );
            treeseeedids.push_back(ancestorid);
            if( i == seeeds[s]->listofancestorids.size() -1 )
            {
               root = ancestorid;
            }
            currid = ancestorid;
         }
         else
            break;
      }
   }

   for( size_t i = 0; i < treeseeeds.size(); ++i )
   {
      SeeedPtr seeed;
      std::stringstream decompfilename;

      seeed = treeseeeds[i];

      decompfilename << workfolder << "/" << getSeeedFolderLatex(seeed);

      seeed->showScatterPlot(conshdlrdata->seeedpool, TRUE, decompfilename.str().c_str(), draft );
   }

   //  finishedSeeeds[0]->showScatterPlot(this, TRUE, "./testdecomp/001.pdf") ;

   firstsibldist = 1. / (childs[root].size() - 1 );
   preambel.precision(2);

   preambel << "\\documentclass[a4paper,landscape]{scrartcl}\n\\usepackage{fancybox}\n\\usepackage{tikz}";
   preambel << "\n\\usetikzlibrary{positioning}\n\\title{Detection Tree}\n\\date{}\n\\begin{document}\n\n";
   preambel << "\\begin{tikzpicture}[level/.style={sibling distance=" << firstsibldist << "\\textwidth/#1}, level distance=12em, ->, dashed]\n\\node";



   /** start writing file */
   ofs.open (filename, std::ofstream::out );
   ofs << preambel.str();

   /** iterate tree and write file */
   curr = root;
   while ( curr != -1 )
   {
      if( !visited[curr] )
      {
         /** write node */
         ofs << writeSeeedIncludeLatex( allrelevantseeeds[curr], workfolder );
         /* set node visited */
         visited[curr] = TRUE;
         if( parents[curr] != -1 )
            finishnextchild(childs[parents[curr]], childsfinished[parents[curr]], curr);

      }
      if ( unfinishedchildexists(childsfinished[curr] ) )
      {
         int unfinishedchild = getfirstunfinishedchild(childsfinished[curr], childs[curr] );
         /* is first child unfinihsed? */
         //         if( unfinishedchild == childs[curr][0] )
         ofs << " child { node " ;
         curr = unfinishedchild;
         ++currheight;
      }
      else
      {
         if ( parents[curr] != -1 ){
            ofs << writeSeeedDetectorChainInfoLatex( allrelevantseeeds[curr], currheight, conshdlrdata->helpvisucounter);
            ++conshdlrdata->helpvisucounter;
         }
         --currheight;
         curr = parents[curr];
         if( curr != -1)
            ofs << " } " ;
      }
   }

   ofs << ";" << std::endl;
   for( size_t i = 0; i < treeseeeds.size(); ++i)
   {
      ofs << writeSeeedInfoLatex( treeseeeds[i] );
   }

   ofs << closing << std::endl;

   ofs.close();

   return SCIP_OKAY;
}




/** write
 *  out all detected or provided decompositions */
/** write family tree **/
SCIP_RETCODE DECwriteFamilyTree(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< filename the output should be written to (including directory) */
   const char*           workfolder,         /**< directory in which should be worked */
   int                   ndecompositions,    /**< the number of (complete) decompositions in order of a certain measure (atm: max white) */
   SCIP_Bool             draft               /**< draft mode will not visualize non-zeros but is faster and takes less memory */
   )
{

	SCIP_CONSHDLR* conshdlr;
	SCIP_CONSHDLRDATA* conshdlrdata;
	assert(scip != NULL);

	conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
	assert(conshdlr != NULL);

	conshdlrdata = SCIPconshdlrGetData(conshdlr);

	std::vector<SeeedPtr> tovisualize(0);
	assert(conshdlrdata != NULL);

	if ( SCIPconshdlrDecompExistsSelected(scip) )
	{
	   for ( size_t i = 0; tovisualize.size() <= (size_t) ndecompositions &&  i < conshdlrdata->selected->size(); ++i )
	   {
	      if( conshdlrdata->listall->at(conshdlrdata->selected->at(i))->isComplete() )
	         tovisualize.push_back( conshdlrdata->listall->at(conshdlrdata->selected->at(i)  ) );
	   }
	}
	else
	{
	   SCIPconshdlrDecompUpdateSeeedlist(scip);
	   for ( size_t i = 0; tovisualize.size() <= (size_t) ndecompositions &&  i < conshdlrdata->listall->size(); ++i )
	   {
	      if( conshdlrdata->listall->at(i)->isComplete() )
	         tovisualize.push_back( conshdlrdata->listall->at(i) );
	   }
	}

	SCIPconshdlrDecompWriteFamilyTreeLatexFile(scip, filename, workfolder, tovisualize, draft);

	return SCIP_OKAY;
}

/** returns the best known decomposition, if available and NULL otherwise */
DEC_DECOMP* DECgetBestDecomp(
   SCIP*                 scip                /**< SCIP data structure */
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


   //seeedpool = ( SCIPconshdlrDecompIsBestCandidateUnpresolved(scip) ? conshdlrdata->seeedpoolunpresolved :  conshdlrdata->seeedpool );

   if( conshdlrdata->seeedpool == NULL )
      conshdlrdata->seeedpool = new gcg::Seeedpool(scip, CONSHDLR_NAME, TRUE);

   seeedpool = conshdlrdata->seeedpool;
   seeedpoolunpresolved = conshdlrdata->seeedpoolunpresolved ;

   DECconshdlrDecompSortDecompositionsByScore(scip);

   seeed = conshdlrdata->candidates->at( 0 ).first;

   SCIPdialogMessage(scip, NULL, "In get bestdecomp\n");

   if( SCIPconshdlrDecompIsBestCandidateUnpresolved(scip) )
   {
      std::vector<SeeedPtr> seeedtotranslate(0);
      std::vector<SeeedPtr> translatedseeeds(0);

      seeedtotranslate.push_back(seeed);
      seeedpool->translateSeeeds(seeedpoolunpresolved, seeedtotranslate, translatedseeeds);
      seeed = translatedseeeds[0];
   }


   seeedpool->createDecompFromSeeed(seeed, &decomp) ;

   return decomp;


   /** OLD HACK! */
   if( conshdlrdata->selectedexists )
   {

      seeedpool->createDecompFromSeeed(conshdlrdata->listall->at( conshdlrdata->selected->at(0) ), &decomp) ;
      return decomp;
   }

   if( conshdlrdata->ndecomps > 0 )
      return conshdlrdata->decdecomps[0];

   else if ( conshdlrdata->createbasicdecomp)
   {
      SCIP_RETCODE retcode;
      DEC_DECOMP* decomp = NULL;
      retcode = DECcreateBasicDecomp(scip, &decomp);
      assert(retcode == SCIP_OKAY);
      assert(decomp != NULL );

      retcode = SCIPconshdlrDecompAddDecdecomp(scip, decomp);
      if( retcode != SCIP_OKAY )
      {
         SCIPerrorMessage("Could not add decomp to cons_decomp!\n");
         return NULL;
      }

      assert(conshdlrdata->ndecomps > 0);
      assert(conshdlrdata->decdecomps[0] != NULL);
      return conshdlrdata->decdecomps[0];
   }

   SCIPdebugMessagePrint(scip, "no decomps out there \n");

   return NULL;
}

/** writes out a list of all detectors */
void DECprintListOfDetectors(
   SCIP*                 scip                /**< SCIP data structure */
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

/** returns whether the detection has been performed */
SCIP_Bool DEChasDetectionRun(
   SCIP*                 scip                /**< SCIP data structure */
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

/** returns the character of the detector */
char DECdetectorGetChar(
   DEC_DETECTOR*         detector            /**< pointer to detector */
)
{
   if( detector == NULL )
     return '0';
   else
      return detector->decchar;
}

/** display statistics about detectors */
SCIP_RETCODE GCGprintDetectorStatistics(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file or NULL for standard output */
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
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "  %-10.10s       :   %8.2f %10d    ", conshdlrdata->detectors[i]->name, conshdlrdata->detectors[i]->dectime, conshdlrdata->detectors[i]->ndecomps );
      for( j = 0; j < conshdlrdata->detectors[i]->ndecomps; ++j )
      {
         SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " %d", DECdecompGetNBlocks(conshdlrdata->detectors[i]->decomps[j]));
      }
      SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
   }
   return SCIP_OKAY;
}

/** resets the parameters to their default value */
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
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if(SCIPgetNVars(scip) + SCIPgetNConss(scip) < DEFAULT_LEVENSHTEIN_MAXMATRIXHALFPERIMETER)
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", FALSE) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIP_Result result;

      char paramname[SCIP_MAXSTRLEN];
      SCIP_Bool paramval;
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPresetParam(scip, paramname) );

      result = SCIP_DIDNOTRUN;
      if( conshdlrdata->detectors[i]->setParamDefault != NULL )
         conshdlrdata->detectors[i]->setParamDefault(scip, conshdlrdata->detectors[i], &result);
      if( !quiet )
      {
         SCIP_Bool written = FALSE;

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
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

/** sets the parameters to aggressive values */
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
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if(SCIPgetNVars(scip) + SCIPgetNConss(scip) < AGGRESSIVE_LEVENSHTEIN_MAXMATRIXHALFPERIMETER)
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

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
            SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
            if( paramval == TRUE )
            {
               SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
               written = TRUE;
            }

            (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
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

/** disables detectors */
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
      (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);

      SCIP_CALL( SCIPsetBoolParam(scip, paramname, FALSE) );
      if( !quiet )
      {
         SCIPinfoMessage(scip, NULL, "%s = FALSE\n", paramname);
      }
   }

   return SCIP_OKAY;
}

/** sets the parameters to fast values */
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

   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/nnonzeros/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/scipconstype/enabled", TRUE) );
   SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamenonumbers/enabled", TRUE) );

   if( SCIPgetNVars(scip) + SCIPgetNConss(scip) < FAST_LEVENSHTEIN_MAXMATRIXHALFPERIMETER )
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", TRUE) );
   else
      SCIP_CALL(SCIPsetBoolParam(scip, "detection/consclassifier/consnamelevenshtein/enabled", FALSE) );

   for( i = 0; i < conshdlrdata->ndetectors; ++i )
   {
      SCIP_Result result;

      result = SCIP_DIDNOTRUN;
      if( conshdlrdata->detectors[i]->setParamFast != NULL )
         conshdlrdata->detectors[i]->setParamFast(scip, conshdlrdata->detectors[i], &result);
      if( !quiet )
      {
         char paramname[SCIP_MAXSTRLEN];
         SCIP_Bool paramval;
         SCIP_Bool written = FALSE;

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/enabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/origenabled", conshdlrdata->detectors[i]->name);
         SCIP_CALL( SCIPgetBoolParam(scip, paramname, &paramval) );
         if( paramval == TRUE )
         {
            SCIPinfoMessage(scip, NULL, "%s = %s\n", paramname, paramval == TRUE ? "TRUE" : "FALSE");
            written = TRUE;
         }

         (void) SCIPsnprintf(paramname, SCIP_MAXSTRLEN, "detectors/%s/finishingenabled", conshdlrdata->detectors[i]->name);
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

/** sets detector parameters values to
 *
 *  - SCIP_PARAMSETTING_DEFAULT which are the default values of all detector parameters
 *  - SCIP_PARAMSETTING_FAST such that the time spend for detection is decreased
 *  - SCIP_PARAMSETTING_AGGRESSIVE such that the detectors produce more decompositions
 *  - SCIP_PARAMSETTING_OFF which turns off all detection
 */
SCIP_RETCODE GCGsetDetection(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PARAMSETTING     paramsetting,       /**< parameter settings */
   SCIP_Bool             quiet               /**< should the parameter be set quiet (no output) */
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
