/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2015 Operations Research, RWTH Aachen University       */
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

/**@file   class_seeedpool.h
 * @brief  class with functions for seeed pool where a seeed is a (potentially incomplete) description of a decomposition (not to confuse with the band from German capital)
 * @author Michael Bastubbe
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEEDPOOL_H__
#define GCG_CLASS_SEEEDPOOL_H__

#include "objscip/objscip.h"
#include <vector>
#include <tr1/unordered_map> //c++ hashmap
#include <unordered_map>
#include <functional>
#include <string>
#include <utility>
#include "gcg.h"

#include "class_seeed.h"
#include "class_consclassifier.h"
#include "class_varclassifier.h"



struct Seeed_Propagation_Data
{
   gcg::Seeedpool* seeedpool;
   gcg::Seeed* seeedToPropagate;
   gcg::Seeed** newSeeeds;
   int nNewSeeeds;
};

namespace gcg {


//typedef boost::shared_ptr<Seeed> SeeedPtr;
typedef Seeed* SeeedPtr;


// Only for pairs of std::hash-able types for simplicity.
// You can of course template this struct to allow other hash functions
struct pair_hash {
    template <class T1, class T2>
    std::size_t operator () (const std::pair<T1,T2> &p) const {
        auto h1 = std::hash<T1>{}(p.first);
        auto h2 = std::hash<T2>{}(p.second);

        // overly simple hash combination
        return h1 ^ h2;
    }
};


class Seeedpool
{   /*lint -esym(1712,Seeedpool)*/

private:
   SCIP*                 						      scip;              	   /**< SCIP data structure */



   int                                          maxndetectionrounds;    /**< maximum number of detection rounds */
   int											         nTotalSeeeds;        	/**< number of created seeeeds, used to give next id */
   std::vector<std::vector<int>> 				   varsForConss; 	   	   /**< stores for every constraint the indices of variables that are contained in the constraint */
   std::vector<std::vector<double>>             valsForConss;           /**< stores for every constraint the coefficients of variables that are contained in the constraint (i.e. have a nonzero coefficient) */
   std::vector<std::vector<int>>    				conssForVars;      		/**< stores for every variable the indices of constraints containing this variable */
   std::vector<SCIP_CONS*> 			   			consToScipCons;	      /**< stores the corresponding scip constraints pointer */
   std::vector<SCIP_VAR*> 				       		varToScipVar;		      /**< stores the corresponding scip variable pointer */
   std::vector<DEC_DETECTOR*> 				   	detectorToScipDetector; /**< stores the corresponding SCIP detector pinter */
   std::vector<DEC_DETECTOR*>                   detectorToFinishingScipDetector; /**< stores the corresponding finishing SCIP detector pointer*/
   std::tr1::unordered_map<SCIP_CONS*, int> 	   scipConsToIndex;	      /**< maps SCIP_CONS* to the corresponding index */
   std::tr1::unordered_map<SCIP_VAR*, int>  	   scipVarToIndex;		   /**< maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map<DEC_DETECTOR*, int>  scipDetectorToIndex;		/**< maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map<DEC_DETECTOR*, int>  scipFinishingDetectorToIndex;    /**< maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map< std::pair<int, int>, SCIP_Real, pair_hash>  valsMap;               /**< maps an entry of the matrix to its value, zeros are omitted */

   int 										         	nVars;                  /**< number of variables */
   int 										         	nConss;                 /**< number of constraints */
   int										         	nDetectors;             /**< number of detectors */

   int                                          nnonzeros;              /**< number of nonzeros */

   int                                          nFinishingDetectors;             /**< number of detectors */


   DEC_DECOMP**                                 decompositions;         /**< decompositions found by the detectors */
   int                                          ndecompositions;        /**< number of decompositions found by the detectors */

   /** oracle data */
   std::vector<std::pair<int,int> >                             candidatesNBlocks;      	/**< candidate for the number of blocks  */


   std::vector<ConsClassifier*>                 consclassescollection; /**< collection of different constraint class distributions  */
   std::vector<VarClassifier*>                  varclassescollection;  /**< collection of different variabale class distributions   */


   SCIP_Bool                                    transformed;            /**< corresponds the matrix datastructure to the transformed problem */

   std::vector<SeeedPtr>                        seeedstopopulate;      /**< seeeds that are translated seeeds from found ones for the original problem */



public:

   std::vector<SeeedPtr>                        incompleteSeeeds;       /**< vector of incomplete seeeds that can be used for initialization */
   std::vector<SeeedPtr>                        currSeeeds;             /**< vector of current (open) seeeds */
   std::vector<SeeedPtr> 						      finishedSeeeds;		   /**< vector of finished seeeds */
   std::vector<SeeedPtr>                        allrelevantseeeds;      /** collection of all relevant seeeds, allrelevaseeeds[i] contains seeed with id i; non relevant seeeds are repepresented by a null pointer */

   /** constructor */
   Seeedpool(
      SCIP*             scip, /**< SCIP data structure */
	  const char*	conshdlrName,
	  SCIP_Bool    transformed
      );

   ~Seeedpool();

   SCIP_RETCODE calcConsClassifierAndNBlockCandidates(
      SCIP*               givenScip        /**< SCIP data structure */
   );

   /** finds seeeds  */
   /*
    * @return user has to free Seeeds
    */
   std::vector<SeeedPtr> findSeeeds(
      );



   /** method to complete a set of incomplete seeeds with the help of all included detectors that implement a finishing method */
   /*
    * @return set of completed decomposition
    * */

   std::vector<SeeedPtr>  finishIncompleteSeeeds(
      std::vector<SeeedPtr> incompleteseeeds
    );

   /** finds decompositions  */
   void findDecompositions(
    );

   /** translates seeeds and classifiers if the index structure of the problem has changed, e.g. due to presolving */
   void translateSeeedData(
      Seeedpool* otherpool,                              /**< old seeedpool */
      std::vector<Seeed*> otherseeeds,                   /**< seeeds to be translated */
      std::vector<Seeed*>& newseeeds,                    /**< translated seeeds (pass empty vector) */
      std::vector<ConsClassifier*> otherconsclassifiers, /**< consclassifiers to be translated */
      std::vector<ConsClassifier*>& newconsclassifiers,  /**< translated consclassifiers (pass empty vector) */
      std::vector<VarClassifier*> othervarclassifiers,   /**< varclassifiers to be translated */
      std::vector<VarClassifier*>& newvarclassifiers     /**< translated varclassifiers (pass empty vector) */
   );

   /** translates seeeds if the index structure of the problem has changed, e.g. due to presolving */
   void translateSeeeds(
      Seeedpool* otherpool,                              /**< old seeedpool */
      std::vector<Seeed*> otherseeeds,                   /**< seeeds to be translated */
      std::vector<Seeed*>& newseeeds                     /**< translated seeeds (pass empty vector) */
   );

   void populate(std::vector<SeeedPtr> seeeds);

   SCIP_RETCODE prepareSeeed( SeeedPtr seeed);

   void freeCurrSeeeds();


   void addSeeedToIncomplete(SeeedPtr seeed);

   void addSeeedToCurr(SeeedPtr seeed);

   void addSeeedToFinished(SeeedPtr seeed);

   void sortAllRelevantSeeeds();

   /** access the variable indices of matrix constraint-wise */
   const  int *  getVarsForCons(int consIndex);

   /** access the coefficients constraint-wise */
    const  SCIP_Real *  getValsForCons(int consIndex);


   /** access coefficient matrix variable-wise */
   const  int * getConssForVar(int varIndex);

   /** returns the number of variables for the given constraint */
   int getNVarsForCons(int consIndex);

   /** returns the number of constraints for the given variable */
   int getNConssForVar(int varIndex);

   SCIP_VAR* getVarForIndex(int varIndex);

   SCIP_CONS* getConsForIndex(int consIndex);

   DEC_DETECTOR* getDetectorForIndex(int detectorIndex);

   DEC_DETECTOR* getFinishingDetectorForIndex(int detectorIndex);

   SCIP_Real getVal(int row, int col);

   int getIndexForVar(SCIP_VAR* var);

   int getIndexForCons(SCIP_CONS* cons);

   int getIndexForDetector(DEC_DETECTOR* detector);

   int getIndexForFinishingDetector(DEC_DETECTOR* detector);

   int getNewIdForSeeed();

   void decrementSeeedcount();

   DEC_DECOMP** getDecompositions();

   int getNDecompositions();

   int getNNonzeros();

   int getNDetectors();

   int getNFinishingDetectors();

   int getNVars();

   int getNConss();

   std::vector<int> getSortedCandidatesNBlocks();

   void addCandidatesNBlocks(
      int                 candidate            /**< candidate for block size */
      );

   void calcCandidatesNBlocks();

   int getNConssClassDistributions();

   int* getConssClassDistribution(int consclassdistr);

   std::vector<int> getConssClassDistributionVector(int consclassdistr);

   int getNClassesOfDistribution(int consclassdistr);

   /** returns number of different constraint classifiers */
   int getNConsClassifiers();

   /** returns pointer to a constraint classifier */
   ConsClassifier* getConsClassifier(
      int classifierIndex                     /**< index of constraint classifier */
   );

   ConsClassifier* createConsClassifierForSCIPConstypes(
      );

   ConsClassifier* createConsClassifierForConsnamesDigitFreeIdentical(
      );

   ConsClassifier* createConsClassifierForConsnamesLevenshteinDistanceConnectivity(
      int connectivity
         );

   ConsClassifier* createConsClassifierForNNonzeros(
      );

   /** adds a constraint classifier if it is no duplicate of an existing constraint classifier */
   void addConsClassifier(
      ConsClassifier*               classifier              /**< consclassifier to be added */
   );

   /** adds constraint classifiers with a reduced number of classes */
   void reduceConsclasses();

   /** returns number of different variable classifiers */
   int getNVarClassifiers();

   /** returns pointer to a variable classifier */
   VarClassifier* getVarClassifier(
      int classifierIndex                     /**< index of variable classifier */
   );

   VarClassifier* createVarClassifierForSCIPVartypes(
   );

   /** adds a variable classifier if it is no duplicate of an existing variable classifier */
   void addVarClassifier(
      VarClassifier*               classifier              /**< varclassifier to be added */
   );

   /** adds variable classifiers with a reduced number of classes */
   void reduceVarclasses();

   std::vector<SeeedPtr> removeSomeOneblockDecomps(
      std::vector<SeeedPtr> givenseeeds);



   /**
    * creates a decomposition for a given seeed
    */
   SCIP_RETCODE createDecompFromSeeed(
      SeeedPtr       seeed,                                 /** seeed the decomposition is created for */
      DEC_DECOMP**   newdecomp                              /** the new decomp created from the seeed */
      );

   /**
    * creates a seeed for a given decomposition, atm dummy method
    */
   SCIP_RETCODE createSeeedFromDecomp(
      DEC_DECOMP* decomp,                                    /** decomposition the seeed is created for */
      SeeedPtr*   newseeed                                   /** the new seeed created from the decomp */
      );


   /**
    * returns transformation information
    */
   SCIP_Bool getTransformedInfo(
      );

private:

   /** calculates necessary data for translating seeeds and classifiers */
   void calcTranslationMapping(
      Seeedpool* origpool,
      std::vector<int>& rowothertothis,
      std::vector<int>& rowthistoother,
      std::vector<int>& colothertothis,
      std::vector<int>& colthistoother,
      std::vector<int>& missingrowinthis
   );

   /** returns translated Seeeds derived from given mapping data */
   std::vector<Seeed*> getTranslatedSeeeds(
      std::vector<Seeed*>& otherseeeds,                     /**< seeeds to be translated */
      std::vector<int>& rowothertothis,
      std::vector<int>& rowthistoother,
      std::vector<int>& colothertothis,
      std::vector<int>& colthistoother
   );

   /** returns translated ConsClassifiers derived from given mapping data */
   std::vector<ConsClassifier*> getTranslatedConsClassifiers(
      std::vector<ConsClassifier*>& otherclassifiers,       /**< consclassifiers to be translated */
      std::vector<int>& rowothertothis,
      std::vector<int>& rowthistoother
   );

   /** returns translated VarClassifiers derived from given mapping data */
   std::vector<VarClassifier*> getTranslatedVarClassifiers(
      std::vector<VarClassifier*>& otherclassifiers,        /**< varclassifiers to be translated */
      std::vector<int>& colothertothis,
      std::vector<int>& colthistoother
   );

};

} /* namespace gcg */
#endif /* GCG_CLASS_SEEEDPOOL_H__ */
