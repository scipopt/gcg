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
   std::vector<SeeedPtr> 						      currSeeeds;				   /**< vector of current (open) seeeds */
   std::vector<SeeedPtr> 						      finishedSeeeds;		   /**< vector of current (open) seeeds */

   int                                          maxndetectionrounds;    /**< maximum number of detection rounds */
   int											         nTotalSeeeds;        	/**< number of created seeeeds, used to give next id */
   std::vector<std::vector<int>> 				   varsForConss; 	   	   /**< stores for every constraint the indices of variables that are contained in the constraint */
   std::vector<std::vector<double>>             valsForConss;           /**< stores for every constraint the coefficients of variables that are contained in the constraint (i.e. have a nonzero coefficient) */
   std::vector<std::vector<int>>    				conssForVars;      		/**< stores for every variable the indices of constraints containing this variable */
   std::vector<SCIP_CONS*> 			   			consToScipCons;	      /**< stores the corresponding scip constraints pointer */
   std::vector<SCIP_VAR*> 				       		varToScipVar;		      /**< stores the corresponding scip variable pointer */
   std::vector<DEC_DETECTOR*> 				   	detectorToScipDetector; /**< stores the corresponding SCIP detector pinter */
   std::vector<DEC_DETECTOR*>                   detectorToFinishingScipDetector; /**< stores the corresponding finishing SCIP detector pinter */
   std::tr1::unordered_map<SCIP_CONS*, int> 	   scipConsToIndex;	      /**< maps SCIP_CONS* to the corresponding index */
   std::tr1::unordered_map<SCIP_VAR*, int>  	   scipVarToIndex;		   /**< maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map<DEC_DETECTOR*, int>  scipDetectorToIndex;		/**< maps SCIP_VAR* to the corresponding index */
   std::tr1::unordered_map<DEC_DETECTOR*, int>  scipFinishingDetectorToIndex;    /**< maps SCIP_VAR* to the corresponding index */


   std::tr1::unordered_map< std::pair<int, int>, SCIP_Real, pair_hash>  valsMap;               /**< maps an entry of the matrix to its value, zeros are omitted */

   int 										         	nVars;                  /**< number of variables */
   int 										         	nConss;                 /**< number of constraints */
   int										         	nDetectors;             /**< number of detectors */

   int                                          nFinishingDetectors;             /**< number of detectors */


   DEC_DECOMP**                                 decompositions;         /**< decompositions found by the detectors */
   int                                          ndecompositions;        /**< number of decompositions found by the detectors */

   /** oracle data */
   std::vector<int>                             candidatesNBlocks;      /**< candidate for the number of blocks  */

   std::vector<std::vector<int> >               consclassescollection;  /**< collection of different constraint class distributions  */
   std::vector<int >                            consclassesnclasses;    /**< number of classes of the corresponding distribution */

   SCIP_Bool                                    transformed;            /**< corresponds the matrix datastructure to the transformed problem */

   std::vector<SeeedPtr>                        translatedOrigSeeeds;   /**< seeeds that are translated seeeds from found ones for the original problem */

public:

   /** constructor */
   Seeedpool(
      SCIP*             scip, /**< SCIP data structure */
	  const char*	conshdlrName,
	  SCIP_Bool    transformed
      );

   ~Seeedpool();

   /** finds seeedss  */
   /*
    * @return user has to free
    */
   std::vector<SeeedPtr> findSeeeds(
      );


   /** finds decompositions  */
   void findDecompositions(
    );

     std::vector<SeeedPtr> translateSeeeds( Seeedpool* otherpool, std::vector<Seeed*> otherseeeds );

   void populate(std::vector<SeeedPtr> seeeds);

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

   int getNDetectors();

   int getNFinishingDetectors();

   int getNVars();

   int getNConss();

   std::vector<int> getCandidatesNBlocks() const;

   void addCandidatesNBlocks(
      int                 candidate            /**< candidate for block size */
      );

   void calcCandidatesNBlocks();

   int getNConssClassDistributions();

   int* getConssClassDistribution(int consclassdistr);

   std::vector<int> getConssClassDistributionVector(int consclassdistr);

   int getNClassesOfDistribution(int consclassdistr);

   void addConssClassesForSCIPConstypes(
      );

   void addConssClassesForConsnamesDigitFreeIdentical(
      );

   void addConssClassesForConsnamesLevenshteinDistanceConnectivity(
      int connectivity
         );

   void addConssClassesForNNonzeros(
      );

   void addConssClassDistribution(
      std::vector<int>              conssClassDistribution,
      std::vector<SCIP_CONS*>       indexToCons
      );

   bool distributionIsNoDuplicateOfDistributions(
      std::vector<int>              compDistribution,
      int                           nClasses,
      std::vector<std::vector<int>> distributions
      );


};

} /* namespace gcg */
#endif /* GCG_CLASS_SEEEDPOOL_H__ */
