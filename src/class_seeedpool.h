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

/**@file   class_seeedpool.h
 * @brief  class with functions for seeed pool where a seeed is a (potentially incomplete) description of a decomposition (not to confuse with the band from German capital)
 * @author Michael Bastubbe
 * @author Julius Hense
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEEDPOOL_H__
#define GCG_CLASS_SEEEDPOOL_H__

#include "objscip/objscip.h"
#include <vector>

#if __cplusplus >= 201103L
#include <unordered_map>
using std::unordered_map;
#else
#include <tr1/unordered_map>
using std::tr1::unordered_map;
#endif

#include <functional>
#include <string>
#include <utility>
#include "gcg.h"

#include "class_seeed.h"
#include "class_consclassifier.h"
#include "class_varclassifier.h"


/**
 * \brief interface data structure for the detector calling methods
 */
struct Seeed_Propagation_Data
{
   gcg::Seeedpool* seeedpool;                /**< current seeedpool to consider */
   gcg::Seeed* seeedToPropagate;             /**< seeed (aka partial decomposition) to be propagted in next detector call */
   gcg::Seeed** newSeeeds;                   /**< array of new seeeds to filled by the detetcor methods */
   int nNewSeeeds;                           /**< number of new neeed, set by the detector methods */
};


namespace gcg{

enum GCG_PROBLEM_TRANSFORMED_STATUS
{
   original,
   transformed,
   native

};



/** forward declaration */
typedef Seeed* SeeedPtr;



/**
 *  @brief combine two hash function of objects of a pair to get a vaulue for the pair
 */
struct pair_hash
{
   template<class T1, class T2>
   std::size_t operator()(
      const std::pair<T1, T2> &p) const
   {
      auto h1 = std::hash<T1>{}( p.first );
      auto h2 = std::hash<T2>{}( p.second );

      /* overly simple hash combination */
      return h1 ^ h2;
   }
};


/**
 * class to manage the detection process and data for one coefficient matrix of a MIP, usually there is one seeedpool for the original and one seeedpool for the transformed problem
 */
class Seeedpool
{ /*lint -esym(1712,Seeedpool)*/

private:
   SCIP* scip;                                                 /**< SCIP data structure */
   std::vector<SeeedPtr> incompleteSeeeds;                     /**< vector of incomplete seeeds that can be used for initialization */
   std::vector<SeeedPtr> currSeeeds;                           /**< vector of current (open) seeeds */
   std::vector<SeeedPtr> finishedSeeeds;                       /**< vector of finished seeeds */
   std::vector<SeeedPtr> ancestorseeeds;                       /**< @todo revise comment! collection of all relevant seeeds, allrelevaseeeds[i] contains seeed with id i; non relevant seeeds are repepresented by a null pointer */
   int maxndetectionrounds;                                    /**< maximum number of detection rounds */
   std::vector<std::vector<int>> varsForConss;                 /**< stores for every constraint the indices of variables
                                                                 *< that are contained in the constraint */
   std::vector<std::vector<double>> valsForConss;              /**< stores for every constraint the coefficients of
                                                                 *< variables that are contained in the constraint (i.e.
                                                                 *< have a nonzero coefficient) */
   std::vector<std::vector<int>> conssForVars;                 /**< stores for every variable the indices of constraints
                                                                 *< containing this variable */
   std::vector<SCIP_CONS*> consToScipCons;                     /**< stores the corresponding scip constraints pointer */
   std::vector<SCIP_VAR*> varToScipVar;                        /**< stores the corresponding scip variable pointer */
   std::vector<DEC_DETECTOR*> detectorToScipDetector;          /**< stores the corresponding SCIP detector pinter */
   std::vector<DEC_DETECTOR*> detectorToFinishingScipDetector; /**< stores the corresponding finishing SCIP detector pointer*/
   std::vector<DEC_DETECTOR*> detectorToPostprocessingScipDetector; /**< stores the corresponding postprocessing SCIP detector pointer*/
   std::vector<std::vector<int>> conssadjacencies;
   unordered_map<SCIP_CONS*, int> scipConsToIndex;   /**< maps SCIP_CONS* to the corresponding index */
   unordered_map<SCIP_VAR*, int> scipVarToIndex;     /**< maps SCIP_VAR* to the corresponding index */
   unordered_map<DEC_DETECTOR*, int> scipDetectorToIndex;        /**< maps SCIP_VAR* to the corresponding index */
   unordered_map<DEC_DETECTOR*, int> scipFinishingDetectorToIndex;     /**< maps SCIP_VAR* to the corresponding
                                                                                   *< index */
   unordered_map<DEC_DETECTOR*, int> scipPostprocessingDetectorToIndex;     /**< maps SCIP_VAR* to the corresponding
                                                                                   *< index */

   unordered_map<std::pair<int, int>, SCIP_Real, pair_hash> valsMap;   /**< maps an entry of the matrix to its
                                                                                   *< value, zeros are omitted */

   std::vector<SCIP_VAR*> unpresolvedfixedtozerovars;                            /**< helping data structure to collect SCIP_VAR* that are fixed to yero in transformed problem*/

   int nVars;                    /**< number of variables */
   int nConss;                   /**< number of constraints */
   int nDetectors;               /**< number of detectors */
   int nFinishingDetectors;      /**< number of finishing detectors */
   int nPostprocessingDetectors; /**< number of postprocessing detectors */
   int nnonzeros;                /**< number of nonzero entries in the coefficient matrix */


   /** oracle data */
   std::vector<int> usercandidatesnblocks;               /**< candidate for the number of blocks that were given by the user and thus will be handled priorized */
   std::vector<std::pair<int, int>> candidatesNBlocks;   /**< candidate for the number of blocks  */

   std::vector<SCIP_Real>        dualvalsrandom;                                      /**< vector of random dual values, used for strong detection scores */
   std::vector<SCIP_Real>        dualvalsoptimaloriglp;                               /**< vector of dual values of the optimal solved original lp */
   SCIP_Bool                     dualvalsrandomset;                                   /**< are the random dual values set */
   SCIP_Bool                     dualvalsoptimaloriglpcalculated;                     /**< are the optimal dual values from original lp calulated? */

   SCIP_Bool transformed;                                /**< corresponds the matrix data structure to the transformed
                                                         *< problem */
   SCIP_Bool benders;                                    /**< is this supposed to find decompositions for benders decomposition */

   std::vector<SeeedPtr> seeedstopopulate;               /**< seeeds that are translated seeeds from found ones for the
                                                         *< original problem */


   /**
    * @brief method to calculate and set the optimal dual values from original lp, used for strong detection score
    * @return scip return code
    */
   SCIP_RETCODE calculateDualvalsOptimalOrigLP();

   /**
    * @brief method that shuffles randomly and set dual variable values, used for strong detection score
    * @return scip return code
    */
   SCIP_RETCODE shuffleDualvalsRandom();

public:

   std::vector<ConsClassifier*> consclassescollection;   /**< collection of different constraint class distributions  */
   std::vector<VarClassifier*> varclassescollection;     /**< collection of different variable class distributions   */

   SCIP_Real classificationtime;                         /**< time that was consumed by the classification of the constraint and variables classifiers */
   SCIP_Real nblockscandidatescalctime;                  /**< time that was used to calulate the candidates of te block number */
   SCIP_Real postprocessingtime;                         /**< time that was spent in postproceesing decomposigtions */
   SCIP_Real scorecalculatingtime;                       /**< time that was spent by calculating scores */
   SCIP_Real translatingtime;                            /**< time that was spent by transforming seeeds between presolved and unpresolved problem */


   /**
    * @brief constructor
    * @param scip SCIP data structure
    * @param conshdlrName  name of the conshandler maintaining the seeedpool, should be "decomp"
    * @param transformed true if the seeedpool is created for the presolved version of the problem
    * @param benders true if the seeedpool is created for for detecting benders decompositions
    */
   /** constructor */
   Seeedpool(
      SCIP* scip,                /**< SCIP data structure */
      const char* conshdlrName,  /**< name of the conshandler maintaining the seeedpool */
      SCIP_Bool transformed,      /**< true if the seeedpool is created for the presolved version of the problem */
      SCIP_Bool benders           /**< true if the seeedpool is created for for detecting benders decompositions*/
      );

   /**
    * destructor
    */
   ~Seeedpool();


   /**
    * @brief creates constraint and variable classifiers, and deduces block number candidates
    * @param givenScip SCIP data structure
    * @return scip return code
    */
   SCIP_RETCODE calcClassifierAndNBlockCandidates(
      SCIP* givenScip /**< SCIP data structure */
      );


   /**
    * @brief create the constraint adjacency datastructure that is used (if created) for some methods to faster access the constarints that have variables in common
    */
   void createConssAdjacency();



   /**
    * @brief  constructs seeeds using the registered detectors
    * @return vector of found seeeds, that has to be freed by the caller of findSeeeds()
    */
   std::vector<SeeedPtr> findSeeeds();



   /**
    * @brief sorts seeeds in finished seeeds data structure according to the scoretype that is currently used
    */
   void sortFinishedForScore();


   /**
    * @brief method to complete a set of incomplete seeeds with the help of all included detectors that implement a finishing method
    * @param incompleteseeeds the set of incompleted seeeds
    * @return  vector of completed decompositions, has to be freed by the caller
    */
   std::vector<SeeedPtr> finishIncompleteSeeeds(
      std::vector<SeeedPtr> incompleteseeeds /**< the set of incompleted seeeds */
      );


   /**
    * @brief calls findSeeeds method and sort the finished decompositions by score
    */
   void findDecompositions();



   /**
    * @brief  returns seeed with the corresponding id or NULL if it is not found in this seeedpool
    * @param seeedid id of seed is looked for
    * @return a seeed with id as id or NULL no such seeed exists
    */
   gcg::Seeed* findFinishedSeeedByID(
      int      seeedid
      );



   /**
    * @brief adds a seeed to ancestor seeeds
    * @param seeed taht is added to the ancesotr seeeds
    */
   void addSeeedToAncestor(
      SeeedPtr seeed
      );

   /**
    * @brief adds a seeed to current seeeds (data structure for seeeds that are goin to processed in the propagation rounds)
    * @param seeed pointer of seeed to be added
    */
   void addSeeedToCurr(
      SeeedPtr seeed
      );


   /**
    * @brief adds a seeed to finished seeeds
    * @param seeed pointer of seeed that is going to be added to the finished seeeds (data structure to carry finished decompositions)
    * @param success pointer that is set to TRUE iff the the seeeds was successfully added (i.e. it is no duplicate of a known seeed)
    * @see addSeeedToFinishedUnchecked()
    */
   void addSeeedToFinished(
      SeeedPtr seeed,
      SCIP_Bool* success
      );


   /**
    * @brief adds a seeed to finished seeeds without checking for duplicates, dev has to check this on his own
    * @param seeed pointer of seeed that is going to be added unchecked to the finished seeeds (data structure to carry finished decompositions)
    * @see addSeeedToFinished()
    */
   void addSeeedToFinishedUnchecked(
      SeeedPtr seeed
      );

   /**
    * @brief adds a seeed to incomplete seeeds
    * @param seeed to be added
    * @param success set in method to TRUE iff the seeeds was added successfully (i.e. it is no duplicate of the existing seeeds)
    */
   void addSeeedToIncomplete(
      SeeedPtr seeed,
      SCIP_Bool* success
      );


   /**
    * @brief check if seeed is a duplicate of an existing incomplete seeed
    * @param seeed seeed to be checked
    * @return TRUE iff seeed is a duplicate of an existing incomplete seeed
    */
   SCIP_Bool isSeeedDuplicateofIncomplete(
      SeeedPtr seeed
      );


   /**
    * @brief check if seeed is a duplicate of an existing finished seeed
    * @param seeed seeed to be checked
    * @return TRUE iff seeed is a duplicate of an existing finished seeed
    */
   SCIP_Bool isSeeedDuplicateofFinished(
      SeeedPtr seeed
      );


   /**
    * @brief calculates the strong decomposition score of a seeed
    * @param seeed pointer of seeed the score is calculated for
    * @param score pointer to store the calculated score
    * @return scip return code
    */
   /** calculates the strong decomposition score of a seeed */
   SCIP_RETCODE calcStrongDecompositionScore(
      SeeedPtr seeed,
      SCIP_Real* score
      );



   /**
    * @brief checks if there are continouos variables
    * @return TRUE iff there are continouos variables
    */
   SCIP_Bool areThereContinuousVars();



   /**
    * @brief clears ancestor seeed data structure,
    * @note does not free the seeeds themselves
    */
   void clearAncestorSeeeds();


   /**
    * @brief clears current seeed data structure
    * @note does not free the seeeds themselves
    */
   void clearCurrentSeeeds();


   /**
    * @brief clears finished seeed data structure
    * @note does not free the seeeds themselves
    */
   void clearFinishedSeeeds();


   /**
    * @brief clears incomplete seeed data structure
    * @note does not free the seeeds themselves
    */
   void clearIncompleteSeeeds();


   /**
    * @brief returns a seeed from ancestor seeed data structure with given index
    * @param seeedindex index of seeed in ancestor seeed data structure
    * @return seeed from ancestor seeed data structure
    */
   SeeedPtr getAncestorSeeed(
      int seeedindex /**< index of seeed in ancestor seeed data structure */
      );


   /**
    * returns a seeed from current (open) seeed data structure
    * @param seeedindex index of seeed in current seeed data structure
    * @return  index of seeed in current (open) seeed data structure
    */
   /**  */
   SeeedPtr getCurrentSeeed(
      int seeedindex /**< index of seeed in current (open) seeed data structure */
      );


   /**
    * @brief returns a seeed from finished seeed data structure
    * @param seeedindex index of seeed in finished seeed data structure
    * @return  seeed from finished seeed data structure
    */
   /** returns a seeed from finished seeed data structure */
   SeeedPtr getFinishedSeeed(
      int seeedindex /**< index of seeed in finished seeed data structure */
      );


   /**
    * @brief returns a seeed from incomplete seeed data structure
    * @param seeedindex index of seeed in incomplete seeed data structure
    * @return  a seeed from incomplete seeed data structure
    */
   SeeedPtr getIncompleteSeeed(
      int seeedindex /**< index of seeed in incomplete seeed data structure */
      );


   /**
    * @brief returns size of ancestor seeed data structure
    * @return size of ancestor seeed data structure
    */
   int getNAncestorSeeeds();



   /**
    * @brief returns size of current (open) seeed data structure
    * @return size of current (open) seeed data structure
    */
   int getNCurrentSeeeds();



   /**
    * returns size of finished seeed data structure
    * @return  size of finished seeed data structure
    */
   int getNFinishedSeeeds();


   /**
    * returns size of incomplete seeed data structure
    * @return size of incomplete seeed data structure
    */
   int getNIncompleteSeeeds();


   /**
    * returns total number of constraints ( ranged constraints (lb \neq ub) are counted twice )
    * @return total number of constraints where ranged constraints are counted twice
    */
   int getNTotalConss(
   );



   /**
    * @brief returns the number of nonzero entries in the coefficient matrix
    * @return the number of nonzero entries in the coefficient matrix
    */
   long getNTotalNonzeros();



   /**
    * @brief returns true if the given seeed is a duplicate of a seeed that is already contained in finished seeeds or current seeeds data structure
    * @param seeed pointer of seeed that is to check
    * @return true iff the given seeed is a duplicate of a seeed that is already contained in
    *  finished seeeds or current seeeds data structure
    */
   bool hasDuplicate(
      SeeedPtr seeed
      );


   /**
    * @brief returns whether or not this seeedpool is for detectiong decompositions for benders
    * @return whether or not this seeedpool is for detectiong decompositions for benders
    */
   SCIP_Bool isForBenders();


   /**
    * @brief  translates seeeds and classifiers if the index structure  of the problem has changed, e.g. due to presolving
    * @param otherpool old seeedpool
    * @param otherseeeds vector of seeeds to be translated
    * @param newseeeds translated seeeds (caller should pass empty vector)
    * @param otherconsclassifiers consclassifiers to be translated
    * @param newconsclassifiers ranslated consclassifiers (pass empty vector)
    * @param othervarclassifiers varclassifiers to be translated
    * @param newvarclassifiers translated varclassifiers (pass empty vector)
    */
    void translateSeeedData(
      Seeedpool* otherpool,                              /**< old seeedpool */
      std::vector<Seeed*> otherseeeds,                   /**< seeeds to be translated */
      std::vector<Seeed*>& newseeeds,                    /**< translated seeeds (pass empty vector) */
      std::vector<ConsClassifier*> otherconsclassifiers, /**< consclassifiers to be translated */
      std::vector<ConsClassifier*>& newconsclassifiers,  /**< translated consclassifiers (pass empty vector) */
      std::vector<VarClassifier*> othervarclassifiers,   /**< varclassifiers to be translated */
      std::vector<VarClassifier*>& newvarclassifiers     /**< translated varclassifiers (pass empty vector) */
      );


    /**
     * @brief translates seeeds if the index structure of the problem has changed, e.g. due to presolving
     * @param otherpool old seeedpool
     * @param otherseeeds seeeds to be translated
     * @param newseeeds  translated seeeds (pass empty vector)
     */
   /** translates seeeds if the index structure of the problem has changed, e.g. due to presolving */
   void translateSeeeds(
      Seeedpool* otherpool,            /**< old seeedpool */
      std::vector<Seeed*> otherseeeds, /**< seeeds to be translated */
      std::vector<Seeed*>& newseeeds   /**< translated seeeds (pass empty vector) */
      );


   /**
    * @brief registers seeeds, e.g. done with translated seeeds from the original problem
    * @param seeeds vector of seeeds that are about to be populated
    */
   void populate(
      std::vector<SeeedPtr> seeeds
      );

   /**
    * @brief registers a seeed for a seeedpool, sets it up by and calculates its hash value
    * @param seeed pointer f seeed to be prepared
    * @return scip return code
    */
   SCIP_RETCODE prepareSeeed(
      SeeedPtr seeed
      );

   /**
    * @brief sorts seeeds in allrelevantseeeds data structure by ascending id
    * @note needed for fast access to ancestors
    */
   void sortAllRelevantSeeeds();


   /**
    * @brief returns whether a constraint is a cardinality constraint, i.e. of the \sum){i} x_i == b
    * @param consindexd index of constraint that is be checked
    * @return returns whether a constraint is a cardinality constraint
    */
   bool isConsCardinalityCons(
         int  consindexd
         );

   /**
    * @brief is cons with specified index partitioning packing, or covering constraint?
    * @param consindex indec of cons to be checked
    * @return  whether a constraint  is partitioning packing, or covering constraint?
    */
   bool isConsSetppc(
      int  consindexd
      );


   /**
    * @brief is cons with specified indec partitioning, or packing covering constraint?
    * @param consindexd index of the given cons
    * @return is cons with specified indec partitioning, or packing covering constraint
    */
   bool isConsSetpp(
      int  consindexd
      );


   /**
    * @brief returns the variable indices of the coefficient matrix for a constraint
    * @param consIndex index of the constraint to be considered
    * @return the variable indices of the coefficient matrix for a constraint
    */
   /** returns the variable indices of the coefficient matrix for a constraint */
   const int* getVarsForCons(
      int consIndex /**< index of the constraint to be considered */
      );

   /**
    * @brief returns the nonzero coefficients of the coefficient matrix for a constraint
    * @param consIndex  index of the constraint to be considered
    * @return array of coefficients of in matrix for constraints
    * @note same order as in @see getVarsForCons()
    */
   const SCIP_Real* getValsForCons(
      int consIndex /**< index of the constraint to be considered */
      );


   /**
    * \brief returns the constraint indices of the coefficient matrix for a variable
    * @param varIndex index of the variable to be considered
    * @return vector of constraint indices that have a nonzero entry with this variable
    */
   const int* getConssForVar(
      int varIndex /**< index of the variable to be considered */
      );


/**
 * @brief returns the value of the optimal lp relaxation dual value of the given constrainr rid correspondoning problem of the seeedpool; if it is not calculated yet it will be calculated
 * @param consindex index of constraint the value is asked for
 * @return the value of the optimal lp relaxation dual value of the given constrainr rid correspondoning problem of the seeedpool
 */
   SCIP_Real  getDualvalOptimalLP(
      int  consindex
      );

   /**
    * @brief return the a random value of the dual variable of the corresponding ; if it is not calculated yet it will be calculated
    * @param consindex  index of constraint the value is asked for
    * @return  the a random value of the dual variable of the corresponding
    */
   SCIP_Real  getDualvalRandom(
      int  consindex
   );


   /**
    * @brief returns the number of variables for a given constraint
    * @param consIndex  index of the constraint to be considered
    * @return the number of variables for a given constraint
    */
   /** returns the number of variables for a given constraint */
   int getNVarsForCons(
      int consIndex /**< index of the constraint to be considered */
      );

   /**
    * returns the number of constraints for a given variable where the var has a nonzero entry in
    * @param varIndex index of the variable to be considered
    * @return the number of constraints for a given variable
    */
   /** returns the number of constraints for a given variable */
   int getNConssForVar(
      int varIndex /**< index of the variable to be considered */
      );

   /**
    * return array of constraint indices that have a common variable with the given constraint
    * @param consIndex index of the constraint the neighboring constraints are
    * @return return array of constraint indices that have a common variable with the given constraint
    * @note constraint adjacency data structure has to initilized
    */
   const int* getConssForCons(
      int consIndex /**< index of the constraint to be considered */
      );


   /**
    * @brief returns the number of constraints for a given constraint
    * @param consIndex index of the constraint to be considered
    * @return the number of constraints for a given constraint
    */
   int getNConssForCons(
      int consIndex /**< index of the constraint to be considered */
      );


   /**
    * @brief returns SCIP variable related to a variable index
    * @param varIndex index of the variable to be considered
    * @return SCIP variable pointer related to a variable index
    */
   SCIP_VAR* getVarForIndex(
      int varIndex /**< index of the variable to be considered */
      );


   /**
    * @brief returns the SCIP constraint related to a constraint index
    * @param consIndex index of the constraint to be considered
    * @return the SCIP constraint related to a constraint index
    */
   SCIP_CONS* getConsForIndex(
      int consIndex /**< index of the constraint to be considered */
      );

   /**
    * @brief returns the detector related to a detector index
    * @param detectorIndex index of the detector to be considered
    * @return pointer to detector related to a detector index
    */
   DEC_DETECTOR* getDetectorForIndex(
      int detectorIndex /**< index of the detector to be considered */
      );

   /**
    * @brief returns the detector related to a finishing detector index
    * @param detectorIndex index of the finishing detector to be considered
    * @return detector pointer related to a finishing detector index
    */
   DEC_DETECTOR* getFinishingDetectorForIndex(
      int detectorIndex /**< index of the finishing detector to be considered */
      );


   /**
    * @brief returns the detector related to a finishing detector index
    * @param detectorIndex index of the postprocessing detector to be considered
    * @return detector pointer related to a postprocessing detector index
    */
   DEC_DETECTOR* getPostprocessingDetectorForIndex(
      int detectorIndex /**< index of the postprocessing detector to be considered */
      );


   /**
    * @brief returns a coefficient from the coefficient matrix
    * @param row index of the constraint to be considered
    * @param col index of the variable to be considered
    * @return a coefficient from the coefficient matrix
    */
   SCIP_Real getVal(
      int row, /**< index of the constraint to be considered */
      int col  /**< index of the variable to be considered */
      );

   /**
    * @brief returns the variable index related to a SCIP variable
    * @param var variable pointer the index is asked for
    * @return the variable index related to a SCIP variable
    */
   int getIndexForVar(
      SCIP_VAR* var
      );


   /**
    * @brief returns the constraint index related to a SCIP constraint
    * @param cons the SCIP constraint pointer the index is asked for
    * @return the constraint index related to a SCIP constraint
    */
   int getIndexForCons(
      SCIP_CONS* cons
      );


   /**
    * @brief the detector index related to a detector
    * @param detector pointer of detector the index is asked for
    * @return the detector index related to a detector
    */
   int getIndexForDetector(
      DEC_DETECTOR* detector
      );


   /**
    * @brief returns the finishing detector index related to a detector
    * @param detector pointer of finishing detector
    * @return the finishing detector index
    */
   int getIndexForFinishingDetector(
      DEC_DETECTOR* detector
      );


   /**
    * @brief returns the postprocessing detector index related to a detector
    * @param detector pointer to the postprocessing detector
    * @return the postprocessing detector index
    */
   int getIndexForPostprocessingDetector(
      DEC_DETECTOR* detector
      );


   /**
    * @brief returns a new unique id for a new seeed
    * @return  a new integer unique id for a seeed
    */
   int getNewIdForSeeed();


   /**
    * @brief returns the number of propagating detectors used in the seeedpool
    * @return number of detectors
    */
   int getNDetectors();

   /**
    * @brief returns the number of nonzero entries in the coefficient matrix
    * @return the number of nonzero entries in the coefficient matrix
    */
   int getNNonzeros();

   /**
    * @brief returns the number of finishing detectors used in the seeedpool
    * @return  the number of finishing detectors used in the seeedpool
    */
   int getNFinishingDetectors();


   /**
    * @brief returns the number of postprocessing detectors used in the seeedpool
    * @return the number of postprocessing detectors used in the seeedpool
    */
   int getNPostprocessingDetectors();


   /**
    * @brief return the number of variables considered in the seeedpool
    * @return the number of variables considered in the seeedpool
    */
   int getNVars();


   /**
    * @brief returns the number of variables considered in the seeedpool
    * @return number of variables considered in the seeedpool
    */
   int getNConss();


   /**
    * @brief returns the corresponding scip data structure
    * @return the corresponding scip data structure
    */
   SCIP* getScip();


   /**
    * @brief returns scip cons for corresponing id
    * @param consid the index of the constraint
    * @return pointer of constraint for the given index
    */
   SCIP_CONS* getScipCons(
      int consid
      );


   /**
    * @brief returns scip var for corresponding id
    * @param varid the index of the variable
    * @return returns scip var for corresponding id
    */
   SCIP_VAR* getScipVar(
      int varid
   );


   /**
    * @brief returns the candidates for number of blocks added by the user followed by the found ones sorted in descending order by how often a candidate was proposed  )
    * @note @see getSortedCandidatesNBlocksFull()
    * @return the candidates for number of blocks sorted in descending order by how often a candidate was added
    */
   std::vector<int> getSortedCandidatesNBlocks();

   /**
    * @brief returns (candidate-nvotes) pair for candidates for number of blocks added by the user followed by the found ones sorted in descending order by how often a candidate was proposed, candidates are first   )
    * @note @see getSortedCandidatesNBlocks()
    * @return vector of (candidate-nvotes) pair  the candidates for block size sorted in descending order by how often a candidate was added with nvotes information*
    */
    std::vector<std::pair<int, int>> getSortedCandidatesNBlocksFull();


    /**
     * @brief adds and counts how often added a candidate for block size and
     * @param candidate number of candidate
     */
   void addCandidatesNBlocks(
      int candidate /**< candidate for block size */
      );


   /**
    * @brief adds a candidate for block size and counts how often a candidate is added
    * @param candidate candidate for block size
    * @param nvotes number of votes this candidates will get
    */
   void addCandidatesNBlocksNVotes(
      int candidate, /**< candidate for block size */
      int nvotes     /**< number of votes this candidates will get */
      );


   /**
    * @brief adds a candidate for block size given by the user
    * @param candidate user candidate for block size
    */
   void addUserCandidatesNBlocks(
      int candidate /**< candidate for block size */
      );


   /**
    * @brief returns number of user-given block size candidates
    * @return number of user-given block size candidates
    */
   int getNUserCandidatesNBlocks();


   /**
    * @brief calculates and adds block size candidates using constraint classifications and variable classifications
    */
   void calcCandidatesNBlocks();

   /**
    * @brief adds a constraint classifier if it is no duplicate of an existing constraint classifier
    * @param classifier pointer to consclassifier to be added
    */
   void addConsClassifier(
      ConsClassifier* classifier /**< consclassifier to be added*/
      );


   /**
    * @brief returns a new constraint classifier
    *  where all constraints with identical SCIP constype are assigned to the same class
    * @return pointer to constraint classifier
    */
   ConsClassifier* createConsClassifierForSCIPConstypes();


   /**
    * @brief returns a new constraint classifier
    *  where all constraints with identical Miplib constype are assigned to the same class
    * @return a new constraint classifier according to miplib 2010 constype
    */
   ConsClassifier* createConsClassifierForMiplibConstypes();


   /**
    * @brief returns a new constraint classifier
    *  where all constraints with identical consname (ignoring digits) are assigned to the same class
    * @return new constraint classifer according to consname (ignoring digits)
    */
   ConsClassifier* createConsClassifierForConsnamesDigitFreeIdentical();


   /**
    * @brief returns a new constraint classifier where all constraints whose consnames do not a have levenshtein distance to each other higher than a given connectivity are assigned to the same class
    * @param connectivity minimum levenshtein distance for two consnames to be incident
    * @return new constraint classifier according to levenshtein distance of consnames
    */
   ConsClassifier* createConsClassifierForConsnamesLevenshteinDistanceConnectivity(
      int connectivity /**< given connectivity */
      );


   /**
    * returns a new constraint classifier where all constraints with identical number of nonzero coefficients are assigned to the same class
    * @return a new constraint classifier where all constraints with identical number of nonzero coefficients are assigned to the same class
    */
   ConsClassifier* createConsClassifierForNNonzeros();


   /**
    * @brief returns pointer to a constraint classifier
    * @param classifierIndex index of constraint classifier
    * @return pointer to a cosntraint classifier with the given index
    */
   ConsClassifier* getConsClassifier(
      int classifierIndex /**< index of constraint classifier */
      );


   /**
    * @brief returns the assignment of constraints to classes of a classifier as integer array
    * @param classifierIndex index of the constraint classifier the assignment i asked for
    * @return array conatining for each constraint id the assigned class of the constraint for the classifier of the given index
    */
   int* getConsClassifierArray(
      int classifierIndex /**< index of constraint classifier */
      );

   /**
    * @brief returns number of different constraint classifiers
    * @return number of different constraint classifiers
    */
   int getNConsClassifiers();


   /**
    * @brief adds constraint classifiers with a reduced number of classes
    */
   void reduceConsclasses();


   /**
    * @brief  adds a variable classifier if it is no duplicate of an existing variable classifier
    * @param classifier varclassifier to be added
    */
   void addVarClassifier(
      VarClassifier* classifier /**< */
      );

   /**
    * @brief returns a new variable classifier where all variables with identical objective function value are assigned to the same class
    * @return var classififier according to objective function coefficient
    */
   VarClassifier* createVarClassifierForObjValues();


   /**
    * @brief returns a new variable classifier where all variables are assigned to class zero, positive or negative according to their objective function value sign all class zero variables are assumed to be only master variables (set via DECOMPINFO)
    * @return var classififier according to objective function coefficient sign
    */
   VarClassifier* createVarClassifierForObjValueSigns();

   /**
    * @brief returns a new variable classifier where all variables with identical SCIP vartype are assigned to the same class
    * @return variable classifier concerning SCIP vartype
    */
   VarClassifier* createVarClassifierForSCIPVartypes();

   /**
    * @brief returns number of different variable classifiers
    * @return  number of different variable classifiers
    */
   int getNVarClassifiers();

   /**
    * @brief returns pointer to a variable classifier with given index
    * @param classifierIndex index of variable classifer
    * @return pointer to a variable classifier with given index
    */
   VarClassifier* getVarClassifier(
      int classifierIndex /**< index of variable classifier */
      );


   /**
    * @brief  returns the assignment of variables to classes of a classifier as integer array
    * @param classifierIndex index of the variables  classifier the assignment vector is requested for
    * @return array of class indices for the variable indices
    */
   int* getVarClassifierArray(
      int classifierIndex /**< index of constraint classifier */
      );

   /**
    * @brief adds variable classifiers with a reduced number of classes
    * this is done by greedliy merge two smallest classes into one new until limit from settings is requested
    */
   void reduceVarclasses();


   /**
    * @brief returns a vector of seeeds such that all seeeds of given a vector of seeeds having only one block are removed except for the two seeeds with the lowest numbers of masterconss
    * @param givenseeeds vector of seeeds to be reduced
    * @return the reduced vector of seeeds
    */
   /** returns a vector of seeeds where all seeeds of given seeeds having only one block are removed
    *  except for the two seeeds with the lowest numbers of masterconss */
   std::vector<SeeedPtr> removeSomeOneblockDecomps(
      std::vector<SeeedPtr> givenseeeds
      );


   /**
    * @brief creates a decomposition DEC_DECOMP structure for a given seeed
    * @param seeed seeed the decomposition is created for
    * @param newdecomp the new decomp created from the seeed
    * @return scip return code
    */
   SCIP_RETCODE createDecompFromSeeed(
      SeeedPtr seeed,         /** seeed the decomposition is created for */
      DEC_DECOMP** newdecomp  /** the new decomp created from the seeed */
      );


   /**
    *  creates a seeed for a given decomposition
    *  the resulting seeed will not have a detectorchaininfo or any ancestor or finishing detector data
    *  only use this method if the seeedpool is for the transformed problem
    *  the resulting seeed may only be added to the seeedpool for the presolved problem
    * @param decomp decomposition the seeed is created for
    * @param newseeed the new seeed created from the decomp
    * @return scip return code
    */
    SCIP_RETCODE createSeeedFromDecomp(
      DEC_DECOMP* decomp, /**< decomposition the seeed is created for */
      SeeedPtr* newseeed /**< the new seeed created from the decomp */
      );

    /**
     * @brief returns true if the matrix structure corresponds to the transformed problem
     * @return TRUE if the matrix structure corresponds to the transformed problem
     */
   SCIP_Bool getTransformedInfo();


   /**
    * @brief output method for json file writer to write block candidate information
    * @param scip SCIP data structure
    * @param file  output file or NULL for standard output
    * @return scip return code
    */
   SCIP_RETCODE printBlockcandidateInformation(
    SCIP*                 scip,               /**< SCIP data structure */
    FILE*                 file                /**< output file or NULL for standard output */
   );


   /**
    * @brief output method for json file writer to write classifier candidate information
    * @param scip SCIP data structure
    * @param file output file or NULL for standard output
    * @return scip return code
    */
   SCIP_RETCODE printClassifierInformation(
    SCIP*                 scip,               /**< SCIP data structure */
    FILE*                 file                /**< output file or NULL for standard output */
   );


private:


   /**
    * @brief calculates necessary data for translating seeeds and classifiers
    * @param origpool original seeedpool
    * @param rowothertothis constraint index mapping from old to new seeedpool
    * @param rowthistoother constraint index mapping new to old seeedpool
    * @param colothertothis variable index mapping from old to new seeedpool
    * @param colthistoother variable index mapping from new to old seeedpool
    * @param missingrowinthis missing constraint indices in new seeedpool
    */
   /** calculates necessary data for translating seeeds and classifiers */
   void calcTranslationMapping(
      Seeedpool* origpool, /** original seeedpool */
      std::vector<int>& rowothertothis,   /** constraint index mapping from old to new seeedpool */
      std::vector<int>& rowthistoother,   /** constraint index mapping new to old seeedpool */
      std::vector<int>& colothertothis,   /** variable index mapping from old to new seeedpool */
      std::vector<int>& colthistoother,   /** variable index mapping from new to old seeedpool */
      std::vector<int>& missingrowinthis  /** missing constraint indices in new seeedpool */
      );


   /**
    * @brief returns translated eeeds derived from given mapping data
    * @param otherseeeds seeeds to be translated
    * @param rowothertothis constraint index mapping from old to new seeedpool
    * @param rowthistoother constraint index mapping new to old seeedpool
    * @param colothertothis variable index mapping from old to new seeedpool
    * @param colthistoother variable index mapping from new to old seeedpool
    * @return vector of translated seeed pointers
    */
    std::vector<Seeed*> getTranslatedSeeeds(
      std::vector<Seeed*>& otherseeeds,   /**< seeeds to be translated */
      std::vector<int>& rowothertothis,   /** constraint index mapping from old to new seeedpool */
      std::vector<int>& rowthistoother,   /** constraint index mapping new to old seeedpool */
      std::vector<int>& colothertothis,   /** variable index mapping from old to new seeedpool */
      std::vector<int>& colthistoother    /** variable index mapping from new to old seeedpool */
      );


    /**
     * @brief returns translated ConsClassifiers derived from given mapping data
     * @param otherclassifiers consclassifiers to be translated
     * @param rowothertothis constraint index mapping from old to new seeedpool
     * @param rowthistoother constraint index mapping new to old seeedpool
     * @return vector of translated ConsClassifier
     */
   std::vector<ConsClassifier*> getTranslatedConsClassifiers(
      std::vector<ConsClassifier*>& otherclassifiers, /**< consclassifiers to be translated */
      std::vector<int>& rowothertothis,   /** constraint index mapping from old to new seeedpool */
      std::vector<int>& rowthistoother    /** constraint index mapping new to old seeedpool */
      );

   /**
    * @brief returns translated VarClassifiers derived from given mapping data
    * @param otherclassifiers varclassifiers to be translated
    * @param colothertothis variable index mapping from old to new seeedpool
    * @param colthistoother variable index mapping from new to old seeedpool
    * @return vector of translated VarClassifier
    */
   std::vector<VarClassifier*> getTranslatedVarClassifiers(
      std::vector<VarClassifier*>& otherclassifiers, /**< varclassifiers to be translated */
      std::vector<int>& colothertothis,   /** variable index mapping from old to new seeedpool */
      std::vector<int>& colthistoother    /** variable index mapping from new to old seeedpool */
      );



};
/* class Seeedpool */



} /* namespace gcg */
#endif /* GCG_CLASS_SEEEDPOOL_H__ */
