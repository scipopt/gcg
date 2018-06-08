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

/**@file   class_Seeed.h
 * @brief  class with functions for seeed
 * @author Michael Bastubbe
 * @author Hannah Hechenrieder
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEED_H__
#define GCG_CLASS_SEEED_H__

#include "objscip/objscip.h"
#include "struct_detector.h"
#include "cons_decomp.h"

#include <vector>
#include <string>

#include "class_consclassifier.h"
#include "class_varclassifier.h"
#include "graph/graph_gcg.h"
#include "graph/graph.h"

namespace gcg {


/**
 * @brief enumeration to display if a decomposition was given by the user and if so, how it was processed after adding
 */
enum USERGIVEN
{
   NOT = 0,                            /**< this seeed was not given by the user */
   PARTIAL = - 1,                      /**< this partial seeed was given by the user as it is*/
   COMPLETE = - 2,                     /**< this complete seeed was given by the user as it is*/
   COMPLETED_CONSTOMASTER = - 3        /**< this seeed was partially given by the user and then completed by setting all missing constraints to the master*/
};

class Seeedpool;


/*!
 * @brief class to manage partial decompositions (aka seeed), each seeed corresponds to one seeedpool which contains the problem information, there is one seeedpool for the original and the transformed problem.
 */
class Seeed
{
private:
   SCIP* scip;                                                 /**< SCIP data structure */
   int id;                                                     /**< unique id of the seeed, unique */
   int nBlocks;                                                /**< number of blocks the partial decomposition currently has */
   int nVars;                                                  /**< number of variables */
   int nConss;                                                 /**< number of constraints */
   std::vector<int> masterConss;                               /**< vector containing indices of master constraints */
   std::vector<int> masterVars;                                /**< vector containing indices of master variables (these variables are supposed to have all nonzero entries in master constraints) */
   std::vector<std::vector<int>> conssForBlocks;               /**< conssForBlocks[k] contains a vector of indices of all
                                                                 *< constraints assigned to block k */
   std::vector<std::vector<int>> varsForBlocks;                /**< varsForBlocks[k] contains a vector of indices of all
                                                                 *< variables assigned to block k */
   std::vector<int> linkingVars;                               /**< vector containing indices of linking variables */
   std::vector<std::vector<int>> stairlinkingVars;             /**< vector containing indices of staircase linking variables
                                                                 *< of the blocks (stair-linking variables are registered only
                                                                 *< in their first block) */
   std::vector<int> openVars;                                  /**< vector containing indices of variables that are not
                                                                 *< assigned yet*/
   std::vector<int> openConss;                                 /**< vector containing indices of constraints that are not
                                                                 *< assigned yet*/
   std::vector<bool> isvaropen;                                /**< help vector for fast query if a variable is still open */
   std::vector<bool> isconsopen;                               /**< help vector for fast query if a constraint is still open */
   std::vector<bool> isvarmaster;                              /**< help vector for fast query if a variable is assigned to be a only-master variable */
   std::vector<bool> isconsmaster;                             /**< help vector for fast query if a constraint is assigned to be a master constraint */

   std::vector<int>  ncoeffsforblock;                          /**< number of coeffs per block */

   SCIP_Bool         calculatedncoeffsforblock;                /**< is the  number of coeff per block already calculated*/
   int               ncoeffsformaster;                         /**< number of master coefficients */
   std::vector<std::vector<int>> ncoeffsforblockformastercons; /**< number of coeffs a block has in a certain master constraint */

   bool varsforblocksorted;                                    /**< bool to store if the varsforblocks datastructure is sorted atm */
   bool stairlinkingvarsforblocksorted;                        /**< bool to store if the stairlinkingvarsforblock datastructure is sorted atm */
   bool conssforblocksorted;                                   /**< bool to store if the conssforblock datastructure is sorted atm */
   bool linkingvarssorted;                                     /**< bool to store if the linkingvars datastructure is sorted atm */
   bool mastervarssorted;                                      /**< bool to store if the mastervars datastructure is sorted atm */
   bool masterconsssorted;                                     /**< bool to store if the masterconsssorted datastructure is sorted atm */


   std::vector<int> bookedAsMasterConss;                       /**< vector containing indices of constraints that are not
                                                                 *< assigned yet but booked as master conss */
   std::vector<std::pair<int, int>> bookedAsBlockConss;        /**< vector containing indices of constraints that are not
                                                                 *< assigned yet but booked as block conss and the block */
   std::vector<int> bookedAsLinkingVars;                       /**< vector containing indices of variables that are not
                                                                 *< assigned yet but booked as linking vars */
   std::vector<int> bookedAsMasterVars;                        /**< vector containing indices of variables that are not
                                                                 *< assigned yet but booked as master vars */
   std::vector<std::pair<int, int>> bookedAsBlockVars;         /**< vector containing indices of variables that are not
                                                                 *< assigned yet but booked as block vars and the block */
   std::vector<std::pair<int, int>> bookedAsStairlinkingVars;  /**< vector containing indices of variables that are not
                                                                 *< assigned yet but booked as stairlinking vars and the
                                                                 *< first block of the stairlinking var */
   long hashvalue;                                             /**< hash value of this partial decomposition, decompositions with same has value are considered to be identical */
   bool changedHashvalue;                                      /**< are there any changes concerning the hash value since it
                                                                 *< was calculated last time */

   bool isselected;                                            /**< is this seeed selected */

   bool isagginfoalreadytoexpensive;                            /**< is agginfo already known to be to expensive to calculate*/

   const static int primes[];                                   /**< an array of prime numbers to calculate the hashvalue */
   const static int nPrimes;                                    /**< size of the array of prime numbers */

   bool isFinishedByFinisher;                                   /**< was this seeed finished by the finishseeed() method of a detector */

   /** aggregation information */
   SCIP_Bool            agginfocalculated;                             /**< is aggregation information for the blocks already calculated */
   int                  nrepblocks;                                    /**< number of block representatives */
   std::vector<std::vector<int>> reptoblocks;                          /**< translation of the block representatives to (old) blocks */
   std::vector<int>     blockstorep;                                   /**< translation of the (old) blocks to the block representatives */
   std::vector<std::vector<std::vector<int> > > pidtopidvarmaptofirst; /**< [nrepblocks][blockstorep[k].size()][nvarsforprob] collection of varmaps of probindices from k-th subproblem to the zeroth block that is represented */

   /** statistic information */
   std::vector<DEC_DETECTOR*> detectorChain;          /**< vector containing detectors that worked on that seeed */
   std::vector<std::string> detectorchaininfo;        /**< vector containing information about the detector call */
   std::vector<SCIP_Bool> detectorChainFinishingUsed; /**< vector containing whether the finishing method of the
                                                        *< corresponding detector was used on that seeed */
   std::vector<SCIP_Real> detectorClockTimes;         /**< vector containing detector times in seconds  */
   std::vector<SCIP_Real> pctVarsToBorder;            /**< vector containing the fraction of variables assigned to the
                                                        *< border for each detector working on that seeed*/
   std::vector<SCIP_Real> pctVarsToBlock;             /**< vector containing the fraction of variables assigned to a block
                                                        *< for each detector working on that seeed*/
   std::vector<SCIP_Real> pctVarsFromFree;            /**< vector containing the fraction of variables that are not longer
                                                        *< open for each detector working on that seeed*/
   std::vector<SCIP_Real> pctConssToBorder;           /**< vector containing the fraction of constraints assigned to the
                                                        *< border for each detector working on that seeed*/
   std::vector<SCIP_Real> pctConssToBlock;            /**< vector containing the fraction of constraints assigned to a block
                                                        *< for each detector working on that seeed*/
   std::vector<SCIP_Real> pctConssFromFree;           /**< vector containing the fraction of constraints that are not longer
                                                         *< open for each detector working on that seeed*/
   std::vector<int> nNewBlocks;                       /**< vector containing information how many new blocks a detector has assigned */

   std::vector<IndexClassifier*> usedClassifier;      /**< vector containing pointer to the cons- or varclassifier
                                                         *< a detector made use of for each detector working on that seeed
                                                         *< (NULL if no classifier was used) */
   std::vector<std::vector<int>> classesToMaster;     /**< vector containing the vector of classindices that were assigned
                                                         *< to master by the classifier used by a detector
                                                         *< (empty vector if no classifier was used) */
   std::vector<std::vector<int>> classesToLinking;    /**< vector containing the vector of classindices that were assigned
                                                         *< to linking by the classifier used by a detector
                                                         *< (empty vector if no classifier was used) */

   std::vector<int> listofancestorids;                /**< vector containing decomposition indices that are ancestors of this seeed */


   USERGIVEN usergiven;                               /**< is this seeed partially or completely given by user */
   bool isfromlegacymode;                             /**< true if this seeed stems from a detector operating in legacy mode */
   SCIP_Real score;                                   /**< classc score to evaluate the partial */
   SCIP_Real maxwhitescore;                           /**< score corresponding to the max white measure */
   SCIP_Real bendersscore;                            /**< score to evaluate the seeeds */
   SCIP_Real benderareascore;                         /**< 1 - fraction of white area in master constraints to complete area */

   SCIP_Real strongdecompositionscore;                /**< strong decomposition score  */

   SCIP_Real borderareascore;                         /**< 1 - fraction of border area to complete area */


   SCIP_Real maxwhitescoreagg;                        /**< score corresponding to the max white measure according to aggregated blocks */

   SCIP_Real blockareascore;                          /**< 1 - fraction of block area to complete area */
   SCIP_Real blockareascoreagg;                       /**< 1 - fraction of aggregated block area to complete area */


   SCIP_Real maxforeseeingwhitescore;                /**< maximum foreseeing  white area score (i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking) */
   SCIP_Real maxforeseeingwhitescoreagg;             /**< maximum foreseeing  white area score with respect to aggregatable blocks  (i.e. maximize fraction of white area score considering problem with copied linking variables and corresponding master constraints; white area is nonblock and nonborder area, stairlinking variables count as linking) */

   SCIP_Real setpartfwhitescore;                      /** setpartitioning maximum foreseeing  white area score (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )*/
   SCIP_Real setpartfwhitescoreagg;                   /** setpartitioning maximum foreseeing  white area score with respect to aggregateable  (i.e. convex combination of maximum foreseeing white area score and a boolean score rewarding a master containing only setppc and cardinality constraints )*/


   char* detectorchainstring;                         /**< string formed by the chars of the detectors involved for this seeed  */

   /** datastructure to store information if this seeed stems from a seeed concerning the unpresolved problem */
   bool stemsFromUnpresolved;             /**< seeed has at least one ancestor that is a seeed from unpresolved problem */
   bool isfromunpresolved;                /**< seeed is from unpresolved problem */
   bool isFinishedByFinisherUnpresolved;  /**< was the ancestor seeed for the unpresolved problem finished by the
                                            *< finishseeed() method of a detector */
   DEC_DETECTOR* finishedUnpresolvedBy;   /**< index of finishing detector of unpresolved ancestor seeed */
   Seeedpool*    seeedpool;               /**< seeedpool for the corresponding problem , corresponds either to the presolved or the unpresolved problem*/

   /**
    *
    * @brief checks blocks for identity by graph automorphism check, done by bliss, identity is only found if variables are in correct order
    * @param b1 block id of first block
    * @param b2 block id of second block
    * @param varmap maps variable indices (corresponding to  seeedpool indices) of block 2 to block 1
    * @param varmap2 maps variable pointers of block 2 to those of block 1 if both blocks (problems) are identical
    * @param identical pointer to store if the subproblems are identical
    */
      void checkIdenticalBlocksBliss(
      int                  b1,                     /**< block id of first block */
      int                  b2,                     /**< block id of second block */
      std::vector<int>&    varmap,                 /**< maps variable indices (corresponding to  seeedpool indices) of block 2 to block 1 */
      SCIP_HASHMAP*        varmap2,                /**< maps variable pointers of block 2 to those of block 1 if both blocks (problems) are identical*/
      SCIP_Bool*           identical               /**< pointer to store if the subproblems are identical  */
      );


      /**
       * @brief checks blocks for identity by brute force, identity is only found if variables are in correct order
       * @param b1 block id of first block
       * @param b2 block id of second block
       * @param varmap maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1
       * @param varmap2 maps variable pointers of block 2 to those of block 1 if both blocks (problems) are identical
       * @param identical pointer to store if the subproblems are identical
       */
   void checkIdenticalBlocksBrute(
      int                  b1,                     /**< block id of first block */
      int                  b2,                     /**< block id of second block */
      std::vector<int>&    varmap,                 /**< maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1 */
      SCIP_HASHMAP*        varmap2,                /**< maps variable pointers of block 2 to those of block 1 if both blocks (problems) are identical*/
      SCIP_Bool*           identical               /**< pointer to store if the subproblems are identical  */
      );


   /**
    * @brief check some necessary conditions for two blocks to be identical
    * @param b1 block id of first block
    * @param b2 block id of second block
    * @param notidentical pointer to store whether or not the non-identity is proven
    * @return scip return code
    */
   SCIP_RETCODE checkIdenticalBlocksTrivial(
      int                  b1,                     /**< block id of first block */
      int                  b2,                     /**< block id of second block */
      SCIP_Bool*           notidentical            /**< pointer to store whether or not the non-identity is proven */
      );

public:

   /**
    * constructor
    * @param scip data structure
    * @param id that is given to this seeed
    * @param seeedpool this seeed is created for
    */
   Seeed(
      SCIP* scip,                                  /**< scip data structure */
      int id,                                      /**< id that is given to this seeed */
      Seeedpool* seeedpool                         /**< seeedpool this seeed is created for */
      );

   /**
    * copy constructor
    * @param seeedToCopy seeed to be copied
    */
   Seeed(
      const Seeed *seeedToCopy /**< seeed to be copied */
      );

   /**
    *  destructor
    */
   ~Seeed();


    /**
     *
     * @brief adds a block, returns the number of the new block
     * */
   int addBlock();


   /**
    * @brief incorporates the needed time of a certain detector in the detector chain
    * @param clocktime time to add
    */
   void addClockTime(
      SCIP_Real clocktime /**< time to be added */
      );

   /**
    * @brief incorporates the changes from ancestor seeed into the statistical data structures
    * @param ancestor seeed whose propagation yielded to the current seeed
    */
  void addDecChangesFromAncestor(
      Seeed* ancestor                  /**< seeed whose propagation yielded to the current seeed */
      );

   /**
    * @brief adds a detectorchain information string to the corresponding vector (that carries information for each detector call)
    * @param decinfo information string (about the detector call) to add
    * */
   void addDetectorChainInfo(
      const char* decinfo              /**< information string (about the detector call) to add  */
      );

   /**
    *
    * @brief bookkeeping information: adds number of new blocks created by a detector added to detector chain
    * @param nnewblocks number of new added blocks by latest detector call
    */
   void addNNewBlocks(
      int nnewblocks                   /**< number of new added blocks by latest detector call */
      );

   /**
    * @brief bookkeeping information: fraction of constraints that are not longer open for a detector added to detector chain
    * @param pct fraction of constraints that are not longer open
    */
   void addPctConssFromFree(
      SCIP_Real pct                    /**< fraction of constraints that are not longer open */
      );

   /**
    *  @brief bookkeeping information: adds fraction of constraints assigned to a block for a detector added to detector chain
    * @param pct fraction of constraints assigned to a block
    *  */
   void addPctConssToBlock(
      SCIP_Real pct                    /**< fraction of constraints assigned to a block */
      );

   /**
    *  @brief bookkeeping information: adds fraction of constraints assigned to the border for a detector added to detector chain
    * @param pct fraction of constraints assigned to the border
    */
   void addPctConssToBorder(
      SCIP_Real pct                    /**< fraction constraints assigned to the border */
      );

   /**
    *  @brief bookkeeping information: adds fraction of variables that are not longer open for a detector added to detector chain
    *  @param pct fraction of variables that are not longer open
    */
   void addPctVarsFromFree(
      SCIP_Real pct                    /**< fraction of variables that are not longer open */
      );


   /**
    *  @brief bookkeeping information: adds fraction of variables assigned to a block for a detector added to detector chain
    *  @param pct fraction of variables assigned to a block
    *  */
   void addPctVarsToBlock(
      SCIP_Real pct                     /**< fraction of variables assigned to a block */
      );

   /**
    * @brief bookkeeping information: adds fraction of variables assigned to the border for a detector added to detector chain
    * @param pct fraction of variables assigned to a block
    */
   void addPctVarsToBorder(
      SCIP_Real pct                    /**< fraction of variables assigned to a block */
      );

   /**
    * @brief method to check if at leas one constraint is assigned to some block
    * @returns true if at least one constraint is assigned to a block
    *  */
   bool alreadyAssignedConssToBlocks();


   /**
    *
    */
   /**
    * @brief assigns open conss to master according to the cons assignment information given in constoblock hashmap,
    * @param constoblock hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
    * @param givenNBlocks number of blocks the hashmap contains
    * @return  scip return code
    * \note for conss assigned to blocks according to constoblock there is no assignment \see assignSeeedFromConstoblock
    * */
   SCIP_RETCODE assignBorderFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks          /**< number of blocks the hashmap contains */
       );


   /**
    * @brief assigns open vars to stairlinking if they can be found in exactly two consecutive blocks, returns
    * @return true iff at least one stairlinkingvar  was assigned
    */
   bool assignCurrentStairlinking(
      );

   /**
    * @brief adds blocks and assigns open conss to such a new block or to master
    *  according to the cons assignment information given in constoblock hashmap
    *  @param constoblock hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
    *  @param additionalNBlocks number of (additional) blocks the hashmap contains
    *  @return scip return code
    *  \see assignSeeedFromConstoblockVector
    *  */
   SCIP_RETCODE assignSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int additionalNBlocks     /**< number of (additional) blocks the hashmap contains */
       );

   /*!
    * @brief adds blocks and assigns open conss to such a new block or to master
    *  according to the cons assignment information given in constoblock vector
    *  @param constoblock vector containing an assignment of conss to a block or to master
    *  @param additionalNBlocks number of (additional) blocks the vector contains
    *  @return scip return code
    *  \see  assignSeeedFromConstoblock()  */
   SCIP_RETCODE assignSeeedFromConstoblockVector(
      std::vector<int> constoblock, /**< vector containing an assignment of conss to a block or to master
                                      *< (master is indicated by assigning cons to index additionalNBlocks) */
      int additionalNBlocks        /**< number of (additional) blocks the vector contains */
      );

   /**
    * @brief books a constraint to be added to the block constraints of the given block (by calling flushBooked all bookings are in fact performed)
    * @param consToBlock constraint index to assign
    * @param block index of block cons is assigned to
    * @return scip return code
    *  \see flushBooked()
    */
   SCIP_RETCODE bookAsBlockCons(
      int consToBlock,
      int block
      );

   /**
    * @brief books a variable to be added to the block constraints of the given block (by calling flushBooked all bookings are in fact performed)
    * @param varToBlock variable index to be booked for block assignment
    * @param block index of block variables is assigned to
    * @return  scip return code
    * \see flushBooked()
    */
   SCIP_RETCODE bookAsBlockVar(
      int varToBlock,
      int block
      );

   /**
    * @brief  books a constraint to be added to the master constraints (by calling flushBooked all bookings are in fact performed)
    * @param consToMaster index of the constraint to be booked for master assignment
    * @return  scip return code
    * \see flushBooked()
    * */
   SCIP_RETCODE bookAsMasterCons(
         int consToMaster /*< this index can be computed by the function Seeedpool::getIndexForCons */
      );

   /**
    * @brief books a variable to be added to the master variables (by calling flushBooked all bookings are in fact performed)
    * @param varToMaster index index of the variable to be booked for master assignment
    * @return  scip return code
    * \see flushBooked()
    */
   SCIP_RETCODE bookAsMasterVar(
      int varToMaster
      );

   /**
    * @brief books a variable to be added to the linking variables (by calling flushBooked all bookings are in fact performed)
    * @param varToLinking index of variable that is booked for assigning to linking
    * @return scip return code
     * \see flushBooked()
    *  */
   SCIP_RETCODE bookAsLinkingVar(
      int varToLinking
      );

   /**
    * @brief books a variable to be added to the stairlinking variables of the given block and the following block (after calling
    *  flushBooked)
    * @param varToStairlinking index of variables to be assigned as stairlinking variable
    * @param firstBlock stairlinking variables hit exactly two consecutive blocks, this is the indwex of the first of these blocks
    * @return scip return code
     * \see flushBooked()
    *  */
   SCIP_RETCODE bookAsStairlinkingVar(
      int varToStairlinking,
      int firstBlock
      );

   /**
    * @brief checks if aggregation of sub problems is possible and stores the corresponding aggregation information
    */
   void calcAggregationInformation( );

   /**
    * @brief calculates the hash value of the seeed for comparing
    */
   void calcHashvalue();

   /**
    * @brief reassigns linking vars to stairlinkingvars if possible
    *  potentially reorders blocks for making a maximum number of linking vars stairlinking
    *  if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
    *  otherwise, the stairlinking assignment is done greedily
    *  precondition: seeed does not have any stairlinking vars
    */
   void calcStairlinkingVars(
        );

   /**
    * @brief counts for each pair of block and master constraint, how many nonzero entries the variables of the blocks have in the master constraint
    */
   void calcNCoeffsForBlockForMastercons(
        );


   /**
    *  @brief changes the block order in a way such that all linking vars that are potentially stairlinking
    *  may be reassigned to stairlinking
    * @param g graph with blocks as nodes and weighted edges for the number of
                         potentially stairlinkingvars connecting two blocks
    * @note precondition: all potentially stairlinking vars have a staircase structure */
   void changeBlockOrderStaircase(
        GraphGCG* g /**< graph with blocks as nodes and weighted edges for the number of
                         potentially stairlinkingvars connecting two blocks */
        );


   /**
    * @brief changes the block order in a way such that some linking vars that are potentially stairlinking
    *  may be reassigned to stairlinking using a greedy method
    *  \param g graph with blocks as nodes and weighted edges for the number of
    *                   potentially stairlinkingvars connecting two blocks
    */
   void changeBlockOrderGreedily(
      GraphGCG* g /**< graph with blocks as nodes and weighted edges for the number of
                       potentially stairlinkingvars connecting two blocks */
        );

   /**
    * @brief changes the order of the blocks according to the given mapping
    * \param oldToNewBlockIndex the mapping from old to new block indices
    * \note precondition: given mapping needs to be an adequately sized permutation */
   void changeBlockOrder(
        std::vector<int> oldToNewBlockIndex /**< the mapping from old to new block indices */
        );

   /**
    * @brief returns true iff all constraints are assigned and deletes the vector open conss if so
    * @return true iff all constraints are assigned
    * */
   bool checkAllConssAssigned();



   /**
    * @brief returns true if the assignments in the seeed are consistent
    * the following checks are performed:
    * 1) check if nblocks is set appropriately
    * 2) check for empty (row- and col-wise) blocks
    * 3) every variable is assigned at most once
    * 4) check if all not assigned variables are open vars
    * 5) check if all open vars are not assigned
    * 6) every constraint is assigned at most once
    * 7) check if all not assigned constraints are open cons
    * 8) check if all open conss are not assigned
    * 9) check if the datastructures are sorted
    * 10) check if variables hitting a cons are either in the cons's block or border or still open
    * @return true if the seeed seems to be consistent
    * */
   bool checkConsistency(
      );

   /**
    * @brief assigns all open constraints and open variables
    *  strategy: assigns all conss and vars to the same block if they are connected
    *  a cons and a var are adjacent if the var appears in the cons
    *  @return scip return code
    */
   SCIP_RETCODE completeByConnected(
      );


   /**
    * @brief computes components corresponding to connectedness of conss and vars as in @see completeByConnectedConssAdjacency
    * and assigns them accordingly but one of largest components
    * \see completeByConnected
    *  @return scip return code
    */
   SCIP_RETCODE assignSmallestComponentsButOneConssAdjacency(
        );


   /**
    * @brief try to reassign each mastercons to one block without inducing conflicts
    * @param success pointer to store whether at least one master constraint was reassigned
    * @return scip return code
    */
   SCIP_RETCODE postprocessMasterToBlocks(
        SCIP_Bool* success
      );


   /**
    * @brief try to reassign each mastercons to one block without inducing conflicts using the cons adjacency data structure
    * @param success pointer to store whether at least one master constraint was reassigned
    * @return scip return code
    */
   SCIP_RETCODE postprocessMasterToBlocksConssAdjacency(
        SCIP_Bool* success
      );

   /**
      * @brief assigns all open constraints and open variables
      *  strategy: assigns all conss and vars to the same block if they are connected
      *  a cons and a var are adjacent if the var appears in the cons
      *  \note   this relies on the consadjacency structure of the seeedpool
      *  hence it cannot be applied in presence of linking variables
      *  @return scip return code
      */
    SCIP_RETCODE completeByConnectedConssAdjacency(
         );

   /**
    * @brief assigns all open constraints and open variables
    *  strategy: assigns a cons (and related vars) to a new block if possible, if not to an existing block if possible (by means of prior var assignments)
    *  and finally to master, if there does not exist such a block
    *  @return scip return code
    */
   SCIP_RETCODE completeGreedily(
        );

   /**
    * @brief returns true if the given detector used a consclassifier
    * @param detectorchainindex index of the detector in the detectorchain
    * @return true iff the given detector used a consclassifier
    */
   bool consClassifierUsed(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /**
    * @brief: assigns every open cons/var in the following manner:
    *  - to the respective block if it hits exactly one blockvar/blockcons and no open vars/conss
    *  - to master/linking if it hits blockvars/blockconss assigned to different blocks
    *  - and every cons to master that hits a master var
    *  - and every var to master if it does not hit any blockcons and has no open cons
    *  - leave the cons/variableopen if nothing from the above holds
    *  @return scip return code
    *  */
   SCIP_RETCODE considerImplicits(
        );


   /**
    * @brief copies the given seeed's classifier statistics
    * @param otherseeed seeed whose classifier statistics are to be copied
    * @return scip return code
    */
   SCIP_RETCODE copyClassifierStatistics(
      const Seeed* otherseeed
      );

   /**
    * @brief deletes empty blocks and sets nblocks accordingly, a block is considered to be empty if no constraint is assigned to it, variables in blocks with no constraints become open
    * @param variables if true, then blocks with no constraints but at least one variable are considered to be nonempty
    * @return scip return code
    */
   SCIP_RETCODE deleteEmptyBlocks(
      bool variables
   );

   /**
    * @brief deletes a cons from list of open conss
    * @param opencons id of the cons that is not considered open anymore
    * @return scip return code
    */
   SCIP_RETCODE deleteOpencons(
      int opencons
      );

   /**
    * @brief deletes a var from the list of open vars
    * @param openvar id of the var that is not considered open anymore
    * @return scip return code
    */
   /** d */
   SCIP_RETCODE deleteOpenvar(
      int openvar
      );

   /**
    * @brief prints out the aggregation information that is calculated yet, i.e. if there has been identified identical blocks
    * @return scip return code
    */
   SCIP_RETCODE displayAggregationInformation();

   /**
    * @brief displays the assignments of the conss to blocks and master
    * @return scip return code
    */
   SCIP_RETCODE displayConss();


   /**
    * @brief displays the relevant information of the seeed
    * @param detailLevel pass a value that indicates how detailed the output should be:
    *                         0: brief overview
    *                         1: block and detector info
    *                         2: cons and var assignments
    * @return scip return code
    */
   SCIP_RETCODE displayInfo(
      int detailLevel
      );

   /**
    * @brief displays the relevant information of the seeed
    * @return scip return code
    */
   SCIP_RETCODE displaySeeed(
      );

   /**
    * @brief displays the assignments of the vars
    * @return scip return code
    */
   SCIP_RETCODE displayVars(
      );

   /**
    *@brief computes and returns the score of the given type of the seeed
    * @param type the scoretype that should be calculated
    * @return the score value (usually in [0,1] with 1 best poosible )
    * \see enum scoretype in cons_decomp.h for a list of scoretypes
    */
   /**  */
   SCIP_Real evaluate(
      SCORETYPE  type
      );

   /**
    * @brief every constraint is either assigned to master or open
    *  according to the cons assignment information given in constoblock hashmap
    *  variables are set accordingly
    * @note precondition: no constraint or variable is already assigned to a block
    * @param constoblock hashmap assigning cons indices (not SCIP_Cons* !!) to block indices (master assignment is indicated by assigning cons to index additionalNBlocks)
    * @param givenNBlocks number of blocks the hashmap contains
    * @return scip return code
    */
   SCIP_RETCODE filloutBorderFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks          /**< number of blocks the hashmap contains */
      );


   /**
    * @brief  assigns all conss to master or a block
    *  according to the cons assignment information given in constoblock hashmap
    * @param constoblock hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks)
    * @param givenNBlocks number of blocks the hashmap contains
    * @return scip return code
  *  calculates implicit variable assignment through cons assignment
    * @note precondition: no cons or var is already assigned to a block and constoblock contains information for every cons */

   SCIP_RETCODE filloutSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks          /**< number of blocks the hashmap contains */
      );


   /**
    * @brief reassigns variables classified as linking to master if the variable only hits master conss
    * @return scip return code
    */
   SCIP_RETCODE findVarsLinkingToMaster(
      );


   /**
    * @brief reassigns variables classified as linking to stairlinking if the variable hits conss in exactly two consecutive
    * blocks
    * @return scip return code
    */
   SCIP_RETCODE findVarsLinkingToStairlinking(
      );

   /**
    * @brief calculates potential stair linking variables with their blocks
    * @return a vector of pairs of var index and vector of (two) block indices
    *  the related linking variable hits exactly these two blocks given in the related vector
    */
   std::vector< std::pair< int, std::vector< int > > > findLinkingVarsPotentiallyStairlinking(
      );

   /**
    * @brief assigns all booked constraints and variables and deletes them from list of open cons and open vars
    * @return scip return code
    */
   SCIP_RETCODE flushBooked();


   /**
    * @brief gets seeed id of given ancestor id
    * @param ancestorindex index of ancestor seeed in ancestor list
    * @return seeed id of given ancestor id
    */
   int getAncestorID(
      int ancestorindex /**< index of ancestor in listofancestorids data structure */
      );



   /**
    * @brief get ancestor ids as vector
    * @return vector of ids of all ancestors id
    */
   std::vector<int> getAncestorList(
      );


   /**
    * set ancestor list directly
    * @param newlist new list of ancestor ids
    */
   void setAncestorList(
      std::vector<int> newlist
      );


   /**
    * adds ancestor id to back of list
    * @param ancestor id of ancestor that is to be added
    */
   void addAncestorID(
      int ancestor
      );


   /**
    * @brief get a vector of block ids that are identical to block with id repid
    * @param repid id of the representative block
    * @return vector of block ids that are identical to block with id repid
    */
   const std::vector<int> & getBlocksForRep(
      int  repid
      );


   /**
    * @brief the detectorchainstring contains the chars of all detectors that worked on this seeed in this order
    * @return detectorchainstring (containing the chars of all detectors that worked on this seeed in this order)
    */
   char* getDetectorChainString();


   /**
    * @brief returns detectorchain info of detector related to given detectorchain index
    * @param detectorchainindex index of the detector in the detectorchain
    * @return detectorchaininfo of the detector with the given index in the detector chain
    */
   std::string getDetectorchainInfo(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns the time that the detector related to the given detectorchainindex needed for detecting
    * @param detectorchainindex index of the detector the time that should be returned
    * @return the clock time for the corresponding detector in the chain
    */
   SCIP_Real getDetectorClockTime(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /**
    * @brief returns a vector of the clock times that each detector needed that was involved in this seeed
    * @return vector of the clock times
    */
   std::vector<SCIP_Real> getDetectorClockTimes();


   /**
    * @brief returns a string containing statistical data of the numbers of constraints and variables in the components:
    * in particular: ncomponents, percentage_min_nconss, percentage_max_nconss, percentage_median_nconss,
    * percentage_mean_nconss , percentage_min_nvars, percentage_max_nvars, percentage_median_nvars, percentage_mean_nvars
    * @return returns string with statistical data
    * @note used for features for miplib 2017
    */
   std::string getComponentInformation(
   );

   /**
    * @brief returns the data of the consclassifier that the given detector made use of
    * @param detectorchainindex index of the detector in the detectorchain
    * @param classifier a pointer to the used consclassifier (set by method)
    * @param consclassesmaster  a vector containing all indices of the consclasses assigned to master (set by method)
    * @return scip return code
    */
   SCIP_RETCODE getConsClassifierData(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      ConsClassifier** classifier, /**< a pointer to the used consclassifier */
      std::vector<int>& consclassesmaster /**< a vector containing all indices of the consclasses assigned to master */
      );

   /**
    * @brief returns array containing constraints assigned to a block
    * @param block id of the block the constraint indices are returned
    * @return array containing constraints assigned to a block
    */
   const int* getConssForBlock(
      int block
      );


   /**
    * @brief returns detector chain as array of detector pointers
    * @return detector chain as array of detector pointers
    */
   DEC_DETECTOR** getDetectorchain();


   /**
    * @brief returns the detectorchain as a vector of detector pointers
    * @return the detectorchain as a vector of detector pointers
    */
   std::vector<DEC_DETECTOR*> getDetectorchainVector();


   /**
    * @brief returns a string displaying all detector-related information, i.e. clock times and assignment data
    * @param detectorchainindex index of the detector in the detectorchain
    * @return string displaying all detector-related information, i.e. clock times and assignment data
    */
   /** returns a string displaying all detector-related information, i.e. clock times and assignment data */
   std::string getDetectorStatistics(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns a string displaying classifier information if such a classifier was used
    * @param detectorchainindex index of the detector in the detectorchain
    * @param displayConssVars pass true if constraints and variables of the respective classes should be displayed
    * @return string displaying classifier information if such a classifier was used
    */
   std::string getDetectorClassifierInfo(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      bool displayConssVars /**< pass true if constraints and variables of the respective classes should be displayed */
      );



   /**
    * @brief returns true iff this seeed was finished by finishSeeed() method of a detector
    * @return true iff this seeed was finished by finishSeeed() method of a detector
    */
   bool getFinishedByFinisher();


   /**
    * @brief returns true if the seeed is finished by a finisher in the unpresolved problem
    * @return true if the seeed is finished by a finisher in the unpresolved problem
    */
   bool getFinishedByFinisherUnpresolved();


   /**
    * @brief returns the detector that finished this seeed in the unpresolved problem if there exists one, NULL otherwise
    * @return the detector that finished this seeed in the unpresolved problem if there exists one, NULL otherwise
    * @note after finihed for unpresolved problem the transformed seeed (for the presolved problem) might be not completed
    */
   DEC_DETECTOR* getFinishedUnpresolvedBy();


   /**
    * @brief returns the calculated hash value of this seeed
    * @return the calculated hash value of this seeed
    */
   /** returns the calculated hash value of this seeed */
   long getHashValue();


   /**
    * @brief returns the unique id of the seeed
    * @return the unique id of the seeed
    */
   int getID();


   /**
    * @brief returns array containing all linking vars indices
    * @return array containing all linking vars indices
    * @note when accessed it is suppossed to be sorted
    */
   const int* getLinkingvars();


   /**
    * returns array containing all master conss indices
    * @return array containing all master conss indices
    * @note when accessed it is suppossed to be sorted
    */
   const int* getMasterconss();


   /**
    * returns array containing all master vars (hitting only constraints in the master, aka static variables) indices
    * @return array containing all master vars (hitting only constraints in the master, aka static variables) indices
    */
   const int* getMastervars();


   /**
    * @brief returns the "maximum white score"
    * @return  returns the "maximum white score"
    * @note "maximum white score" is fraction of the area of the decomposed matrix that is neither block or border
    */
   SCIP_Real getMaxWhiteScore();


   /**
    * @brief returns the experimental benders score
    * in detail:
    * bendersscore = max ( 0., 1 - ( 1 - blockareascore + (1 - borderareascore - bendersareascore ) ) ) with
    * blockareascore = blockarea / totalarea
    * borderareascore = borderarea / totalarea
    * bendersareascore = bendersarea /totalarea with
    * bendersarea = A + B - PENALTY with
    * A = nmasterconshittingonlyblockvars * nblockvarshittngNOmasterconss
    * B = nlinkingvarshittingonlyblockconss * nblockconsshittingonlyblockvars
    * PENALTY = \sum_{b=1}^(nblocks) \sum_{blockvars bv of block b hitting a master constraint} \sum_{all blocks b2 != b} nblockcons(b2)
    * @return experimental benders score
    */
   SCIP_Real getBendersScore();


   /**
    * @brief returns the number of nonzero coeffs in a certain block
    * @param blockid of the block the number of nozerors are requested for
    * @return number of nonzero coeffs in a certain block
    */
   int  getNCoeffsForBlock(
      int blockid
      );


   /**
    * returns the number of nonzero coeffs in master
    * @return the number of nonzero coeffs in master
    */
   int  getNCoeffsForMaster(
      );


   /**
    * @brief returns the score of the seeed (depending on used scoretype)
    * @param type the scoretype
    * @return the score
    * @see enum scoretype in cons_decomp.h
    */
   SCIP_Real getScore(
      SCORETYPE type
      );


   /**
    * @brief checks if all master constraints set partitioning, set packing, set cover, or cardinality constraints
    * @return TRUE iff all master constraints set partitioning, set packing, set cover, or cardinality constraints
    */
   SCIP_Bool hasSetppccardMaster(
   );


   /**
    * @brief checks iff all master constraints set partitioning, set packing, or set cover constraints
    * @return TRUE iff all master constraints set partitioning, set packing, or set cover
    */
   SCIP_Bool hasSetppcMaster(
   );


   /**
    * @brief checks iff all master constraints set partitioning, or set packing constraints
    * @return TRUE iff all master constraints set partitioning, or set packing constraints
    */
   SCIP_Bool hasSetppMaster(
   );


   /**
    * @brief returns the USERGIVEN status of this seeeds
    * @return the USERGIVEN status of this seeeds
    * @see enum USERGIVEN
    */
   USERGIVEN getUsergiven();


   /**
    * @brief returns number of ancestor seeeds
    * @return number of ancestor seeeds
    */
   int getNAncestors();


   /**
    * @brief returns the number of blocks
    * @return number of blocks
    */
   int getNBlocks();


   /**
    * @brief returns the number of constraints
    * @return number of constraints
    */
   int getNConss();


   /**
    * @brief returns size of the vector containing conss assigned to a block
    * @param block id of the block the number of constraints is asked for
    * @return size of the vector containing conss assigned to a block
    */
   int getNConssForBlock(
      int block
      );


   /**
    * @brief returns size of the detectorchain info vector
    * @return size of the detectorchain info vector
    */
   int getNDetectorchainInfo();

   /**
    * @brief returns the number of detectors the seeed is propagated by
    * @return  number of detectors the seeed is propagated by
    */
   int getNDetectors();

   /**
    * @brief returns the number used classifiers
    * @return number used classifiers
    */
   int getNUsedClassifier();

   /**
    * @brief returns size of the vector containing linking vars
    * @return size of the vector containing linking vars
    */
   int getNLinkingvars();


   /**
    * @brief returns size of the vector containing master conss
    * @return
    */
   int getNMasterconss();


   /**
    * @brief returns size of the vector containing master vars (hitting only constraints in the master)
    * @return size of the vector containing master vars (hitting only constraints in the master)
    */
   int getNMastervars();


   /**
    * @brief returns number of blocks a detector added
    * @param detectorchainindex index of the detector in the detectorchain
    * @return number of blocks a detector added
    */
   int getNNewBlocks(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief number of blocks the detectors in the detectorchain added
    * @return number of blocks the detectors in the detectorchain added
    */
   std::vector<int> getNNewBlocksVector();


   /**
    * @brief returns total number of stairlinking vars
    * @return total number of stairlinking vars
    */
   int getNTotalStairlinkingvars();


   /**
    * @brief returns size of vector containing constraints not assigned yet
    * @return returns size of vector containing constraints not assigned yet
    */
   int getNOpenconss();


   /**
    * @brief returns size of vector containing variables not assigned yet
    * @return size of vector containing variables not assigned yet
    */
   int getNOpenvars();


   /**
    * @brief returns the number of blockrepresentatives
    * @return the number of blockrepresentatives
    */
   int getNReps();



   /**
    * @brief returns size of the vector containing stairlinking vars
    * @param block id of the block the size of the stairlinking vector is asked for
    * @return size of the vector containing stairlinking vars
    */
   int getNStairlinkingvars(
      int block
      );


   /**
    * @brief returns number of vars
    * @return number of vars
    */
   int getNVars();



   /**
    * @brief returns size of the vector containing vars assigned to a block
    * @param block id of the block the number of variables is asked for
    * @return size of the vector containing vars assigned to a block
    */
   int getNVarsForBlock(
      int block
      );


   /**
    * @brief returns array containing constraints not assigned yet
    * @return array containing constraints not assigned yet
    */
   const int* getOpenconss();


   /**
    * @brief returns a vector containing constraint ids not assigned yet as vector
    * @return returns a vector containing constraint ids not assigned yet as vector
    */
   std::vector<int> getOpenconssVec();



   /**
    * @brief returns array containing variables not assigned yet
    * @return returns array containing variables not assigned yet
    */
   const int* getOpenvars();

   /**
    * returns array containing variables not assigned yet as vector
    * @return array containing variables not assigned yet as vector
    */
   std::vector<int> getOpenvarsVec();


   /**
    * @brief returns fraction of variables assigned to the border for a detector
    * @param detectorchainindex index of the detector in the detectorchain
    * @return fraction of variables assigned to the border for a detector
    */
   SCIP_Real getPctVarsToBorder(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );



   /**
    * @brief returns fraction of variables assigned to the border for detectors in detectorchain
    * @return vector of fractions of variables assigned to the border for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctVarsToBorderVector();


   /**
    * @brief returns fraction of variables assigned to a block for a detector
    * @param detectorchainindex index of the detector in the detectorchain
    * @return fraction of variables assigned to a block for a detector
    */
   SCIP_Real getPctVarsToBlock(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns fraction of variables assigned to a block for detectors in detectorchain
    * @return vector of fractions of variables assigned to a block for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctVarsToBlockVector();


   /**
    * @brief returns fraction of variables that are not longer open for a detector
    * @param detectorchainindex  index of the detector in the detectorchain
    * @return index of the detector in the detectorchain
    */
   SCIP_Real getPctVarsFromFree(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns fraction of variables that are not longer open for detectors in detectorchain
    * @return vecort or fractions of variables that are not longer open for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctVarsFromFreeVector();


   /**
    * @brief returns fraction of constraints assigned to the border for a detector
    * @param detectorchainindex index of the detector in the detectorchain
    * @return returns fraction of constraints assigned to the border for a detector
    */
   /**  */
   SCIP_Real getPctConssToBorder(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns fraction of constraints assigned to the border for detectors in detectorchain
    * @return vector of fractions of constraints assigned to the border for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctConssToBorderVector();


   /**
    * @brief returns fraction of constraints assigned to a block for a detector
    * @param detectorchainindex  index of the detector in the detectorchain
    * @return fraction of constraints assigned to a block for a detector
    */
   SCIP_Real getPctConssToBlock(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns fraction of constraints assigned to a block for detectors in detectorchain
    * @return vector of fractions of constraints assigned to a block for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctConssToBlockVector();


   /**
    * @brief returns fraction of constraints that are not longer open for a detector
    * @param detectorchainindex index of the detector in the detectorchain
    * @return fraction of constraints that are not longer open for a detector
    */
   SCIP_Real getPctConssFromFree(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief returns fraction of constraints that are not longer open for detectors in detectorchain
    * @return vector of fractions of constraints that are not longer open for detectors in detectorchain
    */
   std::vector<SCIP_Real> getPctConssFromFreeVector();


   /**
    * @brief  returns index of the representative block for a block, this might be blockid itself
    * @param blockid id of the block the representative is asked for
    * @return index of the representative block for a block, this might be blockid itself
    */
   int getRepForBlock(
      int blockid
      );

   /**
    * @brief returns the represenation varmap as vector for represenative repid and the blockrepid-th block that is represented by repid
    * @param repid
    * @param blockrepid
    * @return the represenation varmap as vector for represenative repid and the blockrepid-th block that is represented by repid
    */

   std::vector<int> & getRepVarmap(
      int repid,
      int blockrepid
      );


   /**
    * @brief returns the corresponding seeedpool
    * @return corresponding seeedpool
    */
   Seeedpool* getSeeedpool();


   /**
    * @brief returns array containing stairlinking vars,
    * @note if a stairlinking variable links block i and i+1 it is only stored in vector of block i
    * @param block id of the block the stairlinking variable varctor is asked for
    * @return array containing stairlinking vars,
    */
   const int* getStairlinkingvars(
      int block
      );


   /**
    * @brief returns true if this seeed stems from the unpresolved problem
    * @return  TRUE iff seeed stems from the unpresolved problem
    */
   bool getStemsFromUnpresolved();


   /**
    * @brief returns the data of the varclassifier that the given detector made use of
    * @param detectorchainindex index of the detector in the detectorchain
    * @param classifier a pointer to the used varclassifier, , to be set by method
    * @param varclasseslinking a vector containing all indices of the varclasses assigned to linking, to be filled by method
    * @param varclassesmaster a vector containing all indices of the varclasses assigned to master, to be filled by method
    * @return data of the varclassifier that the given detector made use of
    */
   SCIP_RETCODE getVarClassifierData(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      VarClassifier** classifier, /**< a pointer to the used varclassifier */
      std::vector<int>& varclasseslinking, /**< a vector containing all indices of the varclasses assigned to linking */
      std::vector<int>& varclassesmaster /**< a vector containing all indices of the varclasses assigned to master */
      );


   /**
    * @brief returns array containing vars of a block
    * @param block id of the block the vars are requested for
    * @return returns array containing vars of a block
    */
   const int* getVarsForBlock(
      int block
      );


   /**
    * @brief  returns index in variables array of a block for a variable
    * @param varid the id of the variable the index
    * @param block the corresponding block id
    * @return  returns index in variables array of a block for a variable
    */
   int getVarProbindexForBlock(
      int varid,
      int block
   );


   /**
    * @brief prepare the seeed such that all predecessors have the folloing property:
    * all variables in the master problem are binary variables
    * thus all other variables are assigned to a block
    * requirement: all constraints and variables are open when this method is called
    */
   void initOnlyBinMaster();

   /**
    * @brief checks if calculation of aggregation information is considered to be to expansive
    * @return TRUE iff calculation of aggregation information is considered to be to expansive
    */
   SCIP_Bool isAgginfoToExpensive();



   /**
    * @brief   returns true if this seeed is complete,
    *  i.e. it has no more open constraints and variables
    * @return TRUE iff this seeed is complete
    */
   bool isComplete();


   /**
    * @brief returns true if the cons is a cons of the block
    * @param cons id of constraint to check
    * @param block id of the blocks
    * @return true iff the cons is a cons of the block
    */
   /**  */
   bool isConsBlockconsOfBlock(
      int cons,
      int block
      );

   /**
    * @brief returns true if the cons is a master cons
    * @param cons id of ccons to check if it is master constraint
    * @return true iff the cons is a master cons
    */
   bool isConsMastercons(
      int cons
      );


   /**
    * @brief returns true if the cons is an open cons
    * @param cons id of cons to check
    * @return true iff the cons is an open cons
    */
   bool isConsOpencons(
      int cons
      );


   /**
    * @brief returns true if the seeed is from a detector operating in legacymode
    * @return retruns true iff this partial decomposition is found during legacy mode
    */
   bool isFromLegacymode();


   /**
    * @brief returns true if the seeed is from the unpresolved problem
    * @return true iff the seeed is from the unpresolved problem
    */
   bool isFromUnpresolved();


   /**
    * returns true if the seeed is currently selected in explore menue
    * @return  true iff the seeed is currently selected in explore menue
    */
   bool isSelected();


   /**
    * @brief method to check whether this seeed is equal to a given other seeed ( \see  isEqual(Seeed*))
    * @param otherseeed seeed to check euality with
    * @param isequal pointer to store whether seeeds are identical
    * @param sortseeeds should conss and vars be sorted before comparing the seeeds?
    * @return
    */
   SCIP_RETCODE isEqual(
      Seeed* otherseeed,   /**< other seeed */
      SCIP_Bool* isequal,  /**< pointer to store whether seeeds are identical */
      bool sortseeeds      /**< should conss and vars be sorted before comparing the seeeds? */
      );

   /**
    * @brief method to check whether this seeed is equal to a given other seeed
    * @param other seed to check equality with
    * @return true iff seeeds are equal
    */
   bool isEqual(
      Seeed* other /**< other seeed */
      );

   /**
    * @brief returns true if this seeed was propagated by specified detector
    * @param detectorID pointer to detector to check for
    * @return true iff this seeed was propagated by  detectorID
    */
   bool isPropagatedBy(
      DEC_DETECTOR* detectorID
      );


   /**
    * @brief returns true if this seeed is considered to be trivial,
    *  i.e. all conss are in one block, all conss are in border, all variables linking or mastervars, or all constraints and variables are open
    * @return true iff this seeed is considered to be trivial
    */
   bool isTrivial();


   /**
    * @brief returns true if the var is assigned to the block
    * @param var id of var to check
    * @param block id of block to check
    * @return true iff the var is assigned to the block
    */
   bool isVarBlockvarOfBlock(
      int var,
      int block
      );


   /**
    * @brief returns true if the var is a linking var
    * @param var id of var to check
    * @return true iff the var is a linking var
    */
   bool isVarLinkingvar(
      int var
      );


   /**
    * @brief  returns true if the var is a master var
    * @param var id of var to check
    * @return  true iff the var is a master var
    */
   bool isVarMastervar(
      int var
      );



   /**
    * @brief returns true if the var is an open var
    * @param var id of var to check
    * @return  true iff the var is an open var
    */
   /**  */
   bool isVarOpenvar(
      int var
      );


   /**
    * @brief returns true if the var is a stairlinking var
    * @param var id of var to check
    * @return true if the var is a stairlinking var
    */
   bool isVarStairlinkingvar(
      int var
      );


   /**
    * @brief returns true if the var is a stairlinkingvar of a speciefied block
    * @param var id of var to check if it is a stairlinking variable hitting specified block
    * @param block id of block to check
    * @return true iff the var is a stairlinkingvar of a speciefied block
    */
   bool isVarStairlinkingvarOfBlock(
      int var,
      int block
      );


   /**
    * @brief prints classifier information as described in \see cls reader
    * @param scip scip data structure
    * @param file output file
    * @return scip return code
    */
   SCIP_RETCODE printClassifierInformation(
      SCIP*                scip,
      FILE*                file);



   /**
    * @brief refine seeed with focus on blocks: assigns open conss and vars if they can be found in blocks (without respect to open vars and conss  @see assignHittingOpenconss(), @see assignHittingOpenvars())
    * @note seeed is might be not complete
    * @return scip return code
    */
   SCIP_RETCODE refineToBlocks(
      );

   /**
    * @brief refine seeed with focus on master: do obvious ( @see considerImplicits()) assignments and
    *  assign other conss and vars to master if possible (@see assignOpenPartialHittingToMaster())
    * @return scip return code
    */
   /**  */
   SCIP_RETCODE refineToMaster(
      );


   /**
    * @brief registers statistics for a used consclassifier
    * @param detectorchainindex index of the detector in the detectorchain
    * @param classifier the used consclassifier
    * @param consclassesmaster vector of classindices that were assigned to master
    */
   void setConsClassifierStatistics(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      ConsClassifier* classifier, /**< the used consclassifier */
      std::vector<int> consclassesmaster /**< vector of classindices that were assigned to master */
      );


   /**
    * @brief adds a constraint to a block, does not delete this cons from list of open conss
    * @param consToBlock id of cons to add
    * @param block id of block to add
    * @return scip return code
    */
   SCIP_RETCODE setConsToBlock(
      int consToBlock,
      int block
      );


   /**
    * @brief adds a constraint to the master constraints, does not delete this cons from list of open conss
    * @param consToMaster id of cons to add
    * @return scip return code
    */
   SCIP_RETCODE setConsToMaster(
      int consToMaster
      );


   /**
    *  sets the detectorchain with the given vetcor of detector pointers
    * @param detectorChain vetcor of detector pointers
    */
   void setDetectorchain(
      std::vector<DEC_DETECTOR*> detectorChain
      );


   /**
    * @brief sets seeed to be propagated by a detector
    * @param detector pointer to detector that is registered for this seeed
    * @return scip return code
    */
   SCIP_RETCODE setDetectorPropagated(
      DEC_DETECTOR* detector
      );


   /**
    * @brief sets seeed to be finished by a detector
    * @param detector pointer to detector that has finished this seeeds
    * @return scip return code
    */
   SCIP_RETCODE setFinishingDetectorPropagated(
      DEC_DETECTOR* detector
      );


   /**
    * @brief sets whether this seeed was finished by a finishing detector
    * @param finished is this seeds finished by a finishing detector
    */
   void setFinishedByFinisher(
      bool finished
      );


   /**
    * @brief sets whether this seeed is finished by a finisher in the unpresolved problem
    * @param finishedByFinisherUnpresolved is this seeed finished by a finisher in the unpresolved problem
    *
    */
   void setFinishedByFinisherUnpresolved(
      bool finishedByFinisherUnpresolved
      );

   /**
    * @brief  sets the detector that finished the seeed in the unpresolved problem
    * @param detector pointe of detector that has finished this seeed in unpresolved problem
    */
   void setFinishedUnpresolvedBy(
      DEC_DETECTOR* detector
      );


   /**
    * @brief sets whether this seeed stems from a detector operating in legacymode
    * @param legacymode true iff this seeed stems from legacy mode detection
    */
   void setLegacymode(
      bool legacymode
      );


   /**
    * @brief sets number of blocks, only increasing number allowed
    * @param nBlocks new number of blocks
    * @return scip return code
    */
   SCIP_RETCODE setNBlocks(
      int nBlocks
      );

   /**
    * @brief sets the id of the seeed
    * @param id id to be set
    * @return scip return code
    */
   SCIP_RETCODE setID(
      int id
      );


   /**
    * @brief sets whether this seeed is from the unpresolved problem
    * @param unpresolved is the seeed from unpresolved problem
    */
   void setIsFromUnpresolved(
      bool unpresolved
      );


   /**
    * @brief set if selection status of this seeeds
    * @param selected
    */
   /** sets whether this seeed is selected */
   void setSelected(
      bool selected
      );


   /**
    * @brief set the corresponding seeedpool of this seeeds
    * @param seeedpool pointer seeedpool to be set
    */
   void setSeeedpool(
      Seeedpool* seeedpool
      );


   /**
    * @brief sets whether this seeed stems from an unpresolved problem seeed
    * @param stemsfromunpresolved has this seeed ancestors from the unpresolved probelm
    */
   void setStemsFromUnpresolved(
      bool stemsfromunpresolved
      );


   /**
    * @brief sets whether this seeed is usergiven
    * @param usergiven is this seeed user given
    */
   void setUsergiven(
      USERGIVEN usergiven
      );


   /**
    * @brief registers statistics for a used varclassifier
    * @param detectorchainindex index of the detector in the detectorchain
    * @param classifier the used varclassifier
    * @param varclasseslinking vector of classindices that were assigned to linking
    * @param varclassesmaster  vector of classindices that were assigned to master
    */
   /** registers statistics for a used varclassifier */
   void setVarClassifierStatistics(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      VarClassifier* classifier, /**< the used varclassifier */
      std::vector<int> varclasseslinking, /**< vector of classindices that were assigned to linking */
      std::vector<int> varclassesmaster /**< vector of classindices that were assigned to master */
      );


   /**
    * @brief adds a variable to the linking variables,  does not delete this var from list of open vars
    * @param varToBlock id of var to be added
    * @param block id of block to be added
    * @return scip return code
    */
   SCIP_RETCODE setVarToBlock(
      int varToBlock,
      int block
      );


   /**
    * @brief adds a variable to the linking variables, does not delete this var from list of open vars
    * @param varToLinking if of var to be set to linking
    * @return scip return code
    */
   SCIP_RETCODE setVarToLinking(
      int varToLinking
      );



   /** directly adds a variable to the master variables (hitting only constraints in the master)
    *  does not delete this var from list of open vars */
   SCIP_RETCODE setVarToMaster(
      int varToMaster
      );

   /**
    * @brief adds a variable to the stairlinking variabl, does not delete this var from list of open vars
    * @param varToStairLinking id of variable to be added
    * @param block1 id of block one
    * @param block2 id of block two
    * @note stairlinking variables are only registered in block with smalller index
    * @return scip return code
    */
   SCIP_RETCODE setVarToStairlinking(
      int varToStairLinking,
      int block1,
      int block2
      );


   /**
    * @brief generates and opens a gp visualization of the seeed
    * @see visual/pdfreader and
    */
   void showVisualisation();


   /**
    * @brief returns true if this seeed is a userseeed that should be completed by setting unspecified constraints to master
    * @return TRUE iff this seeed is a userseeed that should be completed by setting unspecified constraints to master
    */
   SCIP_Bool shouldCompletedByConsToMaster();


   /**
    * @brief sorts the vars and conss datat structures  by their indices
    */
   void sort();


   /**
    * @brief returns a short caption for this seeed
    * @return short caption for this seeed
    */
   const char* getShortCaption();



   /**
    * @brief sets the detector chain short string
    * @param detectorchainstring the detector chain string to set
    * @return scip return code
    */
   SCIP_RETCODE setDetectorChainString(
      char* detectorchainstring
      );


   /**
    * @brief set statistical vector of numbers of newly assigned blocks per involved detector
    * @param newvector vector of numbers of newly assigned blocks per involved detector
    */
   void setNNewBlocksVector(
      std::vector<int>  newvector
);


   /**
    * @brief set statistical vector of fractions of constraints set to blocks per involved detector
    * @param newvector vector of fractions of constraints set to blocks per involved detector
    */
   void setPctConssToBlockVector(
      std::vector<SCIP_Real> newvector
      );


   /**
    * @brief set statistical vector of fractions of constraints that are not longer open  per involved detector
    * @param newvector vector of fractions of constraints that are not longer open  per involved detector
    */
   void setPctConssFromFreeVector(
      std::vector<SCIP_Real> newvector
   );

   /**
    * @brief set statistical vector of fractions of constraints assigned to the border per involved detector
    * @param newvector vector of fractions of constraints assigned to the border per involved detector
    */
   void setPctConssToBorderVector(
      std::vector<SCIP_Real> newvector
      );


   /**
    * @brief set statistical vector of fraction sof variables assigned to the border per involved detector
    * @param newvector vector of fractions of variables assigned to the border per involved detector
    */
   void setPctVarsToBorderVector(
      std::vector<SCIP_Real> newvector
      );


   /**
    * @brief set statistical vector of fractions of variables assigned to a block per involved detector
    * @param newvector vector of fractions of variables assigned to a block per involved detector
    */
   void setPctVarsToBlockVector(
      std::vector<SCIP_Real> newvector
   );



   /**
    * @brief set statistical vector of variables that are not longer open per involved detector
    * @param newvector vector of fractions of variables that are not longer open per involved detector
    */
   void setPctVarsFromFreeVector(
      std::vector<SCIP_Real> newvector
      );

   /**
    * @brief set statistical vector of the times that the detectors needed for detecting per involved detector
    * @param newvector vector of the times that the detectors needed for detecting per involved detector
    */
   void setDetectorClockTimes(
      std::vector<SCIP_Real> newvector
      );


   /**
    * @brief returns true if the given detector used a varclassifier
    * @param detectorchainindex index of the detector in the detectorchain
    * @return true if the given detector used a varclassifier
    */
   bool varClassifierUsed(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );


   /**
    * @brief write this seeed to file in dec format for the corresponding seeedpool
    * @param file pointer to file to write to
    * @param seeedpooltowriteto if this is not the seeedpool of the seeed, the seeed is tried to be transformed
    * @param result will be set to SCIP_SUCCESS if writeing (including possible transformation) was successful
    * @return scip return code
    */
   SCIP_RETCODE writeAsDec(
      FILE* file,
      Seeedpool*   seeedpooltowriteto,
      SCIP_RESULT* result
      );



   /**
    * @brief creates and sets a detector chain short string for this seeed, is built from detector chain
    * @return scip return code
    */
   /** creates and sets a detector chain short string for this seeed */
   SCIP_RETCODE buildDecChainString();

private:


   /**
    * @brief  adds empty entries for all classifier statistics for a detector added to the detector chain
    */
   void addEmptyClassifierStatistics();


   /** assigns open cons
    *  - to master if it hits blockvars of different blocks
    *  - to the respective block if it hits a blockvar of exactly one block and no stairlinking var
    *  - to master if it hits a stairlinking var but there is no block the cons may be assigned to
    *  - to the block with the lowest number of conss if it hits a stairlinking var and there are blocks the cons may be
    *    assigned to
    *  - leave it open if it hits no blocks yet
    *  @return true iff some assignment was made by the method
    */
   bool assignHittingOpenconss(
      );



   /** @brief assigns every open var
    *  - to the respective block if it hits blockconss of exactly one block
    *  - to linking if it hits blockconss of more than one different blocks
    *  - leave the var open otherwise
    *  @return true iff there is a var that has been assigned in this call*/
   bool assignHittingOpenvars(
      );


   /**
    * @brief assigns every open cons to master that hits
    *  - exactly one block var and at least one open var or
    *  - a master var
    *  - or leave it open elsewise
    *  @return scip return code
    */
   SCIP_RETCODE assignOpenPartialHittingConsToMaster(
      );



   /**
    * @brief assigns open conss/vars that hit exactly one block and at least one open var/cons to border
    * @return scip return code
    */
   SCIP_RETCODE assignOpenPartialHittingToMaster(
      );



   /**
    * @brief assigns every open var to linking that hits
    *  - exactly one block cons and at least one open cons
    *  - leave it open otherwise
    *  @return scip return code
    */
   SCIP_RETCODE assignOpenPartialHittingVarsToMaster(
      );


   /**
    * @brief calculates the number of nonzero coefficients for the blocks
    * @return scip return code
    */
   SCIP_RETCODE calcNCoeffsForBlocks(
   );



   /**
    * @brief calc maximum white score for this seeed
    * @note "maximum white score" is fraction of the area of the decomposed matrix that is neither block or border
    */
   void calcmaxwhitescore();


   /**
    * @brief calculates the experimental benders score
    * in detail:
    * bendersscore = max ( 0., 1 - ( 1 - blockareascore + (1 - borderareascore - bendersareascore ) ) ) with
    * blockareascore = blockarea / totalarea
    * borderareascore = borderarea / totalarea
    * bendersareascore = bendersarea /totalarea with
    * bendersarea = A + B - PENALTY with
    * A = nmasterconshittingonlyblockvars * nblockvarshittngNOmasterconss
    * B = nlinkingvarshittingonlyblockconss * nblockconsshittingonlyblockvars
    * PENALTY = \sum_{b=1}^(nblocks) \sum_{blockvars bv of block b hitting a master constraint} \sum_{all blocks b2 != b} nblockcons(b2)
    */
   void calcbendersscore();

   /**
    * @brief calculates classical score
    * classical score = (0.6 * ( borderscore ) + 0.2 * ( linkingscore ) + 0.2 * ( densityscore ) ) with
    * borderarescore = 1 - (borderarea / matrixarea)
    * densityscore = min_{b \in blocks} blockdensity(b) , blockdensity = nonzeros(b)/(nvars(b)*nconss(b))
    * linkingscore = 1 - ( 0.5 + 0.5 * varratio ) with
    * varratio = \prod_{b \in blocks} nlinkingvarshittingblock(b) / (ntotalinkingvars + mastervars)
    * @return scip return code
    */
   SCIP_RETCODE calcclassicscore();

   /**
    * @brief calculates borderareascore
    * borderareascore = = 1 - (borderarea / matrixarea)
    */
   void calcborderareascore();

   /**
    * @brief calculates max foreseeeing white score
    * \see calcmaxwhitescore() for matrix with this adaptions: (to apply DW decomposition)
    * linking vars are copied for each pricing prob they occur and new constraints assuring equality of the copies are added to master problem
    */
   void calcmaxforeseeingwhitescore();


   /**
    * @brief calculates max foreseeing white score aggregated
    * \see calcmaxforeseeingwhitescore(), but identical blocks are only considered once for block area
    */
   void calcmaxforeseeingwhitescoreagg();

   /**
    * @brief calculates setpartfwhitescore
    * setpartfwhitescore = 0.5 * setpartindicator + 0.5
    * setpartindicator = 1 if master problem contains only setppc or cardianlity cosntraints = 0, otherwise
    * \see calcmaxforeseeingwhitescore()
    * \see hasSetppccardMaster()
    */
   void calcsetpartfwhitescore();


   /**
    * @brief calculates setpartfwhitescoreagg
    * \same as \see calcsetpartfwhitescore but identical blocks are only considered once for block area
    */
   void calcsetpartfwhitescoreagg();

   /**
    * @brief calculates the bendersarea score
    * \see calcbendersscore()
    * in detail:
    * bendersareascore = bendersarea /totalarea with
    * bendersarea = A + B - PENALTY with
    * A = nmasterconshittingonlyblockvars * nblockvarshittngNOmasterconss
    * B = nlinkingvarshittingonlyblockconss * nblockconsshittingonlyblockvars
    * PENALTY = \sum_{b=1}^(nblocks) \sum_{blockvars bv of block b hitting a master constraint} \sum_{all blocks b2 != b} nblockcons(b2)
    */
   void calcbenderareascore();

   /**
    * @brief calculates block area score
    * blockareascore = (blockarea/totalarea)
    */
   void calcblockareascore();


   /**
    * @brief calculates block area score aggregated
    * blockareascoreagg is the same as @see calcblockareascore but identical blocks are only considered once for blockarea
    */
   void calcblockareascoreagg();

};

} /* namespace gcg */
#endif /* GCG_CLASS_Seeed_H__ */
