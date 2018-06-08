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


enum USERGIVEN
{
   NOT = 0,
   PARTIAL = - 1,
   COMPLETE = - 2,
   COMPLETED_CONSTOMASTER = - 3
};

class Seeedpool;

class Seeed
{
private:
   SCIP* scip;                                                 /**< SCIP data structure */
   int id;                                                     /**< id of the seeed */
   int nBlocks;                                                /**< number of blocks the decomposition currently has */
   int nVars;                                                  /**< number of variables */
   int nConss;                                                 /**< number of constraints */
   std::vector<int> masterConss;                               /**< vector containing indices of master constraints */
   std::vector<int> masterVars;                                /**< vector containing indices of master variables */
   std::vector<std::vector<int>> conssForBlocks;               /**< conssForBlocks[k] contains a vector of indices of all
                                                                 *< constraints assigned to block k */
   std::vector<std::vector<int>> varsForBlocks;                /**< varsForBlocks[k] contains a vector of indices of all
                                                                 *< variables assigned to block k */
   std::vector<int> linkingVars;                               /**< vector containing indices of linking variables */
   std::vector<std::vector<int>> stairlinkingVars;             /**< vector containing indices of staircase linking variables
                                                                 *< of the blocks (stairlinking variables can be found only
                                                                 *< in the vector of their first block) */
   std::vector<int> openVars;                                  /**< vector containing indices of variables that are not
                                                                 *< assigned yet*/
   std::vector<int> openConss;                                 /**< vector containing indices of constraints that are not
                                                                 *< assigned yet*/
   std::vector<bool> isvaropen;
   std::vector<bool> isconsopen;
   std::vector<bool> isvarmaster;
   std::vector<bool> isconsmaster;

   std::vector<int>  ncoeffsforblock;                          /**< number of coeffs per block */

   SCIP_Bool         calculatedncoeffsforblock;                        /**< is the  number of coeff per block already calculated*/
   int               ncoeffsformaster;                                    /**< number of master coefficients */

   bool varsforblocksorted;
   bool stairlinkingvarsforblocksorted;
   bool conssforblocksorted;
   bool linkingvarssorted;
   bool mastervarssorted;
   bool masterconsssorted;


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
   long hashvalue;
   bool changedHashvalue;                                      /**< are there any changes concerning the hash value since it
                                                                 *< was calculated last time */

   bool isselected;                                            /**< is this seeed selected */

   bool isagginfoalreadytoexpensive;                            /** is agginfo already known to be to expensive */

   const static int primes[];
   const static int nPrimes;

   bool isFinishedByFinisher;                         /**< was this seeed finished by the finishseeed() method of a detector */

   std::vector<std::vector<int>>    ncoeffsforblockformastercons;    /**< number of coeffs a block has in a certain master constraint */

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

   std::vector<int> listofancestorids;                /**< vector containing detector indices that worked on that seeed */
   USERGIVEN usergiven;                               /**< is this seeed partially or completely given by user */
   bool isfromlegacymode;                             /**< true if this seeed stems from a detector operating in legacy mode */
   SCIP_Real score;                                   /**< score to evaluate the seeeds */
   SCIP_Real maxwhitescore;                           /**< score corresponding to the max white measure */
   SCIP_Real bendersscore;                           /**< score to evaluate the seeeds */
   SCIP_Real benderareascore;                          /**< 1 - fraction of white area iin master constraints to complete area */

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

   Seeedpool*    seeedpool;               /**< seeedpool for the corresponding problem */

   /** checks blocks for identity by graph automorphism check done by bliss, identity is only found if variables are in correct order */
   void checkIdenticalBlocksBliss(
      Seeedpool*           seeedpool,
      int                  b1,
      int                  b2,
      std::vector<int>&    varmap,         /**< maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1 */
      SCIP_HASHMAP*        varmap2,
      SCIP_Bool*           identical
      );


   /** checks blocks for identity by brute force, identity is only found if variables are in correct order */
   void checkIdenticalBlocksBrute(
      Seeedpool*           seeedpool,
      int                  b1,
      int                  b2,
      std::vector<int>&    varmap,         /**< maps variable indices (corresponding to  seeedpool indices) of prob2 to prob1 */
      SCIP_HASHMAP*        varmap2,
      SCIP_Bool*           identical
      );

   SCIP_RETCODE checkIdenticalBlocksTrivial(
      Seeedpool*           givenseeedpool,
      int                  b1,
      int                  b2,
      SCIP_Bool*           notidentical
      );

public:

   /** constructor
    *  initially, all conss and vars are open */
   Seeed(
      SCIP* scip,       /**< scip data structure */
      int id,           /**< id that is given to this seeed */
      Seeedpool* seeedpool
      );

   /** copy constructor */
   Seeed(
      const Seeed *seeedToCopy /**< seeed to be copied */
      );

   /** destructor */
   ~Seeed();

   SCIP_Bool isconshittingblockca(
      gcg::Seeedpool* seeedpool,
      int masterconsid,
      int b
      );

   /** adds a block, returns the number of the new block */
   int addBlock();

   /** incorporates the needed time of a certain detector in the detector chain */
   void addClockTime(
      SCIP_Real clocktime /**< time to be added */
      );

   /** incorporates the changes from ancestor seeed */
   void addDecChangesFromAncestor(
      Seeed* ancestor
      );

   /** adds a detector chain info */
   void addDetectorChainInfo(
      const char* decinfo
      );

   /** adds number of new blocks created by a detector added to detector chain */
   void addNNewBlocks(
      int nnewblocks
      );

   /** adds fraction of constraints that are not longer open for a detector added to detector chain */
   void addPctConssFromFree(
      SCIP_Real pct
      );

   /** adds fraction of constraints assigned to a block for a detector added to detector chain */
   void addPctConssToBlock(
      SCIP_Real pct
      );

   /** adds fraction of constraints assigned to the border for a detector added to detector chain */
   void addPctConssToBorder(
      SCIP_Real pct
      );

   /** adds fraction of variables that are not longer open for a detector added to detector chain */
   void addPctVarsFromFree(
      SCIP_Real pct
      );

   /** adds fraction of variables assigned to a block for a detector added to detector chain */
   void addPctVarsToBlock(
      SCIP_Real pct
      );

   /** adds fraction of variables assigned to the border for a detector added to detector chain */
   void addPctVarsToBorder(
      SCIP_Real pct
      );

   /** returns true if at least one constraint is assigned to a block */
   bool alreadyAssignedConssToBlocks();

   /** assigns open conss to master according to the cons assignment information given in constoblock hashmap */
   SCIP_RETCODE assignBorderFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks,          /**< number of blocks the hashmap contains */
      Seeedpool* seeedpool       /**< a seeedpool that uses this seeed */
      );

   /** assigns open vars to stairlinking if they can be found in two consecutive blocks, returns true if stairlinkingvars
    *  are assigned */
   bool assignCurrentStairlinking(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** adds blocks and assigns open conss to such a new block or to master
    *  according to the cons assignment information given in constoblock hashmap */
   SCIP_RETCODE assignSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int additionalNBlocks,     /**< number of (additional) blocks the hashmap contains */
      Seeedpool* seeedpool       /**< a seeedpool that uses this seeed */
      );

   /** adds blocks and assigns open conss to such a new block or to master
    *  according to the cons assignment information given in constoblock vector */
   SCIP_RETCODE assignSeeedFromConstoblockVector(
      std::vector<int> constoblock, /**< vector containing an assignment of conss to a block or to master
                                      *< (master is indicated by assigning cons to index additionalNBlocks) */
      int additionalNBlocks,        /**< number of (additional) blocks the vector contains */
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
      );

   /** books a constraint to be added to the block constraints of the given block (after calling flushBooked) */
   SCIP_RETCODE bookAsBlockCons(
      int consToBlock,
      int block
      );

   /** books a variable to be added to the block constraints of the given block (after calling flushBooked) */
   SCIP_RETCODE bookAsBlockVar(
      int varToBlock,
      int block
      );

   /** books a constraint to be added to the master constraints (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterCons(
         int consToMaster /*< this index can be computed by the function Seeedpool::getIndexForCons */
      );

   /** books a variable to be added to the master variables (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterVar(
      int varToMaster
      );

   /** books a variable to be added to the linking variables (after calling flushBooked) */
   SCIP_RETCODE bookAsLinkingVar(
      int varToLinking
      );

   /** books a variable to be added to the stairlinking variables of the given block and the following block (after calling
    *  flushBooked) */
   SCIP_RETCODE bookAsStairlinkingVar(
      int varToStairlinking,
      int firstBlock
      );

   /** checks if aggregation of sub problems is possible and stores the corresponding aggreagtion information */
   void calcAggregationInformation(
      Seeedpool*  seeedpool
      );

   /** calculates the hash value of the seeed for comparing */
   void calcHashvalue();

   /** reassigns linking vars stairlinkingvars if possible
    *  potentially reorders blocks for making a maximum number of linking vars stairlinking
    *  if all vars that connect exactly two blocks have a staircase structure, all of them become stairlinkingvars
    *  otherwise, the stairlinking assignment is done greedily
    *  precondition: seeed does not have any stairlinking vars */
   void calcStairlinkingVars(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );


   void calcNCoeffsForBlockForMastercons(
      Seeedpool*           givenseeedpool
      );


   /** changes the block order in a way such that all linking vars that are potentially stairlinking
    *  may be reassigned to stairlinking
    *  precondition: all potentially stairlinking vars have a staircase structure */
   void changeBlockOrderStaircase(
        GraphGCG* g /**< graph with blocks as nodes and weighted edges for the number of
                         potentially stairlinkingvars connecting two blocks */
        );

   /** changes the block order in a way such that some linking vars that are potentially stairlinking
    *  may be reassigned to stairlinking using a greedy method */
   void changeBlockOrderGreedily(
      GraphGCG* g /**< graph with blocks as nodes and weighted edges for the number of
                       potentially stairlinkingvars connecting two blocks */
        );

   /** changes the order of the blocks according to the given mapping
    *  precondition: given mapping needs to be an adequately sized permutation */
   void changeBlockOrder(
        std::vector<int> oldToNewBlockIndex /**< the mapping from old to new block indices */
        );

   /** returns true if all constraints are assigned and deletes the vector open conss if so */
   bool checkAllConssAssigned();

   /** returns true if the assignments in the seeed are consistent */
   bool checkConsistency(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** assigns all open constraints and open variables
     *  strategy: assigns all conss and vars to the same block if they are connected
     *  a cons and a var are adjacent if the var appears in the cons */
   SCIP_RETCODE completeByConnected(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** assigns all open constraints and open variables
     *  strategy: assigns all conss and vars to the same block if they are connected
     *  a cons and a var are adjacent if the var appears in the cons */
   SCIP_RETCODE assignSmallestComponentsButOneConssAdjacency(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );


   /** try to reassign each  mastercons to one block without inducing conflicts  */
   SCIP_RETCODE postprocessMasterToBlocks(
      Seeedpool* seeedpool, /**< a seeedpool that uses this seeed */
      SCIP_Bool* success
      );


   /** try to reassign each  mastercons to one block without inducing conflicts  */
   SCIP_RETCODE postprocessMasterToBlocksConssAdjacency(
      Seeedpool* seeedpool, /**< a seeedpool that uses this seeed */
      SCIP_Bool* success
      );


   /** assigns all open constraints and open variables
     *  strategy: assigns all conss same block if they are connected
     *  two constraints are adjacent if there is a common variable
     *  this relies on the consadjacency structure of the seeedpool
     *  hence it cannot be applied in presence of linking variables */
    SCIP_RETCODE completeByConnectedConssAdjacency(
       Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
       );



   /** assigns all open constraints and open variables
    *  strategy: assigns a cons (and related vars) to any block if possible by means of prior var assignments
    *  and to master, if there does not exist such a block */
   SCIP_RETCODE completeGreedily(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** returns true if the given detector used a consclassifier */
   bool consClassifierUsed(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** assigns every open cons/var
    *  - to the respective block if it hits exactly one blockvar/blockcons and no open vars/conss
    *  - to master/linking if it hits blockvars/blockconss assigned to different blocks
    *  - and every cons to master that hits a master var
    *  - and every var to master if it does not hit any blockcons and has no open cons */
   SCIP_RETCODE considerImplicits(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** copies the given seeed's classifier statistics */
   SCIP_RETCODE copyClassifierStatistics(
      const Seeed* otherseeed
      );

   /** deletes empty blocks */
   SCIP_RETCODE deleteEmptyBlocks(
      bool variables
   );

   /** deletes a cons from list of open conss */
   SCIP_RETCODE deleteOpencons(
      int opencons
      );

   /** deletes a var from the list of open vars */
   SCIP_RETCODE deleteOpenvar(
      int openvar
      );

   SCIP_RETCODE displayAggregationInformation();

   /** displays the assignments of the conss */
   SCIP_RETCODE displayConss(Seeedpool* seeedpool);

   /** displays the relevant information of the seeed */
   SCIP_RETCODE displayInfo(
      Seeedpool* seeedpool, /**< a seeedpool that uses this seeed */
      int detailLevel /**< pass a value that indicates how detailed the output should be:
                              0: brief overview
                              1: block and detector info
                              2: cons and var assignments */
      );

   /*@todo is initialization in declaration necessary? */
   /** displays the relevant information of the seeed */
   SCIP_RETCODE displaySeeed(
      Seeedpool* seeedpool = NULL /**< a seeedpool that uses this seeed */
      );

   /*@todo is initialization in declaration necessary? */
   /** displays the assignments of the vars */
   SCIP_RETCODE displayVars(
      Seeedpool* seeedpool = NULL /**< a seeedpool that uses this seeed */
      );

   /** computes the score of the given seeed based on the border, the average density score and the ratio of linking
    * variables */
   SCIP_Real evaluate(
      Seeedpool* seeedpool, /**< a seeedpool that uses this seeed */
      SCORETYPE  type
      );

   /** assigns all conss to master or declares them to be open (and declares all vars to be open)
    *  according to the cons assignment information given in constoblock hashmap
    *  precondition: no cons or var is already assigned to a block */
   SCIP_RETCODE filloutBorderFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks,          /**< number of blocks the hashmap contains */
      Seeedpool* seeedpool       /**< a seeedpool that uses this seeed */
      );

   /** assigns all conss to master or a block
    *  according to the cons assignment information given in constoblock hashmap
    *  calculates implicit variable assignment through cons assignment
    *  precondition: no cons or var is already assigned to a block and constoblock contains information for every cons */
   SCIP_RETCODE filloutSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock, /**< hashmap assigning cons indices (not SCIP_Cons* !!) to block indices
                                   *< (master assignment is indicated by assigning cons to index additionalNBlocks) */
      int givenNBlocks,          /**< number of blocks the hashmap contains */
      Seeedpool* seeedpool       /**< a seeedpool that uses this seeed */
      );

   /** reassigns variables classified as linking to master if the variable only hits master conss */
   SCIP_RETCODE findVarsLinkingToMaster(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** reassigns variables classified as linking to stairlinking if the variable hits conss in exactly two consecutive
    * blocks */
   SCIP_RETCODE findVarsLinkingToStairlinking(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** returns a vector of pairs of var indices and vectors of (two) block indices
    *  the related linking variable hits exactly the two blocks given in the related vector */
   std::vector< std::pair< int, std::vector< int > > > findLinkingVarsPotentiallyStairlinking(
      Seeedpool* seeedpool
      );

   /** assigns all booked constraints and variables and deletes them from list of open cons and open vars */
   SCIP_RETCODE flushBooked();

   /** returns ancestor id of given ancestor */
   int getAncestorID(
      int ancestorindex /**< index of ancestor in listofancestorids data structure */
      );

   /** returns ancestor id of given ancestor */
   std::vector<int> getAncestorList(
      );

   void setAncestorList(
      std::vector<int> newlist
      );

   /** adds ancestor id of given ancestor */
   void addAncestorID(
      int ancestor
      );

   const std::vector<int> & getBlocksForRep(int repid);

   /** returns detectorchainstring */
   char* getDetectorChainString();

   /** returns detectorchain info of detector related to given detectorchain index */
   std::string getDetectorchainInfo(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns the time that the detector related to the given detectorchainindex needed for detecting */
   SCIP_Real getDetectorClockTime(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns the time that the detectors needed for detecting */
   std::vector<SCIP_Real> getDetectorClockTimes();


   std::string getComponentInformation(
            );

   /** returns the data of the consclassifier that the given detector made use of */
   SCIP_RETCODE getConsClassifierData(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      ConsClassifier** classifier, /**< a pointer to the used consclassifier */
      std::vector<int>& consclassesmaster /**< a vector containing all indices of the consclasses assigned to master */
      );

   /** returns array containing constraints assigned to a block */
   const int* getConssForBlock(
      int block
      );

   /** returns the detectorchain */
   DEC_DETECTOR** getDetectorchain();

   /** returns the detectorchain as a vector */
   std::vector<DEC_DETECTOR*> getDetectorchainVector();

   /** returns a string displaying all detector-related information, i.e. clock times and assignment data */
   std::string getDetectorStatistics(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns a string displaying classifier information if such a classifier was used */
   std::string getDetectorClassifierInfo(
      Seeedpool* seeedpool, /**< a seeedpool that uses this seeed */
      int detectorchainindex, /**< index of the detector in the detectorchain */
      bool displayConssVars /**< pass true if constraints and variables of the respective classes should be displayed */
      );

   /** returns true if this seeed was finished by finishSeeed() method of a detector */
   bool getFinishedByFinisher();

   /** returns true if the seeed is finished by a finisher in the unpresolved problem */
   bool getFinishedByFinisherUnpresolved();

   /** returns the detector that finished this seeed in the unpresolved problem if there exists one, NULL otherwise */
   DEC_DETECTOR* getFinishedUnpresolvedBy();

   /** returns the calculated hash value of this seeed */
   long getHashValue();

   /** returns the id of the seeed */
   int getID();

   /** returns array containing all linking vars */
   const int* getLinkingvars();

   /** returns array containing all master conss */
   const int* getMasterconss();

   /** returns array containing all master vars (hitting only constraints in the master) */
   const int* getMastervars();

   /** returns the "maximum white score" */
   SCIP_Real getMaxWhiteScore();

   /** returns the "maximum white score" */
   SCIP_Real getBendersScore();


   /** returns the number of nonzero coeffs in a certain block */
   int  getNCoeffsForBlock(
      gcg::Seeedpool* seeedpool,
      int blockid
      );

   /** returns the number of nonzero coeffs in master */
   int  getNCoeffsForMaster(
      gcg::Seeedpool* seeedpool
      );


   /** returns the score of the seeed (depending on used scoretype) */
   SCIP_Real getScore(
      SCORETYPE type
      );


   /* Are all master constraints set partitioning, set packing, set cover, or cardinality constraints */
   SCIP_Bool hasSetppccardMaster(
      gcg::Seeedpool* seeedpool
   );

   /* Are all master constraints set partitioning, set packing, or set cover constraints */
   SCIP_Bool hasSetppcMaster(
      gcg::Seeedpool* seeedpool
   );


   /* Are all master constraints set partitioning, or set packing constraints */
   SCIP_Bool hasSetppMaster(
      gcg::Seeedpool* seeedpool
   );



   /** returns whether this seeed is usergiven */
   USERGIVEN getUsergiven();

   /** returns number of ancestor seeeds */
   int getNAncestors();

   /** returns number of blocks */
   int getNBlocks();

   /** returns number of conss */
   int getNConss();

   /** returns size of the vector containing conss assigned to a block */
   int getNConssForBlock(
      int block
      );

   /** returns size of the detectorchain info vector */
   int getNDetectorchainInfo();

   /** returns the number of detectors the seeed is propagated by */
   int getNDetectors();

   /** returns the number used classifiers */
   int getNUsedClassifier();

   /** returns size of the vector containing linking vars */
   int getNLinkingvars();

   /** returns size of the vector containing master conss */
   int getNMasterconss();

   /** returns size of the vector containing master vars (hitting only constraints in the master) */
   int getNMastervars();

   /** returns number of blocks a detector added */
   int getNNewBlocks(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns number of blocks the detectors in the detectorchain added */
   std::vector<int> getNNewBlocksVector();

   /** returns total number of stairlinking vars */
   int getNTotalStairlinkingvars();

   /** returns size of vector containing constraints not assigned yet */
   int getNOpenconss();

   /** returns size of vector containing variables not assigned yet */
   int getNOpenvars();

   /** returns the number of blockrepresentatives */
   int getNReps();

   /** returns size of the vector containing stairlinking vars */
   int getNStairlinkingvars(
      int block
      );

   /** returns number of vars */
   int getNVars();

   /** returns size of the vector containing vars assigned to a block */
   int getNVarsForBlock(
      int block
      );

   /** returns array containing constraints not assigned yet */
   const int* getOpenconss();

   /** returns array containing constraints not assigned yet  as vector*/
   std::vector<int> getOpenconssVec();

   /** returns array containing variables not assigned yet */
   const int* getOpenvars();

   /** returns array containing variables not assigned yet as vector*/
   std::vector<int> getOpenvarsVec();

   /** returns fraction of variables assigned to the border for a detector */
   SCIP_Real getPctVarsToBorder(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of variables assigned to the border for detectors in detectorchain */
   std::vector<SCIP_Real> getPctVarsToBorderVector();

   /** returns fraction of variables assigned to a block for a detector */
   SCIP_Real getPctVarsToBlock(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of variables assigned to a block for detectors in detectorchain */
   std::vector<SCIP_Real> getPctVarsToBlockVector();

   /** returns fraction of variables that are not longer open for a detector */
   SCIP_Real getPctVarsFromFree(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of variables that are not longer open for detectors in detectorchain */
   std::vector<SCIP_Real> getPctVarsFromFreeVector();

   /** returns fraction of constraints assigned to the border for a detector */
   SCIP_Real getPctConssToBorder(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of constraints assigned to the border for detectors in detectorchain */
   std::vector<SCIP_Real> getPctConssToBorderVector();

   /** returns fraction of constraints assigned to a block for a detector */
   SCIP_Real getPctConssToBlock(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of constraints assigned to a block for detectors in detectorchain */
   std::vector<SCIP_Real> getPctConssToBlockVector();

   /** returns fraction of constraints that are not longer open for a detector */
   SCIP_Real getPctConssFromFree(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   /** returns fraction of constraints that are not longer open for detectors in detectorchain */
   std::vector<SCIP_Real> getPctConssFromFreeVector();

   /** returns index of the representative block */
   int getRepForBlock(
      int blockid
      );

   std::vector<int> & getRepVarmap(
      int repid,
      int blockrepid
      );

   /** returns the corresponding seeedpool */
   Seeedpool* getSeeedpool();

   /** returns array containing stairlinking vars */
   const int* getStairlinkingvars(
      int block
      );

   /** returns true if this seeed stems from the unpresolved problem */
   bool getStemsFromUnpresolved();

   /** returns the data of the varclassifier that the given detector made use of */
   SCIP_RETCODE getVarClassifierData(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      VarClassifier** classifier, /**< a pointer to the used varclassifier */
      std::vector<int>& varclasseslinking, /**< a vector containing all indices of the varclasses assigned to linking */
      std::vector<int>& varclassesmaster /**< a vector containing all indices of the varclasses assigned to master */
      );

   /** returns array containing vars of a block */
   const int* getVarsForBlock(
      int block
      );

   /** returns array containing vars of a block */
   int getVarProbindexForBlock(
      int varid,
      int block
   );

   void initOnlyBinMaster();

   SCIP_Bool isAgginfoToExpensive();


   /** returns true if this seeed is complete,
    *  i.e. it has no more open constraints and variables */
   bool isComplete();

   /** returns true if the cons is a cons of the block */
   bool isConsBlockconsOfBlock(
      int cons,
      int block
      );

   /** returns true if the cons is a master cons */
   bool isConsMastercons(
      int cons
      );

   /** returns true if the cons is an open cons */
   bool isConsOpencons(
      int cons
      );

   /** returns true if the seeed is from a detector operating in legacymode */
   bool isFromLegacymode();

   /** returns true if the seeed is from the unpresolved problem */
   bool isFromUnpresolved();

   /** returns true if the seeed is selected */
   bool isSelected();

   /* method to check whether this seeed is equal to a given other seeed (calls isEqual(Seeed*)) */
   SCIP_RETCODE isEqual(
      Seeed* otherseeed,   /**< other seeed */
      SCIP_Bool* isequal,  /**< pointer to store whether seeeds are identical */
      bool sortseeeds      /**< should conss and vars be sorted before comparing the seeeds? */
      );

   /* method to check whether this seeed is equal to a given other seeed */
   bool isEqual(
      Seeed* other /**< other seeed */
      );

   /** returns true if this seeed was propagated by a detector */
   bool isPropagatedBy(
      DEC_DETECTOR* detectorID
      );

   /** returns true if this seeed is trivial,
    *  i.e. all conss are in one block, all conss are in border, all variables linking or mastervars */
   bool isTrivial();

   /** returns true if the var is assigned to the block */
   bool isVarBlockvarOfBlock(
      int var,
      int block
      );

   /** returns true if the var is a linking var */
   bool isVarLinkingvar(
      int var
      );

   /** returns true if the var is a master var */
   bool isVarMastervar(
      int var
      );

   /** returns true if the var is an open var */
   bool isVarOpenvar(
      int var
      );

   /** returns true if the var is a stairlinking var */
   bool isVarStairlinkingvar(
      int var
      );

   /** returns true if the var is a stairlinkingvar of the block */
   bool isVarStairlinkingvarOfBlock(
      int var,
      int block
      );

   SCIP_RETCODE printClassifierInformation(
      SCIP*                scip,
      gcg::Seeedpool*      seeedpool,
      FILE*                file);



   /** refine seeed with focus on blocks: assigns open conss and vars if they can be found
    *  in blocks without respect to open vars and conss (assignHittingOpenconss(), assignHittingOpenvars()) */
   SCIP_RETCODE refineToBlocks(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** refine seeed with focus on master: do obvious (considerImplicits()) assignments and
    *  assign other conss and vars to master if possible (assignOpenPartialHittingToMaster()) */
   SCIP_RETCODE refineToMaster(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** registers statistics for a used consclassifier */
   void setConsClassifierStatistics(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      ConsClassifier* classifier, /**< the used consclassifier */
      std::vector<int> consclassesmaster /**< vector of classindices that were assigned to master */
      );

   /** directly adds a constraint to a block
    *  does not delete this cons from list of open conss */
   SCIP_RETCODE setConsToBlock(
      int consToBlock,
      int block
      );

   /** directly adds a constraint to the master constraints
    *  does not delete this cons from list of open conss */
   SCIP_RETCODE setConsToMaster(
      int consToMaster
      );

   /** sets the whole detectorchain */
   void setDetectorchain(
      std::vector<DEC_DETECTOR*> detectorChain
      );

   /** sets seeed to be propagated by a detector */
   SCIP_RETCODE setDetectorPropagated(
      DEC_DETECTOR* detectorID
      );

   /** sets seeed to be propagated by a finishing detector */
   SCIP_RETCODE setFinishingDetectorPropagated(
      DEC_DETECTOR* detectorID
      );

   /** sets whether this seeed was finished by a detector */
   void setFinishedByFinisher(
      bool finished
      );

   /** sets whether this seeed is finished by a finisher in the unpresolved problem */
   void setFinishedByFinisherUnpresolved(
      bool finishedByFinisherUnpresolved
      );

   /** sets the detector that finished the seeed in the unpresolved problem */
   void setFinishedUnpresolvedBy(
      DEC_DETECTOR* detector
      );

   /** sets whether this seeed stems from a detector operating in legacymode */
   void setLegacymode(
      bool legacymode
      );

   /** sets number of blocks, only increasing number allowed */
   SCIP_RETCODE setNBlocks(
      int nBlocks
      );

   /** sets the id of this seeed */
   SCIP_RETCODE setID(
      int id
      );

   /** sets whether this seeed is from the unpresolved problem */
   void setIsFromUnpresolved(
      bool unpresolved
      );

   /** sets whether this seeed is selected */
   void setSelected(
      bool selected
      );

   /** set the corresponding seeedpool */
   void setSeeedpool(
      Seeedpool* seeedpool
      );

   /** sets whether this seeed stems from an unpresolved problem seeed */
   void setStemsFromUnpresolved(
      bool stemsfromunpresolved
      );

   /** sets whether this seeed is usergiven */
   void setUsergiven(
      USERGIVEN usergiven
      );

   /** registers statistics for a used varclassifier */
   void setVarClassifierStatistics(
      int detectorchainindex, /**< index of the detector in the detectorchain */
      VarClassifier* classifier, /**< the used varclassifier */
      std::vector<int> varclasseslinking, /**< vector of classindices that were assigned to linking */
      std::vector<int> varclassesmaster /**< vector of classindices that were assigned to master */
      );

   /** directly adds a variable to the linking variables
    *  does not delete this var from list of open vars */
   SCIP_RETCODE setVarToBlock(
      int varToBlock,
      int block
      );

   /** directly adds a variable to the linking variables
    *  does not delete this var from list of open vars */
   SCIP_RETCODE setVarToLinking(
      int varToLinking
      );

   /** directly adds a variable to the master variables (hitting only constraints in the master)
    *  does not delete this var from list of open vars */
   SCIP_RETCODE setVarToMaster(
      int varToMaster
      );

   /** directly adds a variable to the stairlinking variables
    *  does not delete this var from list of open vars */
   SCIP_RETCODE setVarToStairlinking(
      int varToStairLinking,
      int block1,
      int block2
      );

   /** generates and opens a gp visualization of the seeed */
   void showVisualisation();

   /** returns true if this seeed is a userseeed that should be completed by setting unspecified constraints to master */
   SCIP_Bool shouldCompletedByConsToMaster();

   /** sorts the vars and conss by their indices */
   void sort();

   /** returns a short caption for this seeed */
   const char* getShortCaption();

   /** sets the detector chain short string */
   SCIP_RETCODE setDetectorChainString(
      char* detectorchainstring
      );

   void setNNewBlocksVector(
      std::vector<int>  newvector
);

   void setPctConssToBlockVector(
      std::vector<SCIP_Real> newvector
      );

   /** returns fraction of constraints that are not longer open for detectors in detectorchain */
   void setPctConssFromFreeVector(
      std::vector<SCIP_Real> newvector
      );

   /** returns fraction of constraints assigned to the border for detectors in detectorchain */
   void setPctConssToBorderVector(
      std::vector<SCIP_Real> newvector
      );

   /** returns fraction of variables assigned to the border for detectors in detectorchain */
   void setPctVarsToBorderVector(
      std::vector<SCIP_Real> newvector
      );

   /** returns fraction of variables assigned to a block for detectors in detectorchain */
   void setPctVarsToBlockVector(
      std::vector<SCIP_Real> newvector
   );

   /** returns fraction of variables that are not longer open for detectors in detectorchain */
   void setPctVarsFromFreeVector(
      std::vector<SCIP_Real> newvector
      );

   /** returns the time that the detectors needed for detecting */
   void setDetectorClockTimes(
      std::vector<SCIP_Real> newvector
      );

   /** returns true if the given detector used a varclassifier */
   bool varClassifierUsed(
      int detectorchainindex /**< index of the detector in the detectorchain */
      );

   SCIP_RETCODE writeAsDec(
      FILE* file,
      //GCG_PROBLEM_TRANSFORMED_STATUS transformed,
      Seeedpool*   seeedpool,
      SCIP_RESULT* result
      );


   /** creates and sets a detector chain short string for this seeed */
   SCIP_RETCODE buildDecChainString();

private:

   /** adds empty entries for all classifier statistics for a detector added to the detector chain */
   void addEmptyClassifierStatistics();

   /** assigns every open cons
    *  - to master if it hits blockvars of different blocks
    *  - to the respective block if it hits a blockvar of exactly one block and no stairlinking var
    *  - to master if it hits a stairlinking var but there is no block the cons may be assigned to
    *  - to the block with the lowest number of conss if it hits a stairlinking var and there are blocks the cons may be
    *    assigned to
    *  returns true if there is a cons that has been assigned */
   bool assignHittingOpenconss(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** assigns every open var
    *  - to the respective block if it hits blockconss of exactly one block
    *  - to linking if it hits blockconss of more than one different blocks
    *  returns true if there is a var that has been assigned */
   bool assignHittingOpenvars(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** assigns every open cons to master that hits
    *  - exactly one block var and at least one open var or
    *  - a master var */
   SCIP_RETCODE assignOpenPartialHittingConsToMaster(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** assigns open conss/vars that hit exactly one block and at least one open var/cons to border */
   SCIP_RETCODE assignOpenPartialHittingToMaster(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );
   /** assigns every open var to linking that hits
    *  - exactly one block cons and at least one open cons */
   SCIP_RETCODE assignOpenPartialHittingVarsToMaster(
      Seeedpool* seeedpool /**< a seeedpool that uses this seeed */
      );

   /** calculates the number of nonzero coefficients for the blocks */
   SCIP_RETCODE calcNCoeffsForBlocks(
   Seeedpool*   seeedpool
   );


   void calcmaxwhitescore();

   void calcbendersscore();

   SCIP_RETCODE calcclassicscore();

   void calcborderareascore();

   void calcmaxforeseeingwhitescore();

   void calcmaxforeseeingwhitescoreagg();

   void calcsetpartfwhitescore();

   void calcsetpartfwhitescoreagg();

   void calcbenderareascore();

   void calcblockareascore();

   void calcblockareascoreagg();


};

} /* namespace gcg */
#endif /* GCG_CLASS_Seeed_H__ */
