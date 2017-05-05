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

/**@file   class_Seeed.h
 * @brief  class with functions for seeed
 * @author Michael Bastubbe
 * @author Hannah Hechenrieder
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef GCG_CLASS_SEEED_H__
#define GCG_CLASS_SEEED_H__


#include "objscip/objscip.h"
#include <vector>
#include "struct_detector.h"
#include <string>



namespace gcg {

enum USERGIVEN
{
   NOT = 0,
   PARTIAL = -1,
   COMPLETE = -2,
   COMPLETED_CONSTOMASTER = -3
};

class Seeedpool;


class Seeed
{

private:
	SCIP*							         scip;
   int								      id;						         /**< id of the seeed */
   int 								      nBlocks;				            /**< number of blocks the decomposition currently has */
   int 								      nVars;                        /**< number of variables */
   int								      nConss;                       /**< number of constraints */
   std::vector<int>					   masterConss;			         /**< vector containing indices of master constraints */
   std::vector<int>					   masterVars;				         /**< vector containing indices of master variables */
   std::vector<std::vector<int>>    conssForBlocks; 		         /**< conssForBlocks[k] contains a vector of indices of all constraints assigned to block k */
   std::vector<std::vector<int>>    varsForBlocks; 			      /**< varsForBlocks[k] contains a vector of indices of all variables assigned to block k */
   std::vector<int> 				      linkingVars;			         /**< vector containing indices of linking variables */
   std::vector<std::vector<int>>    stairlinkingVars;		         /**< vector containing indices of staircase linking variables of the blocks (stairlinking variables can be found only in the vector of their first block) */
   std::vector<int> 				      openVars;				         /**< vector containing indices of variables that are not assigned yet*/
   std::vector<int> 				      openConss;				         /**< vector containing indices of constraints that are not assigned yet*/
   std::vector<int>                 bookedAsMasterConss;          /**< vector containing indices of constraints that are not assigned yet but booked as master conss */
   std::vector<std::pair<int, int>> bookedAsBlockConss;           /**< vector containing indices of constraints that are not assigned yet but booked as block conss and the block */
   std::vector<int>                 bookedAsLinkingVars;          /**< vector containing indices of variables that are not assigned yet but booked as linking vars */
   std::vector<int>                 bookedAsMasterVars;           /**< vector containing indices of variables that are not assigned yet but booked as master vars */
   std::vector<std::pair<int, int>> bookedAsBlockVars;            /**< vector containing indices of variables that are not assigned yet but booked as block vars and the block */
   std::vector<std::pair<int, int>> bookedAsStairlinkingVars;     /**< vector containing indices of variables that are not assigned yet but booked as stairlinking vars and the first block of the stairlinking var */
   std::vector<bool> 				   propagatedByDetector;	      /**< propagatedByDetector[i] is this seeed propagated by detector i */
   bool 							         openVarsAndConssCalculated;   /**< are the openVars and openCons calculated */
   long                             hashvalue;
   SCIP_Real                        score;                        /**< score to evaluate the seeeds */
   SCIP_Real                        maxwhitescore;                /**< score corresponding to the max white measure */

   bool                             changedHashvalue;             /**< are there any changes concerning the hash value since it was calculated last time */

   const static int                 primes[];
   const static int                 nPrimes;

public:

   bool                             isFinishedByFinisher;         /**< was this seeed finished by the finishseeed() method of a detector */
   /** statistic information */
   std::vector<DEC_DETECTOR*>       detectorChain;              /**< vector containing detectors that worked on that seeed */
   std::vector<std::string>         detectorchaininfo;          /**< vector containing information about the detector call */
   std::vector<SCIP_Bool>           detectorChainFinishingUsed; /**< vector containing whether the finishing method of the corresponding detector was used on that seeed */
   std::vector<SCIP_Real>           detectorClockTimes;         /**< vector containing detector times in seconds  */
   std::vector<SCIP_Real>           pctVarsToBorder;            /**< vector containing the fraction of variables assigned to the border for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctVarsToBlock;             /**< vector containing the fraction of variables assigned to a block for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctVarsFromFree;            /**< vector containing the fraction of variables that are not longer open for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssToBorder;           /**< vector containing the fraction of constraints assigned to the border for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssToBlock;            /**< vector containing the fraction of constraints assigned to a block for each detector working on that seeed*/
   std::vector<SCIP_Real>           pctConssFromFree;           /**< vector containing the fraction of constraints that are not longer open for each detector working on that seeed*/
   std::vector<int>                 nNewBlocks;                 /**< vector containing detector indices that worked on that seeed */

   std::vector<int>                 listofancestorids;          /**< vector containing detector indices that worked on that seeed */

   USERGIVEN                        usergiven;                  /**< is this seeed partially or completely given by user */

   /** datastructure to store information if this seeed stems from a seeed concerning the unpresolved problem */
   bool                             stemsFromUnpresolved;
   bool                             isFinishedByFinisherUnpresolved; /**< was the ancestor seeed for the unpresolved problem finished by the finishseeed() method of a detector */
   DEC_DETECTOR*                    finishedUnpresolvedBy;           /**< index of dinishing detector of unpresolved ancestor seeed */


   /** constructor */
   Seeed(
      SCIP* scip,                   /**< scip data structure */
      int id,      		   	      /**< id that is given to this seeed */
      int nDetectors,               /**< number of detectors */
      int nConss,				         /**< number of constraints */
      int nVars				         /**< number of variables */
   );

   /** copy constructor */
   Seeed(
      const Seeed *seeedToCopy,     /**< seeed to be copied */
      Seeedpool* seeedpool          /**< TODO why? */
   );

   /** destructor */
   ~Seeed(
   );


   /** adds a block, returns the number of the new block */
   int addBlock(
   );

   /** incorporates the the needed time of a certain detector in the detector chain */
   void addClockTime(
      SCIP_Real clocktime           /**< time to be added */
   );

   /** incorporates the changes from ancestor seeed */
   void addDecChangesFromAncestor(
      Seeed* ancestor
   );

   /** adds a detector chain info */
   void addDetectorChainInfo(
      const char* decinfo
   );

   /** returns true if at least one constraint is assigned to a block */
   bool alreadyAssignedConssToBlocks();

   /** returns true if the open vars and conss are calculated */
   bool areOpenVarsAndConssCalculated(
   );

   /** assigns open conss and vars if they can be found in blocks */
   SCIP_RETCODE assignAllDependent(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** fills out the border of the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignBorderFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns open vars to stairlinking if they can be found in two consecutive blocks, returns true if stairlinkingvars are assigned */
   bool assignCurrentStairlinking(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns open conss if they include blockvars, returns true if open conss are assigned */
   bool assignHittingOpenconss(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns open vars if they can be found in one block, returns true if open vars are assigned */
   bool assignHittingOpenvars(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** fills out the seeed with the hashmap constoblock if there are still assigned conss and vars */
   SCIP_RETCODE assignSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** fills out the seeed with the vector constoblock  */
   SCIP_RETCODE assignSeeedFromConstoblockVector(
      std::vector<int> constoblock,
      int additionalNBlocks,
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
      int consToMaster
   );

   /** books a variable to be added to the master variables (after calling flushBooked) */
   SCIP_RETCODE bookAsMasterVar(
      int varToMaster
   );

   /** books a variable to be added to the linking variables (after calling flushBooked) */
   SCIP_RETCODE bookAsLinkingVar(
      int varToLinking
   );

   /** books a variable to be added to the stairlinking variables of the given block and the following block (after calling flushBooked) */
   SCIP_RETCODE bookAsStairlinkingVar(
      int varToStairlinking,
      int firstBlock
   );

   /** calculates the hashvalue of the seeed for comparing */
   void calcHashvalue(
   );

   /** calculates vector containing constraints not assigned yet */
   void  calcOpenconss(
   );

   /** constructs vector containing variables not assigned yet */
   void  calcOpenvars(
   );

   /** returns whether all constraints are assigned and deletes the vector open cons if all are assigned */
   bool checkAllConsAssigned(
   );

   /** TODO checks the consistency of this seeed */
   bool checkConsistency(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** TODO checks the consistency of this seeed */
   bool checkConsistency(
   );

   /** TODO checks the consistency of this seeed partially */
   bool checkVarsAndConssConsistency(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns all open constraints and open variables
    *  strategy: assigns all conss and vars to the same block if they are indirectly connected
    *  a cons and a var are directly connected if the var appears in the cons */
   SCIP_RETCODE completeByConnected(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns all open constraints and open variables
    *  strategy: assigns a cons (and related vars) to any block if possible by means of prior var assignments
    *  and to master, if there does not exist such a block */
   SCIP_RETCODE completeGreedily(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns the open cons and open vars which are implicitly assigned */
   SCIP_RETCODE considerImplicits(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** deletes empty blocks */
   SCIP_RETCODE deleteEmptyBlocks(
   );

   /** deletes a cons from list of open conss */
   SCIP_RETCODE deleteOpencons(
      int opencons
   );

   /** deletes a var from the list of open vars */
   SCIP_RETCODE deleteOpenvar(
      int openvar
   );

   /** displays the assignments of the conss */
   SCIP_RETCODE displayConss(
   );

   /** displays the relevant information of the seeed */
   SCIP_RETCODE displaySeeed(
      Seeedpool* seeedpool = NULL   /**< a seeedpool that uses this seeed */
   );

   /** displays the assignments of the vars */
   SCIP_RETCODE displayVars(
      Seeedpool* seeedpool = NULL   /**< a seeedpool that uses this seeed */
   );

   /** computes the score of the given seeed based on the border, the average density score and the ratio of linking variables*/
   SCIP_Real evaluate(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** fills out the border of the seeed with the hashmap constoblock */
   SCIP_RETCODE filloutBorderFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** fills out a seeed with the hashmap constoblock */
   SCIP_RETCODE filloutSeeedFromConstoblock(
      SCIP_HASHMAP* constoblock,
      int givenNBlocks,
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** finds linking-variables that are actually master-variables, i.e. the variable is adjacent to only master-constraint */
   SCIP_RETCODE findVarsLinkingToMaster(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** finds linking-variables that are actually stairlinking-variables, i.e. the variable is adjacent to constraints in exactly two block */
   SCIP_RETCODE findVarsLinkingToStairlinking(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns all booked constraints and variables and deletes them from list of open cons and open vars */
   SCIP_RETCODE flushBooked(
   );

   /** returns array containing constraints assigned to a block */
   const int* getConssForBlock(
      int block
   );

   /** returns the detectorchain */
   DEC_DETECTOR** getDetectorchain(
   );

   /** returns true if this seeed was finished by finishSeeed() method of a detector */
   bool getFinishedByFinisher(
   );

   /** returns the calculated hash value of this seeed */
   long getHashValue(
   );

   /** returns the id of the seeed */
   int getID(
   );

   /** returns array containing all linking vars */
   const int* getLinkingvars(
   );

   /** returns array containing all master conss */
   const int* getMasterconss(
   );

   /** returns array containing all master vars (hitting only constraints in the master) */
   const int* getMastervars(
   );

   /** returns the "maximum white score" (the smaller the better) */
   SCIP_Real getMaxWhiteScore();

   /** returns number of blocks */
   int getNBlocks(
   );

   /** returns number of conss */
   int getNConss(
   );

   /** returns size of the vector containing conss assigned to a block */
   int getNConssForBlock(
      int block
   );

   /** returns the number of detectors the seeed is propagated by */
   int getNDetectors(
   );

   /** returns size of the vector containing linking vars */
   int getNLinkingvars(
   );

   /** returns size of the vector containing master conss */
   int getNMasterconss(
   );

   /** returns size of the vector containing master vars (hitting only constraints in the master) */
   int getNMastervars(
   );

   /** returns total number of stairlinking vars */
   int getNTotalStairlinkingvars(
   );

   /** returns size of vector containing constraints not assigned yet */
   int getNOpenconss(
   );

   /** returns size of vector containing variables not assigned yet */
   int getNOpenvars(
   );

   /** returns size of the vector containing stairlinking vars */
   int getNStairlinkingvars(
      int block
   );

   /** returns number of vars */
   int getNVars(
   );

   /** returns size of the vector containing vars assigned to a block */
   int getNVarsForBlock(
      int block
   );

   /** returns array containing constraints not assigned yet */
   const int* getOpenconss(
   );

   /** returns array containing variables not assigned yet */
   const int* getOpenvars(
   );

   /** returns array containing stairlinking vars */
   const int* getStairlinkingvars(
      int block
   );

   /** returns array containing vars of a block */
   const int* getVarsForBlock(
      int block
   );

   /** returns true if this seeed is complete,
    *  i.e. it has at no more open constraints and variables */
   bool isComplete(
   );

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

   /* method to check whether this seeed is equal to given other seeed (calls isEqual(Seeed*)) */
   SCIP_RETCODE isEqual(
      Seeed* otherseeed,            /**< other seeed */
      SCIP_Bool* isequal,           /**< pointer to store whether seeeds are identical */
      bool sortseeeds               /**< should conss and vars be sorted before comparing the seeeds? */
   );

   /* method to check whether this seeed is equal to given other seeed */
   bool isEqual(
      Seeed* other                  /**< other seeed */
   );

   /** returns true if this seeed was propagated by a detector */
   bool isPropagatedBy(
      DEC_DETECTOR* detectorID
   );

   /** returns true if this seeed is trivial,
    *  i.e. all conss are in one block, all conss are in border, all variables linking or mastervars */
   bool isTrivial(
   );

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

   /** refine seeed: do obvious (considerImplicits()) and some non-obvious assignments assignOpenPartialHittingToMaster() */
   SCIP_RETCODE refineToMaster(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** add a constraint to a block */
   SCIP_RETCODE setConsToBlock(
      int consToBlock,
      int block
   );

   /** add a constraint to the master constraints */
   SCIP_RETCODE setConsToMaster(
      int consToMaster
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

   /** sets number of blocks, only increasing number allowed */
   SCIP_RETCODE setNBlocks(
      int nBlocks
   );

   /** sets the id of this seeed */
   SCIP_RETCODE setID(
      int id
   );

   /** sets whether open vars and conss are calculated */
   SCIP_RETCODE setOpenVarsAndConssCalculated(
      bool value
   );

   /** adds a variable to a block */
   SCIP_RETCODE setVarToBlock(
      int varToBlock,
      int block
   );

   /** adds a variable to the linking variables */
   SCIP_RETCODE setVarToLinking(
      int varToLinking
   );

   /** adds a variable to the master variables (hitting only constraints in the master) */
   SCIP_RETCODE setVarToMaster(
      int varToMaster
   );

   /** adds a variable to the stairlinking variables */
   SCIP_RETCODE setVarToStairlinking(
      int varToStairLinking,
      int block1,
      int block2
   );

   /** displays seeed data in a scatterplot */
   void showScatterPlot(
      Seeedpool* seeedpool,         /**< a seeedpool that uses this seeed */
      SCIP_Bool writeonly = FALSE,
      const char* filename = NULL,
      SCIP_Bool draft = FALSE,
      SCIP_Bool colored = TRUE
   );

   /** returns true if this seeed is a userseeed that should be completed by setting unspecified constraints to master */
   SCIP_Bool shouldCompletedByConsToMaster(
   );

   /** sorts the vars and conss by their indices */
   void sort(
   );

   /** displays the assignments of the vars */
   SCIP_RETCODE writeScatterPlot(
      Seeedpool* seeedpool,         /**< a seeedpool that uses this seeed */
      const char* filename
   );

   /** returns a short caption for this seeed */
   const char* getShortCaption(
   );

private:

   /** assigns open conss (and vars) that hit a block and other open vars (or conss) that are open to border */
   SCIP_RETCODE assignOpenPartialHittingToMaster(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns open conss that hit a block and other open vars that are open to border */
   SCIP_RETCODE assignOpenPartialHittingConsToMaster(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );

   /** assigns open vars that hit a block and other open conss that are open to border */
   SCIP_RETCODE assignOpenPartialHittingVarsToMaster(
      Seeedpool* seeedpool          /**< a seeedpool that uses this seeed */
   );


   //  bool compare_blocks(int a, int b);

};

} /* namespace gcg */
#endif /* GCG_CLASS_Seeed_H__ */
