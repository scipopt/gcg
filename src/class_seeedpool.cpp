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

/**@file   class_colpool.cpp
 * @brief  class with functions for colpool
 * @author Michael Bastubbe
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg.h"
#include "objscip/objscip.h"
#include "class_seeedpool.h"
#include "struct_detector.h"

#include <exception>

#define SCIP_CALL_EXC(x)   do                                                                                  \
                       {                                                                                      \
                          SCIP_RETCODE _restat_;                                                              \
                          if( (_restat_ = (x)) !=  SCIP_OKAY )                                                \
                          {                                                                                   \
                             SCIPerrorMessage("Error <%d> in function call\n", _restat_);                     \
                             throw std::exception();                                                          \
                           }                                                                                  \
                       }                                                                                      \
                       while( FALSE )

struct Seeed_Propagation_Data
{
	gcg::Seeedpool* seeedpool;
	gcg::Seeed* seeedToPropagate;
	std::vector<gcg::Seeed*> newSeeeds;
};


/** constraint handler data */
struct SCIP_ConshdlrData
{
   DEC_DECOMP**          decdecomps;         /**< array of decomposition structures */
   DEC_DETECTOR**        detectors;          /**< array of structure detectors */
   int*                  priorities;         /**< priorities of the detectors */
   int                   ndetectors;         /**< number of detectors */
   SCIP_CLOCK*           detectorclock;      /**< clock to measure detection time */
   SCIP_Bool             hasrun;             /**< flag to indicate whether we have already detected */
   int                   ndecomps;           /**< number of decomposition structures  */
   SCIP_Bool             createbasicdecomp;  /**< indicates whether to create a decomposition with all constraints in the master if no other specified */
};

namespace gcg {

/** local methods */

SCIP_CONS* consGetRelevantRepr(SCIP* scip, SCIP_CONS* cons){

		return cons;
}

SCIP_VAR* varGetRelevantRepr(SCIP* scip, SCIP_VAR* var){

		return var;
}

SCIP_Bool seeedIsNoDuplicate(SeeedPtr seeed, std::vector<SeeedPtr> const & currSeeeds, std::vector<SeeedPtr> const & finishedSeeeeds){

	return TRUE;
}


/** constructor */
 Seeedpool::Seeedpool(
    SCIP*             	givenScip, /**< SCIP data structure */
	const char*  		conshdlrName
    ):scip(givenScip), currSeeeds(0), nTotalSeeeds(0),nVars(SCIPgetNVars(givenScip) ), nConss(SCIPgetNConss(givenScip) ), nDetectors(0)
 {
	 SCIP_CONS** conss;
	 SCIP_VAR** vars;

	 SCIP_CONSHDLR* conshdlr;  /** cons_decomp to get detectors */
	 SCIP_CONSHDLRDATA* conshdlrdata;


	 int relevantVarCounter = 0;
	 int relevantConsCounter = 0;

	 /** store all enabled detectors */

	 conshdlr = SCIPfindConshdlr(scip, conshdlrName);
	 assert(conshdlr != NULL);
	 conshdlrdata = SCIPconshdlrGetData(conshdlr);
	 assert(conshdlrdata != NULL);

	 for(int d = 0; d < conshdlrdata->ndetectors; ++d )
	 {
		 DEC_DETECTOR *detector;
		 detector = conshdlrdata->detectors[d];
		 assert(detector != NULL);
		 conshdlrdata->priorities[d] = detector->priority;
	 }

	 SCIPdebugMessage("Sorting %i detectors\n", conshdlrdata->ndetectors);
	 SCIPsortIntPtr(conshdlrdata->priorities, (void**)conshdlrdata->detectors, conshdlrdata->ndetectors);

	 SCIPdebugMessage("Trying %d detectors.\n", conshdlrdata->ndetectors);
	 for(int d = 0; d < conshdlrdata->ndetectors; ++d )
	 {
		 DEC_DETECTOR* detector;

		 detector = conshdlrdata->detectors[d];
		 assert(detector != NULL);
		 if( !detector->enabled || detector->propagateSeeed == NULL)
			 continue;

		 scipDetectorToIndex[detector] = nDetectors;
		 detectorToScipDetector.push_back(detector);
		 ++nDetectors;

	 }


	 /** initilize matrix datastructures */
	 conss = SCIPgetConss(scip);
	 vars = SCIPgetVars(scip);

	 /** assign an index to every cons and var
	  * @TODO: are all constraints/variables relevant? (probvars etc)  */
	 for(int i = 0; i < nConss; ++i)
	 {
		 SCIP_CONS* relevantCons;

		 relevantCons = consGetRelevantRepr(scip, conss[i]);
		 if( relevantCons != NULL )
		 {
			 scipConsToIndex[relevantCons] = relevantConsCounter ;
			 consToScipCons.push_back(relevantCons);
			 ++relevantConsCounter;
		 }
	 }

	 for(int i = 0; i < nVars; ++i)
	 {
		 SCIP_VAR* relevantVar;

		 relevantVar = varGetRelevantRepr(scip, vars[i]);

		 if( relevantVar != NULL )
		 {
			 scipVarToIndex[relevantVar] = relevantVarCounter ;
			 varToScipVar.push_back(relevantVar);
			 ++relevantVarCounter;
		 }
	 }

	 /** from here on nVars and nConss represents the relevant numbers */
	 nVars = relevantVarCounter;
	 nConss = relevantConsCounter;
	 varsForConss = std::vector<std::vector<int>>(nConss);
	 conssForVars = std::vector<std::vector<int>>(nVars);

	 assert(varToScipVar.size() == nVars);
	 assert(consToScipCons.size() == nConss);

	 /** assumption: now every relevant constraint and variable has its index and is stored in the corresponding unordered_map */
	 /** find constraint <-> variable relationships and store them in both directions */
	 for(int i = 0; i < (int)consToScipCons.size() ; ++i)
	 {
		 SCIP_CONS* cons;
		 SCIP_VAR** currVars;
		 int 		nCurrVars;
		 SCIP_Bool  success;

		 cons = consToScipCons[i];

		 SCIP_CALL_ABORT( SCIPgetConsNVars(scip, cons, &nCurrVars, &success ) );
		 assert(success);

		 SCIP_CALL_ABORT( SCIPallocBufferArray(scip, &currVars, nCurrVars) );
		 SCIPgetConsVars(scip, cons, currVars, nCurrVars, &success );

		 for(int currVar = 0; currVar < nCurrVars; ++currVar)
		 {
			 int varIndex;

			 std::tr1::unordered_map<SCIP_VAR*, int>::const_iterator iterVar = scipVarToIndex.find(currVars[currVar]);

			 if(iterVar == scipVarToIndex.end() )
				 continue;

			 varIndex = iterVar->second;

			 varsForConss[i].push_back(varIndex);
			 conssForVars[varIndex].push_back(i);

		 }
	 }

	 /* populate seeed vector with empty seeed */

	 currSeeeds.push_back(new Seeed(nTotalSeeeds,nDetectors,nConss,nVars) );

	 nTotalSeeeds++;


 }//end constructor

 Seeedpool::~Seeedpool(){

 }


 /** finds decompositions  */
 DEC_DECOMP**    Seeedpool::findDecompostions(
 ){

	 /** 1) read parameter, as there are: maxrounds
	  *  2) loop rounds
	  *  3) every seeed in seeeds
	  *  4) every detector not registered yet propagetes seeed
	  *  5)  */

	 int maxRounds;
	 SEEED_PROPAGATION_DATA* seeedPropData;
	 DEC_DECOMP** decompostions;

	 maxRounds = 2;
	 seeedPropData = new SEEED_PROPAGATION_DATA();
	 seeedPropData->seeedpool = this;

	 for(int round = 0; round < maxRounds; ++round)
		 for(size_t s = 0; s < currSeeeds.size(); ++s )
		 {
			 SeeedPtr seeedPtr;
			 seeedPtr= currSeeeds[s];

			 for(int d = 0; d < nDetectors; ++d)
			 {
				 std::vector<SeeedPtr>::const_iterator newSIter;
				 std::vector<SeeedPtr>::const_iterator newSIterEnd;

				 SCIP_RESULT* result;

				 if(seeedPtr->isPropagatedBy(d) )
					 continue;

				 seeedPropData->seeedToPropagate = seeedPtr;
				 seeedPropData->newSeeeds = std::vector<SeeedPtr>(0);

				 SCIP_CALL_ABORT(detectorToScipDetector[d]->propagateSeeed(scip, seeedPropData, result) );

				 assert(seeedPtr->isPropagatedBy(d));

				 newSIter = seeedPropData->newSeeeds.begin();
				 newSIterEnd = seeedPropData->newSeeeds.end();

				 for(; newSIter != newSIterEnd; ++newSIter)
				 {
					 if(seeedIsNoDuplicate(*newSIter, currSeeeds, finishedSeeeds ))
					 {
						 currSeeeds.push_back(*newSIter);
					 }
				 }
			 }
		 }

	 delete seeedPropData;

	 return decompostions;
 }

 /** access coefficient matrix constraint-wise */
 std::vector<int> const & Seeedpool::getVarsForCons(int cons){
	 return varsForConss[cons];
 }

 /** access coefficient matrix variable-wise */
 std::vector<int> const & Seeedpool::getConssForVar(int var){
	 return conssForVars[var];
 }




} /* namespace gcg */
