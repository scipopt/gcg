/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Colum Generation                                 */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id$"
//#define SCIP_DEBUG
/**@file   branch_generic.c
 * @ingroup BRANCHINGRULES
 * @brief  branching rule based on vanderbeck's generic branching scheme
 * @author Marcel Schmickerath
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#define SCIP_DEBUG

#include <assert.h>
#include <string.h>

#include "branch_generic.h"
#include "relax_gcg.h"
#include "cons_masterbranch.h"
#include "cons_origbranch.h"
#include "pricer_gcg.h"
#include "scip/cons_linear.h"
#include "type_branchgcg.h"
#include "pub_gcgvar.h"
#include "event_genericbranchvaradd.h"

#include <stdio.h>
#include <stdlib.h>

#define BRANCHRULE_NAME          "generic"
#define BRANCHRULE_DESC          "generic branching rule by Vanderbeck"
#define BRANCHRULE_PRIORITY      999999
#define BRANCHRULE_MAXDEPTH      -1
#define BRANCHRULE_MAXBOUNDDIST  1.0

typedef SCIP_Real ComponentBoundSequence[3];  // [[comp], [sense], [bound]]    sense=1 means >=, sense=0 means <

struct GCG_BranchData
{
	ComponentBoundSequence**   C;             /**< S[k] bound sequence for block k */ //!!! sort of each C[i]=S[i] is important !!!
	int*               sequencesizes;                 /**< number of bounds in S[k] */
	int                Csize;
	ComponentBoundSequence*   S;             /**< component bound sequence which induce the current branching constraint */
	int                Ssize;
	int                blocknr;             /**< number of block branching was performed */
	int                childnr;
	SCIP_Real          lhs;
	int                nchildNodes;
	SCIP_Real*         childlhs;
	SCIP_CONS*         mastercons;          /**< constraint enforcing the branching restriction in the master problem */
	GCG_BRANCHDATA**   childbranchdatas;
};

struct GCG_Strip
{
	SCIP_VAR*          mastervar;             /**< master variable */
	SCIP_Real          mastervarValue;
	int                blocknr;               /**< number of the block in which the strip belong */
	SCIP_Real*         generator;             /**< corresponding generator to the mastervar */
	SCIP_Bool*         compisinteger;         /**< ? is comp with integer origvars? */
	int                generatorsize;
	ComponentBoundSequence**   C;             /**< often NULL, only needed for ptrilocomp */
	int                Csize;
	int*               sequencesizes;
};



/** set of component bounds in separate */
struct GCG_Record
{
	ComponentBoundSequence**   record;             /**< returnvalue of separte function */
	int                recordsize;
	int*               sequencesizes;
};
typedef struct GCG_Record GCG_RECORD;

/*
 * Callback methods
 */

/* define not used callback as NULL*/
#define branchCopyGeneric NULL
#define branchFreeGeneric NULL
#define branchExitGeneric NULL
#define branchInitsolGeneric NULL
#define branchExitsolGeneric NULL

/*
 * branching specific interface methods
 */

/** method for calculating the generator of mastervar*/
static
SCIP_RETCODE getGenerators(SCIP* scip, SCIP_Real** generator, int* generatorsize, SCIP_Bool** compisinteger, int blocknr, SCIP_VAR** mastervars, int nmastervars, SCIP_VAR* mastervar)
{
	int i;
	int j;
	int k;
	SCIP_VAR** origvarsunion;
	SCIP_VAR** origvars;
	SCIP_Real* origvals;
	int norigvars;
	int nvarsinblock;

	i = 0;
	j = 0;
	k = 0;
	*generatorsize = 0;
	nvarsinblock = 0;
	origvarsunion = NULL;
	assert(mastervars != NULL);

	for(i=0; i<nmastervars; ++i)
	{
		origvars = GCGmasterVarGetOrigvars(mastervars[i]);
		norigvars = GCGmasterVarGetNOrigvars(mastervars[i]);

		if(blocknr != GCGvarGetBlock(mastervars[i]))
			continue;
		else
			++nvarsinblock;
		if(*generatorsize==0 && norigvars>0)
		{
			*generatorsize = norigvars;
			SCIP_CALL( SCIPallocMemoryArray(scip, generator, *generatorsize) );
			SCIP_CALL( SCIPallocMemoryArray(scip, compisinteger, *generatorsize) );
			SCIP_CALL( SCIPallocMemoryArray(scip, &origvarsunion, *generatorsize) );
			for(j=0; j<*generatorsize; ++j)
			{
				origvarsunion[j] = origvars[j];
				(*generator)[j] = 0;
				(*compisinteger)[j] = TRUE;
				if(SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_IMPLINT)
					(*compisinteger)[j] = FALSE;
			}
		}
		else
		{
			for(j=0; j<norigvars; ++j)
			{
				int oldgeneratorsize;

				oldgeneratorsize = *generatorsize;

				for(k=0; k<oldgeneratorsize; ++k)
				{
					if(origvarsunion[k] == origvars[j])
					{
						break;
					}
					if(k == norigvars-1)
					{
						++(*generatorsize);
						SCIP_CALL( SCIPreallocMemoryArray(scip, generator, *generatorsize) );
						SCIP_CALL( SCIPreallocMemoryArray(scip, compisinteger, *generatorsize) );
						SCIP_CALL( SCIPreallocMemoryArray(scip, &origvarsunion, *generatorsize) );
						origvarsunion[*generatorsize-1] = origvars[j];
						(*generator)[*generatorsize-1] = 0;
						(*compisinteger)[(*generatorsize)-1] = TRUE;
						if(SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_CONTINUOUS || SCIPvarGetType(origvars[j]) == SCIP_VARTYPE_IMPLINT)
							(*compisinteger)[(*generatorsize)-1] = FALSE;
					}
				}
			}
		}
	}

	origvars = GCGmasterVarGetOrigvars(mastervar);
	norigvars = GCGmasterVarGetNOrigvars(mastervar);
	origvals = GCGmasterVarGetOrigvals(mastervar);

	for(i=0; i<norigvars; ++i)
	{
		for(j=0; j<*generatorsize; ++j)
		{
			if(origvarsunion[j]==origvars[i])
			{
				if(SCIPvarGetType(origvars[i]) == SCIP_VARTYPE_CONTINUOUS)
					(*compisinteger)[j] = FALSE;				 
				if(!SCIPisZero(scip, origvals[i]))
					(*generator)[j] = origvals[i];
				break;
			}
		}
	}

	SCIPfreeMemoryArray(scip, &origvarsunion);

	return SCIP_OKAY;
}

/** method for calculating the median over all fractional components values if its the minimum return ceil(arithm middle)*/
static
SCIP_Real GetMedian(SCIP* scip, SCIP_Real* array, int arraysize, SCIP_Real min)
{  
	SCIP_Real Median;
	SCIP_Real swap;
	int l;
	int r;
	int i;
	int j;
	int MedianIndex;
	SCIP_Real arithmMiddle;

	r = arraysize -1;
	l = 0;
	arithmMiddle = 0;

	if( arraysize & 1)
		MedianIndex = arraysize/2;
	else
		MedianIndex = arraysize/2 -1;

	while(l < r-1)
	{
		Median = array[ MedianIndex ];
		i = l;
		j = r;
		do
		{
			while( SCIPisLT(scip, array[i], Median) ) 
				++i;
			while( SCIPisGT(scip, array[j], Median) )
				--j;
			if( i <=j )
			{
				swap = array[i];
				array[i] = array[j];
				array[j] = swap;
				++i;
				--j;
			}
		} while( i <=j );
		if( j < MedianIndex )
			l = i;
		if( i > MedianIndex )
			r = j;
	}
	Median = array[ MedianIndex ];

	if(  SCIPisEQ(scip, Median, min) )
	{
		SCIPdebugMessage("arithmmiddle \n");
		for(i=0; i<arraysize; ++i)
			arithmMiddle += array[i];

		arithmMiddle /= arraysize;
		Median = SCIPceil(scip, arithmMiddle);
	}

	return Median;
}

// comparefunction for lexicographical sort
static
SCIP_DECL_SORTPTRCOMP(ptrcomp)
{
	struct GCG_Strip* strip1;
	struct GCG_Strip* strip2;
	int i;

	strip1 = (struct GCG_Strip*) elem1;
	strip2 = (struct GCG_Strip*) elem2;

	i = 0;

	assert(strip1->generatorsize == strip2->generatorsize);

	for( i=0; i< strip1->generatorsize; ++i)
	{
		if( strip1->generator[i] > strip2->generator[i] )
			return -1;
		if( strip1->generator[i] < strip2->generator[i] )
			return 1;
	}

	return 0;
}

// lexicographical sort using scipsort
// !!! changes the array
static
SCIP_RETCODE LexicographicSort( struct GCG_Strip** array, int arraysize)
{

	SCIPdebugMessage("Lexicographic sorting\n");

	SCIPsortPtr((void**)array, ptrcomp, arraysize );

	//change array
	return SCIP_OKAY;
}


// compare function for ILO: returns 1 if bd1 < bd2 else -1 
static
int ILOcomp( struct GCG_Strip* strip1, struct GCG_Strip* strip2, ComponentBoundSequence** C, int NBoundsequences, int* sequencesizes, int p) // ComponentBoundSequence* S, int Ssize, int* IndexSet, int indexsetsize)
{
	int i;
	int ivalue;
	int j;
	int k;
	int Nupper;
	int Nlower;
	SCIP_Bool returnvalue;
	ComponentBoundSequence** CopyC;
	int* newsequencesizes;

	i = -1;
	j = 0;
	k = 0;
	Nupper = 0;
	Nlower = 0;

	//lexicographic Order ?
	if( C == NULL || NBoundsequences <= 1 )
		return (*ptrcomp)( strip1, strip2);// == -1);

	assert(C != NULL);
	assert(NBoundsequences > 0);
	assert(p-1 < NBoundsequences);
	
	//find i which is in all S in C on position p (not exactly like pseudocode ?
	while( sequencesizes[k] < p-1 )
	{
		++k;
		assert(k < NBoundsequences);
	}
	i = C[k][p-1][0];
//	SCIPdebugMessage("i = %d\n", i);
	assert(i <= strip1->generatorsize);
	ivalue = C[k][p-1][2];

	assert(i>=0);

	//calculate subset of C
	for(j=0; j< NBoundsequences; ++j)
	{
		if(sequencesizes[j] >= p)
		{
			if(C[j][p-1][1] == 1)
				++Nupper;
			else 
				++Nlower;
		}
	}

	if( strip1->generator[i]>=ivalue && strip2->generator[i]>=ivalue )
	{
		k=0;
		SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nupper) );
		SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nupper) );
		for(j=0; j< NBoundsequences; ++j)
		{
			if( sequencesizes[j] >= p )
				assert(C[j][p-1][0] == i);

			if(sequencesizes[j] >= p && C[j][p-1][1] == 1)
			{
				CopyC[k] = C[j]; 
				newsequencesizes[k] = sequencesizes[j]; 
				++k;
			}
		}

		assert( k <= Nupper+1 );

		returnvalue = ILOcomp( strip1, strip2, CopyC, Nupper, newsequencesizes, p+1); 

		return returnvalue;
	}


	if( strip1->generator[i]<ivalue && strip2->generator[i]<ivalue )
	{
		k=0;
		SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Nlower) );
		SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Nlower) );
		for(j=0; j< NBoundsequences; ++j)
		{
			if( sequencesizes[j] >= p )
				assert(C[j][p-1][0] == i);

			if( sequencesizes[j]>=p && C[j][p-1][1] == 0)
			{
				CopyC[k] = C[j];    
				newsequencesizes[k] = sequencesizes[j];  
				++k;
			}
		}

		assert( k <= Nlower+1);

		returnvalue = ILOcomp( strip1, strip2, CopyC, Nlower, newsequencesizes, p+1);// S, Ssize, IndexSet, indexsetsize, p+1);

		SCIPfreeMemoryArray(scip, &newsequencesizes);
		SCIPfreeMemoryArray(scip, &CopyC);

		return returnvalue; 
	}
	if( strip1->generator[i] > strip2->generator[i])
		return 1;
	else 
		return -1;
}

// comparefunction for induced lexicographical sort
static
SCIP_DECL_SORTPTRCOMP(ptrilocomp)
{
	struct GCG_Strip* strip1;
	struct GCG_Strip* strip2;
	int returnvalue;

	strip1 = (struct GCG_Strip*) elem1;
	strip2 = (struct GCG_Strip*) elem2;

	returnvalue = ILOcomp( strip1, strip2, strip1->C, strip1->Csize, strip1->sequencesizes, 1); //NULL, 0, strip1->IndexSet, strip1->generatorsize, 1);

	return returnvalue;
}

// induced lexicographical sort
static
SCIP_RETCODE InducedLexicographicSort( SCIP* scip, struct GCG_Strip*** array, int arraysize, ComponentBoundSequence** C, int NBoundsequences, int* sequencesizes )
{
	int i;
	
	SCIPdebugMessage("Induced Lexicographic sorting\n");

	if( NBoundsequences == 0 )
		return LexicographicSort( array, arraysize );
	assert( C!= NULL );

<<<<<<< HEAD
	//ILOQSort( scip, array, arraysize, C, sequencesize, 0, arraysize-1 );

	//set data in the strips for the ptrcomp
	//n=array[0]->generatorsize;
	/*
	SCIP_CALL( SCIPallocMemoryArray(scip, &IS, n) );
	for( i=0; i<n; ++i )
		IS[i]=i;
	 */

	for( i=0; i<arraysize; ++i ){
<<<<<<< HEAD
		(*array)[i]->C = C;
		(*array)[i]->Csize = NBoundsequences; 
		(*array)[i]->sequencesizes = sequencesizes;
=======
		array[i]->C = C;
		array[i]->Csize = NBoundsequences;
		array[i]->sequencesizes = sequencesizes;
>>>>>>> ab47181035a7ef9f6e168d2b5114ac46e38f030d
		//array[i]->strip->IndexSet = IS;
		//array[i]->strip->Indexsetsize = n;
	}

<<<<<<< HEAD
	SCIPsortPtr( array, ptrilocomp, arraysize);
=======
=======
	assert(arraysize > 0);
	
	if(arraysize <= 1)
		return SCIP_OKAY;
	
	//SCIPdebugMessage("set Sortdata\n");
	for( i=0; i<arraysize; ++i )
	{
		array[i]->Csize = NBoundsequences;
		array[i]->sequencesizes = sequencesizes;
		array[i]->C = C;
	}
<<<<<<< HEAD
	SCIPdebugMessage("start sorting\n");
>>>>>>> 8b02942... clean up
=======
	//SCIPdebugMessage("start sorting\n");
>>>>>>> dfa55e8... empty mastercons in root branchdata
	SCIPsortPtr((void**)array, ptrilocomp, arraysize);
>>>>>>> ab47181035a7ef9f6e168d2b5114ac46e38f030d

	//free memory
	//SCIPdebugMessage("free Sortdata\n");

	return SCIP_OKAY;
}


// separation at the root node
static
SCIP_RETCODE Separate( SCIP* scip, struct GCG_Strip** F, int Fsize, int* IndexSet, int IndexSetSize, ComponentBoundSequence* S, int Ssize, struct GCG_Record** record )
{
	int i;
	int j;
	//int n;
	int k;
	int l;
	int Jsize;
	int* J;
	SCIP_Real median;
	SCIP_Real min;
	int Fupper;
	int Flower;
	int* priority;
	struct GCG_Strip** copyF;
	ComponentBoundSequence* upperLowerS;
	//ComponentBoundSequence* lowerS;
	SCIP_Real* alpha;
	SCIP_Real* compvalues;
	SCIP_Real  muF;
	SCIP_Real maxPriority;
	SCIP_Bool found;
	//	SCIP_VAR* origvar;
	//	unsigned int seed;

	//	seed = 0;
	i = 0;
	j = 0;
	k = 0;
	l = 0;
	Jsize = 0;
	Fupper = 0;
	Flower = 0;
	muF = 0;
	//copySsize = 0;
	min = INT_MAX;
	maxPriority = INT_MIN;
	found = FALSE;

	SCIPdebugMessage("Separate\n");

	if(Fsize == 0 || IndexSetSize == 0)
		return SCIP_OKAY;
	
	SCIPdebugMessage("Fsize = %d; Ssize = %d, IndesSetSize = %d\n", Fsize, Ssize, IndexSetSize);

	assert( F != NULL ); 
	assert( IndexSet != NULL );

	for(j=0; j<Fsize; ++j)
		muF += F[j]->mastervarValue;
	SCIP_CALL( SCIPallocMemoryArray(scip, &alpha, IndexSetSize) );

	for(k=0; k<IndexSetSize; ++k)
	{
		ComponentBoundSequence* copyS;
		SCIP_Real mu_F;
		SCIP_Bool even;

		even = TRUE;
		mu_F = 0;
		i = IndexSet[k]; 
		alpha[k] = 0;

		if(!F[0]->compisinteger[i])
			continue;

		for(j=0; j<Fsize; ++j)
		{
			alpha[k] += F[j]->generator[i] * F[j]->mastervarValue;
			if(F[j]->generator[i] != 0 && F[j]->mastervarValue != 0)
			{
				SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
				SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
			}
		}
		if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
		{
			++Jsize;
			SCIPdebugMessage("alpha[%d] = %g\n", k, alpha[k]);
		}
		if( SCIPisGT(scip, alpha[k] - SCIPfloor(scip, alpha[k]), 0) )
		{
			found = TRUE;

			SCIPdebugMessage("found a fractional alpha\n");

			//add to record
			++Ssize;
			SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, Ssize) );
			for(l=0; l < Ssize-1; ++l)
				for(j=0; j<3; ++j)
					copyS[l][j] = S[l][j];


			copyS[Ssize-1][0] = i;
			copyS[Ssize-1][1] = 1;

			SCIP_CALL( SCIPallocMemoryArray(scip, &compvalues, Fsize) );
			for(l=0; l<Fsize; ++l)
			{
				compvalues[l] = F[l]->generator[i];
				if( SCIPisLT(scip, compvalues[l], min) )
					min = compvalues[l];
			}
			median = GetMedian(scip, compvalues, Fsize, min);
			SCIPdebugMessage("median = %g \n", median);
			SCIPdebugMessage("compisinteger = %d \n", F[0]->compisinteger[i]);
			SCIPdebugMessage("TRUE = %d \n", TRUE);
			SCIPdebugMessage("FALSE = %d \n", FALSE);
			j = 0;
			do
			{
				mu_F = 0;
				if(even)
				{
					median = median+j;
					even = FALSE;
				}
				else 
				{
					median = median-j;
					even = TRUE;
				}

				for(l=0; l<Fsize; ++l)
				{
					if( SCIPisGE(scip, F[l]->generator[i], median) ) //F[l]->generator[i] >= median )
						mu_F += F[l]->mastervarValue;
				}
				++j;

			}while( SCIPisEQ(scip, mu_F - SCIPfloor(scip, mu_F), 0) );

			SCIPdebugMessage("mu_F = %g \n", mu_F);
			SCIPdebugMessage("newmedian = %g \n", median);
			SCIPdebugMessage("last component = %d \n", i);

			copyS[Ssize-1][2] = median;

			SCIPfreeMemoryArray(scip, &compvalues);

			(*record)->recordsize++;
			if((*record)->recordsize == 1 )
			{
				SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
				SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
			}
			else
			{
				SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
				SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
			}
			(*record)->record[(*record)->recordsize-1] = copyS;
			(*record)->sequencesizes[(*record)->recordsize-1] = Ssize;
			--Ssize;
		}
	}

	if(found)
	{
		SCIPfreeMemoryArray(scip, &alpha);
		SCIPdebugMessage("one S found with size %d\n", (*record)->sequencesizes[(*record)->recordsize-1]);
		//SCIPfreeMemoryArray(scip, &copyS);
		//	return record;
		return SCIP_OKAY;
	}

	//discriminating components
	SCIP_CALL( SCIPallocMemoryArray(scip, &J, Jsize) );
	j=0;
	for(k=0; k<IndexSetSize; ++k)
	{
		if( SCIPisGT(scip, alpha[k], 0) && SCIPisLT(scip, alpha[k], muF) )
		{
			J[j] = IndexSet[k];
			++j;
		}
	}

	//compute priority  (max-min)
	SCIP_CALL( SCIPallocMemoryArray(scip, &priority, Jsize) );
	for(j=0; j<Jsize; ++j)
	{
		int maxcomp;
		int mincomp;

		maxcomp = INT_MIN;
		mincomp = INT_MAX;

		for(l=0; l<Fsize; ++l)
		{
			if(F[l]->generator[J[j]] > maxcomp)
				maxcomp = F[l]->generator[J[j]];
			if(F[l]->generator[J[j]] < mincomp)
				mincomp = F[l]->generator[J[j]];
		}
		priority[J[j]] = maxcomp-mincomp;
	}

	//Partition
	SCIPdebugMessage("Partitioning\n");
	do{
		min = INT_MAX;
		maxPriority = INT_MIN;
	
		//max-min priority
		for(j=0; j<Jsize; ++j)
		{
			if(priority[J[j]] > maxPriority && F[0]->compisinteger[J[j]])
			{
				maxPriority = priority[J[j]];
				i = J[j];
			}
		}
		SCIP_CALL( SCIPallocMemoryArray(scip, &compvalues, Fsize) );
		for(l=0; l<Fsize; ++l)
		{
			compvalues[l] = F[l]->generator[i];
			if( SCIPisLT(scip, compvalues[l], min) )
				min = compvalues[l];
		}
		median = GetMedian(scip, compvalues, Fsize, min); 
		SCIPfreeMemoryArray(scip, &compvalues);

		if(median != 0)
		{
			SCIPdebugMessage("median = %g\n", median);
			SCIPdebugMessage("min = %g\n", min);
		}

		if( SCIPisEQ(scip, median, min) )
		{
			//here with max-min priority
			//	assert(maxPriority == 0);
			for(j=0; j<Jsize; ++j)
			{
				if( i == J[j])
				{
					J[j] = J[Jsize-1];
					break;
				}
			}
			--Jsize;

		}
		assert(Jsize>=0);
	}while( SCIPisEQ(scip, median, min) );

	++Ssize;
	SCIP_CALL( SCIPallocMemoryArray(scip, &upperLowerS, Ssize) );
	for(l=0; l < Ssize-1; ++l)
		for(j=0; j<3; ++j)
			upperLowerS[l][j] = S[l][j];

	upperLowerS[Ssize-1][0] = i;
	upperLowerS[Ssize-1][1] = 1;
	upperLowerS[Ssize-1][2] = median;

	for(k=0; k<Fsize; ++k)
	{
		if( SCIPisGE(scip, F[k]->generator[i], median) )
			++Fupper;
		else 
			--Flower;
	}

	//choose smallest partition
	if( (Flower <= Fupper && Flower > 0) || Fupper <= 0 )
	{
		SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
		j = 0;
		for(k=0; k<Fsize; ++k)
		{
			if( SCIPisLT(scip, F[k]->generator[i], median) )
				copyF[j]=F[k];
		}
		Fsize = Flower;
	}
	else
	{
		upperLowerS[Ssize-1][1] = 0;
		SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
		j = 0;
		for(k=0; k<Fsize; ++k)
		{
			if( SCIPisGE(scip, F[k]->generator[i], median) )
				copyF[j]=F[k];
		}
		Fsize = Fupper;
	}

	//record = 
	Separate( scip, copyF, Fsize, J, Jsize, upperLowerS, Ssize, record );


	SCIPfreeMemoryArray(scip, &copyF);
	SCIPfreeMemoryArray(scip, &upperLowerS);
	SCIPfreeMemoryArray(scip, &priority);
	SCIPfreeMemoryArray(scip, &J);
	SCIPfreeMemoryArray(scip, &alpha);

	//	return record;
	return SCIP_OKAY;
}

// choose a component bound sequence 
static
SCIP_RETCODE ChoseS( SCIP* scip, struct GCG_Record* record, ComponentBoundSequence** S, int* Ssize )
{
	int minSizeOfMaxPriority;  //neede if the last comp priority is euqal to the one in other bound sequences
	int maxPriority;
	int i;
	int j;
	int Index;

	minSizeOfMaxPriority = INT_MAX;
	maxPriority = INT_MIN;
	i = 0;
	j = 0;
	Index = -1;

	SCIPdebugMessage("Chose S \n");

	assert(record->recordsize > 0);

	for( i=0; i< record->recordsize; ++i )
	{
		assert(record->sequencesizes != NULL );
		assert(record->sequencesizes[i] > 0);
		if(maxPriority <= 1 || maxPriority == INT_MIN) // later by pseudocosts e.g.
		{
			if( maxPriority < 1 || maxPriority == INT_MIN) 
			{
				maxPriority = 1; // only choose here first smallest S
				minSizeOfMaxPriority = record->sequencesizes[i];
				Index = i;
			}
			else
				if( record->sequencesizes[i] < minSizeOfMaxPriority )
				{
					minSizeOfMaxPriority = record->sequencesizes[i];
					Index = i;
				}
		}
	}
	assert(maxPriority != INT_MIN);
	assert(minSizeOfMaxPriority != INT_MAX);
	assert(Index >= 0);

	*Ssize = minSizeOfMaxPriority;
	SCIP_CALL( SCIPallocMemoryArray(scip, S, *Ssize) );
	for(i=0; i< *Ssize;++i)
		for(j=0; j<3; ++j)
			(*S)[i][j] = record->record[Index][i][j];

	assert(S!=NULL);
	assert(*S!=NULL);
	//free record
	for( i=0; i< record->recordsize; ++i )
	{
		SCIPfreeMemoryArray(scip, &(record->record[i]) );
	}
	SCIPfreeMemoryArray(scip, &(record->record) );

	SCIPdebugMessage("with size %d \n", *Ssize);

	assert(*S!=NULL);

	return SCIP_OKAY;	
}



// separation at a node other than the root node
static
SCIP_RETCODE Explore( SCIP* scip, ComponentBoundSequence** C, int Csize, int* sequencesizes, int p, struct GCG_Strip** F, int Fsize, int* IndexSet, int IndexSetSize, ComponentBoundSequence** S, int* Ssize, struct GCG_Record** record )
{
	int i;
	int j;
	int k;
	int l;
	SCIP_Real ivalue;
	int isense;
	SCIP_Real median;
	int Fupper;
	int Flower;
	int Cupper;
	int Clower;
	struct GCG_Strip** copyF;
	ComponentBoundSequence* copyS;
	ComponentBoundSequence** CopyC;
	int* newsequencesizes;
	SCIP_Real alpha_i;
	SCIP_Real  muF;
	SCIP_Bool found;
	SCIP_Real mu_F;

	i = 0;
	j = 0;
	k = 0;
	l = 0;
	Fupper = 0;
	Flower = 0;
	Cupper = 0;
	Clower = 0;
	CopyC = NULL;
	newsequencesizes = NULL;
	muF = 0;
	found = FALSE;

	SCIPdebugMessage("Explore\n");
	
	SCIPdebugMessage("with Fsize = %d, Csize = %d, Ssize = %d\n", Fsize, Csize, *Ssize);

	//call separate?
	if( C == NULL || Fsize==0 || IndexSetSize==0 || Csize == 0)
	{
		SCIPdebugMessage("go to Separate\n");
		return Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );
	}

	assert( C!=NULL );
	assert( Csize>0 );
	assert( F != NULL ); 
	assert( IndexSet != NULL );
	assert( sequencesizes != NULL );

	//find i which is in all S in C on position p 
	while( sequencesizes[k] < p )
	{
		SCIPdebugMessage("sequencesizes[%d] = %d\n", k, sequencesizes[k]);
		++k;
		if( k>= Csize )
		{
			SCIPdebugMessage("no %dth element bounded\n", p);
			return Separate( scip, F, Fsize, IndexSet, IndexSetSize, *S, *Ssize, record );
		}
		assert( k < Csize );
	}
	i = C[k][p-1][0];
	isense = C[k][p-1][1];
	ivalue = C[k][p-1][2];
	
	SCIPdebugMessage("i = %d; ivalue = %g\n", i, ivalue);

	for(j=0; j<Fsize; ++j)
		muF += F[j]->generator[i] * F[j]->mastervarValue;

	SCIPdebugMessage("muF = %g\n", muF);

	alpha_i = 0;
	for(j=0; j<Fsize; ++j)
	{
		if(isense == 1)
		{
			if( SCIPisGE(scip, F[j]->generator[i], ivalue) ) //F[j]->generator[i] >= ivalue)
			{
				alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
				if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
				{
					SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
					SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
				}
			}
		}
		else
		{
			if( SCIPisLT(scip, F[j]->generator[i], ivalue) ) //F[j]->generator[i] < ivalue)
			{
				alpha_i += F[j]->generator[i] * F[j]->mastervarValue;
				if(F[j]->mastervarValue != 0 && F[j]->generator[i] != 0)
				{
					SCIPdebugMessage("generator[%d] = %g\n", i, F[j]->generator[i]);
					SCIPdebugMessage("mastervarvalue = %g\n", F[j]->mastervarValue);
				}
			}

		}
	}

	if( SCIPisGT(scip, alpha_i - SCIPfloor(scip, alpha_i), 0) )
	{
		found = TRUE;
		SCIPdebugMessage("alpha[%d] = %g\n", i, alpha_i);

		//add to record
		++(*Ssize);
		SCIP_CALL( SCIPallocMemoryArray(scip, &copyS, *Ssize) );
		for(l=0; l < *Ssize-1; ++l)
			for(j=0; j<3; ++j)
				copyS[l][j] = (*S)[l][j];

		copyS[*Ssize-1][0] = i;
		copyS[*Ssize-1][1] = isense;

		/*
		j = 0;
		do
		{
			mu_F = 0;
			ivalue += j;

			for(l=0; l<Fsize; ++l)
			{
				if( SCIPisGE(scip, F[l]->generator[i], ivalue) ) 
					mu_F += F[l]->mastervarValue;
			}
			++j;
		}while( SCIPisEQ(scip, mu_F - SCIPfloor(scip, mu_F), 0) );
        */
		copyS[*Ssize-1][2] = ivalue;
		
		mu_F = 0;
		for(l=0; l<Fsize; ++l)
		{
			if(isense ==1 )
			{
				if( SCIPisGE(scip, F[l]->generator[i], ivalue) ) 
					mu_F += F[l]->mastervarValue;
			}
			else 
			{
				if( SCIPisLT(scip, F[l]->generator[i], ivalue) ) 
					mu_F += F[l]->mastervarValue;
			}
		}

		if(SCIPisGT(scip, mu_F - SCIPfloor(scip, mu_F), 0))
		{
			(*record)->recordsize++;
			if((*record)->recordsize == 1 )
			{
				SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
				SCIP_CALL( SCIPallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
			}
			else
			{
				SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->record), (*record)->recordsize) );
				SCIP_CALL( SCIPreallocMemoryArray(scip, &((*record)->sequencesizes), (*record)->recordsize) );
			}
			(*record)->record[(*record)->recordsize-1] = copyS;
			(*record)->sequencesizes[(*record)->recordsize-1] = *Ssize;
			--(*Ssize);
		}
		else
		{
			found = FALSE;
		}
	}

	if(found)
	{
		SCIPdebugMessage("found fractional alpha\n");
		return SCIP_OKAY;
	}

	++(*Ssize);
	if(S==NULL || *S==NULL)
	{
		assert(*Ssize == 1);
		SCIP_CALL( SCIPallocMemoryArray(scip, S, *Ssize) );
	}
	else
	{
		assert(*Ssize >1);
		SCIP_CALL( SCIPreallocMemoryArray(scip, S, *Ssize) );
	}
	median = ivalue;
	(*S)[*Ssize-1][0] = i;
	(*S)[*Ssize-1][1] = 1;
	(*S)[*Ssize-1][2] = median;

	for(k=0; k<Fsize; ++k)
	{
		if( SCIPisGE(scip, F[k]->generator[i], median) )
			++Fupper;
		else 
			++Flower;
	}

	//calculate subset of C
	for(j=0; j< Csize; ++j)
	{
		if(sequencesizes[j] >= p)
		{
			if(C[j][p-1][0] == 1)
				++Cupper;
			else 
				++Clower;
		}
	}

	if( SCIPisLE(scip, alpha_i, 0) )
		Fupper = INT_MAX;
	if( SCIPisEQ(scip, alpha_i, muF) )
		Flower = INT_MAX;

	//choose smallest partition
	if( ((Fupper <= Flower && Fupper > 0 ) || Flower <= 0) && Fupper != INT_MAX )
	{
		SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Fupper) );
		j = 0;
		for(k=0; k<Fsize; ++k)
		{
			if( SCIPisGE(scip, F[k]->generator[i], median) )
			{
				copyF[j] = F[k];
				++j;
			}
		}
		Fsize = Fupper;

		//new C
		k=0;
		for(j=0; j< Csize; ++j)
		{
			assert(C[j][p-1][0] == i);

			if(sequencesizes[j] >= p && C[j][p-1][1] == 1)
			{
				C[k] = C[j];   // copyC[k] = C[j];
				sequencesizes[k] = sequencesizes[j];   // newsequencesizes[k] = sequencesizes[j];
				++k;
			}
		}
		Csize = Cupper;
		SCIPdebugMessage("chose upper bound Fupper = %d\n", Fupper);

	}
	else
	{
		SCIPdebugMessage("chose lower bound Flower = %d\n", Flower);
		(*S)[*Ssize-1][1] = 0;
		SCIP_CALL( SCIPallocMemoryArray(scip, &copyF, Flower) );
		j = 0;
		for(k=0; k<Fsize; ++k)
		{
			if( SCIPisLT(scip, F[k]->generator[i], median) )
			{
				copyF[j] = F[k];
				++j;
			}
		}
		Fsize = Flower;

		//new C
		k=0;
		SCIP_CALL( SCIPallocMemoryArray(scip, &CopyC, Clower) );
		SCIP_CALL( SCIPallocMemoryArray(scip, &newsequencesizes, Clower) );
		for(j=0; j< Csize; ++j)
		{
			assert(C[j][p-1][0] == i);

			if(sequencesizes[j] >= p && C[j][p-1][1] == 0)
			{
				CopyC[k] = C[j]; 
				newsequencesizes[k] = sequencesizes[j];   ;
				++k;
			}
		}
		Csize = Clower;
	}
	assert( k <= Csize+1 );
	
	Explore( scip, CopyC, Csize, sequencesizes, p+1, copyF, Fsize, IndexSet, IndexSetSize, S, Ssize, record );

	SCIPfreeMemoryArray(scip, &copyF);
	if(CopyC != NULL)
		SCIPfreeMemoryArray(scip, &CopyC);
	if(newsequencesizes != NULL)
		SCIPfreeMemoryArray(scip, &newsequencesizes);

	return SCIP_OKAY;
}

// callup method for seperate 
static
SCIP_RETCODE CallSeparate( SCIP* scip, struct GCG_Strip** F, int Fsize, ComponentBoundSequence** S, int* Ssize, ComponentBoundSequence** C, int Csize, int* CompSizes )
{
	int i;
	//int n;
	int* IndexSet;
	int IndexSetSize;
	GCG_RECORD* record;
	int exploreSsize;
	ComponentBoundSequence* exploreS;

	assert(Fsize > 0);
	assert(F!=NULL);
	exploreSsize = 0;
	exploreS = NULL;

	SCIPdebugMessage("Calling Separate\n");


	//record = (struct GCG_Record*) malloc(sizeof(struct GCG_Record));
	SCIP_CALL( SCIPallocBuffer(scip, &record) );
	record->recordsize = 0;

	//calculate IndexSet
	IndexSetSize = F[0]->generatorsize;
	assert(IndexSetSize > 0);
	SCIP_CALL( SCIPallocMemoryArray(scip, &IndexSet, IndexSetSize) );
	for( i=0; i<IndexSetSize; ++i )
		IndexSet[i] = i;

	//rootnode?
	if( Csize<=0 )
		Separate( scip, F, Fsize, IndexSet, IndexSetSize, NULL, 0, &record );
	else
	{
		assert( C!=NULL );
		Explore( scip, C, Csize, CompSizes, 1, F, Fsize, IndexSet, IndexSetSize, &exploreS, &exploreSsize, &record);
		if(exploreS != NULL)
			SCIPfreeMemoryArray(scip, &exploreS);
	}	

	ChoseS( scip, record, S, Ssize );
	assert(*S!=NULL);

#ifdef SCIP_DEBUG
	for( i=0; i<*Ssize; ++i)
	{
		SCIPdebugMessage("S[%d][0] = %g\n", i, (*S)[i][0]);
		SCIPdebugMessage("S[%d][1] = %g\n", i, (*S)[i][1]);
		SCIPdebugMessage("S[%d][2] = %g\n", i, (*S)[i][2]);
	}
#endif
	SCIPfreeMemoryArray(scip, &IndexSet);
	SCIPfreeBuffer(scip, &record);

	return SCIP_OKAY;	
}

/** callback deletion method for branching data*/
static
GCG_DECL_BRANCHDATADELETE(branchDataDeleteGeneric)
{
	assert(scip != NULL);
	assert(branchdata != NULL);

	SCIPdebugMessage("branchDataDeleteGeneric:\n");

	if(*branchdata == NULL)
	{
	   SCIPdebugMessage("branchDataDeleteGeneric: cannot delte empty branchdata\n");

	   return SCIP_OKAY;
	}

	if((*branchdata)->mastercons != NULL)
	   SCIPdebugMessage("branchDataDeleteGeneric: child blocknr %d, %s\n", (*branchdata)->blocknr,
			SCIPconsGetName((*branchdata)->mastercons) );
	else
	   SCIPdebugMessage("branchDataDeleteGeneric: child blocknr %d, empty mastercons\n", (*branchdata)->blocknr);

	/* release constraint that enforces the branching decision */
	if( (*branchdata)->mastercons != NULL )
	{
		SCIP_CALL( SCIPreleaseCons(GCGrelaxGetMasterprob(scip), &(*branchdata)->mastercons) );
		(*branchdata)->mastercons = NULL;
	}

	if((*branchdata)->S != NULL && (*branchdata)->Ssize > 0)
	{
		SCIPfreeMemoryArray(scip, &((*branchdata)->S));
		(*branchdata)->S = NULL;
	}

	if( (*branchdata)->nchildNodes > 0)
	{
		if((*branchdata)->childlhs != NULL)
		{
			SCIPfreeMemoryArray(scip, &((*branchdata)->childlhs));
			(*branchdata)->childlhs = NULL;
		}
		if((*branchdata)->childbranchdatas!= NULL) // this should not happen
		{
			SCIPfreeMemoryArray(scip, &((*branchdata)->childbranchdatas));
			(*branchdata)->childbranchdatas = NULL;
		}
	}
	/*
	if((*branchdata)->C != NULL && (*branchdata)->Csize > 0)
	{
		SCIPfreeMemoryArray(scip, &((*branchdata)->C));
		(*branchdata)->C = NULL;
	}
	if((*branchdata)->sequencesizes != NULL && (*branchdata)->Csize > 0 )
	{
		SCIPfreeMemoryArray(scip, &((*branchdata)->sequencesizes));
		(*branchdata)->sequencesizes = NULL;
	}
	*/

	SCIPfreeMemory(scip, branchdata);
	*branchdata = NULL;

	return SCIP_OKAY;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_Bool pruneChildNodeByDominanceGeneric(
		SCIP*                 scip,               /**< SCIP data structure */
		SCIP_Real             lhs,                /**< lhs for childnode which is checkes to be pruned */
		ComponentBoundSequence* childS,           /**< Component Bound Sequence defining the childnode */
		int                   childSsize,
		SCIP_CONS*            masterbranch,
		int                   childBlocknr             /**< number of the block for the childnode */
)
{
	SCIP_CONS* cons;
	int i;

	i = 0;
	SCIPdebugMessage("Prune by dominance\n");
	cons = GCGconsMasterbranchGetParentcons(masterbranch);

	if( cons == NULL)
	{
		SCIPdebugMessage("cons == NULL\n");
		return FALSE;
	}
	while( GCGconsMasterbranchGetParentcons(cons) != NULL )
	{
		GCG_BRANCHDATA* parentdata;

		cons = GCGconsMasterbranchGetParentcons(cons);  //skip the first ancestor node
		parentdata = GCGconsMasterbranchGetBranchdata(cons);

		if(parentdata->blocknr == childBlocknr && parentdata->Ssize >= childSsize)
		{
			for( i=0; i<childSsize+1; ++i)
			{
				if(parentdata->childlhs[i] != lhs)
					continue;
				if( parentdata->S[i][0] == childS[i][0] && parentdata->S[i][2] == childS[i][2] )
				{
					if(parentdata->S[i][1] != childS[i][1] && i < childSsize-1 )
						break; //subset = FALSE;
					if( i == childSsize-1 )
					{
						if(childSsize == parentdata->Ssize)
						{
							SCIPdebugMessage("child pruned\n");
							return TRUE;
						}
						if( parentdata->S[i][1] != childS[i][1] )  //child is induced by parantdata->S
						{
							SCIPdebugMessage("child pruned\n");
							return TRUE;
						}
					}
				}
				else 
					break; //subset = FALSE; no subset or other lhs
			}
		}
	}

	SCIPdebugMessage("child not pruned\n");
	return FALSE;
}

/** for given component bound sequence S, create |S|+1 Vanderbeck branching nodes */
static
SCIP_RETCODE createChildNodesGeneric(
		SCIP*                 scip,               /**< SCIP data structure */
		SCIP_BRANCHRULE*      branchrule,         /**< branching rule */
		ComponentBoundSequence* S,              /**< Component Bound Sequence defining the nodes */
		int                   Ssize,
		int                   blocknr,             /**< number of the block */
		struct GCG_Strip**           F,                   /**< strips with mu>0 */  //for rhs, will be small than
		int                   Fsize,
		SCIP_CONS*            parentcons,
		int*                  nmasternodes,
		SCIP_Bool             createorignodes
)
{
	SCIP*  masterscip;
	int i;
	int p;
	int k;
	SCIP_Real pL;
	SCIP_Real L;
	SCIP_Real lhs;
	int nmastervars;
	int nmastervars2;
	int ncopymastervars;
	int nbranchcands;
	int allnorigvars;
	SCIP_Real mu;  // mu(S)
	SCIP_VAR** mastervars;
	SCIP_VAR** mastervars2;
	SCIP_VAR** branchcands;
	SCIP_VAR** copymastervars;
	SCIP_VAR** allorigvars;
	GCG_BRANCHDATA* parentdata;

	SCIPdebugMessage("Vanderbeck branching rule Node creation for blocknr %d with %d identical blocks \n", blocknr, GCGrelaxGetNIdenticalBlocks(scip, blocknr));

	lhs = 0;
	p = 0;
	k = 0;
	i = 0;
	L = 0;
	pL = GCGrelaxGetNIdenticalBlocks(scip, blocknr);
	mu = 0;
	*nmasternodes = 0;
	if( parentcons != NULL)
	{
		parentdata = GCGconsMasterbranchGetBranchdata(parentcons);
		if(parentdata == NULL)
		{
			SCIP_CALL( SCIPallocMemory(scip, &parentdata) );
			parentdata->Ssize = 0;
		}
		parentdata->nchildNodes = 0;
		if(parentdata->S == NULL || parentdata->Ssize == 0)
		{
			parentdata->S = S;
			parentdata->Ssize = Ssize;
			parentdata->childlhs = NULL;
			parentdata->C = NULL;
			parentdata->sequencesizes = NULL;
		}
	}
	else
		parentdata = NULL;

	SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

	// get variable data of the master problem
	masterscip = GCGrelaxGetMasterprob(scip);
	SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
	nmastervars2 = nmastervars;
	assert(nmastervars >= 0); 
	SCIP_CALL( SCIPallocMemoryArray(scip, &copymastervars, nmastervars) );
	SCIP_CALL( SCIPallocMemoryArray(scip, &mastervars2, nmastervars) );

	for(i=0; i<nmastervars; ++i)
	{
		mastervars2[i] = mastervars[i];
		copymastervars[i] = mastervars[i];
	}

	assert(scip != NULL);
	assert(Ssize > 0);
	assert(S != NULL);
	assert(F != NULL);
	assert(Fsize > 0);

	SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

	SCIPdebugMessage("Vanderbeck branching rule: creating %d nodes\n", Ssize+1);

	for( p=0; p<Ssize+1; ++p )
	{
		SCIP_NODE* origchild;
		SCIP_CONS* origcons;
		GCG_BRANCHDATA* branchchilddata;
		char childname[SCIP_MAXSTRLEN];

		mu = 0;

		/* allocate branchdata for same child and store information */
		SCIP_CALL( SCIPallocMemory(scip, &branchchilddata) );
		branchchilddata->blocknr = blocknr;
		branchchilddata->mastercons = NULL;
		branchchilddata->childnr = p;
		branchchilddata->nchildNodes = 0;
		SCIPdebugMessage("p = %d \n", p);
		if( p == Ssize )
		{
			SCIP_CALL( SCIPallocMemoryArray(scip, &(branchchilddata->S), Ssize) );
			branchchilddata->Ssize = Ssize;
		}
		else
		{
			SCIP_CALL( SCIPallocMemoryArray(scip, &(branchchilddata->S), p+1) );
			branchchilddata->Ssize = p+1;
		}
		for( k=0; k<=p; ++k)
		{
			ComponentBoundSequence compBound;

			if( k == Ssize )
			{
				assert( p == Ssize );
				compBound[0] = S[k-1][0]; //k <-> p ?
				compBound[2] = S[k-1][2];
				compBound[1] = S[k-1][1];
				branchchilddata->S[k-1][0] = compBound[0];
				branchchilddata->S[k-1][1] = compBound[1];
				branchchilddata->S[k-1][2] = compBound[2];
			}
			else
			{
				if( k < p-1 )
				{
					compBound[0] = S[k][0];  //k <-> p ?
					compBound[1] = S[k][1];
					compBound[2] = S[k][2];
				}
				else
				{
					compBound[0] = S[k][0]; //k <-> p ?
					compBound[2] = S[k][2];
					if( S[p][1] == 1)
						compBound[1] = 0;
					else
						compBound[1] = 1;
				}
				branchchilddata->S[k][0] = compBound[0];
				branchchilddata->S[k][1] = compBound[1];
				branchchilddata->S[k][2] = compBound[2];
			}
		}
		SCIPdebugMessage("branchchilddata set\n");
		//last node?
		if( p == Ssize )
		{
			lhs = pL;
		}
		else
		{
			//calculate mu
			for(k=0;k<=p;++k)
			{
				ncopymastervars = nmastervars2;
				for( i=0; i<ncopymastervars; ++i)
				{
					SCIP_VAR* swap;
					SCIP_Real generator_i;
					SCIP_Real* generator;
					int generatorsize;
					SCIP_Bool* compisinteger;
					
					if(GCGvarGetBlock(mastervars2[i]) == blocknr)
					{
						getGenerators(scip, &generator, &generatorsize, &compisinteger, blocknr, mastervars, nmastervars, mastervars2[i]);
						generator_i = generator[ (int) SCIPceil(scip, S[k][0]-0.5)];

						//						SCIPdebugMessage("generator_i = %g \n", generator_i);
						//						SCIPdebugMessage("S[k][0] = %d\n", S[k][0]);
						//						SCIPdebugMessage("S[k][2] = %d\n", S[k][2]);
						//						SCIPdebugMessage("mastvervalue = %g\n", SCIPgetSolVal(masterscip, NULL, mastervars2[i]));

						if(S[k][1] == 1 )
						{
							if( SCIPisGE(scip, generator_i, S[k][2]) ) 
							{
								if( k == p )
									mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
							}
							else
							{
								if(ncopymastervars > 0)
								{
									swap = mastervars2[i];
									mastervars2[i] = mastervars2[ncopymastervars-1];
									mastervars2[ncopymastervars-1] = swap;
									--ncopymastervars;
									--i;
								}
							}
						}
						else  //nested erasing
						{
							if( SCIPisLT(scip, generator_i, S[k][2]) )
							{
								if( k == p )
									mu += SCIPgetSolVal(masterscip, NULL, mastervars2[i]);
							}
							else
							{
								if(ncopymastervars > 0)
								{
									swap = mastervars2[i];
									mastervars2[i] = mastervars2[ncopymastervars-1];
									mastervars2[ncopymastervars-1] = swap;
									--ncopymastervars;
									--i;
								}
							}
						}
					}
					else
					{
						if(ncopymastervars > 0)
						{
							swap = mastervars2[i];
							mastervars2[i] = mastervars2[ncopymastervars-1];
							mastervars2[ncopymastervars-1] = swap;
							--ncopymastervars;
							--i;
						}
					}
				}
				nmastervars2 = ncopymastervars;
			}
			if( p == Ssize-1)
				L = SCIPceil(scip, mu);
			else
			{
				L = mu;
			}
			lhs = pL-L+1;
			SCIPdebugMessage("pL-L+1 = %g \n", pL-L+1);
		}
		pL = L;

		branchchilddata->lhs = lhs;
		branchchilddata->childlhs = NULL;
		branchchilddata->C = NULL;
		branchchilddata->sequencesizes = NULL;
		SCIPdebugMessage("L = %g \n", L);
		SCIPdebugMessage("lhs set to %g \n", lhs);

		if( parentcons == NULL || !pruneChildNodeByDominanceGeneric(scip, lhs, branchchilddata->S, branchchilddata->Ssize, parentcons, blocknr) )
		{
			++(*nmasternodes);
			if( parentcons != NULL)
			{
				++parentdata->nchildNodes;
				if(parentdata->nchildNodes == 1)
				{
					SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
					SCIP_CALL( SCIPallocMemoryArray(scip, &(parentdata->childbranchdatas), parentdata->nchildNodes) );
				}
				else
				{
					SCIP_CALL( SCIPreallocMemoryArray(scip, &(parentdata->childlhs), parentdata->nchildNodes) );
					SCIP_CALL( SCIPreallocMemoryArray(scip, &(parentdata->childbranchdatas), parentdata->nchildNodes) );
				}
				//assert( p == parentdata->nchildNodes-1 ); 
				parentdata->childlhs[parentdata->nchildNodes-1] = lhs;
				parentdata->childbranchdatas[parentdata->nchildNodes-1] = branchchilddata;
			}
			if(createorignodes)
			{
				// define name for constraint 
				(void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", p, lhs);
				SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );
				
				//create cons
				SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
						GCGconsOrigbranchGetActiveCons(scip), branchrule, branchchilddata) );
				SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
				SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
			}
		}
	}

	SCIPfreeMemoryArray(scip, &mastervars2);
	SCIPfreeMemoryArray(scip, &copymastervars);

	return SCIP_OKAY; 
}

/** returns the number of successor nodes needed for branch_master while using the generic branching scheme */
int GCGbranchGenericGetNChildnodes(
		SCIP*            masterscip,
		SCIP_Bool        createorignodes
)
{
	int nmasternodes;
	SCIP* scip;
	SCIP_Bool feasible;
	SCIP_Bool SinC;
	SCIP_VAR** branchcands;
	SCIP_VAR** allorigvars;
	SCIP_VAR** mastervars;
	int nmastervars;
	SCIP_CONS* masterbranchcons;
	SCIP_CONS* parentcons;
	int nbranchcands;
	GCG_BRANCHDATA* branchdata;
	SCIP_VAR* mastervar;
	SCIP_Real mastervarValue;
	ComponentBoundSequence* S;
	ComponentBoundSequence** C;
	struct GCG_Strip** F;
	int Ssize;
	int Csize;
	int Fsize;
	int* sequencesizes;
	int blocknr;
	int i;
	int c;
	int k;
	int allnorigvars;

	blocknr = -2;
	Ssize = 0;
	Csize = 0;
	Fsize = 0;
	i = 0;
	c = 0;
	feasible = TRUE;
	SinC = TRUE;
	S = NULL;
	F = NULL;
	nmasternodes = 0;

	assert(masterscip != NULL);

	SCIPdebugMessage("Get number of childnodes for Vanderbecks generic branching\n");

	if(createorignodes)
	{
		scip = masterscip;
		masterscip = GCGrelaxGetMasterprob(scip);
	}
	else
	{
		SCIPdebugMessage("get masterprob\n");
		scip = GCGpricerGetOrigprob(masterscip);
		SCIPdebugMessage("get branchcands\n");
	}
	assert(scip != NULL);
	SCIP_CALL( SCIPgetLPBranchCands(masterscip, &branchcands, NULL, NULL, &nbranchcands, NULL) );

	SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );
	SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );

	assert(nbranchcands > 0);

	for( i=0; i<nbranchcands; ++i )
	{
		mastervar = branchcands[i];
		assert(GCGvarIsMaster(mastervar));
		blocknr = GCGvarGetBlock(mastervar);
		if(blocknr >= -1)
			break;
	}
	if( blocknr < -1 )
	{
		feasible = TRUE;
		SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");

		return -1;
	}
	else
		feasible = FALSE;

	k = 0;
	masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);

	//calculate F and the strips
	SCIPdebugMessage("calculate F\n");
	for( i=0; i<nbranchcands; ++i )
	{
		mastervar = branchcands[i];
		assert(GCGvarIsMaster(mastervar));

		if(blocknr == GCGvarGetBlock(mastervar))
		{
			mastervarValue = SCIPgetSolVal(masterscip, NULL, mastervar);
			if( SCIPisGT(scip, mastervarValue - SCIPfloor(scip, mastervarValue), 0) )
			{
				struct GCG_Strip* strip;
				SCIP_CALL( SCIPallocBuffer(scip, &strip) );
				if(Fsize == 0)
				{
					++Fsize;
					SCIP_CALL( SCIPallocMemoryArray(scip, &F, Fsize) );
				}
				else
				{
					++Fsize;
					SCIP_CALL( SCIPreallocMemoryArray(scip, &F, Fsize) );
				}
				strip->blocknr = blocknr;
				strip->mastervar = mastervar;
				strip->mastervarValue = mastervarValue;

				getGenerators(scip, &(strip->generator), &(strip->generatorsize), &(strip->compisinteger), blocknr, mastervars, nmastervars, mastervar);
				F[k] = strip;
				++k;
			}
		}
	}

	//old data to regard?
	if( masterbranchcons != NULL )
	{
		SCIPdebugMessage("masterbranchcons != NULL\n");
		//calculate C
		Csize = 1;
		SCIP_CALL( SCIPallocMemoryArray(scip, &C, Csize) );
		SCIP_CALL( SCIPallocMemoryArray(scip, &sequencesizes, Csize) );

		parentcons = GCGconsMasterbranchGetParentcons(masterbranchcons);
		
		if( parentcons != NULL)
		{
			SCIPdebugMessage("parentcons != NULL\n");
			branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
			assert(branchdata != NULL);
			assert(branchdata->Ssize > 0);
			SCIP_CALL( SCIPallocMemoryArray(scip, &(C[0]), branchdata->Ssize) );
			for(i=0; i<branchdata->Ssize; ++i)
			{
				C[0][i][0] = branchdata->S[i][0];
				C[0][i][1] = branchdata->S[i][1];
				C[0][i][2] = branchdata->S[i][2];
			}
			sequencesizes[0] = branchdata->Ssize;

			parentcons = GCGconsMasterbranchGetParentcons(parentcons);
		}
		else 
			Csize = 0;
		while( parentcons != NULL )
		{
			SCIPdebugMessage("while parentcons != NULL\n");
			branchdata = GCGconsMasterbranchGetBranchdata(parentcons);
			assert(branchdata != NULL);
			if( branchdata->S == NULL || branchdata->Ssize == 0)
				break;
			//S not yet in C ?
			for( c=0; c<Csize && SinC; ++c)
			{
				//SinC = TRUE;
				if(branchdata->Ssize == sequencesizes[c])
				{
					for( i=0; i<branchdata->Ssize; ++i)
					{
						if(branchdata->S[i][0] != C[c][i][0] || branchdata->S[i][1] != C[c][i][1] || branchdata->S[i][2] != C[c][i][2] )
						{
							SinC = FALSE;
							break;
						}
					}
				}
				else
					SinC = FALSE;
			}
			if(!SinC)
			{
				SCIPdebugMessage("!SinC\n");
				++Csize;
				SCIP_CALL( SCIPreallocMemoryArray(scip, &C, Csize) );
				SCIP_CALL( SCIPreallocMemoryArray(scip, &sequencesizes, Csize) );
				SCIP_CALL( SCIPallocMemoryArray(scip, &(C[Csize-1]), branchdata->Ssize) );
				
				for(i=0; i<branchdata->Ssize; ++i)
				{
					C[Csize-1][i][0] = branchdata->S[i][0];
					C[Csize-1][i][1] = branchdata->S[i][1];
					C[Csize-1][i][2] = branchdata->S[i][2];
				}
				sequencesizes[Csize-1] = branchdata->Ssize;
			}
			parentcons = GCGconsMasterbranchGetParentcons(parentcons);
		}

		if(C != NULL){
			SCIPdebugMessage("C != NULL\n");
			SCIPdebugMessage("Csize = %d\n", Csize);
			
			for(i=0; i<Csize; ++i)
			{
				for(c=0; c<sequencesizes[i]; ++c)
				{
					SCIPdebugMessage("C[%d][%d][0] = %g\n", i, c, C[i][c][0]);
					SCIPdebugMessage("C[%d][%d][1] = %g\n", i, c, C[i][c][1]);
					SCIPdebugMessage("C[%d][%d][2] = %g\n", i, c, C[i][c][2]);
				}
			}
			
			SCIP_CALL( InducedLexicographicSort(scip, F, Fsize, C, Csize, sequencesizes) );
			SCIP_CALL( CallSeparate(scip, F, Fsize, &S, &Ssize, C, Csize, sequencesizes) );
		}
		else
		{
			SCIPdebugMessage("C == NULL\n");
			SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
			SCIP_CALL( CallSeparate( scip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
		}
		SCIPfreeMemoryArray(scip, &sequencesizes);
		SCIPfreeMemoryArray(scip, &C);
	}
	else
	{
		SCIPdebugMessage("masterbanchcons == NULL\n");
		SCIP_CALL( InducedLexicographicSort( scip, F, Fsize, NULL, 0, NULL ) );
		SCIP_CALL( CallSeparate( scip, F, Fsize, &S, &Ssize, NULL, 0, NULL ) );
	}
	assert(S!=NULL);

	if( feasible )
	{
		SCIPdebugMessage("Vanderbeck generic branching rule could not find variables to branch on!\n");
		return -1;
	}

	/* create the |S|+1 child nodes in the branch-and-bound tree */
	SCIP_CALL( createChildNodesGeneric(scip, NULL, S, Ssize, blocknr, F, Fsize, masterbranchcons, &nmasternodes, createorignodes) );

	SCIPdebugMessage("free F\n");
	for(i=0; i<Fsize; ++i)
	{
		struct GCG_Strip* strip;
		strip = F[i];
		SCIPfreeMemoryArray(scip, &(strip->compisinteger));
		SCIPfreeMemoryArray(scip, &(strip->generator));
		SCIPfreeBuffer(scip, &strip);
	}
	SCIPfreeMemoryArray(scip, &F);

	return nmasternodes;
}

/** callback activation method */
static
GCG_DECL_BRANCHACTIVEMASTER(branchActiveMasterGeneric)
{
	SCIP* masterscip;
	SCIP_VAR** mastervars;
	SCIP_VAR** copymastervars;
	SCIP_VAR** allorigvars;
	int allnorigvars;
	int nmastervars;
	int nnewmastervars;
	int i;
	int p;
	char name[SCIP_MAXSTRLEN];

	i = 0;
	p = 0;
	nmastervars = 0;

	assert(scip != NULL);
	assert(branchdata != NULL);
	assert(branchdata->S != NULL);

	masterscip = scip;//GCGrelaxGetMasterprob(scip);
	assert(masterscip != NULL);
	SCIP_CALL( SCIPgetVarsData(masterscip, &mastervars, &nmastervars, NULL, NULL, NULL, NULL) );
	SCIP_CALL( SCIPgetVarsData(scip, &allorigvars, &allnorigvars, NULL, NULL, NULL, NULL) );

	SCIP_CALL( SCIPallocMemoryArray(scip, &copymastervars, nmastervars) );

	for(i=0; i<nmastervars; ++i)
		copymastervars[i] = mastervars[i];


	SCIPdebugMessage("branchActiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->blocknr,
			branchdata->Ssize);

	/* create corresponding constraint in the master problem, if not yet created */
	if( branchdata->mastercons == NULL )
	{

		(void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "child(%d, %g)", branchdata->childnr, branchdata->lhs);

		// create constraint for child
		SCIP_CALL( SCIPcreateConsLinear(masterscip, &(branchdata->mastercons), name, 0, NULL, NULL,
				branchdata->lhs, SCIPinfinity(scip), TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE) );

		//add mastervars
		for(p=0; p< branchdata->Ssize; ++p)
		{
			nnewmastervars = nmastervars;
			for(i=0; i<nnewmastervars; ++i)
			{
				SCIP_Real* generator;
				SCIP_Bool* compisinteger;
				int generatorsize;
				SCIP_Real generator_i;
				if( GCGvarGetBlock(copymastervars[i]) == branchdata->blocknr )
				{
					getGenerators(scip, &generator, &generatorsize, &compisinteger, branchdata->blocknr, mastervars, nmastervars, copymastervars[i]);

					generator_i = generator[(int) SCIPceil(scip, branchdata->S[p][0]-0.5)];

					if( branchdata->S[p][1] == 1)
					{
						if( SCIPisGE(scip, generator_i, branchdata->S[p][2]) )
						{
							if( p == branchdata->Ssize-1 )
								//add var to constraint
								SCIP_CALL( SCIPaddCoefLinear(masterscip, branchdata->mastercons, copymastervars[i], 1.0) );
						}
						else
						{
							//small down array
							copymastervars[i] = copymastervars[nnewmastervars-1];
							--i;
							--nnewmastervars;
						}
					}
					else
					{
						if( SCIPisLT(scip, generator_i, branchdata->S[p][2]) )
						{
							if( p == branchdata->Ssize-1 )
								//add var to constraint
								SCIP_CALL( SCIPaddCoefLinear(masterscip, branchdata->mastercons, copymastervars[i], 1.0) );
						}
						else
						{
							//small down array
							copymastervars[i] = copymastervars[nnewmastervars-1];
							--i;
							--nnewmastervars;
						}
					}
				}
				else
				{
					//small down array
					copymastervars[i] = copymastervars[nnewmastervars-1];
					--i;
					--nnewmastervars;
				}
			}
			nmastervars = nnewmastervars;
		}
	}
	/* add constraint to the master problem that enforces the branching decision */
	SCIP_CALL( SCIPaddCons(masterscip, branchdata->mastercons) );

	SCIPfreeMemoryArray(scip, &copymastervars);

	return SCIP_OKAY;
}

/** callback deactivation method */
static
GCG_DECL_BRANCHDEACTIVEMASTER(branchDeactiveMasterGeneric)
{
	SCIP* masterscip;

	assert(scip != NULL);
	assert(branchdata != NULL);
	assert(branchdata->mastercons != NULL);

	masterscip = scip;//GCGrelaxGetMasterprob(scip);
	assert(masterscip != NULL);

	SCIPdebugMessage("branchDeactiveMasterGeneric: Block %d, Ssize %d)\n", branchdata->blocknr,
			branchdata->Ssize);

	/* remove constraint from the master problem that enforces the branching decision */
	assert(branchdata->mastercons != NULL);
	SCIP_CALL( SCIPdelCons(masterscip, branchdata->mastercons) );

	return SCIP_OKAY;
}



/** callback propagation method */
static
GCG_DECL_BRANCHPROPMASTER(branchPropMasterGeneric)
{
	assert(scip != NULL);
	assert(branchdata != NULL);
	assert(branchdata->mastercons != NULL);
	assert(branchdata->S != NULL);

	SCIPdebugMessage("branchPropMasterGeneric: Block %d ,Ssize %d)\n", branchdata->blocknr,
			branchdata->Ssize);

	return SCIP_OKAY;
}


/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpGeneric)
{  /*lint --e{715}*/
	SCIPdebugMessage("Execlp method of generic branching\n");

	*result = SCIP_DIDNOTRUN;

	return SCIP_OKAY;
}



/** branching execution method for relaxation solutions */
static
SCIP_DECL_BRANCHEXECEXT(branchExecextGeneric)
{  /*lint --e{715}*/
	SCIP* masterscip;
	SCIP_Bool feasible;
	SCIP_Bool discretization;
	SCIP_CONS* masterbranchcons;
	GCG_BRANCHDATA* branchdata;
	int i;
	int j;
	i = 0;
	j = 0;
	feasible = TRUE;
	branchdata = NULL;

	assert(branchrule != NULL);
	assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
	assert(scip != NULL);
	assert(result != NULL);

	SCIPdebugMessage("Execrel method of Vanderbecks generic branching\n");

	*result = SCIP_DIDNOTRUN;

	/* the branching scheme only works for the discretization approach */
	SCIP_CALL( SCIPgetBoolParam(scip, "relaxing/gcg/discretization", &discretization) );
	if( !discretization )
	{
		SCIPdebugMessage("Generic branching only for discretization approach\n");
		return SCIP_OKAY;
	}

	/* do not perform Ryan & Foster branching if we have neither a set partitioning nor a set covering structure */
	if( GCGrelaxIsMasterSetCovering(scip) || GCGrelaxIsMasterSetPartitioning(scip) )
	{
		SCIPdebugMessage("Generic branching executed on a set covering or set partitioning problem\n");
	}

	/* check whether the current original solution is integral */
#ifdef SCIP_DEBUG
	SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), TRUE, TRUE, TRUE, TRUE, &feasible) );
#else
	SCIP_CALL( SCIPcheckSol(scip, GCGrelaxGetCurrentOrigSol(scip), FALSE, TRUE, TRUE, TRUE, &feasible) );
#endif
	if( feasible )
	{
		SCIPdebugMessage("node cut off, since origsol was feasible, solval = %f\n",
				SCIPgetSolOrigObj(scip, GCGrelaxGetCurrentOrigSol(scip)));

		*result = SCIP_CUTOFF;

		return SCIP_OKAY;
	}

	masterscip = GCGrelaxGetMasterprob(scip);

	masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
	if(masterbranchcons!=NULL)
		branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

	if(branchdata!=NULL)
	{
		for(i=0; i<branchdata->nchildNodes; ++i)
		{
			SCIP_CONS* origcons;
			SCIP_NODE* origchild;
			GCG_BRANCHDATA* nodebranchdata;
			char childname[SCIP_MAXSTRLEN];
			
			//set branchdata
			SCIP_CALL( SCIPallocMemory(scip, &nodebranchdata) );
			nodebranchdata->lhs = branchdata->childbranchdatas[i]->lhs;
			nodebranchdata->nchildNodes = 0;
			nodebranchdata->blocknr = branchdata->childbranchdatas[i]->blocknr;
			nodebranchdata->mastercons = NULL;
			nodebranchdata->childnr = branchdata->childbranchdatas[i]->childnr;
			SCIP_CALL( SCIPallocMemoryArray(scip, &(nodebranchdata->S), branchdata->childbranchdatas[i]->Ssize) );
			nodebranchdata->Ssize =	branchdata->childbranchdatas[i]->Ssize;
			for(j=0; j< nodebranchdata->Ssize; ++j)
			{
				nodebranchdata->S[j][0] = branchdata->childbranchdatas[i]->S[j][0];
				nodebranchdata->S[j][1] = branchdata->childbranchdatas[i]->S[j][1];
				nodebranchdata->S[j][2] = branchdata->childbranchdatas[i]->S[j][2];
			}
			
			// define name for constraint 
			(void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", i, nodebranchdata->lhs);
			SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

			//create cons
			SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
					GCGconsOrigbranchGetActiveCons(scip), branchrule, nodebranchdata) );
			SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
			SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
		}
		if(branchdata->childbranchdatas != NULL)
		{
			SCIPfreeMemoryArray(scip, &(branchdata->childbranchdatas));
			branchdata->childbranchdatas = NULL;
		}
	}
	else
	{
		SCIPdebugMessage("Branchdata is NULL!\n");

		GCGbranchGenericGetNChildnodes(scip, TRUE);
	}

	*result = SCIP_BRANCHED;

	return SCIP_OKAY;
}

/** branching execution method for not completely fixed pseudo solutions */
static
SCIP_DECL_BRANCHEXECPS(branchExecpsGeneric)
{  /*lint --e{715}*/
	SCIP* masterscip;
	SCIP_CONS* masterbranchcons;
	GCG_BRANCHDATA* branchdata;
	int i;
	
	i = 0;

	assert(branchrule != NULL);
	assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
	assert(scip != NULL);
	assert(result != NULL);

	SCIPdebugMessage("Execps method of Vanderbecks generic branching\n");

	return SCIP_OKAY;

	*result = SCIP_DIDNOTRUN;

	masterscip = GCGrelaxGetMasterprob(scip);

	assert(branchrule != NULL);
	assert(strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0);
	assert(scip != NULL);
	assert(result != NULL);

	SCIPdebugMessage("Execps method of Vanderbecks generic branching\n");

	*result = SCIP_DIDNOTRUN;

	if( !GCGrelaxIsMasterSetCovering(scip) || !GCGrelaxIsMasterSetPartitioning(scip) )
	{
		//SCIPdebugMessage("Executing generic branching, where master is neither set covering nor set partitioning \n");
	}

	masterbranchcons = GCGconsMasterbranchGetActiveCons(masterscip);
	if(masterbranchcons!=NULL)
		branchdata = GCGconsMasterbranchGetBranchdata(masterbranchcons);

	if(branchdata!=NULL)
	{
		for(i=0; i<branchdata->nchildNodes; ++i)
		{
			SCIP_CONS* origcons;
			SCIP_NODE* origchild;
			char childname[SCIP_MAXSTRLEN];
			// define name for constraint 
			(void) SCIPsnprintf(childname, SCIP_MAXSTRLEN, "child(%d, %g)", i, branchdata->lhs);
			SCIP_CALL( SCIPcreateChild(scip, &origchild, 0.0, SCIPgetLocalTransEstimate(scip)) );

			//create cons
			SCIP_CALL( GCGcreateConsOrigbranch(scip, &origcons, childname, origchild,
					GCGconsOrigbranchGetActiveCons(scip), branchrule, branchdata->childbranchdatas[i]) );
			SCIP_CALL( SCIPaddConsNode(scip, origchild, origcons, NULL) );
			SCIP_CALL( SCIPreleaseCons(scip, &origcons) );
		}
		SCIPfreeMemoryArray(scip, &(branchdata->childbranchdatas));
	}

	*result = SCIP_BRANCHED;

	return SCIP_OKAY;
}

/** initialization method of branching rule (called after problem was transformed) */
static
SCIP_DECL_BRANCHINIT(branchInitGeneric)
{  
	assert(branchrule != NULL);

	SCIP_CALL( GCGrelaxIncludeBranchrule(scip, branchrule, branchActiveMasterGeneric, 
			branchDeactiveMasterGeneric, branchPropMasterGeneric, NULL, branchDataDeleteGeneric) );

	return SCIP_OKAY;
}

/** creates the most infeasible LP branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleGeneric( 
		SCIP*                 scip                /**< SCIP data structure */
)
{
	SCIP_BRANCHRULEDATA* branchruledata;

	/* create branching rule data */
	branchruledata = NULL;

	SCIPdebugMessage("Include method of Vanderbecks generic branching\n");

	/* include branching rule */
	SCIP_CALL( SCIPincludeBranchrule(scip, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, 
			BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchCopyGeneric,
			branchFreeGeneric, branchInitGeneric, branchExitGeneric, branchInitsolGeneric, 
			branchExitsolGeneric, branchExeclpGeneric, branchExecextGeneric, branchExecpsGeneric,
			branchruledata) );

	/* include event handler for adding generated mastervars to the branching constraints */
	SCIP_CALL( SCIPincludeEventHdlrGenericbranchvaradd(scip) );

	return SCIP_OKAY;
}



