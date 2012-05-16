/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* #define SCIP_DEBUG */
//#define CHECKCONSISTENCY
/**@file stat.c
 * @brief  Prints information about the best decomposition
 * @author Alexander Gross
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#include "scip/scip.h"
#include "stat.h"
#include "scip_misc.h"
#include "pub_decomp.h"
#include "cons_decomp.h"
#include "struct_detector.h"



SCIP_RETCODE writeDecompositionData(SCIP* scip)
{
   DECDECOMP* decomposition=NULL;
   DEC_DETECTOR* detector;
   DEC_DECTYPE type;
   const char* typeName;
   const char* detectorName;
   int i,nBlocks,nLinkingCons,nLinkingVars;
   int* nVarsInBlocks;
   int* nConsInBlocks;
   char* ausgabe;

   SCIPallocBufferArray(scip,&ausgabe,SCIP_MAXSTRLEN);

   decomposition=DECgetBestDecomp(scip);
   type=DECdecdecompGetType(decomposition);
   typeName=DECgetStrType(type);

   detector=DECdecdecompGetDetector(decomposition);
   detectorName=detector->name;

   nBlocks=DECdecdecompGetNBlocks(decomposition);

   nVarsInBlocks=DECdecdecompGetNSubscipvars(decomposition);
   nConsInBlocks=DECdecdecompGetNSubscipconss(decomposition);

   nLinkingVars=DECdecdecompGetNLinkingvars(decomposition);
   nLinkingCons=DECdecdecompGetNLinkingconss(decomposition);

   //Print
   SCIPinfoMessage(scip,NULL,"Decomposition:\n");
   SCIPsnprintf(ausgabe,SCIP_MAXSTRLEN,"Decomposition Type: %s \n",typeName);
   SCIPinfoMessage(scip,NULL,ausgabe);

   SCIPsnprintf(ausgabe,SCIP_MAXSTRLEN,"Decomposition Detector: %s\n",detectorName);
   SCIPinfoMessage(scip,NULL,ausgabe);

   SCIPinfoMessage(scip,NULL,"Number of Blocks: %d \n",nBlocks);
   SCIPinfoMessage(scip,NULL,"Number of LinkingVars: %d\n",nLinkingVars);
   SCIPinfoMessage(scip,NULL,"Number of LinkingCons: %d\n",nLinkingCons);


   for(i=0;i<nBlocks;i++)
   {
      SCIPinfoMessage(scip,NULL,"Vars in Block %d: %d\n", i ,nVarsInBlocks[i]);
      SCIPinfoMessage(scip,NULL,"Cons in Block %d: %d\n",i,nConsInBlocks[i]);
   }

   SCIPfreeBufferArray(scip, &ausgabe);
return SCIP_OKAY;
}
