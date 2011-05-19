/*
 * struct_decomp.h
 *
 *  Created on: Jun 8, 2010
 *      Author: mbergner
 */

#ifndef STRUCT_DECOMP_H_
#define STRUCT_DECOMP_H_
#include "scip/scip.h"
enum DecDecompType
{
   DEC_ARROWHEAD, DEC_STAIRCASE, DEC_DIAGONAL, DEC_BORDERED, DEC_UNKNOWN
};

typedef enum DecDecompType DECDECOMPTYPE;

struct DecDecomp
{
   int            nblocks;
/*   SCIP**         subscips; */
   SCIP_VAR***    subscipvars;
   int*           nsubscipvars;
   SCIP_CONS***   subscipconss;
   int*           nsubscipconss;
   SCIP_CONS**    linkingconss;
   int            nlinkingconss;
   SCIP_CONS**    linkingcuts;
   int            nlinkingcuts;
   SCIP_VAR**     linkingvars;
   int            nlinkingvars;
   SCIP_HASHMAP*  vartoblock;
   SCIP_HASHMAP*  constoblock;
   SCIP_HASHMAP*  varindex;
   SCIP_HASHMAP*  consindex;
   DECDECOMPTYPE  type;
};

typedef struct DecDecomp DECDECOMP;
#endif /* STRUCT_DECOMP_H_ */
