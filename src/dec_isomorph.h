/*
 * dec_isomorph.h
 *
 *  Created on: Jul 1, 2013
 *      Author: peters
 */

#ifndef DEC_ISOMORPH_H__
#define DEC_ISOMORPH_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for connected constraints and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDetectionIsomorphism(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif
#endif /* DEC_ISOMORPH_H__ */
