/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2012 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   event_origdiving.h
 * @ingroup EVENTS 
 * @brief  eventhdlr for origdiving event
 * @author Christian Puchert
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_EVENT_ORIGDIVING_H__
#define __SCIP_EVENT_ORIGDIVING_H__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** informs the event handler that a diving heuristic has been called */
extern
SCIP_RETCODE GCGeventOrigdivingCalled(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< diving heuristic */
   );

/** informs the event handler that a diving heuristic has found a new solution */
extern
SCIP_RETCODE GCGeventOrigdivingNewDivingsol(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_SOL*             sol                 /**< new solution */
   );

/** updates diving loop statistics of a diving heuristic */
extern
SCIP_RETCODE GCGeventOrigdivingDiveround(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur                /**< diving heuristic */   
   );

/** updates LP statistics of a diving heuristic */
extern
SCIP_RETCODE GCGeventOrigdivingUpdateLPstats(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_HEUR*            heur,               /**< diving heuristic */   
   SCIP_Longint          nlpiters,           /**< number of new LP iterations */
   int                   npricerounds        /**< number of new pricing rounds */
   );

/** creates event handler for origdiving event */
extern
SCIP_RETCODE SCIPincludeEventHdlrOrigdiving(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
