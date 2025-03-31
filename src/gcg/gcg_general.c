/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* Copyright (C) 2010-2025 Operations Research, RWTH Aachen University       */
/*                         Zuse Institute Berlin (ZIB)                       */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with GCG; see the file LICENSE. If not visit gcg.or.rwth-aachen.de.*/
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   gcg_general.c
 * @brief  gcg general public methods
 * @author Steffan Schlein
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "gcg/gcg_general.h"
#include "pub_gcg.h"
#include "gcg/gcggithash.h"

#include "scip/pub_message.h"
#include "scip/scip_message.h"


/** returns complete GCG version number in the format "major . minor tech"
 *
 *  @return complete GCG version
 */
SCIP_Real GCGversion(
   void
   )
{
   return (SCIP_Real)(GCG_VERSION)/100.0;
}

/** returns GCG major version
 *
 *  @return major GCG version
 */
int GCGmajorVersion(
   void
   )
{
   return GCG_VERSION/100;
}

/** returns GCG minor version
 *
 *  @return minor GCG version
 */
int GCGminorVersion(
   void
   )
{
   return (GCG_VERSION/10) % 10;
}

/** returns GCG technical version
 *
 *  @return technical GCG version
 */
int GCGtechVersion(
   void
   )
{
   return GCG_VERSION % 10;
}

/** returns GCG sub version number
 *
 *  @return subversion GCG version
 */
int GCGsubversion(
   void
   )
{
   return GCG_SUBVERSION;
}

/** prints out GCG version
 */
void GCGprintVersion(
   GCG*                  gcg,                /**< GCG data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIP* scip;
   assert(gcg != NULL);
   scip = GCGgetOrigprob(gcg);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "GCG version %d.%d.%d",
      GCGmajorVersion(), GCGminorVersion(), GCGtechVersion());
#if GCG_SUBVERSION > 0
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ".%d", GCGsubversion());
#endif
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " [GitHash: %s]", GCGgetGitHash());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "Copyright (C) 2010-2025 Operations Research, RWTH Aachen University\n");
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "                        Konrad-Zuse-Zentrum fuer Informationstechnik Berlin (ZIB)\n\n");
}
