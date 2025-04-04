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

/**@file   xternal.c
 * @brief  documentation page for GCG's API (no other pages)
 */

/**@defgroup PUBLICAPI-GCG Public API of GCG
 * @brief methods and headers of the public API of \GCG
 *
 * The public API of \GCG is separated into a Core API and a Plugin API.
 * The first contains all methods that can be accessed by including the header gcg.h.
 * The Plugin API is a collection of methods that are provided by the default plugins of \GCG.
 * The Plugin API is provided by gcgplugins.c.
 */

/**@defgroup INTERNALAPI-GCG Internal API of GCG
 * @brief internal API methods that should only be used by the core of \GCG
 *
 * This page lists the header files of internal API methods. In contrast to the public API, these internal methods
 * should not be used by user plugins and extensions of \GCG. Please consult
 * \ref PUBLICCOREAPI "the Core API" and \ref PUBLICPLUGINAPI "Plugin API" for the complete API available to user plugins.
 */

/**@defgroup PUBLICCOREAPI Core API
* @ingroup PUBLICAPI-GCG
* @brief methods and headers of the plugin-independent API provided by \GCG.
*
* In order facilitate the navigation through the core API of \GCG, it is structured into different modules.
*/

/**@defgroup DATASTRUCTURES Data Structures
  * @ingroup PUBLICCOREAPI
  * @brief Commonly used data structures
  */

/**@defgroup TYPEDEFINITIONS Type Definitions
  * @ingroup PUBLICCOREAPI
  * @brief Type definitions and callback declarations
  */


/**@defgroup BLISS Bliss
  * @ingroup PUBLICCOREAPI
  * @brief Methods concerning BLISS
  */

/**@defgroup PublicConsClassifierMethods Constraint classifiers
  * @ingroup PUBLICCOREAPI
  * @brief Public methods for constraint classifiers.
  */

/**@defgroup PublicVarClassifierMethods Variable classifiers
  * @ingroup PUBLICCOREAPI
  * @brief Public methods for variable classifiers.
  */

/**@defgroup DECOMP Decomposition
  * @ingroup PUBLICCOREAPI
  * @brief Public methods concerning the decomposition.
  *
  */

/**@defgroup HEURISTICS Heuristics
  * @ingroup PUBLICCOREAPI
  * @brief Public methods concerning heuristics.
  */

/**@defgroup PRICING_PUB Pricing
  * @ingroup PUBLICCOREAPI
  * @brief All pricing-related public functionalities.
  *
  */

/**@defgroup PRICINGJOB Pricing Job
  * @ingroup PRICING_PUB
  */

/**@defgroup PRICINGPROB Pricing Problem
  * @ingroup PRICING_PUB
  */

/**@defgroup PRICING_PRIV Pricing
  * @ingroup INTERNALAPI-GCG
  * @brief All pricing-related internal functionalities.
  */

/**@defgroup PublicScoreMethods Scores
  * @ingroup PUBLICCOREAPI
  * @brief Public methods for scores.
  */

/**@defgroup SEPARATORS_PUB Separators
  * @ingroup PUBLICCOREAPI
  * @brief Public methods for separators.
  */

/**@defgroup MISC Miscellaneous
  * @ingroup PUBLICCOREAPI
  * @brief Public methods from the scip_misc.c file.
  */




/**@defgroup PUBLICPLUGINAPI-GCG Plugin API of GCG
  * @ingroup PUBLICAPI-GCG
  * @brief core API extensions provided by the default plugins of \GCG.
  *
  * All of the modules listed below provide functions that are allowed to be used by user-written extensions of \GCG.
  */

  /**@defgroup BENDERS-GCG Benders' Decomposition
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a description of all methods and files provided by the Benders' decomposition.
   *
   */

  /**@defgroup BRANCHINGRULES-GCG Branching Rules
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all branching rule which are currently available.
   *
   * A detailed description what a branching rule does and how to add a branching rule to \GCG can be found
   * \ref own-branching-rule "here".
   */

  /**@defgroup CONSHDLRS-GCG  Constraint Handler
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all constraint handlers which are currently available.
   *
   * A detailed description what a constraint handler does and how to add a constraint handler to \GCG can be found
   * in the SCIP documentation.
   */

  /**@defgroup DETECTORS Detectors
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all detectors which are currently available.
   *
   * A detailed description what a detector does and how to add a detector to \GCG can be found
   * \ref detection "here".
   */

  /**@defgroup CLASSIFIERS Classifiers
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all classifiers which are currently available.
   *
   * A detailed description what a classifier does can be found \ref classifiers "here"
   * and a guide on how to add a classifier to \GCG can be found \ref own-classifier "here".
   *
   */

  /**@defgroup DIALOGS Dialogs
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all dialogs which are currently available.
   *
   * A detailed description what a dialog does and how to add a dialog to \GCG can be found
   * n the SCIP documentation.
   */

  /**@defgroup DISPLAYS Displays
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all displays (output columns)  which are currently available.
   *
   * A detailed description what a display does and how to add a display to \GCG can be found
   * in the SCIP documentation.
   *
   */

  /**@defgroup FILEREADERS-GCG File Readers
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all file readers which are currently available.
   *
   * A detailed description what a file reader does and how to add a file reader to \GCG can be found
   * in the SCIP documentation.
   */

  /**@defgroup NODESELECTORS Node Selectors
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all node selectors which are currently available.
   *
   * A detailed description what a node selector does and how to add a node selector to \GCG can be found
   * in the SCIP documentation.
   */

  /**@defgroup PRICINGSOLVERS Pricing Solvers
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all pricing solvers which are currently available.
   *
   * A detailed description what a pricing solver does and how to add a pricing solver to \GCG can be found
   * \ref pricing "here".
   */

   /**@defgroup PRIMALHEURISTICS Primal Heuristics
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all primal heuristics which are currently available.
   *
   * A detailed description what a primal heuristic does and how to add a primal heuristic to \GCG can be found
   * \ref own-primal-heuristic "here".
   */

   /**@defgroup DIVINGHEURISTICS Diving Heuristics
   * @ingroup PRIMALHEURISTICS
   * @brief This page contains a list of all diving heuristics which are currently available.
   *
   * A detailed description what a diving heuristic does can be found
   * \ref diving-heuristics "here".
   */

  /**@defgroup RELAXATORS-GCG Relaxators
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all relaxators which are currently available.
   */

  /**@defgroup SCORES Scores
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all scores which are currently available.
   */

  /**@defgroup SEPARATORS Separators
   * @ingroup PUBLICPLUGINAPI-GCG
   * @brief This page contains a list of all separators  which are currently available.
   *
   * A detailed description what a separator does and how to add a separator to \GCG can be found
   * in the SCIP documentation.
   */

  /**@defgroup TYPEDEFINITIONS Type Definitions
   * @ingroup PUBLICCOREAPI
   * This page lists headers containing branch-and-price specific public methods provided by \GCG.
   *
   * All of the headers listed below include functions that are allowed to be called by external users. Besides those
   * functions it is also valid to call methods that are listed in one of the headers of the (default) \GCG plug-ins; in
   * particular, this holds for relax_gcg.h and pricer_gcg.h.
   *
   */

  /**@defgroup PARTITIONS Partitions
   * @ingroup DATASTRUCTURES
   * @brief C++ Data structures that can store partitions of variables and constraints.
   */

  /**@defgroup GCG_EXTENDEDMASTERCONSDATA GCG Extended Master Cons Data
   * @ingroup DATASTRUCTURES
   */
