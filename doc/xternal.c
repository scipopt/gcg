/**@defgroup PUBLICAPI Public API of GCG
 * @brief methods and headers of the public C-API of \GCG
 *
 * \PUBLICAPIDESCRIPTION
 *
 *
 */

/**@defgroup PUBLICCOREAPI Core API
* @ingroup PUBLICAPI
* @brief methods and headers of the plugin-independent C-API provided by the \GCG header file scip.h.
*
* This module comprises methods provided by the header file scip.h. Including this header into a user-written extension
* suffices to have all plugin-independent functionality of \SCIP available. Plugin-independent
* user functionality includes the
*
* - creation of problems that \SCIP should solve
* - fine-grained access to initiate the solving process of \SCIP
* - access to all sorts of solving process statistics
* - commonly used data structures and algorithms
* - the management of plugins
* - ...
*
* In order facilitate the navigation through the core API of \SCIP, it is structured into different modules.
*/

/**@defgroup DATASTRUCTURES Data Structures
  * @ingroup PUBLICCOREAPI
  * @brief all data structures implemented in \GCG.
  */

/**@defgroup TYPEDEFINITIONS Type Definitions
  * @ingroup PUBLICCOREAPI
  */


/**@defgroup BLISS Bliss
  * @ingroup PUBLICCOREAPI
  */

/**@defgroup DECOMP Decomposition
  * @ingroup PUBLICCOREAPI
  * @brief all decomposition algorithmics implemented in \GCG.
  *
  */

/**@defgroup HEURISTICS Heuristics
  * @ingroup PUBLICCOREAPI
  */

/**@defgroup PRICING_PUB Pricing
  * @ingroup PUBLICCOREAPI
  * @brief This page contains all pricing-related public functionalities.
  *
  */

/**@defgroup PRICINGJOB Pricing Job
  * @ingroup PRICING_PUB
  */

/**@defgroup PRICINGPROB Pricing Problem
  * @ingroup PRICING_PUB
  */

/**@defgroup SEPARATORS_PUB Separators
  * @ingroup PUBLICCOREAPI
  */

/**@defgroup MISC Miscellaneous
  * @ingroup PUBLICCOREAPI
  */




/**@defgroup PUBLICPLUGINAPI Plugin API of GCG
  * @ingroup PUBLICAPI
  * @brief core API extensions provided by the default plugins of \SCIP, includable via scipdefplugins.h.
  *
  * All default plugins of \SCIP, especially the default \ref CONSHDLRS "constraint handlers", provide
  * valuable extensions to the \ref PUBLICCOREAPI "core API" of \SCIP. These methods are made available
  * by including scipdefplugins.h to user-written extensions.
  *
  * For a better overview, this page lists all default plugin headers structured into modules based on their individual
  * topic.
  *
  * All of the modules listed below provide functions that are allowed to be used by user-written extensions of \SCIP.
  */
 /**@defgroup INTERNALAPI Internal API of GCG
  * @brief internal API methods that should only be used by the core of \GCG
  *
  * This page lists the header files of internal API methods. In contrast to the public API, these internal methods
  * should not be used by user plugins and extensions of GCG. Please consult
  * \ref PUBLICCOREAPI "the Core API" and \ref PUBLICPLUGINAPI "Plugin API" for the complete API available to user plugins.
  *
  */

 /**@defgroup LPIS LP Solver Interface
  * @ingroup PUBLICPLUGINLPI
  * @brief methods and files provided by the LP solver interface of \GCG
  *
  * \SCIP uses external tools to solve LP relaxations. The communication
  * is realized through an LP interface.
  *
  * This page lists public interface methods that every LP interface provides.
  * Find the concrete implementation for your LP solver
  * under "src/lpi/".
  *
  * @see \ref LPI for a list of available LP solver interfaces
  */

  /**@defgroup BENDERS Benders' Decomposition
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a description of all methods and files provided by the Benders' decomposition.
   *
   */

  /**@defgroup BRANCHINGRULES Branching Rules
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all branching rule which are currently available.
   *
   * A detailed description what a branching rule does and how to add a branching rule to GCG can be found
   * \subpage BRANCH "here".
   */

  /**@defgroup CONSHDLRS  Constraint Handler
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all constraint handlers which are currently available.
   *
   * A detailed description what a constraint handler does and how to add a constraint handler to \SCIP can be found
   * in the SCIP documentation.
   */

  /**@defgroup DETECTORS Detectors
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all detectors which are currently available.
   *
   * A detailed description what a detector does and how to add a detector to GCG can be found
   * \ref DETECT "here".
   */

  /**@defgroup DIALOGS Dialogs
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all dialogs which are currently available.
   *
   * A detailed description what a dialog does and how to add a dialog to \SCIP can be found
   * n the SCIP documentation.
   */

  /**@defgroup DISPLAYS Displays
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all displays (output columns)  which are currently available.
   *
   * A detailed description what a display does and how to add a display to \SCIP can be found
   * in the SCIP documentation.
   *
   */

  /**@defgroup FILEREADERS File Readers
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all file readers which are currently available.
   *
   * A detailed description what a file reader does and how to add a file reader to \SCIP can be found
   * in the SCIP documentation.
   */

  /**@defgroup NODESELECTORS Node Selectors
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all node selectors which are currently available.
   *
   * A detailed description what a node selector does and how to add a node selector to \SCIP can be found
   * in the SCIP documentation.
   */

   /**@defgroup PRICING Pricing
    * @ingroup PUBLICPLUGINAPI
    * @brief This page contains a list of all pricers, pricing solvers and the pricing jobs and problem structures.
    *
    */

  /**@defgroup PRICERS Pricers
   * @ingroup PRICING
   * @brief This page contains a list of all pricers which are currently available.
   *
   * Per default there exist no variable pricer. A detailed description what a variable pricer does and how to add a
   * variable pricer to \SCIP can be found in the SCIP documentation.
   */

  /**@defgroup PRICINGSOLVERS Pricing solvers
   * @ingroup PRICING
   * @brief This page contains a list of all pricing solvers which are currently available.
   *
   * A detailed description what a pricing solver does and how to add a pricing solver to GCG can be found
   * \ref PRICINGSOLVER "here".
   */

   /**@defgroup PRIMALHEURISTICS Primal Heuristics
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all primal heuristics which are currently available.
   *
   * A detailed description what a primal heuristic does and how to add a primal heuristic to \SCIP can be found
   * \ref HEUR "here".
   */

  /**@defgroup RELAXATORS Relaxators
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all relaxators which are currently available.
   */

  /**@defgroup SEPARATORS Separators
   * @ingroup PUBLICPLUGINAPI
   * @brief This page contains a list of all separators  which are currently available.
   *
   * A detailed description what a separator does and how to add a separator to \SCIP can be found
   * in the SCIP documentation.
   */

  /**@defgroup TYPEDEFINITIONS Type Definitions
   * @ingroup PUBLICCOREAPI
   * This page lists headers containing branch-and-price specific public methods provided by GCG.
   *
   * All of the headers listed below include functions that are allowed to be called by external users. Besides those
   * functions it is also valid to call methods that are listed in one of the headers of the (default) GCG plug-ins; in
   * particular, this holds for relax_gcg.h and pricer_gcg.h.
   *
   */
