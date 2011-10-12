/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program                         */
/*          GCG --- Generic Column Generation                                */
/*                  a Dantzig-Wolfe decomposition based extension            */
/*                  of the branch-cut-and-price framework                    */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   main.c
 * @brief  Main file for C compilation
 * @author Gerald Gamrath
 * @author Martin Bergner
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define GCG_VERSION 90
#define GCG_SUBVERSION 2

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "gcgplugins.h"
#include "gcggithash.h"

/** returns GCG major version */
static
int GCGmajorVersion(
   void
   )
{
   return GCG_VERSION/100;
}

/** returns GCG minor version */
static
int GCGminorVersion(
   void
   )
{
   return (GCG_VERSION/10) % 10;
}

/** returns GCG technical version */
static
int GCGtechVersion(
   void
   )
{
   return GCG_VERSION % 10;
}
#if GCG_SUBVERSION > 0
/** returns GCG sub version number */
static
int GCGsubversion(
   void
   )
{
   return GCG_SUBVERSION;
}
#endif

static
void GCGprintVersion(
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   SCIPmessageFPrintInfo(file, "GCG version %d.%d.%d",
      GCGmajorVersion(), GCGminorVersion(), GCGtechVersion());
#if GCG_SUBVERSION > 0
   SCIPmessageFPrintInfo(file, ".%d", GCGsubversion());
#endif
   SCIPmessageFPrintInfo(file, " [GitHash: %s]", GCGgetGitHash());
   SCIPmessageFPrintInfo(file, "\n");
}
static
SCIP_RETCODE SCIPrunGCGShell(
   int                        argc,               /**< number of shell parameters */
   char**                     argv,               /**< array with shell parameters */
   const char*                defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   GCGprintVersion(NULL);

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );

   /* include coloring plugins */
   SCIP_CALL( SCIPincludeGcgPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessShellArguments(scip, argc, argv, defaultsetname) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}


int
main(
   int                        argc,
   char**                     argv
   )
{
  SCIP_RETCODE retcode;

  retcode = SCIPrunGCGShell(argc, argv, "gcg.set");

  if( retcode != SCIP_OKAY )
  {
     SCIPprintError(retcode, stderr);
     return -1;
  }

  return 0;
}
