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
 * @author Christian Puchert
 */

/*--+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define GCG_VERSION 100
#define GCG_SUBVERSION 0

#include <string.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/scipshell.h"
#include "gcgplugins.h"
#include "cons_decomp.h"
#include "gcggithash.h"
#include "relax_gcg.h"

#if SCIP_VERSION < 300
#error GCG can only be compiled with SCIP version 3.0.0 or higher
#endif

/** returns GCG major version */
static
int GCGmajorVersion(void)
{
   return GCG_VERSION/100;
}

/** returns GCG minor version */
static
int GCGminorVersion(void)
{
   return (GCG_VERSION/10) % 10;
}

/** returns GCG technical version */
static
int GCGtechVersion(void)
{
   return GCG_VERSION % 10;
}
#if GCG_SUBVERSION > 0
/** returns GCG sub version number */
static
int GCGsubversion(void)
{
   return GCG_SUBVERSION;
}
#endif

static
void GCGprintVersion(
   SCIP*                 scip,               /**< SCIP data structure */
   FILE*                 file                /**< output file (or NULL for standard output) */
   )
{
   assert(scip != NULL);

   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "GCG version %d.%d.%d",
      GCGmajorVersion(), GCGminorVersion(), GCGtechVersion());
#if GCG_SUBVERSION > 0
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, ".%d", GCGsubversion());
#endif
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, " [GitHash: %s]", GCGgetGitHash());
   SCIPmessageFPrintInfo(SCIPgetMessagehdlr(scip), file, "\n");
}


static
SCIP_RETCODE readParams(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename            /**< parameter file name */
   )
{
   if( SCIPfileExists(filename) )
   {
      SCIPinfoMessage(scip, NULL, "reading user parameter file <%s>\n", filename);
      SCIP_CALL( SCIPreadParams(scip, filename) );
   }
   else
      SCIPinfoMessage(scip, NULL, "user parameter file <%s> not found - using default parameters\n", filename);

   return SCIP_OKAY;
}

static
SCIP_RETCODE fromCommandLine(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           filename,           /**< input file name */
   const char*           decname             /**< decomposition file name (or NULL) */
   )
{
   SCIP_RESULT result;
   /********************
    * Problem Creation *
    ********************/

   SCIPinfoMessage(scip, NULL, "\nread problem <%s>\n", filename);
   SCIPinfoMessage(scip, NULL, "============\n\n");
   SCIP_CALL( SCIPreadProb(scip, filename, NULL) );
   if( decname != NULL )
   {
      SCIPinfoMessage(scip, NULL, "\nread decomposition <%s>\n", decname);
      SCIPinfoMessage(scip, NULL, "==================\n\n");
      SCIP_CALL( SCIPreadProb(scip, decname, NULL) );
   }


   /*******************
    * Problem Solving *
    *******************/

   /* solve problem */
   SCIPinfoMessage(scip, NULL, "\nsolve problem\n");
   SCIPinfoMessage(scip, NULL, "=============\n\n");
   SCIP_CALL( SCIPpresolve(scip) );

   if( decname == NULL )
   {
      SCIP_CALL( DECdetectStructure(scip, &result) );
      if( result != SCIP_SUCCESS )
      {
         SCIPinfoMessage(scip, NULL, "No decomposition exists and could be detected. You need to specify one.\n");
         return SCIP_OKAY;
      }
   }

   SCIP_CALL( SCIPsolve(scip) );

   SCIPinfoMessage(scip, NULL, "\nprimal solution:\n");
   SCIPinfoMessage(scip, NULL, "================\n\n");
   SCIP_CALL( SCIPprintBestSol(scip, NULL, FALSE) );


   /**************
    * Statistics *
    **************/

   SCIPinfoMessage(scip, NULL, "\nStatistics\n");
   SCIPinfoMessage(scip, NULL, "==========\n\n");

   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   return SCIP_OKAY;
}

/** evaluates command line parameters and runs GCG appropriately in the given SCIP instance */
static
SCIP_RETCODE SCIPprocessGCGShellArguments(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{  /*lint --e{850}*/
   char* probname = NULL;
   char* decname = NULL;
   char* settingsname = NULL;
   char* mastersetname = NULL;
   char* logname = NULL;
   SCIP_Bool quiet;
   SCIP_Bool paramerror;
   SCIP_Bool interactive;
   int i;

   /********************
    * Parse parameters *
    ********************/

   quiet = FALSE;
   paramerror = FALSE;
   interactive = FALSE;
   for( i = 1; i < argc; ++i )
   {
      if( strcmp(argv[i], "-l") == 0 )
      {
         i++;
         if( i < argc )
            logname = argv[i];
         else
         {
            printf("missing log filename after parameter '-l'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-q") == 0 )
         quiet = TRUE;
      else if( strcmp(argv[i], "-s") == 0 )
      {
         i++;
         if( i < argc )
            settingsname = argv[i];
         else
         {
            printf("missing settings filename after parameter '-s'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-m") == 0 )
      {
         i++;
         if( i < argc )
            mastersetname = argv[i];
         else
         {
            printf("missing master settings filename after parameter '-m'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-f") == 0 )
      {
         i++;
         if( i < argc )
            probname = argv[i];
         else
         {
            printf("missing problem filename after parameter '-f'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-d") == 0 )
      {
         i++;
         if( i < argc )
            decname = argv[i];
         else
         {
            printf("missing decomposition filename after parameter '-d'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-c") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_CALL( SCIPaddDialogInputLine(scip, argv[i]) );
            interactive = TRUE;
         }
         else
         {
            printf("missing command line after parameter '-c'\n");
            paramerror = TRUE;
         }
      }
      else if( strcmp(argv[i], "-b") == 0 )
      {
         i++;
         if( i < argc )
         {
            SCIP_FILE* file;

            file = SCIPfopen(argv[i], "r");
            if( file == NULL )
            {
               printf("cannot read command batch file <%s>\n", argv[i]);
               SCIPprintSysError(argv[i]);
               paramerror = TRUE;
            }
            else
            {
               while( !SCIPfeof(file) )
               {
                  char buffer[SCIP_MAXSTRLEN];

                  (void)SCIPfgets(buffer, sizeof(buffer), file);
                  if( buffer[0] != '\0' )
                  {
                     SCIP_CALL( SCIPaddDialogInputLine(scip, buffer) );
                  }
               }
               SCIPfclose(file);
               interactive = TRUE;
            }
         }
         else
         {
            printf("missing command batch filename after parameter '-b'\n");
            paramerror = TRUE;
         }
      }
      else
      {
         printf("invalid parameter <%s>\n", argv[i]);
         paramerror = TRUE;
      }
   }
   if( interactive && probname != NULL )
   {
      printf("cannot mix batch mode '-c' and '-b' with file mode '-f'\n");
      paramerror = TRUE;
   }
   if( probname == NULL && decname != NULL )
   {
      printf("cannot read decomposition file without given problem\n");
      paramerror = TRUE;
   }

   if( !paramerror )
   {

      /***********************************
       * create log file message handler *
       ***********************************/

      if( quiet )
      {
         SCIPsetMessagehdlrQuiet(scip, quiet);
      }

      if( logname != NULL )
      {
         SCIPsetMessagehdlrLogfile(scip, logname);
      }


      /***********************************
       * Version and library information *
       ***********************************/

      SCIPprintVersion(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      SCIPprintExternalCodes(scip, NULL);
      SCIPinfoMessage(scip, NULL, "\n");

      /*****************
       * Load settings *
       *****************/

      if( settingsname != NULL )
      {
         SCIP_CALL( readParams(scip, settingsname) );
      }
      else if( defaultsetname != NULL )
      {
         SCIP_CALL( readParams(scip, defaultsetname) );
      }

      if( mastersetname != NULL )
      {
         SCIP_CALL( readParams(GCGrelaxGetMasterprob(scip), mastersetname) );
      }

      /**************
       * Start SCIP *
       **************/

      if( probname != NULL )
      {
         SCIP_CALL( fromCommandLine(scip, probname, decname) );
      }
      else
      {
         SCIPinfoMessage(scip, NULL, "\n");
         SCIP_CALL( SCIPstartInteraction(scip) );
      }
   }
   else
   {
      printf("\nsyntax: %s [-l <logfile>] [-q] [-s <settings>] [-f <problem>] [-m <mastersettings>] [-d <decomposition>] [-b <batchfile>] [-c \"command\"]\n"
         "  -l <logfile>        : copy output into log file\n"
         "  -q                  : suppress screen messages\n"
         "  -s <settings>       : load parameter settings (.set) file\n"
         "  -m <mastersettings> : load master parameter settings (.set) file\n"
         "  -f <problem>        : load and solve problem file\n"
         "  -d <decomposition>  : load decomposition file\n"
         "  -b <batchfile>      : load and execute dialog command batch file (can be used multiple times)\n"
         "  -c \"command\"        : execute single line of dialog commands (can be used multiple times)\n\n",
         argv[0]);
   }

   return SCIP_OKAY;
}

static
SCIP_RETCODE SCIPrunGCGShell(
   int                   argc,               /**< number of shell parameters */
   char**                argv,               /**< array with shell parameters */
   const char*           defaultsetname      /**< name of default settings file */
   )
{
   SCIP* scip = NULL;

   /*********
    * Setup *
    *********/

   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(&scip) );
   GCGprintVersion(scip, NULL);


   /* include coloring plugins */
   SCIP_CALL( SCIPincludeGcgPlugins(scip) );

   /**********************************
    * Process command line arguments *
    **********************************/
   SCIP_CALL( SCIPprocessGCGShellArguments(scip, argc, argv, defaultsetname) );


   /********************
    * Deinitialization *
    ********************/

   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   return SCIP_OKAY;
}

/** main function called first */
int
main(
   int                   argc,               /**< number of arguments */
   char**                argv                /**< array of arguments */
   )
{
  SCIP_RETCODE retcode;

  retcode = SCIPrunGCGShell(argc, argv, "gcg.set");

  if( retcode != SCIP_OKAY )
  {
     SCIPprintError(retcode);
     return -1;
  }

  return 0;
}
