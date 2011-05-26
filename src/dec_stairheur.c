/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2010 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not email to scip@zib.de.      */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
#pragma ident "@(#) $Id: presol_stairheur.c,v 1.24 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   presol_stairheur.c
 * @ingroup PRESOLVERS
 * @brief  stairheur presolver
 * @author Martin Bergner
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG
#include <assert.h>

#include "dec_stairheur.h"
#include "scip_misc.h"

#define PRESOL_NAME            "stairheur"
#define PRESOL_DESC            "presolver template"
#define PRESOL_PRIORITY               -9000001 /**< priority of the presolver (>= 0: before, < 0: after constraint handlers) */
#define PRESOL_MAXROUNDS              -1 /**< maximal number of presolving rounds the presolver participates in (-1: no limit) */
#define PRESOL_DELAY                TRUE /**< should presolver be delayed, if other presolvers found reductions? */
#define DEFAULT_PASSES                10 /**< Number of passes for the whole algorithm */
#define DEFAULT_ROCPASSES            100 /**< Number of passes for the roc algorithm */
#define LISTSIZE_FACTOR                2 /**< the factor that the ordered* arrays are bigger then the scip arrays (INTEGRAL)*/
/*
 * Data structures
 */
typedef struct lentry ListEntry;
struct lentry {
   ListEntry *next;
   ListEntry *prev;
   void*      ptr;
   int index;
   int oldindex;
};

/** presolver data */
struct SCIP_StairheurData
{
   SCIP_HASHMAP* consentry;
   SCIP_HASHMAP* varentry;
   SCIP_HASHMAP* varcons;
   ListEntry*    conssfront;
   ListEntry*    varsfront;
   ListEntry*    conssback;
   ListEntry*    varsback;
   int           rocpasses;
   int           passes;
   int           nvars;
   int           nconss;
};

struct consArray
{
   SCIP_CONS** conss;
   int         nconss;
};
typedef struct consArray STAIRHEUR_ConsArray;

/*
 * Local methods
 */
/* put your local methods here, and declare them static */

/** deletes an entry from the list */
static
SCIP_RETCODE deleteListEntry(
   ListEntry*  entry,      /**< the element to be deleted */
   ListEntry** front,      /**< pointer to the beginning of the list */
   ListEntry** back        /**< pointer to the end of the list */
   )
{
   assert(entry != NULL);
   assert((*front)->prev == NULL);
   assert((*back)->next == NULL);
   if(entry->next != NULL)
   {
      entry->next->prev = entry->prev;
   }
   if(entry->prev != NULL)
   {
      entry->prev->next = entry->next;
   }

   /* special cases of back and front*/
   if(entry == *back)
   {
      *back = entry->prev;
   }
   else if(entry == *front)
   {
      *front = entry->next;
   }

   entry->next = NULL;
   entry->prev = NULL;

   assert((*front)->prev == NULL);
   assert((*back)->next == NULL);

   return SCIP_OKAY;
}

/** puts an entry to the front of the list */
static
SCIP_RETCODE pushEntryFront(
   ListEntry** list,    /**< pointer to the beginning of list */
   ListEntry*  entry    /**< the element to be prepended */
   )
{
   assert(list != NULL);
   assert((*list)->prev == NULL);
   assert((*list)->next != NULL); /* this may in bad instances, but its not critical */

   assert(entry->next == NULL);
   assert(entry->prev == NULL);
   entry->next = *list;
   (*list)->prev = entry;
   entry->index = (*list)->index-1;
   *list = entry;

   assert((*list)->prev == NULL);
   return SCIP_OKAY;
}

/** enumerates the entries of the list, signals if the indices have changed since the last run */
static
SCIP_Bool calculateIndices(
   ListEntry* list   /**< pointer to the list */
   )
{
   int i;
   ListEntry *current;
   SCIP_Bool changed;
   assert(list != NULL);

   changed = FALSE;
   i = 0;
   for(current = list; current != NULL; current = current->next)
   {
      assert(current->next != list);
      current->index = i;
      if(current->oldindex != current->index)
         changed = TRUE;
      current->oldindex = current->index;
      ++i;
   }
   return changed;
}

/** initialize the list for variables */
static
SCIP_RETCODE initialiseVarList(
   SCIP*          scip,    /**< SCIP data structure */
   ListEntry**    list,    /**< Pointer to the list start */
   ListEntry**    back,    /**< Pointer to the back of the list */
   SCIP_HASHMAP*  varmap   /**< Map mapping variables to the list element */
   )
{
   ListEntry *curr;
   ListEntry *prev;
   SCIP_VAR** vars;
   int nvars;
   int i;

   assert(scip != NULL);
   assert(list != NULL);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   prev = NULL;
   assert(nvars > 0);

   /* do the first one separately*/
   assert(SCIPvarIsActive(vars[0]));
   SCIP_CALL(SCIPallocMemory(scip, list));
   (*list)->index = 0;
   (*list)->oldindex = 0;
   (*list)->ptr = vars[0];
   (*list)->next = NULL;
   (*list)->prev = NULL;
   prev = *list;
   *back = *list;
   SCIP_CALL(SCIPhashmapSetImage(varmap, vars[0], *list ));
   for(i = 1; i < nvars; ++i)
   {
      assert(SCIPvarIsActive(vars[i]));
      SCIP_CALL(SCIPallocBlockMemory(scip, &curr));
      curr->index = i;
      curr->oldindex = i;
      curr->ptr = vars[i];
      curr->next = NULL;
      curr->prev = prev;
      prev->next = curr;
      prev = curr;
      *back = curr;
      SCIP_CALL(SCIPhashmapSetImage(varmap, vars[i], curr ));
   }

   assert((*list)->prev == NULL);
   assert((*back)->next == NULL);
   return SCIP_OKAY;
}

/** initialize the list for constraints */
static
SCIP_RETCODE initialiseConsList(
   SCIP*          scip,    /**< SCIP data structure */
   ListEntry**    list,    /**< Pointer to the beginning of the list */
   ListEntry**    back,    /**< Pointer to the end of the list */
   SCIP_HASHMAP*  consmap  /**< Hashmap mapping constrainst to list entries */
   )
{
   ListEntry *curr;
   ListEntry *prev;
   SCIP_CONS** conss;
   int nconss;
   int i;
   int j;

   assert(scip != NULL);
   assert(list != NULL);

   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);
   prev = NULL;
   assert(nconss > 0);

   /* do the first one separately*/
   assert(SCIPconsIsActive(conss[0]));
   SCIP_CALL(SCIPallocMemory(scip, list));
   (*list)->index = 0;
   (*list)->oldindex = 0;
   (*list)->ptr = conss[0];
   (*list)->next = NULL;
   (*list)->prev = NULL;
   prev = *list;
   *back = *list;
   SCIP_CALL(SCIPhashmapSetImage(consmap, conss[0], *list));
   for(i = 1; i < nconss; ++i)
   {
      SCIP_Bool hasvars;
      int nvars;
      SCIP_VAR** vars;
      assert(SCIPconsIsActive(conss[i]));
      hasvars = FALSE;
      vars = SCIPgetVarsXXX(scip, conss[i]);
      nvars = SCIPgetNVarsXXX(scip, conss[i]);

      for(j = 0; j < nvars && !hasvars; ++j)
      {
         hasvars = SCIPvarIsActive(vars[j]);
      }
      SCIPfreeMemoryArray(scip, &vars);
      if(!hasvars)
         continue;
      SCIP_CALL(SCIPallocBlockMemory(scip, &curr));
      curr->index = i;
      curr->oldindex = i;
      curr->ptr = conss[i];
      curr->next = NULL;
      curr->prev = prev;
      prev->next = curr;
      prev = curr;
      *back = curr;
      SCIP_CALL(SCIPhashmapSetImage(consmap, conss[i], curr));
   }
   assert((*list)->prev == NULL);
   assert((*back)->next == NULL);
   return SCIP_OKAY;
}

/** walks through the list and frees every element */
static
void freeList(
   SCIP*       scip,    /**< SCIP data structure */
   ListEntry** list     /**< Pointer to the list */
   )
{
   ListEntry *curr;
   ListEntry *next;
   assert(scip != NULL);
   assert(list != NULL);
   for(curr = *list; curr != NULL; curr = next)
   {
      next = curr->next;
      SCIPfreeBlockMemoryNull(scip, &curr);
   }
}

/** returns the constraints a variable belongs to */
static
SCIP_CONS** SCIPgetVarConss(
   SCIP_STAIRHEURDATA*  stairheurdata,    /**< presolver data data structure */
   SCIP_VAR*         var            /**< the variable for which the constrains should be returned */
   )
{
   STAIRHEUR_ConsArray* consarray;
   assert(stairheurdata != NULL);
   assert(var != NULL);
   assert(SCIPvarIsActive(var));
   consarray = (STAIRHEUR_ConsArray*) SCIPhashmapGetImage(stairheurdata->varcons, var);
   if(consarray == NULL)
      return NULL;
   else
      return consarray->conss;
}

/** returns the number of constraints a variable belongs to */
static
int SCIPgetNVarConss(
   SCIP_STAIRHEURDATA*  stairheurdata,    /**< presolver data data structure */
   SCIP_VAR*         var            /**< the variable for which the number should be returned */
   )
{
   STAIRHEUR_ConsArray* consarray;
   assert(stairheurdata != NULL);
   assert(var != NULL);
   assert(SCIPvarIsActive(var));
   consarray = (STAIRHEUR_ConsArray*) SCIPhashmapGetImage(stairheurdata->varcons, var);
   if(consarray == NULL)
      return 0;
   else
      return consarray->nconss;
}

/** runs the roc algorithm from J.R. King,  and V. Nakornchai */
static
SCIP_RETCODE runROC(
   SCIP*             scip,       /**< SCIP data structure */
   SCIP_STAIRHEURDATA*  stairheurdata  /**< presolver data data structure */
   )
{
   SCIP_CONS** conss;
   SCIP_VAR** vars;

   int nconss;
   int nvars;
   ListEntry *curr;
   int j;
   SCIP_Bool changed;
   int*        indices;
   ListEntry **listentries;

   /*
    * will hold preeliminary indices, the higher the index the earlier the
    * variable in the list, serves as a fast indicator in order to find out the
    * relative position of two variables and constraints in the array
    */
   SCIPdebugMessage("Starting one round of the ROC algorithm.\n");
   assert(scip != NULL);
   assert(stairheurdata != NULL);
   conss = SCIPgetConss(scip);
   nconss = SCIPgetNConss(scip);

   vars = SCIPgetVars(scip);
   nvars = SCIPgetNVars(scip);
   changed = FALSE;

   SCIP_CALL(SCIPallocMemoryArray(scip, &indices, MAX(nconss, nvars)));
   SCIP_CALL(SCIPallocMemoryArray(scip, &listentries, MAX(nconss, nvars)));

   assert(indices != NULL);
   assert(listentries != NULL);
   while(changed == FALSE)
   {

#ifndef NDEBUG /* visualisation of the matrix */
      char *array;
#endif

      /* go backwards through the variables and sort the constraints */
      for( curr = stairheurdata->varsback; curr != NULL ; curr = curr->prev)
      {
         SCIP_VAR*   var;
         SCIP_CONS** curconss;
         int         ncurcons;

         var = (SCIP_VAR*)curr->ptr;
         assert(SCIPvarIsActive(var));
         curconss = SCIPgetVarConss(stairheurdata, var);
         ncurcons = SCIPgetNVarConss(stairheurdata, var);

         for( j = 0; j < ncurcons; ++j) /* Loop over the constraints of each variable */
         {
            ListEntry *entry;
            /* add them to a local array */
            entry = (ListEntry*)SCIPhashmapGetImage(stairheurdata->consentry, curconss[j]);
            assert(entry != NULL);

            assert(listentries != NULL);
            listentries[j] = entry;

            /* add the index array in the same order as the pointers */
            assert(indices != NULL);
            indices[j] = listentries[j]->index;

            /* delete the old ones from the list */
            SCIP_CALL(deleteListEntry(listentries[j], &stairheurdata->conssfront, &stairheurdata->conssback));
            assert(listentries[j]->prev == NULL);
            assert(listentries[j]->next == NULL);
            assert(stairheurdata->conssfront->prev == NULL);
            assert(stairheurdata->conssback->next == NULL);
         }
         /* sort both the pointers and the index array according to the index array */
         SCIPsortIntPtr(indices, (void**)listentries, ncurcons);

         /* put the old ones to the front of the list */
         for( j = ncurcons-1; j >= 0; --j)
         {
            SCIP_CALL(pushEntryFront(&stairheurdata->conssfront, listentries[j]));
            assert(stairheurdata->conssfront->prev == NULL);
            assert(stairheurdata->conssfront == listentries[j]);
         }

      }
      /* go backwards through the constraints and sort the variables */
      for( curr = stairheurdata->conssback; curr != NULL; curr = curr->prev)
      {
         SCIP_CONS* cons;
         SCIP_VAR** curvars;
         int        ncurvars;
         int        nactivevars;

         cons = (SCIP_CONS*) curr->ptr;
         assert(SCIPconsIsActive(cons));
         curvars = SCIPgetVarsXXX(scip, cons);
         ncurvars = SCIPgetNVarsXXX(scip, cons);

         nactivevars = 0;
         for( j = 0; j < ncurvars; ++j) /* Loop over variable of each constraint */
         {
            if(SCIPvarIsActive(curvars[j]))
            {
               ListEntry* entry;
               /* add them to a local array */
               entry = (ListEntry*)SCIPhashmapGetImage(stairheurdata->varentry, curvars[j]);
               assert(entry != NULL);

               assert(listentries != NULL);
               listentries[nactivevars] = entry;
               assert(SCIPvarIsActive(curvars[j]));

               /* add the index array in the same order as the pointers */
               assert(indices != NULL);
               indices[nactivevars] = listentries[nactivevars]->index;

               /* delete them from the ordered constraints */
               SCIP_CALL(deleteListEntry(listentries[nactivevars], &stairheurdata->varsfront, &stairheurdata->varsback));
               assert(listentries[nactivevars]->prev == NULL);
               assert(listentries[nactivevars]->next == NULL);
               assert(stairheurdata->varsfront->prev == NULL);
               assert(stairheurdata->varsback->next == NULL);
               nactivevars++;
            }
         }
         SCIPfreeMemoryArray(scip, &curvars);
         /* sort both the pointers and the index array according to the index array */
         SCIPsortIntPtr(indices, (void**)listentries, nactivevars);

         /* put the old ones to the front of the list */
         for( j = nactivevars-1; j >=0 ; --j)
         {
            SCIP_CALL(pushEntryFront(&stairheurdata->varsfront, listentries[j]));
            assert(stairheurdata->varsfront->prev == NULL);
            assert(stairheurdata->varsfront == listentries[j]);
         }
      }

      changed = !(calculateIndices(stairheurdata->varsfront) || calculateIndices(stairheurdata->conssfront));
#if 0 /* visualisation of the matrix */
      SCIP_CALL(SCIPallocMemoryArray(scip, &array, nvars));
      for(curr = stairheurdata->conssfront; curr != NULL; curr = curr->next)
      {
         SCIP_VAR** curvars;
         int ncurvars;
         SCIP_CONS* cons = (SCIP_CONS*)curr->ptr;
         ListEntry *entry;
         for(j = 0; j < nvars; ++j)
         {
            array[j] = ' ';
         }
         curvars = SCIPgetVarsXXX(scip, cons);
         ncurvars = SCIPgetNVarsXXX(scip, cons);
         entry = (ListEntry*) SCIPhashmapGetImage(stairheurdata->consentry, cons);
         assert(listentries != NULL);

         SCIPinfoMessage(scip, NULL, "%d\t", entry->oldindex );
         for(j = 0; j < ncurvars; ++j)
         {
            if(!SCIPvarIsActive(curvars[j]))
               continue;
            entry = (ListEntry*) SCIPhashmapGetImage(stairheurdata->varentry, curvars[j]);
            assert(listentries != NULL);
            array[entry->oldindex] = '-';
            array[entry->index] = 'x';
         }
         for(j = 0; j < nvars; ++j)
         {
            SCIPinfoMessage(scip, NULL, "%c", array[j]);
         }
         SCIPinfoMessage(scip, NULL, "\n");
      }
      SCIPinfoMessage(scip, NULL, "\n\n");
      SCIPfreeMemoryArray(scip, &array);
#endif
   }

   SCIPfreeMemoryArray(scip, &indices);
   SCIPfreeMemoryArray(scip, &listentries);
   return SCIP_OKAY;
}

/** tries to detect a staircase structure signals the number of blocks */
static
SCIP_RETCODE checkStaircase(
   SCIP*                scip,           /**< SCIP data structure */
   SCIP_STAIRHEURDATA*  stairheurdata,  /**< presolver data data structure */
   int*                 jend,           /**< the last nonzero row in a variable */
   int*                 jbegin,         /**< the first nonzero row in a variable */
   int*                 iend,           /**< the last nonzero variable in a row */
   int*                 ibegin,         /**< the first nonzero variable in a row */
   int*                 result          /**< signals the number of detected blocks */
   )
{
   int nconss;
   int nvars;
   int *jmin;
   int *jmax;
   int *vmin;
   int i;
   int max;
   int min;
   int n;
   int v;
   FILE *fp;
   ListEntry *curr;

   SCIPdebugMessage("Checking staircase structure.\n");
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   SCIP_CALL(SCIPallocMemoryArray(scip, &jmin, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &jmax, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &vmin, nconss));

   /* initialize data arrays*/
   for( i = 0; i < nconss; ++i)
   {
      jmin[i] = -1;
      jmax[i] = 0;
      vmin[i] = -1;
   }

   max = 0;
   n = 0;
   v = SCIPgetNVars(scip) + 100;
   min = SCIPgetNVars(scip) + 100;
   for( curr = stairheurdata->conssfront; curr != NULL; curr = curr->next)
   {
      n = MAX(n, iend[curr->index]-ibegin[curr->index]);
      v = MIN(v, iend[curr->index]-ibegin[curr->index]);
      assert(v >= 0);

      /* test the conformity of the indices array, namely
       *  IBEG ( i ) < IBEG(i + k) for all positive k (paper page 234)
       */
      max = MAX(max, iend[curr->index]);
      jmin[curr->index] = ibegin[curr->index];
      jmax[curr->index] = max;
      if( curr->index != 0 )
      {
         vmin[curr->index-1] =1+ jmax[curr->index-1] - jmin[curr->index];
         min = MIN(min, vmin[curr->index-1]);
         assert(min >= 0);
      }
   }
   vmin[nconss-1] = 0;

   if( min > nvars / 2 )
      *result = -1;
   else
   {
      assert(n - v > 0);
      *result = (nvars - v) / (n - v);
   }
   SCIPdebugMessage("Writing vmin information to file.\n");
   fp = fopen("pvmin.gp", "w");
   if(fp != NULL)
   {
      SCIPinfoMessage(scip, fp, "plot \"-\" with p ps 0.5\n");
      for(i = 0; i < nconss; ++i)
      {
         SCIPinfoMessage(scip, fp, "%d\n", vmin[i]);
      }
      fclose(fp);
   }

   SCIPfreeMemoryArray(scip, &jmin);
   SCIPfreeMemoryArray(scip, &jmax);
   SCIPfreeMemoryArray(scip, &vmin);

   return SCIP_OKAY;
}

/** runs all parts of the clustering algorithm */
static
SCIP_RETCODE runClusteringHeur(
   SCIP*                scip,           /**< SCIP data structure */
   SCIP_STAIRHEURDATA*  stairheurdata   /**< presolver data data structure */
   )
{
   /* initialize structures*/
   int nconss;
   int nvars;

   /* endless loop*/
   int pass;
   int i;
   int j;

   int* ibegin;
   int* iend;
   int result;
   int* jbegin;
   int* jend;
   ListEntry *curr;
   assert(scip != 0);
   assert(stairheurdata != 0);
   /* compute indices */
   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   SCIP_CALL(SCIPallocMemoryArray(scip, &ibegin, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &iend, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &jbegin, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &jend, nvars));
   assert(ibegin != NULL);
   assert(iend != NULL);
   for(i = 0; i < nconss; ++i)
   {
      iend[i] = 0;
      ibegin[i] = nvars;
   }
   for( pass = 0; pass < stairheurdata->passes; ++pass )
   {
      SCIPdebugMessage("Clustering pass %d/%d.\n", pass+1 ,stairheurdata->passes );
      /* run rank ordering clustering */
      SCIP_CALL(runROC(scip, stairheurdata));

      SCIPdebugMessage("Building ibegin and iend.\n");
      for( curr = stairheurdata->conssfront; curr != NULL; curr = curr->next )
      {
         SCIP_CONS* cons;
         SCIP_VAR** curvars;

         int ncurvars;
         int max = 0;
         int min;

         min = SCIPgetNVars(scip) + 1;
         cons = (SCIP_CONS*)curr->ptr;
         assert(cons != NULL);

         ncurvars = SCIPgetNVarsXXX(scip, cons);
         curvars = SCIPgetVarsXXX(scip, cons);
         /* find out max and min index */

         for( j = 0; j < ncurvars; ++j )
         {
            ListEntry *entry;
            if(!SCIPvarIsActive(curvars[j]))
               continue;

            entry = (ListEntry*) SCIPhashmapGetImage(stairheurdata->varentry, curvars[j]);
            assert(entry != NULL);
            assert(entry->index < nvars);

            max = MAX(max, entry->index);
            min = MIN(min, entry->index);
         }
         SCIPfreeMemoryArray(scip, &curvars);
         assert(curr->index >= 0);
         assert(curr->index < nconss);
         assert(iend != NULL);
         iend[curr->index] = max;
         assert(ibegin != NULL);
         ibegin[curr->index] = min;
         if(curr->index > 0)
         {
            assert(ibegin[curr->index] >= ibegin[curr->index-1]);
         }
      }
#ifndef NDEBUG
      for( i = 0; i < nconss; ++i)
      {
         for(j = 1; j+i < nconss; ++j)
         {
            /* check conformity of the indices array */
            assert(ibegin[i] <= ibegin[i+j] );
         }
      }
#endif
      SCIPdebugMessage("Builing jbegin and jend\n");
      for( curr = stairheurdata->varsfront; curr != NULL; curr = curr->next )
      {
         SCIP_VAR* var;
         int max = 0;
         int min = SCIPgetNConss(scip) + 1;

         SCIP_CONS** curconss;
         int         ncurconss;

         var = (SCIP_VAR*) curr->ptr;
         assert(var != NULL);
         assert(SCIPvarIsActive(var));
         assert(curr->index >= 0);
         assert(curr->index < nvars);

         curconss = SCIPgetVarConss(stairheurdata, var);
         ncurconss = SCIPgetNVarConss(stairheurdata, var);

         for( j = 0; j < ncurconss; ++j)
         {
            ListEntry* entry;
            entry = (ListEntry*)SCIPhashmapGetImage(stairheurdata->consentry, curconss[j]);
            assert(entry != NULL);
            assert(entry->index < nconss);

            /* compute minimal and maximal row index for the variable */
            max = MAX(max, entry->index);
            min = MIN(min, entry->index);
         }
         assert(jend != NULL);
         jend[curr->index] = max;
         assert(jbegin != NULL);
         jbegin[curr->index] = min;
         if(curr->index > 0)
         {
            assert(jbegin[curr->index] >= jbegin[curr->index-1]);
         }
         /* check conformity of the indices array */
#ifndef NDEBUG
         for( i = 0; i < nvars; ++i)
         {
            for(j = 1; j+i < nvars; ++j)
            {
               /* check conformity of the indices array */
               assert(jbegin[i] <= jbegin[i+j] );
            }
         }
#endif
      }

      SCIP_CALL(checkStaircase(scip, stairheurdata, jend, jbegin, iend, ibegin, &result));

   }
   /* check for staircase structure */
   SCIP_CALL(checkStaircase(scip, stairheurdata, jend, jbegin, iend, ibegin, &result));

   SCIPfreeMemoryArray(scip, &ibegin);
   SCIPfreeMemoryArray(scip, &iend);
   SCIPfreeMemoryArray(scip, &jbegin);
   SCIPfreeMemoryArray(scip, &jend);

   /* TODO: Blocking */

   if( result == -1 )
      return SCIP_ERROR;

   SCIPinfoMessage(scip, NULL, "Blocks: %d\n", result);
   return SCIP_OKAY;
}

/** presolving initialization method of presolver (called when presolving is about to begin) */
#if 0
static
SCIP_DECL_PRESOLINITPRE(presolInitpreStairheur)
{
   SCIP_STAIRHEURDATA* stairheurdata;
   int i;
   int j;
   int nconss;
   int nvars;

   assert(scip != NULL);
   assert(presol != NULL);

   stairheurdata = SCIPpresolGetData(presol);
   assert(stairheurdata != NULL);

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);

   SCIP_CALL(SCIPallocMemoryArray(scip, &stairheurdata->orderedconss, LISTSIZE_FACTOR*nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &stairheurdata->orderedvars, LISTSIZE_FACTOR*nvars));
   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->consindex, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->varindex, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->varcons, SCIPblkmem(scip), nvars));
   for( i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var;
      STAIRHEUR_ConsArray *consarray;
      var = SCIPgetVars(scip)[i];
      SCIP_CALL(SCIPallocMemory(scip, &consarray));
      SCIP_CALL(SCIPallocMemoryArray(scip, &consarray->conss, nconss)); /* TODO: this is a lot of memory needed here, can be made more clever! */
      consarray->nconss = 0;
      SCIP_CALL(SCIPhashmapSetImage(stairheurdata->varcons, var, consarray));
   }
   stairheurdata->frontindconss = LISTSIZE_FACTOR*nconss-1;
   stairheurdata->backindconss = LISTSIZE_FACTOR*nconss-1;
   stairheurdata->frontindvars = LISTSIZE_FACTOR*nvars-1;
   stairheurdata->backindvars = LISTSIZE_FACTOR*nvars-1;


   /* sort the constraints and variables to the end of the arrays */
   for(i = 0; i < nconss; ++i)
   {
      SCIP_CONS* cons;
      SCIP_VAR** curvars;
      int        ncurvars;

      cons = SCIPgetConss(scip)[i];
      /* set the image of the current constraint and add it to the end of the
       * cons index array */
      assert(SCIPhashmapGetImage(stairheurdata->consindex, cons) == NULL);
      stairheurdata->orderedconss[stairheurdata->frontindconss] = cons;
      SCIP_CALL(SCIPhashmapSetImage(stairheurdata->consindex, cons, (void*)stairheurdata->frontindconss));
      --(stairheurdata->frontindconss);
      assert(stairheurdata->frontindconss >= 0);

      curvars = SCIPgetVarsXXX(scip, cons);
      ncurvars = SCIPgetNVarsXXX(scip, cons);

      for( j = 0; j < ncurvars; ++j)
      {
         STAIRHEUR_ConsArray *consarray;
         /* if the variable is not active or it has already been handled, skip it*/
         if(!SCIPvarIsActive(curvars[j]) || SCIPhashmapGetImage(stairheurdata->varindex, curvars[j]) != NULL)
         {
            continue;
         }
         /* add the constraint to the variables cons array */
         consarray = (STAIRHEUR_ConsArray*)SCIPhashmapGetImage(stairheurdata->varcons, curvars[j]);
         assert(consarray != NULL);
         consarray->conss[consarray->nconss] = cons;
         ++(consarray->nconss);

         assert(SCIPhashmapGetImage(stairheurdata->varindex, curvars[j]) == NULL);
         stairheurdata->orderedvars[stairheurdata->frontindvars] = curvars[j];
         SCIP_CALL(SCIPhashmapSetImage(stairheurdata->varindex, curvars[j], (void*)stairheurdata->frontindvars));
         --(stairheurdata->frontindvars);
         assert(stairheurdata->frontindvars >= 0);
      }
   }
   for(i = stairheurdata->frontindvars; i >= 0; --i)
   {
      stairheurdata->orderedvars[i] = NULL;
   }
   for(i = stairheurdata->frontindconss; i >= 0; --i)
   {
      stairheurdata->orderedconss[i] = NULL;
   }
   return SCIP_OKAY;
}
#endif



/** execution method of presolver */
extern
SCIP_RETCODE detectStructureStairheur(
      SCIP*                 scip,           /**< SCIP data structure       */
      SCIP_STAIRHEURDATA*   stairheurdata,  /**< stairheur data structure  */
      SCIP_RESULT*          result          /**< pointer to hold result    */
      )
{

   int i;
   int j;
   int nconss;
   int nvars;

   assert(scip != NULL);
   assert(stairheurdata != NULL);

   nconss = SCIPgetNConss(scip);
   nvars = SCIPgetNVars(scip);
   stairheurdata->nconss = nconss;
   stairheurdata->nvars = nvars;

   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->consentry, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->varentry, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&stairheurdata->varcons, SCIPblkmem(scip), nvars));

   for( i = 0; i < nvars; ++i)
   {
      SCIP_VAR* var;
      STAIRHEUR_ConsArray *consarray;
      var = SCIPgetVars(scip)[i];
      SCIP_CALL(SCIPallocMemory(scip, &consarray));
      SCIP_CALL(SCIPallocMemoryArray(scip, &consarray->conss, nconss)); /* TODO: this is a lot of memory needed here, can be made more clever! */
      consarray->nconss = 0;
      SCIP_CALL(SCIPhashmapSetImage(stairheurdata->varcons, var, consarray));
   }

   /* sort the constraints and variables to the end of the arrays */
   for(i = 0; i < nconss; ++i)
   {
      SCIP_CONS* cons;
      SCIP_VAR** curvars;
      int        ncurvars;

      cons = SCIPgetConss(scip)[i];
      /* set the image of the current constraint and add it to the end of the
       * cons index array */

      curvars = SCIPgetVarsXXX(scip, cons);
      ncurvars = SCIPgetNVarsXXX(scip, cons);

      for( j = 0; j < ncurvars; ++j)
      {
         STAIRHEUR_ConsArray *consarray;
         /* if the variable is not active, skip it*/
         if(!SCIPvarIsActive(curvars[j]))
         {
            continue;
         }
         /* add the constraint to the variables cons array */
         consarray = (STAIRHEUR_ConsArray*)SCIPhashmapGetImage(stairheurdata->varcons, curvars[j]);
         assert(consarray != NULL);
         consarray->conss[consarray->nconss] = cons;
         ++(consarray->nconss);
      }
      SCIPfreeMemoryArray(scip, &curvars);
   }
   SCIP_CALL(initialiseVarList(scip, &stairheurdata->varsfront, &stairheurdata->varsback, stairheurdata->varentry));
   SCIP_CALL(initialiseConsList(scip, &stairheurdata->conssfront, &stairheurdata->conssback, stairheurdata->consentry));

   SCIP_CALL(runClusteringHeur(scip, stairheurdata));

   freeList(scip, &stairheurdata->varsfront);
   freeList(scip, &stairheurdata->conssfront);
   SCIPhashmapFree(&stairheurdata->consentry);
   SCIPhashmapFree(&stairheurdata->varentry);
   for( i = 0; i < SCIPgetNVars(scip); ++i)
   {
      SCIP_VAR* var;
      STAIRHEUR_ConsArray *consarray;
      var = SCIPgetVars(scip)[i];
      consarray = (STAIRHEUR_ConsArray*)SCIPhashmapGetImage(stairheurdata->varcons, var);
      assert(consarray != NULL);
      SCIPfreeMemoryArray(scip, &consarray->conss);
      SCIPfreeMemory(scip, &consarray);
      SCIP_CALL(SCIPhashmapRemove(stairheurdata->varcons, var));
   }
   SCIPhashmapFree(&stairheurdata->varcons);
   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** allocates and initializes the stairheurdata */
extern
SCIP_RETCODE createStairheurData(
      SCIP*                scip,           /**< SCIP data structure */
      SCIP_STAIRHEURDATA** stairheurdata   /**< stairheur data structure */
   )
{
   assert(scip != NULL);
   assert(stairheurdata != NULL);

   SCIP_CALL(SCIPallocMemory(scip, stairheurdata));
   return SCIP_OKAY;
}

/** frees the stairheurdata */
extern
void freeStairheurData(
      SCIP*                scip,           /**< SCIP data structure */
      SCIP_STAIRHEURDATA** stairheurdata   /**< stairheur data structure */
   )
{
   assert(scip != NULL);
   assert(stairheurdata != NULL);
   SCIPfreeMemory(scip, stairheurdata)
}

/** creates the stairheur detection and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludeDetectionStairheur(
      SCIP*                 scip, /**< SCIP data structure */
      SCIP_STAIRHEURDATA*   stairheurdata
)
{


   assert(stairheurdata != NULL);

   /* add stairheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "stairheur/passes", "Maximal number of passes of the wholealgorithm", &stairheurdata->passes, FALSE, DEFAULT_PASSES, 1, INT_MAX, NULL, NULL ) );
   SCIP_CALL( SCIPaddIntParam(scip, "stairheur/rocpasses", "Maximal number of passes of the roc algorithm", &stairheurdata->rocpasses, FALSE, DEFAULT_ROCPASSES, 1, INT_MAX, NULL, NULL ) );
   return SCIP_OKAY;
}
