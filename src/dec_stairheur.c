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

/**@file   dec_stairheur.c
 * @brief  stairheur presolver
 * @author Martin Bergner
 * @author Mathias Luers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
/* #define SCIP_DEBUG */

#include <assert.h>
#include <string.h>
#include "dec_stairheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "pub_decomp.h"
#include "scip_misc.h"
#include "scip/pub_misc.h"
#include "scip/struct_var.h"

#define DEC_DETECTORNAME      "stairheur"    /**< name of the detector */
#define DEC_DESC              "detects staircase matrices via matrix reordering" /**< detector description */
#define DEC_PRIORITY          1200           /**< priority of the detector */
#define DEC_DECCHAR           's'            /**< display character of detector */
#define DEC_ENABLED           FALSE         /**< should detector be called by default */

/* Default parameter settings*/
#define DEFAULT_MAXBLOCKS                       20       /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                       2        /**< value for the minimum number of blocks to be considered */
#define DEFAULT_PRIORITY                        DEC_PRIORITY
#define DEFAULT_DESIREDBLOCKS                   0       /**< value for the desired number of blocks (for all blocking types). Set to zero for self determination of block number */
#define DEFAULT_ENABLEBLOCKINGDYNAMIC           TRUE     /**< Enable blocking type 'dynamic' */
#define DEFAULT_ENABLEBLOCKINGSTATIC            TRUE     /**< Enable blocking type 'static' */
#define DEFAULT_ENABLEBLOCKINGASSOONASPOSSIBLE  TRUE     /**< Enable blocking type 'as soon as possible' */
#define DEFAULT_ENABLEMULTIPLEDECOMPS           TRUE     /**< Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks' */
#define DEFAULT_MAXITERATIONSROC                1000000  /**< The maximum of iterations of the ROC-algorithm. -1 for no iteration limit */

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)

/** TODO:
 * currently, all vars from the first column where a linking var appears until the end of the block are considered as linking vars, although there might be empty columns. This could be changed so that these empty columns are considered as subscipvars and not linking vars.
 *
 * In some cases a block can consist of linking vars exclusively. This makes no real sense.
 *
 * For some instances the assertion regarding the consistency of the arrays ibegin and jbegin fails
 * */

struct node {
 void* data;
 struct node* prev;
 struct node* next;
};
typedef struct node NODE;


struct list {
   /* a doubly linked list */
   NODE* nil;
   int size;
};

/** a doubly linked list */
typedef struct list LIST;

struct iterator {
   NODE* node;
   LIST* list;
};

/** an iterator for traversing lists */
typedef struct iterator ITERATOR;

/*lists*/
LIST* SCIPlistCreate(SCIP* scip);
LIST* SCIPlistCreateInt(SCIP* scip, int from, int to);
LIST* SCIPlistCopyShallow(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistPushFront(SCIP* scip, LIST* list, void* data);
SCIP_Bool SCIPlistPopFront(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistPushBack(SCIP* scip, LIST* list, void* data);
SCIP_Bool SCIPlistPopBack(SCIP* scip, LIST* list);
NODE* SCIPlistInsert(SCIP* scip, ITERATOR* it, void* data);
SCIP_Bool SCIPlistAssign(ITERATOR* it, void* data);
SCIP_Bool SCIPlistAssignList(LIST* list1, LIST* list2);
SCIP_Bool SCIPlistIsEmpty(LIST* list);
SCIP_Bool SCIPlistErase(SCIP* scip, ITERATOR* it);
SCIP_Bool SCIPlistDelete(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistReset(SCIP* scip, LIST* list, SCIP_Bool delete_data);
SCIP_Bool SCIPlistDeleteData(SCIP* scip, LIST* list);
SCIP_Bool SCIPlistResetNested(SCIP* scip, LIST* list, SCIP_Bool delete_data);
SCIP_Bool SCIPlistDeleteNested(SCIP* scip, LIST* list);
ITERATOR SCIPlistFind(ITERATOR first, ITERATOR last, void* value, SCIP_Bool (*comp_func)(void* a, void* b));
SCIP_Bool SCIPlistMove(ITERATOR position, ITERATOR it);
void SCIPlistMoveFront(ITERATOR position);
void SCIPlistMoveBack(ITERATOR position);
SCIP_RETCODE SCIPlistRearrange(SCIP* scip, LIST* list, LIST* order);
SCIP_Bool SCIPlistForeach(LIST* list, int(*func)(void*));
void SCIPlistPrint(LIST* list, int(*printfunc)(void*));

/*node*/
NODE* SCIPnodeCreate(SCIP* scip, void* data);

/*iteratos*/
ITERATOR SCIPiteratorBegin(LIST* list);
ITERATOR SCIPiteratorEnd(LIST* list);
SCIP_Bool SCIPiteratorNext(ITERATOR* it);
SCIP_Bool SCIPiteratorPrev(ITERATOR* it);
void SCIPiteratorEquals(ITERATOR* it1, ITERATOR* it2);
SCIP_Bool SCIPiteratorIsEqual(ITERATOR it1, ITERATOR it2);

/*print, compare callback*/
int printstring(void *s);
int printint(void *i);
static int compare (const void * a, const void * b);
static SCIP_Bool compare_int(void* a, void * b);



/** creates an empty lists and returns a pointer to that list */
LIST* SCIPlistCreate(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   LIST* list;
   NODE* node;

   assert( scip != NULL);

   if( SCIPallocBlockMemory(scip, &list) == SCIP_NOMEMORY )
      return NULL;

   node = SCIPnodeCreate(scip, NULL);
   list->nil = node;

   /* list->nil->next contains a pointer to the last node */
   list->nil->next = list->nil;

   /* list->nil->prev contains a pointer to the first node */
   list->nil->prev = list->nil;
   list->nil->data = node->data;
   list->size = 0;
   return list;
}

/** creates a list with integers running from 'from' to 'to'. */
LIST* SCIPlistCreateInt(
   SCIP*                 scip,               /**< SCIP data structure */
   int                   from,               /**< Start index */
   int                   to                  /**< End index */
   )
{
   LIST* list;
   int* data;
   int i;
   SCIP_Bool success;

   list = SCIPlistCreate(scip);
   for( i = from; i <= to; ++i )
   {
      if( SCIPallocMemory(scip, &data) == SCIP_NOMEMORY )
         return NULL;

      *data=i;
      success = SCIPlistPushBack(scip, list, (void*) data);
      if( !success )
         return NULL;
   }
   return list;
}

/** creates a copy of 'list', but does not copy the data stored in list->data. */
LIST* SCIPlistCopyShallow(SCIP* scip, LIST* list)
{
   LIST* copy;
   ITERATOR it;
   //case: list exists
   if( list )
   {
      copy = SCIPlistCreate(scip);
      for( it = SCIPiteratorBegin(list); ! SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)); SCIPiteratorNext(&it) )
      {
         SCIPlistPushBack(scip, copy, it.node->data);
      }
      return copy;
   }
   //case: list does not exist
   else
   {
      return NULL;
   }
}

/** inserts a new node with a pointer to 'data' at the front of the list. */
SCIP_Bool SCIPlistPushFront(SCIP* scip, LIST* list, void* data)
{
   /** adds a node at the front of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return FALSE;
   //it points to the beginning of list
   it = SCIPiteratorBegin(list);
   SCIPlistInsert(scip, &it, data);
   return TRUE;
}

/** removes the node at the front of the list. */
SCIP_Bool SCIPlistPopFront(SCIP* scip, LIST* list)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (SCIPlistIsEmpty(list)) ) return FALSE;

   it = SCIPiteratorBegin(list);
   SCIPlistErase(scip, &it);
   return TRUE;
}

/** inserts a new node with a pointer to 'data' at the back of the list. */
SCIP_Bool SCIPlistPushBack(SCIP* scip, LIST* list, void* data)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return FALSE;
   //it points to the end of list
   it = SCIPiteratorEnd(list);
   SCIPlistInsert(scip, &it, data);
   return TRUE;
}

/** removes the node at the back of the list. */
SCIP_Bool SCIPlistPopBack(SCIP* scip, LIST* list)
{
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (SCIPlistIsEmpty(list)) ) return FALSE;

   it = SCIPiteratorEnd(list);
   SCIPiteratorPrev(&it);
   SCIPlistErase(scip, &it);
   return TRUE;
}

/** creates a new node before 'it->node' and returns the new node. */
NODE* SCIPlistInsert(SCIP* scip, ITERATOR* it, void* data)
{
   NODE* newnode;
   NODE* successor;
   NODE* predecessor;
   //case: 'it' points to no node or list
   if(it->node == NULL && it->list == NULL) return NULL;
   newnode=SCIPnodeCreate(scip, data);
   successor = it->node;
   predecessor = it->node->prev;
   //adjust all pointers
   //note: due to sentinel node no distinction of cases is necessary(begin, middle, end)
   predecessor->next = newnode;
   newnode->prev = predecessor;
   newnode->next = successor;
   successor->prev = newnode;

   ++it->list->size;
   return newnode;
}

/** assigns new data to the node the iterator 'it' points to. */
SCIP_Bool SCIPlistAssign(ITERATOR* it, void* data)
{
   if( it && it->node )
   {
      it->node->data=data;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** if target and origin have the same size, this function assigns the data pointers of origin to target. */
SCIP_Bool SCIPlistAssignList(LIST* target, LIST* origin)
{
   ITERATOR it1;
   ITERATOR it2;
   if( target && origin && target->size == origin->size )
   {
      for( it1 = SCIPiteratorBegin(target), it2 = SCIPiteratorBegin(origin); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(target)); SCIPiteratorNext(&it1), SCIPiteratorNext(&it2) )
      {
         SCIPlistAssign(&it1, it2.node->data);
      }
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** creates a node with a void pointer to data and returns the node. */
NODE* SCIPnodeCreate(SCIP* scip, void* data)
{
   NODE* node;
   if( SCIPallocBlockMemory(scip, &node) == SCIP_NOMEMORY )
      return NULL;
   node->data = data;
   node->prev = NULL;
   node->next = NULL;
   return node;
}

/** removes the node the iterator 'it' points to from the list.
 *
 * Returns TRUE if 'node' is in 'list'. Returns FALSE if 'node' is not in 'list'.
 *
 * The memory the data pointer points to is not deallocated.*/
SCIP_Bool SCIPlistErase(SCIP* scip, ITERATOR* it)
{
   NODE* successor;
   NODE* predecessor;
   //case: iterator, node or list do not exist or empty list
   if( (! it) || (! it->node) || (! it->list) || SCIPlistIsEmpty(it->list) )
   {
      return FALSE;
   }
   successor = it->node->next;
   predecessor = it->node->prev;
   predecessor->next = successor;
   successor->prev = predecessor;
   SCIPfreeBlockMemory(scip, &it->node);
   it->node = successor;
   --it->list->size;
   return TRUE;
}

/** deletes the entire list, but not the data stored in the list.
 *
 *  For deallocating memory of the list data, call 'list_deleta_data' before. */
SCIP_Bool SCIPlistDelete(SCIP* scip, LIST* list)
{
   if( ! list )
   {
      return FALSE;
   }
   else
   {
      while( ! SCIPlistIsEmpty(list) )
      {
         SCIPlistPopBack(scip, list);
      }
      //deallocate sentinel node
//      SCIPfreeMemory(scip, &list->nil->data);
      SCIPfreeBlockMemory(scip, &list->nil);
      SCIPfreeBlockMemory(scip, &list);
      list = NULL;
      return TRUE;
   }
}

/** resets the list, effectively removing all nodes, resulting in an empty list.
 *
 *  For deallocating memory of the list data, set delete_data = TRUE. */
SCIP_Bool SCIPlistReset(SCIP* scip, LIST* list, SCIP_Bool delete_data)
{
   if( ! list )
   {
      return FALSE;
   }
   else
   {
      if( delete_data )
      {
         SCIPlistDeleteData(scip, list);
      }
      while( ! SCIPlistIsEmpty(list) )
      {
         SCIPlistPopBack(scip, list);
      }
      return TRUE;
   }
}


/** deallocates the memory the data pointers of the list points to. */
SCIP_Bool SCIPlistDeleteData(SCIP* scip, LIST* list)
{
   ITERATOR it;
   if( ! list )
   {
      return FALSE;
   }
   else
   {
      for( it = SCIPiteratorBegin(list); ! ( SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it) )
      {
         SCIPfreeMemory(scip, &it.node->data);
      }
      return TRUE;
   }
}

/** deallocates all memory for a nested list including the data. */
SCIP_Bool SCIPlistDeleteNested(SCIP* scip, LIST* list)
{
   ITERATOR it1;
   if( ! list )
   {
      return FALSE;
   }
   else
   {
      for( it1 = SCIPiteratorBegin(list);  ! ( SCIPiteratorIsEqual(it1, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it1) )
         {
            SCIPlistDeleteData(scip, it1.node->data);
            SCIPlistDelete(scip, it1.node->data);
         }
      SCIPlistDelete(scip, list);
      return TRUE;
   }
}

/** resets the nested list, effectively removing all nodes, resulting in an empty list.
 *
 *  For deallocating memory of the list data, set delete_data = TRUE. */
SCIP_Bool SCIPlistResetNested(SCIP* scip, LIST* list, SCIP_Bool delete_data)
{
   ITERATOR it1;
   if( ! list )
   {
      return FALSE;
   }
   else
   {
      for( it1 = SCIPiteratorBegin(list);  ! ( SCIPiteratorIsEqual(it1, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it1) )
         {
            if( delete_data )
            {
               SCIPlistDeleteData(scip, it1.node->data);
            }
            SCIPlistDelete(scip, it1.node->data);
         }
      SCIPlistReset(scip, list, FALSE);
      return TRUE;
   }
}

/** returns TRUE iff list exists and is empty. */
SCIP_Bool SCIPlistIsEmpty(LIST* list)
{
   //list is empty if the next pointer of the sentinel points to itself
   if( list && list->nil->next == list->nil )
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** runs the function 'func' for every node of the list. */
SCIP_Bool SCIPlistForeach(LIST* list, int(*func)(void*))
{
   ITERATOR it;
   for( it = SCIPiteratorBegin(list); ! ( SCIPiteratorIsEqual(it, SCIPiteratorEnd(list)) ); SCIPiteratorNext(&it) )
   {
      if( func(it.node->data) != 0 )
      {
         return FALSE;
      }
   }
   return TRUE;
}

/** prints the contents of the list. */
void SCIPlistPrint(
      LIST* list,             /**< list to print */
      int(*printfunc)(void*)  /**< function to print contents of the list */
      )
{
   printf("( ");
   SCIPlistForeach(list, printfunc);
   printf(")\n");
}

/** searches in the interval [first, last) for 'value'.
 * @param first position to start the search
 * @param last position before the search is stopped
 * @param value value which shall be found
 * @param comp_func a predicate function
 * @return Iterator pointing to the position of the first appearance of 'value'. Iterator points to the position after the end of the list, if 'value' is not found. */
ITERATOR SCIPlistFind(ITERATOR first, ITERATOR last, void* value, SCIP_Bool (*comp_func)(void* a, void* b))
{
   for (;
        //first != last AND end of list is not reached yet
        ! (SCIPiteratorIsEqual(first, last)) && ! (SCIPiteratorIsEqual(first, SCIPiteratorEnd(first.list)));
        SCIPiteratorNext(&first)
        )
   {
      if ( comp_func(value, first.node->data) )
      {
         break;
      }
   }
   return first;
}

/** moves the element at 'position' to the location in front of element 'it'
 *
 * Both elements must be in the same list. */
SCIP_Bool SCIPlistMove(ITERATOR position, ITERATOR it)
{
   NODE* old_successor;
   NODE* old_predecessor;
   NODE* new_successor;
   NODE* new_predecessor;
   if( position.node && it.node && position.list == it.list )
   {
      //case 'it' and 'position' point to the same node or 'it' points to the succeeding node of 'position': no moving
      if( position.node == it.node || position.node->next == it.node )
      {
         return TRUE;
      }
      //all other cases
      //adjust pointers at the old position
      old_successor = position.node->next;
      old_predecessor = position.node->prev;
      old_predecessor->next = old_successor;
      old_successor->prev = old_predecessor;

      //adjust pointers at the new position
      new_successor = it.node;
      new_predecessor = it.node->prev;
      new_predecessor->next = position.node;
      position.node->prev = new_predecessor;
      position.node->next = new_successor;
      new_successor->prev = position.node;

      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** moves the element at 'position' to the front of the list. */
void SCIPlistMoveFront(ITERATOR position)
{
   SCIPlistMove(position, SCIPiteratorBegin(position.list));
}

/** moves the element at 'position' to the end of the list. */
void SCIPlistMoveBack(ITERATOR position)
{
   SCIPlistMove(position, SCIPiteratorEnd(position.list));
}

/** rearranges elements of list according to the ordering of order.
 *
 * example: list = (a b c d); order = (3 2 4 1)
 * after calling SCIPlistRearrange(list, order): list = (c b d a)
 * both lists must have the same size
 * order must have elements from 1 to list->size */
SCIP_RETCODE SCIPlistRearrange(SCIP* scip, LIST* list, LIST* order)
{
   LIST* new_list;
   ITERATOR it1, it2;
   int i;
   if( list && order && list->size == order->size )
   {
      new_list = SCIPlistCreate(scip);
      for( it1 = SCIPiteratorBegin(order); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(order)); SCIPiteratorNext(&it1) )
      {
         for( it2 = SCIPiteratorBegin(list), i = 1; i < *(int*)it1.node->data; ++i )
         {
            SCIPiteratorNext(&it2);
         }
         SCIPlistPushBack(scip, new_list, it2.node->data);
      }
      SCIPlistAssignList(list, new_list);
      SCIPlistDelete(scip, new_list);
      return SCIP_OKAY;
   }
   else
   {
      return SCIP_ERROR;
   }
}

/** print function for strings.
 *
 * Can be used as a parameter in SCIPlistPrint(LIST* list, int(*printfunc)(void*)). */
int printstring(void *s)
{
   printf("%s\n", (char *)s);
   return 0;
}

/** print function for integers.
 *
 * Can be used as a parameter in SCIPlistPrint(LIST* list, int(*printfunc)(void*)). */
int printint(void *i)
{
   printf("%i ", *(int*)i);
   return 0;
}

/** returns an iterator pointing to the first element of list. */
ITERATOR SCIPiteratorBegin(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if( list )
   {
      it.list = list;
      it.node = list->nil->next;
   }
   return it;
}

/** returns an iterator pointing to the past-the-end element of list. */
ITERATOR SCIPiteratorEnd(LIST* list)
{
   //the past-the-end element of list is the sentinal node 'nil'
   ITERATOR it = {NULL, NULL};
   if( list )
   {
      it.list = list;
      it.node = list->nil;
   }
   return it;
}

/** sets the iterator 'it' to the next element. */
SCIP_Bool SCIPiteratorNext(ITERATOR* it)
{
   if( it && it->list && it->node )
   {
      it->node = it->node->next;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** sets the iterator 'it' to the previous element. */
SCIP_Bool SCIPiteratorPrev(ITERATOR* it)
{
   if( it && it->list && it->node )
   {
      it->node = it->node->prev;
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/** assigns it1 = it2. */
void SCIPiteratorEquals(ITERATOR* it1, ITERATOR* it2)
{
   it1->list = it2->list;
   it1->node = it2->node;
}

/** returns TRUE if it1 points to the same element as it2. */
SCIP_Bool SCIPiteratorIsEqual(ITERATOR it1, ITERATOR it2)
{
   if( it1.list == it2.list && it1.node == it2.node )
   {
      return TRUE;
   }
   else
   {
      return FALSE;
   }
}

/*
 * Data structures
 */

/** A struct that contains 4 hashmaps, which maps variables and constraints to their position in the constraint matrix (Ax<=b) and vice versa */
struct IndexMap
{
   /** index in problem -> constraint */
   SCIP_HASHMAP* indexcons;
   /** constraint -> index in problem */
   SCIP_HASHMAP* consindex;
   /** index in problem -> variable */
   SCIP_HASHMAP* indexvar;
   /** variable -> index in problem */
   SCIP_HASHMAP* varindex;
};
typedef struct IndexMap INDEXMAP;

/** detector data */
struct DEC_DetectorData
{
//   DEC_DECOMP* decdecomp;
   SCIP_VAR*** varsperblock;
   int* nvarsperblock;
   SCIP_CONS*** consperblock;
   int *nconsperblock;
   SCIP_VAR** linkingvars;
   int nlinkingvars;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   SCIP_HASHMAP* vartoblock;
   SCIP_HASHMAP* constoblock;
   int blocks;
   int maxblocks;
   int minblocks;
   SCIP_CONS** relevantConss; //array with all non-empty constraints
   int nRelevantConss;        //number of relevants constraints
   INDEXMAP* indexmap;
   int* ibegin; //array, ibegin[i]: index of first nonzero entry in row i
   int* iend;   //array, iend[i]: index of last nonzero entry in row i
   int* jbegin; //array, jbegin[j]: index of first nonzero entry in column j
   int* jend;   //array, jend[j]: index of last nonzero entry in column j
   int* jmin;   //array, jmin[i]: index of first nonzero column of the i-th row
   int* jmax;   //array, jmax[i]: the last nonzero entry among all rows prior to and including the i-th row
   int* minV;   //array, minV[i]: number of linking variables corresponding to a partitioning after the i-th row
   int* width;  //array, width[i]: width of the band (of nonzero entries after ROC) at row i
   int* hashmapindices;  //array with integers running from 0 to maximum(nvars, ncons)+1 (for usage of hash maps)
   LIST* rowsWithConstrictions;
   LIST* blockedAfterrow;
   SCIP_CLOCK* clock;
   SCIP_Bool found;
   int desiredblocks;
   SCIP_Bool enableblockingdynamic;  //Enable blocking type 'dynamic'
   SCIP_Bool enableblockingstatic;  //Enable blocking type 'static'
   SCIP_Bool enableblockingassoonaspossible;  //Enable blocking type 'as soon as possible'
   SCIP_Bool enablemultipledecomps;  //Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks
   int maxiterationsROC;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

//debugging methods
/** prints out detailed information on the contents of detectordata*/
#ifdef SCIP_DEBUG
static
void PrintDetectordata(
   SCIP*      scip,                 /**< SCIP data structure */
   DEC_DETECTORDATA* detectordata   /**< detectordata instance */
   )
{
   int i;
   int j;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   SCIPinfoMessage(scip, NULL, "================DETECTORDATA============\n");
   SCIPinfoMessage(scip, NULL, "# blocks: %i\n", detectordata->blocks);
   for( i = 0; i < detectordata->blocks; ++i )
   {
      SCIPinfoMessage(scip, NULL, "Block #%i (#vars: %i, #conss: %i):\n", i+1, detectordata->nvarsperblock[i], detectordata->nconsperblock[i]);
      SCIPinfoMessage(scip, NULL, "Variables (block, index):\n");
      for( j = 0; j < detectordata->nvarsperblock[i]; ++j )
      {
         var = detectordata->varsperblock[i][j];
         SCIPinfoMessage(scip, NULL, "\t%s (%i, %i)\n", SCIPvarGetName(var), *(int*) SCIPhashmapGetImage(detectordata->vartoblock, (void*) var), *(int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) var));
      }
      SCIPinfoMessage(scip, NULL, "Constraints:\n");
      for( j = 0; j < detectordata->nconsperblock[i]; ++j )
      {
         cons = detectordata->consperblock[i][j];
         SCIPinfoMessage(scip, NULL, "\t%s (%i, %i)\n", SCIPconsGetName(cons), *(int*) SCIPhashmapGetImage(detectordata->constoblock, (void*) cons), *(int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons));
      }
      SCIPinfoMessage(scip, NULL, "========================================\n");
   }
   SCIPinfoMessage(scip, NULL, "Linking variables #%i (varindex) :\n", detectordata->nlinkingvars);
   for( j = 0; j < detectordata->nlinkingvars; ++j )
   {
      var = detectordata->linkingvars[j];
      SCIPinfoMessage(scip, NULL, "\t%s (%i)\n", SCIPvarGetName(var), *(int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) var));
   }
   SCIPinfoMessage(scip, NULL, "========================================\n");
   SCIPinfoMessage(scip, NULL, "Linking constraints #%i (consindex) :\n", detectordata->nlinkingconss);
   for( j = 0; j < detectordata->nlinkingconss; ++j )
   {
      cons = detectordata->linkingconss[j];
      SCIPinfoMessage(scip, NULL, "\t%s (%i)\n", SCIPconsGetName(cons), *(int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons));
   }
   SCIPinfoMessage(scip, NULL, "========================================\n");
}

static void printArray(int* array, int size, const char* name)
{
   int i;
   printf("%s=[ ", name);
   for( i = 0; i < size; ++i )
   {
      printf("%i ", array[i]);
   }
   printf("]\n");
}

static void printNested(LIST* list, const char* name)
{
   ITERATOR it2;
   printf("%s=( ", name);
   for( it2 = SCIPiteratorBegin(list); ! SCIPiteratorIsEqual(SCIPiteratorEnd(list), it2); SCIPiteratorNext(&it2) )
   {
      SCIPlistPrint(it2.node->data, printint);
   }
   printf(")\n");
}
#endif

/** allocates memory for an indexmap. */
static
SCIP_RETCODE indexmapCreate(
      SCIP* scip,          /**< SCIP data structure  */
      INDEXMAP** indexmap, /**< address of the pointer which shall store the index map*/
      int nconss,          /**< number of constraints */
      int nvars            /**< number of variables */
      )
{
   INDEXMAP* imap;
   assert(scip != NULL);
   assert(nconss > 0);
   assert(nvars > 0);
   SCIP_CALL( SCIPallocMemory(scip, &imap) );
   assert(imap != NULL);

   SCIP_CALL( SCIPhashmapCreate(&imap->indexvar, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&imap->varindex, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&imap->indexcons, SCIPblkmem(scip), nconss) );
   SCIP_CALL( SCIPhashmapCreate(&imap->consindex, SCIPblkmem(scip), nconss) );

   *indexmap = imap;
   return SCIP_OKAY;
}

/** deallocates memory of indexmap. */
static
void indexmapFree(SCIP* scip, INDEXMAP* indexmap)
{
   SCIPhashmapFree(&indexmap->indexvar);
   SCIPhashmapFree(&indexmap->varindex);
   SCIPhashmapFree(&indexmap->indexcons);
   SCIPhashmapFree(&indexmap->consindex);
   SCIPfreeMemory(scip, &indexmap);
}

static
void indexmapInit(INDEXMAP* indexmap, SCIP_VAR** vars, int nvars, SCIP_CONS** conss, int nconss, int* hashmapindices)
{
   int i;
   int* hashmapindex;
   SCIP_VAR* var;
   SCIP_CONS* cons;
   for( i = 0; i < nvars; ++i )
   {
      var = vars[i];
      //careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = hashmapindices + i+1;
      assert( ! SCIPhashmapExists(indexmap->indexvar, (void*) hashmapindex));
      SCIPhashmapInsert(indexmap->indexvar, (void*) hashmapindex, (void*) var);
      assert( ! SCIPhashmapExists(indexmap->varindex, (void*) var));
      SCIPhashmapInsert(indexmap->varindex, (void*) var, (void*) hashmapindex);
   }
   for( i = 0; i < nconss; ++i )
   {
      cons = conss[i];
      //careful: hashmapindex+1, because '0' is treated as an empty hashmap entry, which causes an error
      hashmapindex = hashmapindices + i+1;
      assert( ! SCIPhashmapExists(indexmap->indexcons, (void*) hashmapindex));
      SCIPhashmapInsert(indexmap->indexcons, (void*) hashmapindex, (void*) cons);
      assert( ! SCIPhashmapExists(indexmap->consindex, (void*) cons));
      SCIPhashmapInsert(indexmap->consindex, (void*) cons, (void*) hashmapindex);
   }
}

/** predicate function for sorting arrays.
 *
 * Returns the value of a - b (after casting into integers). */
static
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

/** predicate function to compare two integers.
 *
 * Returns true if the value of a is equal to the value of b (a==b). */
static
SCIP_Bool compare_int(void* a, void * b)
{
   return ( *(int*)a == *(int*)b ? TRUE : FALSE );
}

/** returns the maximum of a and b: max(a,b). */
static
int maximum(int a, int b)
{
   return (a > b ? a : b);
}

/** returns the value of the maximum in the array a.
 *
 *  Or 0 if a is empty or invalid.*/
static
int maxArray(int* a, int num_elements)
{
   int i;
   int max;
   if( num_elements > 0 && a != NULL )
   {
      max = a[0];
      for (i = 1; i<num_elements; i++)
      {
         if (a[i] > max)
         {
            max = a[i];
         }
      }
      return(max);
   }
   //case: empty array
   else
   {
      return 0;
   }
}

/** returns the value of the minimum in the array a.
 *
 *  Or 0 if a is empty or invalid. */
static
int minArray(int* a, int num_elements)
{
   int i;
   int min;
   if( num_elements > 0 && a != NULL )
   {
      min = a[0];
      for (i = 1; i<num_elements; i++)
      {
         if (a[i] < min)
         {
            min = a[i];
         }
      }
      return(min);
   }
   //case: empty array
   else
   {
      return 0;
   }
}

#ifndef NDEBUG
#ifdef SCIP_DEBUG
/** returns the value of the minimum in the list between the iterators it1 and it2
 *
 *  Or -1 if a is empty or invalid. */
static
int minList(ITERATOR first, ITERATOR last)
{
   int min;
   //first is valid
   if( first.list && first.node && ! SCIPlistIsEmpty(first.list) )
   {
      min = *(int*)first.node->data;
      SCIPiteratorNext(&first);
      for (;
           //first != last AND end of list is not reached yet
           ! (SCIPiteratorIsEqual(first, last)) && ! (SCIPiteratorIsEqual(first, SCIPiteratorEnd(first.list)));
           SCIPiteratorNext(&first)
           )
      {
         if ( *(int*)first.node->data < min )
         {
            min = *(int*)first.node->data;
         }
      }
      return min;
   }
   //first invalid
   else
   {
      return -1;
   }
}
#endif
#endif

/** switches the data the pointers p1 and p2 points to. */
static
void switchPointers(void** p1, void** p2)
{
   void* p3; /* local for swap */
    p3 = *p2;
    *p2 = *p1;
    *p1= p3;
}

//debug ?
#ifndef NDEBUG
/** returns the problem name without the path */
static const char* getProbNameWithoutPath(SCIP* scip)
{
   const char* pname;
   //remove '/' from problem name
   pname = strrchr(SCIPgetProbName(scip), '/');
   if( pname == NULL )
   {
      pname = SCIPgetProbName(scip);
   }
   else
   {
      pname = pname+1;
   }
   return pname;
}


static void checkConsistencyOfIndexarrays(DEC_DETECTORDATA* detectordata, int nvars)
{
   int i;
   for( i = 0; i < detectordata->nRelevantConss - 1; ++i )
   {
      assert(detectordata->ibegin[i] <= detectordata->ibegin[i+1]);
   }
   for( i = 0; i < nvars - 1; ++i )
   {
      assert(detectordata->jbegin[i] <= detectordata->jbegin[i+1]);
   }
}


//debug ?
/** creates a data and a gnuplot file for the initial problem.
 * @param scip < SCIP data structure
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
static
SCIP_RETCODE plotInitialProblem(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   int j;
   int* varindex;
   int* consindex;
   SCIP_VAR* var;
   SCIP_VAR** vars;
   int nvars;
   SCIP_CONS* cons;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotProblem!\n");
   }
   else
   {
      for( i = 0; i < detectordata->nRelevantConss; ++i )
      {
         cons = detectordata->relevantConss[i];
         consindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons);
         assert(consindex != NULL);
         //Get array of variables from constraint
         nvars = SCIPgetNVarsXXX(scip, cons);
         SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
         SCIP_CALL( SCIPgetVarsXXX(scip, cons, vars, nvars) );
         for( j = 0; j < nvars; ++j )
         {
            var = vars[j];
            varindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) var);
            assert(varindex != NULL);
            fprintf(output, "%i %i\n", *varindex, *consindex);
         }
         SCIPfreeBufferArray(scip, &vars);
      }
   }
   fclose(output);

   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nunset xtics\nunset ytics\nunset border\nset pointsize 0.05\nset xrange [0:%i]\nset yrange[%i:0]\nplot '%s' lt 0 pt 5 notitle", pdffile, SCIPgetNVars(scip), detectordata->nRelevantConss, datafile);
   fclose(output);
   return SCIP_OKAY;
}

//debug ?
/** creates a data and a gnuplot file for the blocked problem.
 * @param scip < SCIP data structure
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
static SCIP_RETCODE plotBlocking(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   int j;
   int k;
   int nvars;
   int* varindex;
   int* consindex;
   SCIP_VAR** vars;
   SCIP_CONS* cons;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotBlocking!\n");
   }
   else
   {
      //loop over all blocks
      for( i = 0; i < detectordata->blocks; ++i )
      {
         //loop over all constraints in block
         for( j = 0; j < detectordata->nconsperblock[i]; ++j )
         {
            cons = detectordata->consperblock[i][j];
            consindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->consindex, (void*) cons);
            assert(consindex != NULL);
            nvars = SCIPgetNVarsXXX(scip, cons);
            SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
            SCIP_CALL( SCIPgetVarsXXX(scip, cons, vars, nvars) );
            //loop over all vars in constraint
            for( k = 0; k < nvars; ++k )
            {
               varindex = (int*) SCIPhashmapGetImage(detectordata->indexmap->varindex, (void*) vars[k]);
               assert(varindex != NULL);
               fprintf(output, "%i %i\n", *varindex, *consindex);
            }//loop over all vars in constraint
            SCIPfreeBufferArray(scip, &vars);
         }//loop over all constraints in block
         fprintf(output, "\n");
      }//loop over all blocks
   }
   fclose(output);

   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nunset xtics\nunset ytics\nunset border\nset style line 1 lt 0 lw 1 pt 5\nset style line 2 lt 9 lw 1 pt 5\nset pointsize 0.05\nset xrange [0:%i]\nset yrange[%i:0]\nplot for [i=0:%i:1] '%s' every :::i::(i+1) linestyle (i%%2+1) notitle", pdffile, SCIPgetNVars(scip), detectordata->nRelevantConss, detectordata->blocks-1, datafile);
   fclose(output);
   return SCIP_OKAY;
}


/** creates a data and a gnuplot file for the graph representing the array minV (number of linking variables).
 * @param detectordata < presolver data data structure
 * @param filename name of the output files (without any filename extension) */
static
void plotMinV(SCIP* scip, DEC_DETECTORDATA* detectordata, char* filename)
{
   FILE* output;
   char datafile[256];
   char blockingfile[256];
   char gpfile[256];
   char pdffile[256];
   int i;
   ITERATOR it1;
   //filenames
   sprintf(datafile, "%s.dat", filename);
   sprintf(blockingfile, "%s_blocked_at.dat", filename);
   sprintf(gpfile, "%s.gp", filename);
   sprintf(pdffile, "%s.pdf", filename);

   //datafile
   output = fopen(datafile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotMinV!\n");
   }
   else
   {
      //write data to datafile
      for( i = 0; i < detectordata->nRelevantConss -1; ++i )
      {
         fprintf(output, "%i\n", detectordata->minV[i]);
      }
   }
   fclose(output);

   //blocking points
   output = fopen(blockingfile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for blocking output in plotMinV!\n");
   }
   else
   {
      //write data to blockingfile
      for( it1 = SCIPiteratorBegin(detectordata->blockedAfterrow); ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(detectordata->blockedAfterrow)); SCIPiteratorNext(&it1) )
      {
         fprintf(output, "%i %i\n", *(int*)it1.node->data - 1, detectordata->minV[*(int*)it1.node->data-1]);
      }
   }
   fclose(output);
   //write Gnuplot file
   output = fopen(gpfile, "w");
   fprintf(output, "set terminal pdf\nset output \"%s\"\nset style line 1 lt 1 lc rgb \"black\"\nplot '%s' title '# verb. Variablen' ls 1 with lines, \\\n '%s' lt 0 pt 4 with points title \"Blockgrenze\"", pdffile, datafile, blockingfile);
   fclose(output);
}

#ifdef SCIP_DEBUG
static
void writeParams(SCIP* scip, DEC_DETECTORDATA* detectordata, char* paramfile, int ROC_iterations, int tau, double time)
{
   FILE* output;
   int i;
   int nvars;
   int ncons;
   int nonzeros;
   int zeros;
   float sparsity;
   int minimum_linking_vars;
   output=fopen(paramfile, "w");
   if (output == NULL)
   {
      SCIPinfoMessage(scip, NULL, "Can't open file for output in plotMinV!\n");
   }
   else
   {
      nvars = SCIPgetNOrigVars(scip);
//      ncons = SCIPgetNConss(scip);
      ncons = detectordata->nRelevantConss;
      nonzeros = 0;
      for( i = 0; i < ncons; ++i )
      {
         nonzeros += SCIPgetNVarsXXX(scip,  SCIPgetConss(scip)[i]);
      }
      zeros = nvars*ncons - nonzeros;
      sparsity = (float) nonzeros / (nvars*ncons);
      minimum_linking_vars = minList(SCIPiteratorBegin(detectordata->rowsWithConstrictions), SCIPiteratorEnd(detectordata->rowsWithConstrictions));
      SCIPdebugMessage("minList.\n");
      fprintf(output, "# of rows\n%i\n", ncons);
      fprintf(output, "# of columns\n%i\n", nvars);
      fprintf(output, "# of nonzeros\n%i\n", nonzeros);
      fprintf(output, "# of zeros\n%i\n", zeros);
      fprintf(output, "# sparsity\n%f\n", sparsity);
      fprintf(output, "# detection time in seconds\n%f\n", time);
      fprintf(output, "# tau\n%i\n", tau);
      fprintf(output, "# of blocks\n%i\n", detectordata->blocks);
      fprintf(output, "# of iterations\n%i\n", ROC_iterations);
      fprintf(output, "# of minimum linking vars\n%i\n", minimum_linking_vars);
      fprintf(output, "# of linking vars\n%i\n", detectordata->nlinkingvars);
      for( i = 0; i < detectordata->blocks; ++i )
      {
         fprintf(output, "block # %i\n", i+1);
         fprintf(output, "# nonlinking vars\n%i\n", detectordata->nvarsperblock[i]);
         fprintf(output, "# cons per block\n%i\n", detectordata->nconsperblock[i]);
      }
      fclose(output);
   }
}
#endif
#endif

/** scans all constraints of the constraint array of the scip object,
 * and stores pointers to all constraints that have at least one variable in detectordata->relevantConss.
 * Thus it removes all empty constraints.
 */
static
SCIP_RETCODE findRelevantConss(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   SCIP_CONS** cons_array;
   LIST* relevantConssIndices;
   int i;
   int* data;
   ITERATOR it1;
   cons_array = SCIPgetConss(scip);
   relevantConssIndices = SCIPlistCreate(scip);
   for( i = 0; i < SCIPgetNConss(scip); ++i )
   {
      if( SCIPgetNVarsXXX(scip, cons_array[i]) > 0 )
      {
         SCIP_CALL( SCIPallocMemory(scip, &data) );
         *data = i;
         SCIPlistPushBack(scip, relevantConssIndices, data);
      }
   }
   //debug
//   SCIPlistPrint(relevantConssIndices, printint);
   //allocate memory for detectordata->relevantConss and store pointers of relevant conss
   detectordata->nRelevantConss = relevantConssIndices->size;
   SCIPdebugMessage("nRelevantConss: %i \n", detectordata->nRelevantConss);
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->relevantConss, detectordata->nRelevantConss) );
   for( i = 0, it1 = SCIPiteratorBegin(relevantConssIndices); i < detectordata->nRelevantConss; ++i, SCIPiteratorNext(&it1) )
   {
      detectordata->relevantConss[i] = cons_array[*(int*)it1.node->data];
   }
   SCIPlistDeleteData(scip, relevantConssIndices);
   SCIPlistDelete(scip, relevantConssIndices);

   return SCIP_OKAY;
}

/** creates a nested list with the indices of the nonzero entries of each row.
 *
 * example:
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting list:
 *  ( (1 2 4)
 *    (2 3)
 *    (5)    )
 */
static
SCIP_RETCODE rowindices_list(
      SCIP* scip,                      /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      SCIP_HASHMAP* indexcons,         /**< hashmap index -> constraint */
      SCIP_HASHMAP* varindex,          /**< hashmap variable -> index*/
      LIST** rowindices                 /**< list to store the row indices list*/
      )
{
   //create the rowindices list
   int i;
   int j;
   int* data;
   LIST* rowindices_row;
   int ncons; //number of constraints of the problem
   int nvars; //number of variables in a constraint
   int* probindices;
   int* hashmapindex;
   SCIP_CONS* cons; //one constraint of the problem
   SCIP_VAR** vars; //array of variables that occur in a constraint (unequal zero)

   *rowindices = SCIPlistCreate(scip);
   ncons = detectordata->nRelevantConss;
   for( i = 0; i < ncons; ++i )
   {
      hashmapindex = &detectordata->hashmapindices[i+1];
      cons = (SCIP_CONS*) SCIPhashmapGetImage(indexcons, (void*) hashmapindex);
      nvars = SCIPgetNVarsXXX(scip, cons);
      SCIP_CALL( SCIPallocBufferArray(scip, &vars, nvars) );
      SCIPgetVarsXXX(scip, cons, vars, nvars);
      //allocate memory for the array of probindices
      SCIP_CALL( SCIPallocMemoryArray(scip, &probindices, nvars) );
      //fill the array with the indices of the variables of the current constraint
      for( j = 0; j < nvars; ++j )
      {
         probindices[j] = *(int*) SCIPhashmapGetImage(varindex, vars[j]);
      }
      //sort the elements of probindices ('<')
      qsort(probindices, nvars, sizeof(int), compare);
      //store a copy of the elements of probindices in the list rowindices_row
      rowindices_row = SCIPlistCreate(scip);
      for( j = 0; j < nvars; ++j )
      {
         SCIP_CALL( SCIPallocMemory(scip, &data) );
         *data = probindices[j];
         SCIPlistPushBack(scip, rowindices_row, data);
      }
      //deallocate memory
      SCIPfreeMemoryArray(scip, &probindices);
      SCIPfreeBufferArray(scip, &vars);
      //add rowindices_row to the list rowindices
      SCIPlistPushBack(scip, *rowindices, rowindices_row);
   }
   //debug
//   printNested(*rowindices, "rowindices in rowindices_list");
   return SCIP_OKAY;
}

/** creates a nested list with the indices of the nonzero entries of each column.
 *
 * example:
 *
 * constraint matrix:
 *
 *  1 1 0 1 0
 *
 *  0 1 1 0 0
 *
 *  0 0 0 0 1
 *
 *  resulting list:
 *  ( (1)
 *    (1 2)
 *    (2)
 *    (1)
 *    (3)    )
 */
static
SCIP_RETCODE columnindices_list(
      SCIP* scip,                      /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata,  /**< detector data data structure */
      LIST* rowindices,                /**< A list with the row indices (achieved from calling rowindices_list() ) */
      LIST** columnindices              /**< list to store the column indices list*/
      )
{
   LIST** columnindices_array;
   LIST* rowindices_row;
   int* data;
   int position;
   int nvars;
   int i;
   ITERATOR it1;
   ITERATOR it2;
   nvars = SCIPgetNVars(scip);
   //create the columnindices_array with pointers to empty lists
   SCIP_CALL( SCIPallocMemoryArray(scip, &columnindices_array, nvars) );
   for( i = 0; i < nvars; ++i )
   {
      columnindices_array[i] = SCIPlistCreate(scip);
   }

   for( it1 = SCIPiteratorBegin(rowindices), i = 0; ! ( SCIPiteratorIsEqual(it1, SCIPiteratorEnd(rowindices)) ); SCIPiteratorNext(&it1), ++i )
   {
      rowindices_row = it1.node->data;
      for( it2 = SCIPiteratorBegin(rowindices_row); ! ( SCIPiteratorIsEqual(it2, SCIPiteratorEnd(rowindices_row)) ); SCIPiteratorNext(&it2) )
      {
         SCIP_CALL( SCIPallocMemory(scip, &data) );
         *data = i+1;
         position = *(int*)(it2.node->data)-1;
         SCIPlistPushBack(scip, columnindices_array[position], data);
      }
   }
   //create a columnindices list instead of an array
   *columnindices = SCIPlistCreate(scip);
   for( i = 0; i < nvars; ++i )
   {
      SCIPlistPushBack(scip, *columnindices, columnindices_array[i]);
   }
   //deallocate memory
   SCIPfreeMemoryArray(scip, &columnindices_array);
   return SCIP_OKAY;
}

/** does the row ordering of the ROC2 algorithm.
 *
 * It also works for the column ordering. In this case the terms row<->column have to be exchanged.
 *
 * @param columnindices A list of the nonzero entries in each column.
 * @param nrows The number of rows of the constraint matrix (=number of relevant constraints)
 * @return A list with the new row order. E.g. (2 3 1) means the second row comes first now, and so on. */
static
LIST* rowOrdering(SCIP* scip, LIST* columnindices, int nrows)
{
   LIST* roworder;
   LIST* new_roworder;
   ITERATOR it1;
   ITERATOR it2;
   ITERATOR it3;
   ITERATOR it4;

   //create a list for the order of the rows ( 1 2 3 ... nrows )
   roworder = SCIPlistCreateInt(scip, 1, nrows);
   new_roworder = SCIPlistCopyShallow(scip, roworder);

   for( it1 = SCIPiteratorEnd(columnindices), SCIPiteratorPrev(&it1); it1.node != it1.list->nil; SCIPiteratorPrev(&it1) )
   {
      for( it2 = SCIPiteratorEnd(roworder), SCIPiteratorPrev(&it2); it2.node != it2.list->nil; SCIPiteratorPrev(&it2) )
      {
         it3 = SCIPlistFind(SCIPiteratorBegin(it1.node->data), SCIPiteratorEnd(it1.node->data), it2.node->data, compare_int);
         if( it3.node->data != NULL )
         {
            it4 = SCIPlistFind(SCIPiteratorBegin(new_roworder), SCIPiteratorEnd(new_roworder), it2.node->data, compare_int);
            SCIPlistMoveFront(it4);
         }
         else
         {
            ;
         }
      }
      SCIPlistAssignList(roworder, new_roworder);
   }
   //deallocate memory
   SCIPlistDelete(scip, new_roworder);
   return roworder;
}

/** stores the first and last entry of the i-th column(row) in begin[i] and end[i] respectively.
 *
 * @param begin Array to store the first nonzero entry of the i-th column (row)
 * @param end Array to store the last nonzero entry of the i-th column (row)
 * @param indices columnindices list (rowindices list) */
static
SCIP_RETCODE formIndexArray(int* begin, int* end, LIST* indices)
{
   ITERATOR it1;
   ITERATOR it2;
   int i;
   assert(begin != NULL && end != NULL && indices != NULL);
   for( it1 = SCIPiteratorBegin(indices), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(indices)); SCIPiteratorNext(&it1), ++i )
   {
      //case: list not empty
      if( ! SCIPlistIsEmpty(it1.node->data) )
      {
         it2 = SCIPiteratorBegin(it1.node->data);
         begin[i] = *(int*)it2.node->data;
         it2 = SCIPiteratorEnd(it1.node->data);
         SCIPiteratorPrev(&it2);
         end[i] = *(int*)it2.node->data;
      }
      //case: list empty
      else
      {
         begin[i] = 0;
         end[i] = 0;
      }
   }
   return SCIP_OKAY;
}


/**returns FALSE if at least one entry of new_array and old_array are different.*/
static
SCIP_Bool arraysAreEqual(int* new_array, int* old_array, int num_elements)
{
   int i;
   for( i = 0; i < num_elements; ++i )
   {
      if( new_array[i] != old_array[i] )
      {
         return FALSE;
      }
   }
   //case: all entries of old and new are equal
   return TRUE;
}

/**permutes the order of rows and columns in inputmap and stores the result in outputmap.
 *
 *  One call of this function is equivalent to one iteration of the ROC2-algortihm. */
static
SCIP_RETCODE rankOrderClusteringIteration(
      SCIP*             scip,          /**< SCIP data structure */
      DEC_DETECTORDATA*  detectordata, /**< presolver data data structure */
      INDEXMAP* inputmap,              /**< indexmap for input */
      INDEXMAP* outputmap              /**< indexmap for output */
      )
{
   LIST* roworder;
   LIST* columnorder;
   LIST* rowindices;
   LIST* columnindices;
   ITERATOR it1;
   int nvars;
   int ncons;
   int i;
   int position;
   int* hashmapindex;
   SCIP_CONS* cons;
   SCIP_VAR* var;

   SCIPdebugMessage("Entering rankOrderClusteringIteration\n");

   assert(scip != NULL);
   assert(detectordata != NULL);
   nvars = SCIPgetNVars(scip);
   ncons = detectordata->nRelevantConss;
   //create the lists containing the positions of nonzero entries; row and column ordering
   rowindices = SCIPlistCreate(scip);
   rowindices_list(scip, detectordata, inputmap->indexcons, inputmap->varindex, &rowindices);
   columnindices = SCIPlistCreate(scip);
   columnindices_list(scip, detectordata, rowindices, &columnindices);
   roworder = rowOrdering(scip, columnindices, ncons);
   SCIPlistRearrange(scip, rowindices, roworder);
   columnorder = rowOrdering(scip, rowindices, nvars);

   //consindex and indexcons
   for( it1 = SCIPiteratorBegin(roworder), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(roworder)) && i < ncons; ++i, SCIPiteratorNext(&it1) )
   {
      position = (*(int*)it1.node->data);
      hashmapindex = &detectordata->hashmapindices[position];
      cons = SCIPhashmapGetImage(inputmap->indexcons, (void*) hashmapindex);
      assert ( cons != NULL);
      //consindex
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->consindex, (void*) cons));
      SCIPhashmapSetImage(outputmap->consindex, (void*) cons, (void*) hashmapindex);
      //indexcons
      assert( SCIPhashmapExists(outputmap->indexcons, (void*) hashmapindex ));
      SCIPhashmapSetImage(outputmap->indexcons, (void*) hashmapindex, cons);
   }
   //varindex and indexvar
   for( it1 = SCIPiteratorBegin(columnorder), i = 0; ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(columnorder)) &&i < nvars; ++i, SCIPiteratorNext(&it1) )
   {
      position = (*(int*)it1.node->data);
      hashmapindex = &detectordata->hashmapindices[position];
      var = (SCIP_VAR*) SCIPhashmapGetImage(inputmap->indexvar, (void*) hashmapindex);
      assert ( var != NULL);

      //varindex
      hashmapindex = &detectordata->hashmapindices[i+1];
      assert( SCIPhashmapExists(outputmap->varindex, (void*) var) );
      SCIPhashmapSetImage(outputmap->varindex, (void*) var, (void*) hashmapindex);
      //indexvar
      assert( SCIPhashmapExists(outputmap->indexvar, (void*) hashmapindex ));
      SCIPhashmapSetImage(outputmap->indexvar, (void*) hashmapindex, var);
   }
   //deallocate memory
   SCIPlistDeleteData(scip, roworder);
   SCIPlistDelete(scip, roworder);
   SCIPlistDeleteData(scip, columnorder);
   SCIPlistDelete(scip, columnorder);
   SCIPlistDeleteNested(scip, rowindices);
   SCIPlistDeleteNested(scip, columnindices);
   return SCIP_OKAY;
}

static
int rankOrderClustering(SCIP* scip, DEC_DETECTORDATA* detectordata, int max_iterations)
{
   int i;
   int nvars;
   int ncons;
   INDEXMAP* indexmap_permuted;
   LIST* rowindices;
   LIST* columnindices;
   int* ibegin_permuted;
   int* iend_permuted;
   int* jbegin_permuted;
   int* jend_permuted;
   assert(scip != NULL);
   assert(detectordata != NULL);

   if( max_iterations <= 0 )
   {
      return max_iterations;
   }
   else
   {
      nvars = SCIPgetNVars(scip);
      ncons = detectordata->nRelevantConss;
      indexmapCreate(scip, &indexmap_permuted, ncons, nvars);
      SCIP_CALL( SCIPallocMemoryArray(scip, &ibegin_permuted, ncons) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &iend_permuted, ncons) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &jbegin_permuted, nvars) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &jend_permuted, nvars) );
      indexmapInit(indexmap_permuted, SCIPgetVars(scip), nvars, detectordata->relevantConss, ncons, detectordata->hashmapindices);
      rowindices = SCIPlistCreate(scip);
      columnindices = SCIPlistCreate(scip);
      i = 0;
      do
      {
         ++i;
         //not more than max_iterations loops. no iteration limit for max_iterations == -1
         if(i > max_iterations && max_iterations != -1) {break;}
         SCIPdebugMessage("Iteration # %i of ROC2\n", i);
         rankOrderClusteringIteration(scip, detectordata, detectordata->indexmap, indexmap_permuted);
         //form the new index arrays after the permutation
         rowindices_list(scip, detectordata, indexmap_permuted->indexcons, indexmap_permuted->varindex, &rowindices);
         columnindices_list(scip, detectordata, rowindices, &columnindices);
         formIndexArray(ibegin_permuted, iend_permuted, rowindices);
         formIndexArray(jbegin_permuted, jend_permuted, columnindices);
         SCIPlistResetNested(scip, rowindices, TRUE);
         SCIPlistResetNested(scip, columnindices, TRUE);
         //switch between index arrays containing new and old indices
         switchPointers( (void*) &detectordata->ibegin, (void*) &ibegin_permuted);
         switchPointers( (void*) &detectordata->iend, (void*) &iend_permuted);
         switchPointers( (void*) &detectordata->jbegin, (void*) &jbegin_permuted);
         switchPointers( (void*) &detectordata->jend, (void*) &jend_permuted);
         //switch between hash maps containing new and old indices
         switchPointers( (void*) &detectordata->indexmap, (void*) &indexmap_permuted);
      }
      //while Index Arrays change
      while( ! (arraysAreEqual(detectordata->ibegin, ibegin_permuted, ncons )
             && arraysAreEqual(detectordata->iend, iend_permuted, ncons)
             && arraysAreEqual(detectordata->jbegin, jbegin_permuted, nvars)
             && arraysAreEqual(detectordata->jend, jend_permuted, nvars)));

      indexmapFree(scip, indexmap_permuted);
      SCIPfreeMemoryArray(scip, &ibegin_permuted);
      SCIPfreeMemoryArray(scip, &iend_permuted);
      SCIPfreeMemoryArray(scip, &jbegin_permuted);
      SCIPfreeMemoryArray(scip, &jend_permuted);
   }
   return (i-1);
}

/** finds rows with local minima regarding the number of linking variables and stores them in detectordata->rowsWithConstrictions */
static
SCIP_RETCODE rowsWithConstriction(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   //if blocking is performed after row i+1; local minima
   int i;
   int* data;
   for( i = 1; i < detectordata->nRelevantConss - 2; ++i )
   {
      //is minV[i] a local minimum?    < or <=   ? What does make more sense?
      if( detectordata->minV[i] < detectordata->minV[i-1] && detectordata->minV[i] < detectordata->minV[i+1] )
      {
         SCIP_CALL( SCIPallocMemory(scip, &data) );
         *data = i+1;
         SCIPlistPushBack(scip, detectordata->rowsWithConstrictions, data);
      }
   }
   return SCIP_OKAY;
}

/** assigns variables to a block, divided into linking variables and nonlinking variables.*/
static
SCIP_RETCODE assignVarsToBlock(
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int block,                       /**< number of current block, first block = 1 */
      int first_var,                   /**< index of first variable in the block*/
      int last_var,                    /**< index of last variable in the block*/
      int first_linkingvar             /**< index of first linking variable*/
      )
{
   int i;
   int j;
   int* hashmapindex;
   SCIP_VAR* var;
   //assign the subscipvars (=nonlinking vars)
   //debug assertion not useful in case of constant block size?
//   assert(first_linkingvar >= first_var);
   detectordata->nvarsperblock[block-1] = maximum(first_linkingvar - first_var, 0);
   for( i = first_var, j = 0; j < detectordata->nvarsperblock[block-1]; ++i, ++j )
   {
      hashmapindex = &detectordata->hashmapindices[i];
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->indexmap->indexvar, (void*) hashmapindex);
      assert(var != NULL);
      detectordata->varsperblock[block-1][j] = var;
      /* insert var into hash map vartoblock */
      assert(!SCIPhashmapExists(detectordata->vartoblock, var));
      SCIP_CALL( SCIPhashmapInsert(detectordata->vartoblock, var, (void*)&detectordata->hashmapindices[block]) );
   }
   //assign linking vars
   for( i = first_linkingvar; i <= last_var; ++i )
   {
      hashmapindex = &detectordata->hashmapindices[i];
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->indexmap->indexvar, (void*) hashmapindex);
      assert(var != NULL);
      detectordata->linkingvars[detectordata->nlinkingvars] = var;
      ++detectordata->nlinkingvars;
   }
   return SCIP_OKAY;
}

/** assigns constraints in the interval [first_cons, last_cons] to 'block'. */
static
SCIP_RETCODE assignConsToBlock(SCIP* scip, DEC_DETECTORDATA* detectordata, int block, int first_cons, int last_cons)
{
   int i;
   int j;
   int* hashmapindex;
   SCIP_CONS* cons;
   //assign the constraints to the current block
   detectordata->nconsperblock[block-1] = last_cons - first_cons + 1;
   for( i = first_cons, j = 0; i <= last_cons; ++i, ++j )
   {
      hashmapindex = &detectordata->hashmapindices[i];
      cons = (SCIP_CONS*) SCIPhashmapGetImage(detectordata->indexmap->indexcons, (void*) hashmapindex);
      assert(cons != NULL);
      detectordata->consperblock[block-1][j] = cons;
      /* insert cons into hash map vartoblock */
      assert(!SCIPhashmapExists(detectordata->constoblock, cons));
      SCIP_CALL( SCIPhashmapInsert(detectordata->constoblock, cons, (void*)&detectordata->hashmapindices[block]) );
   }
   SCIPlistPushBack(scip, detectordata->blockedAfterrow, &detectordata->hashmapindices[last_cons]);
   return SCIP_OKAY;
}

/** returns the largest column index of a nonzero entry between rows [from_row, to_row] */
static
int getMaxColIndex(DEC_DETECTORDATA* detectordata, int from_row, int to_row)
{
   //some pointer arithmetic
   return maxArray(detectordata->iend + (from_row -1), to_row - from_row + 1);
}

/** returns the column index of the first nonzero entry in 'row'. Rows start counting at 1, not 0. */
static
int getMinColIndex(DEC_DETECTORDATA* detectordata, int row)
{
   return detectordata->ibegin[row-1];
}

/** determines if a blocking at 'block_at_row' is a valid blocking
 *
 * @param detectordata detectordata data structure
 * @param prev_block_first_row first row of the previous block
 * @param prev_block_last_row last row of the previous block
 * @param block_at_row the row for which you want to determine if the blocking is valid
 * @return TRUE if blocking is valid, else FALSE
 */
static
SCIP_Bool isValidBlocking(DEC_DETECTORDATA* detectordata, int prev_block_first_row, int prev_block_last_row, int block_at_row)
{
   int last_column_prev_block;
   int first_column_current_block;

   //if the function is called for the first block, the blocking is always valid
   if( prev_block_last_row == 0 )
   {
      return TRUE;
   }
   last_column_prev_block = getMaxColIndex(detectordata, prev_block_first_row, prev_block_last_row);
   first_column_current_block = getMinColIndex(detectordata, block_at_row);
   return ( first_column_current_block > last_column_prev_block ? TRUE : FALSE);
}

/** this functions looks for rows to block at, which creates block of size min_block_size or bigger
 *
 * @param it_constrictions Iterator pointing to a list of constraints (detectordata->rowsWithConstrictions)
 * @param min_block_size minimum number of rows to be in a block
 * @param prev_block_last_row the last row of the preceding block
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no candidate was found
 */
static
ITERATOR findBlockingCandidate(ITERATOR it_constrictions, int min_block_size, int prev_block_last_row)
{
   for( ;; )
   {
      //end of the list?
      if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)) )
      {
         return it_constrictions;
      }
      //does a blocking to the next row forming a constriction comprise more rows than min_block_size
      if( (* (int*) (it_constrictions.node->data) - prev_block_last_row) >= min_block_size )
      {
         return it_constrictions;
      }
      //advance iterator to next element
      SCIPiteratorNext(&it_constrictions);
   }
}

/** this functions determines the next row to block at
 *
 * @param detectordata detectordata data structure
 * @param it_constrictions Iterator pointing to a list of constraints (detectordata->rowsWithConstrictions)
 * @param min_block_size minimum number of rows to be in a block
 * @param prev_block_first_row the first row of the preceding block
 * @param prev_block_last_row the last row of the preceding block
 * @return Iterator pointing to a node which contains a suitable row for blocking; If the iterator points after the last element, no row was found
 */
static
ITERATOR nextRowToBlockAt(DEC_DETECTORDATA* detectordata, ITERATOR it_constrictions, int min_block_size, int prev_block_first_row, int prev_block_last_row)
{
   assert(it_constrictions.list != NULL);
   assert(it_constrictions.node != NULL);

   //end of the constriction list?
   if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)) )
   {
      return it_constrictions;
   }

   for( ;; )
   {
      //find a blocking candidate
      it_constrictions = findBlockingCandidate(it_constrictions, min_block_size, prev_block_last_row);
      //case: no candidate found
      if( SCIPiteratorIsEqual(it_constrictions, SCIPiteratorEnd(it_constrictions.list)) )
      {
         break;
      }
      //case: candidate found
      else
      {
         //valid blocking
         if( isValidBlocking(detectordata, prev_block_first_row, prev_block_last_row, *(int*) it_constrictions.node->data) )
         {
            break;
         }
         //invalid blocking
         else
         {
            SCIPiteratorNext(&it_constrictions);
         }
      }
   }
   return it_constrictions;
}

static
int calculateNdecompositions(DEC_DETECTORDATA* detectordata)
{
   int nblockingtypes;
   int nblockingspertype;

   nblockingtypes = 0;
   //get the number of enabled blocking types
   if( detectordata->enableblockingdynamic )
   {
      ++nblockingtypes;
   }
   if( detectordata->enableblockingstatic )
   {
      ++nblockingtypes;
   }
   if( detectordata->enableblockingassoonaspossible )
   {
      ++nblockingtypes;
   }

   //get the number of blockings per blocking type
   if( detectordata->enablemultipledecomps )
   {
      nblockingspertype = detectordata->maxblocks - detectordata->minblocks + 1;
   }
   else
   {
      nblockingspertype = 1;
   }

   return nblockingtypes * nblockingspertype;
}

static
void checkParameterConsistency(DEC_DETECTORDATA* detectordata, SCIP_RESULT* result)
{
   //maxblocks < nRelevantsCons?

   //desired blocks <= maxblocks?

   //is  minblocks <= maxblocks?
   if( detectordata->enablemultipledecomps )
   {
      if( detectordata->minblocks > detectordata->maxblocks )
      {
         SCIPerrorMessage("minblocks > maxblocks. Setting minblocks = maxblocks.\n");
         detectordata->minblocks = detectordata->maxblocks;
      }
   }

   //is at least one blocking type enabled?
   if( ! detectordata->enableblockingassoonaspossible && ! detectordata->enableblockingstatic &&! detectordata->enableblockingdynamic )
   {
      SCIPerrorMessage("No blocking type enabled, cannot perform blocking.\n");
      *result = SCIP_DIDNOTRUN;
   }
}

/** tries to dynamically divide the problem into subproblems (blocks)*/
static
SCIP_RETCODE blockingDynamic(
      SCIP* scip,                      /**< scip object */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int tau,                         /**< desired number of blocks */
      int nvars                        /**< number of variables in the problem*/
      )
{
   int block;
   int prev_block_first_row;
   int prev_block_last_row;
   int current_row;
   int min_block_size;
   //notation: i=current block; im1=i-1=previous block; ip1=i+1=next block
   int max_col_index_im1;
   int min_col_index_ip1;
   int max_col_index_i;
   ITERATOR it1;
   //debug
   SCIPdebugMessage("Starting Blocking...\n");
   SCIPdebugMessage("Max blocks: %i\n", detectordata->maxblocks);
   block = 1;
   prev_block_first_row = 0;
   prev_block_last_row = 0;
   max_col_index_im1 = 0;
   min_col_index_ip1 = 1;
   min_block_size = round( (float) detectordata->nRelevantConss / (2 * tau ));
   it1 = SCIPiteratorBegin(detectordata->rowsWithConstrictions);

   for( it1 = nextRowToBlockAt(detectordata, it1, min_block_size, prev_block_first_row, prev_block_last_row);
         ! SCIPiteratorIsEqual(it1, SCIPiteratorEnd(it1.list)) && block < detectordata->maxblocks;
         it1 = nextRowToBlockAt(detectordata, it1, min_block_size, prev_block_first_row, prev_block_last_row) )
   {
      current_row = * (int*) (it1.node->data);
      max_col_index_i = getMaxColIndex(detectordata, prev_block_last_row + 1, current_row);
      min_col_index_ip1 = getMinColIndex(detectordata, current_row + 1);
      SCIPdebugMessage("assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, current_row);
      SCIPdebugMessage("vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, max_col_index_i, min_col_index_ip1, max_col_index_i);
      //assign the variables and constraints to block
      assignVarsToBlock(detectordata, block, max_col_index_im1+1, max_col_index_i, min_col_index_ip1);
      assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, current_row);
      //update variables in the while loop
      max_col_index_im1 = max_col_index_i;
      prev_block_first_row = prev_block_last_row + 1;
      prev_block_last_row = current_row;
      ++block;
   }
   //assign the remaining (< M/2tau) cons and vars to the last block; no new linking vars are added
   //debug
   SCIPdebugMessage("last time: assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPdebugMessage("last time: vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, nvars, nvars+1, nvars);
   assignVarsToBlock(detectordata, block, max_col_index_im1+1, nvars, nvars+1);
   assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPlistPopBack(scip, detectordata->blockedAfterrow);

   detectordata->blocks = block;
   detectordata->found = TRUE;
   //debug plot the blocking  plot for [i=1:2:1] 'test.dat' every :::i::i lt i pt 5
#ifndef NDEBUG
   {
   char filename1[256];
   char filename2[256];
   char paramfile[256];

   sprintf(filename1, "%s_dynamic_blocking", getProbNameWithoutPath(scip));
   sprintf(filename2, "%s_dynamic_minV", getProbNameWithoutPath(scip));
   sprintf(paramfile, "%s_dynamic.params", getProbNameWithoutPath(scip));
   plotBlocking(scip, detectordata, filename1);
   plotMinV(scip, detectordata, filename2);
   //debug
#ifdef SCIP_DEBUG
   PrintDetectordata(scip, detectordata);
   SCIP_CALL( SCIPstopClock(scip, detectordata->clock) );
   writeParams(scip, detectordata, paramfile, ROC_iterations, tau, SCIPgetClockTime(scip, detectordata->clock));
#endif
   }
#endif

   return SCIP_OKAY;
}

/** returns the number of rows in a block in order to distribute the number of rows evenly across the blocks */
static
int rowsInConstantBlock(int block, int desired_blocks, int nrows)
{
   if( block <= desired_blocks - (nrows % desired_blocks) )
   {
      return (nrows / desired_blocks);
   }
   else
   {
      return ((nrows / desired_blocks) + 1);
   }
}

/** creates blocks with the same number of rows*/
static
SCIP_RETCODE blockingStatic(
      SCIP* scip,                      /**< scip object */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int desired_blocks,              /**< desired number of blocks */
      int nvars                        /**< number of variables in the problem*/
      )
{
   int block;
   int prev_block_last_row;
   int current_row;
   //notation: i=current block; im1=i-1=previous block; ip1=i+1=next block
   int max_col_index_im1;
   int min_col_index_ip1;
   int max_col_index_i;

   block = 1;
   prev_block_last_row = 0;
   max_col_index_im1 = 0;
   min_col_index_ip1 = 1;
   current_row = 0;
   //blocks 1 to (desired_blocks-1)
   for( block = 1; block < desired_blocks; ++block )
   {
      current_row += rowsInConstantBlock(block, desired_blocks, detectordata->nRelevantConss);
      max_col_index_i = getMaxColIndex(detectordata, prev_block_last_row + 1, current_row);
      min_col_index_ip1 = getMinColIndex(detectordata, current_row + 1);

      //first check if three adjacent blocks overlap; in this case all variables are linking
      if( min_col_index_ip1 <= max_col_index_im1 )
      {
         SCIPdebugMessage("assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, current_row);
         SCIPdebugMessage("vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, max_col_index_i, max_col_index_im1+1, max_col_index_i);
         assignVarsToBlock(detectordata, block, max_col_index_im1 + 1, max_col_index_i, max_col_index_im1 + 1);
      }
      else //no overlap of three adjacent blocks, only some vars are linking
      {
         //assign the variables and constraints to block
         SCIPdebugMessage("assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, current_row);
         SCIPdebugMessage("vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, max_col_index_i, min_col_index_ip1, max_col_index_i);
         assignVarsToBlock(detectordata, block, max_col_index_im1 + 1, max_col_index_i, min_col_index_ip1);
      }
      assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, current_row);
      //update variables in the while loop
      max_col_index_im1 = max_col_index_i;
      prev_block_last_row = current_row;
   }
   //last block
   //assign the remaining cons and vars to the last block; no new linking vars are added
   //debug
   SCIPdebugMessage("last time: assignVarsToBlock: block, from_row, to_row: %i, %i, %i\n", block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPdebugMessage("last time: vars in block: %i - %i, linking vars: %i - %i\n", max_col_index_im1+1, nvars, nvars+1, nvars);
   assignVarsToBlock(detectordata, block, max_col_index_im1+1, nvars, nvars+1);
   assignConsToBlock(scip, detectordata, block, prev_block_last_row + 1, detectordata->nRelevantConss);
   SCIPlistPopBack(scip, detectordata->blockedAfterrow);

   detectordata->blocks = block;
   detectordata->found = TRUE;
#ifndef NDEBUG
   {
   char filename1[256];
   char filename2[256];
   char paramfile[256];

   //debug
//   PrintDetectordata(scip, detectordata);
   sprintf(filename1, "%s_static_blocking_%i", getProbNameWithoutPath(scip), detectordata->blocks);
   sprintf(filename2, "%s_static_minV_%i", getProbNameWithoutPath(scip), detectordata->blocks);
   sprintf(paramfile, "%s_static.params", getProbNameWithoutPath(scip));
   plotBlocking(scip, detectordata, filename1);
   plotMinV(scip, detectordata, filename2);
//   SCIP_CALL( SCIPstopClock(scip, detectordata->clock) );
//   writeParams(scip, detectordata, paramfile, ROC_iterations, tau, SCIPgetClockTime(scip, detectordata->clock));
   }
#endif

   return SCIP_OKAY;
}

static
SCIP_RETCODE blockingAsSoonAsPossible(
      SCIP* scip,                      /**< scip object */
      DEC_DETECTORDATA* detectordata,  /**< presolver data data structure */
      int desired_blocks,              /**< desired number of blocks */
      int nvars                        /**< number of variables in the problem*/
      )
{
   int block;
   block = 0;
   detectordata->blocks = block;
   detectordata->found = TRUE;
   return SCIP_OKAY;
}

/** copies the variable and block information to the decomp structure */
static
SCIP_RETCODE copyDetectorDataToDecomp(
      SCIP*             scip,         /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata, /**< presolver data data structure */
      DEC_DECOMP*        decdecomp        /**< DECOMP data structure */
      )
{

   assert(scip != 0);
   assert(detectordata != 0);
   assert(decdecomp != 0);

   DECdecompSetNBlocks(decdecomp, detectordata->blocks);
   SCIP_CALL( DECdecompSetType(decdecomp, DEC_DECTYPE_STAIRCASE) );
   SCIP_CALL( DECdecompSetSubscipvars(scip, decdecomp, detectordata->varsperblock, detectordata->nvarsperblock) );
   SCIP_CALL( DECdecompSetSubscipconss(scip, decdecomp, detectordata->consperblock, detectordata->nconsperblock) );
   SCIP_CALL( DECdecompSetLinkingvars(scip, decdecomp, detectordata->linkingvars, detectordata->nlinkingvars) );
   SCIP_CALL( DECdecompSetLinkingconss(scip, decdecomp, detectordata->linkingconss, detectordata->nlinkingconss) );

   //hashmaps: shallow copy
   DECdecompSetVarindex(decdecomp, detectordata->indexmap->varindex);
   DECdecompSetConsindex(decdecomp, detectordata->indexmap->consindex);
   DECdecompSetVartoblock(decdecomp, detectordata->vartoblock);
   DECdecompSetConstoblock(decdecomp, detectordata->constoblock);
   //debug
//   PrintDetectordata(scip, detectordata);
//   DECdecompPrintDecomp(scip, decdecomp);
   return SCIP_OKAY;
}

/** resets detectordata such that it can be used for the next decomposition */
static
void resetDetectordata(DEC_DETECTORDATA* detectordata)
{
   SCIPhashmapRemoveAll(detectordata->vartoblock);
   SCIPhashmapRemoveAll(detectordata->constoblock);
   detectordata->nlinkingvars = 0;
   detectordata->nlinkingconss = 0;
}

static
SCIP_RETCODE blocking(
      SCIP* scip,
      DEC_DETECTORDATA* detectordata,
      DEC_DECOMP*** decdecomps,
      int* ndecdecomps,
      int nvars,
      int ncons,
      SCIP_RESULT* result
      )
{
   int n;   // maximum width of the band after ROC
   int v;   // minimum width of the band after ROC
   int tau; // desired number of blocks

   tau = 0;

   assert(*ndecdecomps == 0);
   //debug
   SCIPdebugMessage("Entering Blocking\n");
   //if multiple decompositions disabled
   if( detectordata->enablemultipledecomps == FALSE )
   {
      //if desiredblocks == 0 let the algorithm determine the desired number of blocks
      if( detectordata->desiredblocks == 0 )
      {
         n = maxArray(detectordata->width, ncons);
         v = minArray(detectordata->width, ncons);
         tau = round((nvars - v)/(n - v));
         SCIPdebugMessage("<n><v><tau>: <%i><%i><%i>\n", n, v, tau);
         if( tau > detectordata->maxblocks )
         {
            tau = detectordata->maxblocks;
         }
         //debug
         SCIPdebugMessage("detectordata->enablemultipledecomps == FALSE. detectordata->desiredblocks == 0. Calculating tau = %i\n", tau);
         //continue only if tau >= 2
         if( tau < 2 )
         {
            *result = SCIP_DIDNOTFIND;
            return SCIP_OKAY;
         }
      }
      else
      {
         tau = detectordata->desiredblocks;
      }
   }

   //variant 1
//   if( detectordata->enablemultipledecomps )
//   {
//      for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
//      {
//         if( detectordata->enableblockingdynamic )
//         {
//            blockingDynamic(scip, detectordata, tau, nvars);
//         }
//         if( detectordata->enableblockingstatic )
//         {
//            blockingStatic(scip, detectordata, tau, nvars);
//         }
//         if( detectordata->enableblockingassoonaspossible )
//         {
//            blockingAssoonaspossible(scip, detectordata, tau, nvars);
//         }
//      }
//   }
//   else
//   {
//      if( detectordata->enableblockingdynamic )
//      {
//         blockingDynamic(scip, detectordata, tau, nvars);
//      }
//      if( detectordata->enableblockingstatic )
//      {
//         blockingStatic(scip, detectordata, tau, nvars);
//      }
//      if( detectordata->enableblockingassoonaspossible )
//      {
//         blockingAssoonaspossible(scip, detectordata, tau, nvars);
//      }
//   }


   //variant 2
   //dynamic blocking
   if( detectordata->enableblockingdynamic )
   {
      //debug
      SCIPdebugMessage("detectordata->enableblockingdynamic == TRUE. \n");
      SCIP_CALL( rowsWithConstriction(scip, detectordata) );
      if( detectordata->enablemultipledecomps )
      {
         //debug
         SCIPdebugMessage("detectordata->enablemultipledecomps == TRUE. \n");
         for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
         {
            //debug
            SCIPdebugMessage("tau = %i \n", tau);
            resetDetectordata(detectordata);
            blockingDynamic(scip, detectordata, tau, nvars);
            //debug
            SCIPdebugMessage("dynamic blocking: copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[%i]);\n", *ndecdecomps);
            copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
            //debug
            DECdecompPrintDecomp(scip, (*decdecomps)[*ndecdecomps]);
            *ndecdecomps += 1;
         }
      }
      else
      {
         //debug
         SCIPdebugMessage("detectordata->enablemultipledecomps == FALSE. \n");
         resetDetectordata(detectordata);
         //debug
         SCIPdebugMessage("tau = %i \n", tau);
         blockingDynamic(scip, detectordata, tau, nvars);
         //debug
         SCIPdebugMessage("dynamic blocking: copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[%i]);\n", *ndecdecomps);
         copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
         //debug
         DECdecompPrintDecomp(scip, (*decdecomps)[*ndecdecomps]);
         *ndecdecomps += 1;
      }
   }

   //static blocking
   if( detectordata->enableblockingstatic )
   {
      //debug
      SCIPdebugMessage("detectordata->enableblockingstatic == TRUE. \n");
      if( detectordata->enablemultipledecomps )
      {
         for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
         {
            //debug
            SCIPdebugMessage("tau = %i \n", tau);
            resetDetectordata(detectordata);
            blockingStatic(scip, detectordata, tau, nvars);
            copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
            //debug
            DECdecompPrintDecomp(scip, (*decdecomps)[*ndecdecomps]);
            *ndecdecomps += 1;
         }
      }
      else
      {
         //debug
         SCIPdebugMessage("detectordata->enablemultipledecomps == FALSE. \n");
         resetDetectordata(detectordata);
         //debug
         SCIPdebugMessage("tau = %i \n", tau);
         blockingStatic(scip, detectordata, tau, nvars);
         //debug
         SCIPdebugMessage("static blocking: copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[%i]);\n", *ndecdecomps);
         copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
         //debug
         DECdecompPrintDecomp(scip, (*decdecomps)[*ndecdecomps]);
         *ndecdecomps += 1;
      }
   }

   //blocking ASAP
   if( detectordata->enableblockingassoonaspossible )
   {
      if( detectordata->enablemultipledecomps )
      {
         for( tau = detectordata->minblocks; tau <= detectordata->maxblocks; ++tau )
         {
            resetDetectordata(detectordata);
            blockingAsSoonAsPossible(scip, detectordata, tau, nvars);
            copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
            *ndecdecomps += 1;
         }
      }
      else
      {
         resetDetectordata(detectordata);
         blockingAsSoonAsPossible(scip, detectordata, tau, nvars);
         copyDetectorDataToDecomp(scip, detectordata, (*decdecomps)[*ndecdecomps]);
         *ndecdecomps += 1;
      }
   }
   return SCIP_OKAY;
}


static
DEC_DECL_INITDETECTOR(initStairheur)
{

   int i;
   int nvars;
   int nconss;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   SCIP_CALL( SCIPcreateWallClock(scip, &detectordata->clock) );
   SCIP_CALL( SCIPstartClock(scip, detectordata->clock) );
   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->consperblock, detectordata->maxblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varsperblock, detectordata->maxblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nconsperblock, detectordata->maxblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->nvarsperblock, detectordata->maxblocks) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->linkingvars, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->linkingconss, nconss) );
   for( i = 0; i < detectordata->maxblocks; ++i )
   {
      detectordata->nvarsperblock[i] = 0;
      detectordata->nconsperblock[i] = 0;
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->consperblock[i], nconss) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->varsperblock[i], nvars) );
   }

   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->ibegin, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->iend, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jbegin, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jend, nvars) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jmin, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->jmax, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->minV, nconss-1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->width, nconss) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &detectordata->hashmapindices, maximum(nvars, nconss) + 1) );
   for( i = 0; i < maximum(nvars, nconss)+1; ++i )
   {
      detectordata->hashmapindices[i] = i;
   }
   detectordata->rowsWithConstrictions = SCIPlistCreate(scip);
   detectordata->blockedAfterrow = SCIPlistCreate(scip);

   detectordata->nlinkingvars = 0;
   detectordata->nlinkingconss = 0;
   /* create hash tables */
   indexmapCreate(scip, &detectordata->indexmap, nconss, nvars);
//   SCIP_CALL( SCIPhashmapCreate(&detectordata->indexvar, SCIPblkmem(scip), nvars) );
//   SCIP_CALL( SCIPhashmapCreate(&detectordata->varindex, SCIPblkmem(scip), nvars) );
//   SCIP_CALL( SCIPhashmapCreate(&detectordata->indexcons, SCIPblkmem(scip), nconss) );
//   SCIP_CALL( SCIPhashmapCreate(&detectordata->consindex, SCIPblkmem(scip), nconss) );
   return SCIP_OKAY;
}

/** presolving deinitialization method of presolver (called after presolving has been finished) */
static
DEC_DECL_EXITDETECTOR(exitStairheur)
{
   int i;
   DEC_DETECTORDATA* detectordata;

   assert(scip != NULL);

   detectordata = DECdetectorGetData(detector);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(detector), DEC_DETECTORNAME) == 0);

   if( !detectordata->found )
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   /* free presolver data */
   for( i = 0; i < detectordata->maxblocks; ++i )
   {
      SCIPfreeMemoryArray(scip, &detectordata->consperblock[i]);
      SCIPfreeMemoryArray(scip, &detectordata->varsperblock[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->varsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nvarsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->consperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nconsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->linkingvars);
   SCIPfreeMemoryArray(scip, &detectordata->linkingconss);
   SCIPfreeMemoryArray(scip, &detectordata->relevantConss);

   SCIPfreeMemoryArray(scip, &detectordata->ibegin);
   SCIPfreeMemoryArray(scip, &detectordata->iend);
   SCIPfreeMemoryArray(scip, &detectordata->jbegin);
   SCIPfreeMemoryArray(scip, &detectordata->jend);
   SCIPfreeMemoryArray(scip, &detectordata->jmin);
   SCIPfreeMemoryArray(scip, &detectordata->jmax);
   SCIPfreeMemoryArray(scip, &detectordata->minV);
   SCIPfreeMemoryArray(scip, &detectordata->width);
   SCIPfreeMemoryArray(scip, &detectordata->hashmapindices);
   SCIP_CALL( SCIPfreeClock(scip, &detectordata->clock) );
   //delete lists
   //data had to be deleted before because of SCIP memory management
   SCIPlistDelete(scip, detectordata->rowsWithConstrictions);
   SCIPlistDelete(scip, detectordata->blockedAfterrow);
   //free deep copied hash maps
   //DO NOT FREE varindex and consindex because they are only shallow copied and contain the final permutation
   //debug
   SCIPhashmapFree(&detectordata->indexmap->indexvar);
   SCIPhashmapFree(&detectordata->indexmap->indexcons);
   SCIPfreeMemory(scip, &detectordata->indexmap);
   SCIPfreeMemory(scip, &detectordata);
   return SCIP_OKAY;
}

static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildStair)
{
   int i;
   int ncons; //number of constraints in the problem
   int nvars; //number of variables in the problem
   int ndecs;
   SCIP_VAR** vars_array;
   SCIP_CONS** cons_array;
   LIST* rowindices;
   LIST* columnindices;
#ifndef NDEBUG
   int ROC_iterations;
#endif

   assert(scip != NULL);
   assert(detectordata != NULL);
   assert(decdecomps != NULL);
   assert(ndecdecomps != NULL);

   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   SCIPwriteParams(scip, NULL, TRUE, TRUE);
   checkParameterConsistency(detectordata, result);
   ndecs = calculateNdecompositions(detectordata);
   SCIPdebugMessage("%i decompositions will be created\n", ndecs);
   *ndecdecomps = 0;

   /* allocate space for output data */
   SCIP_CALL( SCIPallocMemoryArray(scip, decdecomps, ndecs) );
   for( i = 0; i < ndecs; ++i )
   {
      SCIP_CALL_ABORT( DECdecompCreate(scip, &(*decdecomps)[i]) );
   }
   //remove empty constraints
   SCIP_CALL( findRelevantConss(scip, detectordata) );
   //SCIP_CALL( DECdecompCreate(scip, &detectordata->decdecomp) );
   nvars = SCIPgetNVars(scip);
   vars_array = SCIPgetVars(scip);
   ncons = detectordata->nRelevantConss;
   cons_array = detectordata->relevantConss;
   //initialize hash maps for keeping track of variables and constraints and their corresponding indices after being permuted by the ROC2-algorithm
   indexmapInit(detectordata->indexmap, vars_array, nvars, cons_array, ncons, detectordata->hashmapindices);
#ifndef NDEBUG
   {
      char filename[256];
      sprintf(filename, "%s_initial_problem", getProbNameWithoutPath(scip));
      plotInitialProblem(scip, detectordata, filename);
   }
#endif
   //initialize index arrays ibegin, iend, jbegin, jend
   rowindices = SCIPlistCreate(scip);
   columnindices = SCIPlistCreate(scip);
   rowindices_list(scip, detectordata, detectordata->indexmap->indexcons, detectordata->indexmap->varindex, &rowindices);
   //debug
//   printNested(rowindices, "rowindices in detectAndBuildStair");
   columnindices_list(scip, detectordata, rowindices, &columnindices);
   formIndexArray(detectordata->ibegin, detectordata->iend, rowindices);
   formIndexArray(detectordata->jbegin, detectordata->jend, columnindices);

   //debug
//   printArray(detectordata->ibegin, ncons, "ibegin");
//   printArray(detectordata->iend, ncons, "iend");

   //====================
   //===ROC2 algorithm===
   //====================
   SCIPdebugMessage("starting ROC2 algorithm\n");


#ifndef NDEBUG
   ROC_iterations = rankOrderClustering(scip, detectordata, detectordata->maxiterationsROC);
   {
      char filename[256];
      sprintf(filename, "%s_ROC", getProbNameWithoutPath(scip));
      plotInitialProblem(scip, detectordata, filename);
   }
   //check conditions for arrays ibegin and jbegin: ibegin[i]<=ibegin[i+k] for all positive k
   if( ROC_iterations < detectordata->maxiterationsROC || detectordata->maxiterationsROC  == -1 )
   {
      checkConsistencyOfIndexarrays(detectordata, nvars);
   }
#else
   (void) rankOrderClustering(scip, detectordata, detectordata->maxiterationsROC);
#endif
   //arrays jmin, jmax and minV
   SCIPdebugMessage("calculating index arrays\n");
   detectordata->jmin[0] = detectordata->ibegin[0];
   detectordata->jmax[0] = detectordata->iend[0];
   detectordata->width[0] = detectordata->iend[0] - detectordata->ibegin[0];
   for( i = 1; i < ncons; ++i )
   {
      detectordata->width[i] = detectordata->iend[i] - detectordata->ibegin[i];
      detectordata->jmin[i] = detectordata->ibegin[i];
      detectordata->jmax[i] = maximum(detectordata->iend[i], detectordata->jmax[i-1]);
      detectordata->minV[i-1]=1 + (detectordata->jmax[i-1] - detectordata->jmin[i]);
   }
   //====================
   //=====BLOCKING=======
   //====================
   //create the hashmaps constoblock and vartoblock
   SCIP_CALL( SCIPhashmapCreate(&detectordata->vartoblock, SCIPblkmem(scip), nvars) );
   SCIP_CALL( SCIPhashmapCreate(&detectordata->constoblock, SCIPblkmem(scip), detectordata->nRelevantConss) );

   blocking(scip, detectordata, decdecomps, ndecdecomps, nvars, ncons, result);
   //debug
   SCIPdebugMessage("Detected %i decompositions. Block sizes are ", *ndecdecomps);
   for( i = 0; i < *ndecdecomps; ++i )
   {
      SCIPinfoMessage(scip, NULL, "%i ", DECdecompGetNBlocks( (*decdecomps)[i] ));
   }
   SCIPinfoMessage(scip, NULL, "\n");

   //deallocate memory
   SCIPlistDeleteData(scip, detectordata->rowsWithConstrictions);
   SCIPlistDeleteNested(scip, rowindices);
   SCIPlistDeleteNested(scip, columnindices);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}


/** creates the stairheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionStairheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL( SCIPallocMemory(scip, &detectordata) );

   assert(detectordata != NULL);
   detectordata->found = FALSE;
//   detectordata->decdecomp = NULL;
   SCIP_CALL( DECincludeDetector(scip, DEC_DETECTORNAME, DEC_DECCHAR, DEC_DESC, DEC_PRIORITY, DEC_ENABLED, detectordata, detectAndBuildStair, initStairheur, exitStairheur) );

   /* add stairheur presolver parameters */
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/desiredblocks", "The desired number of blocks. 0 means automatic determination of the number of blocks.", &detectordata->desiredblocks, FALSE, DEFAULT_DESIREDBLOCKS, 0, 1000000, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/enableblockingdynamic", "Enable blocking type 'dynamic'", &detectordata->enableblockingdynamic, FALSE, DEFAULT_ENABLEBLOCKINGDYNAMIC, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/enableblockingstatic", "Enable blocking type 'static'", &detectordata->enableblockingstatic, FALSE, DEFAULT_ENABLEBLOCKINGSTATIC, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/enableblockingassoonaspossible", "Enable blocking type 'as soon as possible", &detectordata->enableblockingassoonaspossible, FALSE, DEFAULT_ENABLEBLOCKINGASSOONASPOSSIBLE, NULL, NULL) );
   SCIP_CALL( SCIPaddBoolParam(scip, "detectors/stairheur/enablemultipledecomps", "Enables multiple decompositions for all enabled blocking types. Ranging from minblocks to maxblocks", &detectordata->enablemultipledecomps, FALSE, DEFAULT_ENABLEMULTIPLEDECOMPS, NULL, NULL) );
   SCIP_CALL( SCIPaddIntParam(scip, "detectors/stairheur/maxiterationsROC", "The maximum number of iterations of the ROC-algorithm. -1 for no limit", &detectordata->maxiterationsROC, FALSE, DEFAULT_MAXITERATIONSROC, -1, 1000000, NULL, NULL) );
   return SCIP_OKAY;
}
