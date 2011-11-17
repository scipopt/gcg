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
#pragma ident "@(#) $Id: dec_stairheur.c,v 1.24 2010/01/04 20:35:45 bzfheinz Exp $"

/**@file   dec_stairheur.c
 * @ingroup DETECTORS
 * @brief  stairheur presolver
 * @author Mathias Luers
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/
#define SCIP_DEBUG

#include <assert.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>

#include "dec_stairheur.h"

#include "cons_decomp.h"
#include "struct_decomp.h"
#include "scip_misc.h"
#include "scip/struct_var.h"

#define DEC_DETECTORNAME      "stairheur"   /**< name of the detector */
#define DEC_PRIORITY          -100              /**< priority of the detector */

/* Default parameter settings*/
#define DEFAULT_BLOCKS                    2     /**< number of blocks */
#define DEFAULT_CONSWEIGHT                5     /**< weight for constraint hyperedges */
#define DEFAULT_RANDSEED                  1     /**< random seed for the hmetis call */
#define DEFAULT_TIDY                      TRUE  /**< whether to clean up afterwards */
#define DEFAULT_DUMMYNODES	              0.2   /**< percentage of dummy vertices*/

#define DEFAULT_MAXBLOCKS                 10    /**< value for the maximum number of blocks to be considered */
#define DEFAULT_MINBLOCKS                 2     /**< value for the minimum number of blocks to be considered */

#define DEFAULT_METIS_UBFACTOR            5.0   /**< default unbalance factor given to metis on the commandline */
#define DEFAULT_METIS_VERBOSE             FALSE /**< should metis be verbose */
#define DEFAULT_METISUSEPTYPE_RB          TRUE  /**< Should metis use the rb or kway partitioning algorithm */
#define DEFAULT_PRIORITY                  DEC_PRIORITY

#define DWSOLVER_REFNAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_ref.txt", (name), (blocks), (cons), (dummy)

#define GP_NAME(name, blocks, cons, dummy) "%s_%d_%d_%.1f_%d.gp", (name), (blocks), (cons), (dummy)


/*include the code for doubly linked lists */
#include <stdio.h>
#include <stdlib.h>

/** TODO:
 * deallocate all memory (lists, iterators, data stored in lists...)
 *
 * */

typedef enum { false, true } bool;

struct node {
 void* data;
 struct node* prev;
 struct node* next;
};
typedef struct node NODE;

struct list {
   /** a doubly linked list
    * the first and last node are referred to by nil->next and nil->prev, respectively.
    */
   NODE* nil;
   int size;
};
typedef struct list LIST;

struct iterator {
   NODE* node;
   LIST* list;
};
typedef struct iterator ITERATOR;

LIST* list_create_empty(void);
LIST* list_create(void* data);
LIST* list_int_create(int from, int to);
LIST* list_flat_copy(LIST* list);
bool list_push_front(LIST* list, void* data);
bool list_pop_front(LIST* list);
bool list_push_back(LIST* list, void* data);
bool list_pop_back(LIST* list);
NODE* list_insert(ITERATOR* it, void* data);
bool list_assign(ITERATOR* it, void* data);
bool list_assign_list(LIST* list1, LIST* list2);
bool list_empty(LIST* list);
bool list_erase(ITERATOR* it);
bool list_delete(LIST* list);
bool list_delete_data(LIST* list);
bool list_delete_nested(LIST* list);
ITERATOR list_find(ITERATOR first, ITERATOR last, void* value, bool (*comp_func)(void* a, void* b));
bool list_move(ITERATOR position, ITERATOR it);
void list_move_front(ITERATOR position);
void list_move_back(ITERATOR position);
SCIP_RETCODE list_rearrange(LIST* list, LIST* order);
bool list_foreach(LIST* list, int(*func)(void*));

NODE* node_create(void* data);
int printstring(void *s);
int printint(void *i);
void list_prin(LIST* list, int(*func)(void*));
static int compare (const void * a, const void * b);
static bool compare_int(void* a, void * b);

ITERATOR* iterator_create(void);
ITERATOR iterator_begin(LIST* list);
ITERATOR iterator_end(LIST* list);
bool iterator_next(ITERATOR* it);
bool iterator_prev(ITERATOR* it);
void iterator_equals(ITERATOR it1, ITERATOR it2);
bool iterator_is_equal(ITERATOR it1, ITERATOR it2);

LIST* list_create_empty()
{
   LIST* list;
   NODE* node;
   if(! (list = malloc(sizeof(LIST)))) return NULL;
   //debugging
//   int* i;
//   i = malloc(sizeof(int));
//   *i = 0;
   //node = node_create(i);
   node = node_create(NULL);
   list->nil = node;
   list->nil->next = list->nil;
   list->nil->prev = list->nil;
   list->nil->data = node->data;
   list->size = 0;
   return list;
}

LIST* list_create(void* data)
{
   LIST* list;
   list = list_create_empty();
   list_push_front(list, data);
   return list;
}

LIST* list_int_create(int from, int to)
{
   LIST* list;
   int* data;
   int i;
   list = list_create_empty();
   for(i = from; i <= to; ++i)
   {
      data = (int*) malloc(sizeof(int));
      *data=i;
      list_push_back(list, (void*) data);
   }
   return list;
}

LIST* list_flat_copy(LIST* list)
{
   /**creates a copy of list, but does not copy the data stored in list->data */
   LIST* copy;
   ITERATOR it;
   //case: list exists
   if(list)
   {
      copy = list_create_empty();
      for(it = iterator_begin(list); ! iterator_is_equal(it, iterator_end(list)); iterator_next(&it))
      {
         list_push_back(copy, it.node->data);
      }
      return copy;
   }
   //case: list does not exist
   else
   {
      return NULL;
   }
}

bool list_push_front(LIST* list, void* data)
{
   /** adds a node at the front of the list */
   NODE* newnode;
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return false;
   newnode = node_create(data);
   //it points to the beginning of list
   it = iterator_begin(list);
   list_insert(&it, data);
   return true;
}

bool list_pop_front(LIST* list)
{
   /** removes the node at the front of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (list_empty(list)) ) return false;

   it = iterator_begin(list);
   list_erase(&it);
   return true;
}

bool list_push_back(LIST* list, void* data)
{
   /** adds a node at the back of the list */
   NODE* newnode;
   ITERATOR it = {NULL, NULL};
   //case: list does not exist
   if(! list) return false;
   newnode = node_create(data);
   //it points to the beginning of list
   it = iterator_end(list);
   list_insert(&it, data);
   return true;
}

bool list_pop_back(LIST* list)
{
   /** removes the node at the back of the list */
   ITERATOR it = {NULL, NULL};
   //case: list does not exist or is empty
   if( (! list) || (list_empty(list)) ) return false;

   it = iterator_end(list);
   iterator_prev(&it);
   list_erase(&it);
   return true;
}

NODE* list_insert(ITERATOR* it, void* data)
{
   /** creates a new node before 'it->node' and returns the new node*/
   NODE* newnode;
   NODE* successor;
   NODE* predecessor;
   //case: 'it' points to no node or list
   if(it->node == NULL && it->list == NULL) return NULL;
   newnode=node_create(data);
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

bool list_assign(ITERATOR* it, void* data)
{
   /** assign new data to 'it->node' */
   if(it && it->node)
   {
      it->node->data=data;
      return true;
   }
   else
   {
      return false;
   }
}

bool list_assign_list(LIST* target, LIST* origin)
{
   /**if target and origin have the same size, this function assigns the data pointers of origin to target */
   ITERATOR it1;
   ITERATOR it2;
   if(target && origin && target->size == origin->size)
   {
      for(it1 = iterator_begin(target), it2 = iterator_begin(origin); ! iterator_is_equal(it1, iterator_end(target)); iterator_next(&it1), iterator_next(&it2))
      {
         list_assign(&it1, it2.node->data);
      }
      return true;
   }
   else
   {
      return false;
   }
}

NODE* node_create(void* data)
{
   /** creates a node with a void pointer to data and returns the node */
   NODE* node;
   if(! (node = malloc(sizeof(NODE)))) return NULL;
   node->data = data;
   node->prev = NULL;
   node->next = NULL;
   return node;
}

bool list_erase(ITERATOR* it)
{
   /** removes the node 'it->node' from the list 'it->list'
    * returns true if 'node' is in 'list'
    * returns false if 'node' is not in 'list'
    * the memory the data pointer points to is not deallocated*/
   NODE* successor;
   NODE* predecessor;
   //case: iterator, node or list do not exist or empty list
   if( (! it) || (! it->node) || (! it->list) || list_empty(it->list))
   {
      return false;
   }
   successor = it->node->next;
   predecessor = it->node->prev;
   predecessor->next = successor;
   successor->prev = predecessor;
   free(it->node);
   it->node = successor;
   --it->list->size;
   return true;
}

bool list_delete(LIST* list)
{
   /** deletes the entire list
    *  for deallocating memory of the list data, call 'list_deleta_data' before */
   if(! list)
   {
      return false;
   }
   else
   {
      while(! list_empty(list))
      {
         list_pop_back(list);
      }
      //deallocate sentinel node
      free(list->nil);
      free(list);
      list = NULL;
      return true;
   }
}

bool list_delete_data(LIST* list)
{
   /** deallocates the memory the data pointers of the list points to */
   ITERATOR it;
   if(! list)
   {
      return false;
   }
   else
   {
      for( it = iterator_begin(list); ! ( iterator_is_equal(it, iterator_end(list)) ); iterator_next(&it))
      {
         free(it.node->data);
      }
      return true;
   }
}
bool list_delete_nested(LIST* list)
{
   /** deallocates all memory for a nested list including the data */
   ITERATOR it1;
   if(! list)
   {
      return false;
   }
   else
   {
      for(it1 = iterator_begin(list);  ! ( iterator_is_equal(it1, iterator_end(list)) ); iterator_next(&it1))
         {
            list_delete_data(it1.node->data);
            list_delete(it1.node->data);
         }
      list_delete(list);
      return true;
   }
}

bool list_empty(LIST* list)
{
   /** returns true iff list exists and has 0 nodes (only sentinel node)*/
//   if( list && list->size == 0)
   //list is empty if the next pointer of the sentinel points to itself
   if( list && list->nil->next == list->nil)
   {
      return true;
   }
   else
   {
      return false;
   }
}

bool list_foreach(LIST* list, int(*func)(void*))
{
   /** runs the function 'func' for every node of the list 'node' */
   ITERATOR it;
   for(it = iterator_begin(list); ! ( iterator_is_equal(it, iterator_end(list)) ); iterator_next(&it))
   {
      if(func(it.node->data) != 0)
      {
         return false;
      }
   }
   return true;
}

void list_prin(LIST* list, int(*func)(void*))
{
   printf("( ");
   list_foreach(list, func);
   printf(")\n");
}

ITERATOR list_find(ITERATOR first, ITERATOR last, void* value, bool (*comp_func)(void* a, void* b))
{
   for (; ! (iterator_is_equal(first, last)); iterator_next(&first))
   {
      if ( comp_func(value, first.node->data) )
      {
         break;
      }
   }
   return first;
}

bool list_move(ITERATOR position, ITERATOR it)
{
   /** moves the element at 'position' to the location in front of element 'it'
    * both elements must be in the same list
    */
   NODE* old_successor;
   NODE* old_predecessor;
   NODE* new_successor;
   NODE* new_predecessor;
   if(position.node && it.node && position.list == it.list)
   {
      //case 'it' and 'position' point to the same node or 'it' points to the succeeding node of 'position': no moving
      if(position.node == it.node || position.node->next == it.node)
      {
         return true;
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

      return true;
   }
   else
   {
      return false;
   }
}

void list_move_front(ITERATOR position)
{
   /** moves the element at 'position' to the front of the list */
   list_move(position, iterator_begin(position.list));
}

void list_move_back(ITERATOR position)
{
   /** moves the element at 'position' to the end of the list */
   list_move(position, iterator_end(position.list));
}

SCIP_RETCODE list_rearrange(LIST* list, LIST* order)
{
   /** rearranges elements of list according to the ordering of order
    * example: list = (a b c d); order = (3 2 4 1)
    * after calling list_rearrange(list, order): list = (c b d a)
    * both lists must have same size!
    * order must have elements from 1 to list->size
    */
   LIST* new_list;
   ITERATOR it1, it2;
   int i;
   if( list && order && list->size == order->size)
   {
      new_list = list_create_empty();
      for(it1 = iterator_begin(order); ! iterator_is_equal(it1, iterator_end(order)); iterator_next(&it1))
      {
         for( it2 = iterator_begin(list), i = 1; i < *(int*)it1.node->data; ++i)
         {
            iterator_next(&it2);
         }
         list_push_back(new_list, it2.node->data);
      }
      list_assign_list(list, new_list);
      list_delete(new_list);
      return SCIP_OKAY;
   }
   else
   {
      return SCIP_ERROR;
   }
}

int printstring(void *s)
{
   printf("%s\n", (char *)s);
   return 0;
}

int printint(void *i)
{
//   int* n;
//   n=(int*) i;
//   printf("%i\n", *n);
//   return 0;
   printf("%i ", *(int*)i);
   return 0;
}

//forward iterators
ITERATOR* iterator_create()
{
   ITERATOR* it;
   if(! (it = malloc(sizeof(ITERATOR)))) return NULL;
   it->list = NULL;
   it->node = NULL;
   return it;
}

ITERATOR iterator_begin(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil->next;
   }
   return it;
}

ITERATOR iterator_end(LIST* list)
{
   ITERATOR it = {NULL, NULL};
   if(list)
   {
      it.list = list;
      it.node = list->nil;
   }
   return it;
}


bool iterator_next(ITERATOR* it)
{
   /** sets it to the next element of the list
    */
   if(it && it->list && it->node)
   {
      it->node = it->node->next;
      return true;
   }
   else
   {
      return false;
   }
}

bool iterator_prev(ITERATOR* it)
{
   /** sets it to the previous element of the list
    */
   if(it && it->list && it->node)
   {
      it->node = it->node->prev;
      return true;
   }
   else
   {
      return false;
   }
}

void iterator_equals(ITERATOR it1, ITERATOR it2)
{
   /** assigns it1 = it2 */
   it1.list = it2.list;
   it1.node = it2.node;
}

bool iterator_is_equal(ITERATOR it1, ITERATOR it2)
{
   if(it1.list == it2.list && it1.node == it2.node)
   {
      return true;
   }
   else
   {
      return false;
   }
}

/*
 * Data structures
 */
///** score data structure **/
//struct Dec_StairheurScores
//{
//   SCIP_Real borderscore;
//   SCIP_Real minkequicutscore;
//   SCIP_Real equicutscorenormalized;
//   SCIP_Real densityscore;
//   SCIP_Real linkingscore;
//};
typedef struct Dec_StairheurScores DEC_STAIRHEURSCORES;

/** detector data */
struct DEC_DetectorData
{
   DECDECOMP* decdecomp;
   SCIP_VAR*** varsperblock;
   int* nvarsperblock;
   SCIP_CONS*** consperblock;
   int *nconsperblock;
   SCIP_CONS** linkingconss;
   int nlinkingconss;
   int blocks;
   int maxblocks;
   int minblocks;
   SCIP_CONS** relevantConss;
   int nRelevantConss;
   SCIP_HASHMAP* indexvar_old;  // index in problem -> variable
   SCIP_HASHMAP* indexvar_new;  // index in permuted problem -> variable
   SCIP_HASHMAP* varindex_old;  // variable -> index in problem
   SCIP_HASHMAP* varindex_new;  // variable -> index in permuted problem
   SCIP_HASHMAP* indexcons_old; // index in problem -> constraint
   SCIP_HASHMAP* indexcons_new; // index in permuted problem -> constraint
   SCIP_HASHMAP* consindex_old; // constraint -> index in problem
   SCIP_HASHMAP* consindex_new; // constraint -> index in permuted problem
   int* ibegin_old; //array, ibegin[i]: index of first nonzero entry in row i before permutation
   int* iend_old;   //array, iend[i]: index of last nonzero entry in row i before permutation
   int* jbegin_old; //array, jbegin[j]: index of first nonzero entry in column j before permutation
   int* jend_old;   //array, jend[j]: index of last nonzero entry in column j before permutation
   int* ibegin_new; //array, ibegin[i]: index of first nonzero entry in row i after permutation
   int* iend_new;   //array, iend[i]: index of last nonzero entry in row i after permutation
   int* jbegin_new; //array, jbegin[j]: index of first nonzero entry in column j after permutation
   int* jend_new;   //array, jend[j]: index of last nonzero entry in column j after permutation
   int* jmin;   //array, jmin[i]: index of first nonzero column of the i-th row
   int* jmax;   //array, jmax[i]: the last nonzero entry among all rows prior to and including the i-th row
   int* minV;   //array, minV[i]: number of linking variables corresponding to a partitioning after the i-th row
   int* width;  //array, width[i]: width of the band (of nonzero entries after ROC) at row i
   int n;   //maximum width of the band (of nonzero entries after ROC)
   int v;   //minimum width of the band
   int tau; //approximate number of blocks
   LIST* rowsWithConstrictions;
   SCIP_Bool found;
   int priority;
};

/*
 * Local methods
 */

/* put your local methods here, and declare them static */

static
int compare (const void * a, const void * b)
{
  return ( *(int*)a - *(int*)b );
}

static
bool compare_int(void* a, void * b)
{
   return ( *(int*)a == *(int*)b ? true : false );
}

//int round(double number)
//{
//    return (number >= 0) ? (int)(number + 0.5) : (int)(number - 0.5);
//}
//
static
int maximum(int a, int b)
{
   return (a > b ? a : b);
}

static
int max_array(int* a, int num_elements)
{
   int i;
   int max;
   if(num_elements > 0 && a != NULL)
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

static
int min_array(int* a, int num_elements)
{
   int i;
   int min;
   if(num_elements > 0 && a != NULL)
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

static
void ChangeHashMapPointers(SCIP_HASHMAP** hm1, SCIP_HASHMAP** hm2)
{
   SCIP_HASHMAP* hm3; /* local for swap */
    hm3 = *hm2;
    *hm2 = *hm1;
    *hm1= hm3;
}

static
void ChangeIntPointers(int** i1, int** i2)
{
   int* i3; /* local for swap */
    i3 = *i2;
    *i2 = *i1;
    *i1= i3;
}

static
void findRelevantConss(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   /** scans all constraints of the constraint array of the scip object,
    * and stores pointers to all constraints that have at least one variable in detectordata->relevantConss.
    * Thus it removes all empty constraints.
    */
   SCIP_CONS** cons_array;
   LIST* relevantConssIndices;
   int i;
   int* data;
   ITERATOR it1;
   cons_array = SCIPgetConss(scip);
   relevantConssIndices = list_create_empty();
   for(i = 0; i < SCIPgetNConss(scip); ++i)
   {
      if(SCIPgetNVarsXXX(scip, cons_array[i]) > 0)
      {
         data = malloc(sizeof(int));
         *data = i;
         list_push_back(relevantConssIndices, data);
      }
   }
   list_prin(relevantConssIndices, printint);
   //allocate memory for detectordata->relevantConss and store pointers of relevant conss
   detectordata->nRelevantConss = relevantConssIndices->size;
   printf("nRelevantConss: %i \n", detectordata->nRelevantConss);
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->relevantConss, detectordata->nRelevantConss));
   for(i = 0, it1 = iterator_begin(relevantConssIndices); i < detectordata->nRelevantConss; ++i, iterator_next(&it1))
   {
      printf("i, *(int*)it1.node->data: %i, %i \n", i, *(int*)it1.node->data);
      detectordata->relevantConss[i] = cons_array[*(int*)it1.node->data];
   }
   list_delete_data(relevantConssIndices);
   list_delete(relevantConssIndices);
}

static
LIST* rowindices_list(SCIP* scip, DEC_DETECTORDATA* detectordata, SCIP_HASHMAP* indexcons, SCIP_HASHMAP* varindex)
{
   //create the rowindices list
   int i;
   int j;
   int* data;
   LIST* rowindices;
   LIST* rowindices_row;
   int ncons; //number of constraints of the problem
   int nconsvars; //number of variables in a constraint
   int* probindices;
   SCIP_CONS* cons; //one constraint of the problem
   SCIP_VAR** consvars; //array of variables that occur in a constraint (unequal zero)

   rowindices = list_create_empty();
   ncons = detectordata->nRelevantConss;
   for(i = 0; i < ncons; ++i)
   {
      cons = (SCIP_CONS*) SCIPhashmapGetImage(indexcons, (void*) (i+1));
      nconsvars = SCIPgetNVarsXXX(scip, cons);
      consvars = SCIPgetVarsXXX(scip, cons);
      //allocate memory for the array of probindices
      probindices=(int*)malloc(nconsvars*sizeof(int));
      //fill the array with the indices of the variables of the current constraint
      for(j = 0; j < nconsvars; ++j)
      {
//         probindices[j] = SCIPvarGetProbindex(consvars[j])+1;
         probindices[j] = (int*) SCIPhashmapGetImage(varindex, consvars[j]);
      }
      //sort the elements of probindices ('<')
      qsort(probindices, nconsvars, sizeof(int), compare);
      //store a copy of the elements of probindices in the list rowindices_row
      rowindices_row = list_create_empty();
      for(j = 0; j < nconsvars; ++j)
      {
         data = (int*) malloc(sizeof(int));
         *data = probindices[j];
         list_push_back(rowindices_row, data);
      }
      //deallocate memory
      free(probindices);
      //add rowindices_row to the list rowindices
      list_push_back(rowindices, rowindices_row);
   }
   return rowindices;
}

static
LIST* columnindices_list(SCIP* scip, DEC_DETECTORDATA* detectordata, LIST* rowindices)
{
   /** creates a nested list which contains the position of non-zero entries of the boundary constraints */
   LIST** columnindices_array;
   LIST* columnindices;
   LIST* rowindices_row;
   int* data;
   int position;
   int nvars;
   int ncons;
   int i;
   ITERATOR it1;
   ITERATOR it2;
   nvars = SCIPgetNVars(scip);
   ncons = detectordata->nRelevantConss;
   //create the columnindices_array with pointers to empty lists
   columnindices_array = malloc(nvars * sizeof(LIST*));
   for(i = 0; i < nvars; ++i)
   {
      columnindices_array[i] = list_create_empty();
   }

   for(it1 = iterator_begin(rowindices), i = 0; ! ( iterator_is_equal(it1, iterator_end(rowindices)) ); iterator_next(&it1), ++i)
   {
      rowindices_row = it1.node->data;
      for(it2 = iterator_begin(rowindices_row); ! ( iterator_is_equal(it2, iterator_end(rowindices_row)) ); iterator_next(&it2))
      {
         data = (int*) malloc(sizeof(int));
         *data = i+1;
         position = *(int*)(it2.node->data)-1;
         list_push_back(columnindices_array[position], data);
      }
   }
   //create a columnindices list instead of an array
   columnindices = list_create_empty();
   for(i = 0; i < nvars; ++i)
   {
      list_push_back(columnindices, columnindices_array[i]);
   }
   //deallocate memory
   free(columnindices_array);
   return columnindices;
}

static
LIST* row_ordering(LIST* roworder, LIST* columnindices)
{
   /** does the row ordering of the ROC2 algortihm */
   LIST* new_roworder;
   ITERATOR it1;
   ITERATOR it2;
   ITERATOR it3;
   ITERATOR it4;
   new_roworder = list_flat_copy(roworder);

   for( it1 = iterator_end(columnindices), iterator_prev(&it1); it1.node != it1.list->nil; iterator_prev(&it1))
   {
      for( it2 = iterator_end(roworder), iterator_prev(&it2); it2.node != it2.list->nil; iterator_prev(&it2) )
      {
         it3 = list_find(iterator_begin(it1.node->data), iterator_end(it1.node->data), it2.node->data, compare_int);
         if(it3.node->data != NULL)
         {
            it4 = list_find(iterator_begin(new_roworder), iterator_end(new_roworder), it2.node->data, compare_int);
            list_move_front(it4);
         }
         else
         {
            ;
         }
      }
      list_assign_list(roworder, new_roworder);
   }
   //deallocate memory
   list_delete(new_roworder);
   return roworder;
}

static
SCIP_RETCODE formIndexArray(int* begin, int* end, LIST* indices)
{
   ITERATOR it1;
   ITERATOR it2;
   int i;
   assert(begin != NULL && end != NULL && indices != NULL);
   for(it1 = iterator_begin(indices), i = 0; ! iterator_is_equal(it1, iterator_end(indices)); iterator_next(&it1), ++i)
   {
      //case: list not empty
      if(! list_empty(it1.node->data))
      {
         it2 = iterator_begin(it1.node->data);
         begin[i] = *(int*)it2.node->data;
         it2 = iterator_end(it1.node->data);
         iterator_prev(&it2);
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

static
bool IndexArrayChanges(SCIP* scip, DEC_DETECTORDATA* detectordata)
{
   int i;
   //any change in ibegin or iend?
   for(i = 0; i < detectordata->nRelevantConss; ++i)
   {
      if(   detectordata->ibegin_old[i] != detectordata->ibegin_new[i]
         || detectordata->iend_old[i] != detectordata->iend_new[i] )
      {
         return true;
      }
   }
   //any change in jbegin or jend?
   for(i = 0; i < SCIPgetNVars(scip); ++i)
   {
      if(   detectordata->jbegin_old[i] != detectordata->jbegin_new[i]
         || detectordata->jend_old[i] != detectordata->jend_new[i] )
      {
         return true;
      }
   }
   //case: all entries of old and new are equal
   return false;
}

static
SCIP_RETCODE RankOrderClustering(
      SCIP*             scip,       /**< SCIP data structure */
      DEC_DETECTORDATA*  detectordata /**< presolver data data structure */
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
   SCIP_CONS* cons;
   SCIP_VAR* var;

   assert(scip != NULL);
   assert(detectordata != NULL);

   nvars = SCIPgetNVars(scip);
   ncons = detectordata->nRelevantConss;

   //create a list for the order of the rows ( 1 2 3 ... ncons ) and columns ( 1 2 3 ... nvars )
   roworder = list_int_create(1, ncons);
   columnorder = list_int_create(1, nvars);

   //create the lists containing the positions of nonzero entries
   rowindices = rowindices_list(scip, detectordata, detectordata->indexcons_old, detectordata->varindex_old);
   columnindices = columnindices_list(scip, detectordata, rowindices);

   //row ordering
   roworder = row_ordering(roworder, columnindices);
   //rearrange the rowindices list according to the new rowordering
   list_rearrange(rowindices, roworder);
   //column ordering
   columnorder = row_ordering(columnorder, rowindices);

   printf("row ordering:");
   list_prin(roworder, printint);
   printf("column ordering:");
   list_prin(columnorder, printint);

   //consindex and indexcons
   for(it1 = iterator_begin(roworder), i = 0; ! iterator_is_equal(it1, iterator_end(roworder)) && i < ncons; ++i, iterator_next(&it1))
   {
      position = (*(int*)it1.node->data);
      cons = SCIPhashmapGetImage(detectordata->indexcons_old, (void*) position);
      assert ( cons != NULL);

      //consindex
      assert( SCIPhashmapExists(detectordata->consindex_new, (void*) cons));
      SCIPhashmapSetImage(detectordata->consindex_new, (void*) cons, (void*) (i+1));
      //indexcons
      assert( SCIPhashmapExists(detectordata->indexcons_new, (void*) (i+1) ));
      SCIPhashmapSetImage(detectordata->indexcons_new, (void*) (i+1), cons);
   }
   //varindex and indexvar
   for(it1 = iterator_begin(columnorder), i = 0; ! iterator_is_equal(it1, iterator_end(columnorder)) &&i < nvars; ++i, iterator_next(&it1))
   {
      position = (*(int*)it1.node->data);
      var = (SCIP_VAR*) SCIPhashmapGetImage(detectordata->indexvar_old, (void*) position);
      assert ( var != NULL);

      //varindex
      assert( SCIPhashmapExists(detectordata->varindex_new, (void*) var) );
      SCIPhashmapSetImage(detectordata->varindex_new, (void*) var, (void*) (i+1));
      //indexvar
      assert( SCIPhashmapExists(detectordata->indexvar_new, (void*) (i+1) ));
      SCIPhashmapSetImage(detectordata->indexvar_new, (void*) (i+1), var);
   }

   //form the new index arrays after the permutation
   list_delete_nested(rowindices);
   list_delete_nested(columnindices);
   rowindices = rowindices_list(scip, detectordata, detectordata->indexcons_new, detectordata->varindex_new);
   columnindices = columnindices_list(scip, detectordata, rowindices);
   formIndexArray(detectordata->ibegin_new, detectordata->iend_new, rowindices);
   formIndexArray(detectordata->jbegin_new, detectordata->jend_new, columnindices);
//   //debug
//   for(i = 0; i < ncons; ++i)
//   {
//      printf("ibegin_old, ibegin_new: [%i] [%i]\n", detectordata->ibegin_old[i], detectordata->ibegin_new[i]);
//      printf("iend_old, iend_new: [%i] [%i]\n", detectordata->iend_old[i], detectordata->iend_new[i]);
//   }
   //switch between index arrays containing new and old indices
   ChangeIntPointers(&detectordata->ibegin_old, &detectordata->ibegin_new);
   ChangeIntPointers(&detectordata->iend_old, &detectordata->iend_new);
   ChangeIntPointers(&detectordata->jbegin_old, &detectordata->jbegin_new);
   ChangeIntPointers(&detectordata->jend_old, &detectordata->jend_new);
   //switch between hash maps containing new and old indices
   ChangeHashMapPointers(&detectordata->indexvar_old, &detectordata->indexvar_new);
   ChangeHashMapPointers(&detectordata->varindex_old, &detectordata->varindex_new);
   ChangeHashMapPointers(&detectordata->indexcons_old, &detectordata->indexcons_new);
   ChangeHashMapPointers(&detectordata->consindex_old, &detectordata->consindex_new);

   //deallocate memory
   list_delete_data(roworder);
   list_delete(roworder);
   list_delete_data(columnorder);
   list_delete(columnorder);
   list_delete_nested(rowindices);
   list_delete_nested(columnindices);
   return SCIP_OKAY;
}

static
void rowsWithConstriction(DEC_DETECTORDATA* detectordata)
{
   int i;
   int* data;
   for(i = 1; i < detectordata->nRelevantConss - 2; ++i)
   {
      //is minV[i] a local minimum?    < or <=   ? What does make more sense?
      if(detectordata->minV[i] < detectordata->minV[i-1] && detectordata->minV[i] < detectordata->minV[i+1])
      {
         data = (int*) malloc(sizeof(int));
         *data = i;
         list_push_back(detectordata->rowsWithConstrictions, data);
      }
   }
}

/** copies the variable and block information to the decomp structure */
static
SCIP_RETCODE copyDetectorDataToDecomp(
      SCIP*             scip,         /**< SCIP data structure */
      DEC_DETECTORDATA* detectordata, /**< presolver data data structure */
      DECDECOMP*        decomp        /**< DECOMP data structure */
      )
{
   int i;
   assert(scip != 0);
   assert(detectordata != 0);
   assert(decomp != 0);

   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipvars, detectordata->blocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &decomp->subscipconss, detectordata->blocks));

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->linkingconss, detectordata->linkingconss, detectordata->nlinkingconss));
   decomp->nlinkingconss = detectordata->nlinkingconss;

   for( i = 0; i < detectordata->blocks; ++i )
   {
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipconss[i], detectordata->consperblock[i], detectordata->nconsperblock[i]));
      SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->subscipvars[i], detectordata->varsperblock[i], detectordata->nvarsperblock[i]));
   }

   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipconss, detectordata->nconsperblock, detectordata->blocks));
   SCIP_CALL(SCIPduplicateMemoryArray(scip, &decomp->nsubscipvars, detectordata->nvarsperblock, detectordata->blocks));

   decomp->varindex = detectordata->varindex_old;
   decomp->consindex = detectordata->consindex_old;
   decomp->nblocks = detectordata->blocks;
   decomp->type = DEC_STAIRCASE;
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

   nvars = SCIPgetNVars(scip);
   nconss = SCIPgetNConss(scip);
   detectordata->maxblocks = MIN(nconss, detectordata->maxblocks);
   /* initialize variables and constraints per block structures*/
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nconsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->nvarsperblock, detectordata->maxblocks));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->linkingconss, nconss));
   for( i = 0; i < detectordata->maxblocks; ++i )
   {
      detectordata->nvarsperblock[i] = 0;
      detectordata->nconsperblock[i] = 0;
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->consperblock[i], nconss));
      SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->varsperblock[i], nvars));
   }

   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->ibegin_old, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->iend_old, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jbegin_old, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jend_old, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->ibegin_new, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->iend_new, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jbegin_new, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jend_new, nvars));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jmin, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->jmax, nconss));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->minV, nconss-1));
   SCIP_CALL(SCIPallocMemoryArray(scip, &detectordata->width, nconss));

   detectordata->rowsWithConstrictions = list_create_empty();

   detectordata->nlinkingconss = 0;
   /* create hash tables */
   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexvar_old, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexvar_new, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->varindex_old, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->varindex_new, SCIPblkmem(scip), nvars));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexcons_old, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->indexcons_new, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->consindex_old, SCIPblkmem(scip), nconss));
   SCIP_CALL(SCIPhashmapCreate(&detectordata->consindex_new, SCIPblkmem(scip), nconss));
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

   if( !detectordata->found)
   {
      SCIPfreeMemory(scip, &detectordata);
      return SCIP_OKAY;
   }

   /* free presolver data */
   for( i = 0; i < detectordata->maxblocks; ++i)
   {
      SCIPfreeMemoryArray(scip, &detectordata->consperblock[i]);
      SCIPfreeMemoryArray(scip, &detectordata->varsperblock[i]);
   }
   SCIPfreeMemoryArray(scip, &detectordata->varsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nvarsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->consperblock);
   SCIPfreeMemoryArray(scip, &detectordata->nconsperblock);
   SCIPfreeMemoryArray(scip, &detectordata->linkingconss);
   SCIPfreeMemoryArray(scip, &detectordata->relevantConss);
   SCIPfreeMemory(scip, &detectordata);

   SCIPfreeMemoryArray(scip, &detectordata->ibegin_new);
   SCIPfreeMemoryArray(scip, &detectordata->iend_new);
   SCIPfreeMemoryArray(scip, &detectordata->jbegin_new);
   SCIPfreeMemoryArray(scip, &detectordata->jend_new);
   SCIPfreeMemoryArray(scip, &detectordata->ibegin_old);
   SCIPfreeMemoryArray(scip, &detectordata->iend_old);
   SCIPfreeMemoryArray(scip, &detectordata->jbegin_old);
   SCIPfreeMemoryArray(scip, &detectordata->jend_old);
   SCIPfreeMemoryArray(scip, &detectordata->jmin);
   SCIPfreeMemoryArray(scip, &detectordata->jmax);
   SCIPfreeMemoryArray(scip, &detectordata->jmax);
   SCIPfreeMemoryArray(scip, &detectordata->minV);
   SCIPfreeMemoryArray(scip, &detectordata->width);
   //delete lists
   list_delete_data(detectordata->rowsWithConstrictions);
   list_delete(detectordata->rowsWithConstrictions);
   //free deep copied hash maps
   //DO NOT FREE varindex_old and consindex_old because they are only shallow copied and contain the final permuation
   SCIPhashmapFree(&detectordata->indexvar_old);
   SCIPhashmapFree(&detectordata->indexvar_new);
   SCIPhashmapFree(&detectordata->varindex_new);
   SCIPhashmapFree(&detectordata->indexcons_old);
   SCIPhashmapFree(&detectordata->indexcons_new);
   SCIPhashmapFree(&detectordata->consindex_new);
   return SCIP_OKAY;
}

static
DEC_DECL_DETECTSTRUCTURE(detectAndBuildStair)
{
   int i;
   int probindex;
   int ncons; //number of constraints in the problem
   int nvars; //number of variables in the problem
   SCIP_VAR** vars_array;
   SCIP_VAR* var;
   SCIP_CONS** cons_array;
   SCIP_CONS* cons;
   DEC_DETECTOR* stairheur;
   DEC_DETECTORDATA* detectordata;
   LIST* rowindices;
   LIST* columnindices;

   assert(scip != NULL);
   stairheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(stairheur);
   assert(detectordata != NULL);
   assert(strcmp(DECdetectorGetName(stairheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Detecting structure from %s\n", DEC_DETECTORNAME);
   //remove empty constraints
   findRelevantConss(scip, detectordata);

   nvars = SCIPgetNVars(scip);
   vars_array = SCIPgetVars(scip);
   ncons = detectordata->nRelevantConss;
   cons_array = detectordata->relevantConss;
   printf("ncons: %i \n", ncons); //debug
   printf("nvars: %i \n", nvars); //debug
   //initialize hash maps for variables: indexvar_old, indexvar_new, varindex_old, varindex_new
   for(i = 0; i < nvars; ++i)
   {
      var = vars_array[i];
      //careful: probindex+1, because '0' is treated as an empty hashmap entry, which causes an error
//      probindex = (SCIPvarGetProbindex(var)+1);
      probindex = i+1;
      assert( ! SCIPhashmapExists(detectordata->indexvar_old, (void*) probindex));
      SCIPhashmapInsert(detectordata->indexvar_old, (void*) probindex, (void*) var);
      assert( ! SCIPhashmapExists(detectordata->indexvar_new, (void*) probindex));
      SCIPhashmapInsert(detectordata->indexvar_new, (void*) probindex, (void*) var);
      assert( ! SCIPhashmapExists(detectordata->varindex_old, (void*) var));
      SCIPhashmapInsert(detectordata->varindex_old, (void*) var, (void*) probindex);
      assert( ! SCIPhashmapExists(detectordata->varindex_new, (void*) var));
      SCIPhashmapInsert(detectordata->varindex_new, (void*) var, (void*) probindex);
   }

   //initialize hash maps for constraints: indexcons_old, indexcons_new, consindex_old, consindex_new
   for(i = 0; i < ncons; ++i)
   {
      cons = cons_array[i];
      //careful: i+1, because '0' is treated as an empty hashmap entry, which causes an error
      assert( ! SCIPhashmapExists(detectordata->indexcons_old, (void*) (i+1)));
      SCIPhashmapInsert(detectordata->indexcons_old, (void*) (i+1), (void*) cons);
      assert( ! SCIPhashmapExists(detectordata->indexcons_new, (void*) (i+1)));
      SCIPhashmapInsert(detectordata->indexcons_new, (void*) (i+1), (void*) cons);
      assert( ! SCIPhashmapExists(detectordata->consindex_old, (void*) cons));
      SCIPhashmapInsert(detectordata->consindex_old, (void*) cons, (void*) (i+1));
      assert( ! SCIPhashmapExists(detectordata->consindex_new, (void*) cons));
      SCIPhashmapInsert(detectordata->consindex_new, (void*) cons, (void*) (i+1));
   }

   //initialize index arrays ibegin_old, iend_old, jbegin_old, jend_old
   printf("1\n");
   rowindices = rowindices_list(scip, detectordata, detectordata->indexcons_old, detectordata->varindex_old);
   columnindices = columnindices_list(scip, detectordata, rowindices);
   printf("2\n");
   formIndexArray(detectordata->ibegin_old, detectordata->iend_old, rowindices);
   formIndexArray(detectordata->jbegin_old, detectordata->jend_old, columnindices);
   printf("3\n");
   //ROC2 algorithm
   do
   {
      RankOrderClustering(scip, detectordata);
   }
   while(IndexArrayChanges(scip, detectordata));
   //arrays jmin, jmax and minV
   detectordata->jmin[0] = detectordata->ibegin_new[0];
   detectordata->jmax[0] = detectordata->iend_new[0];
   detectordata->width[0] = detectordata->iend_new[0] - detectordata->ibegin_new[0];
   for(i = 1; i < ncons; ++i)
   {
      detectordata->width[i] = detectordata->iend_new[i] - detectordata->ibegin_new[i];
      detectordata->jmin[i] = detectordata->ibegin_new[i];
      detectordata->jmax[i] = maximum(detectordata->iend_new[i], detectordata->jmax[i-1]);
      detectordata->minV[i-1]=1 + (detectordata->jmax[i-1] - detectordata->jmin[i]);
   }
   detectordata->n = max_array(detectordata->width, ncons);
   detectordata->v = min_array(detectordata->width, ncons);
   detectordata->tau = round((nvars - detectordata->v)/(detectordata->n - detectordata->v));
   //debug
   printf("<N> <n> <v> <tau>: <%i> <%i> <%i> <%i>\n", nvars, detectordata->n, detectordata->v, detectordata->tau);
   printf("minV = [ ");
   for(i = 0; i < ncons-1; ++i)
   {
      printf("%i ", detectordata->minV[i]);
   }
   printf("]\n");
   printf("jmin = [ ");
   for(i = 0; i < ncons; ++i)
   {
      printf("%i ", detectordata->jmin[i]);
   }
   printf("]\n");
   printf("jmax = [ ");
   for(i = 0; i < ncons; ++i)
   {
      printf("%i ", detectordata->jmax[i]);
   }
   printf("]\n");

   rowsWithConstriction(detectordata);
   list_prin(detectordata->rowsWithConstrictions, printint); //debug

   //fill detectordata for a single block for testing
   detectordata->blocks = 1;
   detectordata->varsperblock[0] = vars_array;
   detectordata->nvarsperblock[0] = nvars;
   detectordata->consperblock[0] = SCIPgetConss(scip);
   detectordata->nconsperblock[0] = ncons;
   detectordata->linkingconss = NULL;
   detectordata->nlinkingconss = 0;
   detectordata->found = TRUE;
   SCIP_CALL(copyDetectorDataToDecomp(scip, detectordata, detectordata->decdecomp));

   //deallocate memory
   list_delete_nested(rowindices);
   list_delete_nested(columnindices);

   *result = SCIP_SUCCESS;
   return SCIP_OKAY;
}

/** set the decomp structure */
static
DEC_DECL_SETSTRUCTDECOMP(StairheurSetDecomp)
{
   DEC_DETECTOR* stairheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   assert(decdecomp != NULL);
   stairheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(stairheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(stairheur), DEC_DETECTORNAME) == 0);
   SCIPdebugMessage("Setting decdecomp\n");
   detectordata->decdecomp = decdecomp;
}

static
DEC_DECL_GETPRIORITY(getPriority)
{
   DEC_DETECTOR* arrowheur;
   DEC_DETECTORDATA* detectordata;
   assert(scip != NULL);
   arrowheur = DECfindDetector(scip, DEC_DETECTORNAME);
   detectordata = DECdetectorGetData(arrowheur);
   assert(detectordata != NULL);

   assert(strcmp(DECdetectorGetName(arrowheur), DEC_DETECTORNAME) == 0);
   return detectordata->priority;
}

/** creates the stairheur presolver and includes it in SCIP */
SCIP_RETCODE SCIPincludeDetectionStairheur(
   SCIP*                 scip              /**< SCIP data structure */

   )
{
   DEC_DETECTORDATA *detectordata;
   assert(scip != NULL);

   SCIP_CALL(SCIPallocMemory(scip, &detectordata));

   assert(detectordata != NULL);
   detectordata->found = FALSE;
   detectordata->decdecomp = NULL;
   SCIP_CALL(DECincludeDetector(scip, DEC_DETECTORNAME, detectordata, detectAndBuildStair, StairheurSetDecomp, initStairheur, exitStairheur, getPriority));

   /* add stairheur presolver parameters */
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/maxblocks", "The maximal number of blocks", &detectordata->maxblocks, FALSE, DEFAULT_MAXBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/minblocks", "The minimal number of blocks", &detectordata->minblocks, FALSE, DEFAULT_MINBLOCKS, 2, 1000000, NULL, NULL));
   SCIP_CALL(SCIPaddIntParam(scip, "stairheur/priority", "Priority of detector", &detectordata->priority, FALSE, DEFAULT_PRIORITY, INT_MIN, INT_MAX, NULL, NULL));
   return SCIP_OKAY;
}
