/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         List.h
* AUTORS:       Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#ifndef LIST_H
#define LIST_H

/**************** declaration of variables and structs ***********************/

/** entry of dynamic list */
typedef struct sListEntry{
  struct sListEntry *next;
  void *entry;
} ListEntry;

/** dynamic list */
typedef struct {
  long numEntries;
  ListEntry *first;
  ListEntry *current;
  ListEntry *last;
} List;

/**************** declaration of global functions ***************************/
void ListInit(List *list);
void ListFree(List *list);
void ListDel(List *list, void *entry);
void ListAdd(List *list, void *entry);
void ListAddAtPos(List *list, void *entry, long pos);
void *ListGetFirst(List *list);
void *ListGetNext(List *list);
void *ListGetLast(List *list);
long ListGetLength(List *list);

#endif
