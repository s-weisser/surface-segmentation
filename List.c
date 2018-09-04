/**************************************************************
*
* PROJECT:      Surface Segmentation
* FILE:         List.c
* AUTORS:       Steffen Weisser
* AFFILIATION:  Saarland University, Germany
* DATE:         2018
*
* Copyright (C) 2018 Steffen Weisser - All Rights Reserved
* You may use, distribute and modify this code under the
* terms of the MIT license.
*
**************************************************************/

#include <assert.h>
#include <stdlib.h>
#include "List.h"

void ListInit(List *list) {
  list->numEntries = 0;
  list->first = NULL;
  list->current = NULL;
  list->last = NULL;
}

void ListFree(List *list) {
  ListEntry *entry, *next;

  entry = list->first;
  while (entry != NULL) {
    next = entry->next;
    free(entry);
    entry = next;
  }

  list->numEntries = 0;
  list->first = NULL;
  list->current = NULL;
  list->last = NULL;
}

void ListDel(List *list, void *entry) {
  ListEntry *prev, *next;

  if (list->first == NULL)
    return;

  next = list->first;
  if (next->entry == entry) {     // delete first entry?
    list->first = list->first->next;
    list->numEntries--;
    free(next);
    return;
  }

  do {
    prev = next;
    next = next->next;
    if (next == NULL)    // entry not in list?
      return;
  } while (next->entry != entry);

  // entry found, delete it
  prev->next = next->next;
  list->numEntries--;
  free(next);
}

void ListAdd(List *list, void *entry) {

  ListEntry *newEntry = (ListEntry*) calloc(1, sizeof(ListEntry));
  assert(newEntry!=NULL);

  (list->numEntries)++;

  newEntry->entry = entry;

  if (list->first == NULL) {
    list->first = newEntry;
    list->last = newEntry;
  }
  else {
    list->last->next = newEntry;
    list->last = newEntry;
  }
}

void ListAddAtPos(List *list, void *entry, long pos) {
  long i;

  if (pos >= list->numEntries) 
    ListAdd(list, entry);             // add at the end of the list
  else {
    ListEntry *curEntry, *nextEntry;
    ListEntry *newEntry = (ListEntry*) calloc(1, sizeof(ListEntry));
    assert(newEntry!=NULL);
  
    (list->numEntries)++;
  
    newEntry->entry = entry;
  
    if (pos==0) {                     // add at the beginning of the list
      newEntry->next = list->first;
      list->first = newEntry;
    }
    else {                            // add somewhere in the middle
      i=1;
      curEntry = list->first;
      nextEntry = curEntry->next;
      while(i!=pos) {
        curEntry = nextEntry;
        nextEntry = curEntry->next;
        i++;
      }
      newEntry->next = curEntry->next;
      curEntry->next = newEntry;
    }
  }
}

void *ListGetFirst(List *list) {
  list->current = list->first;

  if (list->current != NULL)
    return list->current->entry;
  else
    return NULL;
}

void *ListGetLast(List *list) {
  list->current = list->last;

  if (list->current != NULL)
    return list->current->entry;
  else
    return NULL;
}

void *ListGetNext(List *list) {

  if (list->current != NULL) 
    list->current = list->current->next;
  else 
    return NULL;

  if (list->current != NULL) 
    return list->current->entry;
  else 
    return NULL;
}

long ListGetLength(List *list) {
  return list->numEntries;
}
