/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* nscale.c                   S. Paine rev. 2019 September 26 
*
* Functions to support the Nscale facility.  Each Nscale
* statement in a configuration file generates an entry in
* a linked list maintained and accessed by the functions in
* this file.
************************************************************/

#include <stdio.h>
#include <stdlib.h>

#include "nscale.h"

static Nscale_list_t *head = NULL;


/***********************************************************
* Nscale_list_t *create_Nscale_list_entry(
*         const int col_typenum,
*         const int tagnum,
*         const double Nscale)
*
* Purpose:
*   Creates a new entry in the Nscale list.  If an entry
*   with the same col_typenum and tagnum already exists,
*   its Nscale parameter is overwritten.
*
* Arguments:
*   int col_typenum - a column type number
*   int tagnum      - index number of a tag string
*
* Return:
*   Pointer to the newly created or revised entry, or
*   NULL on failure.
************************************************************/

Nscale_list_t *create_Nscale_list_entry(
        const int col_typenum,
        const int tagnum,
        const double Nscale)
{
    Nscale_list_t *ptr, *tail = NULL;

    for (ptr = head; ptr != NULL; ptr = ptr->next) {
        if (col_typenum == ptr->col_typenum && tagnum == ptr->tagnum) {
            /* Found duplicate entry; replace Nscale parameter. */
            ptr->Nscale = Nscale;
            return ptr;
        }
        tail = ptr;
    }
    /* Create a new entry at the end of the list. */
    if ((ptr = (Nscale_list_t*)malloc(sizeof(Nscale_list_t))) == NULL)
        return NULL;
    ptr->next = NULL;
    ptr->col_typenum = col_typenum;
    ptr->tagnum = tagnum;
    ptr->Nscale = Nscale;
    if (head == NULL)
        head = ptr;
    else
        tail->next = ptr;
    return ptr;
}   /* create_Nscale_list_entry() */


/***********************************************************
* Nscale_list_t *find_Nscale_list_entry(
*         const int col_typenum,
*         const int tagnum)
*
* Purpose:
*   Find the address of an entry in the Nscale list.
*
* Arguments:
*   int col_typenum - a column type number
*   int tagnum - index number of a tag string
*
* Return:
*   Pointer to an existing entry, or NULL if none exists
************************************************************/

Nscale_list_t *find_Nscale_list_entry(
        const int col_typenum,
        const int tagnum)
{
    Nscale_list_t *ptr;

    for (ptr = head; ptr != NULL; ptr = ptr->next) {
        if (col_typenum == ptr->col_typenum && tagnum == ptr->tagnum)
            return ptr;
    }
    return NULL;
}   /* find_Nscale_list_entry() */


/***********************************************************
* double lookup_Nscale(
*         const int col_typenum,
*         const int tagnum)
*
* Purpose:
*   Gets the Nscale factor for a given column type and layer
*   tag.  This is done by multiplying the scale factors for
*   all matching entries in the list, so that, for example,
*   a pair of configuration file statements such as
*
*       Nscale h2o 0.8
*       Nscale troposphere h2o 0.8
*
*   will have the effect of scaling h2o by a factor 0.64 on
*   all layers tagged "troposphere", and by a factor 0.8 on
*   all other layers.
*
* Arguments:
*   int col_typenum - a column type number
*   int tagnum - index number of a tag string
*
* Return:
*   Nscale factor.  A value of 1.0 is returned if no match
*   other than the default is found in the list.
************************************************************/

double lookup_Nscale(
        const int col_typenum,
        const int tagnum)
{
    Nscale_list_t *ptr;
    double Nscale = 1.0;

    for (ptr = head; ptr != NULL; ptr = ptr->next) {
        if ((ptr->tagnum == 0 || ptr->tagnum == tagnum) &&
            (ptr->col_typenum == 0 || ptr->col_typenum == col_typenum))
            Nscale *= ptr->Nscale;
    }
    return Nscale;
}   /* lookup_Nscale() */


/***********************************************************
* void free_Nscale_list(void)
*
* Purpose:
*   Frees the Nscale list
************************************************************/

void free_Nscale_list(void)
{
    Nscale_list_t *ptr, *lptr;

    ptr = head;
    while (ptr != NULL) {
        lptr = ptr;
        ptr = ptr->next;
        free(lptr);
    }
    return;
}   /* free_Nscale_list() */


/***********************************************************
* Nscale_list_t *Nscale_list_head(void)
*
* Purpose:
*   Returns a pointer to the head of the Nscale list, or
*   NULL if the list is empty.
************************************************************/

Nscale_list_t *Nscale_list_head(void)
{
    return head;
}   /* Nscale_list_head() */


#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Tests the functions in this file.
************************************************************/

int main(int argc, char **argv) {
    printf("lookup_Nscale(5, 5) = %g\n", lookup_Nscale(5, 5));
    printf("lookup_Nscale(5, 3) = %g\n", lookup_Nscale(5, 3));
    printf("lookup_Nscale(5, 0) = %g\n", lookup_Nscale(5, 0));
    printf("lookup_Nscale(3, 0) = %g\n", lookup_Nscale(3, 0));
    printf("lookup_Nscale(3, 1) = %g\n", lookup_Nscale(3, 1));
    printf("\n");
    printf("create_Nscale_list_entry(5, 3, 0.9) = %p\n",
        create_Nscale_list_entry(5, 3, 0.9));
    printf("create_Nscale_list_entry(5, 5, 0.8) = %p\n",
        create_Nscale_list_entry(5, 5, 0.8));
    printf("create_Nscale_list_entry(3, 0, 0.9) = %p\n",
        create_Nscale_list_entry(3, 0, 0.9));
    printf("\n");
    printf("lookup_Nscale(5, 5) = %g\n", lookup_Nscale(5, 5));
    printf("lookup_Nscale(5, 3) = %g\n", lookup_Nscale(5, 3));
    printf("lookup_Nscale(5, 0) = %g\n", lookup_Nscale(5, 0));
    printf("lookup_Nscale(3, 0) = %g\n", lookup_Nscale(3, 0));
    printf("lookup_Nscale(3, 1) = %g\n", lookup_Nscale(3, 1));
    printf("\n");
    printf("create_Nscale_list_entry(5, 0, 0.7) = %p\n",
        create_Nscale_list_entry(5, 0, 0.7));
    printf("\n");
    printf("lookup_Nscale(5, 5) = %g\n", lookup_Nscale(5, 5));
    printf("lookup_Nscale(5, 3) = %g\n", lookup_Nscale(5, 3));
    printf("lookup_Nscale(5, 0) = %g\n", lookup_Nscale(5, 0));
    printf("lookup_Nscale(3, 0) = %g\n", lookup_Nscale(3, 0));
    printf("lookup_Nscale(3, 1) = %g\n", lookup_Nscale(3, 1));
    printf("\n");
    printf("create_Nscale_list_entry(5, 5, 0.6) = %p\n",
        create_Nscale_list_entry(5, 5, 0.6));
    printf("lookup_Nscale(5, 5) = %g\n", lookup_Nscale(5, 5));
    printf("lookup_Nscale(5, 3) = %g\n", lookup_Nscale(5, 3));
    printf("lookup_Nscale(5, 0) = %g\n", lookup_Nscale(5, 0));
    printf("lookup_Nscale(3, 0) = %g\n", lookup_Nscale(3, 0));
    printf("lookup_Nscale(3, 1) = %g\n", lookup_Nscale(3, 1));
    free_Nscale_list();
    return 0;
} /* main() */

#endif /* UNIT_TEST */
