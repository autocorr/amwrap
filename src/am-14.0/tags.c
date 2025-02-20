/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* tags.c                         S. Paine rev. 2024 April 16
*
* Tags are strings which can be used to identify groups of
* layers.  The functions here maintain a table of tag
* strings.
************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static int num_tag_strings = 0;
static char **tag_strings  = NULL;

static int initialize_tag_string_table(void);

/***********************************************************
* int initialize_tag_string_table(void)
*
* Purpose:
*   Initializes the tag string table, setting the first
*   entry to point at an empty string.
*
* Return:
*   0 on success, 1 on failure
************************************************************/

static int initialize_tag_string_table(void)
{
    if ((tag_strings = (char**)malloc(sizeof(char*))) == NULL)
        return 1;
    if ((tag_strings[0] = (char*)malloc(sizeof(char))) == NULL) {
        free(tag_strings);
        return 1;
    }
    tag_strings[0][0] = '\0';
    num_tag_strings = 1;
    return 0;
}   /* initialize_tag_string_table() */


/***********************************************************
* void free_tag_string_table(void)
*
* Purpose:
*   Frees the tag string table
************************************************************/

void free_tag_string_table(void)
{
    if (tag_strings == NULL)
        return;
    while (num_tag_strings)
        free(tag_strings[--num_tag_strings]);
    free(tag_strings);
    tag_strings = NULL;
    return;
}   /* free_tag_string_table() */


/***********************************************************
* int get_num_tag_strings(void)
*
* Purpose:
*   Returns the number of tag strings in the tag string
*   table
************************************************************/

int get_num_tag_strings(void)
{
    if (tag_strings == NULL) {
        if (initialize_tag_string_table())
            return 0;
    }
    return num_tag_strings;
}   /* get_num_tag_strings() */


/***********************************************************
* char* tag_string(const int n)
*
* Purpose:
*   Returns a pointer to the tag string with index n.
*
* Arguments:
*   int n - table index of tag string
*
* Return:
*   pointer to tag string, or NULL if n is outside the table
*   range or an error occurs
************************************************************/

char *tag_string(const int n)
{
    if (tag_strings == NULL) {
        if (initialize_tag_string_table())
            return NULL;
    }
    if (n < 0 || n >= num_tag_strings || tag_strings == NULL)
        return NULL;
    return tag_strings[n];
}   /* tag_string() */


/***********************************************************
* int tag_string_index(const char *s)
*
* Purpose:
*   Looks for a string s in the tag string table.  If found,
*   the index of the string is returned.  If s is not
*   already in the table, it is assigned the next available
*   index number, which is returned.
*
* Arguments:
*   char *s - pointer to a tag string to look up in the
*             table
*
* Return:
*   table index of tag string, or -1 on error.
************************************************************/

int tag_string_index(const char *s)
{
    void *tptr;
    int i;

    if (tag_strings == NULL) {
        if (initialize_tag_string_table())
            return -1;
    }
    for (i = 0; i < num_tag_strings; ++i) {
        if (!strcmp(s, tag_strings[i]))
            return i;
    }
    if ((tptr = realloc(tag_strings,
                    (num_tag_strings + 1) * sizeof(char*))) == NULL)
        return -1;
    tag_strings = (char**)tptr;
    if ((tag_strings[num_tag_strings] = malloc(
                    (1 + strlen(s)) * sizeof(char))) == NULL)
        return -1;
    strcpy(tag_strings[num_tag_strings], s);
    return num_tag_strings++;
}   /* tag_string_index() */


#ifdef UNIT_TEST

/***********************************************************
* int main(int argc, char **argv)
*
* Purpose:
*   Simple stand-alone test for these functions, which
*   inserts strings from the command line and writes them
*   back.
*
* Example:
*   $ gcc -D UNIT_TEST tags.c
*   $ ./a.out one times one equals one
*   1
*   2
*   1
*   3
*   1
*      0 
*      1 one
*      2 times
*      3 equals
*   $ 
************************************************************/

int main(int argc, char **argv) {
    int i;

    for (i = 1; i < argc; ++i)
        printf("%d\n", tag_string_index(argv[i]));
    for (i = 0; i < num_tag_strings; ++i)
        printf("%4d %s\n", i, tag_string(i));
    free_tag_string_table();
    return 0;
} /* main() */

#endif /* UNIT_TEST */
