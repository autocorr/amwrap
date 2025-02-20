/***********************************************************
* Smithsonian Astrophysical Observatory
* Submillimeter Receiver Laboratory
* am
*
* tags.h                      S. Paine rev. 2012 December 12
*
* Declarations for tags.c
************************************************************/

#ifndef AM_TAGS_H
#define AM_TAGS_H


int free_tag_string_table(void);
int get_num_tag_strings(void);
char *tag_string(const int);
int tag_string_index(const char *s);

#endif /* AM_TAGS_H */
