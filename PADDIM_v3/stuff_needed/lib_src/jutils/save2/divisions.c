#include <stdio.h>
#include <stdlib.h>

#include "inout.h"



struct nlist *
calcdivisions(int xdim, int ydim, int zdim)
{
int subimg = 1,i, laenge, checksum=0, chunk;
struct nlist *list, *olist;

     olist =(struct nlist*)  malloc(sizeof(struct nlist)*1);

  list = olist; // save start of list

  laenge = ydim*zdim;

  while((ydim)*(zdim)/subimg > 256*256-1) subimg++;

   chunk = ydim*zdim/subimg;

  while(laenge > 256*256-1) {
   
     laenge -= chunk;

     list->num = chunk;

     list->next =(struct nlist*)  malloc(sizeof(struct nlist)*1);

     list = list->next;

  };

     list->num = laenge;
     list->next = NULL;
     printf("last  chunk %i \n",laenge);

   list = olist; //reset pointer to start

   while(list != NULL) {

        checksum +=list->num;
        printf(" mylist %i \n", list->num);
        list = list->next;
    }  

     printf(" checksum %i : orgi %i\n", checksum, ydim*zdim);

   return(olist);
}

void
freemylist(struct nlist* olist)
{
struct nlist *list;

olist = olist->next;

   while(olist != NULL) {

        list = olist;
        olist = olist->next;
        free(list);
    }

}


void
test_main()
{

struct nlist *mylist, *plist;


   mylist = calcdivisions(512, 380, 421);

 plist = mylist;

   while(plist != NULL) {
        printf(" plist %i \n", plist->num);
        plist = plist->next;
    }  

  freemylist(mylist);
   
}
