#include <extUtil.h>


int checkIfNumber( char *string )
{
  int i;
  for(i=0; i<strlen(string); i++) if((string[i]<48)||(string[i]>57)) return(0);
  return(1);
}

