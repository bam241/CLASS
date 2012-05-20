#include "NextName.hxx"
#include <stdio.h>

char *NextName() {
  static long NameNo;
  static char Name[15];

  snprintf(Name,15,"Name_%ld",NameNo);
  NameNo++;
  return Name;
}
