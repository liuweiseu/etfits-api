#include <stdio.h>
#include <time.h>
#include <vector>
#include "s6fits.h"

int main(int argc, char *argv[]) 
{
  s6dataspec_t s6dataspec; 
  s6dataspec.filename = argv[1];
  s6dataspec.sortby_bors = 0;
  s6dataspec.threshold = 0; 
  get_s6data(&s6dataspec);
  return 0;
}
