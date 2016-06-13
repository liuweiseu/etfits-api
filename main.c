#include <stdio.h>
#include <time.h>
#include <vector>
#include "s6fits.h"

void get_num_hits(std::vector<s6hits_t> s6hits)
{
  int num_hits = 0;
  for (std::vector<s6hits_t>::iterator it = s6hits.begin(); it != s6hits.end(); ++it)
  {
    num_hits += 1;
  }
  printf("number of hits: %d\n", num_hits); 
}

int main(int argc, char *argv[]) 
{
  s6dataspec_t s6dataspec; 
  char * filename = argv[1];
  //std::vector<s6hits_t> hits;
  s6dataspec.filename = filename;
  //s6dataspec.s6hits = hits;
  //s6dataspec.sortby_freq = 0;
  s6dataspec.sortby_time = 0;
  //std::vector<int> bors;
  //bors.push_back(3);
  //bors.push_back(6);
  //s6dataspec.bors = bors;
  s6dataspec.sortby_bors = 0;
  s6dataspec.threshold = 20.0; 
  get_s6data(&s6dataspec);
  print_hits_table(s6dataspec.s6hits);
  //if (argc == 3) get_num_hits(s6dataspec.s6hits);
  //time_t time = get_time_over_file (filename);
  //int num_hits = get_hits_over_file (filename);
  //printf("seconds over file: %d\n", (int) time);
  //printf("number of hits: %d\n", num_hits); 
  return 0;
}
