#ifndef S6FITS_H_INCLUDED
#define S6FITS_H_INCLUDED

#include <stdio.h>
#include <time.h>
#include <vector>
#include <stdint.h>

typedef struct 
{
  double julian_date;
  time_t unix_time;
  int bors;
  double ra;
  double dec;
  //int missedpk;
  double detected_power;
  double mean_power;
  int32_t fine_channel_bin;
  unsigned short coarse_channel_bin;
  int user_flag;
  double ifreq;
  double rfreq;
} s6hits_t;

typedef struct
{
  double julian_date;
  time_t unix_time;
  double ra;
  double dec;
  int beam;
  double coarse_channel_bin;
  double power_x;
  double power_y;
  //double ifreq;
  //double rfreq;
} s6ccpowers_t;

typedef struct
{
  double julian_date;
  time_t unix_time;
  double ra;
  int bors;
  double dec;
  int nhits;
  int missedpk;  
} s6hitsheader_t;

typedef struct s6dataspec_t
{
  char * filename;
  //vector where hits requested will be added onto.
  std::vector<s6hits_t> s6hits;
  std::vector<s6ccpowers_t> s6ccpowers;
  std::vector<s6hitsheader_t> s6hitsheaders;
  //sortby_ifreq and sorty_rfreq can't both be > 0, if they are maybe it should
  //default to one or the other (it'll mostly be the same order regardless)
  int sortby_ifreq;
  int sortby_rfreq;
  int sortby_time;
  int sortby_bors;
  int total_missedpk;
  float threshold;
  std::vector<int> bors;
  //same as above but with channels
  std::vector<int> channels;
  int errorcode;
  //status in cfitsio
  int ciftsio_error;
  /*s6dataspec_t() : sortby_ifreq(0),
                   sortby_rfreq(0),
                   sortby_time(0),
                   sortby_bors(0) {}*/
} s6dataspec_t;



int get_s6data (s6dataspec_t * s6dataspec);

int get_s6ccpowers (s6dataspec_t * s6dataspec);

int get_s6hitsheaders (s6dataspec_t * s6dataspec);

time_t get_time_over_file (char * filename);

int get_hits_over_file (char * filename);

void print_hits_structure (std::vector<s6hits_t> s6hits);

void print_hits_table (std::vector<s6hits_t> s6hits);

void print_hits_header_table (std::vector<s6hitsheader_t> s6hitsheaders);

void print_ccpowers_table (std::vector<s6ccpowers_t> s6ccpowers);

#endif
