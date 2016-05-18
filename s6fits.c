#include <string.h>
#include <stdlib.h>
#include "fitsio.h"
#include "s6fits.h"
#include <vector>

/*function declarations*/
//int is_metadata (fitsfile * fptr, int hdupos);

time_t get_time (fitsfile * fptr, int * status);

double get_RA (fitsfile * fptr, int * status);

int get_bors (fitsfile * fptr, int * status);

double get_dec (fitsfile * fptr, int * status);

int get_nhits (fitsfile * fptr, int * status);

int get_missedpk (fitsfile * fptr, int * status);

int is_metadata (time_t time);

/*get_s6data will populate a s6hits_t data type according to the data specs
 * given by the user in the s6dataspec argument. For more details please see
 * s6hits.h, which contains the two data types used here.*/

int get_s6data(s6dataspec_t * s6dataspec) 
{
  fitsfile *fptr;
  /*status must be initialized to zero!*/
  int status = 0;
  int single = 0, hdupos, nkeys, ii;
  int hdunum, hdutype, ncols;
  char * filename = s6dataspec->filename;

  if (!fits_open_file(&fptr, filename, READONLY, &status))
  {
    //fits_get_hdu_num(fptr, &hdupos);
    
    for (; !status; hdupos++) 
    {
      /* calculated and saved before we make the hit to save us the hassle of
 * grabbing it for each new hit in the binary table*/
      time_t time = get_time(fptr, &status);
      double ra = get_RA(fptr, &status);
      int bors = get_bors(fptr, &status);
      double dec = get_dec(fptr, &status);
      int nhits = get_nhits(fptr, &status);
      int missedpk = get_missedpk(fptr, &status); 
      /*enter bintable loop, declare s6fits in here, adding in time, ra, bors,
      * channel, etc inside. Maybe make a function to add metadata*/
      //assume for now loop has been entered
      if (is_metadata(time))
      {
        /*start of tablist.c, for figuring out bintable stuff ONLY!*/ 
        int ncols, anynul;
        long nrows; 

        fits_get_num_rows(fptr, &nrows, &status);
        fits_get_num_cols(fptr, &ncols, &status); 
        
        float detpow_array[nrows];
        fits_read_col (fptr, TFLOAT, 1, 1, 1, nrows, NULL,
            &detpow_array, &anynul, &status); 

        float meanpow_array[nrows];
        fits_read_col (fptr, TFLOAT, 2, 1, 1, nrows, NULL,
            &meanpow_array, &anynul, &status); 

        unsigned short coarch_array[nrows];
        fits_read_col (fptr, TUSHORT, 3, 1, 1, nrows, NULL,
            &coarch_array, &anynul, &status); 
        
        unsigned long finechan_array[nrows];
        fits_read_col (fptr, TULONG, 4, 1, 1, nrows, NULL,
            &finechan_array, &anynul, &status); 

        for (int i = 0; i < nhits; i++) 
        {
          s6hits_t hit;
          hit.time = time;
          hit.ra = ra;
          hit.bors = bors;
          hit.dec = dec;
          hit.nhits = nhits;
          hit.missedpk = missedpk;
          hit.detected_power = detpow_array[i];
          hit.mean_power = meanpow_array[i];
          hit.fine_channel_bin = finechan_array[i];
          hit.coarse_channel_bin = coarch_array[i];
          s6dataspec->s6hits.push_back(hit);
        } 
      }
      
      /*moves to next HDU*/
      fits_get_hdrspace(fptr, &nkeys, NULL, &status); 
      fits_movrel_hdu(fptr, 1, NULL, &status);
    }  
  }    
  status = 0;
  fits_close_file(fptr, &status);
 
  if (status) fits_report_error(stderr, status);
  
  s6dataspec->errorcode = status; 
  
  return (status);
}

/*gets time found in the binary table header block*/
time_t get_time (fitsfile * fptr, int * status)
{
  time_t time;
  fits_read_key(fptr, TINT, "TIME", &time, NULL, status);
  if (*status == KEY_NO_EXIST) time = -1, *status=0;
  return (time);
}

/*gets RA found in the binary table header block*/
double get_RA (fitsfile * fptr, int * status) 
{
  double ra;
  fits_read_key(fptr, TDOUBLE, "RA", &ra, NULL, status);
  if (*status == KEY_NO_EXIST) ra = -1, *status=0;
  return (ra);
}

/*gets beampol found in the binary table header block*/
int get_bors (fitsfile * fptr, int * status)
{
  int bors;
  fits_read_key(fptr, TINT, "BEAMPOL", &bors, NULL, status);
  if (*status == KEY_NO_EXIST) bors = -1, *status=0;
  return (bors);
}

/*gets declanation found in the binary table header block*/
double get_dec (fitsfile * fptr, int * status)
{
  double dec;
  fits_read_key(fptr, TDOUBLE, "DEC", &dec, NULL, status);
  if (*status == KEY_NO_EXIST) dec = -1, *status=0;
  return (dec);
}

/*gets nhits found in the binary table header block. The individual hit will
 * include the number of hits in the 'hit block' it's in. Kinda meta huh?*/
int get_nhits (fitsfile * fptr, int * status)
{
  int nhits;
  fits_read_key(fptr, TINT, "NHITS", &nhits, NULL, status);
  if (*status == KEY_NO_EXIST) nhits = -1, *status=0;
  return (nhits);
}


/*gets missed packets found in the binary table header block*/
int get_missedpk (fitsfile * fptr, int * status)
{
  int missedpk;
  fits_read_key(fptr, TINT, "MISSEDPK", &missedpk, NULL, status);
  if (*status == KEY_NO_EXIST) missedpk = -1, *status=0;
  return (missedpk);
}

int is_metadata (time_t time)
{
  if (time != -1) return 1;
  else return 0; 
}

void print_hits_structure (std::vector<s6hits_t> s6hits) 
{
    for (std::vector<s6hits_t>::iterator hit = s6hits.begin() ; 
         hit != s6hits.end(); ++hit)
    {
      /*header metadata*/
      printf("time: %d\n", hit->time);
      printf("RA: %g\n", hit->ra);
      printf("beampol: %d\n", hit->bors);
      printf("dec: %g\n", hit->dec);
      //printf("channel: %d\n", hit->channel);
      //printf("nhits: %d\n", hit->nhits);
      printf("missedpk: %d\n", hit->missedpk);
      /*binary table data*/
      printf("detected power: %f\n", hit->detected_power);
      printf("mean power: %f\n", hit->mean_power);
      printf("fine channel bin: %lu\n", hit->fine_channel_bin);
      printf("coarse channel bin: %hu\n", hit->coarse_channel_bin);
      /*calculated frequencies*/
      /*unimplemented 
      printf("ifreq: %g\n", hit->ifreq);
      printf("rfreq: %g\n", hit->rfreq); 
      */
      printf("\n");
    }
}

void print_hits_table (std::vector<s6hits_t> s6hits)
{
  printf("%15s", "Time");
  printf("%10s", "RA");
  printf("%10s", "BorS");
  printf("%15s", "DEC");
  printf("%10s", "MSDPK");
  printf("%25s", "DETPOW");
  printf("%20s", "MEANPOW");
  printf("%15s", "FINECHAN");
  printf("%15s\n", "COARCHID"); 
  for (std::vector<s6hits_t>::iterator hit = s6hits.begin() ; 
         hit != s6hits.end(); ++hit)
  {
    printf("%15d", hit->time);
    printf("%10g", hit->ra);
    printf("%10d", hit->bors);
    printf("%15g", hit->dec);
    printf("%10d", hit->missedpk);
    printf("%25f", hit->detected_power);
    printf("%20f", hit->mean_power);
    printf("%15lu", hit->fine_channel_bin);
    printf("%15hu\n", hit->coarse_channel_bin);
  }
}

