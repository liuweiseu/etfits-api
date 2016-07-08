"""
This is an example of Python code using the s6fits API to create list of hits,
where each element is a dictionary of every keyword in the hit. The function 
itself is relatively useless, but it should provide a fairly good idea of how 
to use the (somewhat basic) functions and classes that have been implemented. 

s6hits_t and s6dataspec_t are C structs that contain, respectively, the data for
one individual hit, and the specs for the hits that you want to extract. 
The structs and variables are listed below. Commented out variables have not 
been implemented and cannot be used.

Also note that there are no constructors implemented at the moment, so you must
assign each variable individually.

As you can see, s6hits is of type std::vector<s6hits_t>, or as known to Python,
an s6Vec. s6Vec functions are the same as those in std::vector, so if you are
unsure of a certain function, you can look at the std::vector doc page:

http://en.cppreference.com/w/cpp/container/vector

typedef struct 
{
  time_t time;
  double ra;
  int bors;
  double dec;
  int missedpk;
  float detected_power;
  float mean_power;
  unsigned long fine_channel_bin;
  unsigned short coarse_channel_bin;
  double ifreq;
  double rfreq;
} s6hits_t;

typedef struct
{
  char * filename;
  std::vector<s6hits_t> s6hits;
  int sortby_freq;
  int sortby_time;
  int sortby_bors;
  float threshold;
  std::vector<int> bors;
  std::vector<int> channels;
  int errorcode;
} s6dataspec_t;

Upcoming Python specific API features:
-throwing an error when a user goes out of bounds on the s6Vec class. At the
moment, users are able to continue past the end of the vector until the program
segfaults. There won't be any segfaults should it be used as a proper iterable,
but in some cases (such as using indexes), it will continue past the vector's
end.
"""

import s6fits
import sys

def make_hits_list(hits):
  """
  Creates the main list of hits and navigates through the fits file
  """
  hits_list = []
  
 
  """
  so we have a saved s6hits vector, which functions as an iterator.
  """ 
  for hit in hits:
    hit_dict = make_dictionary(hit)   
    hits_list.append(hit_dict)
  return hits_list

def make_dictionary(hit):
  hit_dict = {}
  hdu_keys = [
          "time", 
          "ra", 
          "bors", 
          "dec", 
          "missedpk"
          ]
  bin_table_keys = [
          "coarse_channel_bin", 
          "fine_channel_bin", 
          "detected_power", 
          "mean_power"
          ]
  keys = hdu_keys + bin_table_keys
  for key in keys:
    val = getattr(hit, key)    
    hit_dict[key] = val
  return hit_dict

"""okay, well that's simple enough, but now I really want to manipulate the
vector more. For example, what if I just want to get a subset of the hits that I
just got? Or I really just hate one hit in particular, how should I get rid of 
it? 

Find out next time when I implement a better looking main.""" 

def remove_single_hit(s6hits, index):
  """removes a hit at a specific index."""

  # .begin() returns an iterator beginning at the first element.
  # note if you want to delete the first element, you still need to add 0
  head = s6hits.begin()
  s6hits.erase(head + index)
   
def remove_range_of_hits(s6hits, index_range):
  head = s6hits.begin()
  """this code seems  a little unintuitive, but much like python's range 
function, std::vectors erase is not inclusive of the last element you include. 
So if you are erasing [1, 5), then you'll need to make sure that erase is 
receiving the numbers 1 and 5"""
  s6hits.erase(head + index_range[0], head + index_range[-1] + 1)

def get_single_hit(s6hits, index):
  """vectors are also indexable, so getting a single hit is an easy task"""
  return s6hits[index]

def get_range_of_hits(s6hits, index_range):
  """gets a range of hits from s6hits and puts it into a new vector"""
  vec = s6fits.s6Vector()
  for i in index_range:
    vec.push_back(s6hits[i]) 
  return vec

def combine_vectors(hitlist1, hitlist2):
  """combines two vectors into a new one and returns. I wouldn't try merging the
two vectors together using one as the return vector for now."""
  vec = s6fits.s6Vector()
  vec.insert(hitlist1.begin(), s6Vector)
  return vec

def main():
  if sys.argv is None:
    print "provide a valid fits file to run this"
    exit()
  """
  initialize your dataspec_t object. For now you must initialize it with no
  arguments and fill them in subsequent lines. If you only specify filename, the
  program will assume you want all the hits in the file you provide.
  """ 
  dataspec = s6fits.s6dataspec_t()
  dataspec.filename = sys.argv[1]
  dataspec.sortby_bors = 0
  dataspec.sortby_time = 0
  dataspec.sortby_ifreq = 0
  datspec.sortby_rfreq = 0
  """
  get_s6data() will either create a new s6fits vector, or if your dataspec
  includes a vector, it will add hits to it. Right now we have no vector of hits
  we want to add to so we'll just leave it empty.
  """  
  s6fits.get_s6data(dataspec)
  s6fits.print_hits_table(dataspec.s6hits)
  #make_hits_list(dataspec.s6hits)
"""
  remove_single_hit(dataspec.s6hits, 10)   
  remove_range_of_hits(dataspec.s6hits, range(0,5))
  index_range = range(0, 100) 
  v = get_range_of_hits(dataspec.s6hits, index_range)
  dataspec2 = s6fits.s6dataspec_t() 
  dataspec2.filename = sys.argv[2]
  s6fits.get_s6data(dataspec2)
  u = s6fits.s6Vector()
  combine_vectors(dataspec.s6hits, dataspec2.s6hits)
"""

if __name__ == "__main__":
  main()  

