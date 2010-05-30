// Perform correlations on FRET data.
// See "usage" message below, or try ./favia -h
//

using namespace std;

#include <sys/types.h>          /* for stat */
#include <sys/stat.h>           /* for stat */
#include <unistd.h>             /* for stat */
#include "Getopt.h"
#include <stdlib.h>             /* atoi, getenv */
#include <values.h>             /* INT_MAX */
#include <math.h>               /* log10() */

#include <iostream>
#include <fstream>
#include "pakdot.h"
#include <boost/format.hpp>

#include <assert.h>

void usage(ostream& foo) {
  foo << 
"Perform correlations on FRET data.\n"
"For conceptual foundations, see\n"
"  http://www.av8n.com/physics/correlation-norm.htm\n"
"\n"
"Typical usage:\n"
"  ./favia [options] > myfile.csv\n"
"\n"
"Long-form options include the following.\n"
"Defaults are shown in square brackets.\n"
"  --verbose            increase verbosity of debugging messages\n" 
"  --xfile=sss          first input file [data007.dat]\n"
"  --yfile=sss          second input file [clone of xfile]\n"
"  --jiffy=ddd          bin size (in s) for input timestamps [4e-12]\n"
"  --max_events=iii     pretend input file has at most this many events\n"
"                       [0 means use all of the input data]\n"
"  --fineness=iii       inverse grain size (in grains per octave) [8]\n"
"  --long_lag=ddd       longest lag (in seconds) to be considered [10]\n"
"  --short_lag=ddd      smallest lag (in seconds) to be considered [1e-8]\n"
"  --zerospike          include point at lag=0 in the output\n"
"  --help               print this message\n"
"\n"
"Note: sss=some string;  ddd=some double;  iii=some int.\n"
"Note: for each long-form option there exists the corresponding\n"
"short-form option, e.g. --jiffy=4e-12  <--> -j 4e-12\n"
"Note: -L is an easier-on-the-eyes synonym for -l\n"
"\n" 
"Input file formats:\n"
"  Binary: 64 bit integers;  each integer is a timestamp,\n"
"  representing an event, i.e. the detection of a photon\n"
"  in units of jiffies.\n"
"  Ascii: One ascii decimal number per line, representing\n"
"  a timestamp, as above.\n"
"Note: we auto-detect the input file format.\n"
"\n"
"Output file column headers are:\n"
"  lag, loglag, dot, dotnormed, bar, residual\n"
"where:\n"
"  lag       is measured in seconds\n"
"  loglag    is log10(lag)\n"
"  dot       is the raw dot product, i.e. the correlation at this lag\n"
"  dotnormed is the normalized dot product\n"
"  bar       is an error bar, i.e. a measure of the uncertainty\n"
"  residual  is measured relative to an estimate of the large-lag asymptote\n"
<< endl;

// Not implemented:
// "  PERIODIC=1   use periodic boundary conditions\n"
// "               [0 means pad with zeros]\n"

}

const char* Getenv(const char* key, const char* dflt){
  char* get = getenv(key);
  if (!get) return dflt;
  return get;
}


vector<daton> readfile_bin(const string fn,
        const int max__events){
  struct stat stats;
  int rslt = stat(fn.c_str(), &stats);
  if (rslt) {
    cerr << "Cannot stat data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  off_t size_b = stats.st_size;
  pdx_t tot__events = size_b / sizeof(cdp_t);
  if (max__events && tot__events > max__events) tot__events=max__events;
  cerr << "total_events: " << tot__events << endl;

  ifstream inch;
  inch.open(fn.c_str(), ios::in | ios::binary);
  if (!inch.good()) {
    cerr << "Favia cannot open binary data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  vector<daton> raw__data(tot__events);

// read raw data into memory:
  {                     
    cdp_t temp;
    pdx_t ii;
    for (ii=0; ii < tot__events; ii++) {
    if (!inch.good()) break;
      inch.read((char*)(&temp), sizeof(temp));
      if (inch.eof()) break;
#ifdef SLOW
      raw__data.push_back(daton(temp, 1));
#else
      raw__data[ii].abscissa = temp;
      raw__data[ii].ordinate = 1;
#endif
    }
    inch.close();
    if (ii != tot__events) {
      cerr << "Error reading file '" << fn << "' ..." << endl;
      cerr << "Bad count: " << ii << "  " << tot__events << endl;
      exit(1);
    }
  }
  return raw__data;
}

vector<daton> readfile_asc(const string fn,
        const int max__events){

  ifstream inch;
  inch.open(fn.c_str(), ios::in);
  if (!inch.good()) {
    cerr << "Favia cannot open ascii data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  vector<daton> raw__data;

// atoll returns a long long int;  make sure that is
// the right thing on this machine:
  assert (sizeof(long long int) == sizeof(cdp_t));

// read raw data into memory:
  {                     
    cdp_t temp;
    string buf;
    pdx_t ii;
    for (ii=0; ; ii++) {
    if (!inch.good()) break;
      getline(inch, buf);
      if (inch.eof()) break;
      temp = atoll(buf.c_str());        // see assertion above
// Alas push_back is slow, but we are obliged to use it
// (rather than indexing to element [ii]) because we 
// cannot determine the size of the array in advance.
// Besides, if you were interested in efficiency, you
// wouldn't be using ascii format anyway.
      raw__data.push_back(daton(temp, 1));
    }
    inch.close();
  }
  cout << "Input file: " << raw__data.size() 
       << " lines" << endl;
  return raw__data;
}

//////////
// Auto-detect the input file format,
// then call the appropriate handler.
vector<daton> readfile(const string fn, const int max__events) {

  ifstream inch;
  inch.open(fn.c_str(), ios::in | ios::binary);
  if (!inch.good()) {
    cerr << "Favia cannot open data file '" << fn << "': ";
    perror(0);
    exit(1);
  }

  cdp_t temp(-1);

  inch.read((char*)(&temp), sizeof(temp));
  inch.close();
// Ascii file should not have any nulls in it,
// certainly not a null in the 8th byte of the file.
// Conversely, the first timestamp in a binary file
// is restricted to be less than 2^56 units.
  if (temp >> 56) {
    return readfile_asc(fn, max__events);
  } else {
    return readfile_bin(fn, max__events);
  }
}

int main(int argc, char** argv) {

  int verbose(0);
  string xfn;
  string yfn;
  int fineness(8);
  pdx_t max_events(0);

  double long_lag(10);
  double short_lag(1e-8);
  double jiffy(4e-12);
  int spikeme(0);
  ostream* ouch(&cout);

// Process commandline options 
  const int ALT(128);
  static struct option long_options[] = {
    {"help",            0, NULL, 'h'},
    {"fineness",        1, NULL, 'f'},
    {"jiffy",           1, NULL, 'j'},
    {"long_lag",        1, NULL, 'l'},
// synonym, since lowercase "l" looks too much like digit "1":
    {"Long_lag",        1, NULL, 'L'},
    {"max_events",      1, NULL, 'm'},
    {"short_lag",       1, NULL, 's'},
    {"verbose",         0, NULL, 'v'},
    {"xfile",           1, NULL, 'x'},
    {"yfile",           1, NULL, 'y'},
    {"zerospike",       0, NULL, 'z'},
  };

  while(1) {
    extern char* optarg;
    char ch = getopt_long (argc, argv, long_options, NULL);
    if (ch == -1)
      break;
    
    if (optarg) if (*optarg == ':' || *optarg == '=') optarg++;
    switch(ch) {
      case 'h':
        usage(cout);
        exit(0);
      case 'f':
        fineness = atoi(optarg);
        break;
      case 'j':
        jiffy = atof(optarg);
        break;
      case 'l':
      case 'L':
        long_lag = atof(optarg);
        break;
      case 'm':
         max_events = atoi(optarg);
         break;
      case 's':
        short_lag = atof(optarg);
        break;
      case 'v':
        verbose++;
        break;
      case 'x':
        xfn = optarg;
        break;
      case 'y':
        yfn = optarg;
        break;
      case 'z':
        spikeme++;
        break;
      default:
        int chx(ch&~ALT);
	fprintf(stderr, "Sorry, option %s'%c' not yet implemented.\n", 
		ch!=chx? "ALT+" : "", chx);
	exit(1);
    }

  }
  if (xfn.length()==0) xfn = "data007.dat";
  if (yfn.length()==0) yfn = xfn;
  cerr << "xfn: " << xfn << endl;
  cerr << "yfn: " << yfn << endl;
  cerr << "jiffy: " << jiffy << endl;
  cerr << "fineness: " << fineness << endl;
  cerr << "max_events: " << max_events << endl;
  cerr << "short_lag: " << short_lag << endl;
  cerr << "long_lag: " << long_lag << endl;
  cerr << "verbose: " << verbose << endl;

// In this program, _mi_ means "... measured in units of ..."

// The "floor" here might make the smallest grain
// slightly smaller than what was requested by -s.
  double shortgrain_mi_bins = floor(short_lag / jiffy);

// Be a little bit defensive:
  if (shortgrain_mi_bins < 1) shortgrain_mi_bins = 1;

// Picture of the grains near the diagonal,
// as multiples of the shortest grain,
// assuming fineness=4:
// ||||||||| | | | |   |   |   |   |

// 000000000 1 1 1 1   2   2   2   3  (bin
// 012345678 0 2 4 6   0   4   8   2   number)
// ----0000111111112222222222222222  (octave number)

// --001111222222223333333333333333 ????

  double long_req_mi_short = ceil(long_lag / (shortgrain_mi_bins * jiffy));
  int octno(floor(log2(long_req_mi_short / fineness)));
  int octnox(octno > 0 ? octno : 0);

  cerr << "shortgrain_mi_bins: " << shortgrain_mi_bins 
        << "  long_req_mi_short: " << long_req_mi_short
        << "  log2: "  << log2(long_req_mi_short)
        << "  octno: "  << octno
        << "  octnox: "  << octnox
        << endl;

  double longgrain_mi_bins = shortgrain_mi_bins * (1 << octnox);
  cerr << "longgrain_mi_bins: " << longgrain_mi_bins
       << "  longgrain:  " << longgrain_mi_bins * jiffy
       << "  long/short: " << longgrain_mi_bins / shortgrain_mi_bins
        << endl;

//??     << "  ratio: " << double(longgrain_mi_bins*jiffy)/long_lag
//??     << "  ratio-1: " << double(longgrain_mi_bins*jiffy)/long_lag - 1.


class thingy{
public:
  vector<daton> raw_data;  
  pdx_t tot_events;
  vector<daton> cg_data;
  cdp_t front_event;
  cdp_t back_event;
  cdp_t aspan_mi_bins;

// constructor:
  thingy(const string fn, const int max__events,
        const double jiffy){
    raw_data = readfile(fn, max__events);
    tot_events = raw_data.size();

// Allocate space for coarse-grained data;
// For safety, overestimate the amount of space needed:
// the size of coarse-grained array can never exceed the
// size of the raw-data array.
    cg_data.resize(tot_events);
    front_event = raw_data[0].abscissa;
    back_event = raw_data[tot_events-1].abscissa;
/*****************************

Later, we will do a normalization calculation.  We will need to feed
it an unbiased estimate of the Poisson rate.  If we had a definite
number of points in a definite interval, this would be easy ... but we
don't have a definite interval.  All we have is the timestamp of the
last event.

Rather than talking about the rate, let's talk about inverse rate,
i.e. the time between pulses.  If we have N events, there are N+1
consecutive time intervals of interest:
 -- the leadin, i.e. the time before the front event;
 -- N-1 inter-event intervals; and
 -- the leadout.

Alas we do not know the leadout.  So we need a method that
(implicitly or explicitly) estimates the leadout, and then
estimates the rate.  So we are using the data to estimate
two things.  Our estimate of the rate will of course depend
on our estimate of the leadout.

As a first attempt, let's just use the N-1 inter-event intervals, then
the average interval is (t_back - t_front) / (N-1).  That is a
disaster when N=1, but we can improve it as follows.

As a second (and final) attempt, we make use of the leadin interval.
This is useful information, but it is not an unbiased estimate of the
inter-event interval, since the logger presumably turned on in the
middle of an interval, and this samples intervals in a biased way,
favoring longer intervals.  So let's count it as half of an
interval.  So our estimate of the average interval is then
(t_back - 0) / (N-1/2).  

Note that when N=1, this implicitly estimates the leadout as being
equal to the leadin.  More generally, this estimates the leadout as
t_back / (2N-1).  For large N this converges to half of the average
inter-event interval.  All this seems entirely reasonable.

***************************/

// span of actual data, 
// i.e. the sum of the N-1 inter-event intervals
// ... not used ...
//  cdp_t dspan_mi_bins = 1 + back_event - front_event;

// augmented span:
// includes the lead-in (before the first bin)
// and includes all the way to the /end/ of the last bin
// but does not include any lead-out (after last bin)
    aspan_mi_bins = 1 + back_event - 0;
    cerr  << "  front_event: " << front_event
      << "  back_event: "  << back_event 
      << "  aspan_mi_bins: "  << aspan_mi_bins
      << endl;
    cerr   << "  start-time of front_event bin: " << front_event*jiffy
           << "  start-time of back_event bin: "  << back_event*jiffy 
           << "  aug span time: "    << aspan_mi_bins*jiffy
           << endl;
    {     // spacing is just FYI;  not used elsewhere
     double spacing = aspan_mi_bins / tot_events;
     cerr << boost::format("approx avg spacing: "
        "%10.0f == %6.2e bins per event\n") 
           % spacing % spacing;
    }
  }
};      /* end class thingy */

  thingy x(xfn, max_events, jiffy);
  thingy y(xfn, max_events, jiffy);

// there are three requirements on the zone size:
// *) Zonesize must be a multiple of the largest grain size.
//    Not largest octave, but largest grain within that octave.
// *) Zonesize + longest lag <= X data size
// *) Zonesize <= Y data size
//
// The treatment of X and Y is not symmetric, because 
// only X gets lagged.

  double aspan_mi_longs = x.aspan_mi_bins / longgrain_mi_bins;
  const int reserved(1);
  double zone_mi_longs = floor(aspan_mi_longs) - reserved;
  double zone_mi_bins = zone_mi_longs * longgrain_mi_bins;

  cerr << "aspan_mi_longs: " << aspan_mi_longs
       << "  zone_mi_longs: " << zone_mi_longs
       << "  zone_mi_bins: " << zone_mi_bins
       << endl;

// There are fineness cells in an ordinary octave, but
// the startup requirement requires twice that many.
// Also, the typical octave uses curlag_mi_grains values in 
// the range [fineness, 2*fineness-1] ... so we will 
// have many uses for fin2:
  pdx_t fin2(2*fineness);

// curlag_mi_grains is the #1 crucial loop variable in all that follows.
// It is manipulated by the outer loop *and* by the inner loop.
// It is an integer, representing a number of grains.
// It is a number in the range [0, fin2]
// ... usually in the range [fineness, fin2-1]
// Note that curlag_mi_grains * grain_mi_bins * jiffy == 
//  start time (in seconds) of the grain we are working on.
  cdp_t curlag_mi_grains(1);
  if (spikeme) curlag_mi_grains=0;         // output will include zero spike

  dot_t zerospike(0);
  dot_t hits(0);

// curgrain_mi_bins is the #2 crucial loop variable.
// It changes in the outer loop, 
// but is constant with respect to the inner loop.

// Meanwhile, shortgrain_mi_bins never changes;  
// it is the size of the shortest grain.
// Also, the bin size never changes; 
// a bin is always one jiffy long.
  cdp_t curgrain_mi_bins = shortgrain_mi_bins;

// main loop over all resolutions
  for (pdx_t kk=0; ; kk++) {

// here with the following loop variables already set up
//   curlag_mi_grains = ...
//   curgrain_mi_bins = ...

    double bogus_zone_mi_grains = double(x.aspan_mi_bins) / curgrain_mi_bins;

    double norm_denom = double(x.tot_events)*double(y.tot_events)
                / bogus_zone_mi_grains;

    double curgrain_mi_sec(jiffy * curgrain_mi_bins);
    double lag = curlag_mi_grains * curgrain_mi_sec;
    cerr << "top of loop:  lag_mi_grains: " << curlag_mi_grains
      << "  curgrain_mi_bins: " << curgrain_mi_bins
      << "  iniial lag: " << lag << " sec"
      << "  denom: " << norm_denom
      << endl;

    if (lag > long_lag) break;

// remap so data starts at zero, with coarse graining:

    pdx_t x_small = tighten(x.front_event, curgrain_mi_bins, 
                &x.raw_data[0], x.raw_data.size(), 
                &x.cg_data[0]);
    pakvec pv1(0, 1+x.back_event/curgrain_mi_bins, &x.cg_data[0], x_small);

    pdx_t y_small = tighten(y.front_event, curgrain_mi_bins, 
                &y.raw_data[0], y.raw_data.size(), 
                &y.cg_data[0]);
    pakvec pv2(0, 1+y.back_event/curgrain_mi_bins, &y.cg_data[0], y_small);

// inner loop over all lags, stepping grain by grain:
    for ( ; curlag_mi_grains < fin2 ; curlag_mi_grains++) {
      lag = curlag_mi_grains * curgrain_mi_sec;
      if (lag > long_lag) goto main_done;

      if (0) cerr << "curgrain_mi_bins: " << curgrain_mi_bins
        << "  curlag_mi_grains: " << curlag_mi_grains
        << "  lag: " << curlag_mi_grains*curgrain_mi_sec << endl;
      pv1.set_bin0time(curlag_mi_grains);
      dot_t dot = pakdot(pv1, pv2);
      if (curlag_mi_grains == 0) zerospike = dot;
      hits += dot;
      double loglag(-9999);
      if (lag) loglag = log10(lag);
      
      double dotnormed = double(dot) / norm_denom;
      double bar(dot ? dotnormed / sqrt(dot) : 0);
      double model = 1 - double(curlag_mi_grains) / bogus_zone_mi_grains;
      double residual = dotnormed - model;
// If you change this output statement, be sure to
// change the usage() message to match:
      *ouch << boost::format
        ("%14.8e, %14.8e, %10Ld, %14.8e, %14.8e, %14.8e,\n")
        % lag % loglag % dot % dotnormed % bar % residual;
      if (dot) if (verbose)  cerr << boost::format
           ("curgrain_mi_bins:%12Ld curlag_mi_grains: %12Ld"
              "  lag: %14.8e  dot: %10Ld  dot/x: %14.8e\n")
                % curgrain_mi_bins % curlag_mi_grains 
                % lag % dot % (double(dot)/double(norm_denom));
    }

// prepare for the next iteration:
    curlag_mi_grains /= 2;
    curgrain_mi_bins *= 2;
  }

main_done:;;;
  
  if (zerospike) {
    dot_t other = hits - zerospike;
    cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
       % zerospike % other % (zerospike + 2*other);
  }

  return 0;
}
