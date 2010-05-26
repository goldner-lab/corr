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

void usage() {
  cout << 
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
"  --jiffy=ddd          unit of time (in s) for input timestamps [4e-12]\n"
"  --max_events=iii     pretend input file has at most this many events\n"
"                       [0 means use all of the input data]\n"
"  --blocksize=iii      number of grains per octave [8]\n"
"  --long_lag=ddd       longest lag (in seconds) to be considered [10]\n"
"  --short_lag=ddd      smallest lag (in seconds) to be considered [1e-8]\n"
"  --zerospike          include point at lag=0 in the output\n"
"  --help               print this message\n"
"\n"
"Note: sss=some string;  ddd=some double;  iii=some int.\n"
"Note: for each long-form option there exists the corresponding\n"
"short-form option, e.g. -j 4e-12\n"
"\n" 
"Input file format: 64 bit integers;  each integer is\n"
"a timestamp, representing an event, i.e. the detection of a photon\n"
"in units of jiffies.\n"
"\n"
"Output file column headers are:\n"
"  lag, loglag, dot, dotnormed, bar, residual\n"
"where:\n"
"  lag      is measured in seconds\n"
"  loglag   is log10(lag)\n"
"  dot      is the dot product, i.e. the correlation at this lag\n"
"  bar      is an error bar, i.e. a measure of the uncertainty\n"
"  residual is measured relative to an estimate of the large-lag asymptote\n"
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


vector<daton> readfile(const string fn, const int max__events) {
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
    cerr << "Cannot open data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  cdp_t temp;
  vector<daton> raw__data(tot__events);

// read raw data into memory:
  {                     
    pdx_t ii;
    for (ii=0; ii < tot__events; ii++) {
    if (!inch.good()) break;
      inch.read((char*)(&temp), sizeof(temp));
      if (inch.eof()) break;
      raw__data[ii].abscissa = temp;
      raw__data[ii].ordinate = 1;
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

int main(int argc, char** argv) {

  int verbose(0);
  string xfn;
  string yfn;
  int blocksize(8);
  pdx_t max_events(0);

  double long_lag(10);
  double short_lag(1e-8);
  double jiffy(4e-12);
  int spikeme(0);

// Process commandline options 
  const int ALT(128);
  static struct option long_options[] = {
    {"help",            0, NULL, 'h'},
    {"blocksize",       1, NULL, 'b'},
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
        usage();
        exit(0);
      case 'b':
        blocksize = atoi(optarg);
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
  cerr << "blocksize: " << blocksize << endl;
  cerr << "max_events: " << max_events << endl;
  cerr << "short_lag: " << short_lag << endl;
  cerr << "long_lag: " << long_lag << endl;
  cerr << "verbose: " << verbose << endl;

// The "floor" here might make the smallest grain
// slightly smaller than what was requested by -s.
  double bins_per_short = floor(short_lag / jiffy);
  double bins_per_long_req = ceil(long_lag / jiffy);
  int doublings = ceil(log2(bins_per_long_req / bins_per_short));
  cerr << "bins_per_short: " << bins_per_short 
        << "  bins_per_long_req: " << bins_per_long_req
        << "  ratio: " << bins_per_long_req / bins_per_short
        << "  log2: "  << log2(bins_per_long_req / bins_per_short)
        << "  doublings: "  << doublings
        << endl;

  double bins_per_long = bins_per_short * (1 << doublings);
  cerr << "bins_per_long: " << bins_per_long
        << "  ratio: " << double(bins_per_long)/bins_per_long_req
        << "  ratio-1: " << double(bins_per_long)/bins_per_long_req - 1.
        << endl;

  vector<daton> x_raw_data = readfile(xfn, max_events);
  pdx_t x_tot_events = x_raw_data.size();
// Allocate array for coarse-grained data, leaving plenty of room:
  vector<daton> cg_data(x_tot_events);

  cdp_t first_event = x_raw_data[0].abscissa;
  cdp_t last_event = x_raw_data[x_tot_events-1].abscissa;

// span of actual data:  not used
// cdp_t dspan_bins = 1 + last_event - first_event;

// augmented span:
// includes the lead-in (time before the first event)
// but does not include any lead-out (time after last event).
  cdp_t aspan_bins = 1 + last_event - 0;

  double longs_per_aspan = aspan_bins / bins_per_long;
  const int reserved(1);
  double longs_per_zone = floor(longs_per_aspan) - reserved;
  double bins_per_zone = longs_per_zone * bins_per_long;
  cerr << "longs_per_aspan: " << longs_per_aspan
       << "  longs_per_zone: " << longs_per_zone
       << "  bins_per_zone: " << bins_per_zone
       << endl;

  cerr  << "  first_event: " << first_event
    << "  last_event: "  << last_event 
    << "  aspan_bins: "  << aspan_bins
    << endl;

  cerr   << "  start-time of first_event bin: " << first_event*jiffy
         << "  start-time of last_event bin: "  << last_event*jiffy 
         << "  aug span time: "    << aspan_bins*jiffy
         << endl;

  {     // spacing is just FYI;  not used elsewhere
   double spacing = aspan_bins / x_tot_events;
   cerr << boost::format("approx avg # of bins per event: %10.0f == %6.2e\n") 
         % spacing % spacing;
  }

// There are blocksize cells in an ordinary octave, but
// the startup requirement requires twice that many.
// Also, the typical octave uses shifts in the range
// [blocksize, 2*blocksize-1] ... so we will have many
// uses for blksiz2:
  pdx_t blksiz2(2*blocksize);

// 'shift' is a number in the range [0, blksiz2]
// ... usually in the range [blocksize, blksiz2-1]
// A shift is essentially a coarse-grained version
// of a lag.
// Note that shift * grainsize * jiffy == 
//  start time (in seconds) of the grain at this shift
  cdp_t shift(1);
  if (spikeme) shift=0;         // output will include zero spike

  dot_t zerospike(0);
  dot_t hits(0);
  cdp_t grainsize = cdp_t(short_lag / jiffy + .5);

// main loop over all resolutions
  for (pdx_t kk=0; ; kk++) {
    if (kk) {
      shift /= 2;
      grainsize *= 2;
    }

    double coarse_bins = double(aspan_bins) / double(grainsize);
    double norm_denom = double(x_tot_events)*double(x_tot_events)
                / coarse_bins;

    double dt(jiffy * grainsize);
    double lag = shift * dt;
    cerr << "top of loop:  shift: " << shift
      << "  grainsize: " << grainsize
      << "  lag: " << lag
      << "  denom: " << norm_denom
      << endl;

    if (lag > long_lag) break;

// remap so data starts at zero, with coarse graining:

    pdx_t small = tighten(first_event, grainsize, 
                &x_raw_data[0], x_raw_data.size(), 
                &cg_data[0]);
    pakvec pv1(0, 1+last_event/grainsize, &cg_data[0], small);
    pakvec pv2(0, 1+last_event/grainsize, &cg_data[0], small);

// inner loop over all shifts at this resolution:
    for ( ; shift < blksiz2 ; shift++) {
      lag = shift * dt;
      if (0) cerr << "grainsize: " << grainsize
        << "  shift: " << shift
        << "  lag: " << shift*dt << endl;
      pv1.set_bin0time(shift);
      dot_t dot = pakdot(pv1, pv2);
      if (shift == 0) zerospike = dot;
      hits += dot;
      double loglag(-9999);
      if (lag) loglag = log10(lag);
      
      double dotnormed = double(dot) / norm_denom;
      double bar(0);
      if (dot) bar = dotnormed / sqrt(dot);
      double model = 1 - double(shift) / coarse_bins;
      double residual = dotnormed - model;
// If you change this output statement, be sure to
// change the usage() message to match:
      cout << boost::format
        ("%14.8e, %14.8e, %10Ld, %14.8e, %14.8e, %14.8e,\n")
        % lag % loglag % dot % dotnormed % bar % residual;
      if (dot) if (verbose)  cerr << boost::format
        ("grainsize:%12Ld shift: %12Ld  lag: %14.8e  dot: %10Ld"
                "  dot/x: %14.8e\n")
        % grainsize % shift % lag % dot % (double(dot)/double(norm_denom));
    }
  }
  
  if (zerospike) {
    dot_t other = hits - zerospike;
    cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
       % zerospike % other % (zerospike + 2*other);
  }

  return 0;
}
