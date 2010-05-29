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
"  --jiffy=ddd          unit of time (in s) for input timestamps [4e-12]\n"
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
    cerr << "Cannot open data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  cdp_t temp;
  string buf;
  vector<daton> raw__data;

// atoll returns a long long int;  make sure that is
// the right thing on this machine:
  assert (sizeof(long long int) == sizeof(cdp_t));

// read raw data into memory:
  {                     
    pdx_t ii;
    for (ii=0; ; ii++) {
    if (!inch.good()) break;
      getline(inch, buf);
      if (inch.eof()) break;
      temp = atoll(buf.c_str());
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
    cerr << "Favia auto-detect cannot open data file '" << fn << "': ";
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


  vector<daton> x_raw_data = readfile(xfn, max_events);
  pdx_t x_tot_events = x_raw_data.size();
// Allocate array for coarse-grained data, leaving plenty of room:
  vector<daton> x_cg_data(x_tot_events);
  cdp_t x_first_event = x_raw_data[0].abscissa;
  cdp_t x_last_event = x_raw_data[x_tot_events-1].abscissa;

  vector<daton> y_raw_data = readfile(yfn, max_events);
  pdx_t y_tot_events = y_raw_data.size();
  vector<daton> y_cg_data(y_tot_events);
  cdp_t y_first_event = y_raw_data[0].abscissa;
  cdp_t y_last_event = y_raw_data[y_tot_events-1].abscissa;

// span of actual data:  not used
// cdp_t dspan_mi_bins = 1 + x_last_event - x_first_event;

// augmented span:
// includes the lead-in (time before the first event)
// but does not include any lead-out (time after last event).
  cdp_t aspan_mi_bins = 1 + x_last_event - 0;

// there are three requirements on the zone size:
// *) Zonesize must be a multiple of the largest grain size.
//    Not largest octave, but largest grain within that octave.
// *) Zonesize + longest lag <= X data size
// *) Zonesize <= Y data size
//
// The treatment of X and Y is not symmetric, because 
// only X gets lagged.

  double aspan_mi_longs = aspan_mi_bins / longgrain_mi_bins;
  const int reserved(1);
  double zone_mi_longs = floor(aspan_mi_longs) - reserved;
  double zone_mi_bins = zone_mi_longs * longgrain_mi_bins;
  cerr << "aspan_mi_longs: " << aspan_mi_longs
       << "  zone_mi_longs: " << zone_mi_longs
       << "  zone_mi_bins: " << zone_mi_bins
       << endl;

  cerr  << "  x_first_event: " << x_first_event
    << "  x_last_event: "  << x_last_event 
    << "  aspan_mi_bins: "  << aspan_mi_bins
    << endl;

  cerr   << "  start-time of x_first_event bin: " << x_first_event*jiffy
         << "  start-time of x_last_event bin: "  << x_last_event*jiffy 
         << "  aug span time: "    << aspan_mi_bins*jiffy
         << endl;

  {     // spacing is just FYI;  not used elsewhere
   double spacing = aspan_mi_bins / x_tot_events;
   cerr << boost::format("approx avg spacing: "
      "%10.0f == %6.2e bins per event\n") 
         % spacing % spacing;
  }

// There are fineness cells in an ordinary octave, but
// the startup requirement requires twice that many.
// Also, the typical octave uses shifts in the range
// [fineness, 2*fineness-1] ... so we will have many
// uses for fin2:
  pdx_t fin2(2*fineness);

// 'shift' is a number in the range [0, fin2]
// ... usually in the range [fineness, fin2-1]
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

    double coarse_bins = double(aspan_mi_bins) / double(grainsize);
    double norm_denom = double(x_tot_events)*double(x_tot_events)
                / coarse_bins;

    double dt(jiffy * grainsize);
    double lag = shift * dt;
    cerr << "top of loop:  shift: " << shift
      << "  grainsize: " << grainsize
      << "  iniial lag: " << lag
      << "  denom: " << norm_denom
      << endl;

    if (lag > long_lag) break;

// remap so data starts at zero, with coarse graining:

    pdx_t x_small = tighten(x_first_event, grainsize, 
                &x_raw_data[0], x_raw_data.size(), 
                &x_cg_data[0]);
    pakvec pv1(0, 1+x_last_event/grainsize, &x_cg_data[0], x_small);

    pdx_t y_small = tighten(y_first_event, grainsize, 
                &y_raw_data[0], y_raw_data.size(), 
                &y_cg_data[0]);
    pakvec pv2(0, 1+y_last_event/grainsize, &y_cg_data[0], y_small);

// inner loop over all shifts at this resolution:
    for ( ; shift < fin2 ; shift++) {
      lag = shift * dt;
      if (lag > long_lag) goto main_done;

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
      *ouch << boost::format
        ("%14.8e, %14.8e, %10Ld, %14.8e, %14.8e, %14.8e,\n")
        % lag % loglag % dot % dotnormed % bar % residual;
      if (dot) if (verbose)  cerr << boost::format
        ("grainsize:%12Ld shift: %12Ld  lag: %14.8e  dot: %10Ld"
                "  dot/x: %14.8e\n")
        % grainsize % shift % lag % dot % (double(dot)/double(norm_denom));
    }
  }

main_done:;;;
  
  if (zerospike) {
    dot_t other = hits - zerospike;
    cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
       % zerospike % other % (zerospike + 2*other);
  }

  return 0;
}
