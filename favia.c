// Perform correlations on FRET data.
// See "usage" message below, or try ./favia -h
//

using namespace std;

#include <sys/types.h>          /* for stat */
#include <sys/stat.h>           /* for stat */
#include <unistd.h>             /* for stat */
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
"  ./favia myfile.dat > myfile.csv\n"
"\n"
"Useful environment variables include the following.\n"
"Defaults are shown in square brackets.\n"
"  VERBOSE=1    turn on informational and debugging messages [0]\n" 
"  MAX_EVENTS=1000  pretend input file has at most this many events\n"
"               [0 means use all of the input file]\n"
"  BLOCKSIZE=16 each new block is this many times bigger than previous [8]\n"
"  MAX_LAG=??  largest lag (in seconds) to be considered [10]\n"
"  MIN_LAG=??  smallest lag (in seconds) to be considered [1e-8]\n"
"\n"
"Input file format: 64 bit integers;  each integer is\n"
"a timestamp, representing an event, i.e. the detection of a photon\n"
"in units of 4 picoseconds\n"
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

int main(int argc, char** argv) {

  double blocksize = atof(Getenv("BLOCKSIZE",  "8"));
  double max_lag = atof(Getenv("MAX_LAG",  "10"));
  double min_lag  = atof(Getenv("MIN_LAG",  "1e-8"));
  pdx_t max_events = atoi(Getenv("MAX_EVENTS", "0"));
  int verbose =  atoi(Getenv("VERBOSE",  "0"));

// If you change the following, be sure to change
// the usage() message accordingly:
  const double raw_dt(4e-12);
  cdp_t grain = cdp_t(min_lag / raw_dt + .5);

  string fn("data007.dat");
  if (argc > 1) {
    string arg = argv[1];
    if (arg == "-h" || arg == "--h" || arg == "--help") {
      usage();
      exit(0);
    }
    fn = arg;
  }
  struct stat stats;
  int rslt = stat(fn.c_str(), &stats);
  if (rslt) {
    cerr << "Cannot stat data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  off_t size_b = stats.st_size;
  pdx_t tot_events = size_b / sizeof(cdp_t);
  if (max_events && tot_events > max_events) tot_events=max_events;
  cerr << "total_events: " << tot_events << endl;

  ifstream inch;
  inch.open(fn.c_str(), ios::in | ios::binary);
  if (!inch.good()) {
    cerr << "Cannot open data file '" << fn << "': ";
    perror(0);
    exit(1);
  }
  cdp_t temp;
  vector<daton> raw_data(tot_events);
  vector<daton> cg_data(tot_events);          // coarse-grained data

// read raw data into memory:
  {                     
    pdx_t ii;
    for (ii=0; ii < tot_events; ii++) {
    if (!inch.good()) break;
      inch.read((char*)(&temp), sizeof(temp));
      if (inch.eof()) break;
      raw_data[ii].abscissa = temp;
      raw_data[ii].ordinate = 1;
    }
    inch.close();
    if (ii != tot_events) {
      cerr << "Bad count: " << ii << "  " << tot_events << endl;
      exit(1);
    }
  }

  pdx_t blksiz2(2*blocksize);
  pdx_t zerospike(0);
  dot_t hits(0);

// 'shift' is a number in the range [0, blksiz2]
// a shift is essentially a coarse-grained version
// of a lag.
  cdp_t shift(1);
// note that shift * grain * raw_dt == realtime
  
  cdp_t first_event = raw_data[0].abscissa;
  cdp_t last_event = raw_data[tot_events-1].abscissa;
  cdp_t lead_out = first_event;     // heuristic construction
  cdp_t span_bins = 1 + last_event + lead_out;

// Here we have a conceptual vector with
// 9 bins labeled 0 through 8 inclusive
//   +----+----+----+----+----+----+----+----+----+
//   |0   |1   |2   |3 * |4   |5 * |6   |7   |8   |9
//   +----+----+----+----+----+----+----+----+----+
//   [-- lead_in  --]              [-- lead_out --]
//
//  record_first   = 0
//  record_last    = 8 = (N-1) = not very interesting
//  N == span_bins = 9 = 1 + last_event + lead_out
//  span_time      = N*dt
//  first_event    = 3
//  last_event     = 5
//  lead_in        = 3
//  lead_out       = 3 by heuristic construction;  not recorded by instrument

  cerr  << "  first_event: " << first_event
    << "  last_event: "  << last_event 
    << "  span_bins: "  << span_bins
    << endl;

  cerr   << "  start-time of first_event bin: " << first_event*raw_dt
         << "  start-time of last_event bin: "  << last_event*raw_dt 
         << "  span time: "    << span_bins*raw_dt
         << endl;

  double spacing = span_bins / tot_events;
  cerr << boost::format("average # of bins per event: %10.0f == %6.2e\n") 
        % spacing % spacing;

// main loop over all resolutions
  for (pdx_t kk=0; ; kk++) {
    if (kk) {
      shift /= 2;
      grain *= 2;
    }

    double coarse_bins = double(span_bins) / double(grain);
    double norm_denom = double(tot_events)*double(tot_events)
                / coarse_bins;

    double dt(raw_dt * grain);
    double lag = shift * dt;
    cerr << "top of loop:  shift: " << shift
      << "  grain: " << grain
      << "  lag: " << lag
      << "  denom: " << norm_denom
      << endl;

    if (lag > max_lag) break;

// remap so data starts at zero, with coarse graining:

    pdx_t small = tighten(first_event, grain, 
                &raw_data[0], raw_data.size(), 
                &cg_data[0]);
    pakvec pv1(0, 1+last_event, &cg_data[0], small);
    pakvec pv2(0, 1+last_event, &cg_data[0], small);

// inner loop over all shifts at this resolution:
    for ( ; shift < blksiz2 ; shift++) {
      lag = shift * dt;
      if (0) cerr << "grain: " << grain
        << "  shift: " << shift
        << "  lag: " << shift*dt << endl;
      pv1.earlyize(shift);
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
        ("grain:%12Ld shift: %12Ld  lag: %14.8e  dot: %10Ld"
                "  dot/x: %14.8e\n")
        % grain % shift % lag % dot % (double(dot)/double(norm_denom));
    }
  }
  
  if (zerospike) {
    dot_t other = hits - zerospike;
    cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
       % zerospike % other % (zerospike + 2*other);
  }

  return 0;
}
