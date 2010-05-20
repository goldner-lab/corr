// Driver program for testing sparse_dot.c
//
// Typical usage: VERBOSE=1 PERIODIC=1 WIDTH=100 ./paktest

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

const char* Getenv(const char* key, const char* dflt){
  char* get = getenv(key);
  if (!get) return dflt;
  return get;
}

int main(int argc, char** argv) {

  double blocksize = atof(Getenv("BLOCKSIZE",  "8"));
  double max_time = atof(Getenv("MAX_TIME",  "10"));
  double min_time  = atof(Getenv("MIN_TIME",  "1e-8"));
  pdx_t max_events = atoi(Getenv("MAX_EVENTS", "0"));
//  int width =    atoi(Getenv("WIDTH",    "100"));
  int verbose =  atoi(Getenv("VERBOSE",  "0"));

  double raw_dt(4e-12);
  cdp_t grain = cdp_t(min_time / raw_dt + .5);

  string fn("data007.dat");
  if (argc > 1) fn = argv[1];
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

  {
    pdx_t ii;
    for (ii=0; ii < tot_events; ii++) {
    if (!inch.good()) break;
      inch.read((char*)(&temp), sizeof(temp));
      if (inch.eof()) break;
      raw_data[ii].abscissa = temp;
      raw_data[ii].ordinate = 1;
    }
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
    double time = shift * dt;
    cerr << "top of loop:  shift: " << shift
      << "  grain: " << grain
      << "  time: " << time
      << "  denom: " << norm_denom
      << endl;

    if (time > max_time) break;

// remap so data starts at zero, with coarse graining:

    pdx_t small = tighten(first_event, grain, 
                &raw_data[0], raw_data.size(), 
                &cg_data[0]);
    pakvec pv1(0, 1+last_event, &cg_data[0], small);
    pakvec pv2(0, 1+last_event, &cg_data[0], small);

// inner loop over all shifts at this resolution:
    for ( ; shift < blksiz2 ; shift++) {
      time = shift * dt;
      if (0) cerr << "grain: " << grain
        << "  shift: " << shift
        << "  time: " << shift*dt << endl;
      pv1.earlyize(shift);
      dot_t dot = pakdot(pv1, pv2);
      if (shift == 0) zerospike = dot;
      hits += dot;
      double logtime(-9999);
      if (time) logtime = log10(time);
      
      double dotnormed = double(dot) / norm_denom;
      double bar(0);
      if (dot) bar = dotnormed / sqrt(dot);
      double model = 1 - double(shift) / coarse_bins;
      double residual = dotnormed - model;
      cout << boost::format
        ("%14.8e, %14.8e, %10Ld, %14.8e, %14.8e, %14.8e,\n")
        % time % logtime % dot % dotnormed % bar % residual;
      if (dot) if (verbose)  cerr << boost::format
        ("grain:%12Ld shift: %12Ld  time: %14.8e  dot: %10Ld"
                "  dot/x: %14.8e\n")
        % grain % shift % time % dot % (double(dot)/double(norm_denom));
    }
  }
  
  if (zerospike) {
    dot_t other = hits - zerospike;
    cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
       % zerospike % other % (zerospike + 2*other);
  }

  return 0;
}
