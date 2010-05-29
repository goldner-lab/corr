// For additional information on strategy, terminology,
// and representation issues, see
//   http://www.av8n.com/physics/correlation-norm.htm
//
// Basic notions:
// The conceptual vector (cv) is a histogram with
// cvnp bins, numbered 0 through cvnp-1.

// The conceptual vector is binary and very sparse, i.e.
// the overwhelming majority of bins contain zero events;
// only a very few of the bins contain one or more events.
// In the leading application, each event is the detection
// of a photon.
//
// A pakvec is a packed representation of such a vector.
// We can efficiently compute the dot product of two
// pakvecs, without ever unpacking them.
// 
// A pakvec is constructed from a list of event times.
// (absolute times, not inter-event intervals)
// plus some additional information.
// The additional information is needed because we 
// cannot tell from the list of event-times how many
// zeros there are in the conceptual vector before
// the first event or after the last event.
// Such details may or may not be important.
// They are definitely important if you are doing
// convolutions with periodic boundary conditions.
//
// Also, even for non-periodic boundary conditions,
// it is sometimes convenient to make the conceptual
// vector rather short, so that some of the data at
// the beginning and/or end of the event-list is
// ignored i.e. not needed to populate the cv.
//
// Furthermore, the pakvec maintains some additional
// derived information, behind your back, to increase 
// efficiency.
//
// Useful feature:  You can shift the start-time of
// the conceptual vector (i.e. the time corresponding
// to bin 0 of the cv) after the pakvec is created,
// by calling set_bin0time(new0time).  This is efficient
// if successive calls have nondecreasing new0time values.
// It is also efficient if the new0time value corresponds 
// to an early time, not beyond the first few events in
// the packed vector.
// Note the orientation(s) here: you would add cv0time to
// the cv bin-time to find the timestamp of an event in 
// the packed vector, whereas conversely you would subtract
// cv0time from the timestamp of an event in the packed
// vector to find the cv bin-time where the event belongs.
//
// Never mess with the value of cv0time directly;
// always call set_cv0time() instead.
// I could make cv0time private, and force everybody to use an
// accessor method, but this would incur a performance penalty.
//
// In contrast, it is OK to change cvnp by direct assignment.

//!! Restriction:  Beware:  Do not do anything that might
//!! disturb the packed data (i.e. the list of event-times)
//!! after the moment when a pointer to said data is
//!! passed to the constructor for the pakvec,
//!! unless/until you are through using the pakvec.
//!! For example, if the list was in an I/O buffer, don't
//!! re-use that buffer.  When in doubt, make a private
//!! copy of the list of event-times and pass the private
//!! copy to the constructor.
//!! 
//!! Using the stdlib vector<daton> has advantages and
//!! disadvantages.  Beware that the data is likely to
//!! get copied to a new place if you resize() the 
//!! vector to a larger size;  this invalidates any
//!! old pointers, including the pointers stored within 
//!! the pakvec.

// The dot product routine assumes the pakvec has
// _strictly_ increasing abscissas.  If in doubt,
// use the tighten() routine to convert nondecreasing
// abscissas to strictly increasing abscissas.

#include <stdint.h>

// cdp_t should be able to represent very large
// numbers, corresponding to the number of bins
// in the conceptual vector.  
// (cdp stands for conceptual data point)
typedef int64_t cdp_t;

// pdx_t should be able to represent moderately large
// numbers, corresponding to the number of events
// in the record.  This number serves as an upper bound
// on the number of points in the packed vector, and
// also as an upper bound on the ordinate of any
// particular daton.
// (pdx stands for packed index)
// Presumably 32-bit ints suffice for this, unless 
// you've got many gigabytes of virtual memory, and a
// compiler than implements an array index operator []
// than can handle huge indexes.
typedef int32_t pdx_t;

// dot_t should be able to represent the _square_ of
// the number of events in the record
typedef int64_t dot_t;

class daton{
public:
  cdp_t abscissa;
  pdx_t ordinate;
// trivial default constructor:
  daton(){
    abscissa = 0;
    ordinate = 0;
  }
// basic constructor:
  daton(const cdp_t abs, const pdx_t ord) {
    abscissa = abs;
    ordinate = ord;
  }
};

class pakvec{
public:
  cdp_t cv0time;	// time (in jiffies) corresponding to
			// bin 0 of the conceptual vector
  cdp_t cvnp;		// number of bins in the cv
  const daton* spdata;	// pointer to static packed data
  pdx_t spnp;		// number of points in static packed data
  const daton* updata;	// pointer to useful packed data
  pdx_t	upnp;		// number of useful points 
  const char* verbose;

inline cdp_t get_cv0time() const {
  return cv0time;
}

// directly using public member variable in innermost loop:
// real    0m13.838s
// user    0m13.545s
// sys     0m0.124s
//
// using inline accessor method, no optimization:
// real    0m18.296s
// user    0m18.201s
// sys     0m0.056s


// Returns the number of items that can be skipped
// when we are not using periodic boundary conditions.
// Some points are skippable because packed data items that 
// correspond to negative bin numbers don't matter.
inline pdx_t skippable(){
  for (pdx_t p1 = 0; p1 < upnp; p1++){
    if (updata[p1].abscissa >= cv0time) return p1;
  }
  return upnp;		// didn't find anything in range
}


inline void set_bin0time(const cdp_t new0time){
// For efficiency, new0time should not be less than
// the old cv0time.
// To say the same thing the other way,
// calling us with new0time less than cv0time
// is grudgingly permitted but is inefficient, because
// it forces us -- via skippable() -- to recalculate 
// upnp and updata _ab initio_.
// This recalculation is inefficient but necessary,
// assuming you want to get the right answer.
  if (new0time < cv0time) {	
    upnp = spnp;        // inefficient: reset to safe value
    updata = spdata;    // ditto
    cerr << "inefficient retrograde set_bin0time" << endl;
  }
  cv0time = new0time;
  pdx_t nnn = skippable();
  updata += nnn;  upnp -= nnn;
}

// Returns zero if ndx is within the valid range of the data,
// nonzero otherwise.  
//   In more detail:
//   Returns 2 if ran out of data without reaching end of window.
//   Returns 1 if indicated data is outside the window.
//   Returns 0 if within window.
//
// Note: when doing convolutions, you probably want to
// arrange (or at least check) that neither vector ever
// returns 2, or neither vector ever returns 1.
// Otherwise, you are using a non-constant fraction of
// the packed data to populate the conceptual vector.
//
// Specifically:  If you have 2 hours worth of raw data,
// and your conceptual vector is 1 minute long, the maximum
// lag in your autocorrelation should be 119 minutes, not 120.
inline pdx_t outwin(const pdx_t ndx) const{
  if (ndx >= upnp) {
      if (verbose) 
        cout << "\t" << verbose << " outwin returns 2" << endl;
      return 2;
  }
  if (updata[ndx].abscissa - cv0time >= cvnp) {
    if (verbose) 
      cout << "\t" << verbose << " outwin returns 1: " 
      	<< updata[ndx].abscissa << " " << cvnp << endl;
    return 1;
  }
  return 0;
}

// Constructor.  Used to create a new pakvec.
inline pakvec(const cdp_t cv0time_, const cdp_t cvnp_, 
	const daton* spdata_, pdx_t size_
){
  cv0time = cv0time_;
  cvnp = cvnp_;
  spdata = spdata_;
  updata = spdata_;
  upnp = spnp = size_;
  set_bin0time(cv0time);         // do skipping if needed.
  verbose = 0;
}

inline pdx_t total(){
  pdx_t it;
  pdx_t rslt(0);
  for (it = 0; it < spnp; it++) {
    rslt += spdata[it].ordinate;
  }
  return rslt;
}

};	// end of class pakvec

////////////////////////////
//
// Interface definitions for routines provided by pakdot.c

// Take the dot product of two packed vectors.
dot_t pakdot(const pakvec pv1, const pakvec pv2, cdp_t & gs1, cdp_t & gs2);
dot_t pakdot(const pakvec pv1, const pakvec pv2);
dot_t pbc_pakdot(const pakvec pv1, const pakvec pv2);

pdx_t tighten(
  const cdp_t first, const cdp_t grain,
  const daton* loose, const pdx_t size, daton* tight
);
