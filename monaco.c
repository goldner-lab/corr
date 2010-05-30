// Generate simple model data for FRET
// See "usage" message below, or try ./monaco -h
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
#include <boost/random.hpp>

#include <assert.h>

void usage(ostream& foo) {
  foo << 
"Generate simple model data for FRET\n"
"\n"
"Typical usage:\n"
"  ./monaco [options] > myfile.dat\n"
"\n"
"Long-form options include the following.\n"
"Defaults are shown in square brackets.\n"
"  --verbose            increase verbosity of debugging messages\n" 
"  --maxdraw=iii        number of random draws [1000000]\n"
"  --ofile=sss          output file [cal.dat]\n"
"  --poisson=ddd        poisson rate (in probability per bin) [.01]\n"
"  --scatter            output ascii decimal csv, two per line\n"
"  --help               print this message and exit\n"
"\n" 
"Output file format:\n"
"  Binary: 64 bit integers;  each integer is a timestamp,\n"
"  representing an event, i.e. the detection of a photon\n"
"\n"
"Note: we don't know and don't care how much 'time'\n"
"there is per bin;  we work on a bin-by-bin basis.\n"
"\n"
"Typical next step:\n"
"  ./favia -z -x cal.dat -j 1e-6 -s 1e-6 -L 100000e-6 > cal.csv\n"
<< endl;

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
  cerr << "  total_events: " << tot__events << endl;

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
  cerr << "Input file: " << fn << endl;

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
  int scatter(0);
  pdx_t maxdraw(1000000);
  double poisson(.01);
  string ofn("cal.dat");

// Process commandline options 
  const int ALT(128);
  static struct option long_options[] = {
    {"help",            0, NULL, 'h'},
    {"maxdraw",         1, NULL, 'm'},
    {"ofile",           1, NULL, 'o'},
    {"poisson",         1, NULL, 'p'},
    {"scatter",         0, NULL, 's'},
    {"verbose",         0, NULL, 'v'},
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
      case 'm':
        maxdraw = atoll(optarg);
        break;
      case 'o':
        ofn = optarg;
        break;
      case 'p':
        poisson = atof(optarg);
        break;
      case 's':
        scatter++;
        break;
      case '?':         // optarg() uses this for any unrecognized 
                        //   option, and has already complained about it.
        cerr << "For help, try\n  " << argv[0]
             << " --help" << endl;
        exit(1);
      default:
        int chx(ch&~ALT);
	fprintf(stderr, "Sorry, option %s'%c' not yet implemented.\n", 
		ch!=chx? "ALT+" : "", chx);
	exit(1);
    }

  }
  cerr << "ofn: " << ofn << endl;
  cerr << "maxdraw: " << maxdraw << " bins" << endl;
  cerr << "poisson: " << poisson << endl;
  cerr << "verbose: " << verbose << endl;

// Define a uniform random number distribution which produces "double"
// values between 0 and 1 (0 inclusive, 1 exclusive).
  boost::lagged_fibonacci607 rng;
  rng.seed(static_cast<boost::uint32_t> (std::time(0)));
  boost::uniform_real<> uni(0.0,1.0);

 // B. Produce Uniform (0, 1)
  boost::variate_generator<boost::lagged_fibonacci607&, 
        boost::uniform_real<> >  jrandom(rng, uni);

  if (scatter) {
    for (pdx_t ii(0); ii < maxdraw; ii++) {
      cout << jrandom() 
        << ", "  << jrandom()
        << endl;
    }
  } else {
    ofstream ouch;
    ouch.open(ofn.c_str(), ios::out | ios::binary);
    if (!ouch.good()) {
      cerr << "monaco: cannot create data file '" << ofn << "': ";
      perror(0);
      exit(1);
    }

    for (cdp_t ii(0); ii < maxdraw; ii++) {
      if (jrandom() < poisson) {
        ouch.write((char*)(&ii), sizeof(ii));
        if (!ouch.good()) {
          cerr << "monaco: write error: ";
          perror(0);
          exit(1);
        }
      }
    }
    ouch.close();
  }

  return 0;
}
