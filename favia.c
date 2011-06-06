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

void usage(ostream& foo)
{
        foo << 
                "Perform correlations on FRET data.\n"
                "For conceptual foundations, see\n"
                "  http://www.av8n.com/physics/correlation-norm.htm\n"
                "\n"
                "Typical usage:\n"
                "  ./favia [options] > myfile.corr\n"
                "\n"
                "Long-form options include the following.\n"
                "Defaults are shown in square brackets.\n"
                "  --xfile=sss          first input file\n"
                "  --yfile=sss          second input file [clone of xfile]\n"
                "  --jiffy=ddd          bin size (in s) for input timestamps [4e-12]\n"
                "  --max_events=iii     pretend each input file has at most this many events\n"
                "                       [0 means use all of the input data]\n"
                "                       (probably --tspan models the physics better)\n"
                "  --fineness=iii       inverse grain size (in grains per octave) [8]\n"
                "  --long_lag=ddd       longest lag (in seconds) to be considered [10]\n"
                "  --short_grain=ddd    smallest grain (in seconds) to be used [1e-8]\n"
                "  --tspan=ddd          exact span of both input files (in seconds)\n"
                "                       [0 means estimate span from last event + leadout]\n"
                "  --zerospike          start output at lag=0 [otherwise at lag=1 grain]\n"
                "  --help               print this message and exit immediately\n"
                "  --verbose            increase verbosity of debugging messages\n" 
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
                "  lag\\tloglag\\tdot\\tdotnormed\\tbar^2\\tresidual\n"
                "where \\t represents a tab and,:\n"
                "  lag       is measured in seconds\n"
                "  loglag    is log10(lag)\n"
                "  dot       is the raw dot product, i.e. the correlation at this lag\n"
                "  dotnormed is the normalized dot product\n"
                "  bar^2     is a measure of uncertainty as detailed by wohland"
                "  residual  is measured relative to an estimate of the large-lag asymptote\n"
                "\n"
                "Useful simple checks:\n"
                "  ./favia -z -j 1 -s 1 -L 12 -t 112 -x count10.dat -y count1.dat\n"
                "  ./favia -z -j 1 -s 1 -L 12 -t 112 -y count10.dat -x count1.dat\n"
                "  ./favia -z -j 1 -s 1 -L 140 -t 300 -x count100.dat -y count1.dat\n"
                "  # where -L 140 is the same as -L 128\n"
                "  # and   -t 300 is the same as -t 288\n"
                << endl;

        // Not implemented:
        // "  PERIODIC=1   use periodic boundary conditions\n"
        // "               [0 means pad with zeros]\n"

}

const char* Getenv(const char* key, const char* dflt)
{
        char* get = getenv(key);
        if (!get) return dflt;
        return get;
}


vector<daton> readfile_bin(const string fn, const int max__events)
{
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

vector<daton> readfile_asc(const string fn, const int max__events)
{
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
        cerr << "Ascii input file: " << raw__data.size() 
             << " lines" << endl;
        return raw__data;
}

//////////
// Auto-detect the input file format,
// then call the appropriate handler.
vector<daton> readfile(const string fn, const int max__events)
{
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

class channel
{
public:
        vector<daton> raw_data;
        pdx_t tot_events;
        vector<daton> cg_data;
        cdp_t front_event;
        cdp_t back_event;
        cdp_t bspan_mi_bins;
        cdp_t qspan_mi_bins;
        double spacing_mi_bins;       // bins per event, average spacing
        // roughly speaking, inverse of the average rate

        // constructor:
        channel(const string fn, const int max__events,
                const double jiffy,
                const cdp_t long_grain_mi_bins,
                const double tspan)
        {
                raw_data = readfile(fn, max__events);
                tot_events = raw_data.size();

                // Allocate space for coarse-grained data;
                // For safety, overestimate the amount of space needed:
                // the size of coarse-grained array can never exceed the
                // size of the raw-data array.
                cg_data.resize(tot_events);
                front_event = raw_data[0].abscissa;
                back_event = raw_data[tot_events-1].abscissa;


                // span of actual data, 
                // i.e. the sum of the N-1 inter-event intervals
                // ... not used ...
                //  cdp_t dspan_mi_bins = 1 + back_event - front_event;

                cerr << boost::format(
                                "  front_event: %16LLf == %13.6g bins == %13.6g s")
                        % front_event % double(front_event)
                        % (front_event*jiffy)
                     << endl;
                cerr << boost::format(
                                "  back_event:  %16LLf == %13.6g bins == %13.6g s")
                        % back_event % double(back_event)
                        % (back_event*jiffy)
                     << endl;


                /*****************************
                 *
                 * Estimating the leadout......
                 *
                 * Note: If tspan has been specified, you can skip this entire
                 * discussion.
                 *
                 * Later, we will do a normalization calculation.  We will need
                 * to feed it an unbiased estimate of the Poisson rate for this
                 * channel.  If we had a definite number of points in a
                 * definite interval, this would be easy ... but unless tspan
                 * has been specified, we don't have a definite interval.  All
                 * we have is the timestamp of the last event.
                 *
                 * Rather than talking about the rate, let's talk about inverse
                 * rate, i.e. the time between pulses.  If we have N events,
                 * there are N+1 consecutive time intervals of interest:
                 *
                 * -- the leadin, i.e. the time before the front event; -- N-1
                 * inter-event intervals; and -- the leadout.
                 *
                 * Alas we do not know the leadout.  So we need a method that
                 * (implicitly or explicitly) estimates the leadout, and then
                 * estimates the rate.  So we are using the data to estimate
                 * two things.  Our estimate of the rate will of course depend
                 * on our estimate of the leadout.
                 *
                 * As a first attempt, let's just use the N-1 inter-event
                 * intervals, then the average interval is (back_event -
                 * front_event) / (N-1).  That is a disaster when N=1, but we
                 * can improve it as follows.
                 *
                 * As a second (and almost final) attempt, we make use of the
                 * leadin interval.  This is useful information, but it is not
                 * an unbiased estimate of the inter-event interval, since the
                 * logger presumably turned on in the middle of an interval,
                 * and this samples intervals in a biased way, favoring longer
                 * intervals.  So let's count it as half of an interval.  So
                 * our estimate of the average interval is then (back_event -
                 * 0) / (N-1/2).
                 *
                 * Note that when N=1, this implicitly estimates the leadout as
                 * being equal to the leadin.  More generally, this estimates
                 * the leadout as back_event / (2N-1).  For large N this
                 * converges to half of the average inter-event interval.  All
                 * this seems entirely reasonable.
                 *
                 * Finally, we add 1 to the span_mi_bins, because we want the
                 * span-time to cover all the time from the start of the front
                 * bin to the *end* of the back bin.
                 *
                 */ 

                if (tspan) {
                        bspan_mi_bins = floor(0.5 + tspan / jiffy);
                } else {
                        // integer arithmetic here; remainder (if any) gets thrown away:
                        cdp_t leadout = back_event / (2*tot_events - 1);
                        bspan_mi_bins = (back_event - 0) + leadout + 1;
                        cerr << boost::format(
                                        "  leadout:     %16LLf == %13.6g bins == %13.6g s")
                                % leadout % double(leadout)
                                % (leadout*jiffy)
                                << endl;
                }

                qspan_mi_bins = long_grain_mi_bins * (bspan_mi_bins / long_grain_mi_bins);
                cerr << boost::format(
                                "  qspan:       %16LLf == %13.6g bins == %13.6g s")
                        % qspan_mi_bins % double(qspan_mi_bins)
                        % (qspan_mi_bins*jiffy)
                        << endl;

                spacing_mi_bins = bspan_mi_bins / tot_events;
                cerr << boost::format(
                                "  avg spacing: %16.0f == %13.6g bins == %13.6g s")
                        % spacing_mi_bins % double(spacing_mi_bins) % (spacing_mi_bins*jiffy)
                        << endl;
        }     /* end constructor */
};      /* end class channel */

int main(int argc, char** argv)
{
        int verbose(0);
        string xfn;
        string yfn;
        int fineness(8);
        pdx_t max_events(0);

        double long_lag_req(10.0);
        double short_grain(1e-8);
        double jiffy(4e-12);
        double tspan(0.0);
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
                {"short_grain",     1, NULL, 's'},
                {"tspan",           1, NULL, 't'},
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
                                long_lag_req = atof(optarg);
                                break;
                        case 'm':
                                max_events = atoi(optarg);
                                break;
                        case 's':
                                short_grain = atof(optarg);
                                break;
                        case 't':
                                tspan = atof(optarg);
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
                        case '?':       // optarg() uses this for any unrecognized 
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
        if (xfn.length()==0) {
                fprintf(stderr, "Expected --xfile\n");
                exit(1);
        }
        if (yfn.length()==0) yfn = xfn;

        // Using floats for loop control is bad practice;
        // round to integral # of bins at first opportunity.

        // In this program, _mi_ means "... measured in units of ..."

        cdp_t shortgrain_mi_bins = floor(0.5 + short_grain / jiffy);

        cerr << "verbose: " << verbose << endl;
        cerr << "xfn: " << xfn << endl;
        cerr << "yfn: " << yfn << endl;
        cerr << "jiffy: " << jiffy << " s" << endl;
        cerr << "fineness: " << fineness << " grains per octave"
                " == " << fineness*log2(10) << " grains per decade"
             << endl;
        cerr << "max_events: " << max_events << endl;
        cerr << "short_grain: " << short_grain << " s" 
             << " --> " << shortgrain_mi_bins << " bins"
             << " == " << shortgrain_mi_bins*jiffy << " s"
             << endl;

        // Be a little bit defensive:
        if (shortgrain_mi_bins < 1) {
                cerr << "Ooops: grain size less than bin size." << endl;
                exit(1);
        }

        // Now need to calculate the appropriate longest grain
        //  and longest lag.

        // Picture of the grains near the diagonal,
        // as multiples of the shortest grain,
        // assuming fineness=4:
        // ||||||||| | | | |   |   |   |   |

        // 000000000 1 1 1 1   2   2   2   3  (bin
        // 012345678 0 2 4 6   0   4   8   2   number)
        // ----0000111111112222222222222222  (octave number)
        // 00000000111111112222222222222222  (octnox)

        // Round to nearest bin, 
        // to avoid floating point representation problems:
        cdp_t longlag_req_mi_bins = floor(0.5 + long_lag_req / jiffy);

        // This is the "requested" long lag, quantized in short grains.
        // (In contrast, the real long lag will be quantized in long grains,
        // but we can't do that until we figure out what the long grain is.)
        // Integer arithmetic rounds down:
        cdp_t longlag_req_mi_shortg = longlag_req_mi_bins / shortgrain_mi_bins;

        // Octave number, valid in the geometric region:
        int octno(floor(log2(double(longlag_req_mi_shortg) / fineness)));
        // Treat the initial linear (non-geometric region) as
        // part of octave zero:
        int octnox(octno > 0 ? octno : 0);
        cdp_t longgrain_mi_bins = shortgrain_mi_bins * (1 << octnox);

        cdp_t longlag_mi_longg = longlag_req_mi_bins / longgrain_mi_bins;
        cdp_t longlag_mi_bins = longlag_mi_longg * longgrain_mi_bins;

        cerr << "long_lag: "  << long_lag_req  << " s" 
             << " --> " 
             << longlag_mi_bins << " bins" 
             << " == " << longlag_mi_bins / shortgrain_mi_bins << " shorts"
             << " == " << longlag_mi_bins / longgrain_mi_bins << " longs"
             << " == " << longlag_mi_bins * jiffy << " s"
             << endl;

        cerr  << "  longlag_req_mi_shortg: " << longlag_req_mi_shortg
              << "  log2: "  << log2(longlag_req_mi_shortg)
              << "  octno: "  << octno
              << "  octnox: "  << octnox
              << endl;

        cerr << "longgrain: " << longgrain_mi_bins << " bins"
             << " == " << longgrain_mi_bins / shortgrain_mi_bins << " shorts"
             << " == " << longgrain_mi_bins * jiffy << " s"
             << endl;

        // sanity check:
        if (longlag_mi_bins * jiffy <= long_lag_req
                        && (longlag_mi_bins + longgrain_mi_bins) * jiffy > long_lag_req){
                /* OK */
        } else {
                cerr << "Ooops: inconsistency involving"
                        " longlag_mi_bins, longgrain_mi_bins, and long_lag_req"
                     << endl;
                exit(1);
        }

        channel x(xfn, max_events, jiffy, longgrain_mi_bins, tspan);
        channel y(yfn, max_events, jiffy, longgrain_mi_bins, tspan);

        // there are three requirements on the zone size:
        // *) Zonesize must be a multiple of the largest grain size.
        //    Not largest octave, but largest grain within that octave.
        // *) Zonesize + longest lag <= X data size
        // *) Zonesize <= Y data size
        //
        // The treatment of X and Y is not symmetric, because 
        // only X gets lagged.

        cdp_t zone_mi_bins = x.qspan_mi_bins;
        // Here is where we actually implement extended boundary conditions:
        // Reserve one long_lag worth of room at the end:
        zone_mi_bins -= longlag_mi_bins;
        if (zone_mi_bins > y.qspan_mi_bins) zone_mi_bins = y.qspan_mi_bins;
        cerr << "zone: " << zone_mi_bins << " bins" 
             << " == " << zone_mi_bins / shortgrain_mi_bins << " shorts"
             << " == " << zone_mi_bins / longgrain_mi_bins << " longs"
             << " == " << zone_mi_bins * jiffy << " s"
             << endl;

        cerr << "longlag: " << longlag_mi_bins << " bins" 
             << " == " << longlag_mi_bins / shortgrain_mi_bins << " shorts"
             << " == " << longlag_mi_bins / longgrain_mi_bins << " longs"
             << " == " << longlag_mi_bins * jiffy << " s"
             << endl;

        if (zone_mi_bins <= 0) {
                cerr << "Oops, no data in zone.\n" << endl;
                exit(1);
        }

        // The fineness is the number of grains in an ordinary octave, but
        // the startup transient is a linear progression, not a geometric
        // or pseudogeometric progression, so it is not an octave at all,
        // and it requires 2*fineness grains.
        // Also, in the typical octave, curlag_mi_grains takes on 
        // values in the range [fineness, 2*fineness-1] ... so we will 
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
        cdp_t didlag_mi_bins;

        // main loop over all resolutions
        for (pdx_t kk=0; ; kk++) {
                // here with the following loop variables already set up
                //   curlag_mi_grains = ...
                //   curgrain_mi_bins = ...

                double zone_mi_grains = zone_mi_bins / curgrain_mi_bins;
                double curgrain_mi_sec(jiffy * curgrain_mi_bins);

                double lag = curlag_mi_grains * curgrain_mi_sec;

                if (curlag_mi_grains * curgrain_mi_bins > longlag_mi_bins) break;

                // remap so data starts at zero, with coarse graining:

                pdx_t x_small = tighten(0, curgrain_mi_bins, 
                                &x.raw_data[0], x.raw_data.size(), 
                                &x.cg_data[0]);
                pakvec pv1(0, zone_mi_bins/curgrain_mi_bins, &x.cg_data[0], x_small);

                pdx_t y_small = tighten(0, curgrain_mi_bins, 
                                &y.raw_data[0], y.raw_data.size(), 
                                &y.cg_data[0]);
                pakvec pv2(0, zone_mi_bins/curgrain_mi_bins, &y.cg_data[0], y_small);

                double old_norm_denom = zone_mi_grains
                        * (curgrain_mi_sec / (x.spacing_mi_bins*jiffy))
                        * (curgrain_mi_sec / (y.spacing_mi_bins*jiffy));

                cerr << "top of loop:  lag_mi_grains: " << curlag_mi_grains
                     << "  curgrain_mi_bins: " << curgrain_mi_bins
                     << "  iniial lag: " << lag << " sec"
                     << "  old_denom: " << old_norm_denom
                     << endl;
                cerr << "zone_mi_grains: " << zone_mi_grains
                     << "  upnp1: " << pv1.upnp
                     << "  upnp2: " << pv2.upnp << endl;

                cerr.flush();
                // inner loop over all lags, stepping grain by grain:
                for (; ; curlag_mi_grains++) {
                        // normally we break at the end of an octave ...
                        if (curlag_mi_grains >= fin2) break;
                        // ... but we bail out in mid-octave if we reach the end of the zone:
                        if (curlag_mi_grains * curgrain_mi_bins > longlag_mi_bins)
                                goto main_done;

                        lag = curlag_mi_grains * curgrain_mi_sec;
                        if (0) cerr << "curgrain_mi_bins: " << curgrain_mi_bins
                                    << "  curlag_mi_grains: " << curlag_mi_grains
                                    << "  lag: " << curlag_mi_grains*curgrain_mi_sec << endl;
                        pv1.set_bin0time(curlag_mi_grains);
                        pakdot_result_t dot = pakdot(pv1, pv2);
                        double a = 0, b = 0;
                        for (pdx_t i=0; pv1.outwin(i) == 0; i++)
                                a += pv1.updata[i].ordinate;
                        for (pdx_t i=0; pv2.outwin(i) == 0; i++)
                                b += pv2.updata[i].ordinate;
                        double norm_denom = zone_mi_grains * (a / zone_mi_grains) * (b / zone_mi_grains);
                        didlag_mi_bins = curlag_mi_grains * curgrain_mi_bins;
                        if (curlag_mi_grains == 0) zerospike = dot.dot;
                        hits += dot.dot;
                        double fakelag(lag);
                        if (fakelag == 0.0) fakelag = jiffy/10.;
                        else if (fakelag < 0.0) fakelag = jiffy/100.;
                        double loglag(log10(fakelag));

                        double dotnormed = double(dot.dot) / norm_denom;
                        double squaresumnormed = double(dot.sum_squares) / pow(norm_denom, 2);

                        double bar2 = (zone_mi_grains * squaresumnormed - pow(dotnormed, 2)) / zone_mi_grains;
                        double model = 1.0;
                        double residual = dotnormed - model;
                        // If you change this output statement, be sure to
                        // change the usage() message to match:
                        *ouch << boost::format
                                ("%15.8e\t%15.8e\t%10Ld\t%15.8e\t%15.8e\n")
                                % fakelag % loglag % dot.dot % dotnormed % bar2;

                        //residual;
                        if (dot.dot) if (verbose)  cerr << boost::format
                                ("curgrain_mi_bins:%12Ld curlag_mi_grains: %12Ld"
                                 "  fakelag: %15.8e  dot: %10Ld  dot/x: %15.8e\n")
                                        % curgrain_mi_bins % curlag_mi_grains 
                                        % fakelag % dot.dot % (double(dot.dot)/double(norm_denom));
                }

                // prepare for the next iteration:
                curlag_mi_grains /= 2;
                curgrain_mi_bins *= 2;
        }

main_done:
          if (didlag_mi_bins != longlag_mi_bins) {
                  cerr << "Ooops: wrong loop exit condition:\n" 
                       << "longlag: " << longlag_mi_bins << " bins" 
                       << " == " << longlag_mi_bins / shortgrain_mi_bins << " shorts"
                       << " == " << longlag_mi_bins / longgrain_mi_bins << " longs"
                       << " == " << longlag_mi_bins * jiffy << " s"
                       << "\n"

                       << " didlag: " <<  didlag_mi_bins << " bins" 
                       << " == " <<  didlag_mi_bins / shortgrain_mi_bins << " shorts"
                       << " == " <<  didlag_mi_bins / longgrain_mi_bins << " longs"
                       << " == " <<  didlag_mi_bins * jiffy << " s"
                       << endl;
                  exit(1);
          } else {
                  if (0) cerr << "OK: didlag == longlag == " << longlag_mi_bins
                              << endl;
          }

          if (zerospike) {
                  dot_t other = hits - zerospike;
                  cerr << boost::format("Zero spike: %d  other: 2*%d  total hits: %d\n")
                          % zerospike % other % (zerospike + 2*other);
          }

          return 0;
}
