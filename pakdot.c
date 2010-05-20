// Perform a dot product on the packed representation
// of two vectors.
//
// See pakdot.h for lots of important comments.

using namespace std;
#include <iostream>
#include "pakdot.h"

// Take the dot product of two packed vectors.
//
//   The type here is pdx_t, because the largest number
//   of coincidences between pv1 and pv2 is the number
//   of slots in a packed vector.
 
dot_t pakdot(const pakvec pv1, const pakvec pv2){

  dot_t rslt = 0;
  pdx_t p1 = 0;
  pdx_t p2 = 0;
  
  for(;;){	// loop over all slots in conceptual vector

// Should check to make sure that outwin returns 1 in both
// cases or returns 2 in both cases, but I'm lazy:
    if (pv1.outwin(p1)) break;	
    if (pv2.outwin(p2)) break;	

// difference in conceptual times:
    cdp_t diff =  (pv1.updata[p1].abscissa - pv1.cv0time)
		- (pv2.updata[p2].abscissa - pv2.cv0time);
    if (pv1.verbose) 
      cout << "\t"  " " << p1 << " " << p2 << " diff: " << diff << endl;

// Note that the following three if statements are intentionally 
// not mutually exclusive.  If diff is 0, both p1 and p2 get
// incremented.
    if (diff == 0) rslt +=   pv1.updata[p1].ordinate 
                           * pv2.updata[p2].ordinate;
    if (diff <= 0) p1++;
    if (diff >= 0) p2++;

  }
  return rslt;
}

// Same as the above, but assuming periodic boundary conditions.
// We use do not use updata or outwin().
// Instead we use spdata and spnp.
dot_t pbc_pakdot(const pakvec pv1, const pakvec pv2){

  dot_t rslt = 0;
  pdx_t p1 = 0;
  pdx_t p2 = 0;
  
  for(pdx_t ii=0; ;ii++) {	// loop over all slots in conceptual vector
    if (p1 >= pv1.spnp) break;	// quit if we've used every data point
    if (p2 >= pv2.spnp) break;  // ditto
// difference in conceptual times:
    cdp_t diff =  (pv1.spdata[p1].abscissa - pv1.cv0time) % pv1.cvnp
		- (pv2.spdata[p2].abscissa - pv2.cv0time) % pv2.cvnp;
    if (0 && pv1.verbose)
      cout << "\t"  " " << p1 << "/" << pv1.cvnp << " " 
      			<< p2 << "/" << pv2.cvnp << " diff: " << diff << endl;

// Note that the following three if statements are intentionally 
// not mutually exclusive.  If diff is 0, both p1 and p2 get
// incremented.
    if (diff == 0) rslt +=   pv1.spdata[p1].ordinate
        		   * pv2.spdata[p2].ordinate;
    if (diff <= 0) p1++;
    if (diff >= 0) p2++;
  }


  return rslt;
}

// The input vector must have non-decreasing
// abscissas, but not necessarily strictly
// increasing.  After coarse-graining, the
// abscissas are quite likely to be "loose"
// i.e. to have duplicates.
// The output vector will be tight, i.e.
// strictly increasing abscissas.
pdx_t tighten(
  const cdp_t first, const cdp_t grain,
  const daton* loose, const pdx_t size, daton* tight
){
  if (!size) return 0;
  pdx_t jj(0);
  tight[jj].abscissa = (loose[0].abscissa - first) / grain;
  tight[jj].ordinate = loose[0].ordinate;
  for (pdx_t ii=1; ii < size; ii++) {
    cdp_t abs = (loose[ii].abscissa - first) / grain;
    if (tight[jj].abscissa == abs) {
      tight[jj].ordinate += loose[ii].ordinate;
    } else {
      jj++;
      tight[jj].abscissa = abs;
      tight[jj].ordinate = loose[ii].ordinate;
    }
  }

  return 1+jj;
}
