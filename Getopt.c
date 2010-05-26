//////////////////////////////////////////////////////////////////////
// Provide a C++ interface to getopt.
// Calculate things that can be calculated, to ensure consistency,
// and so the programmer doesn't need to duplicate effort.

using namespace std;

#include <iostream>     /* needed by some compilers */
#include "Getopt.h"
#include <iostream>
#include <string>

int getopt_long(int argc, char * const argv[],
                  const struct option *longopts, int *longindex){
  
  string optstring;
  const struct option * pp;
  for (pp = longopts; pp->name; pp++){
    if (pp->val < 0 || pp->val > 255) {
      cerr << "Cannot handle val " <<  pp->val
	<<  " for option '" << pp->name << "'" << endl;
    } else {
//      cerr << pp->name << "'" << char(pp->val) << "'" << endl;
      optstring += char(pp->val);
      if (pp->has_arg == 1) optstring += ':';
      if (pp->has_arg == 2) optstring += "::";
    }
  }
//  cerr << optstring << endl;
  return getopt_long(argc, argv, optstring.c_str(), longopts, longindex);
}
