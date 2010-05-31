//////////////

using namespace std;

#include <stdlib.h>     /* for atoll */
#include <iostream>
#include <stdint.h>     /* for int64_t */
#include <assert.h>

int64_t foo;
int64_t fudge(0);

int main(int argc, char** argv){
  assert(sizeof(long long int) == sizeof(int64_t));
  if (argc > 1) {
    fudge = atoll(argv[1]);
    cerr << fudge << endl;
  }
  for(;;) {
    cin.read((char*)(&foo), sizeof(foo));
    if (!cin.good()) break;
    foo += fudge;
    cout.write((char*)(&foo), sizeof(foo));
  }
}
