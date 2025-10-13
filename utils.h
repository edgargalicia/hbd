#include "Math/Vectors.h"
#include <cmath>

static constexpr float rc2 = 3.5*3.5;
static constexpr double Pi = 3.14159265358979323846264338327950288;
static constexpr double DEG2RAD = ( Pi / 180.0 );
static constexpr float acut = 30;
static const float ccut = std::cos( acut * DEG2RAD );

static inline double invsqrt( double x ) { return 1.0 / std::sqrt( x ); }

static inline float cos_angle( const Math::Vec3 a, const Math::Vec3 b ) {
  /*
   *                  ax*bx + ay*by + az*bz
   * cos-vec (a,b) =  ---------------------
   *                      ||a|| * ||b||
   */
  float cosval;
  int m;
  double aa, bb, ip, ipa, ipb, ipab; /* For accuracy these must be double! */

  ip = ipa = ipb = 0.0;
  for (m = 0; ( m < 3 ); m++) {
    aa = a[ m ];
    bb = b[ m ];
    ip += aa * bb;
    ipa += aa * aa;
    ipb += bb * bb;
  }
  ipab = ipa * ipb;
  if (ipab > 0) {
    cosval = ip * invsqrt( ipab ); /*  7 */
  } else {
    cosval = 1;
  }
  /* 25 TOTAL */
  if (cosval > 1.0) {
    return 1.0;
  }
  if (cosval < -1.0) {
    return -1.0;
  }

  return cosval;
}

static inline std::vector<int> Mapping(const std::vector<int> &list, int natoms) {
  std::vector<int> map(natoms,-1);;
  for(size_t i = 0; i != list.size(); ++i) {
    map[ list[i] ] = i;
  }
  return map;
}

static inline std::vector<int> Match(const std::vector<int> &atomList, const std::vector<std::string> &AtomName, const std::string &str) {
  std::vector<int> result;
  for(size_t i = 0; i != atomList.size(); ++i) {
    if (AtomName[i].compare(str) == 0) {
      result.push_back(i);
    }
  }
  return result;
}

