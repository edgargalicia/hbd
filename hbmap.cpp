#include "hbmap.h"
#include "topology.h"
#include <Math/Matrix.h>
#include <Math/Vectors.h>
#include <iostream>
#include <unordered_set>

bool isProton(const std::string& name) {
    return name[0] == 'h';
}

bool isAcceptor(const std::string& name) {
    return name == "n" || name == "o"; // extend if needed
}

float rc2 = 3.5*3.5;
float ccut = 30;

double invsqrt( double x ) { return 1.0 / std::sqrt( x ); }

float cos_angle( const Math::Vec3 a, const Math::Vec3 b ) {
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

HBMap::HBMap(const std::string &str, const Topology &topo) {
  auto preList = String2IntList(str);
  for(auto var : preList) {
    if (isAcceptor(topo.AtomName()[var])) {
      Acceptors.push_back(var);
    }
  }
  MapAcceptors = Mapping(Acceptors, topo.AtomName().size());
  FindDonors(topo);
}

int HBMap::TotalProtons() const {
  int sum = 0;
  for(auto Protons : Donors) {
      sum += Protons.second.size();
  }
  return sum;
}

void HBMap::FindDonors(const Topology &topo) {
  std::unordered_set<int> lookup(Acceptors.begin(), Acceptors.end());
  for(const auto &bond : topo.Bonds()) {
    int a = bond.first;
    int b = bond.second;

    if (lookup.count(a) == 0 && lookup.count(b) == 0) continue;

    if ( isProton(topo.AtomName()[b]) ) {
      Donors[a].push_back(b);
    } else if ( isProton(topo.AtomName()[a]) ) {
      Donors[b].push_back(a);
    }
  }
}

std::vector<int> Mapping(const std::vector<int> &list, int natoms) {
  std::vector<int> map(natoms,-1);;
  for(size_t i = 0; i != list.size(); ++i) {
    map[ list[i] ] = i;
  }
  return map;
}

std::vector<int> Match(const std::vector<int> &atomList, const Topology &topology, const std::string &str) {
  std::vector<int> result;
  for(size_t i = 0; i != atomList.size(); ++i) {
    if (topology.AtomName()[i].compare(str) == 0) {
      result.push_back(i);
    }
  }
  return result;
}

void HBMap::Print() {
  std::cout << "Acceptors size: " << Acceptors.size() << '\n';
  std::cout << "Donors size: " << Donors.size() << '\n';
  std::cout << "Protons size: " << TotalProtons() << '\n';

  // int mono = 0;
  // int di = 0;
  // int tri = 0;
  // for(const auto &var : donorToProtons) {
  //   if (var.second.size() == 1) {
  //     mono++;
  //   } else if (var.second.size() == 2) {
  //     di++;
  //   } else {
  //     tri++;
  //   }
  // }
  // std::cout << "Mono-Protons size: " << mono << '\n';
  // std::cout << "Di-Protons size: " << di << '\n';
  // std::cout << "Tri-Protons size: " << tri << '\n';
}

bool isHBonded( int d, int h, int a, Math::Vec3 *x[], const Box &box,
              float &d_ha, float &ang ) {
  float rda2;
  float ca;
  Math::Vec3 r_da, r_dh;

  if (d == a) {
    return false;
  }

  r_da = *x[d] - *x[a];
/* Apply PBC */
  box.Pbc( &r_da );

  rda2 = r_da.norm2();
  if (rda2 > rc2) {
    return false;
  }

  r_dh = *x[ d ] - *x[ h ];
  box.Pbc( &r_dh );
  ca = cos_angle( r_dh, r_da );
  /* if angle is smaller, cos is larger */
  if (ca >= ccut) {
    d_ha = std::sqrt( rda2 );
    ang = std::acos( ca );
    return true;
  }
  return false;
}

// void read_acceptors( const std::string &str, HBData &hb ) {
//   std::stringstream ss( str );
//   std::vector<int> arr;
//   std::string token;
//
//   while (ss >> token) {
//     if (token.back( ) == ',') token.pop_back( );
//     if (token[ 0 ] == '-') {
//       ss >> token;
//       if (token.back( ) == ',') token.pop_back( );
//       for (int i = arr.back( ) + 1; i <= std::stoi( token ); i++) {
//         arr.push_back( i );
//       }
//     } else {
//       arr.push_back( std::stoi( token ) );
//     }
//   }
//   hb.nacc = static_cast<int>( arr.size( ) );
//   hb.a = new int[ hb.nacc ];
//   for (int i = 0; i < hb.nacc; i++) {
//     hb.a[ i ] = arr[ i ];
//   }
//   arr.clear( );
//   std::cout << "Number of acceptors: " << hb.nacc << std::endl;
// }
//
// void search_donors( const std::string &str, const Box &box, rvec *x[],
//                     HBData &hb ) {
//   std::stringstream ss( str );
//   std::string token;
//
//   float rah2;
//   float rc2 = dOH * dOH;
//   rvec r_ah;
//   std::vector<int> arr;
//   std::vector<int> vprot;
//   std::vector<int> vnh;
//
//   while (ss >> token) {
//     if (token.back( ) == ',') token.pop_back( );
//     if (token[ 0 ] == '-') {
//       ss >> token;
//       if (token.back( ) == ',') token.pop_back( );
//       for (int i = arr.back( ) + 1; i <= std::stoi( token ); i++) {
//         arr.push_back( i );
//       }
//     } else {
//       arr.push_back( std::stoi( token ) );
//     }
//   }
//
//   const int nprot = static_cast<int>( arr.size( ) );
//   std::vector<int> hh;
//   std::copy( arr.begin( ), arr.end( ), std::back_inserter( hh ) );
//   arr.clear( );
//
//   /* Search for donors */
//   for (int ia = 0; ia < hb.nacc; ++ia) {
//     int nppd = 0;
//     for (int ih = 0; ih < nprot; ++ih) {
//       rvec_sub( *x[ hb.a[ ia ] ], *x[ hh[ ih ] ], r_ah );
//       pbc( &r_ah, box );
//
//       rah2 = norm2( r_ah );
//       if (rah2 < rc2) {
//         vprot.push_back( hh[ ih ] );
//         nppd++;
//       }
//     }
//     if (nppd > 0) {
//       arr.push_back( hb.a[ ia ] );
//       vnh.push_back( nppd );
//     }
//   }
//
//   hb.ndon = static_cast<int>( arr.size( ) );
//   hb.d = new Donor[ hb.ndon ];
//   for (int i = 0, k = 0; i < hb.ndon; ++i) {
//     hb.d[ i ].don = arr[ i ];
//     hb.d[ i ].nh = vnh[ i ];
//     hb.d[ i ].hh = new int[ vnh[ i ] ];
//     for (int j = 0; j < vnh[ i ]; ++j) {
//       hb.d[ i ].hh[ j ] = vprot[ k++ ];
//     }
//   }
//
//   arr.clear( );
//   vprot.clear( );
//   vnh.clear( );
//
//   std::cout << "Number of donors: " << hb.ndon << std::endl;
//   std::cout << "Number of protons: " << nprot << std::endl << std::endl;
// }
//
