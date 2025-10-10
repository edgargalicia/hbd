#include "hbmap.h"
#include "topology.h"
#include <Math/Matrix.h>
#include <Math/Vectors.h>
#include <iostream>
#include <utility>

bool isProton(const std::string& name) {
    return name[0] == 'h';
}

bool isAcceptor(const std::string& name) {
    return name == "n" || name == "o"; // extend if needed
}

constexpr float rc2 = 3.5*3.5;
constexpr double Pi = 3.14159265358979323846264338327950288;
constexpr double DEG2RAD = ( Pi / 180.0 );
constexpr float acut = 30;
const float ccut = std::cos( acut * DEG2RAD );

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

HBMap::HBMap(const std::vector<int> &list, const Topology &topo) {
  for(auto var : list) {
    if (isAcceptor(topo.AtomName()[var])) {
      acceptors.push_back(var);
    }
  }
  mapAcceptors = Mapping(acceptors, topo.AtomName().size());
  FindDonors(topo);
  mapDonors = Mapping(donorOrder, topo.AtomName().size());
}

int HBMap::TotalProtons() const {
  int sum = 0;
  for(auto Protons : donors) {
      sum += Protons.second.size();
  }
  return sum;
}

void HBMap::FindDonors(const Topology &topo) {
  std::unordered_set<int> lookup(acceptors.begin(), acceptors.end());
  for(const auto &bond : topo.Bonds()) {
    int a = bond.first;
    int b = bond.second;

    if (lookup.count(a) == 0 && lookup.count(b) == 0) continue;

    if ( isProton(topo.AtomName()[b]) ) {
      if (donors.find(a) == donors.end()) donorOrder.push_back(a);
      donors[a].push_back(b);
    } else if ( isProton(topo.AtomName()[a]) ) {
      if (donors.find(b) == donors.end()) donorOrder.push_back(b);
      donors[b].push_back(a);
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
  std::cout << "Acceptors size: " << acceptors.size() << '\n';
  std::cout << "Donors size: " << donors.size() << '\n';
  std::cout << "Protons size: " << TotalProtons() << '\n';
}

bool isHBonded( int don, int prot, int acc, std::vector<Math::Vec3> &x, const Box &box,
              float &d_ha, float &ang ) {
  float rda2;
  float ca;
  Math::Vec3 r_da, r_dh;

  if (don == acc) {
    return false;
  }

  r_da = x[don] - x[acc];
  /* Apply PBC */
  box.Pbc( r_da );

  rda2 = r_da.norm2();
  if (rda2 > rc2) {
    return false;
  }

  r_dh = x[ don ] - x[ prot ];
  box.Pbc( r_dh );
  ca = cos_angle( r_dh, r_da );
  /* if angle is smaller, cos is larger */
  if (ca >= ccut) {
    d_ha = std::sqrt( rda2 );
    ang = std::acos( ca );
    return true;
  }
  return false;
}

void PrintTimeSeries( std::ofstream &fp, int acc, int don, int prot, float d_ha, float ang ) {
  fp << acc << " - " << don << " : " << prot << " dist: " << d_ha << " ang: " << ang/3.1416*180 << '\n';
}

GroupActivity CountGroupActivity(const std::unordered_set<int> &activeAcceptors, const std::unordered_set<int> &activeDonors, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set) {
  GroupActivity g;
  for (int a : activeAcceptors) {
    if (group1Set.count(a)) ++g.g1Acc;
    if (group2Set.count(a)) ++g.g2Acc;
  }

  for (int d : activeDonors) {
    if (group1Set.count(d)) ++g.g1Don;
    if (group2Set.count(d)) ++g.g2Don;
  }

  return g;
}

void PrintStats( std::ofstream &fp, int hbond, int intra1, int intra2, int inter, const std::unordered_set<int> &activeAcceptors, const std::unordered_set<int> &activeDonors, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set ) {

  auto g = CountGroupActivity(activeAcceptors, activeDonors, group1Set, group2Set);

  fp << std::left
     << std::setw(12) << hbond
     << std::setw(12) << intra1
     << std::setw(12) << intra2
     << std::setw(12) << inter
     << std::setw(12) << activeAcceptors.size()
     << std::setw(12) << activeDonors.size()
     << std::setw(12) << g.g1Don
     << std::setw(12) << g.g1Acc
     << std::setw(12) << g.g2Don
     << std::setw(12) << g.g2Acc
     << '\n';
}

void addBondPresence(std::unordered_map<BondKey, std::vector<int>, BondKeyHash> &bondPresence,
                     std::unordered_map<BondKey, BondType, BondKeyHash> &bondTypeMap,
                     int acc, int don, int frame,
                     const std::unordered_set<int> &group1Set,
                     const std::unordered_set<int> &group2Set) {
  BondKey key {acc, don};
  bondPresence[key].push_back(frame);

  if (bondTypeMap.find(key) == bondTypeMap.end()) {
    if (group1Set.count(acc) && group1Set.count(don))
      bondTypeMap[key] = BondType::Intra1;
    if (group2Set.count(acc) && group2Set.count(don))
      bondTypeMap[key] = BondType::Intra2;
    else
      bondTypeMap[key] = BondType::Inter;
  }
}

std::unordered_map<BondKey, double, BondKeyHash> computeBondPresencePercentages(
  const std::unordered_map<BondKey, std::vector<int>, BondKeyHash> &bonds, int totalFrames) {

  std::unordered_map<BondKey, double, BondKeyHash> percentages;

  for (const auto &entry : bonds) {
    const BondKey &key = entry.first;
    const auto &frames = entry.second;

    double percentage = (100.0 * frames.size()) / static_cast<int>(totalFrames);
    percentages[key] = percentage;
  }

  return percentages;
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
