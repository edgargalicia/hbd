#include "hbmap.h"
#include "utils.h"
#include "topology.h"
#include "Math/Vectors.h"
#include <iostream>
#include <utility>
#include <iomanip>

static inline bool isProton(const std::string& name) {
    return name[0] == 'h';
}

static inline bool isAcceptor(const std::string& name) {
    return name == "n" || name == "o"; // extend if needed
}


HBMap::HBMap(const std::vector<int> &list, const Topology &topo) {
  for(auto var : list) {
    if (isAcceptor(topo.AtomName()[var])) {
      acceptors.push_back(var);
    }
  }
  mapAcceptors = Mapping(acceptors, topo.AtomName().size());
  findDonors(topo);
  mapDonors = Mapping(donorOrder, topo.AtomName().size());
}

int HBMap::TotalProtons() const {
  int sum = 0;
  for(auto Protons : donors) {
      sum += Protons.second.size();
  }
  return sum;
}

void HBMap::findDonors(const Topology &topo) {
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

void HBMap::Print() {
  std::cout << "Acceptors size: " << acceptors.size() << '\n';
  std::cout << "Donors size: " << donors.size() << '\n';
  std::cout << "Protons size: " << TotalProtons() << '\n';
}

bool isHBonded( const BondAtoms &at, const HBondContext &ctx, HBGeometry &hbg ) {
  float rda2;
  float ca;
  Math::Vec3 r_da, r_dh;

  if (at.bondpair.donor == at.bondpair.acceptor) {
    return false;
  }

  // const std::vector<Math::Vec3> &x;
  // const Box &box;
  r_da = ctx.coords[at.bondpair.donor] - ctx.coords[at.bondpair.acceptor];
  /* Apply PBC */
  ctx.box.Pbc( r_da );

  rda2 = r_da.norm2();
  if (rda2 > rc2) {
    return false;
  }

  r_dh = ctx.coords[ at.bondpair.donor ] - ctx.coords[ at.proton ];
  ctx.box.Pbc( r_dh );
  ca = cos_angle( r_dh, r_da );
  /* if angle is smaller, cos is larger */
  if (ca >= ccut) {
    hbg.d_da = std::sqrt( rda2 );
    auto r_ha = ctx.coords[ at.proton ] - ctx.coords[ at.bondpair.acceptor ];
    ctx.box.Pbc( r_ha );
    hbg.d_ha = r_ha.norm();
    hbg.ang = std::acos( ca )/Pi*180.0f;
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

void PrintStats( std::ofstream &fp, const float &dt, const HBondStats &hbs, const std::unordered_set<int> &group1Set, const std::unordered_set<int> &group2Set) {
  auto g = CountGroupActivity(hbs.activeAcceptors, hbs.activeDonors, group1Set, group2Set);

  fp << std::left
     << std::setw(12) << std::fixed << std::setprecision(5) << dt
     << std::setw(12) << hbs.hbond
     << std::setw(12) << hbs.intra1
     << std::setw(12) << hbs.intra2
     << std::setw(12) << hbs.inter
     << std::setw(12) << hbs.activeAcceptors.size()
     << std::setw(12) << hbs.activeDonors.size()
     << std::setw(12) << g.g1Don
     << std::setw(12) << g.g1Acc
     << std::setw(12) << g.g2Don
     << std::setw(12) << g.g2Acc
     << '\n';
}

void addBondPresence(std::unordered_map<BondPair, std::vector<int>, BondKeyHash> &bondPresence,
                     std::unordered_map<BondPair, BondType, BondKeyHash> &bondTypeMap,
                     const BondAtoms &at, int frame,
                     const std::unordered_set<int> &group1Set,
                     const std::unordered_set<int> &group2Set) {
  auto &key = at.bondpair;
  bondPresence[key].push_back(frame);

  if (bondTypeMap.find(key) == bondTypeMap.end()) {
    if (group1Set.count(key.acceptor) && group1Set.count(key.donor)) {
      bondTypeMap[key] = BondType::Intra1;
    } else if (group2Set.count(key.acceptor) && group2Set.count(key.donor)) {
      bondTypeMap[key] = BondType::Intra2;
    } else {
      bondTypeMap[key] = BondType::Inter;
    }
  }
}

std::unordered_map<BondPair, double, BondKeyHash> computeBondPresencePercentages(
  const std::unordered_map<BondPair, std::vector<int>, BondKeyHash> &bonds, int totalFrames) {

  std::unordered_map<BondPair, double, BondKeyHash> percentages;

  for (const auto &entry : bonds) {
    const BondPair &key = entry.first;
    const auto &frames = entry.second;

    double percentage = (100.0 * frames.size()) / static_cast<int>(totalFrames);
    percentages[key] = percentage;
  }

  return percentages;
}

void HBondStats::reset() {
  hbond = 0;
  intra1 = 0;
  intra2 = 0;
  inter = 0;
  activeAcceptors.clear();
  activeDonors.clear();
}
