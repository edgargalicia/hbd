#include "logic.h"
#include <Math/Matrix.h>
#include <fstream>
#include <iostream>
#include <params.h>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <vector>


Logic::Logic() : filename{"input.dat"}, conf{filename} {}

void Logic::Init() {
  std::cout << "===========================\n";
  std::cout << "Read settings of trajectory\n";
  conf.Print();

  box = InitBox(conf.Box());

  std::cout << "===========================\n";
  std::cout << "Read topology\n";
  topo.Read(conf, box);
  topo.Print();

  frame.Init(conf.Size());

  hbmap = HBMap(String2IntList( conf.GetGroup1() ), topo);
  hbmap.Print();
  hbmap = HBMap(String2IntList( conf.GetGroup2() ), topo);
  hbmap.Print();
}

void Logic::Run() {
  std::cout << "===========================\n";
  std::cout << "CYCLE:\n";
  infile.open(conf.TrajName());
  if (!infile.is_open()) {
    throw std::runtime_error( "Failed to open trajectory file: " + conf.TrajName() );
  }

  int hbond;
  int intra1;
  int intra2;
  int inter;
  std::unordered_set<int> activeAcceptors;
  std::unordered_set<int> activeDonors;
  float d_ha;
  float ang;
  std::cout << "\tAssess presence of bonds\n";
  std::cout << "\tGather statistics\n";
  std::vector<int> sel;
  // while (frame.Read(infile)) {

  auto group1 = String2IntList( conf.GetGroup1() );
  auto group2 = String2IntList( conf.GetGroup2() );
  std::unordered_set<int> group1Set(group1.begin(), group1.end());
  std::unordered_set<int> group2Set(group2.begin(), group2.end());
  std::ofstream fpTime("timeseries.txt");
  std::ofstream fpStats("stats.txt");
  std::unordered_map<BondKey, std::vector<int>, BondKeyHash> bond_presence;
  for (size_t i = 0; i != 1000; ++i) {
    hbond = 0;
    intra1 = 0;
    intra2 = 0;
    inter = 0;
    activeAcceptors.clear();
    activeDonors.clear();
    frame.Read(infile);
    if (frame.Step % 1000 == 0) {
      std::cout << "Frame step: " << frame.Step << '\r' << std::flush;
    }

    // std::cout << "MAP\n";
    // std::cout << hbmap.getAccMap()[140] << '\n';
    // std::cout << hbmap.getAcceptors()[28] << '\n';
    // std::cout << hbmap.getDonMap()[115] << '\n';
    // std::cout << hbmap.getDonors()[3] << '\n';
    // std::cout << hbmap.getDonUno().find(115)->first << '\n';
    // std::cout << hbmap.getDonUno().find(115)->second[0] << '\n';
    // std::cout << hbmap.getDonUno().find(115)->second[1] << '\n';
    // std::cout << isHBonded(115, 179, 140, frame.Coords, box, d_ha, ang) << '\n';

    Select(sel, frame);
    hbmap = HBMap(sel, topo);
    // hbmap.Print();

    for(auto acc : hbmap.getAcceptors()) {
      for(const auto &donorKV : hbmap.getDonUno()) {
        const int &don = donorKV.first;
        for(const auto &prot : donorKV.second) {
          if ( isHBonded(don, prot, acc, frame.Coords, box, d_ha, ang)  ) {
            PrintTimeSeries(fpTime, acc, don, prot, d_ha, ang);
            addBondPresence(bond_presence, acc, don, frame.Step);
            // std:: cout << acc << " - " << don << " : " << prot << " dist: " << d_ha << " ang: " << ang/3.1416*180 << '\n';

            if (group1Set.count(acc) && group1Set.count(don)) {
              ++intra1;
              // std::cout << "INTRA 1" << ' ';
            } else if (group2Set.count(acc) && group2Set.count(don)) {
              ++intra2;
              // std::cout << "INTRA 2" << ' ';
            } else {
              ++inter;
              // std::cout << "INTER" << ' ';
            }
            // std::cout << '\n';

            activeAcceptors.insert(acc);
            activeDonors.insert(don);
            ++hbond;
          }
        }
      }
    }

    PrintStats(fpStats, hbond, intra1, intra2, inter, activeAcceptors, activeDonors, group1Set, group2Set);
    // std::cout << "No. H-Bonds: " << hbond << '\n';
    // std::cout << "No. intra 1 H-Bonds: " << intra1 << '\n';
    // std::cout << "No. intra 2 H-Bonds: " << intra2 << '\n';
    // std::cout << "No. inter H-Bonds: " << inter << '\n';
    // std::cout << "Active acceptors: " << activeAcceptors.size() << '\n';
    // std::cout << "Active donors: "   << activeDonors.size() << '\n';

  }

  auto percentages = computeBondPresencePercentages(bond_presence, frame.Step);

  std::ofstream fpBondStats("bond_percentages.txt");
    for (const auto &entry : percentages) {
        fpBondStats << entry.first.acceptor << " -> " << entry.first.donor
                    << " : " << entry.second << "%\n";
    }
    fpBondStats.close();
  fpTime.close();
  fpStats.close();
}

void Logic::Finish() {
  std::cout << std::endl;
  std::cout << "Last step: " << frame.Step << '\n';

  // Density(frame, "o", nslabs);
  infile.close();
  std::cout << "===========================\n";
  std::cout << "Print information\n";
  std::cout << std::endl;
}
