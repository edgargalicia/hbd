#include "logic.h"
#include <Math/Matrix.h>
#include <chrono>
#include <fstream>
#include <iomanip>
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

  size_t totalLines = 0;
  std::string tmp;
  std::ifstream infile_count(conf.TrajName());
  int dotCounter = 0;
  while (std::getline(infile_count, tmp)) {
    ++totalLines;

    if (totalLines % 500000 == 0) {
      std::cout << "\rEstimating size of the file ";
      for (size_t i = 0; i != dotCounter % 3 + 1; ++i) std::cout << ".";
      std::cout << "   " << std::flush;
      ++dotCounter;
    }
  }
  infile_count.close();
  std::cout << std::endl;

  size_t linesPerFrame = conf.Size()+2;
  size_t totalFrames = totalLines / linesPerFrame;

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

  auto group1 = String2IntList( conf.GetGroup1() );
  auto group2 = String2IntList( conf.GetGroup2() );
  std::unordered_set<int> group1Set(group1.begin(), group1.end());
  std::unordered_set<int> group2Set(group2.begin(), group2.end());

  std::ofstream fpStats("stats.txt");
  fpStats << std::left
    << std::setw(12) << "# HBs"
    << std::setw(12) << "Intra HBs1"
    << std::setw(12) << "Intra HBs2"
    << std::setw(12) << "Inter HBs"
    << std::setw(12) << "Active acc"
    << std::setw(12) << "Active don"
    << std::setw(12) << "Group1 don"
    << std::setw(12) << "Group1 acc"
    << std::setw(12) << "Group2 don"
    << std::setw(12) << "Group2 acc"
    << '\n';

  std::unordered_map<BondKey, std::vector<int>, BondKeyHash> bond_presence;
  int total_hbond = 0;
  int total_intra1 = 0;
  int total_intra2 = 0;
  int total_inter = 0;
  int totalActiveAcceptors = 0;
  int totalActiveDonors = 0;
  int totalGroup1Donors = 0;
  int totalGroup1Acceptors = 0;
  int totalGroup2Donors = 0;
  int totalGroup2Acceptors = 0;
  std::unordered_map<BondKey, BondType, BondKeyHash> bondTypeMap;

  size_t processedFrames = 0;
  auto startTime = std::chrono::steady_clock::now();

  while (frame.Read(infile)) {
    ++processedFrames;
  // for (size_t i = 0; i != 10; ++i) {
  //   frame.Read(infile);
    hbond = 0;
    intra1 = 0;
    intra2 = 0;
    inter = 0;
    activeAcceptors.clear();
    activeDonors.clear();
    
    Select(sel, frame);
    hbmap = HBMap(sel, topo);

    for(auto acc : hbmap.getAcceptors()) {
      for(const auto &donorKV : hbmap.getDonUno()) {
        const int &don = donorKV.first;
        for(const auto &prot : donorKV.second) {
          if ( isHBonded(don, prot, acc, frame.Coords, box, d_ha, ang)  ) {
            addBondPresence(bond_presence, bondTypeMap, acc, don, frame.Step, group1Set, group2Set);

            if (group1Set.count(acc) && group1Set.count(don)) {
              ++intra1;
              total_intra1++;
            } else if (group2Set.count(acc) && group2Set.count(don)) {
              ++intra2;
              total_intra2++;
            } else {
              ++inter;
              ++total_inter;
            }

            activeAcceptors.insert(acc);
            activeDonors.insert(don);
            ++hbond;
            total_hbond++;
          }
        }
      }
    }

    auto g = CountGroupActivity(activeDonors, activeAcceptors, group1Set, group2Set);
    totalGroup1Donors += g.g1Don;
    totalGroup1Acceptors += g.g1Acc;
    totalGroup2Donors += g.g2Don;
    totalGroup2Acceptors += g.g2Acc;
    totalActiveAcceptors += activeAcceptors.size();
    totalActiveDonors += activeDonors.size();

    PrintStats(fpStats, hbond, intra1, intra2, inter, activeAcceptors, activeDonors, group1Set, group2Set);

    if (frame.Step % 100 == 0) {
      // std::cout << "Frame step: " << frame.Step << '\r' << std::flush;

      auto now = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = now - startTime;
      double percent = 100.0 * processedFrames / totalFrames;
      double estimatedTotalTime = elapsed.count() / (processedFrames / static_cast<double>(totalFrames));
      double etaSeconds = estimatedTotalTime - elapsed.count();
      int etaMin = static_cast<int>(etaSeconds) / 60;
      int etaSec = static_cast<int>(etaSeconds) % 60;

      std::cout << "\rProcessed " << processedFrames << '/' << totalFrames
                << " (" << std::fixed << std::setprecision(1) << percent << "%)"
                << ", ETA: " << etaMin << "m " << etaSec << "s "
                << std::flush;
    }
  }
  std::cout << std::endl;

  auto percentages = computeBondPresencePercentages(bond_presence, frame.Step);
  std::ofstream fpTimeSeries("bond_timeseries.txt");
  fpTimeSeries << std::left
               << std::setw(12) << "# Frame"
               << std::setw(12) << "BondID"
               << '\n';

  std::unordered_map<BondKey, int, BondKeyHash> bondIDMap;
  int bondID = 0;
  for (const auto &entry : percentages) {
    bondIDMap[entry.first] = ++bondID;
  }

  for (const auto &entry : bond_presence) {
    const BondKey &bond = entry.first;
    auto it = bondIDMap.find(bond);
    if (it == bondIDMap.end()) continue;

    int id = it->second;
    for (int frameNum : entry.second) {
      fpTimeSeries << std::setw(12) << frameNum
                   << std::setw(12) << id
                   << '\n';
    }
    fpTimeSeries << "\n\n";
  }

  fpTimeSeries.close();

  std::ofstream fpBondStats("bond_percentages.txt");
  fpBondStats << std::left
              << std::setw(6) << "# ID"
              << std::setw(12) << "Acceptor"
              << std::setw(12) << "Donor"
              << std::setw(12) << "Type"
              << std::setw(15) << "% of Frames"
              << '\n';

  int id = 0;
  for (const auto &entry : percentages) {
    BondType bType = bondTypeMap[entry.first];
    std::string typeStr;
    switch (bType) {
      case BondType::Intra1: typeStr = "Intra1"; break;
      case BondType::Intra2: typeStr = "Intra2"; break;
      case BondType::Inter: typeStr = "Inter"; break;
    }

    fpBondStats << std::right << std::setw(4) << ++id << "  "
                << std::left<< std::setw(12) << entry.first.acceptor
                << std::setw(12) << entry.first.donor
                << std::setw(12) << typeStr
                << std::setw(15) << std::fixed << std::setprecision(2)
                << entry.second
                << '\n';
  }
  // fpBondStats << std::left
  //             << std::setw(12) << "Hbonds"
  fpBondStats << "\n\n";
  fpBondStats << "# Hbonds: " << total_hbond / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "# Intra1: " << total_intra1 / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "# Intra2: " << total_intra2 / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "# Inter: " << total_inter / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "# Active Acceptors: " << totalActiveAcceptors / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "# Active Donors: " << totalActiveDonors / static_cast<double>(frame.Step) << '\n';
  fpBondStats << '\n';
  fpBondStats << "# Average active donors per group per frame:\n";
  fpBondStats << "#   Group1 acceptors: " << totalGroup1Acceptors / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "#   Group1 donors: " << totalGroup1Donors / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "#   Group2 acceptors: " << totalGroup2Acceptors / static_cast<double>(frame.Step) << '\n';
  fpBondStats << "#   Group2 donors: " << totalGroup2Donors / static_cast<double>(frame.Step) << '\n';
  fpBondStats.close();
  fpStats.close();
}

void Logic::Finish() {
  std::cout << "Last step: " << frame.Step << '\n';

  infile.close();
  std::cout << "===========================\n";
  std::cout << "Print information\n";
}
