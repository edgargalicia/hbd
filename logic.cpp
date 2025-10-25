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

void Logic::Init() {
  std::cout << "===========================\n";
  std::cout << "Read settings of trajectory\n";
  conf.Print();

  box = InitBox(conf.Box());

  std::cout << "===========================\n";
  std::cout << "Read topology\n";
  topo = ReadTopology(conf, box);
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

  // Skip frames before Begin()
  if (conf.Begin() != 0) {
    std::cout << "Skipping " << conf.Begin() << " frames...\n";
    for (int i = 0; i < conf.Begin(); ++i) {
      if (!frame.Read(infile)) {
        throw std::runtime_error("Not enough frames in file to skip to Begin()");
      }
    }
  }

  if (conf.End() == -1) {
    size_t totalLines = 0;
    std::string tmp;
    std::ifstream infile_count(conf.TrajName());
    size_t dotCounter = 0;
    int lineString = 500000;
    while (std::getline(infile_count, tmp)) {
      ++totalLines;

      if (totalLines % lineString == 0) {
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
    conf.SeteFrames(totalFrames);
  }

  HBGeometry hbg;
  std::cout << "\tAssess presence of bonds\n";
  std::cout << "\tGather statistics\n";
  std::vector<int> sel;

  auto group1 = String2IntList( conf.GetGroup1() );
  auto group2 = String2IntList( conf.GetGroup2() );
  std::unordered_set<int> group1Set{std::unordered_set<int> (group1.begin(), group1.end())};
  std::unordered_set<int> group2Set{std::unordered_set<int> (group2.begin(), group2.end())};

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

  std::unordered_map<BondPair, std::vector<int>, BondKeyHash> bond_presence;
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
  std::unordered_map<BondPair, BondType, BondKeyHash> bondTypeMap;

  HBondStats hbs;
  size_t processedFrames = 0;
  auto startTime = std::chrono::steady_clock::now();

  HBBinner hbbin;

  // while (frame.Read(infile)) {
  for (int iframe = conf.Begin(); iframe != conf.End(); ++iframe) {
    frame.Read(infile);
    ++processedFrames;
    topo.UpdateBonds(frame, box);
    hbs.reset();
    Select(sel, frame, conf.getZa(), conf.getZb(), conf.GetAxis());
    hbmap = HBMap(sel, topo);

    for(auto acc : hbmap.getAcceptors()) {
      for(const auto &donorKV : hbmap.getDonUno()) {
        const int &don = donorKV.first;
        for(const auto &prot : donorKV.second) {
          BondAtoms atoms{{acc, don}, prot};
          HBondContext ctx{frame.Coords(), box};
          if ( isHBonded(atoms, ctx, hbg) ) {
            hbbin.add(hbg);
            addBondPresence(bond_presence, bondTypeMap, atoms, frame.Step(), group1Set, group2Set);

            if (group1Set.count(acc) && group1Set.count(don)) {
              ++hbs.intra1;
              total_intra1++;
            } else if (group2Set.count(acc) && group2Set.count(don)) {
              ++hbs.intra2;
              total_intra2++;
            } else {
              ++hbs.inter;
              ++total_inter;
            }

            hbs.activeAcceptors.insert(acc);
            hbs.activeDonors.insert(don);
            ++hbs.hbond;
            total_hbond++;
          }
        }
      }
    }

    auto g = CountGroupActivity(hbs.activeDonors, hbs.activeAcceptors, group1Set, group2Set);
    totalGroup1Donors += g.g1Don;
    totalGroup1Acceptors += g.g1Acc;
    totalGroup2Donors += g.g2Don;
    totalGroup2Acceptors += g.g2Acc;
    totalActiveAcceptors += hbs.activeAcceptors.size();
    totalActiveDonors += hbs.activeDonors.size();

    PrintStats(fpStats,hbs,group1Set,group2Set);

    if (frame.Step() % 100 == 0) {
      // std::cout << "Frame step: " << frame.Step << '\r' << std::flush;

      auto now = std::chrono::steady_clock::now();
      std::chrono::duration<double> elapsed = now - startTime;
      double percent = 100.0 * processedFrames / ( conf.End() - conf.Begin() );
      double estimatedTotalTime = elapsed.count() / (processedFrames / static_cast<double>(conf.End() - conf.Begin()));
      double etaSeconds = estimatedTotalTime - elapsed.count();
      int etaMin = static_cast<int>(etaSeconds) / 60;
      int etaSec = static_cast<int>(etaSeconds) % 60;

      std::cout << "\rProcessed " << processedFrames << '/' << conf.End() - conf.Begin()
                << " (" << std::fixed << std::setprecision(1) << percent << "%)"
                << ", ETA: " << etaMin << "m " << etaSec << "s "
                << std::flush;
    }
  }
  std::cout << std::endl;

  auto percentages = computeBondPresencePercentages(bond_presence, frame.Step());
  std::ofstream fpTimeSeries("bond_timeseries.txt");
  fpTimeSeries << std::left
               << std::setw(12) << "# Frame"
               << std::setw(12) << "BondID"
               << '\n';

  std::unordered_map<BondPair, int, BondKeyHash> bondIDMap;
  int bondID = 0;
  for (const auto &entry : percentages) {
    bondIDMap[entry.first] = ++bondID;
  }

  for (const auto &entry : bond_presence) {
    const BondPair &bond = entry.first;
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
  fpBondStats << "# Hbonds: " << total_hbond / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "# Intra1: " << total_intra1 / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "# Intra2: " << total_intra2 / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "# Inter: " << total_inter / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "# Active Acceptors: " << totalActiveAcceptors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "# Active Donors: " << totalActiveDonors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << '\n';
  fpBondStats << "# Average active donors per group per frame:\n";
  fpBondStats << "#   Group1 acceptors: " << totalGroup1Acceptors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "#   Group1 donors: " << totalGroup1Donors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "#   Group2 acceptors: " << totalGroup2Acceptors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats << "#   Group2 donors: " << totalGroup2Donors / static_cast<double>(frame.Step()) << '\n';
  fpBondStats.close();
  fpStats.close();
  std::ofstream fpBin("bins.txt");
  hbbin.print(fpBin);
  fpBin.close();
}

void Logic::Finish() {
  std::cout << "Last Step: " << frame.Step() << '\n';

  infile.close();
  std::cout << "===========================\n";
  std::cout << "Print information\n";
}
