#include "logic.h"
#include <Math/Matrix.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include <params.h>
#include <stdexcept>


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

  // std::vector<int> grp = String2IntList(conf.GetGroup1());
  // std::cout << "Group 1 size: " << grp.size() << std::endl;
  //
  // grp = String2IntList(conf.GetGroup2());
  // std::cout << "Group 2 size: " << grp.size() << std::endl;

  hbmap = HBMap(conf.GetGroup1(), topo);
  hbmap.Print();
  hbmap = HBMap{conf.GetGroup2(), topo};
  hbmap.Print();
}

void Logic::Run() {
  std::cout << "===========================\n";
  std::cout << "CYCLE:\n";
  infile.open(conf.TrajName());
  if (!infile.is_open()) {
    throw std::runtime_error( "Failed to open trajectory file: " + conf.TrajName() );
  }

  int hbond = 0;
  float d_ha;
  float ang;
  std::cout << "\tAssess presence of bonds\n";
  std::cout << "\tGather statistics\n";
  // while (frame.Read(infile)) {
  for (size_t i = 0; i != 1; ++i) {
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

    for(auto acc : hbmap.getAcceptors()) {
      for(auto &don : hbmap.getDonUno()) {
        for(const auto &prot : don.second) {
          if (isHBonded(don.first, prot, acc, frame.Coords, box, d_ha, ang) ) {
            std:: cout << acc << " - " << don.first << " : " << prot << " dist: " << d_ha << " ang: " << ang/3.1416*180 << '\n';
            ++hbond;
          }
        }
      }
    }

  }

  std::cout << std::endl;
  std::cout << "Intra hbonds: " << hbond << '\n';
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
