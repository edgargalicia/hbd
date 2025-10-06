#include "hbmap.h"
#include "logic.h"
#include <Math/Matrix.h>
#include <fstream>
#include <iostream>
#include <params.h>
#include <stdexcept>


Logic::Logic() : filename{"input.dat"}, conf{filename} {}

void Logic::Init() {
  std::cout << "===========================\n";
  std::cout << "Read settings of trajectory\n";
  conf.Print();

  Box box = InitBox(conf.Box());

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

  HBMap hbmap(conf.GetGroup1(), topo);
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

  std::cout << "\tAssess presence of bonds\n";
  std::cout << "\tGather statistics\n";
  while (frame.Read(infile)) {
    if (frame.Step % 1000 == 0) {
      std::cout << "Frame step: " << frame.Step << '\r' << std::flush;
    }
  }

  std::cout << std::endl;
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
