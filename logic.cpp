#include "logic.h"
#include <fstream>
#include <iostream>
#include <vector>


Logic::Logic() : filename{"input.dat"}, conf{filename} {}

void Logic::Init() {
  std::cout << "===========================\n";
  std::cout << "Read settings of trajectory\n";
  conf.Print();

  std::cout << "===========================\n";
  std::cout << "Read topology\n";
  topo.Read(conf);
  topo.Print();

  frame.Init(conf.Size());

  std::cout << "Group: " + conf.GetGroup1() << '\n';
  std::vector<int> grp = String2IntList(conf.GetGroup1());

}

void Logic::Run() {
  std::cout << "===========================\n";
  std::cout << "CYCLE:\n";
  infile.open(conf.TrajName());
  if (!infile.is_open()) {
    std::cerr << "Failed to open " << conf.TrajName() << '\n';
  }

  std::cout << "\tAssess presence of bonds\n";
  std::cout << "\tGather statistics\n";
  while (!infile.eof()) {

    frame.Read(infile);

    if (frame.Step % 1000 == 0) {
      std::cout << "Frame step: " << frame.Step << '\r' << std::flush;
    }
  }

}

void Logic::Finish() {
  std::cout << std::endl;
  std::cout << "Last step: " << frame.Step << '\n';

  // Density(frame, "o", nslabs);
  infile.close();
  std::cout << "===========================\n";
  std::cout << "Print information\n";
}
