#include "topology.h"
#include "Math/Matrix.h"
#include "Math/Vectors.h"
#include <algorithm>
#include <cctype>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

void Topology::Print() {
  std::cout << "Topology size: " << atomNames.size() << '\n';
  std::cout << "Bonds size: " << bonds.size() << '\n';

  std::unordered_map<std::string, int>  freqAtoms;
  for (const auto &atom : atomNames) {
    ++freqAtoms[atom];
  }

  for (const auto & frqAt : freqAtoms) {
    std::cout << frqAt.first << ": " << frqAt.second << '\n';
  }
}

void Topology::ComputeBonds(const std::vector<Math::Vec3> &coords, const Box &box) {
  ClearBonds(); // You need to implement this
  float rcBond2 = 1.3f * 1.3f;
  Math::Vec3 dist;

  for (size_t i = 0; i != coords.size() - 1; ++i) {
    for (size_t j = i + 1; j != coords.size(); ++j) {
      dist = coords[i] - coords[j];
      box.Pbc(dist);
      if (dist.norm2() < rcBond2) {
        PushBond(i, j);
      }
    }
  }
}

Topology ReadTopology(const Config &config, const Box &box) {
  Topology top;
  std::ifstream fp(config.TrajName());
  if ( !fp.is_open() ) {
    throw std::runtime_error( "Failed to open topology file: " + config.TrajName() );
  }

  std::string line;
  std::getline(fp, line);
  std::getline(fp, line);

  std::string atomname;
  Math::Vec3 coord;
  std::vector<Math::Vec3> coords;
  for (size_t i = 0; i != config.Size(); ++i) {
    std::getline(fp, line);
    std::istringstream iss(line);
    iss >> atomname >> coord;
    std::transform(atomname.begin(), atomname.end(), atomname.begin(),
                   [](unsigned char c) {
                     return std::tolower(c);
                   });
    top.PushAtom(atomname);
    coords.push_back(coord);
  }

  top.ComputeBonds(coords, box);
  fp.close();
  return top;
}

void Topology::UpdateBonds(const Frame &frame, const Box &box) {
  // Store previous bonds
  std::vector<std::pair<int, int>> oldBonds = Bonds();

  // Recompute bonds
  ComputeBonds(frame.Coords(), box);

  auto newBonds = Bonds();
}

void Frame::Init(int natoms) {
  step = 0;
  coords.resize(natoms);
}

bool Frame::Read(std::ifstream &fp) {
  std::string line;
  if ( !std::getline(fp, line) ) {
    return false;
  }
  if ( !std::getline(fp, line) ) {
    throw std::runtime_error("Truncated frame header");
  }

  std::string dummy;
  for (auto &atom : coords) {
    if ( !std::getline(fp, line) ) {
      throw std::runtime_error("Unexpected EOF inside frame");
    }
    std::istringstream iss(line);
    if ( !(iss >> dummy >> atom[0] >> atom[1] >> atom[2]) ) {
      throw std::runtime_error("Malformed line in frame");
    }
  }
  ++step;
  return true;
}

Box InitBox( const Math::Matrix33 &mbox ) {
  Box box;
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      box.mbox[ i ][ j ] = mbox[ i ][ j ];
    }
    box.fbox[ i ] = mbox[i].norm();
    box.hbox[ i ] = 0.5 * box.fbox[ i ];
    box.mhbox[ i ] = -( box.hbox[ i ] );
  }
  return box;
}

void Box::Pbc( Math::Vec3 &dx ) const {
  for (int i = 2; i >= 0; i--) {
    while (( dx )[ i ] > hbox[ i ]) {
      dx -= mbox[i];
    }
    while (( dx )[ i ] <= mhbox[ i ]) {
      dx += mbox[i];
    }
  }
}

void Select( std::vector<int> &list, const Frame &frame, float za, float zb, Axis axis) {
  list.clear();
  for (size_t i = 0; i != frame.Size(); ++i) {

    if (frame.Coords()[i][axis] > za && frame.Coords()[i][axis] < zb) {
      list.push_back(i);
    }
  }
}
