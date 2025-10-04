void read_acceptors( const std::string &str, HBData &hb ) {
  std::stringstream ss( str );
  std::vector<int> arr;
  std::string token;

  while (ss >> token) {
    if (token.back( ) == ',') token.pop_back( );
    if (token[ 0 ] == '-') {
      ss >> token;
      if (token.back( ) == ',') token.pop_back( );
      for (int i = arr.back( ) + 1; i <= std::stoi( token ); i++) {
        arr.push_back( i );
      }
    } else {
      arr.push_back( std::stoi( token ) );
    }
  }
  hb.nacc = static_cast<int>( arr.size( ) );
  hb.a = new int[ hb.nacc ];
  for (int i = 0; i < hb.nacc; i++) {
    hb.a[ i ] = arr[ i ];
  }
  arr.clear( );
  std::cout << "Number of acceptors: " << hb.nacc << std::endl;
}

void search_donors( const std::string &str, const Box &box, rvec *x[],
                    HBData &hb ) {
  std::stringstream ss( str );
  std::string token;

  float rah2;
  float rc2 = dOH * dOH;
  rvec r_ah;
  std::vector<int> arr;
  std::vector<int> vprot;
  std::vector<int> vnh;

  while (ss >> token) {
    if (token.back( ) == ',') token.pop_back( );
    if (token[ 0 ] == '-') {
      ss >> token;
      if (token.back( ) == ',') token.pop_back( );
      for (int i = arr.back( ) + 1; i <= std::stoi( token ); i++) {
        arr.push_back( i );
      }
    } else {
      arr.push_back( std::stoi( token ) );
    }
  }

  const int nprot = static_cast<int>( arr.size( ) );
  std::vector<int> hh;
  std::copy( arr.begin( ), arr.end( ), std::back_inserter( hh ) );
  arr.clear( );

  /* Search for donors */
  for (int ia = 0; ia < hb.nacc; ++ia) {
    int nppd = 0;
    for (int ih = 0; ih < nprot; ++ih) {
      rvec_sub( *x[ hb.a[ ia ] ], *x[ hh[ ih ] ], r_ah );
      pbc( &r_ah, box );

      rah2 = norm2( r_ah );
      if (rah2 < rc2) {
        vprot.push_back( hh[ ih ] );
        nppd++;
      }
    }
    if (nppd > 0) {
      arr.push_back( hb.a[ ia ] );
      vnh.push_back( nppd );
    }
  }

  hb.ndon = static_cast<int>( arr.size( ) );
  hb.d = new Donor[ hb.ndon ];
  for (int i = 0, k = 0; i < hb.ndon; ++i) {
    hb.d[ i ].don = arr[ i ];
    hb.d[ i ].nh = vnh[ i ];
    hb.d[ i ].hh = new int[ vnh[ i ] ];
    for (int j = 0; j < vnh[ i ]; ++j) {
      hb.d[ i ].hh[ j ] = vprot[ k++ ];
    }
  }

  arr.clear( );
  vprot.clear( );
  vnh.clear( );

  std::cout << "Number of donors: " << hb.ndon << std::endl;
  std::cout << "Number of protons: " << nprot << std::endl << std::endl;
}

