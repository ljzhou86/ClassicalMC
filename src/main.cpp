#include<fstream>
#include<iostream>
#include<string>
#include<sstream>
#include<vector>

#include"Metropolis.hpp"
#include "mpi.h"

using namespace std;

int main( int argc, char *argv[] ){
  // some constant
  double kb=8.6173324e-2; // meV per K
  
  // initialize MPI environment
  int id;
  int p;
  MPI_Init ( &argc, &argv );
  MPI_Comm_rank ( MPI_COMM_WORLD, &id );
  MPI_Comm_size ( MPI_COMM_WORLD, &p );

  // read temper file
  double Tmax, Tmin;
  int nT;
  double * T;
  ifstream temper;
  temper.open("temper");
  temper >> Tmin >> Tmax >> nT;
  T = new double [nT];
  for(int i=0;i<nT;i++) {
    T[i] = (Tmax-Tmin)/(nT-1)*i + Tmin;
  }
  temper.close();
    
  // create and initialize metropolis objects
  // set temperature for each metropolis object
  Metropolis* metro;
  vector< Metropolis* > metro_list;
  for( int iT=0; iT<nT; iT++ ) {
    if ( iT % p == id ) {
      metro = new Metropolis;
      ostringstream ss;
      ss << "output." << iT;
      metro->output.open(ss.str());
      metro->init();
      metro->beta = 1.0/ ( kb * T[iT] ); // inverse temperature in meV^-1
      metro->output << "# temperature T= ";
      metro->output << T[iT] << " (Kelvin) ";
      metro->output << "beta= " << metro->beta << " (1/meV) "<< endl;;
      metro_list.push_back(metro);
    }
  }

  //metro.one_step();

  for( int i=0; i<metro_list.size(); i++ )
    metro_list[i]->run();

  // MPI finalize
  MPI_Finalize();
}
