#include<iostream>
#include<fstream>
#include<string>
#include<cstdlib>
#include<cmath>
#include<ctime>

#include"Metropolis.hpp"

using namespace std;

Metropolis::~Metropolis() {
  delete [] state;
  output.close();
}

//===================================
// initialize from input file
int Metropolis::init() {

  // initialize the random generator
  srand (time(NULL));
  
  ifstream input;
  input.open("input");

  string line;
  string::size_type sz;
  
  // MC steps
  getline(input, line);
  steps_eq = stoi (line);  // equilibrium steps
  getline(input, line);
  steps_pr = stoi (line);  // production steps
  
  // number of spins per unit cell
  int n;
  getline(input, line);
  n = stoi(line, &sz);
  output << "# spin per unit cell= " << n << endl;

  // size of super cell
  int n1, n2, n3;
  getline(input, line);
  n1 = stoi(line, &sz);
  line = line.substr(sz);
  n2 = stoi(line, &sz);
  line = line.substr(sz);
  n3 = stoi(line, &sz);
  output << "# size of super cell= " << n1 << ' ' << n2 << ' ' << n3 << endl;

  // allocate space for state
  // each spin is described by (x,y,z)
  nions = n * n1 * n2 * n3 ;
  state = new double [ nions*3 ];

  M[0]=0.0; M[1]=0.0; M[2]=0.0;
  for( int i=0; i<nions; i++ ) {
    // initial magnetization along z-direction
    state[ i*3+0 ] = 0.0; //(float) rand() / RAND_MAX ;  
    state[ i*3+1 ] = 0.0; //(float) rand() / RAND_MAX ;
    state[ i*3+2 ] = 1.0; //(float) rand() / RAND_MAX ;
    M[0] += state [ i*3+0 ] / nions;
    M[1] += state [ i*3+1 ] / nions;
    M[2] += state [ i*3+2 ] / nions;
  }
  output << "# initial magnetization: ";
  output << M[0] << ' ';
  output << M[1] << ' ';
  output << M[2] << endl;

  // setup the hamiltonian
  pp_neigh = new int* [nions];
  pp_J     = new double* [nions];
  p_nneigh = new int [nions];
  nneigh_max = 20;
  for(int i=0; i<nions; i++) {
    pp_neigh[i] = new int    [nneigh_max];
    pp_J[i]     = new double [nneigh_max];
    p_nneigh[i] = 0;  // initially zero neighbors for each site
  }

  int nh; // number of hamiltonian terms in input file
  getline(input, line);
  nh = stoi ( line );

  int off1, off2, off3, ion1, ion2, index1, index2;
  double coupling;
  int i1,j1,k1,i2,j2,k2;
  for(int ih=0; ih < nh; ih++ ) {
    //output << "ih,nh " << ih << ' ' << nh << endl;
    getline(input, line);
    off1 = stoi(line, &sz);
    line = line.substr(sz);
    off2 = stoi(line, &sz);
    line = line.substr(sz);
    off3 = stoi(line, &sz);
    line = line.substr(sz);
    ion1 = stoi(line, &sz)-1;
    line = line.substr(sz);
    ion2 = stoi(line, &sz)-1;
    line = line.substr(sz);
    coupling = stof( line );
    
    for(i1=0; i1<n1; i1++ ) {
      for(j1=0; j1<n2; j1++ ) {
	for(k1=0; k1<n3; k1++ ) {
	  i2=(i1+off1) % n1;
	  j2=(j1+off2) % n2;
	  k2=(k1+off3) % n3;
	  i2 = i2<0 ? i2 + n1 : i2;
	  j2 = j2<0 ? j2 + n2 : j2;
	  k2 = k2<0 ? k2 + n3 : k2;
	  index1 = ( i1*n2*n3 + j1*n3 + k1 ) * n + ion1;
	  index2 = ( i2*n2*n3 + j2*n3 + k2 ) * n + ion2;

	  //output << "#i1,j1,k1 " << i1 << ' ' << j1 << ' ' << k1 ;
	  //output << "#i2,j2,k2 " << i2 << ' ' << j2 << ' ' << k2 ;
          //output << "#coupling " << coupling << endl;
	  //output << "#ion1,ion2 " << ion1 << ' ' << ion2 << endl;
	  //output << "# index1,index2 " << index1 << ' ' << index2 << endl;
	  
	  // index1 has a neighbor of index2
	  p_nneigh[index1]+=1;
	  pp_neigh[index1][ p_nneigh[index1]-1 ] = index2;
	  pp_J    [index1][ p_nneigh[index1]-1 ] = coupling;
	  // index2 also has a neighbor of index1
	  //p_nneigh[index2]+=1;
	  //pp_neigh[index2][ p_nneigh[index2]-1 ] = index1;
	  //pp_J    [index2][ p_nneigh[index2]-1 ] = coupling;
	}
      }
    }
  }
  //for(int i=0; i<nions; i++)
  //  output << "# neighbors of site " << i << " : " << p_nneigh[i] << endl;

  // initial energy
  energy = 0.0;
  for( int ion1=0; ion1 < nions; ion1++ ) { // loop over sites
    for ( int i=0; i < p_nneigh[ion1]; i++ ) { // loop over its neighbors
      ion2 = pp_neigh[ion1][i];
      energy += state[ 3*ion1 + 0 ] * state[ 3*ion2 + 0 ] * pp_J[ion1][i] / 2 ;
      energy += state[ 3*ion1 + 1 ] * state[ 3*ion2 + 1 ] * pp_J[ion1][i] / 2 ;
      energy += state[ 3*ion1 + 2 ] * state[ 3*ion2 + 2 ] * pp_J[ion1][i] / 2 ;
    }
  }
  output << "# initial energy= " << energy << endl;
  return 1;
}

//=============================================================================
// choose randomly an ion and try to update its direction
int Metropolis::one_step() {
  //output << "# beta " << beta << endl;
  int ion_selected;
  ion_selected = (int) ( (float) rand() / RAND_MAX * nions + 1 ) % nions;

  //output << "ion_selected " << ion_selected << endl;

  int ion1, ion2, coupling;
  // a random vector as the new direction of ion_selected
  double x1,x2;
  double new_x, new_y, new_z;
  while (true) {
    x1 = ( (float) rand() / RAND_MAX - 0.5 ) *2.0;
    x2 = ( (float) rand() / RAND_MAX - 0.5 ) *2.0;
    if ( x1*x1 + x2*x2 < 1.0 ) break;
  }
  new_x = 2.0 * x1 * sqrt( 1.0 - x1*x1 - x2*x2 );
  new_y = 2.0 * x2 * sqrt( 1.0 - x1*x1 - x2*x2 );
  new_z = 1.0 - 2.0 * ( 1.0 - x1*x1 - x2*x2 );
  
  // energy related to the selected ion in current state
  double old_E = 0.0;
  double new_E = 0.0;

  for ( int i=0; i < p_nneigh[ion_selected]; i++ ) { 
    ion2 = pp_neigh[ion_selected][i];

    old_E += state[ ion_selected*3+0 ] * state[ ion2*3+0 ] * pp_J[ion_selected][i];
    old_E += state[ ion_selected*3+1 ] * state[ ion2*3+1 ] * pp_J[ion_selected][i];
    old_E += state[ ion_selected*3+2 ] * state[ ion2*3+2 ] * pp_J[ion_selected][i];

    new_E += new_x * state[ ion2*3+0 ] * pp_J[ion_selected][i];
    new_E += new_y * state[ ion2*3+1 ] * pp_J[ion_selected][i];
    new_E += new_z * state[ ion2*3+2 ] * pp_J[ion_selected][i];
  }

  // probability to accept the new direction
  double W;
  W = min ( exp( -beta * ( new_E - old_E ) ), 1.0 );

  // to accept the new direction
  if ( W > (float) rand() / RAND_MAX ) {
    energy += (new_E - old_E);
    M[0] += ( new_x - state[ ion_selected*3+0] ) / nions;
    M[1] += ( new_y - state[ ion_selected*3+1] ) / nions;
    M[2] += ( new_z - state[ ion_selected*3+2] ) / nions;    
    state[ ion_selected*3 + 0 ] = new_x;
    state[ ion_selected*3 + 1 ] = new_y;
    state[ ion_selected*3 + 2 ] = new_z;
  }

  //cout << "old_E,new_E,energy= " << old_E << " " <<new_E << " " << energy<<endl;
}

//====================================================================
int Metropolis::run() {
  for( int i=0; i<steps_eq*nions; i++ )
    one_step();

  int nstat=0;
  double e_avg, e2_avg, m_avg, m2_avg;
   e_avg=0.0;
  e2_avg=0.0;
   m_avg=0.0;
  m2_avg=0.0;
  for( int i=0; i<steps_pr; i++ ) {
    for(int j=0; j<nions; j++ ) {
      one_step();
    }
    
    // statistics
    output << energy<<' ';
    output << M[0] <<' ';
    output << M[1] <<' ';
    output << M[2] <<' '<<endl;

    e_avg += energy;
    e2_avg += energy * energy;
    m_avg += sqrt( M[0]*M[0] + M[1]*M[1] + M[2]*M[2] );
    m2_avg +=      M[0]*M[0] + M[1]*M[1] + M[2]*M[2];
    nstat += 1;
  }
  
  e_avg /= nstat;
  e2_avg /= nstat;
  m_avg /= nstat;
  m2_avg /= nstat;

  output<<"T,beta,<e>,<e2>,<m>,<m2> ";
  output<<1.0/beta/8.6173324e-2<<' ';
  output<<beta<<' ';
  output<<e_avg<<' ';
  output<<e2_avg<<' ';
  output<<m_avg<<' ';
  output<<m2_avg<<endl;
}
