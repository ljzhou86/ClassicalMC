#include<fstream>
#include<map>
//
//
using namespace std;

class Metropolis {

private:
  int nions;       // number of magnetic ions
  double * state;  // the state of spin system 

  int nneigh_max;  //
  int * p_nneigh;  // pointer to the number of neighbors
  int** pp_neigh;  // pointer for the list of neighbors
  double** pp_J;   // pointer for the list of J

  int steps_eq;    // equilibrium steps 
  int steps_pr;    // production steps

  double energy;   // internal energy per ion
  double M[3];     // average magnetization
  double pair_energy(double,double,double,double,double); // pair contribution to energy
  
public:
  double beta;     // inverse temperature (1/meV)
  ofstream output; // each object has its unique output file
  
  int init();  // initialize from input file
  int run();
  int one_step();

  ~Metropolis();
    
};
