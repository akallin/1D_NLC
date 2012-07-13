/*******************************************************************
Linked Cluster Expansion program for the 1D TFIM model
Roger Melko, Ann Kallin, Katie Hyatt June 2012
baesd on a Lanczos code from November 2007
********************************************************************/

#include <utility>
#include <fstream> 
#include <vector> 
#include <math.h>
using namespace std;

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <iomanip>
#include <blitz/array.h>

BZ_USING_NAMESPACE(blitz)

#include "Lanczos_07.h"
#include "GenHam.h"
#include "lapack.h"
#include "simparam.h"
#include "graphs.h"
#include "entropy.h"
  

int main(){


    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    Array<l_double,1> eVec;

    J=prm.JJ_;
    h=prm.hh_;


    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightHigh;

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");
    //ReadGraphsFromFile(fileGraphs,"L16pbc.dat");

    ofstream fout("output_1D.dat");
    fout.precision(10);
    cout.precision(10);
    
    J=1;           
    double hvals[24] = {.2,.3,.4,.5,.6,.7,.75,.875,.9375,.96875,1,1.03125,1.0625,1.125,1.25,1.5,1.75,2.0,2.5,3,3.5,4,4.5,5};
    
    for(int hh=20; hh<24; hh++){
      h=hvals[hh];  
      cout <<  h <<" ";
      WeightHigh.push_back(-h); //Weight for site zero
      double RunningSumHigh = WeightHigh[0];      
      
      for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph
  	
	//---Generate the Hamiltonian---
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,fileGraphs.at(i).LowField); 
	
	
	LANCZOS lancz(HV.Vdim);  //dimension of reduced Hilbert space (Sz sector)
	HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	//---Diagonalize and get Eigenvector---
	energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	
	//---Entropy Calculation---
	Entropy1D(eVec);
	
	//---Energy NLC Calculation---
	WeightHigh.push_back(energy);
	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++)
	  WeightHigh.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightHigh[fileGraphs.at(i).SubgraphList[j].first];
	
	// cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<" energy "<<setprecision(12)<<energy<<endl;
	// cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumHigh += WeightHigh.back();
	// cout <<"RunningSumHigh = "<< RunningSumHigh <<endl;
      }
      
      fout<<"h= "<<h<<" J= "<<J<<" Energy= "<<RunningSumHigh<<endl;
      
      WeightHigh.clear();
      RunningSumHigh=0;
    }
    
    fout.close();
    return 0;
}
