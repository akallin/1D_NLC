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
    Array<l_double,1> entVec;

    J=prm.JJ_;
    h=prm.hh_;


    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightEnergy;
    vector<double> WeightEntropy;
    double RunningSumEnergy(0);
    double RunningSumEntropy(0);

    ReadGraphsFromFile(fileGraphs, "lineargraphs.dat");
    //ReadGraphsFromFile(fileGraphs,"L16pbc.dat");

    ofstream fout("output_1D.dat");
    fout.precision(10);
    cout.precision(10);
    
    J=1;           
    double hvals[24] = {.2,.3,.4,.5,.6,.7,.75,.875,.9375,.96875,1,1.03125,1.0625,1.125,1.25,1.5,1.75,2.0,2.5,3,3.5,4,4.5,5};
    
    for(int hh=0; hh<24; hh++){
      h=hvals[hh];  
      cout <<  "h= " <<h <<" ";

      //One Site Graph
      WeightEnergy.push_back(-h); //Energy weight for zero graph (one site)
      WeightEntropy.push_back(0);
      double RunningSumEnergy = WeightEnergy[0];      
      double RunningSumEntropy = 0;

      //Two Site Graph
      WeightEnergy.push_back(-sqrt(1+4*h*h)+2*h);
      WeightEntropy.push_back(-log(1+8*h*h)+log(2+8*h*h));
      //End of 2 site system stuff!!

      RunningSumEnergy+=WeightEnergy.back();
      RunningSumEntropy+=WeightEntropy.back();

      for (int i=2; i<fileGraphs.size()-5; i++){ //skip the zeroth graph
  	
	//---Generate the Hamiltonian---
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,fileGraphs.at(i).LowField); 
	
	LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
	HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	
	//---Diagonalize and get Eigenvector---
	energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	//cout << "Graph " << i <<" energy: " << energy << endl;
	
	//---Energy/Entropy NLC Calculation---
	WeightEnergy.push_back(energy);
	Entropy1D(eVec, entVec);
	WeightEntropy.push_back(entVec(1));
	//cout<<"Entropy "<<i<<" = "<<WeightEntropy.back()<<endl;

	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightEntropy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEntropy[fileGraphs.at(i).SubgraphList[j].first];
	}

	// cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<" energy "<<setprecision(12)<<energy<<endl;
	// cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumEnergy += WeightEnergy.back();
	RunningSumEntropy += WeightEntropy.back();
	//cout << "RunningSumEnergy " << i << " = "<< RunningSumEnergy <<endl;
      }
      
      cout<<" Energy= "<<RunningSumEnergy<<" Entropy= "<<RunningSumEntropy<<endl;
      
      WeightEnergy.clear();
      WeightEntropy.clear();
      RunningSumEnergy=0;
      RunningSumEntropy=0;
    }
    
    fout.close();
    return 0;
}
