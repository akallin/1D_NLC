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
  

int main(int argc, char** argv){

  int CurrentArg = 1;
  string InputFile;
  string OutputFile = "output_2d.dat";
  bool LF = false;
  // flags to set the input file (need to do that), output file (not used), and low field
  while (CurrentArg < argc)
    {
      if (argv[ CurrentArg ] == string("-i") || argv[ CurrentArg ] == string("--input"))
        {
	  InputFile = string(argv[ CurrentArg + 1 ]);
        }
      if (argv[ CurrentArg ] == string("-o") || argv[ CurrentArg ] == string("--output"))
        {
	  OutputFile = string(argv[ CurrentArg + 1 ]);
        }
      if (argv[ CurrentArg ] == string("-LF") || argv[ CurrentArg ] == string("--lowfield"))
	{
	  LF = true;
	}
      CurrentArg++;
    }


    double energy;

    PARAMS prm;  //Read parameters from param.dat  : see simparam.h
    double J;
    double h;

    Array<l_double,1> eVec;
    Array<long double,1> entVec;

    J=prm.JJ_;
    h=prm.hh_;


    vector< graph > fileGraphs; //graph objects
    
    vector<double> WeightEnergy;
    vector<double> WeightEntropy, WeightMagnetization;
    double RunningSumEnergy(0);
    double RunningSumEntropy(0), RunningSumMagnetization;
    double mag;

    ReadGraphsFromFile(fileGraphs, InputFile);

    ofstream fout(OutputFile.c_str());
    fout.precision(10);
    cout.precision(10);
    
    J=1;     

    const int numhVals = 1;
    //double hvals[numhVals] = {0.2,0.5,1,2,2.5,3.04,3.043,3.044,3.04405,3.0441,3.04415,3.0442,3.045,3.05,4,6,8,10,20,2000};
    double hvals[numhVals] = {3.0441};
    //    double alphas[6] = {0.5,0.75,1.0,1.5,2.0,2.5};
    double alphas[1]={2};
    double alpha = 2.0;
    for(int q1=0;q1<1;q1++){ alpha = alphas[q1];
      //cout << "-----------------------------------" << endl << "S_"<<alpha<<endl;
      
    for(int hh=0; hh<numhVals; hh++){
      //h=hvals[hh];  
      h = hvals[hh];
      //cout <<  "h= " <<h <<" " << endl;

      //One Site Graph
      WeightEnergy.push_back(-h); //Energy weight for zero graph (one site)
      WeightEntropy.push_back(0);
      WeightMagnetization.push_back(1.);
      RunningSumEnergy = WeightEnergy.back();      
      RunningSumEntropy = 0;
      RunningSumMagnetization = WeightMagnetization.back();
      
      //Two Site Graph
      //WeightEnergy.push_back(-sqrt(1+4*h*h)+2*h);
      //WeightEntropy.push_back(TwoSiteEntropy(h,alpha));
      //End of 2 site system stuff!!

      //RunningSumEnergy+=WeightEnergy.back();
      //RunningSumEntropy+=WeightEntropy.back();

      for (int i=1; i<fileGraphs.size(); i++){ //skip the zeroth graph
  	
	//---Generate the Hamiltonian---
	//	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,fileGraphs.at(i).LowField); 
	//	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,h<1.0); 
	//      GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,0); 
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,LF); 
	

	LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
	HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	
	//---Diagonalize and get Eigenvector---
	energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	//cout << "Graph " << i <<" energy: " << energy << endl;
	
	//---Energy/Entropy NLC Calculation---
	WeightEnergy.push_back(energy);
	//	Entropy1D(alpha, eVec, entVec, mag);
	Entropy2D(alpha, eVec, entVec, mag, fileGraphs.at(i).RealSpaceCoordinates);
	
	WeightEntropy.push_back(entVec(0));
	WeightMagnetization.push_back(mag);
	//cout<<"Entropy "<<i<<" = "<<WeightEntropy.back()<<endl;

	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightEntropy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEntropy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightMagnetization.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightMagnetization[fileGraphs.at(i).SubgraphList[j].first];
	}

	// cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<" energy "<<setprecision(12)<<energy<<endl;
	// cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumEnergy += WeightEnergy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumEntropy += WeightEntropy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumMagnetization += WeightMagnetization.back()*fileGraphs.at(i).LatticeConstant;
	//	cout <<"S_ " <<alpha <<" h= "<< h<< " RunningSumEnergy " << i <<" "<< RunningSumEnergy << " Entropy= " << RunningSumEntropy 
	//  << " Magnetization= " << RunningSumMagnetization << endl;
      }
      
      cout<<"S_"<<alpha<<" h= " <<h<<" Energy= "<<RunningSumEnergy<<" Entropy= "<<RunningSumEntropy<<" Magnetization= "<<RunningSumMagnetization<<endl;
      
      WeightEnergy.clear();
      WeightEntropy.clear();
      WeightMagnetization.clear();
      RunningSumEnergy=0;
      RunningSumEntropy=0;
      RunningSumMagnetization=0;
    }
    }
    fout.close();
    return 0;
}
