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
#include <sstream>

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
  double alpha = 2.0;
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
      if (argv [ CurrentArg ] == string("-s") || argv[ CurrentArg ] == string("-S"))
	{
	  alpha = atof( argv [ CurrentArg + 1]);
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
    vector<double> WeightLineEntropy, WeightCornerEntropy, WeightMagnetization;
    double RunningSumEnergy(0);
    double RunningSumLineEntropy(0), RunningSumCornerEntropy(0), RunningSumMagnetization;
    double mag;

    ReadGraphsFromFile(fileGraphs, InputFile);

    ofstream fout(OutputFile.c_str());
    fout.precision(10);
    cout.precision(10);
    
    J=1;     


    const int numhVals = 23;
    //23 values
    double hvals[numhVals] = {0.2,0.5,1.0,1.5,2.0,2.5,3.0,3.0441,3.5,4.0,4.5,5.0,5.5,6.0,6.5,7.0,7.5,8.0,8.5,9.0,9.5,10,2000};
    //double hvals[numhVals] = {3.0441};

    
    string magFile;
    double magOne;

    for(int hh=0; hh<numhVals; hh++){
      //h=hvals[hh];  
      h = hvals[hh];
      //cout <<  "h= " <<h <<" " << endl;
      ostringstream s;
      s<<"./MagFRiles/mag"<<h<<".input";
      magFile = s.str();
      s.clear();
      
      ifstream magIn(magFile.c_str());
      if(magIn){
	magIn >> magOne;
      }
      else{ 	magOne=1.0;      }
      magIn.close();

      //One Site Graph
      WeightEnergy.push_back(-h); //Energy weight for zero graph (one site)
      WeightLineEntropy.push_back(0);
      WeightCornerEntropy.push_back(0);
      WeightMagnetization.push_back(0);
      RunningSumEnergy = WeightEnergy.back();      
      RunningSumLineEntropy = 0;
      RunningSumCornerEntropy = 0;
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
	GENHAM HV(fileGraphs.at(i).NumberSites,J,h,fileGraphs.at(i).AdjacencyList,LF,magOne); 
	

	LANCZOS lancz(HV.Vdim);  //dimension of Hilbert space
	HV.SparseHamJQ();  //generates sparse matrix Hamiltonian for Lanczos

	
	//---Diagonalize and get Eigenvector---
	energy = lancz.Diag(HV, 1, prm.valvec_, eVec); // Hamiltonian, # of eigenvalues to converge, 1 for -values only, 2 for vals AND vectors
	//cout << "Graph " << i <<" energy: " << energy << endl;
	
	//---Energy/Entropy NLC Calculation---
	WeightEnergy.push_back(energy);
	//	Entropy1D(alpha, eVec, entVec, mag);
	Entropy2D(alpha, eVec, entVec, fileGraphs.at(i).RealSpaceCoordinates);
	
	WeightLineEntropy.push_back(entVec(1));
	WeightCornerEntropy.push_back(entVec(0));
	WeightMagnetization.push_back(Magnetization(eVec));
	//cout<<"Entropy "<<i<<" = "<<WeightEntropy.back()<<endl;

	for (int j = 0; j<fileGraphs.at(i).SubgraphList.size(); j++){
	  WeightEnergy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightEnergy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightLineEntropy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightLineEntropy[fileGraphs.at(i).SubgraphList[j].first];
	  WeightCornerEntropy.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightCornerEntropy[fileGraphs.at(i).SubgraphList[j].first];	  
	  WeightMagnetization.back() -= fileGraphs.at(i).SubgraphList[j].second * WeightMagnetization[fileGraphs.at(i).SubgraphList[j].first];
	}

	// cout<<"h="<<h<<" J="<<J<<" graph #"<<i<<" energy "<<setprecision(12)<<energy<<endl;
	// cout<<"WeightHigh["<<i<<"] = "<<WeightHigh.back()<<endl;
	RunningSumEnergy += WeightEnergy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumLineEntropy += WeightLineEntropy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumCornerEntropy += WeightCornerEntropy.back()*fileGraphs.at(i).LatticeConstant;
	RunningSumMagnetization += WeightMagnetization.back()*fileGraphs.at(i).LatticeConstant;
	//	cout <<"S_ " <<alpha <<" h= "<< h<< " RunningSumEnergy " << i <<" "<< RunningSumEnergy << " Entropy= " << RunningSumEntropy 
	//  << " Magnetization= " << RunningSumMagnetization << endl;
      }
      
      cout<<"S_"<<setw(4)<< alpha<<" h= " <<setw(6)<<h<<" Energy= "<<setw(15)<<RunningSumEnergy<<" LineEnt= "<<setw(15)<<RunningSumLineEntropy
	  <<" CornerEnt= "<<setw(15)<<RunningSumCornerEntropy<<" Magnetization= "<<setw(15)<<RunningSumMagnetization<<endl;
      
      ofstream magOut(magFile.c_str());
      magOut << RunningSumMagnetization;
      magOut.close();
      

      WeightEnergy.clear();
      WeightLineEntropy.clear();
      WeightCornerEntropy.clear();
      WeightMagnetization.clear();
      RunningSumEnergy=0;
      RunningSumLineEntropy=0;
      RunningSumCornerEntropy=0;
      RunningSumMagnetization=0;
    }

    fout.close();
    return 0;
}
