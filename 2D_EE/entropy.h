//function to calculate the RDM and entropies for a 1d 
#ifndef entropy_H
#define entropy_H

inline double TwoSiteEntropy(double h, double alpha)
{
  double CommonEnt;
  double DiffEnt;
  double unLog;
  CommonEnt = 0.5 + (1. + sqrt(1. + 4.*h*h))/(8.*h*h);
  DiffEnt = h*sqrt(1.+2.*h*h+sqrt(1.+4.*h*h))/2./sqrt(2.0);
  
  unLog = pow(CommonEnt + DiffEnt,alpha) + pow(abs(CommonEnt - DiffEnt),alpha);
  
  // cout << CommonEnt << "  " << DiffEnt << endl;
  if(alpha==1.0){
    // cout << CommonEnt - DiffEnt << endl;
    return -(CommonEnt + DiffEnt)*log(CommonEnt + DiffEnt) 
      - abs(CommonEnt - DiffEnt)*log(abs(CommonEnt - DiffEnt));
  }
  else{
    return (1./(1.-alpha))*log(unLog);
  }
}

inline void Entropy1D(double alpha, Array<l_double,1>& eigs, Array<l_double,1>& ents, double& mag)
{
  // The dimension is number of eigenvalues
  long int Dim = eigs.size();

  // Get number of sites from the dimension
  int Nsite = log2(Dim); //cout << "Nsite = " << Nsite << endl;
 
  // The starting dimensions of region A and region B
  int Adim=1;
  int Bdim=Dim;

  // A rectangular matrix containing the eigenvalues, used to get the RDM
  Array<long double,2> SuperMat;

  // The RDM (SuperMat squared) and RDM squared
  Array<double,2> DM(Adim,Adim);
  Array<double,2> DMsq(Adim,Adim);

  // Some temp variables
  long double temp2(0), temp3(0), temp4(0), temp5(0), temp6(0);
  int a(0),b(0);

  // Eigenvalues of the RDM get put in dd
  vector<double> dd;
  long double vN;

  //Hardcoded to contain only 2 ents (though I may only do *one*!)
  ents.resize(2,0);
  ents = 0;

  // Initializing some things
  long double magnetization(0);
  long double renyi(0);
  l_double norm(0);

  // measure the magnetization
  for(int i=0; i<Dim; i++){ 
    for (int sp=0; sp<Nsite; sp++){
      temp3 += (i>>sp)&1; 
    }
    magnetization += abs(temp3*2-Nsite)*eigs(i)*eigs(i);
    norm += eigs(i)*eigs(i);
    temp3=0;
  }
  // mag was passed by ref to the function, so this is what it returns
  mag = magnetization;
  magnetization = 0;

  // Measure Renyi up to half the graph size
  for(int Asite=1; Asite<(Nsite+2)/2; Asite++){
    // Dimension of region A doubles (from last one) and B gets cut in half
    Adim*=2;
    Bdim/=2;

    // sqrt of RDM gets the properdimensions and initialized
    SuperMat.resize(Adim,Bdim); 
    SuperMat=0;
    
    //More initialization
    a=0;b=0;
    temp3=0;
    norm=0;
    magnetization=0;
    
    for(int i=0; i<Dim; i++){ 
      // extractifying the region A and region B states
      //b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
      //c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
      a = i&(Adim-1);
      b = (i>>Asite)&(Bdim-1);
      
      SuperMat(a,b) = eigs(i);

    
    }
    norm = 0;
    
    DM.resize(Adim,Adim);
    DM=0;
    temp2=0;

    //multiplying the supermat by its transpose to get the RDM
    DM=0;
    for(int i=0; i<Adim; i++){
      for(int j=0; j<Adim; j++){
	temp2=0;
	for(int k=0; k<Bdim; k++){
	  temp2 += SuperMat(i,k)*SuperMat(j,k);
	}
	DM(i,j) = temp2;
      }
    }
    
    //Diagonalizing the RDM
    while(dd.size()>0){dd.erase(dd.begin());}
    diagWithLapack_R(DM,dd); 
    renyi=0; vN=0; temp5=0;
    
    for(int s=0;s<dd.size();s++){
      if(dd[s]<0){dd[s]=0;}
      temp5=log(dd[s]);
      if(!(temp5>-1000000000)){temp5=0;}
      vN+=-dd[s]*temp5;
      if(abs(dd[s])<1e-15){dd[s]=0;}
      renyi+=pow(dd[s],alpha);
    }    

    if(alpha==1.0){temp6 = vN;}
    else{temp6 = 1./(1.-alpha)*log(renyi);}
    ents(1)+=temp6;
    if(Asite<(Nsite+1)/2){ ents(1)+=temp6;}

    //if(Asite<(Nsite+1)/2){ ents(0)+=vN; ents(1)+=-log(renyi);}// ents(1)+=renyi; }
    //else{cout << "Asite:"<<Asite<<" Nsite:"<<Nsite<<endl;}

    //    cout <<" "<< setprecision(15) << -log(renyi) <<" "<< setprecision(15) << vN;    
    
    /*------Commenting out the other S_2 calculation because it's unnecessary---------
    DM.resize(Adim,Adim);
    DM=0;
    temp2=0;
    //multiplying the supermat by its transpose to get the RDM
    DM=0;
    for(int i=0; i<Adim; i++){
      for(int j=0; j<Adim; j++){
	temp2=0;
	for(int k=0; k<Bdim; k++){
	  temp2 += SuperMat(i,k)*SuperMat(j,k);
	}
	DM(i,j) = temp2;
      }
    }

    //Square the DM
    DMsq.resize(Adim,Adim);
    DMsq=0;
    for(int i=0; i<Adim; i++){
      for(int j=0; j<Adim; j++){
	temp2=0;
	for(int k=0; k<Adim; k++){
	  temp2 += DM(i,k)*DM(j,k);
	}
	DMsq(i,j) = temp2; 
      }
    }

    //cout << DMsq << endl;
    
    renyi=0;
    norm=0;
    for(int s=0;s<Adim;s++){
      norm += DM(s,s);
      // cout << "Norm: " << norm << "   ";
      renyi += DMsq(s,s);
      // cout << "Renyi: "<< renyi << endl;
    }
    ---------------------------------------------------------------------*/

    //  cout << "Norm"  << "     " << setprecision(15) << norm << endl;
    //cout << "Asites = " << Asite << "   Renyi"  << "  " << 
    //    cout <<" "<< setprecision(15) << -log(renyi/norm) << endl;   
    //  cout << endl;

  }
  //cout << endl;
}

inline void Entropy2D(double alpha, Array<l_double,1>& eigs, Array<l_double,1>& ents, double& mag, int xMax, int yMax)
{
  // The dimension is number of eigenvalues
  long int Dim = eigs.size();

  // Get number of sites from the dimension
  int Nsite = log2(Dim); 

  // The dimensions of region A
  int xSize(0), ySize(0), Adim(0), Bdim(0);

  // A rectangular matrix containing the eigenvalues, used to get the RDM
  Array<long double,2> SuperMat;

  // The RDM (SuperMat squared) and RDM squared
  Array<double,2> DM(Adim,Adim);
  Array<double,2> DMsq(Adim,Adim);

  // ------ Line Terms!! ------
  // ---- Horizontal ----
  xSize = xMax;
  // Iterate over the horizontal cuts
  for(int ySize=1; ySize<yMax; ySize++){
    // Get the dimensions of region A and B;
    Adim = 1<<xSize*ySize; 
    Bdim = Dim/Adim; 

    // Initialize the matrix of eigenvalues
    SuperMat.resize(Adim,Bdim); 
    SuperMat=0;
  }
  // Make the SuperMat of Eigs (as in 1d example)
  // Figure out dimensions of region A,B (xSize * ySize)
  // If it's > xMax * yMax/2 then switch regions A and B
  // In the future we can just multiply all renyis by 2 except the middle one for an even system.




  // ---- Vertical ----

  // ------ Corner Terms!! ------

}
#endif
