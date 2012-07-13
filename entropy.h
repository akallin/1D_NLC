//function to calculate the RDM and entropies for a 1d 
#ifndef entropy_H
#define entropy_H

inline void Entropy1D(Array<l_double,1>& eigs)
{
  long int Dim = eigs.size();
  int Nsite = log2(Dim); //cout << "Nsite = " << Nsite << endl;
  // int Nsite = 16; //This should be an input parameter
  // long int Dim = (int)pow(2.0,Nsite);
  // if(Dim%2==1){cout << "DIMENSION ERRRRRRRROR" << endl;}
  //  vector <long double> eigs; //This will be the vector passed from the main program
  //can just calculate the number of sites from the vector length actually!
  //  eigs.resize(Dim,-1);

  // Create original basis
  //  for (unsigned long i1=0; i1<Dim; i1++){
  //      Basis.push_back(i1);
  //   BasPos.at(i1)=Basis.size()-1;
  //   Vdim++;
  //   }
  // }
  //  int Bsize = 9;
  //  int Dim9 = 512;
  //  int states=0;
  //  BasPosMat.resize(Dim9,-1);
  //Make a basis position vector for only region B
  //  for (unsigned long i1=0; i1<Dim9; i1++){
  //      BasPosMat.at(i1)=states; //then this is a possible state in region B
  //      states++;
  //    }
  //  }
  
  int Adim=1;
  int Bdim=Dim;
  


  Array<long double,2> SuperMat;
  int a(0),b(0);
  long double renyi(0);
  l_double norm(0);
  Array<double,2> DM(Adim,Adim);
  Array<double,2> DMsq(Adim,Adim);
  long double temp2(0);
  long double temp3(0);
  long double magnetization(0);
  

  vector<double> dd;
  long double vN;
  long double temp5;
  
  //for(int Asite=1; Asite<Nsite; Asite++){
  for(int Asite=1; Asite<Nsite/2+1; Asite++){
    Adim*=2;
    Bdim/=2;
    SuperMat.resize(Adim,Bdim); 
    SuperMat=0;
    a=0;b=0;
    temp3=0;
    norm=0;
    magnetization=0;
    // cout << "Adim = " << Adim << "  Bdim = " << Bdim << endl;
    //if(Asite==1){
    //  cout << eigs << endl;}
    
    for(int i=0; i<Dim; i++){ 
      // extractifying the region A and region B states
      //b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
      //c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
      a = i&(Adim-1);
      b = (i>>Asite)&(Bdim-1);
      
      SuperMat(a,b) = eigs(i);//sqrt(Dim)/Asite*eigs(i); 
      //    cout << "Adim = " << Adim << "  Bdim = " << Bdim << endl;
      //cout << "a " << a << "   b " << b << endl;
      temp3=0;
      if(Asite==1){
	for (int sp=0; sp<Nsite; sp++){
	  temp3 += (i>>sp)&1; 
	}
	magnetization += abs(temp3*2-Nsite)*eigs(i)*eigs(i);
	norm += eigs(i)*eigs(i);
	//cout << "i: " << i << "  mag: " << temp3*2-Nsite << endl;
	temp3=0;
      }
    }
    
    if(Asite==1){
      //  cout << "Mag = ";
      cout << magnetization/Nsite/norm << "   ";
      // cout << "Dim = " << Dim << "    Norm = " << norm << endl<< endl;
    }
    magnetization = 0;
    norm = 0;
    //   if(Adim==2||Bdim==2) {cout << SuperMat << endl;}
    
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
    
    //Diagonalizing the RDM (don't need this for Renyi, right??)
    while(dd.size()>0){dd.erase(dd.begin());}
    diagWithLapack_R(DM,dd); 
    renyi=0; vN=0; temp5=0;
    
    for(int s=0;s<dd.size();s++){
      renyi+=dd[s]*dd[s];
      temp5=log(dd[s]);
      if(!(temp5>-1000000000)){temp5=0;}
      vN+=-dd[s]*temp5;
    }
    
    cout <<" "<< setprecision(15) << -log(renyi) <<" "<< setprecision(15) << vN;    
    
    
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
    
    //  cout << "Norm"  << "     " << setprecision(15) << norm << endl;
    //cout << "Asites = " << Asite << "   Renyi"  << "  " << 
    cout <<" "<< setprecision(15) << -log(renyi/norm) << endl;
    
    //  cout << endl;
  }
  cout << endl;
}
#endif
