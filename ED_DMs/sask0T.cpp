#include"header.h"
#include"lapack.h"
#include<vector>
#include<fstream>
#include <blitz/array.h>
BZ_USING_NAMESPACE(blitz)
using namespace std;

int main()
{
  long int Dim = 65536, Vdim = 12870;
  int Nsite = 16;
  int temp;
  vector <long int> Basis, BasPosMat;
  vector <long double> eigs;

  // Create original basis
  for (unsigned long i1=0; i1<Dim; i1++){
    temp = 0;
    for (int sp =0; sp<Nsite; sp++){
      temp += (i1>>sp)&1; 
    }
    if (temp==(Nsite/2) ){ 
      Basis.push_back(i1);
      //   BasPos.at(i1)=Basis.size()-1;
      //   Vdim++;
    }
  }

  int Bsize = 9;
  int Dim9 = 512;
  int states=0;
  BasPosMat.resize(Dim9,-1);
  eigs.resize(Vdim,-1);
    
  //Make a basis position vector for only region B
  for (unsigned long i1=0; i1<Dim9; i1++){
    temp = 0;
    for (int sp =0; sp<Bsize; sp++){
      temp += (i1>>sp)&1; 
    }
    if (temp<=8){ //if half or less of the spins are up
      BasPosMat.at(i1)=states; //then this is a possible state in region B
      states++;
    }
  }
  
  int Adim = 128;
  int Bdim = 512;
  long double Z = 0;  // The normalization?
  Array<long double,2> supermat(Adim,states); 
  int a(0),b(0),c(0),d(0);
  long double temp3=0;
  double temp2=0;
  Array<double,2> DM(Adim,Adim);
  vector<double> dd;
  long double renyi(0),vN(0),temp5(0);

  ifstream fin("XYzerovector.dat");
  supermat=0;

  a=0;b=0;c=0;d=0;
  temp3=0;
  for(int i=0; i<Vdim; i++){ 
	  a = Basis[i];
	  //somehow extractying the region A and region B states
	  b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
	  c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
	  d = BasPosMat[c];
	  fin >> temp3;
	  //this matrix is like the sqrt of the RDM (sorta)
	  supermat(b,d) = temp3; 
  }



  //multiplying the supermat by its transpose to get the RDM
  temp2=0;
  DM=0;
  for(int i=0; i<Adim; i++){
    for(int j=0; j<Adim; j++){
      temp2=0;
      for(int k=0; k<states; k++){
	temp2 += supermat(i,k)*supermat(j,k);
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
	  //	cout << dd[s] <<endl;
  }

  cout << "XY T=0 Renyi"  << "     " << setprecision(15) << -log(renyi)*2 << endl;
  cout << "XY T=0 von Neumann"  << "     " << setprecision(15) << vN*2 << endl;

  
  // Doing it all again for the Heisenberg model and then for XXZ

  ifstream fin1("HEISzerovector.dat");
  supermat=0;

  a=0;b=0;c=0;d=0;
  temp3=0;
  for(int i=0; i<Vdim; i++){ 
	  a = Basis[i];
	  b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
	  c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
	  d = BasPosMat[c];
	  fin1 >> temp3;
	  supermat(b,d) = temp3; 
  }


  temp2=0;
  DM=0;
  for(int i=0; i<Adim; i++){
    for(int j=0; j<Adim; j++){
      temp2=0;
      for(int k=0; k<states; k++){
	temp2 += supermat(i,k)*supermat(j,k);
      }
      DM(i,j) = temp2;
    }
  }
  

  while(dd.size()>0){dd.erase(dd.begin());}
  diagWithLapack_R(DM,dd); 
  renyi=0; vN=0; temp5=0;

  for(int s=0;s<dd.size();s++){
	  renyi+=dd[s]*dd[s];
	  temp5=log(dd[s]);
	  if(!(temp5>-1000000000)){temp5=0;}
	  vN+=-dd[s]*temp5;
	  //	cout << dd[s] <<endl;
  }

  cout << "Heis T=0 Renyi"  << "     " << setprecision(15) << -log(renyi)*2 << endl;
  cout << "Heis T=0 von Neumann"  << "     " << setprecision(15) << vN*2 << endl;


  ifstream fin2("XXZzerovector.dat");
  supermat=0;

  a=0;b=0;c=0;d=0;
  temp3=0;
  for(int i=0; i<Vdim; i++){ 
	  a = Basis[i];
	  b = (((a>>11)&31)<<2)+(((a>>7)&1)<<1)+((a>>3)&1);
	  c = (((a>>8)&7)<<6)+(((a>>4)&7)<<3)+(a&7);
	  d = BasPosMat[c];
	  fin2 >> temp3;
	  supermat(b,d) = temp3; 
  }


  temp2=0;
  DM=0;
  for(int i=0; i<Adim; i++){
    for(int j=0; j<Adim; j++){
      temp2=0;
      for(int k=0; k<states; k++){
	temp2 += supermat(i,k)*supermat(j,k);
      }
      DM(i,j) = temp2;
    }
  }
  

  while(dd.size()>0){dd.erase(dd.begin());}
  diagWithLapack_R(DM,dd); 
  renyi=0; vN=0; temp5=0;

  for(int s=0;s<dd.size();s++){
	  renyi+=dd[s]*dd[s];
	  temp5=log(dd[s]);
	  if(!(temp5>-1000000000)){temp5=0;}
	  vN+=-dd[s]*temp5;
	  //	cout << dd[s] <<endl;
  }

  cout << "XXZ T=0 Renyi"  << "     " << setprecision(15) << -log(renyi)*2 << endl;
  cout << "XXZ T=0 von Neumann"  << "     " << setprecision(15) << vN*2 << endl;



  return 0;
}
