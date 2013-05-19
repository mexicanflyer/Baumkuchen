#include<sys/stat.h>
#include<iostream>
#include<sstream>
#include<fstream>
#include<math.h>
#include<complex>
#include<stdlib.h>
#define I_SQRTPI    0.56418958354776   /* 1/sqrt(pi) */

using namespace std;
// lapack
extern "C" void zgetrf_( int*, int* , complex<double>* , int*, int* , int* );
extern "C" void zgetri_( int*, complex<double>* , int*, int* , complex<double>*, int* , int* );
extern "C" void zgemm_( char*, char*, int*, int*, int*, complex<double>*, complex<double>*, int*, complex<double>*, int*, complex<double>*,complex<double>*, int* );

//matrix operations
void zgeTranspose( complex<double> *Transposed, complex<double> *M ,int n)
{
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) Transposed[i+n*j] = M[i*n+j];
}

void zgeDagger( complex<double> *Daggered, complex<double> *M ,int n)
{
  for(int i=0;i<n;i++)for(int j=0;j<n;j++) Daggered[i+n*j] = conj(M[i*n+j]);
}

void zgeProduct(complex<double> *AB,complex<double> *A, complex<double> *B, int N)
{
  char TRANSA='N',TRANSB='N';
  complex<double> alph=1.0,beta=0.0;
  zgemm_(&TRANSA,&TRANSB,&N,&N,&N,&alph,A,&N,B,&N,&beta,AB,&N);
}

void MatxVec(complex<double> *Tmp2,complex<double> *K,complex<double> *Tmp,int N)
{
  for(int j=0;j<N;j++){
    Tmp2[j]=0.0;
    for(int i=0;i<N;i++) Tmp2[j]+=K[i+j*N]*Tmp[i];} 
}

void zgeInverse(complex<double> *invA, complex<double> *A, int N)
{
  int LWORK=10*N;
  int *permutations;
  complex<double> *WORK, *tempA;
  
  tempA = new complex<double>[N*N];
  permutations =  new int[2*N];
  WORK = new complex<double>[N*N];
  int INFO;
  
  zgeTranspose(tempA,A,N);
  zgetrf_( &N, &N, tempA , &N, permutations , &INFO );
  if (INFO != 0) {
    cout<<"zgeInverse: Error at zgetrf  \n"; exit(0);
  }
  zgetri_( &N, tempA , &N, permutations , WORK, &LWORK, &INFO );
  if (INFO != 0) {
    cout<<"zgeInverse: Error at zgetri  \n"; exit(0);
  }
  zgeTranspose(invA,tempA,N);
  delete[] WORK;
  delete[] tempA;
  delete[] permutations;
}

//

void CreateCayleyOperator(complex<double> *Cayley,complex<double> a,complex<double> delx,complex<double> delt,int N)
{
  complex<double> One[N*N],S[N*N],InvS[N*N],D[N*N],ISD[N*N],Deno[N*N],Nume[N*N];
  complex<double> invdelx;
  int i,j;
  
  complex<double> I(0.0,1.0);
  
  for(i=0;i<N;i++)One[i+N*i]=1.0;
  
  for(i=0;i<N;i++){
    if(i>0)D[i+N*i-1]=1.0;
    D[i+N*i]=-2.0;
    if(i<N)D[i+N*i+1]=1.0;
  }
  for(i=0;i<N*N;i++)S[i]=One[i]+D[i]/10.0;
  zgeInverse(InvS,S,N);
  zgeProduct(ISD,InvS,D,N);
  invdelx=1.0/delx;
  
  for(i=0;i<N*N;i++){
    Deno[i]=One[i]-a*0.5*(0.0 + I)* delt*invdelx*invdelx * ISD[i];
    Nume[i]=One[i]+a*0.5*(0.0 + I)* delt*invdelx*invdelx * ISD[i];
  }
  
  zgeInverse(ISD,Deno,N);
  zgeProduct(Cayley,ISD,Nume,N);
}

void WriteIt(complex<double> *Mesh,double width,int N,string FileName)
{
  double x,y,p;
  int i,j;
  ofstream ofs;
  ofs.open(FileName.c_str());
  ofs.setf(ios::scientific);
  for(i=0;i<N;i++){
    x=(double)i*width/(double)N;
    for(j=0;j<N;j++){
      y=(double)j*width/(double)N;
      p=abs(Mesh[i+N*j]);
      ofs << x << " " << y << " " << p*p << endl;
    }
    ofs << endl;
  }
  ofs.close();
  ofs.clear();
}

string int2str(int i){
   stringstream ss;
   ss << i;
   return ss.str();
}

void SetGaussian(complex<double> *Mesh, double *r0,double *p0,double vX,double width,int N,complex<double> b)
{
  complex<double>I(0.,1.);
  double x,y,dx,dy,vP;
  complex<double> psi_x,psi_y;
  
  vP=1.0/vX;
  for(int j=0;j<N;j++){
    y=(double)j*width/(double)N;
    dy = y - r0[1];
    psi_y = exp(-.5*vP*vP*dy*dy+I*p0[1]*dy);
    for(int i=0;i<N;i++){
      x=(double)i*width/(double)N;
      dx = x - r0[0];
      psi_x = exp(-.5*vP*vP*dx*dx+I*p0[0]*dx);
      Mesh[i+N*j] = vP*I_SQRTPI*psi_y*psi_x *(cos(b*(x-width/2)*y) + I*sin(b*(x-width/2)*y));
    }
  }
}

void EvolveIt(complex<double> *Mesh,complex<double> *K,int N,int xy)
{
  complex<double> Tmp[N],Tmp2[N];
  if(xy==0){
    //x
    for(int j=0;j<N;j++){
      for(int i=0;i<N;i++)Tmp[i]=Mesh[i+N*j];
      MatxVec(Tmp2,K,Tmp,N);
      for(int i=0;i<N;i++)Mesh[i+N*j]=Tmp2[i];
    }
  }else{
    //y
    for(int i=0;i<N;i++){
      for(int j=0;j<N;j++)Tmp[j]=Mesh[i+N*j];
      MatxVec(Tmp2,K,Tmp,N);
      for(int j=0;j<N;j++)Mesh[i+N*j]=Tmp2[j];
    }
  }
}

void EvolveBxy(complex<double> *Mesh,complex<double> b,int N, double width, int PM)
{
  complex<double> I(0.,1.);
  double x,y;
  double db=real(b);
  for(int i=0;i<N;i++){for(int j=0;j<N;j++){
      x=(double)i*width/(double)N;
      y=(double)j*width/(double)N;
      Mesh[i+j*N]=(cos(db*x*y) +(double)PM*I*sin(db*x*y))*Mesh[i+j*N];
    }}
}

///////////////////////////////////////////////////////////////////////////////
int main()
{
  int i,j;
  const int N = 150;
  
  complex<double> I(0.,1.);
  complex<double> a,delx;
  complex<double> Mesh[N*N],Kb[N*N],Kx[N*N],Ky[N*N];
  
  double width=10.0,vX=0.4;
  complex<double> delt=0.01;
  
  // case 1 with momentum
  //complex<double> b=0.0;  
  //double r0[2]={5.0,5.0};
  //double p0[2]={10.0,0.0};

  // case 2 cyclotron
  complex<double> b=4.0;  
  double r0[2]={6.0,5.0};
  double p0[2]={0.0,0.0};
  
  SetGaussian(Mesh,r0,p0,vX,width,N,b);
  delx=(complex<double>)(width/(double)N);
  a=(complex<double>)(0.5*(1.0/6.0 - b*b*delt*delt/72.0));
  CreateCayleyOperator(Kb,a,delx,delt,N);
  a=(complex<double>)(0.25);
  CreateCayleyOperator(Kx,a,delx,delt,N);
  a=(complex<double>)(1.0/3.0);
  CreateCayleyOperator(Ky,a,delx,delt,N);
  
  
  mkdir("./data/",0777);
  complex<double>x,y;
  for(i=0;i<100;i++){
    string FileName="./data/"+int2str(i);
    WriteIt(Mesh,width,N,FileName);
    EvolveIt(Mesh,Kb,N,1);
    EvolveBxy(Mesh,b,N,width,-1);
    EvolveIt(Mesh,Kx,N,0);
    EvolveBxy(Mesh,b,N,width,+1);
    EvolveIt(Mesh,Ky,N,1);
    EvolveBxy(Mesh,b,N,width,-1);
    EvolveIt(Mesh,Kx,N,0);
    EvolveBxy(Mesh,b,N,width,+1);
    EvolveIt(Mesh,Kb,N,1);
    cout << i<<endl;
  }

  return 0; 
}
