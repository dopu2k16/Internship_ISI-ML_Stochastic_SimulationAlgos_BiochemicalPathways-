#include<iostream>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<vector>
#include <random>
#include<chrono>
#include<algorithm>
#include<windows.h>
#include<array>
#include<iomanip>
#include<cstring>
#include<random>
#define RAND_UNIFORME (double)rand()/(double)RAND_MAX

#define n 4 //no. of species
#define m 3 //no. of reactions
#define tend 1
#define cr 0.1
#define epsilon 0.03
#define N 10
int sto[n][m], x[n],kj[m];

using namespace std;

void inti(int sto[n][m]){
sto[0][0]=-1;
sto[0][1]=1;
sto[0][2]=1;
sto[1][0]=-1;
sto[1][1]=1;
sto[1][2]=0;
sto[2][0]=1;

sto[2][1]=-1;
sto[2][2]=-1;
sto[3][0]=0;
sto[3][1]=0;
sto[3][2]=1;

}

//Calculate propensity values for each
void calpropen(int x[],double p[]){
p[0]=cr*x[0]*x[1];
p[1]=cr*x[2];
p[2]=cr*x[2];

}
double sumpropen(double a[]){
double sp=0;
int i;
for( i=0;i<m;i++)
    sp+=a[i];
    return sp;
}

  int select_rxn(double p[], int pn, double sum_propencity, double r){
	int i,reaction = -1;
	double sp = 0.0;
	//int i;
	r = r * sum_propencity;
	for( i=0; i<pn; i++){
		sp += p[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}
void update_count(int x[],int rxn){
int i;
for(i=0;i<n;i++){
    if(sto[i][rxn]!=0)
        x[i]+=sto[i][rxn];
}
}
void display(int x[],double time,int rec){
if(rec!=-1){
if(time<tend){
cout<<"\n";
cout<<"\tAT Time= "<<time;
cout<<"\n";
if(time!=0.0)
cout<<"\t Reaction No : "<<rec+1<<endl;

cout<<"\n";
cout<<"\t Species\n";
cout<<"\n";
cout<<"E\tS\tES\tP"<<endl;
cout<<"---------------------------------------------------------------------------"<<endl;
cout<< x[0]<<"\t"<< x[1]<<"\t"<< x[2]<<"\t"<< x[3]<<"\t";
cout<<"\n\n\n";
}
}
}


void computeMue_sigma(double mean[], double sigma[],double p[],int sto[n][m]) {

double tempfloat=0.0;
for(int nor=0;nor<m;++nor){
for(int nos=0;nos<n;++nos){
        tempfloat=sto[nos][nor]*p[nor];
mean[nor]+=tempfloat;
sigma[nor]+=sto[nos][nor]*tempfloat;

}
}
}


double computeTimeStep(int x[], double p[],int sto[n][m]){
double tauPrime;
double mean[m],sigma[m];
long double var[m];


memset(mean,0,sizeof(mean));
memset(sigma,0,sizeof(sigma));

computeMue_sigma(mean,sigma,p,sto);
/*for(int is=0;is<n;is++)

{
     cout<<"Mean :"<<mean[is]<<endl;
     //cout<<"Sigma"<<sigma[is]<<endl;
}
*/

double tau,tau1,tau2,epsi,epsixi;

tau=HUGE_VALF;
double sumpropensity=sumpropen(p);

for(int is=0;is<m;++is)
{
    var[is]=fabs(sigma[is]-((long double)(1.0/sumpropensity)*pow(mean[is],2)));
     cout<<"VAR :"<<var[is]<<endl;
    }

     epsi = epsilon*sumpropensity;
     epsixi = epsi *epsi;

    for(int j=0;j<m;j++){


				tau1 = (epsi/fabs(mean[j]));

				tau2=epsixi/fabs(var[j]);
				tau = min(tau,tau1);
                tau=min(tau,tau2);

               }


  tauPrime=tau;

  return tauPrime;
}

 unsigned long random_binomial(unsigned long no, double p)
{
double t=p/(1-p);
double u=RAND_UNIFORME;
double p0=pow((1-p),no);
double g=p0;
unsigned int k=0;
while (u>g)
    {
    p0*=t*(no-k)/(k+1);
    g+=p0;
    k++;
    };
return k;
}



void update_leap(int k,long int nooftimes,int x[]){
for(int i=0;i<n;i++)
x[i]+=nooftimes*sto[i][k];
}


double DirectMethod(int x[],double p[],double tim){
    double t2=tim,tau,sum_propensity=0.0;
    int reac=0;
int counter=0;
srand(time(NULL));
while(counter<10 && t2<tend){
        //while(c<100){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;
double r2=(double)rand()/RAND_MAX;
calpropen(x,p);
sum_propensity=sumpropen(p);
//display(x,t2,reac);
tau=((-log(r1))/sum_propensity);

reac=select_rxn(p,m,sum_propensity,r2);
update_count(x,reac);
t2+=tau;
display(x,t2,reac);
//t2+=tau;
counter++;
Sleep(2000);
}
return t2;
}


void print(int x[], double ti){

cout<<"Tau Leaping"<<endl;
cout<<"At Time : "<<ti<<endl;
cout<<"\t Species\n";
cout<<"\n";
cout<<"E\tS\tES\tP"<<endl;
cout<<"---------------------------------------------------------------------------"<<endl;
cout<< x[0]<<"\t"<< x[1]<<"\t"<< x[2]<<"\t"<< x[3]<<"\t";
cout<<"\n\n\n";
}

void update_lj(int l[],int x[],int rec){
  if(rec==1)
        l[0]=min(x[0],x[1]);
else if(rec==2)
   l[1]=x[2];
    else
    l[2]=x[2];
}


int main(){
   cout<<"\n";
cout<<"Reaction 1: E + S->ES "<<endl;
cout<<"Reaction 2: ES->E + S"<<endl;
cout<<"Reaction 3: ES->E + P "<<endl;
cout<<"\n\n";
cout<<"Enter the initial count for species E"<<endl;
cout<<"Enter the initial count for species S"<<endl;
cout<<"Enter the initial count for species ES"<<endl;
cout<<"Enter the initial count for species P"<<endl;;
int l[m];
double t=0.0;
 int kj;
 double p[m];
	double a0 = 0.0;
	double timepoint;
for(int i=0;i<n;i++)
cin>>x[i];

inti(sto);

/*for(int j=0;j<m;j++){
for(int i=0;i<n;i++)
    {
        if(sto[i][j]<0)
l[j]=min(l[i],x[i]/abs(sto[i][j]))
    }
} */


     srand(time(NULL));
l[0]=min(x[0],x[1]);
l[1]=x[2];
l[2]=x[2];
     while(t<tend)
     {
         timepoint=0.0;
         calpropen(x,p);
         a0=sumpropen(p);
         timepoint=computeTimeStep(x,p,sto);
          cout<<timepoint<<endl;
         if(timepoint>(N/a0))
         {

              for(int j=0;j<m;j++){
                if(p[j]*timepoint>l[j])
                   {

                    kj=random_binomial(l[j],1);
                     update_leap(j,kj,x);   //updating species count after executing each reaction
                      update_lj(l,x,j+1);  //updating the value of Lj after each Reaction Rj fires
                cout<<"Reaction No : "<<j+1<<"occurs "<<kj<<" times"<<endl;
                   }
                else {
                 kj=random_binomial(l[j],p[j]/l[j]);
                  update_leap(j,kj,x);
                cout<<"Reaction No : "<<j+1<<" occurs "<<kj<<" times"<<endl;
                  update_lj(l,x,j+1);
                }
              }
               t+=timepoint;
         print(x,t);
         }

     else{

        double t1=DirectMethod(x,p,t);
        t+=t1;

     }
    Sleep(2000);
     }
return 0;
 }
