#include<iostream>
#include<cmath>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<vector>
#include<chrono>
#include<random>
#include<math.h>
#include<algorithm>
#include<windows.h>
#define n 4 //no. of species
#define m 3 //no. of reactions
#define tend 1
#define cr 0.1
#define epsilon 0.03
#define nc 10
#define N 10

int sto[n][m], x[n],kj[m],j;
double p[n],mean[n],var[n];
int gi[4]={2,2,1,0};

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

double sumprop_critical(double p[],vector<int>& critical){
double sp=0.0;
for(int i=0;i<critical.size();i++){
    sp+=p[critical[i]];
}
return sp;
}

void critical_noncritical(int x[],double p[],int sto[n][m],vector<int>& critical, vector<int>& noncritical)
{
   double b[n],l[m];
    for(int i=0;i<m;i++){
for(int j=0;j<n;j++){
    b[j]=x[j]/abs(sto[j][i]);
}
l[i]=*min_element(b,b+n);
}

for(int j=0;j<m;j++){
    if(l[j]<nc && p[j]>0)
        //flag[j]=1;
        critical.push_back(j);  //if critical
    else
        //flag[j]=0;
        noncritical.push_back(j); //if non critical

}

}



 double sum_critical(double a){
double sp=0;
return sp+a;
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
if(time<tend){
cout<<"\n";
cout<<"\tAT Time= "<<time;
cout<<"\n";
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

void computeMue_sigma(double mean[], double sigma[],double p[],int sto[n][m],vector<int>&noncritical) {

double tempfloat=0.0;
for(int ir=0;ir<m;++ir){
        //if(flag[ir]==)
for(int is=0;is<n;++is){
        tempfloat+=sto[is][ir]*p[noncritical[ir]];
mean[is]+=tempfloat;
sigma[is]+=sto[is][ir]*tempfloat;

}
}
}


void print(int rec,int kj,int x[]){

    cout<<"Reaction :"<<rec+1<<" \t"<<"occurs "<< kj<<endl;
    cout<<"\t Species\n";
cout<<"\n";
cout<<"E\tS\tES\tP"<<endl;
cout<<"---------------------------------------------------------------------------"<<endl;
cout<< x[0]<<"\t"<< x[1]<<"\t"<< x[2]<<"\t"<< x[3]<<"\t";
cout<<"\n\n\n";
}


double computeTimeStep(int x[], double p[],int sto[n][m],vector<int> & noncritical){
double tauPrime;
double mean[n],sigma[n],var[n];

computeMue_sigma(mean,sigma,p,sto,noncritical);

double tau,epsi,epsixi,epsixisq;
double xi;
tau=HUGE_VAL;
double sumpropensity=sumpropen(p);

for(int is=0;is<n;++is)
{
    var[is]=fabs(sigma[is]-(1.0/sumpropensity)*mean[is]*mean[is]);

    }

for(int is=0;is<n;++is){

    //taup=(HUGE_VALF*0.5);
    xi=x[is];
     switch(gi[is]){
    case 0:
          break;
    case 1:
         epsi = epsilon;
				epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabsf(mean[is]));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/var[is]);
				break;
    case 2:
           epsi = 0.5*epsilon;
           epsixi = epsi * xi;
				epsixi = max(epsixi,1.0);
				tau = min(tau,epsixi/fabs(mean[is]));
				epsixisq = epsixi*epsixi;
				tau = min(tau,epsixisq/var[is]);
				break;
    default:
				break;



          }

  }
  tauPrime=tau;
  //cout<<tau<<endl;
  return tauPrime;
}



int getPoisson(double lambda) {
  //srand(time(NULL));
  double L = exp(-lambda);
  double p = 1.0;
  int k = 0;

  do {
    k++;
    p *= (double)rand()/RAND_MAX;
  } while (p > L);

  return k - 1;
}


/*
void update_leap(int ,long int nooftimes,int x[]){
for(int i=0;i<n;i++)
x[i]+=nooftimes*sto[i][k];
}
*/

double DirectMethod(int x[],double p[],double tim){
    double t2=tim,tau,sum_propensity=0.0;
    int reac=0;
int counter=0;
srand(time(NULL));
while(counter<10 && t2<tend){

        tau=0.0;
        double r1=(double)rand()/RAND_MAX;
double r2=(double)rand()/RAND_MAX;
calpropen(x,p);
sum_propensity=sumpropen(p);
display(x,t2,reac);
tau=((-log(r1))/sum_propensity);

reac=select_rxn(p,m,sum_propensity,r2);
update_count(x,reac);
t2+=tau;
counter++;
Sleep(2000);
}
return t2;
}



 int main(){

     vector<int>noncritical;
     vector<int>critical;
double t=0.0,t1,t2,ptau,ptau1,ptaudash,timepoint;
int b[n],l[m],flag[m];

   cout<<"\n";
cout<<"Reaction 1: E + S->ES "<<endl;
cout<<"Reaction 2: ES->E + S"<<endl;
cout<<"Reaction 3: ES->E + P "<<endl;
cout<<"\n\n";
cout<<"Enter the initial count for species E"<<endl;
cout<<"Enter the initial count for species S"<<endl;
cout<<"Enter the initial count for species ES"<<endl;
cout<<"Enter the initial count for species P"<<endl;;
for( int i=0;i<n;i++)
cin>>x[i];

inti(sto);
calpropen(x,p);

critical_noncritical(x,p,sto,critical,noncritical);


while(t< tend){


    calpropen(x,p);
  double a0=sumpropen(p);



  //  critical.resize(0);
    //noncritical.resize(0);
    //memset(critical,0,sizeof(critical));
    //memset(noncritical,0,sizeof(noncritical));

    critical_noncritical(x,p,sto,critical,noncritical);

    double ac=sumprop_critical(p,critical);


       if(noncritical.empty())
        timepoint=HUGE_VALF;
    else
    timepoint=computeTimeStep(x,p,sto,noncritical);

    label: if(timepoint<(N/a0))
    {

    //call direct method & execute it for 100 times by selrxn & tau calculation

    t+=DirectMethod(x,p,t);

    }
    else{
    double r3=(double)rand()/RAND_MAX;
    ptaudash=(1/ac*(1/r3));


    if(timepoint<ptaudash){
        ptau=timepoint;

      for(int j=0;j<critical.size();j++)
           kj[critical[j]]=0;

            for(int j=0;j<noncritical.size();j++)
                   kj[noncritical[j]]=getPoisson(p[j]*ptau);


      /*for(j=0;j<m;j++){
            if(flag[j]==1 &&p[j]>0)
            kj[j]=0;//for critical reactions
    else
    kj[j]= getPoisson(p[j]*ptau);
} */

   }

    else //if(ptaudash<=timepoint)
    {
        ptau=ptaudash;

       int rex=select_rxn(p,critical.size(),ac,(double)rand()/RAND_MAX);
            kj[rex]=1;  //for critical reactions

     for(int j=0;j<critical.size()&& j!=rex;j++)
        kj[critical[j]]=0;

    for(int j=0;j<noncritical.size();j++)
        kj[noncritical[j]]=getPoisson(p[noncritical[j]]*ptau);


      /*  rex=select_rxn(p,sizeof(critical),ac,(double)rand()/RAND_MAX);
            kj[rex]=1;  //for critical reactions
        for(int j=0;j<m && j!=rex;j++){
        if(flag[j]==1 && p[j]>0){
           //rex=select_rxn(p,nocr,ac, ran);
            kj[j]=0;  //for critical reactions
            }
else
    kj[j]= getPoisson(p[j]*ptau);
} */


}
for(int i=0;i<n;i++){
        double alpha=0.0;
    for(int j=0;j<m;j++){
        alpha+=x[i]+kj[j]*sto[i][j];

        }
        if(alpha<0.0)
        {
            timepoint=0.5*timepoint;
            goto label;
        }
        else
            x[i]+=alpha;
    }

t+=ptau;
}


}

return 0;

 }












