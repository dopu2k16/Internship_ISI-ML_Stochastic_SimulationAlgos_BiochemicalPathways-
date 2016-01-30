#include<iostream>
#include<fstream>
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#include<iomanip>
#define cr 0.1//np. of reactions

#define m 999
#define n 1000
using namespace std;
ofstream outFile;



//int m,n;
int x[n];
double prop[m],sum_prop=0.0;


void update_species_prop(int rxn,int *x,double sum_prop)
{
    x[rxn]-=1;
    x[rxn+1]+=1;
    prop[rxn]=cr*x[rxn];
    prop[rxn+1]=cr*x[rxn+1];
    //sum_prop+=prop[rxn+1]-prop[rxn];
}
void update_sumpropen(int rxn,double *prop,double sum_prop)
{
    sum_prop+=prop[rxn+1]-prop[rxn];

}

  int select_rxn(double *prop, int pn, double sum_propencity, double r){
	int i,reaction = -1;
	double sp = 0.0;

	r = r * sum_propencity;
	for( i=0; i<pn; i++){
		sp += prop[i];
		if(r < sp){
			reaction = i;
			break;
		}
	}
	return reaction;
}

/*
void display(int x[],double time,int rec){
//ofstream outFile;
outFile.open("Chain.csv");
outFile<<"Reaction No : "<<rec+1<<endl;
outFile<<"Time "<<endl;

for(int i=0;i<n;i++)
    outFile<<x[i]<<",";
outFile<<endl;
outFile.close();
}
*/
int main(){
//ofstream outFile;
outFile.open("Chain.csv");
     double  t=0.0,tau;//sum_propensity=0.0;
    int reac=0,counter =0;

   //ofstream outFile;


for(int i=0;i<n;i++)
    x[i]=100;

for(int i=0;i<m;i++)
    prop[i]=cr*x[i];

sum_prop=prop[0]*m;
//cout<<sum_prop;

//for(int i=0;i<n;i++)
  //  cout<<x[i];


srand(time(NULL));
while(counter<10000000){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;

double r2=(double)rand()/RAND_MAX;

//display(x,t,reac);
tau=((-log(r1))/sum_prop);

reac=select_rxn(prop,m,sum_prop,r2);
//cout<<reac+1<<endl;
update_species_prop(reac,x,sum_prop);

t+=tau;
//outFile<<t<<","<<reac+1<<","<<r2<<endl;

counter++;

}
for(int i=0;i<n;i++)
outFile<<x[i]<<",";

outFile<<","<<endl;

outFile.close();
return 0;
}








