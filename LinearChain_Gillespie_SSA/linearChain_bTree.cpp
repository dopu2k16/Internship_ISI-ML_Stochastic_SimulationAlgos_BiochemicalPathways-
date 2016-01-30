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
double prop[m],sum_propensity=0.0;
double arr[2*m-1];

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


void BuildTree()
{
    int j=0;

    for(int i=m-1;i<2*m-1;i++)
    {
        arr[i]=prop[j];     //filling the tree array with reaction propensities
        j++;
    }

    for(int i=m-2;i>=0;i--)
        arr[i]=arr[2*i+1]+arr[2*i+2];   //creating the array complete tree with sum of left and right
//for(int i=0;i<2*m-1;i++)
  //  cout<<arr[i]<<endl;
}

int searchY(int pos,double val)   //search function for next reaction firing
{
    int p=pos;

    int left=2*p+1;
   // cout<<left<<endl;
    int right=2*p+2;
    //cout<<right<<endl;

         if(left>=m-1 && right>=m)
             return pos;

        if(val<=arr[left])
        {

         return  searchY(left,val);

        }
        else{
            val=arr[p]-val;
          return   searchY(right,val);
        //cout<<right<<endl;

        }


    }


int main(){
    //int selrxn;
    outFile.open("chain_Btree.csv");
double sump,tau;

    double t=0.0;
     //array to store old propensities values of each rxn channels

for(int i=0;i<n;i++)
    x[i]=100;

for(int i=0;i<m;i++)
    prop[i]=cr*x[i];
int counter=0;

srand(time(NULL));

while(counter<1000000){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;

double r2=(double)rand()/RAND_MAX;


BuildTree();

sum_propensity=arr[0];
//cout<<"Sum propensity"<<sum_propensity<<endl;
tau=((-log(r1))/sum_propensity);
sump=r2*sum_propensity;

int selrxn=searchY(0,sump);

//cout<<"Sel rxn:"<<selrxn+1;


update_species_prop(selrxn,x,sum_propensity);


t+=tau;
counter++;
//outFile<<t<<","<<selrxn+1<<","<<r2<<endl;

}

for(int i=0;i<n;i++)
outFile<<x[i]<<",";

outFile<<","<<endl;

outFile.close();

return 0;
}


