#include <algorithm>
#include <iostream>
#include <vector>
#include<conio.h>
#include <iterator>
#include<cstdlib>
#include<cmath>
#include<ctime>
#include <functional>
#include<climits>
#include<windows.h>
#include<time.h>
#include<stdlib.h>
#include<cstring>
#include<map>
#include<list>
#include<utility>
#include<new>
#define n 8 //no. of species
#define m 5   //no. of reactions
#define tend 1    //time to end
#define cr 0.1  //rate constant for all reaction channels
using namespace std;

int sto[n][m],x[n];  //x[] is species array, sto[] is stochiometric matrix , dep[][] is dependency graph matrix
double p[m];   //p[] is propensity array , t[] is time array
double arr[2*m-1];


void intia(int sto[n][m]){

      int i,j;
for(i=0;i<n;i++){  //Setting all elements of sto[][] & dep[][] to zero
    for(j=0;j<m;j++){
    sto[i][j]=0;
    //dep[i][j]=0;
    }
    }

sto[0][0]=-1;   //Initialization of stochiometric matrix sto[][] for a given set of 7 species & 5 reaction channels
sto[0][4]=1;

sto[1][0]=-1;
sto[1][4]=1;

sto[2][0]=1;
sto[2][1]=-1;
sto[2][2]=1;
sto[2][3]=-2;
sto[3][1]=1;

sto[5][5]=-1;


sto[6][3]=1;
sto[7][4]=-1;


 }

//Calculate propensity values for first time by calling this function from main()
void callpropen(int x[]){
p[0]=cr*x[0]*x[1];
p[1]=cr*x[1]*x[2];
p[2]=cr*x[2]*x[4];
p[3]=cr*x[2]*x[2];
p[4]=cr*x[5]*x[7];

}
//Recalculate the propensity functions for each dependent reactions
void recalpropen(int x[],int rec){
switch(rec)
{
    case 1: p[0]=cr*x[0]*x[1];
            break;
            case 2:  p[1]=cr*x[1]*x[2];
                        break;
                        case 3: p[2]=cr*x[2]*x[4];
                               break;
                               case 4: p[3]=cr*x[2]*x[2];
                                      break;
                                      case 5: p[4]=cr*x[5]*x[7];
                                              break;
}
/*
if(rec==0)   //for rxn R1
    p[0]=cr*x[0]*x[1];
else if(rec==1)   //for rxn R2
    p[1]=cr*x[1]*x[2];
else if(rec==2)  //for rxn R3
p[2]=cr*x[2]*x[4];
else if(rec==3)   //for rxn R4
p[3]=cr*x[2]*x[2];
else
p[4]=cr*x[5]*x[7];  //for rxn R5
*/
}
//function to update the species count
void speciesupdate_count(int x[],int rxn){
int i;
switch(rxn)
{
case 1:
    x[0]--,x[1]++,x[2]--;
    break;
case 2:
    x[2]--,x[1]--,x[3]++;
    break;

case 3: x[2]++,x[4]--,x[5]++;
        break;

case 4: x[6]++,x[2]-=2;
       break;

case 5: x[5]--,x[7]--,x[1]++;
      break;
}
/*
for(i=0;i<n;i++){
    if(sto[i][rxn]!=0)
        x[i]+=sto[i][rxn];
}
*/
}

//function to print species count
void display_count(int x[],double time,int reac){
if(time<tend){
cout<<"\n";
cout<<"\tAT Time= "<<time;
cout<<"\n";

cout<<"\t Reaction No : "<<reac<<endl;

cout<<"\n";
cout<<"\t Species\n";
cout<<"\n";
cout<<"A\tB\tC\tD\tE\tF\tG"<<endl;
cout<<"---------------------------------------------------------------------------"<<endl;
cout<< x[0]<<"\t"<< x[1]<<"\t"<< x[2]<<"\t"<< x[3]<<"\t"<< x[4]<<"\t"<< x[5]<<"\t"<< x[6];
cout<<"\n\n\n";

}}


void BuildTree()
{
    int j=0;

    for(int i=m-1;i<2*m-1;i++)
    {
        arr[i]=p[j];     //filling the tree array with reaction propensities
        j++;
    }

    for(int i=m-2;i>=0;i--)
        arr[i]=arr[2*i+1]+arr[2*i+2];   //creating the array complete tree with sum of left and right
for(int i=0;i<2*m-1;i++)
    cout<<arr[i]<<endl;
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
double sum_propensity,sump,tau;
    map<int,list<int>> mymap{{1,{2,3,4}},{2,{1,3,4}},{3,{2,4,5}},{4,{2,3}},{5,{1,2}}};
//Dependency graph plotted as map Rxn No,{dependent reactions}
    map<int,list<int>>:: iterator it;
    list<int>:: iterator li;
//vector<double> tree(2*m-1);
    double t=0.0;
     //array to store old propensities values of each rxn channels
  cout<<"\n";
cout<<"Reaction 1: A + B->C \n";
cout<<"Reaction 2: C+ B->D \n";
cout<<"Reaction 3:  E+ C->2C+ F\n";
cout<<"Reaction 4: 2C-> G\n";
cout<<"Reaction 5: F + H->B\n";
cout<<endl;
cout<<"Enter the initial count for species A\n";
cout<<"Enter the initial count for species B\n";
cout<<"Enter the initial count for species C\n";
cout<<"Enter the initial count for species D\n";
cout<<"Enter the initial count for species E\n";
cout<<"Enter the initial count for species F\n";
cout<<"Enter the initial count for species G\n";
cout<<"Enter the initial count for species H\n";
for(int i=0;i<n;i++)
cin>>x[i];  //user input for no. of initial species

intia(sto);  //calling intia() function to initialise the stochiometric  matrix
 //callpropen(x);   //calculating propensities of each rxn channels for first time
 srand(time(NULL));

callpropen(x);
//build_array(p);

while(t<tend){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;
       // printf("\n");
        //printf("%f\n",r1);
double r2=(double)rand()/RAND_MAX;


BuildTree();
//sum_propensity=build_tree(p,tree);
sum_propensity=arr[0];
cout<<"Sum propensity"<<sum_propensity<<endl;
tau=((-log(r1))/sum_propensity);
sump=r2*sum_propensity;

int selrxn=searchY(0,sump);

cout<<"Sel rxn:"<<selrxn+1;
it=mymap.find(selrxn);
    //cout<<(*it).first<<"{";
speciesupdate_count(x,selrxn);

        for(li=(*it).second.begin();li!=(*it).second.end();li++)
        {
          // cout<<(*li);
           recalpropen(x,(*li));
           //update_proparray(arr,p,(*li));
        }

t+=tau;
display_count(x,t,selrxn);
//ch=getch();
Sleep(3000);
}
return 0;

}
