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
#define n 7 //no. of species
#define m 5   //no. of reactions
#define tend 1    //time to end
#define cr 0.1  //rate constant for all reaction channels
using namespace std;

int sto[n][m],x[n],dep[m][m];  //x[] is species array, sto[] is stochiometric matrix , dep[][] is dependency graph matrix
double t[m],p[m];   //p[] is propensity array , t[] is time array

void intia(int sto[n][m], int dep[m][m]){

      int i,j;
for(i=0;i<n;i++){  //Setting all elements of sto[][] & dep[][] to zero
    for(j=0;j<m;j++){
    sto[i][j]=0;
    dep[i][j]=0;
    }
    }

sto[0][0]=-1;   //Initialization of stochiometric matrix sto[][] for a given set of 7 species & 5 reaction channels
sto[0][4]=1;

sto[1][0]=-1;
sto[1][1]=-1;

sto[2][0]=1;
sto[2][1]=-1;

sto[3][1]=1;
sto[3][2]=-1;
sto[3][3]=1;

sto[4][4]=-1;

sto[5][2]=1;
sto[5][3]=-1;

sto[6][3]=1;
sto[6][4]=-1;

dep[0][0]=-1;   //Initialization of dependency graph matrix dep[][]
dep[1][0]=1;    //for self loop dep[][]= -1
dep[1][1]=-1;          //for other dependent reactions i.e other edges dep[][]= 1
dep[1][2]=1;
dep[2][2]=-1;
dep[2][3]=1;
dep[3][2]=1;
dep[3][3]=-1;
dep[3][4]=1;
dep[4][0]=1;
dep[4][2]=1;
dep[4][4]=-1;
}

//Calculate propensity values for first time by calling this function from main()
void callpropen(int x[]){
p[0]=cr*x[0]*x[1];
p[1]=cr*x[1]*x[2];
p[2]=cr*x[3]*x[4];
p[3]=cr*x[5];
p[4]=cr*x[4]*x[6];


}
//Recalculate the propensity functions for each dependent reactions
void recalpropen(int x[],int rec){
    if(rec==0)   //for rxn R1
    p[0]=cr*x[0]*x[1];
else if(rec==1)   //for rxn R2
    p[1]=cr*x[1]*x[2];
else if(rec==2)  //for rxn R3
p[2]=cr*x[3]*x[4];
else if(rec==3)   //for rxn R4
p[3]=cr*x[5];
else
p[4]=cr*x[4]*x[6];  //for rxn R5


}
//function to update the species count
void speciesupdate_count(int x[],int rxn){
int i;
for(i=0;i<n;i++){
    if(sto[i][rxn]!=0)
        x[i]+=sto[i][rxn];
}
}
//unnecessary function to print heap structure for debugging option only. Not called explicitly
void print(const vector < double>& v)
{
    vector < double>::const_iterator k;
    for(k = v.begin(); k!= v.end(); k++)
    {
        cout << *k << " ";
    }
    cout << endl;
}
//function to print species count
void display_count(int x[],double time,int reac){
if(time<tend){
cout<<"\n";
cout<<"\tAT Time= "<<time;
cout<<"\n";

cout<<"\t Reaction No : "<<reac+1<<endl;

cout<<"\n";
cout<<"\t Species\n";
cout<<"\n";
cout<<"A\tB\tC\tD\tE\tF\tG"<<endl;
cout<<"---------------------------------------------------------------------------"<<endl;
cout<< x[0]<<"\t"<< x[1]<<"\t"<< x[2]<<"\t"<< x[3]<<"\t"<< x[4]<<"\t"<< x[5]<<"\t"<< x[6];
cout<<"\n\n\n";

}}
int main(){
    int selrxn,i,j,key;

    double currenttime=0.0,newt,ti;
    double oldp[m]; //array to store old propensities values of each rxn channels
  cout<<"\n";
cout<<"Reaction 1: A + B->C \n";
cout<<"Reaction 2: B + C->D \n";
cout<<"Reaction 3: D + E->E + F\n";
cout<<"Reaction 4: F->D + G\n";
cout<<"Reaction 5: E+ G->A\n";
cout<<endl;
cout<<"Enter the initial count for species A\n";
cout<<"Enter the initial count for species B\n";
cout<<"Enter the initial count for species C\n";
cout<<"Enter the initial count for species D\n";
cout<<"Enter the initial count for species E\n";
cout<<"Enter the initial count for species F\n";
cout<<"Enter the initial count for species G\n";

for(i=0;i<n;i++)
cin>>x[i];  //user input for no. of initial species

intia(sto,dep);  //calling intia() function to initialise the stochiometric & dependency graph matrix
 callpropen(x);   //calculating propensities of each rxn channels for first time
 srand(time(NULL));
for(i=0;i<m;i++){

    t[i]=(-log(( double)rand()/RAND_MAX))*(1/p[i]);   //calculating initial putative time for each rxn channel
  // cout<<-log(( double)rand()/RAND_MAX)<<endl;
   //cout<<t[i]<<endl;
   }
    vector < double> v(t, t + sizeof(t)/sizeof( double));  //creating a vector v of time
    vector <int>::iterator k;
    make_heap(v.begin(), v.end(), greater < double>());//creating a min heap structure for time
   // print(v);
    for(i=0;i<m;i++)
 oldp[i]=p[i];    //storing the initial propensities to oldp[]

while(currenttime<tend){

    for(i=0;i<m;i++)
    if(v.front()==t[i])
        selrxn=i;    //key of top element(min time) of min heap

        ti=v.front();

    speciesupdate_count(x,selrxn);
     currenttime=ti;

   display_count(x,currenttime,selrxn);

    for(i=0;i<m;i++)  //for each edge(mu,alpha) in dependency graph matrix
      {


        if(dep[selrxn][i]==1) //checking of other dependent reaction on a strike of a reaction
        {
            oldp[i]=p[i];
           // speciesupdate_count(x,selrxn);
           recalpropen(x,i);


         for(j=0;j<m;j++)
            if(t[i]==v[j])
              key=j;


            newt=oldp[i]*(t[i]-currenttime);
            newt=newt/p[i] + currenttime;
            t[i]=newt;

        v[key]=t[i];

        make_heap(v.begin(), v.end(), greater < double>());  //updating the min heap with reusing time calculation
        //print(v);
        }
   // }
   else if(dep[selrxn][i]==-1)
   {
     //t[i]=((-log((double)rand()/RAND_MAX))/p[selrxn])+currenttime;
     oldp[selrxn]=p[i];
     recalpropen(x,i);

      for(j=0;j<m;j++)
            if(t[i]==v[j])
              key=j;


     t[selrxn]=((-log((double)rand()/RAND_MAX))/p[selrxn])+currenttime;

    v[key]=t[selrxn];

    make_heap(v.begin(), v.end(), greater < double>());  //updating heap
    //print(v);
   }

    }

Sleep(3000);
//char ch=cin.get();
}
return 0;
}





