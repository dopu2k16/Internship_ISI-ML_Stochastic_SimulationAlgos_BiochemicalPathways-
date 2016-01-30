#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#define n 4//no. of species
#define m 3//no. of reactions
#define tend 1
#define cr 0.1
int sto[n][m], x[n]; //x[] is species array, sto[] is stochiometric matrix
double p[n],t[m];  //p[] is propensity array , t[] is time array

void intia(int sto[n][m]){ //Initialization of Stochiometric matrix
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
void callpropen(int x[]){
p[0]=cr*x[0]*x[1];
p[1]=cr*x[2];
p[2]=cr*x[2];
//p[3]=c[0]*x[o]*x[1];
}

void caltime(double p[]){ //function to calculate time
int i;
srand(time(NULL));
for(i=0;i<m;i++){

    t[i]=((-log((double)rand()/RAND_MAX))/p[i]);

}
 /*for(i=0;i<m;i++)
    printf("%f\n",t[i]); */ //function to print time(random numbers) generated at each iteration
    }

    void speciesupdate_count(int x[],int rxn,int sto[n][m]){
int i;
for(i=0;i<n;i++){
    if(sto[i][rxn]!=0)
        x[i]+=sto[i][rxn];
}
}

void display_count(int x[],double time,int reac){
printf("\n");
printf("\tAT Time=%f\n ",time);
printf("\n");
//if(time!=0.0)
printf("\t Reaction No : %d\n",reac+1);

printf("\n");
printf("\t Species\n");
printf("\n");
printf(" Species E  Species S  Species ES Species P\n");
printf("-----------------------------------------------------\n");
printf(" %d\t %5d\t %10d\t %d\t  ",x[0],x[1],x[2],x[3]);
printf("\n\n\n");

}
int main(){
    int selrxn,i;
    double currenttime=0.0, mintime;
    //double t[m];
  printf("\n");
printf("Reaction 1: E + S->ES \n");
printf("Reaction 2: ES->E + S\n");
printf("Reaction 3: ES->E + P\n");
printf("\n\n");
printf("Enter the initial count for species E\n");
printf("Enter the initial count for species S\n");
printf("Enter the initial count for species ES\n");
printf("Enter the initial count for species P\n");

for( i=0;i<n;i++)
scanf("%d",&x[i]);

intia(sto);
 //callpropen(x);
//double currenttime=0.0, mintime;

while(currenttime<tend){
        callpropen(x); //calculating propensities values
        caltime(p);
 //srand(time(NULL));
    selrxn=0;
    mintime=t[0];  //let mintime be the first time & then comparing it with others
    for(i=1;i<m;i++){
    if(mintime>t[i])   //searching for the minimum time & reaction mu to occur next
        {
            mintime=t[i];
            selrxn=i;
         }
    }

currenttime=currenttime+mintime;
speciesupdate_count(x,selrxn,sto);
display_count(x,currenttime,selrxn);
//caltime(p);
//callpropen(x);
sleep(3);
}
return 0;
}



