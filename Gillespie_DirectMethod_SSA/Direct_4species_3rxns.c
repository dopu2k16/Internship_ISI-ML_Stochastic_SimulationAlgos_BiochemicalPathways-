#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#define n 4//no. of species
#define m 3//no. of reactions
#define rc 0.1//Rate constant
#define tend 1//End Time
int sto[n][m], x[n];//Species matrix
double  p[n];//propensity vector

void inti(int sto[n][m]){//Initialization of Stochiometric  Matrix
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
void calpropen(int x[]){
p[0]=rc*x[0]*x[1];
p[1]=rc*x[2];
p[2]=rc*x[2];
//p[3]=c[0]*x[o]*x[1];
}
double sumpropen(double a[]){ //calculation of a0
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
void update_count(int x[],int rxn,int sto[n][m]){
int i;
for(i=0;i<n;i++){
    if(sto[i][rxn]!=0)
        x[i]+=sto[i][rxn];
}
}
void display(int x[],double time,int rec){
printf("\n");
printf("\tAT Time=%f\n ",time);
printf("\n");
if(time!=0.0)
printf("\t Reaction No : %d\n",rec+1);

printf("\n");
printf("\t Species\n");
printf("\n");
printf(" Species E  Species S  Species ES Species P\n");
printf(" %d\t %5d\t %10d\t %d\t  ",x[0],x[1],x[2],x[3]);
printf("\n\n\n");

}

int main(){
    double t=0.0,tau,sum_propensity=0.0;
    int reac=0,i;
    //char ch;
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


inti(sto);
srand(time(NULL));
while(t<tend){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;
double r2=(double)rand()/RAND_MAX;
calpropen(x);
sum_propensity=sumpropen(p);
display(x,t,reac);
tau=((-log(r1))/sum_propensity);

reac=select_rxn(p,m,sum_propensity,r2);
update_count(x,reac,sto);
t+=tau;
sleep(3);
//ch=getch();
}
return 0;
}







