#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#define n 7//no. of species
#define m 4//no. of reactions
#define rc 0.1 //Rate constant
#define tend 1 //Time to end is 1 sec
int sto[n][m], x[n];
double p[n];//propensity vector

void inti(int sto[n][m]){
sto[0][0]=-1;//Stochiometrix Matrix
sto[0][1]=0;
sto[0][2]=0;
sto[0][3]=0;

sto[1][0]=-1;
sto[1][1]=0;
sto[1][2]=-1;
sto[1][3]=0;

sto[2][0]=1;
sto[2][1]=-1;
sto[2][2]=0;
sto[2][3]=-1;

sto[3][0]=0;
sto[3][1]=-1;
sto[3][2]=0;
sto[3][3]=0;

sto[4][0]=0;
sto[4][1]=1;
sto[4][2]=-1;
sto[4][3]=0;

sto[5][0]=0;
sto[5][1]=0;
sto[5][2]=1;
sto[5][3]=-1;

sto[6][0]=0;
sto[6][1]=0;
sto[6][2]=0;
sto[6][3]=1;

}

//Calculate propensity values for each
void calpropen(int x[]){
p[0]=rc*x[0]*x[1];

p[1]=rc*x[2]*x[3];
p[2]=rc*x[1]*x[4];
p[3]=rc*x[2]*x[5];

}
double sumpropen(double a[]){ //Calculation of a0
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
    printf("\n");
printf("\tAT Time=%f\n ",time);
printf("\n");
if(time!=0.0)
printf("\t Reaction No : %d\n",rec);

printf("\n");
printf("\t Species\n");
printf("\n");
printf(" Species A  Species B  Species C Species D  Species E  Species F  Species G\n");
printf(" %d\t %5d\t %10d\t %d\t %5d\t %10d\t %5d\n ",x[0],x[1],x[2],x[3],x[4],x[5],x[6]);
printf("\n\n\n");
}

int main(){
    //char ch;
    double t=0.0,tau,sum_propensity=0.0;
    int reac=0,i;
    printf("\n");
printf("Reaction 1: A + B-> C \n");
printf("Reaction 2: C + D-> E\n");
printf("Reaction 3: B + E-> F \n");
printf("Reaction 4: C + F-> G\n");
printf("\n\n");
printf("Enter the initial count for species A\n");
printf("Enter the initial count for species B\n");
printf("Enter the initial count for species C\n");
printf("Enter the initial count for species D\n");
printf("Enter the initial count for species E\n");
printf("Enter the initial count for species F\n");
printf("Enter the initial count for species G\n");
//printf("Enter the initial count for species H\n");

for( i=0;i<n;i++)
scanf("%d",&x[i]);

inti(sto);
srand(time(NULL));
while(t<tend){
        tau=0.0;
        double r1=(double)rand()/RAND_MAX;
       // printf("\n");
        //printf("%f\n",r1);
double r2=(double)rand()/RAND_MAX;
//printf("%f\n",r2);
calpropen(x);
sum_propensity=sumpropen(p);
display(x,t,reac+1);
tau=((-log(r1))/sum_propensity);

reac=select_rxn(p,m,sum_propensity,r2);
update_count(x,reac);
t+=tau;
//ch=getch();
sleep(3);
}
return 0;
}








