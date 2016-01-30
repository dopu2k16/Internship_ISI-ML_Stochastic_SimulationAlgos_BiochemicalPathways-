#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<limits.h>
#include<time.h>
#define n 20//no. of species
#define m 18
#define rc 0.1//np. of reactions
#define tend 1
int sto[n][m], x[n]; //Stochiometric matrix & species matrix
double p[n];



void inti(int sto[n][m]){
    int i,j;
for(i=0;i<n;i++){
    for(j=0;j<m;j++){
    sto[i][j]=0;}}

sto[0][0]=-1;
sto[0][1]=-1;

sto[1][0]=-1;
sto[1][6]=-1;

sto[2][0]=1;
sto[2][1]=-1;
sto[2][2]=-1;
sto[2][14]=-1;

sto[3][1]=1;
sto[3][2]=-1;
sto[3][4]=-1;
sto[3][5]=-1;

sto[4][2]=1;
sto[4][3]=-1;

sto[5][3]=-1;
sto[5][6]=-1;
sto[5][7]=-1;

sto[6][3]=1;
sto[6][4]=-1;
sto[6][5]=-1;
sto[6][7]=-1;
sto[6][8]=-1;
sto[7][4]=1;
sto[7][10]=-1;

sto[8][7]=1;
sto[8][8]=-1;
sto[8][10]=-1;
sto[8][12]=-1;

sto[9][10]=1;
sto[9][11]=-1;
sto[9][12]=-1;
sto[9][13]=-1;
sto[10][5]=1;
sto[10][13]=-1;
sto[11][12]=1;
sto[11][17]=-1;

sto[12][6]=1;
sto[12][9]=-1;
sto[12][14]=-1;

sto[13][9]=-1;
sto[13][11]=-1;
sto[14][8]=1;
sto[14][11]=1;
sto[14][16]=-1;
sto[14][17]=-1;

sto[15][9]=1;
sto[15][15]=-1;
sto[15][16]=-1;

sto[16][15]=1;

sto[17][16]=1;

sto[18][14]=1;
sto[18][17]=1;

sto[19][13]=1;

}

//Calculate propensity values for each
void calpropen(int x[]){
p[0]=rc*x[0]*x[1];
p[1]=rc*x[0]*x[2];
p[2]=rc*x[2]*x[3];
p[3]=rc*x[4]*x[5];
p[4]=rc*x[3]*x[6];
p[5]=rc*x[3]*x[6];
p[6]=rc*x[1]*x[5];
p[7]=rc*x[5]*x[6];
p[8]=rc*x[5]*x[8];
p[9]=rc*x[12]*x[13];
p[10]=rc*x[7]*x[8];
p[11]=rc*x[13]*x[9];
p[12]=rc*x[8]*x[9];
p[13]=rc*x[10]*x[9];
p[14]=rc*x[2]*x[12];
p[15]=rc*x[0]*x[15];
p[16]=rc*x[15]*x[14];
p[17]=rc*x[11]*x[14];

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
printf("\nAT Time=%f ",time);
printf("\n");
if(time!=0.0)
printf("\t Reaction No : %d\n",rec+1);

printf("\n");
printf("\t Individual Species Count:\n");
printf("\n");
printf("A\tB\tC\tD\tE\tF\tG\tH\tI\tJ\tK\tL\tM\tN\tO\tP\tQ\tR\tS\tT\n");
printf("\n");
printf("%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n ",x[0],x[1],x[2],x[3],x[4],x[5],x[6],x[7],x[8],x[9],x[10],x[11],x[12],x[13],x[14],x[15],x[16],x[17],x[18],x[19]);
printf("\n\n\n");
}

int main(){
    double t=0.0,tau,sum_propensity=0.0;
    int reac=0,i;
//char ch;
printf("\n");
printf("Reaction 1: A + B-> C\n");
printf("Reaction 2: A + C-> D\n");
printf("Reaction 3: C + D-> E\n");
printf("Reaction 4: E + F-> G\n");
printf("Reaction 5: D + G-> H\n");
printf("Reaction 6: D + G-> K\n");
printf("Reaction 7: B + F-> M\n");
printf("Reaction 8: F + G-> J\n");
printf("Reaction 9: F + I-> N\n");
printf("Reaction 10: M + N-> P\n");
printf("Reaction 11: H + I-> J\n");
printf("Reaction 12: N + J-> O\n");
printf("Reaction 13: I + J-> L\n");
printf("Reaction 14: K + J-> T\n");
printf("Reaction 15: C + M-> S\n");
printf("Reaction 16: A + P-> Q\n");
printf("Reaction 17: P + O-> R\n");
printf("Reaction 18: L + O-> S\n");
printf("\n");
printf("Enter the initial count for species A\n");
printf("Enter the initial count for species B\n");
printf("Enter the initial count for species C\n");
printf("Enter the initial count for species D\n");
printf("Enter the initial count for species E\n");
printf("Enter the initial count for species F\n");
printf("Enter the initial count for species G\n");
printf("Enter the initial count for species H\n");
printf("Enter the initial count for species I\n");
printf("Enter the initial count for species J\n");
printf("Enter the initial count for species K\n");
printf("Enter the initial count for species L\n");
printf("Enter the initial count for species M\n");
printf("Enter the initial count for species N\n");
printf("Enter the initial count for species O\n");
printf("Enter the initial count for species P\n");
printf("Enter the initial count for species Q\n");
printf("Enter the initial count for species R\n");
printf("Enter the initial count for species S\n");
printf("Enter the initial count for species T\n");
//printf("\n\n");

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
display(x,t,reac);
tau=((-log(r1))/sum_propensity);

reac=select_rxn(p,m,sum_propensity,r2);
update_count(x,reac,sto);
t+=tau;
sleep(5);
//ch=getch();
}
return 0;
}








