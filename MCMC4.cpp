#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <vector>
#include <iostream>
#include<time.h> 

//g++-10 MCMC4.cpp -o MCMC4.e -lm
using namespace std;

//Function that give back a random number drawn in function of a gaussian with a mean and a sigma given by the input. 
double random_gaussienne(double mean, double sigma){
	bool condi=false;
	double min_x=mean-10*sigma;
	while (condi==false){
		double x= min_x+rand()/(RAND_MAX+1.0)*20*sigma; 
		double y=rand()/(RAND_MAX+1.0);
		double f_x=1/(sigma*sqrt(2*M_PI))*exp(-pow((x-mean),2)/(2*pow(sigma,2)));
		if (y<f_x){
			condi=true;
			return x; 
		}
	}
}

//Function that give back the factorial of the number given in input
double facto(int a){ 
	double fact=1;
	for (int c = 1; c <= a; c++){fact = fact * c;}
	return fact;
}


//Function that gives back the likelihood of the data given in input (pourcentage) given a defined model parametrised by the input vector param and the number of dwarf ingested in the survey to obtain the data Ngalax. 
double poisson(vector< vector< double> >pourcentage,vector< double> param, int Ngalax, double (&facto)(int)){
	double L=1;
	double mv[16]= {-4.5,-4.75,-5,-5.25,-5.5,-5.75,-6,-6.25,-6.5,-6.75,-7,-7.25,-7.5,-7.75,-8,-8.25};// grid of magnitude in the V band
	double rh[12]={1.8,1.9,2,2.1,2.2,2.3,2.4,2.5,2.6,2.7,2.8,2.9};// grid of log(rh)
	vector <vector <double> > model(12, vector<double>(16));
	double rho=1.8-param[0]*param[1];
	//For each mv bin and rh bin the model is calculated:
	for (int m=0;m<12;m++){
		for (int n=0;n<16;n++){
			float eff=0;
			// In order to be more precise the recovery of dwarf galaxies, for each mv-log(rh) bin is calculated thanks to an integral over the bin:
			for (int o=0;o<10;o++){
				for (int p=0;p<10;p++){
					double mv_limm=((rh[m]+p*0.01-rho)/param[1]);
					double beta=erfc((mv[n]-o*0.025-mv_limm)/(sqrt(2)*param[2]));
					eff=eff+0.5*Ngalax*beta;
				}
			}
		
			model[m][n]=eff/100.0;
			L*=pow(model[m][n],pourcentage[m][n])/facto(pourcentage[m][n])*exp(-model[m][n]);
			
		}
	}
	return L;
}

int main(int argc, const char * argv[]){

int champ=1; //the number of the field in which we want to fit the recovery fractions
char line1[10000];
//Obtention of the calculated recovery fractions for the given filed
vector <vector <double> > data(12, vector<double>(16));
char entree[1000];
sprintf(entree,"DL/DL_champ%d.dat", champ); //Input file with the calculated recovery fraction such as in Doliva-Dolinsky and al, 2022
FILE* file_percentage=fopen(entree,"r");
for(int i=0; i<12; i++){
	fgets(line1,10000,file_percentage);
	char *ligne=line1;
	char *strToken =strtok(ligne," ");
	data[i][0]=atof(strToken);
	for (int j=1; j<16;j++){
		strToken = strtok ( NULL, " " );
		data[i][j]=atof(strToken);
	}

}
fclose(file_percentage);
//Output file:
FILE * file_sortie;
char sortie[1000];
sprintf(sortie,"resultat/sortie_%d.txt", champ);
file_sortie=fopen(sortie,"w");

srand(time(NULL)); 
double Ngal=5; // Number of dwarf galaxies added to the survey to calculate the recovery fractions
vector< double> erreur={0.1,0.05,0.05}; // Vector that contains the error
vector< vector <double>> parametres={{-5.4,-0.5,0.4}}; // Vector that contains a guess of the parameters value
vector<double> L;
L.push_back(poisson(data,parametres[0],Ngal,facto));
int number_points_accepted=0; //Number of accepted points
int number_points=0; //Total number of points
int g;
double L_tmp;
double u_rand;
bool condition;
vector <double> p_tmp;

for (int i=1;i<20000;i++){
	number_points+=1;
	g=parametres.size()-1;
	p_tmp={random_gaussienne(parametres[i-1][0],erreur[0]),random_gaussienne(parametres[i-1][1],erreur[1]),random_gaussienne(parametres[i-1][2],erreur[2])}; //We draw new parameters thanks following a gaussian function
	L_tmp=poisson(data,p_tmp,Ngal,facto); //We calculate the new likelihood
	u_rand=rand()/(RAND_MAX+1.0);
	if (u_rand<L_tmp/L[g]){ //If the point is accepted 
		parametres.push_back(p_tmp);
		L.push_back(L_tmp);
		number_points_accepted+=1;
		fprintf(file_sortie,"%f %f %f %e \n",parametres[g+1][0],parametres[g+1][1],parametres[g+1][2],L_tmp);
		condition=true;
		cout<<number_points_accepted<<" "<<number_points<<" "<<g<<" "<<u_rand<<" "<<L_tmp<<" "<<L[g]<<endl;
		cout<<parametres[g+1][0]<<" "<<parametres[g+1][1]<<" "<<parametres[g+1][2]<<endl;
		
	}
	else{//If the point is not accepted 
		parametres.push_back(parametres[g]);
		L.push_back(L[g]);
		cout<<number_points_accepted<<" "<<number_points<<" "<<g<<" "<<u_rand<<" "<<L_tmp<<" "<<L[g]<<endl;
		cout<<parametres[g][0]<<" "<<parametres[g][1]<<" "<<parametres[g][2]<<endl;
		fprintf(file_sortie,"%f %f %f %e \n",parametres[g][0],parametres[g][1],parametres[g][2],L[g]);
	}
}

}
