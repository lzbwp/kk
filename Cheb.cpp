#include <cstdlib>
#include <iostream>
#include <cmath>
#include <iomanip>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <complex>


#define GAP 0.05


#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define NR_END 1
#define FREE_ARG char*


static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)

using namespace std;
typedef complex<double> COMPLEX;
typedef vector<COMPLEX> dvec;
vector<vector<int> > basis;
double eps = pow (2.0, -49.0);
double PI=4.0*atan(1.0), QTF=0.5;
double ALPHA, T_X = 1, T_Y = 1, K_X, K_Y ;
int NUM=100,B=0.0;

double A=6.0;
vector<double> a, b;
const COMPLEX II(0.0,1.0);


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}




double pythag(double a, double b)
{
	double absa,absb;
	absa=fabs(a);
	absb=fabs(b);
	if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void read_basis(int n){
	int size, entry, i, j;
	vector<int> temp;
	basis.clear();
	ostringstream file;
	string s;
	ifstream myfile;
	file << "basis/basis" << n << ".txt";
	s = file.str();
	myfile.open(s.c_str());
	myfile >> size;
	for(i = 0; i < size; i++){
		for(j = 0; j < 5; j++){
			myfile >> entry;
			temp.push_back(entry);
		}
		basis.push_back(temp);
		temp.clear();
	}
	myfile.close();
}


double norm(dvec phi){
	double norm = 0;
	int i;
	for(i = 0; i < phi.size(); i++){
		norm += abs(phi[i])*abs(phi[i]);
	}
	
	return norm;
}

dvec hphi(dvec phi){
	dvec hphi;
	COMPLEX z(0.0,0.0);
	int i,ev_odd,Ri;
	double distance;
	for(i = 0; i < basis.size(); i++)
	{
		
		ev_odd = (basis[i][3]+basis[i][4])&1;
		if(ev_odd==1) //B site
		{
			
			distance= sqrt((basis[i][3]*1.5+0.5)*(basis[i][3]*1.5+0.5)+(basis[i][4]*0.5*sqrt(3.0))*(basis[i][4]*0.5*sqrt(3.0)));
		   /* if(i<4)
			{
				Ri=1;
			}
			else
			{
				Ri=0;
			}*/
		
		//z = (-GAP-Ri*0.5*ALPHA)*phi[i];
		z = (0.0,0.0)*phi[i];
		if(basis[i][0]+1){z-=phi[basis[i][0]]*T_X;}
		if(basis[i][1]+1){z-=phi[basis[i][1]]*T_Y;}
		if(basis[i][2]+1){z-=phi[basis[i][2]]*T_Y;}
		}
		else // A site
		{	
			distance= sqrt((basis[i][3]*1.5+1)*(basis[i][3]*1.5+1)+(basis[i][4]*0.5*sqrt(3.0))*(basis[i][4]*0.5*sqrt(3.0)));
		   /* if(i<1)
			{
				Ri=1;
			}
			else
			{
				Ri=0;
			}
			if(i==0)
			{
				Ri=1;
			}
			else
			{
				Ri=0;
			}*/
	//	z = (GAP-ALPHA/distance*exp(-QTF*distance))*phi[i];
	
		z = (0.0,0.0)*phi[i];
		if(basis[i][0]+1){z-=phi[basis[i][0]]*T_X;}
		if(basis[i][1]+1){z-=phi[basis[i][1]]*T_Y;}
		if(basis[i][2]+1){z-=phi[basis[i][2]]*T_Y;}
		}	

		hphi.push_back(z/A);
	}
	return hphi;
}

double inner(dvec phi1, dvec phi2){
	COMPLEX c1(0,0);
	double c2;
	int i;
	for(i = 0; i < phi1.size(); i++){
		c1 += conj(phi1[i]) * phi2[i];
	}
	c2 = real(c1);
	return c2;
}



void expand(dvec &phi)
{
	double t;
	COMPLEX z;
        while(phi.size() < basis.size())
		{
				//t = (double)rand()/32768.;
				
				z = COMPLEX(0.0,0.0);
			
                phi.push_back(z);
        }
}

void initialize(dvec &phin, dvec &phinm)
{
	int i;
	double inne,nm,nn;
	ofstream file;
	phinm.clear(); 
	a.clear(); b.clear();
	phinm.push_back(1);
	//phinm.push_back(1);
	expand(phinm);
	nm = norm(phinm);
	if(nm!=1)
	{
	for(i = 0; i < phinm.size(); i++){
		phinm[i] = phinm[i]/sqrt(nm);
	}
	}
    a.push_back(1);
	phin = hphi(phinm);
	
	inne=inner(phinm, phin);
    a.push_back(inne);
	nn = norm(phin);
	a.push_back(2.0*nn-a[0]);

	
}


int main(){
	dvec phinm, phin, phinp,coef;
	double nn, inne,energy,G,x,epsilon=0.001;
	int  i,j,Nt;
		

	ostringstream filename;
	string sfile;
	ofstream outfile;
	outfile.precision(16);
	cout.precision(16);
	 	
	phinm.clear(); phin.clear(); phinp.clear(); 

	read_basis(NUM);

	for(ALPHA = 0.0; ALPHA < 0.5; ALPHA += 0.5)
	{
		filename << "data/cheb=" << NUM+3000 << "_ALPHA = " << ALPHA << "_ab.txt";
		sfile = filename.str();
		filename.str("");
		outfile.open(sfile.c_str(), ios::app);

		K_X = 0.0; 	K_Y = 0.0;
			
		initialize(phin, phinm);
		
		for(Nt = 1; Nt<NUM+3000; Nt++)
		{
		
						
			phinp = hphi(phin);
		
			
			for(i = 0; i < basis.size(); i++)
			{
			
				phinp[i] = 2.0*phinp[i] - phinm[i];
				phinm[i] = phin[i];
				phin[i] = phinp[i];
				 
			}
			inne = inner(phin,phinm);
            nn=norm(phin);

			a.push_back(2.0*inne-a[1]);
			a.push_back(2.0*nn-a[0]);
		outfile << "a[" << Nt-1 << "]=" << a[Nt-1] << endl;
    			
		}
		outfile.close();
		filename << "data/2010/cheb=" << NUM+3000 << "_ALPHA = " << ALPHA << ".txt";
		sfile = filename.str();
		filename.str("");
		outfile.open(sfile.c_str(), ios::app);

		

	for(energy = -3.5; energy < 3.5; energy+=0.001)
	{	
		x=(energy-B)/A;
	    
        G=a[0];
		for (j=1;j< a.size();j++)
		{
			G+= 2*a[j]*(exp(-epsilon*j)-exp(epsilon*j)*exp(-2*epsilon*a.size()))/(1-exp(-2*epsilon*a.size()))*cos(j*acos(x));
			 //G+= 2*a[j]*sinh(2*(1-j/a.size()))/sinh(2)*cos(j*acos(x)); //Lorentz kernel
			 //G+= 2*cos(j*acos(x))*a[j]*((a.size()-j+1)*cos(PI*j/(a.size()+1))+sin(PI*j/(a.size()+1))/tan(PI/(a.size()+1)))/(a.size()+1); //Jackson kernel
		}
		
	  G=G/sqrt(1-x*x)/PI/A;
	  outfile << energy << " " << G << endl;
		
	}

	outfile.close();

	}
 
	
	return 0;
}
