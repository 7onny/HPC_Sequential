#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;
#define N 4
#define MAXITER 1000000
#define TOL 0.000001
int main(){
	//Initializing data structures
	double const PI=3.14159265;
	double phi[N+1][N+1],old[N+1][N+1],x[N+1],y[N+1],sol[N+1][N+1];
	for(int i=0; i<=N; i++){
		x[i]=(double)i/N;
		y[i]=(double)i/N;
	}
	for(int i=0; i<=N; i++){
		for(int j=0; j<=N; j++){
			if(i==0){
				phi[i][j]=(sin(PI*y[j])*sin(PI*y[j]));
			}
			else{
				phi[i][j]=0;
			}
			old[i][j]=phi[i][j];
		}
	}

	//File outputs
	ofstream errorLog("error_log.txt",ios_base::trunc);
	ofstream plot("plot_data", ios_base::trunc);

	//
	int k=0;
	bool iterate=true;
	while(iterate && k<MAXITER){
		double iteration_error=0.0;
		for(int i=1; i<N; i++)
			for(int j=1; j<N; j++){
				phi[i][j]=0.25*(old[i-1][j]+old[i+1][j]+old[i][j-1]+old[i][j+1]);
				iteration_error+=fabs(phi[i][j]-old[i][j]);
			}
		for(int i=1;i<N;i++)
			for(int j=1;j<N;j++)
				old[i][j]=phi[i][j];
		if(iteration_error<=TOL)
			iterate=false;
		k++;
	}
	cout<<x[0]<<","<<y[1]<<" = "<<phi[0][1]<<endl;
	system("pause");
	return 0;
}