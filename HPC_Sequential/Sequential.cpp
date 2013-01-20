#include<iostream>
#include<cmath>
#include<fstream>
using namespace std;
#define N 32
#define MAXITER 1000000
#define TOL 0.000001

void printPlotData(double phi[N+1][N+1],double sol[N+1][N+1],double x[N+1],double y[N+1], int iteration, double error, ofstream& plotData);
double getAnalyticalSolution(double x, double y, int precision, const double PI);
double getFourierCoeficient(int n, const double PI);

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
	ofstream plotData("plot_data.txt", ios_base::trunc);

	//Analytical Solution
	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
			sol[i][j]=getAnalyticalSolution(x[i],y[j],100,PI);

	//Jacobi
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
		errorLog<<"Iteration "<<k<<": "<<iteration_error<<endl;

		if(iteration_error<=TOL){
			printPlotData(phi,sol,x,y,k,iteration_error,plotData);
			iterate=false;
		}
		k++;
	}
	
	//system("pause");
	return 0;
}

double getAnalyticalSolution(double x, double y, int precision, const double PI){
	//first summation element
	double alpha=getFourierCoeficient(1,PI);
	double u=alpha*sin(PI*y)*(cosh(PI*x)-sinh(PI*x));
	if(precision>2){	//element 2 introduces a NaN (division by zero)
		for(int n=3; n<=precision; n++){
			alpha=getFourierCoeficient(n,PI);
			u+=alpha*sin(n*PI*y)*(cosh(n*PI*x)-sinh(n*PI*x));
		}
	}
	return u;
}

double getFourierCoeficient(int n, const double PI){
	double num=4*(-1+cos(n*PI));
	double den=(-4*n+pow(n,3))*PI;
	return num/den;
}

void printPlotData(double phi[N+1][N+1],double sol[N+1][N+1],double x[N+1],double y[N+1], int iteration, double error, ofstream& out){
	out<<"N = "<<N<<"\nRequired iterations: "<<iteration<<endl;
	for(int i=0; i<=N; i++)
		for(int j=0; j<=N; j++)
			out<<"("<<x[i]<<","<<y[j]<<") - "<<"Estimated solution: "<<phi[i][j]<<" Analytic solution: "<<sol[i][j]<<endl;
	

}