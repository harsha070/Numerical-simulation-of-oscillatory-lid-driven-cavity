#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;
void update_coefficients();
void correct_pressure_and_velocities();
void initialize();
void solvePressure();
void solveVelocities();
void update_coefficientsP();

const int N=72;

double L=0.6;
double dx=L/N;

double omega=M_PI/6;
double dt=0.024;

double H=0.3;
double dy=L/N;

double mu=0.001;
double density=1000;

double u_star[2*N+1][2*N+1];
double u_prime[2*N+1][2*N+1];
double u_prev[2*N+1][2*N+1];
double u_prev_time[2*N+1][2*N+1];

double v_star[2*N+1][2*N+1];
double v_prime[2*N+1][2*N+1];
double v_prev[2*N+1][2*N+1];
double v_prev_time[2*N+1][2*N+1];

double P_star[2*N+1][2*N+1];
double P_prime[2*N+1][2*N+1];

double temp[2*N][2*N];

double Fe;
double Fn;
double Fw;
double Fs;

double De;
double Dn;
double Dw;
double Ds;

double aP[2*N+1][2*N+1];
double aE[2*N+1][2*N+1];
double aW[2*N+1][2*N+1];
double aN[2*N+1][2*N+1];
double aS[2*N+1][2*N+1];

double S[2*N+1][2*N+1];
double Sdc[2*N+1][2*N+1];

double PSI[2*N+1][2*N+1];

double error1,error2;

double vel=0.01;


void update_coefficients(){
	int i,j;
	for(i=0;i<2*N+1;i++){
		for(j=0;j<2*N+1;j++){
			PSI[i][j]=1;
		}
	}
	//////////////////////////////// FOR X MOMENTUM INTERNAL NODES ////////////////////////////////
	for(i=3;i<=2*N-3;i=i+2)
	{
		for(j=4;j<=2*N-4;j=j+2)
		{
			Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
			Fn= (v_star[i+1][j+1]+v_star[i+1][j-1])*density/2;
			Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
			Fs= (v_star[i-1][j+1]+v_star[i-1][j-1])*density/2;

			De=mu/dx;
			Dw=mu/dx;
			Dn=mu/dy;
			Ds=mu/dy;

			if(Fn>=0 && Fs>=0){
				if(Fe>=0 && Fw>=0){
					aP[i][j]=Fe*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aP[i][j]=-1*Fw*dy+Fe*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aP[i][j]=Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
				}
				else if(Fe<0 && Fw<0){
					aP[i][j]=-Fw*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
				}
			}
			else if(Fn<0 && Fs>=0){
				if(Fe>=0 && Fw>=0){
					aP[i][j]=Fe*dy+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aP[i][j]=Fe*dy-Fw*dy+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aP[i][j]=De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
				}
				else if(Fe<0 && Fw<0){
					aP[i][j]=-1*Fw*dy+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])+Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
				}
			}
			else if(Fn>=0 && Fs<0){
				if(Fe>=0 && Fw>=0){
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Fn+Dn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
			}
			else if(Fn<0 && Fs<0){
				if(Fe>=0 && Fw>=0){
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
				else if(Fe<0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
				}
			}

		}
	}
	//////////////////////////////// FOR X MOMENTUM NODE (1,2) ////////////////////////////////
	i=1;
	j=2;
	Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
	Fn= (v_star[i+1][j+1]+u_star[i+1][j-1])*density/2;
	Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
	Fs= 0.0;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fe>=0){
		if(Fw>=0 && Fn>=0){
			aP[i][j]=Fe*dy+Fn*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw>=0 && Fn<0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy - (Fn*dx);
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fn>=0){
			aP[i][j]=Fn*dx+Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fn<0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
	}
	else if(Fe<0){
		if(Fw>=0 && Fn>=0){
			aP[i][j]=Fn*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=-1*Fe*dy+mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw>=0 && Fn<0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=-1*Fe*dy+mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fn>=0){
			aP[i][j]=Fn*dx+Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fn<0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
	}
	//////////////////////////////// FOR X MOMENTUM NODE (1,2*N-2) ////////////////////////////////
	i=1;
	j=2*N-2;
	Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
	Fn= (v_star[i+1][j+1]+u_star[i+1][j-1])*density/2;
	Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
	Fs= 0.0;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fe>=0){
		if(Fw>=0 && Fn>=0){
			aP[i][j]=Fe*dy+Fn*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aN[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw>=0 && Fn<0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=Fw*dy+mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fn>=0){
			aP[i][j]=Fn*dx+Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fn<0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
	}
	else if(Fe<0){
		if(Fw>=0 && Fn>=0){
			aP[i][j]=Fn*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aN[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw>=0 && Fn<0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=Fw*dy+mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fn>=0){
			aP[i][j]=Fn*dx+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fn<0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aN[i][j]=mu*dx/dy - Fn*dx;
			aE[i][j]=0;
			aS[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])-Fn*PSI[i+1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
	}
	//////////////////////////////// FOR X MOMENTUM NODE (2*N-1,2) ////////////////////////////////
	i=2*N-1;
	j=2;
	Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
	Fn= 0.0;
	Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
	Fs= (v_star[i-1][j+1]+u_star[i-1][j-1])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fe>=0){
		if(Fw>=0 && Fs>=0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aS[i][j]=mu*dx/dy+Fs*dx;
			aE[i][j]=mu*dy/dx;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw>=0 && Fs<0){
			aP[i][j]=Fe*dy-Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aS[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i-2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fs>=0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aS[i][j]=Fs*dx+mu*dx/dy;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fs<0){
			aP[i][j]=Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx;
			aS[i][j]=mu*dx/dy;
			aN[i][j]=0;
			aW[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
	}
	else if(Fe<0){
		if(Fw>=0 && Fs>=0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aS[i][j]=mu*dx/dy+Fs*dx;
			aE[i][j]=mu*dy/dx - Fe*dy;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw>=0 && Fs<0){
			aP[i][j]=-1*Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aE[i][j]=-1*Fe*dy+mu*dy/dx;
			aS[i][j]=mu*dx/dy;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i-2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fs>=0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx-Fe*dy;
			aS[i][j]=Fs*dx+mu*dx/dy;
			aW[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fs<0){
			aP[i][j]=-1*Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aE[i][j]=mu*dy/dx - Fe*dy;
			aS[i][j]=mu*dx/dy;
			aN[i][j]=0;
			aW[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
	}
	//////////////////////////////// FOR X MOMENTUM NODE (2*N-1,2*N-2) ////////////////////////////////
	i=2*N-1;
	j=2*N-2;
	Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
	Fn= 0.0;
	Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
	Fs= (v_star[i-1][j+1]+u_star[i-1][j-1])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fe>=0){
		if(Fw>=0 && Fs>=0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aS[i][j]=mu*dx/dy+Fs*dx;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw>=0 && Fs<0){
			aP[i][j]=Fe*dy-Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aS[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fs>=0){
			aP[i][j]=Fe*dy+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aS[i][j]=Fs*dx+mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
		else if(Fw<0 && Fs<0){
			aP[i][j]=-1*Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+Fe*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aS[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
		}
	}
	else if(Fe<0){
		if(Fw>=0 && Fs>=0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aS[i][j]=mu*dx/dy+Fs*dx;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw>=0 && Fs<0){
			aP[i][j]=-1*Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx+Fw*dy;
			aS[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j]-u_star[i][j-2])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fs>=0){
			aP[i][j]=(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aS[i][j]=Fs*dx+mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i-2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
		else if(Fw<0 && Fs<0){
			aP[i][j]=-1*Fs*dx+(3*mu*dx/dy)+(2*mu*dy/dx)-Fw*dy+density/dt*dx*dy;
			aW[i][j]=mu*dy/dx;
			aS[i][j]=mu*dx/dy;
			aE[i][j]=0;
			aN[i][j]=0;
			S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*vel*dx)+density/dt*dx*dy*u_prev_time[i][j];
			Sdc[i][j]=0.5*(Fw*PSI[i][j-1]*dy*(u_star[i][j-2]-u_star[i][j])+Fs*PSI[i-1][j]*dx*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j]-u_star[i][j+2]));
		}
	}
	//////////////////////////////// FOR X MOMENTUM LEFT EDGE (........,2) ////////////////////////////////
	j=2;
	for(i=3;i<=2*N-3;i=i+2){
		Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
		Fn= (v_star[i+1][j+1]+v_star[i+1][j-1])*density/2;
		Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
		Fs= (u_star[i-1][j+1]+u_star[i-1][j-1])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		// left edge
		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR X MOMENTUM RIGHT EDGE (........,2*N-2) ////////////////////////////////
	j=2*N-2;
	for(i=3;i<=2*N-3;i=i+2){
		Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
		Fn= (v_star[i+1][j+1]+u_star[i+1][j-1])*density/2;
		Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
		Fs= (u_star[i-1][j+1]+u_star[i-1][j-1])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR X MOMENTUM BOTTOM EDGE (1,........) ////////////////////////////////
	i=1;
	for(j=4;j<=2*N-4;j=j+2){
		Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
		Fn= (v_star[i+1][j+1]+u_star[i+1][j-1])*density/2;
		Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
		Fs=0.0;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
		        aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i][j]-u_star[i+2][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fn*dx*PSI[i+1][j]*(u_star[i+2][j]-u_star[i][j])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR X MOMENTUM TOP EDGE (2*N-1,........) ////////////////////////////////
	i=2*N-1;
	for(j=4;j<=2*N-4;j=j+2){
		Fe= (u_star[i][j]+u_star[i][j+2])*density/2;
		Fn=0.0;
		Fw= (u_star[i][j-2]+u_star[i][j])*density/2;
		Fs= (v_star[i-1][j-1]+u_star[i-1][j+1])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		// top edge
		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])-Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])+Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i][j-1]-P_star[i][j+1])*dy+(2*Dn*dx-Fn*dx)*vel+density/dt*dx*dy*u_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(u_star[i][j]-u_star[i-2][j])-Fw*dy*PSI[i][j-1]*(u_star[i][j]-u_star[i][j-2])+Fe*dy*PSI[i][j+1]*(u_star[i][j+2]-u_star[i][j]));
			}
		}
	}



	//////////////////////////////// Y MOMENTUM COEFFICIENTS BEGIN //////////////////////////////////////////

	//////////////////////////////// FOR Y MOMENTUM INTERNAL NODES ////////////////////////////////
	for(i=4;i<=2*N-4;i=i+2)
	{
		for(j=3;j<=2*N-3;j=j+2)
		{
			Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
			Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
			Fw= (u_star[i+1][j-1]+u_star[i-1][j-1])*density/2;
			Fs= (v_star[i-2][j]+v_star[i][j])*density/2;

			De=mu/dx;
			Dw=mu/dx;
			Dn=mu/dy;
			Ds=mu/dy;

			if(Fn>=0 && Fs>=0){
				if(Fe>=0 && Fw>=0){
					aP[i][j]=Fe*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j]-v_star[i][j-2])-Fn*PSI[i+1][j]*dx*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aP[i][j]=-1*Fw*dy+Fe*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j-2]-v_star[i][j])-Fn*PSI[i+1][j]*dx*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aP[i][j]=Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j]-v_star[i][j-2])-Fn*PSI[i+1][j]*dx*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j]-v_star[i][j+2]));
				}
				else if(Fe<0 && Fw<0){
					aP[i][j]=-Fw*dy+Fn*dx+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j-2]-v_star[i][j])-Fn*PSI[i+1][j]*dx*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j]-v_star[i][j+2]));
				}
			}
			else if(Fn<0 && Fs>=0){
				if(Fe>=0 && Fw>=0){
					aP[i][j]=Fe*dy+De*dy+Dw*dy+Dn*dx+Ds*dx;
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j]-v_star[i][j-2])-Fn*PSI[i+1][j]*dx*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aP[i][j]=Fe*dy-Fw*dy+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j-2]-v_star[i][j])-Fn*PSI[i+1][j]*dx*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aP[i][j]=De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j]-v_star[i][j-2])-Fn*PSI[i+1][j]*dx*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j]-v_star[i][j+2]));
				}
				else if(Fe<0 && Fw<0){
					aP[i][j]=-1*Fw*dy+De*dy+Dw*dy+Dn*dx+Ds*dx+density/dt*dx*dy;
					aW[i][j]=Dw*dy;
					aS[i][j]=(Fs+Ds)*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(Fs*PSI[i-1][j]*dx*(v_star[i][j]-v_star[i-2][j])+Fw*PSI[i][j-1]*dy*(v_star[i][j-2]-v_star[i][j])-Fn*PSI[i+1][j]*dx*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j]-v_star[i][j+2]));
				}
			}
			else if(Fn>=0 && Fs<0){
				if(Fe>=0 && Fw>=0){
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=Dn*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
			}
			else if(Fn<0 && Fs<0){
				if(Fe>=0 && Fw>=0){
					aW[i][j]=(Fw+Dw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe>=0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=De*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw>=0){
					aW[i][j]=(Dw+Fw)*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
				else if(Fe<0 && Fw<0){
					aW[i][j]=Dw*dy;
					aS[i][j]=Ds*dx;
					aE[i][j]=(De-Fe)*dy;
					aN[i][j]=(Dn-Fn)*dx;
					aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
					S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
					Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
				}
			}

		}
	}
	//////////////////////////////// FOR Y MOMENTUM NODE (2,1) ////////////////////////////////
	i=2;
	j=1;
	Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
	Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
	Fw= 0.0;
	Fs=(v_star[i][j]+v_star[i-2][j])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=0.0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	//////////////////////////////// FOR Y MOMENTUM NODE (2,2*N-1) ////////////////////////////////
	i=2;
	j=2*N-1;
	Fe= 0.0;
	Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
	Fw= (u_star[i+1][j-1]+u_star[i-1][j-1])*density/2;
	Fs= (v_star[i-2][j]+v_star[i][j])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0.0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0.0;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}


	//////////////////////////////// FOR Y MOMENTUM NODE (2*N-2,1) ////////////////////////////////
	i=2*N-2;
	j=1;
	Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
	Fn= (v_star[i][j]+v_star[i+2][j])*density/2;
	Fw= 0.0;
	Fs= (v_star[i-2][j]+v_star[i][j])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
	}
	if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
	}
	if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
	}
	if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
	}

	//////////////////////////////// FOR Y MOMENTUM NODE (2*N-2,2*N-1) ////////////////////////////////
	i=2*N-2;
	j=2*N-1;

	Fe= 0.0;
	Fn= (v_star[i][j]+v_star[i+2][j])*density/2;
	Fw= (u_star[i+1][j-1]+u_star[i-1][j-1])*density/2;
	Fs= (v_star[i-2][j]+v_star[i][j])*density/2;

	De=mu/dx;
	Dw=mu/dx;
	Dn=mu/dy;
	Ds=mu/dy;

	if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=0.0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}

	//////////////////////////////// FOR Y MOMENTUM LEFT EDGE (........,1) ////////////////////////////////
	j=1;
	for(i=4;i<=2*N-4;i=i+2){
		Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
		Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
		Fw= 0.0;
		Fs= (v_star[i-2][j]+v_star[i][j])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		// left edge
		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=0;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*Dw*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR Y MOMENTUM RIGHT EDGE (........,2*N-1) ////////////////////////////////
	j=2*N-1;
	for(i=4;i<=2*N-4;i=i+2){
		Fe= 0.0;
		Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
		Fw= (u_star[i+1][j-1]+u_star[i-1][j-1])*density/2;
		Fs= (v_star[i][j]+v_star[i-2][j])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy+Fe*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=Dw*dy+Fw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=0;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+2*De*dy-Fw*dy+Fn*dx-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR Y MOMENTUM BOTTOM EDGE (2,........) ////////////////////////////////
	i=2;
	for(j=3;j<=2*N-3;j=j+2){
		Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
		Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
		Fw= (v_star[i-1][j-1]+v_star[i+1][j-1])*density/2;
		Fs=(v_star[i][j]+v_star[i-2][j])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
		        aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Fn*dx+Dn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=Dn*dx;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i][j]-v_star[i+2][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=De*dy;
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=0;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=(Dn*dx-Fn*dx);
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Ds*dx+Fe*dy-Fw*dy+Fn*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fn*dx*PSI[i+1][j]*(v_star[i+2][j]-v_star[i][j])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	}
	//////////////////////////////// FOR Y MOMENTUM TOP EDGE (2*N-2,........) ////////////////////////////////
	i=2*N-2;
	for(j=3;j<=2*N-3;j=j+2){
		Fe= (u_star[i+1][j+1]+u_star[i-1][j+1])*density/2;
		Fn= (v_star[i+2][j]+v_star[i][j])*density/2;
		Fw= (v_star[i-1][j-1]+v_star[i+1][j-1])*density/2;
		Fs=(v_star[i][j]+v_star[i-2][j])*density/2;

		De=mu/dx;
		Dw=mu/dx;
		Dn=mu/dy;
		Ds=mu/dy;

		// top edge
		if(Fn*dx>=0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx>=0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=(Fs*dx+Ds*dx);
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx>=0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
		if(Fn*dx<0 && Fs*dx<0){
			if(Fe*dy>=0 && Fw*dy>=0){
				aW[i][j]=(Fw*dy+Dw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy>=0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=De*dy;
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])-Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy>=0){
				aW[i][j]=(Dw*dy+Fw*dy);
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])+Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
			else if(Fe*dy<0 && Fw*dy<0){
				aW[i][j]=Dw*dy;
				aS[i][j]=Ds*dx;
				aE[i][j]=(De*dy-Fe*dy);
				aN[i][j]=0;
				aP[i][j]=aE[i][j]+aW[i][j]+aN[i][j]+aS[i][j]+Dn*dx+Fe*dy-Fw*dy-Fs*dx+density/dt*dx*dy;
				S[i][j]=(P_star[i-1][j]-P_star[i+1][j])*dx+density/dt*dx*dy*v_prev_time[i][j];
				Sdc[i][j]=0.5*(-1*Fs*dx*PSI[i-1][j]*(v_star[i][j]-v_star[i-2][j])-Fw*dy*PSI[i][j-1]*(v_star[i][j]-v_star[i][j-2])+Fe*dy*PSI[i][j+1]*(v_star[i][j+2]-v_star[i][j]));
			}
		}
	}




	}

	void update_coefficientsP()
{int i,j;
	//////////////////////////////////////	PRESSURE CORRECTION BEGINS ///////////////////////////////////////////

	// Internal Nodes //
	for(i=3;i<=2*N-3;i=i+2){
		for(j=3;j<=2*N-3;j=j+2){
			aW[i][j]=dy*dy/aP[i][j-1];
			aE[i][j]=dy*dy/aP[i][j+1];
			aS[i][j]=dx*dx/aP[i-1][j];
			aN[i][j]=dx*dx/aP[i+1][j];
			aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
			S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;
		}
	}

	// Bottom Edge //
	i=1;
	for(j=3;j<=2*N-3;j=j+2){
		aW[i][j]=dy*dy/aP[i][j-1];
		aE[i][j]=dy*dy/aP[i][j+1];
		aS[i][j]=0.0;
		aN[i][j]=dx*dx/aP[i+1][j];
		aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
		S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;
	}

	// Top Edge //
	i=2*N-1;
	for(j=3;j<=2*N-3;j=j+2){
		aW[i][j]=dy*dy/aP[i][j-1];
		aE[i][j]=dy*dy/aP[i][j+1];
		aN[i][j]=0.0;
		aS[i][j]=dx*dx/aP[i-1][j];
		aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
		S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;
	}

	// Left Edge //
	j=1;
	for(i=3;i<=2*N-3;i=i+2){
		aN[i][j]=dx*dx/aP[i+1][j];
		aE[i][j]=dy*dy/aP[i][j+1];
		aS[i][j]=dx*dx/aP[i-1][j];
		aW[i][j]=0.0;
		aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
		S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;
	}

	// Right Edge //
	j=2*N-1;
	for(i=3;i<=2*N-3;i=i+2){
		aW[i][j]=dy*dy/aP[i][j-1];
		aN[i][j]=dy*dy/aP[i+1][j];
		aS[i][j]=dx*dx/aP[i-1][j];
		aE[i][j]=0.0;
		aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
		S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;
	}

	// Left Bottom //
	i=1;
	j=1;
	aW[i][j]=0.0;
	aE[i][j]=dy*dy/aP[i][j+1];
	aS[i][j]=0.0;
	aN[i][j]=dx*dx/aP[i+1][j];
	aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
	S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;

	// Top Left //
	i=2*N-1;
	j=1;
	aW[i][j]=0.0;
	aE[i][j]=dy*dy/aP[i][j+1];
	aN[i][j]=0.0;
	aS[i][j]=dx*dx/aP[i-1][j];
	aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
	S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;

	// Right Bottom //
	j=2*N-1;
	i=1;
	aN[i][j]=dx*dx/aP[i+1][j];
	aE[i][j]=0.0;
	aS[i][j]=dx*dx/aP[i-1][j];
	aW[i][j]=0.0;
	aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
	S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;

	// Top Right //
	j=2*N-1;
	i=2*N-1;
	aW[i][j]=dy*dy/aP[i][j-1];
	aN[i][j]=0.0;
	aS[i][j]=dx*dx/aP[i-1][j];
	aE[i][j]=0.0;
	aP[i][j]=aN[i][j]+aS[i][j]+aE[i][j]+aW[i][j];
	S[i][j]=(u_star[i][j-1]-u_star[i][j+1])*dy+(v_star[i-1][j]-v_star[i+1][j])*dx;

	/////////////////////////////// PRESSURE CORRECTION ENDS /////////////////////////////

}

void correct_pressure_and_velocities(){
	int i,j;
	double alpha_P=0.4;
	double undex_relax=0.4;
	double sum1,sum2,sum3,sum4;
	// correcting pressure
	for(int i=1;i<=2*N-1;i=i+2){
		for(int j=1;j<=2*N-1;j=j+2){
			P_star[i][j]=P_star[i][j]+alpha_P*P_prime[i][j];
		}
	}

	// correcting u
	for(int i=1;i<=2*N-1;i=i+2){
		for(int j=2;j<=2*N-2;j=j+2){

			u_star[i][j]=u_star[i][j]+(P_prime[i][j-1]-P_prime[i][j+1])*(dy/aP[i][j]);
		}
	}

	// correcting v
	for(int i=2;i<=2*N-2;i=i+2){
		for(int j=1;j<=2*N-1;j=j+2){

			v_star[i][j]=v_star[i][j]+(P_prime[i-1][j]-P_prime[i+1][j])*(dx/aP[i][j]);
		}
	}
    sum1=0;
    sum2=0;
    for (i=1;i<=2*N-1;i=i+2)
        {
            for (j=2;j<=2*N-2;j=j+2)
            {
                sum1=fabs(u_star[i][j]-u_prev[i][j])+sum1;
                sum2=fabs(u_prev[i][j])+sum2;
            }
        }
    error1=sum1;
    sum3=0;
    sum4=0;
    for (i=2;i<=2*N-2;i=i+2)
        {
            for (j=1;j<=2*N-1;j=j+2)
            {
                sum3=fabs(v_star[i][j]-v_prev[i][j])+sum1;
                sum4=fabs(v_prev[i][j])+sum2;
            }
        }
    error2=sum3;

    printf("%lf %lf %lf \n",u_star[2*N-1][N],error1,error2);
}

void initialize(){

	for(int j=1;j<=2*N;j=j+2){
		for(int i=1;i<=2*N;i=i+2){
		P_star[i][j]=0.0;
		P_star[i][1]=0.0;
		P_star[i][2*N-1]=0.0;

		P_prime[i][j]=0.0;
	}
	P_star[1][j]=0.0;
	P_star[2*N-1][j]=0.0;
	}


	//Boundary u and v
	for(int j=1;j<=2*N-1;j=j+1){
		u_star[0][j]=0.0;
		v_star[0][j]=0.0;

		u_star[2*N][j]=vel;
		v_star[2*N][j]=0.0;
	}

	for(int i=1;i<=2*N-1;i=i+1){
		u_star[i][0]=0.0;
		v_star[i][0]=0.0;

		u_star[i][2*N]=0.0;
		v_star[i][2*N]=0.0;
	}

	// Initialize u and v
	for(int i=1;i<=2*N-1;i=i+2){
		for(int j=2;j<=2*N-2;j=j+2){
			u_star[i][j]=vel;
			u_prev[i][j]=vel;
		}
	}
	for(int i=2;i<=2*N-2;i=i+2){
		for(int j=1;j<=2*N-1;j=j+2){
			v_star[i][j]=vel/4;
			v_prev[i][j]=vel/4;
		}
	}
}

void solveVelocities()
{
	int i,j;
	double f;
	double error1,error2,sum1,sum2,sum3,sum4;
	double alpha_u=0.4;
	double alpha_v=0.4;
	double a[2*N][2*N],b[2*N][2*N],c[2*N][2*N],d[2*N][2*N],c_dash[2*N][2*N],d_dash[2*N][2*N];

	for(int i=1;i<=2*N-1;i=i+2){
		for(int j=2;j<=2*N-2;j=j+2){
			u_prev[i][j]=u_star[i][j];

		}
	}


	for(int i=2;i<=2*N-2;i=i+2){
		for(int j=1;j<=2*N-1;j=j+2){
			v_prev[i][j]=v_star[i][j];

		}
	}

do{
	for(j=1;j<=2*N-1;)
	{
		if (j==1)
		{
			a[2][j]=-1*alpha_v*aS[2][j];
        	b[2][j]=aP[2][j];
        	c[2][j]=-1*alpha_v*aN[2][j];
        	d[2][j]=alpha_v*(S[2][j]+Sdc[2][j]+aE[2][j]*v_star[2][j+2])+(1-alpha_v)*aP[2][j]*v_prev[2][j];
			c_dash[2][j]=c[2][j]/b[2][j];
	    	d_dash[2][j]=d[2][j]/b[2][j];
			for(i=4;i<=2*N-2;i=i+2)
			{
            	a[i][j]=-1*alpha_v*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-0.5*aN[i][j];
		    	d[i][j]=alpha_v*(S[i][j]+Sdc[i][j]+aE[i][j]*v_star[i][j+2])+(1-alpha_v)*aP[i][j]*v_prev[i][j];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-2][j]=v_star[2*N-2][j];
			v_star[2*N-2][j]=d_dash[2*N-2][j];

			for(i=2*N-4;i>=2;i=i-2)
			{
				temp[i][j]=v_star[i][j];// to store previous iteration values
				v_star[i][j]=d_dash[i][j]-c_dash[i][j]*v_star[i+2][j];
			}
		}
		else if(j==2*N-1)
		{
			a[2][j]=-1*alpha_v*aS[2][j];
        	b[2][j]=aP[2][j];
        	c[2][j]=-1*alpha_v*aN[2][j];
        	d[2][j]=alpha_v*(S[2][j]+Sdc[2][j]+aW[2][j]*v_star[2][j-2])+(1-alpha_v)*aP[2][j]*v_prev[2][j];
			c_dash[2][j]=c[2][j]/b[2][j];
	    	d_dash[2][j]=d[2][j]/b[2][j];
			for(i=4;i<=2*N-2;i=i+2)
			{
           		a[i][j]=-1*alpha_v*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*alpha_v*aN[i][j];
		    	d[i][j]=alpha_v*(S[i][j]+Sdc[i][j]+aW[i][j]*v_star[i][j-2])+(1-alpha_v)*aP[i][j]*v_prev[i][j];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
		    }

			temp[2*N-2][j]=v_star[2*N-2][j];
			v_star[2*N-2][j]=d_dash[2*N-2][j];

			for(i=2*N-4;i>=2;i=i-2)
			{
				temp[i][j]=v_star[i][j];// to store previous iteration values
				v_star[i][j]=d_dash[i][j]-c_dash[i][j]*v_star[i+2][j];
			}
		}
		else
		{
			a[2][j]=-1*alpha_v*aS[2][j];
        	b[2][j]=aP[2][j];
        	c[2][j]=-1*alpha_v*aN[2][j];
        	d[2][j]=alpha_v*(S[2][j]+Sdc[2][j]+aW[2][j]*v_star[2][j-2]+aE[2][j]*v_star[2][j+2])+(1-alpha_v)*aP[2][j]*v_prev[2][j];
			c_dash[2][j]=c[2][j]/b[2][j];
	    	d_dash[2][j]=d[2][j]/b[2][j];
			for(i=4;i<=2*N-2;i=i+2)
			{
           		a[i][j]=-1*alpha_v*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*alpha_v*aN[i][j];
		    	d[i][j]=alpha_v*(S[i][j]+Sdc[i][j]+aW[i][j]*v_star[i][j-2]+aE[i][j]*v_star[i][j+2])+(1-alpha_v)*aP[i][j]*v_prev[i][j];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-2][j]=v_star[2*N-2][j];
			v_star[2*N-2][j]=d_dash[2*N-2][j];

			for(i=2*N-4;i>=2;i=i-2)
			{
				temp[i][j]=v_star[i][j];// to store previous iteration values
				v_star[i][j]=d_dash[i][j]-c_dash[i][j]*v_star[i+2][j];
			}
		}
		j++;

        if (j<=2*N-2)
		{
			a[1][j]=-1*alpha_u*aS[1][j];
        	b[1][j]=aP[1][j];
        	c[1][j]=-1*alpha_u*aN[1][j];
        	d[1][j]=1*alpha_u*(S[1][j]+Sdc[1][j]+aW[1][j]*u_star[1][j-2]+aE[1][j]*u_star[1][j+2])+(1-alpha_u)*aP[1][j]*u_prev[1][j];
			c_dash[1][j]=c[1][j]/b[1][j];
	    	d_dash[1][j]=d[1][j]/b[1][j];
			for(i=3;i<=2*N-1;i=i+2)
			{
           		a[i][j]=-1*alpha_u*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*alpha_u*aN[i][j];
		    	d[i][j]=1*alpha_u*(S[i][j]+Sdc[i][j]+aW[i][j]*u_star[i][j-2]+aE[i][j]*u_star[i][j+2])+(1-alpha_u)*aP[i][j]*u_prev[i][j];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-1][j]=u_star[2*N-1][j];
			u_star[2*N-1][j]=d_dash[2*N-1][j];

			for(i=2*N-3;i>=1;i=i-2)
			{
				temp[i][j]=u_star[i][j];// to store previous iteration values
				u_star[i][j]=d_dash[i][j]-c_dash[i][j]*u_star[i+2][j];
			}
		}

		j++;
    }
    sum1=0;
    sum2=0;
    sum3=0;
    sum4=0;
    for(i=2;i<=2*N-2;i=i+2)
    {
    	for(j=1;j<=2*N-1;j=j+2)
    	{
    		sum1=sum1+fabs(v_star[i][j]-temp[i][j]);
    		sum2=sum2+fabs(temp[i][j]);
    	}
    }

for(i=1;i<=2*N-1;i=i+2)
    {
    	for(j=2;j<=2*N-2;j=j+2)
    	{
    		sum3=sum3+fabs(u_star[i][j]-temp[i][j]);
    		sum4=sum4+fabs(temp[i][j]);
    	}
    }
    error1=sum1/sum2;
    error2=sum3/sum4;

}while(!(error1<0.00001&&error2<0.00001));

}

void solvePressure()
{
	int i,j;
	double f;
	double error1,error2,sum1,sum2,sum3,sum4;
	double a[2*N][2*N],b[2*N][2*N],c[2*N][2*N],d[2*N][2*N],c_dash[2*N][2*N],d_dash[2*N][2*N];
do{
	for(j=1;j<=2*N-1;j=j+2)
	{
		if (j==1)
		{
			a[1][j]=-1*aS[1][j];
        	b[1][j]=aP[1][j];
        	c[1][j]=-1*aN[1][j];
        	d[1][j]=S[1][j]+aE[1][j]*P_prime[1][j+2];
			c_dash[1][j]=c[1][j]/b[1][j];
	    	d_dash[1][j]=d[1][j]/b[1][j];
			for(i=3;i<=2*N-1;i=i+2)
			{
            	a[i][j]=-1*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*aN[i][j];
		    	d[i][j]=S[i][j]+aE[i][j]*P_prime[i][j+2];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-1][j]=P_prime[2*N-1][j];
			P_prime[2*N-1][j]=d_dash[2*N-1][j];

			for(i=2*N-3;i>=3;i=i-2)
			{
				temp[i][j]=P_prime[i][j];// to store previous iteration values
				P_prime[i][j]=d_dash[i][j]-c_dash[i][j]*P_prime[i+2][j];
			}
		}
		else if(j==2*N-1)
		{
			a[1][j]=-1*aS[1][j];
        	b[1][j]=aP[1][j];
        	c[1][j]=-1*aN[1][j];
        	d[1][j]=S[1][j]+aW[1][j]*P_prime[1][j-2];
			c_dash[1][j]=c[1][j]/b[1][j];
	    	d_dash[1][j]=d[1][j]/b[1][j];
			for(i=3;i<=2*N-1;i=i+2)
			{
            	a[i][j]=-1*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*aN[i][j];
		    	d[i][j]=S[i][j]+aW[i][j]*P_prime[i][j-2];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-1][j]=P_prime[2*N-1][j];
			P_prime[2*N-1][j]=d_dash[2*N-1][j];

			for(i=2*N-3;i>=3;i=i-2)
			{
				temp[i][j]=P_prime[i][j];// to store previous iteration values
				P_prime[i][j]=d_dash[i][j]-c_dash[i][j]*P_prime[i+2][j];
			}
		}
		else
		{
			a[1][j]=-1*aS[1][j];
        	b[1][j]=aP[1][j];
        	c[1][j]=-1*aN[1][j];
        	d[1][j]=S[1][j]+aE[1][j]*P_prime[1][j+2]+aW[1][j]*P_prime[1][j-2];
			c_dash[1][j]=c[1][j]/b[1][j];
	    	d_dash[1][j]=d[1][j]/b[1][j];
			for(i=3;i<=2*N-1;i=i+2)
			{
            	a[i][j]=-1*aS[i][j];
            	b[i][j]=aP[i][j];
            	c[i][j]=-1*aN[i][j];
		    	d[i][j]=S[i][j]+aE[i][j]*P_prime[i][j+2]+aW[i][j]*P_prime[i][j-2];
            	f=1.0/(b[i][j]-(a[i][j]*c_dash[i-2][j]));
            	c_dash[i][j]=c[i][j]*f;
            	d_dash[i][j]=(d[i][j]-(a[i][j]*d_dash[i-2][j]))*f;
			}

			temp[2*N-1][j]=P_prime[2*N-1][j];
			P_prime[2*N-1][j]=d_dash[2*N-1][j];

			for(i=2*N-3;i>=3;i=i-2)
			{
				temp[i][j]=P_prime[i][j];// to store previous iteration values
				P_prime[i][j]=d_dash[i][j]-c_dash[i][j]*P_prime[i+2][j];
			}
		}
	}
	sum1=0;sum2=0;
	for(i=1;i<=2*N-1;i=i+2)
    {
    	for(j=1;j<=2*N-1;j=j+2)
    	{
    		sum1=sum1+fabs(P_prime[i][j]-temp[i][j]);
    		sum2=sum2+fabs(temp[i][j]);
    	}
    }
    error1=sum1/sum2;


   }while(error1>0.0001);
}

string IntToString (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}

int main(){

	int i,j;
	double sum1,sum2,sum3,sum4;
	double allowed_error=0.0005;
	initialize();
	int count=0;
	bool chuck=false;
	int z;
	int steady=0;
    for(z=0;z<=500;z++)
    {
            string num1=IntToString(z);
            string end=".txt";
            num1.append(end);
            cout << num1 << endl;

            count=0;

            vel=0.04*cos(omega*z*dt);
            initialize();
            cout << num1 << " velocity = " << vel << " " << u_star[2*N][1] << "\n";
            do
            {
                update_coefficients();
                solveVelocities();
                update_coefficientsP();
                solvePressure();
                correct_pressure_and_velocities();
                count+=1;

            } //while (error1>allowed_error || error2>allowed_error);
            while(error1>0.0000005 && count<150);

            for(i=2;i<=2*N-2;i=i+2)
                {
                    for(j=1;j<=2*N-1;j=j+2)
                    {
                        v_prev_time[i][j]=v_star[i][j];
                    }
                }

            for(i=1;i<=2*N-1;i=i+2)
                {
                    for(j=2;j<=2*N-2;j=j+2)
                    {
                        u_prev_time[i][j]=u_star[i][j];

                    }
                }

                fstream myfiley;
                //myfilex.open("uvalues.txt");
                myfiley.open(num1.c_str());
                myfiley << vel << endl;
                myfiley << "new row" << endl;
                for(j=0;j<=2*N-2;j=j+2){
                    myfiley << "0.000,0.000" << endl;
                }
                for(i=2;i<=2*N-2;i=i+2){
                    myfiley << "new row" << endl;
                    myfiley << "0.000,0.000" << endl;
                    for(j=2;j<=2*N-2;j=j+2){
                        myfiley << (u_star[i-1][j]+u_star[i+1][j])/2 << ",";
                        myfiley << (v_star[i][j+1]+v_star[i][j-1])/2 << endl;
                    }
                    //cout << "copying... \n";
                }

                myfiley.close();
    }
}

