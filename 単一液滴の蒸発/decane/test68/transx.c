#include <stdio.h>
#include <math.h>
#include "path.h"
#include INC_CC_DEF
#include INC_3_DEF

extern double sigma[spmax],ebyk[spmax],M[spmax],nasaCp[spmax][12];


double Cpi_cal(int sp, double T)
{
	if(T<1000.0)
	{
		return RR*(nasaCp[sp][0]+(nasaCp[sp][1]+(nasaCp[sp][2]+
							(nasaCp[sp][3]+nasaCp[sp][4]*T)*T)*T)*T);
	}
	else
	{
		return RR*(nasaCp[sp][6]+(nasaCp[sp][7]+(nasaCp[sp][8]+
							(nasaCp[sp][9]+nasaCp[sp][10]*T)*T)*T)*T);
	}
}

double Hgi_cal(int sp, double T)
{
	if(T<1000.0)
	{
		return RR*((nasaCp[sp][0]+(nasaCp[sp][1]*0.5+(nasaCp[sp][2]/3.0+
							(nasaCp[sp][3]*0.25+nasaCp[sp][4]*0.2*T)*T)*T)*T)*T+nasaCp[sp][5]);
	}
	else
	{
		return RR*((nasaCp[sp][6]+(nasaCp[sp][7]*0.5+(nasaCp[sp][8]/3.0+
							(nasaCp[sp][9]*0.25+nasaCp[sp][10]*0.2*T)*T)*T)*T)*T+nasaCp[sp][11]);
	}
}

double viscogical(int sp, double T)
{
	double Lsigma,Tas;

	Tas = T/ebyk[sp];
	Lsigma = 1.16145*pow(Tas,-0.14874)+0.52487*exp(-0.7732*Tas)+2.16178*exp(-2.43787*Tas);

	return 26.69*sqrt(M[sp]*1000.0*T)/(sigma[sp]*sigma[sp]*Lsigma);
}

double phical(int sp1, int sp2, double *mui)
{

	return pow(1.0+sqrt(mui[sp1]/mui[sp2])*pow(M[sp2]/M[sp1],0.25),2)/sqrt(8.0*(1.0+M[sp1]/M[sp2]));
}

double viscogcal(double T, double *X)
{
	int sp,j;
	double mui[spmax],mum,BB;

	/* Method of Wilke PGL9-5.13 */

	for(sp=0;sp<spmax;sp++)
	{
		mui[sp] = viscogical(sp,T);
	}

	mum = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		BB = 0.0;
		for(j=0;j<spmax;j++)
		{
			BB += X[j]*phical(sp,j,mui);
		}
		mum += X[sp]*mui[sp]/BB;
	}

	return mum*pow(10.0,-7.0);
}

double lamdagcal(double T, double *X)
{
	/* Wassiljewa Equation PGL10-6.1 */
	int sp,j;
	double lamdai[spmax];
	double mui[spmax],lamdail[spmax],Cvi;
	double lamdam,BB;

	for(sp=0;sp<spmax;sp++)
	{
		mui[sp] = viscogical(sp,T)*pow(10.0,-7.0);
		Cvi = Cpi_cal(sp,T)-RR;
		lamdai[sp] = mui[sp]/M[sp]*(1.32*Cvi+1.77*RR);
	}	
	
	lamdam = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		BB = 0.0;
		for(j=0;j<spmax;j++)
		{
			if(sp==j)
			{
				BB+= X[j]*phical(sp,j,mui);
			}
			else
			{
				BB += 1.065*X[j]*phical(sp,j,mui);
			}
		}
		lamdam += X[sp]*lamdai[sp]/BB;
	}

	return lamdam;
}

double Dij_cal(int i, int j, double T, double P)
{
	double Lsigma,sigmaij,eijbyk,Tab,Mij;

	Mij = 2.0*M[i]*M[j]/(M[i]+M[j]);
	sigmaij = 0.5*(sigma[i]+sigma[j]);
	eijbyk = sqrt(ebyk[i]*ebyk[j]);
	Tab = T/eijbyk;
	Lsigma = 1.06036*pow(Tab,-0.15610)+0.193*exp(-0.47635*Tab)+1.03587*exp(-1.52996*Tab)+1.76474*exp(-3.89411*Tab);

	return 266.0*pow(T,1.5)/(P*sqrt(Mij*1000)*sigmaij*sigmaij*Lsigma);
}

void Diffusion_cal(double *Dsp, double T, double P, double *X)
{
	int sp,j;
	double Dij;

	for(sp=0;sp<spmax;sp++)
	{
		Dsp[sp] = 0.0;
		for(j=0;j<spmax;j++)
		{
			if(sp!=j)
			{
				Dij = Dij_cal(sp,j,T,P);
				Dsp[sp] += X[j]/Dij;
			}
		}
		Dsp[sp] = 1.0/Dsp[sp]*pow(10.0,-4.0);
	}
}

double dHv_cal(double T)
{
	double Tr = T/Tcf;

	return Lf*pow((1.0-Tr)/(1.0-Tb0/Tcf),0.375);
}

/* Liquid phase */
double lamdalcal(double T)
{
	double Astar,alpha,beta,gamma;
	double Tr=T/Tcf;

	/* Latini method */
	Astar = 0.0035;
	alpha = 1.2;
	beta = 0.5;
	gamma = 0.167;

	Astar *= pow(Tb0,alpha)/(pow(M[nFUEL]*1000.0,beta)*pow(Tcf,gamma));

	return	Astar*pow(1.0-Tr,0.38)/pow(Tr,1.0/6.0);
}

double Cplcal(double T)
{
	double Tr=T/Tcf,Cpl;

	Cpl = 1.586+0.49/(1.0-Tr)+ACF*(4.2775+6.3*pow(1.0-Tr,1.0/3.0)/Tr+0.4355/(1.0-Tr));
	Cpl*= RR;
	Cpl+= Cpi_cal(nFUEL,T);
	Cpl/= Mf;

	return Cpl;

	//return 2249.8;
}