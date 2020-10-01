#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "path.h"
#include INC_CC_DEF
#include INC_ALM_H
#include INC_3_DEF
#include INC_KIN
#include INC_TRAN

void output(int k);
void output2(int k);
void output3(int k);
void output4(int k);

/*Å@basic unit is [m], [kg], [K], [mol] */

#define Kmax 5000000
#define Mmax 1000000
#define erefZ 1.0e-8
#define erefT 1.0e-3
#define dtdef 4.0e-6
#define dtdefini  1.0e-6
#define erefCGT 1.0e-6
#define MCGmax 500
#define kchange 200     
#define dtchange 200
#define tmax 0.05

static int Nl=10,Ng=20,NQ=26,NX=100,NZ=200,OUP=0,ELM,NODE;
static double Rs0=0.5e-3,Rsmin=0.005e-3,timestop;
static double Lmax=20.0e-3,Xmax=10.0e-3,rmax=1.6e-3,ld;
static double Tl0=300.0,Pa,Ta,dt=8.0E-6;
static double Rs,Rdot,Rsold,*mdoti,mdot;
static double hrl,byhrl,hrg,byhrg,hs,byhs,hx,byhx,hz,byhz;
static double tref=0.0;
static double tfreq = 5e-4;
static double **Es,**Fs,**Gs,Hs,**Is,**IIs,Js,Ks,Ls,**Ms,Ns,**Os,Us,Vs,**Ws,**cPflgs,**cPflgr,**sPflgx,**sPflgz;
static int **cPflg,**cPflgj,**cPflgi,**sPflg,**sPflgj,**sPflgi;
static char filename[50],filenamea[50],filenameb[50],filenamec[50];

static double **Tgc,**rogc,***Yc,***Xc,**Pc,**dPc,**dPDIFFc,**ux,**uz;
static double **qTgc,**qrogc,***qYc,**qux,**quz;
static double *z,*zm,*x,*xm;
static double **Mmc,**lgc,**Cpgc,**muc,***Dgc,***omegac,***Cpgci,**CgTc;
static double **hgc,***hgci,***mcci;
static double **Gnp,**Gnm;

static int **cflg,**czflg,**cxflg;
static int **cflgA,cflgAmax;

static double **Tgs,**rogs,***Ys,***Xs,**Ps,**dPs,**dPDIFFs,**ur,**uq;
static double **qTgs,**qrogs,***qYs,**qur,**quq;
static double *eta,*etam,*rgm,*rg,*rgm2,*rg2;
static double **Mms,**lgs,**Cpgs,**mus,***Dgs,***omegas,***Cpgsi,**CgTs;
static double **hgs,***hgsi,***mcsi;
static double **Dnp,**Dnm;

static double **Tls,**eTls,**qTls;
static double *xi,*xim,*rlm,*rl;
static double *ssm,*ss,*scm,*sc;
static double **ll,**Cpl;

extern double Y0[spmax],M[spmax];

void alm()
{
	#define vec(Vec,Nv)			(Vec)=(dvec((Nv)))
	#define mat(Mat,Nm,Nv)		(Mat)=(dmat((Nm),(Nv)))
	#define ten(Ten,Nt,Nm,Nv)	(Ten)=(dten((Nt),(Nm),(Nv)))

	mat(Tgc,NZ+1,NX+1);mat(rogc,NZ+1,NX+1);ten(Yc,NZ+1,NX+1,spmax);mat(Pc,NZ+1,NX+1);mat(dPc,NZ+1,NX+1);mat(dPDIFFc,NZ+1,NX+1);
	mat(ux,NZ+1,NX+1);mat(uz,NZ+1,NX+1);
	mat(qTgc,NZ+1,NX+1);mat(qrogc,NZ+1,NX+1);ten(qYc,NZ+1,NX+1,spmax);mat(qux,NZ+1,NX+1);mat(quz,NZ+1,NX+1);ten(Xc,NZ+1,NX+1,spmax);
	mat(Gnp,NZ+1,NX+1);mat(Gnm,NZ+1,NX+1);
	vec(z,NZ+1);vec(zm,NZ);vec(x,NX+1);vec(xm,NX);
	mat(Mmc,NZ+1,NX+1);mat(lgc,NZ+1,NX+1);mat(Cpgc,NZ+1,NX+1);mat(muc,NZ+1,NX+1);ten(Dgc,NZ+1,NX+1,spmax);ten(omegac,NZ+1,NX+1,spmax);
	ten(Cpgci,NZ+1,NX+1,spmax);mat(hgc,NZ+1,NX+1);ten(hgci,NZ+1,NX+1,spmax);ten(mcci,NZ+1,NX+1,spmax);mat(CgTc,NZ+1,NX+1);
	mat(cPflgr,NZ+1,NX+1);mat(cPflgs,NZ+1,NX+1);

	cflg = imat(NZ+1,NX+1);
	cxflg = imat(NZ+1,NX+1);
	czflg = imat(NZ+1,NX+1);
	cPflg = imat(NZ+1,NX+1);
	sPflg = imat(NQ,Ng+1);
	cPflgi = imat(NZ+1,NX+1);
	cPflgj = imat(NZ+1,NX+1);
	sPflgi = imat(NQ,Ng+1);
	sPflgj = imat(NQ,Ng+1);
    
	mat(Es,NZ+NQ+1,NX+1);mat(Fs,NZ+NQ+1,NX+1);mat(Gs,NZ+NQ+1,NX+1);mat(Is,NZ+NQ+1,NX+1);mat(IIs,NZ+NQ+1,NX+1);mat(Ms,NZ+NQ+1,NX+1);mat(Os,NZ+NQ+1,NX+1);mat(Ws,NZ+NQ+1,NX+1);mat(sPflgx,NQ,Ng+1);mat(sPflgz,NQ,Ng+1);
	mat(Tgs,NQ,Ng+1);mat(rogs,NQ,Ng+1);ten(Ys,NQ,Ng+1,spmax);mat(Ps,NQ,Ng+1);mat(dPs,NQ,Ng+1);mat(dPDIFFs,NQ,Ng+1);
	mat(uq,NQ+1,Ng+1);mat(ur,NQ,Ng+1);
	mat(qTgs,NQ,Ng+1);mat(qrogs,NQ,Ng+1);ten(qYs,NQ,Ng+1,spmax);mat(quq,NQ+1,Ng+1);mat(qur,NQ,Ng+1);ten(Xs,NQ,Ng+1,spmax);
	mat(Dnp,NQ,Ng+1);mat(Dnm,NQ,Ng+1);
	vec(eta,Ng+1);vec(etam,Ng);vec(rgm,Ng);vec(rg,Ng+1);vec(rgm2,Ng);vec(rg2,Ng+1);
	mat(Mms,NQ,Ng+1);mat(lgs,NQ,Ng+1);mat(Cpgs,NQ,Ng+1);mat(mus,NQ,Ng+1);ten(Dgs,NQ,Ng+1,spmax);ten(omegas,NQ,Ng+1,spmax);
	ten(Cpgsi,NQ,Ng+1,spmax);mat(hgs,NQ,Ng+1);ten(hgsi,NQ,Ng+1,spmax);ten(mcsi,NQ,Ng+1,spmax);mat(CgTs,NQ,Ng+1);

	mat(Tls,NQ,Nl+1);mat(qTls,NQ,Nl);
	vec(xi,Nl+1);vec(xim,Nl);vec(rl,Nl+1);vec(rlm,Nl+1);
	vec(ss,NQ+1);vec(ssm,NQ);vec(sc,NQ+1);vec(scm,NQ);
	mat(ll,NQ,Nl+1);mat(Cpl,NQ,Nl+1);

	vec(mdoti,NQ);

	#undef vec
	#undef mat
	#undef ten
}

void decimalconvert(char *array)
{
	int i=0;

	do{
		if (array[i]=='.') array[i]='_';
	}while(array[i++]!=0);
}

void filedef()
{
	char charPa[20],charTa[20],charld[20],chard[20];

	sprintf(charPa,"%.1lf",Pa);
	sprintf(charTa,"%.1lf",Ta);
	sprintf(charld,"%.1lf",ld);
	sprintf(chard,"%.1lf",2.0*Rs0);

	decimalconvert(charPa);decimalconvert(charTa);decimalconvert(charld);decimalconvert(chard);

	strcpy(filename,charPa);
	strcat(filename,"p");
	strcat(filename,charTa);
	strcat(filename,"t");
	strcat(filename,charld);
	strcat(filename,"d");
	strcat(filename,chard);
	strcpy(filenamea,"figure");
	strcpy(filenameb,"vector");
	strcpy(filenamec,"record");

}

double Mmcal(double *Y)
{
	int sp;
	double Mm;

	Mm = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Mm += Y[sp]/M[sp];
	}
	return 1.0/Mm;
}

double Cpg_cal(double T, double *Cpgi, double *Y)
{
	int sp;
	double Cpg;
	
	Cpg = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Cpgi[sp] = Cpi_cal(sp,T)/M[sp];
		Cpg += Cpgi[sp]*Y[sp];
	}

	return Cpg;
}

double Hg_cal(double T, double *hgi, double *Y)
{
	int sp;
	double Hg;

	Hg = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		hgi[sp] = Hgi_cal(sp,T)/M[sp];
		Hg += hgi[sp]*Y[sp];
	}

	return Hg;
}

double CgT_cal(double *omega, double *hgi)
{
	double CgT=0.0;
	int sp;

	for(sp=0;sp<spmax;sp++)
	{
		CgT+= omega[sp]*hgi[sp];
	}

	return CgT;
}

void YtoX(double *X, double *Y, double Mm)
{
	int sp;
	double Am=0.0;

	for(sp=0;sp<spmax-1;sp++)
	{
		X[sp] = Y[sp]*Mm/M[sp];
		Am+= X[sp];
	}
	X[nN2] = 1.0 - Am;
}


void propertygs_cal(int i, int j)
{
	Mms[i][j] = Mmcal(Ys[i][j]);
	YtoX(Xs[i][j],Ys[i][j],Mms[i][j]);
	lgs[i][j] = lamdagcal(Tgs[i][j],Xs[i][j]);
	mus[i][j] = viscogcal(Tgs[i][j],Xs[i][j]);
	Diffusion_cal(Dgs[i][j], Tgs[i][j], Pa, Xs[i][j]);
	Cpgs[i][j] = Cpg_cal(Tgs[i][j],Cpgsi[i][j],Ys[i][j]);
	hgs[i][j] = Hg_cal(Tgs[i][j],hgsi[i][j],Ys[i][j]);
	Compute_Chemical_Kinetics(omegas[i][j],Tgs[i][j],Ys[i][j],mcsi[i][j],rogs[i][j]);
	CgTs[i][j] = CgT_cal(omegas[i][j],hgsi[i][j]);
}

void propertygc_cal(int i, int j)
{
	Mmc[i][j] = Mmcal(Yc[i][j]);
	YtoX(Xc[i][j],Yc[i][j],Mmc[i][j]);
	lgc[i][j] = lamdagcal(Tgc[i][j],Yc[i][j]);
	muc[i][j] = viscogcal(Tgc[i][j],Yc[i][j]);
	Cpgc[i][j] = Cpg_cal(Tgc[i][j],Cpgci[i][j],Yc[i][j]);
	Diffusion_cal(Dgc[i][j],Tgc[i][j],Pa,Yc[i][j]);
	hgc[i][j] = Hg_cal(Tgc[i][j],hgci[i][j],Yc[i][j]);
	Compute_Chemical_Kinetics(omegac[i][j],Tgc[i][j],Yc[i][j],mcci[i][j],rogc[i][j]);
}

void propertyl_cal(int i, int j)
{
	ll[i][j] = lamdalcal(Tls[i][j]);
	Cpl[i][j] = Cplcal(Tls[i][j]);
}

double rogcal(double Mm, double T, double P)
{
	return P*Mm/(RR*T);
}

double Pgcal(double Mm, double T, double rog)
{
	return rog*RR*T/Mm;
}

double N2sum(double *Y)
{
	int sp;
	double AA;

	AA = 1.0;
	for(sp=0;sp<spmax-1;sp++)
	{
		AA -= Y[sp];
	}

	return AA;

}

void rlcal()
{
	int i;

	for(i=0;i<=Nl;i++)
	{
		rl[i] = Rs*(double)i/byhrl;
	}
	for(i=0;i<Nl;i++)
	{
		rlm[i] = Rs*((double)i+0.5)/byhrl;
	}
}

void rgcal()
{
	int i;

	for(i=0;i<=Ng;i++)
	{
		rg2[i] = rg[i];
		rg[i] = Rs*pow(rmax/Rs,(double)i/byhrg);
	}
	for(i=0;i<Ng;i++)
	{
		rgm2[i] = rgm[i];
		rgm[i] = Rs*pow(rmax/Rs,((double)i+0.5)/byhrg);
	}
}

void Search()
{
	int i,j,k,flg,ii,jj,sp;
	double x0,z0,aa,bb;
	double /***A,*/**b;
	double u1,u2;
	double v1,v2,v3,v4;
	double r0,s0;

	//A = dmat(3,3);
	b = dmat(spmax+1,3);
	/* Spherical coordinate system */
	/* scalar */
	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{
			x0 = rgm[i-1]*ssm[j];
			z0 = rgm[i-1]*scm[j] + ld;
			flg = 0;
			k = 0;
			do
			{
				if(z0<zm[k])
				{
					flg = 1;
					ii = k;
				}
				k += 1;
			}while(flg==0&&k<=NZ);
			flg = 0;
			k = 0;
			do
			{
				if(x0<xm[k])
				{
					flg = 1;
					jj = k;
				}
				k+=1;
			}while(flg==0&&k<=NX);



			if(jj==0)
			{
				aa = 0.5*xm[jj];
				bb = 0.5*(zm[ii-1]+zm[ii]);
				if(x0>aa)
				{ if(z0>bb){
					b[0][0] = (Tgc[ii][0]-Tgc[ii-1][0])*byhz;
					b[0][1] = 0.25*(Tgc[ii][1]-Tgc[ii][0])*byhx;
					b[0][2] = Tgc[ii][0] - b[0][0]*zm[ii] - b[0][1]*xm[0];
					b[1][0] = (Pc[ii][0]-Pc[ii-1][0])*byhz;
					b[1][1] = 0.25*(Pc[ii][1]-Pc[ii][0])*byhx;
					b[1][2] = Pc[ii][0] - b[1][0]*zm[ii] - b[1][1]*xm[0];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][0][sp]-Yc[ii-1][0][sp])*byhz;
						b[sp+2][1] = 0.25*(Yc[ii][1][sp]-Yc[ii][0][sp])*byhx;
						b[sp+2][2] = Yc[ii][0][sp] - b[sp+2][0]*zm[ii] - b[sp+2][1]*xm[0];
					}
				}
				else
				{
					b[0][0] = (Tgc[ii][0]-Tgc[ii-1][0])*byhz;
					b[0][1] = 0.25*(Tgc[ii-1][1]-Tgc[ii-1][0])*byhx;
					b[0][2] = Tgc[ii-1][0] - b[0][0]*zm[ii-1] - b[0][1]*xm[0];
					b[1][0] = (Pc[ii][0]-Pc[ii-1][0])*byhz;
					b[1][1] = 0.25*(Pc[ii-1][1]-Pc[ii-1][0])*byhx;
					b[1][2] = Pc[ii-1][0] - b[1][0]*zm[ii-1] - b[1][1]*xm[0];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][0][sp]-Yc[ii-1][0][sp])*byhz;
						b[sp+2][1] = 0.25*(Yc[ii-1][1][sp]-Yc[ii-1][0][sp])*byhx;
						b[sp+2][2] = Yc[ii-1][0][sp] - b[sp+2][0]*zm[ii-1] - b[sp+2][1]*xm[0];
					}
				}}
				else
				{if(z0>bb){
					b[0][0] = 0.125*(9.0*(Tgc[ii][0]-Tgc[ii-1][0])-(Tgc[ii][1]-Tgc[ii-1][1]))*byhz;
					b[0][1] = 0.25*(Tgc[ii][1]-Tgc[ii][0])*byhx;
					b[0][2] = Tgc[ii][0] - b[0][0]*zm[ii] - b[0][1]*xm[0];
					b[1][0] = 0.125*(9.0*(Pc[ii][0]-Pc[ii-1][0])-(Pc[ii][1]-Pc[ii-1][1]))*byhz;
					b[1][1] = 0.25*(Pc[ii][1]-Pc[ii][0])*byhx;
					b[1][2] = Pc[ii][0] - b[1][0]*zm[ii] - b[1][1]*xm[0];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = 0.125*(9.0*(Yc[ii][0][sp]-Yc[ii-1][0][sp])-(Yc[ii][1][sp]-Yc[ii-1][1][sp]))*byhz;
						b[sp+2][1] = 0.25*(Yc[ii][1][sp]-Yc[ii][0][sp])*byhx;
						b[sp+2][2] = Yc[ii][0][sp] - b[sp+2][0]*zm[ii] - b[sp+2][1]*xm[0];
					}
				}				
				else
				{
					b[0][0] = 0.125*(9.0*(Tgc[ii][0]-Tgc[ii-1][0])-(Tgc[ii][1]-Tgc[ii-1][1]))*byhz;
					b[0][1] = 0.25*(Tgc[ii-1][1]-Tgc[ii-1][0])*byhx;
					b[0][2] = Tgc[ii-1][0] - b[0][0]*zm[ii-1] - b[0][1]*xm[0];
					b[1][0] = 0.125*(9.0*(Pc[ii][0]-Pc[ii-1][0])-(Pc[ii][1]-Pc[ii-1][1]))*byhz;
					b[1][1] = 0.25*(Pc[ii-1][1]-Pc[ii-1][0])*byhx;
					b[1][2] = Pc[ii-1][0] - b[1][0]*zm[ii-1] - b[1][1]*xm[0];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = 0.125*(9.0*(Yc[ii][0][sp]-Yc[ii-1][0][sp])-(Yc[ii][1][sp]-Yc[ii-1][1][sp]))*byhz;
						b[sp+2][1] = 0.25*(Yc[ii-1][1][sp]-Yc[ii-1][0][sp])*byhx;
						b[sp+2][2] = Yc[ii-1][0][sp] - b[sp+2][0]*zm[ii-1] - b[sp+2][1]*xm[0];
					}
				}
				}}
			else
			{
				aa = (xm[jj-1]+xm[jj])/2;
				bb = (zm[ii-1]+zm[ii])/2;
				if(x0>aa)
				{if(z0>bb){
					b[0][0] = (Tgc[ii][jj]-Tgc[ii-1][jj])*byhz;
					b[0][1] = (Tgc[ii][jj]-Tgc[ii][jj-1])*byhx;
					b[0][2] = Tgc[ii][jj] - b[0][0]*zm[ii] - b[0][1]*xm[jj];
					b[1][0] = (Pc[ii][jj]-Pc[ii-1][jj])*byhz;
					b[1][1] = (Pc[ii][jj]-Pc[ii][jj-1])*byhx;
					b[1][2] = Pc[ii][jj] - b[1][0]*zm[ii] - b[1][1]*xm[jj];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][jj][sp]-Yc[ii-1][jj][sp])*byhz;
						b[sp+2][1] = (Yc[ii][jj][sp]-Yc[ii][jj-1][sp])*byhx;
						b[sp+2][2] = Yc[ii][jj][sp] - b[sp+2][0]*zm[ii] - b[sp+2][1]*xm[jj];
					}
				}
				else
				{
					b[0][0] = (Tgc[ii][jj]-Tgc[ii-1][jj])*byhz;
					b[0][1] = (Tgc[ii-1][jj]-Tgc[ii-1][jj-1])*byhx;
					b[0][2] = Tgc[ii-1][jj] - b[0][0]*zm[ii-1] - b[0][1]*xm[jj];
					b[1][0] = (Pc[ii][jj]-Pc[ii-1][jj])*byhz;
					b[1][1] = (Pc[ii-1][jj]-Pc[ii-1][jj-1])*byhx;
					b[1][2] = Pc[ii-1][jj] - b[1][0]*zm[ii-1] - b[1][1]*xm[jj];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][jj][sp]-Yc[ii-1][jj][sp])*byhz;
						b[sp+2][1] = (Yc[ii-1][jj][sp]-Yc[ii-1][jj-1][sp])*byhx;
						b[sp+2][2] = Yc[ii-1][jj][sp] - b[sp+2][0]*zm[ii-1] - b[sp+2][1]*xm[jj];
					}
				}}
				else
				{if(z0>bb){
					b[0][0] = (Tgc[ii][jj-1]-Tgc[ii-1][jj-1])*byhz;
					b[0][1] = (Tgc[ii][jj]-Tgc[ii][jj-1])*byhx;
					b[0][2] = Tgc[ii][jj-1] - b[0][0]*zm[ii] - b[0][1]*xm[jj-1];
					b[1][0] = (Pc[ii][jj-1]-Pc[ii-1][jj-1])*byhz;
					b[1][1] = (Pc[ii][jj]-Pc[ii][jj-1])*byhx;
					b[1][2] = Pc[ii][jj-1] - b[1][0]*zm[ii] - b[1][1]*xm[jj-1];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][jj-1][sp]-Yc[ii-1][jj-1][sp])*byhz;
						b[sp+2][1] = (Yc[ii][jj][sp]-Yc[ii][jj-1][sp])*byhx;
						b[sp+2][2] = Yc[ii][jj-1][sp] - b[sp+2][0]*zm[ii] - b[sp+2][1]*xm[jj-1];
					}
				}
				else
				{
					b[0][0] = (Tgc[ii][jj-1]-Tgc[ii-1][jj-1])*byhz;
					b[0][1] = (Tgc[ii-1][jj]-Tgc[ii-1][jj-1])*byhx;
					b[0][2] = Tgc[ii-1][jj-1] - b[0][0]*zm[ii-1] - b[0][1]*xm[jj-1];
					b[1][0] = (Pc[ii][jj-1]-Pc[ii-1][jj-1])*byhz;
					b[1][1] = (Pc[ii-1][jj]-Pc[ii-1][jj-1])*byhx;
					b[1][2] = Pc[ii-1][jj-1] - b[1][0]*zm[ii-1] - b[1][1]*xm[jj-1];
					for(sp=0;sp<spmax-1;sp++)
					{
						b[sp+2][0] = (Yc[ii][jj-1][sp]-Yc[ii-1][jj-1][sp])*byhz;
						b[sp+2][1] = (Yc[ii-1][jj][sp]-Yc[ii-1][jj-1][sp])*byhx;
						b[sp+2][2] = Yc[ii-1][jj-1][sp] - b[sp+2][0]*zm[ii-1] - b[sp+2][1]*xm[jj-1];
					}}}
			}


			Tgs[j][i] = b[0][0]*z0+b[0][1]*x0+b[0][2];
			Ps[j][i] = b[1][0]*z0+b[1][1]*x0+b[1][2];
			for(sp=0;sp<spmax-1;sp++)
			{
				Ys[j][i][sp] = b[sp+2][0]*z0+b[sp+2][1]*x0+b[sp+2][2];
			}
			Ys[j][i][nN2] = N2sum(Ys[j][i]);
			propertygs_cal(j,i);
			rogs[j][i] = rogcal(Mms[j][i],Tgs[j][i],Pa);
		}
	}
	/* r Velocity */
	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{
			x0 = rg[i]*ssm[j];
			z0 = rg[i]*scm[j] + ld;
			flg = 0;
			k = 0;
			do
			{
				if(z0<zm[k])
				{
					flg = 1;
					ii = k;
				}
				k += 1;
			}while(flg==0&&k<=NZ);
			flg = 0;
			k = 0;
			do
			{
				if(x0<xm[k])
				{
					flg = 1;
					jj = k;
				}
				k+=1;
			}while(flg==0&&k<=NX);
			if(jj==0)
			{   aa = 0.5*xm[jj];
				bb = 0.5*(zm[ii-1]+zm[ii]);
            if(x0>aa)
				{
				u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
				u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
				v1 = u1*scm[j]+u2*ssm[j];
				u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
				u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
				v2 = u1*scm[j]+u2*ssm[j];			
				if(z0>bb){u1 =0.0625*(9.0*uz[ii+1][0]-uz[ii+1][1]+9.0*uz[ii][0]-uz[ii][1]);
					u2 =0.0;
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v1-v3)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii]- b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.0625*(9.0*uz[ii][0]-uz[ii][1]+9.0*uz[ii-1][0]-uz[ii-1][1]);
					u2 = 0.0;
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v2-v3)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1]- b[0][1]*xm[jj];
				}
			}
			else
			{
				u1 = 0.0625*(9.0*uz[ii+1][0]-uz[ii+1][1]+9.0*uz[ii][0]-uz[ii][1]);
				u2 = 0.0;
				v1 = u1*scm[j]+u2*ssm[j];
				u1 = 0.0625*(9.0*uz[ii-1][0]-uz[ii-1][1]+9.0*uz[ii][0]-uz[ii][1]);
				u2 = 0.0;
				v2 = u1*scm[j]+u2*ssm[j];			
					if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v3-v1)*byhx;
					b[0][2] = v3 - b[0][0]*zm[ii]- b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.5*(uz[ii-1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v3-v2)*byhx;
					b[0][2] = v3- b[0][0]*zm[ii-1]- b[0][1]*xm[jj];
				}
			}}
			else
			{	aa = (xm[jj-1]+xm[jj])/2;
				bb = (zm[ii-1]+zm[ii])/2;
				 if(x0>aa)
				{
					
				u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
				u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
				v1 = u1*scm[j]+u2*ssm[j];
				u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
				u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
				v2 = u1*scm[j]+u2*ssm[j];
				
				if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj-1]+uz[ii][jj-1]);
					u2 = 0.5*(ux[ii][jj]+ux[ii][jj-1]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v1-v3)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii] - b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.5*(uz[ii][jj-1]+uz[ii-1][jj-1]);
					u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v2-v3)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1] - b[0][1]*xm[jj];
				}
			}
				 else
				{

				u1 = 0.5*(uz[ii+1][jj-1]+uz[ii][jj-1]);
				u2 = 0.5*(ux[ii][jj]+ux[ii][jj-1]);
				v1 = u1*scm[j]+u2*ssm[j];
				u1 = 0.5*(uz[ii][jj-1]+uz[ii-1][jj-1]);
				u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
				v2 = u1*scm[j]+u2*ssm[j];
				
				if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v3-v1)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii] - b[0][1]*xm[jj-1];
				}
				else
				{
					u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
					u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
					v3 = u1*scm[j]+u2*ssm[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v3-v2)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1] - b[0][1]*xm[jj-1];
				}
				 }}
			ur[j][i] = b[0][0]*z0+b[0][1]*x0+b[0][2];
		}
	}
	/* shita Velocity */
	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=1;j<NQ;j++)
		{
			x0 = rgm[i-1]*ss[j];
			z0 = rgm[i-1]*sc[j] + ld;
			flg = 0;
			k = 0;
			do
			{
				if(z0<zm[k])
				{
					flg = 1;
					ii = k;
				}
				k += 1;
			}while(flg==0&&k<=NZ);
			flg = 0;
			k = 0;
			do
			{
				if(x0<xm[k])
				{
					flg = 1;
					jj = k;
				}
				k+=1;
			}while(flg==0&&k<=NX);
			if(jj==0)
			{   aa = 0.5*xm[jj];
				bb = 0.5*(zm[ii-1]+zm[ii]);
            if(x0>aa)
				{
				u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
				u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
				v1 = u2*sc[j]-u1*ss[j];
				u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
				u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
				v2 = u2*sc[j]-u1*ss[j];	
				if(z0>bb){u1 =0.0625*(9.0*uz[ii+1][0]-uz[ii+1][1]+9.0*uz[ii][0]-uz[ii][1]);
					u2 =0.0;
					v3 =u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v1-v3)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii]- b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.0625*(9.0*uz[ii][0]-uz[ii][1]+9.0*uz[ii-1][0]-uz[ii-1][1]);
					u2 = 0.0;
					v3 =u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v2-v3)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1]- b[0][1]*xm[jj];
				}
			}
			else
			{
				u1 = 0.0625*(9.0*uz[ii+1][0]-uz[ii+1][1]+9.0*uz[ii][0]-uz[ii][1]);
				u2 = 0.0;
				v1 = u2*sc[j]-u1*ss[j];
				u1 = 0.0625*(9.0*uz[ii-1][0]-uz[ii-1][1]+9.0*uz[ii][0]-uz[ii][1]);
				u2 = 0.0;
				v2 = u2*sc[j]-u1*ss[j];		
					if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
					v3 = u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v3-v1)*byhx;
					b[0][2] = v3 - b[0][0]*zm[ii]- b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.5*(uz[ii-1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
					v3 = u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = 2.0*(v3-v2)*byhx;
					b[0][2] = v3- b[0][0]*zm[ii-1]- b[0][1]*xm[jj];
				}
			}}
			else
			{	aa = (xm[jj-1]+xm[jj])/2;
				bb = (zm[ii-1]+zm[ii])/2;
				 if(x0>aa)
				{
					
				u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
				u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
				v1 =u2*sc[j]-u1*ss[j];
				u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
				u2 = 0.5*(ux[ii-1][jj+1]+ux[ii-1][jj]);
				v2 =u2*sc[j]-u1*ss[j];
				
				if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj-1]+uz[ii][jj-1]);
					u2 = 0.5*(ux[ii][jj]+ux[ii][jj-1]);
					v3 =u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v1-v3)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii] - b[0][1]*xm[jj];
				}
				else
				{
					u1 = 0.5*(uz[ii][jj-1]+uz[ii-1][jj-1]);
					u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
					v3 = u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v2-v3)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1] - b[0][1]*xm[jj];
				}
			}
				 else
				{

				u1 = 0.5*(uz[ii+1][jj-1]+uz[ii][jj-1]);
				u2 = 0.5*(ux[ii][jj]+ux[ii][jj-1]);
				v1 = u2*sc[j]-u1*ss[j];
				u1 = 0.5*(uz[ii][jj-1]+uz[ii-1][jj-1]);
				u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
				v2 =u2*sc[j]-u1*ss[j];
				
				if(z0>bb){
					u1 = 0.5*(uz[ii+1][jj]+uz[ii][jj]);
					u2 = 0.5*(ux[ii][jj+1]+ux[ii][jj]);
					v3 = u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v3-v1)*byhx;
					b[0][2] = v1 - b[0][0]*zm[ii] - b[0][1]*xm[jj-1];
				}
				else
				{
					u1 = 0.5*(uz[ii][jj]+uz[ii-1][jj]);
					u2 = 0.5*(ux[ii-1][jj]+ux[ii-1][jj-1]);
					v3 = u2*sc[j]-u1*ss[j];
					b[0][0] = (v1-v2)*byhz;
					b[0][1] = (v3-v2)*byhx;
					b[0][2] = v2 - b[0][0]*zm[ii-1] - b[0][1]*xm[jj-1];
				}
				 }}
			
			
		
			uq[j][i] = (b[0][0]*z0+b[0][1]*x0+b[0][2]);

		}
	}



	/* cylindrical coordinate system */
	/* schalar */
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			cflg[i][j] = 0;
			cxflg[i][j] = 0;
			czflg[i][j] = 0;
		}
	}

	cflgAmax = 0;
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			r0 = sqrt(xm[j]*xm[j]+(zm[i]-ld)*(zm[i]-ld));
			if(r0<rgm[Ng-3])
			{
				cflg[i][j] = 1;
				if(r0>rgm[0])
				{
					cflg[i][j] = 2;
					cflgAmax+= 1;
					s0 = acos((zm[i]-ld)/r0);
					flg = 0;
					k = 0;
					do
					{
						if(rgm[k]>r0)
						{
							flg = 1;
							jj = k;
						}
						k+=1;
					}while(flg==0);
					flg = 0;
					k = 0;
					do
					{
						if(((double)k+0.5)*hs>s0)
						{
							flg = 1;
							ii = k;
						}
						k+=1;
					}while(flg==0&&k<NQ);

						if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.75)*hs;
					if(aa>r0)
					{   if(bb>s0){b[0][0] = 0.25*(Tgs[ii][jj]-Tgs[ii-1][jj])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii][jj]-Ps[ii-1][jj])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj][sp]-Ys[ii-1][jj][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}
					    else{b[0][0] = 0.25*(Tgs[ii][jj]-Tgs[ii-1][jj])*byhs;
							b[0][1] =  0.125*(9.0*(Tgs[ii][jj+1]-Tgs[ii][jj])-(Tgs[ii-1][jj+1]-Tgs[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii][jj]-Ps[ii-1][jj])*byhs;
							b[1][1] =  0.125*(9.0*(Ps[ii][jj+1]-Ps[ii][jj])-(Ps[ii-1][jj+1]-Ps[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj][sp]-Ys[ii-1][jj][sp])*byhs;
								b[sp+2][1] = 0.125*(9.0*(Ys[ii][jj+1][sp]-Ys[ii][jj][sp])-(Ys[ii-1][jj+1][sp]-Ys[ii-1][jj][sp]))/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}}
					else
					{  if(bb>s0){b[0][0] = 0.25*(Tgs[ii][jj+1]-Tgs[ii-1][jj+1])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj+1] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii][jj+1]-Ps[ii-1][jj+1])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj+1][sp]-Ys[ii-1][jj+1][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj+1][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj];
							}}
					     else{b[0][0] = 0.25*(Tgs[ii][jj+1]-Tgs[ii-1][jj+1])*byhs;
							b[0][1] =  0.125*(9.0*(Tgs[ii][jj+1]-Tgs[ii][jj])-(Tgs[ii-1][jj+1]-Tgs[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj+1] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii][jj+1]-Ps[ii-1][jj+1])*byhs;
							b[1][1] =  0.125*(9.0*(Ps[ii][jj+1]-Ps[ii][jj])-(Ps[ii-1][jj+1]-Ps[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj+1][sp]-Ys[ii-1][jj+1][sp])*byhs;
								b[sp+2][1] = 0.125*(9.0*(Ys[ii][jj+1][sp]-Ys[ii][jj][sp])-(Ys[ii-1][jj+1][sp]-Ys[ii-1][jj][sp]))/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj+1][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj];
							}}}
				}


				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.25*hs;
					if(aa>r0)
					{   if(bb>s0){b[0][0] = 0.25*(Tgs[ii+1][jj]-Tgs[ii][jj])*byhs;
							b[0][1] = 0.125*(9.0*(Tgs[ii][jj+1]-Tgs[ii][jj])-(Tgs[ii+1][jj+1]-Tgs[ii+1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii][jj]-Ps[ii+1][jj])*byhs;
							b[1][1] =0.125*(9.0*(Ps[ii][jj+1]-Ps[ii][jj])-(Ps[ii+1][jj+1]-Ps[ii+1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii+1][jj][sp]-Ys[ii][jj][sp])*byhs;
								b[sp+2][1] = 0.125*(9.0*(Ys[ii][jj+1][sp]-Ys[ii][jj][sp])-(Ys[ii+1][jj+1][sp]-Ys[ii+1][jj][sp]))/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}
					    else{b[0][0] = 0.25*(Tgs[ii+1][jj]-Tgs[ii][jj])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii+1][jj]-Ps[ii][jj])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii+1][jj][sp]-Ys[ii][jj][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}}
					else
					{  if(bb>s0){b[0][0] = 0.25*(Tgs[ii+1][jj+1]-Tgs[ii][jj+1])*byhs;
							b[0][1] =  0.125*(9.0*(Tgs[ii][jj+1]-Tgs[ii][jj])-(Tgs[ii+1][jj+1]-Tgs[ii+1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj+1] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii+1][jj+1]-Ps[ii][jj+1])*byhs;
							b[1][1] = 0.125*(9.0*(Ps[ii][jj+1]-Ps[ii][jj])-(Ps[ii+1][jj+1]-Ps[ii+1][jj]))/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii+1][jj+1][sp]-Ys[ii][jj+1][sp])*byhs;
								b[sp+2][1] = 0.125*(9.0*(Ys[ii][jj+1][sp]-Ys[ii][jj][sp])-(Ys[ii+1][jj+1][sp]-Ys[ii+1][jj][sp]))/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj+1][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj];
							}}
					     else{b[0][0] = 0.25*(Tgs[ii+1][jj+1]-Tgs[ii][jj+1])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj+1] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii+1][jj+1]-Ps[ii][jj+1])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] =  0.25*(Ys[ii+1][jj+1][sp]-Ys[ii][jj+1][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj+1][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj];
							}}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa>r0)
					{   if(bb>s0){b[0][0] = 0.25*(Tgs[ii][jj]-Tgs[ii-1][jj])*byhs;
							b[0][1] = (Tgs[ii-1][jj+1]-Tgs[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii-1][jj] - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii][jj]-Ps[ii-1][jj])*byhs;
							b[1][1] = (Ps[ii-1][jj+1]-Ps[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii-1][jj] - b[1][0]*((double)ii-0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj][sp]-Ys[ii-1][jj][sp])*byhs;
								b[sp+2][1] = (Ys[ii-1][jj+1][sp]-Ys[ii-1][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii-1][jj][sp] - b[sp+2][0]*((double)ii-0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}
					    else{b[0][0] = 0.25*(Tgs[ii][jj]-Tgs[ii-1][jj])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];
							b[1][0] = 0.25*(Ps[ii][jj]-Ps[ii-1][jj])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj-1];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj][sp]-Ys[ii-1][jj][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj-1];
							}}}
					else
					{  if(bb>s0){b[0][0] = 0.25*(Tgs[ii][jj+1]-Tgs[ii-1][jj+1])*byhs;
							b[0][1] = (Tgs[ii-1][jj+1]-Tgs[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii-1][jj+1] - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii][jj+1]-Ps[ii-1][jj+1])*byhs;
							b[1][1] =  (Ps[ii-1][jj+1]-Ps[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii-1][jj+1] - b[1][0]*((double)ii-0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj+1][sp]-Ys[ii-1][jj+1][sp])*byhs;
								b[sp+2][1] = (Ys[ii-1][jj+1][sp]-Ys[ii-1][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii-1][jj+1][sp] - b[sp+2][0]*((double)ii-0.5)*hs - b[sp+2][1]*rgm[jj];
							}}
					     else{b[0][0] = 0.25*(Tgs[ii][jj+1]-Tgs[ii-1][jj+1])*byhs;
							b[0][1] = (Tgs[ii][jj+1]-Tgs[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[0][2] = Tgs[ii][jj+1] - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];
							b[1][0] = 0.25*(Ps[ii][jj+1]-Ps[ii-1][jj+1])*byhs;
							b[1][1] = (Ps[ii][jj+1]-Ps[ii][jj])/(rgm[jj]-rgm[jj-1]);
							b[1][2] = Ps[ii][jj+1] - b[1][0]*((double)ii+0.5)*hs - b[1][1]*rgm[jj];
							for(sp=0;sp<spmax-1;sp++)
							{
								b[sp+2][0] = 0.25*(Ys[ii][jj+1][sp]-Ys[ii-1][jj+1][sp])*byhs;
								b[sp+2][1] = (Ys[ii][jj+1][sp]-Ys[ii][jj][sp])/(rgm[jj]-rgm[jj-1]);
								b[sp+2][2] = Ys[ii][jj+1][sp] - b[sp+2][0]*((double)ii+0.5)*hs - b[sp+2][1]*rgm[jj];
							}}}
				}


					
					Tgc[i][j] = b[0][0]*s0+b[0][1]*r0+b[0][2];
					Pc[i][j] = b[1][0]*s0+b[1][1]*r0+b[1][2];
					for(sp=0;sp<spmax-1;sp++)
					{
						Yc[i][j][sp] = b[sp+2][0]*s0+b[sp+2][1]*r0+b[sp+2][2];
					}
					Yc[i][j][nN2] = N2sum(Yc[i][j]);
					propertygc_cal(i,j);
					rogc[i][j] = rogcal(Mmc[i][j],Tgc[i][j],Pa);
				}
			}
		}
	}
	/* z Vector */
	for(i=1;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			r0 = sqrt(xm[j]*xm[j]+(z[i]-ld)*(z[i]-ld));
			if(r0<rgm[Ng-3])
			{
				czflg[i][j] = 1;
				if(r0>rgm[0])
				{
					czflg[i][j] = 2;
					s0 = acos((z[i]-ld)/r0);
					flg = 0;
					k = 0;
					do
					{
						if(rgm[k]>r0)
						{
							flg = 1;
							jj = k;
						}
						k+=1;
					}while(flg==0);
					flg = 0;
					k = 0;
					do
					{
						if(((double)k+0.5)*hs>s0)
						{
							flg = 1;
							ii = k;
						}
						k+=1;
					}while(flg==0&&k<NQ);
					if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.75)*hs;

					   

					if(aa>r0)
					{    u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
						v1 = u2*scm[ii]-u1*ssm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj-1]-ur[ii-1][jj-1]);
						v3 = u2*sc[NQ]-u1*ss[NQ];
							
						
						if(bb>s0){u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
							v2 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];}

					    else{u1 = 0.0;
							u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj+1]-ur[ii-1][jj+1]);
							v4 = u2*sc[NQ]-u1*ss[NQ];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
						u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
						v2 = u2*scm[ii]-u1*ssm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[ii][jj+1]-ur[ii-1][jj+1]+9.0*ur[ii][jj]-ur[ii-1][jj]);
						v4 = u2*sc[NQ]-u1*ss[NQ];
												
						if(bb>s0){u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
							v1 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];}
					     else{u1 = 0.0;
							u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj-1]-ur[ii-1][jj-1]);
							v3 = u2*sc[NQ]-u1*ss[NQ];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj];}}
				}
				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.25*hs;
					if(aa>r0)
					{  
						
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj-1]-ur[1][jj-1]);
						v1 = u2*sc[0]-u1*ss[0];
						u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
						v3 = u2*scm[ii]-u1*ssm[ii];
												
						if(bb>s0){u1 = 0.0;
							u2 = 0.0625*(9.0*ur[0][jj+1]-ur[1][jj+1]+9.0*ur[0][jj]-ur[1][jj]);
							v2 = u2*sc[0]-u1*ss[0];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj-1];}
					    else{u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*0.5*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.0;
						u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj+1]-ur[1][jj+1]);
						v2 = u2*sc[0]-u1*ss[0];
						u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
						u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
						v4 = u2*scm[ii]-u1*ssm[ii];
						if(bb>s0){
							u1 = 0.0;
							u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj-1]-ur[1][jj-1]);
							v1 = u2*sc[0]-u1*ss[0];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*0.5*hs - b[0][1]*rgm[jj];}
					     else{u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
							v3 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*0.5*hs - b[0][1]*rgm[jj];}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa>r0)
					{   u1 = 0.5*(uq[ii][jj]+uq[ii-1][jj]);
						u2 = 0.5*(ur[ii-1][jj]+ur[ii-1][jj-1]);
						v1 = u2*scm[ii-1]-u1*ssm[ii-1];
						u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj-1]+ur[ii][jj]);
						v3 = u2*scm[ii]-u1*ssm[ii];
						
						
						if(bb>s0){u1 = 0.5*(uq[ii][jj+1]+uq[ii-1][jj+1]);
							u2 = 0.5*(ur[ii-1][jj+1]+ur[ii-1][jj]);
							v2 = u2*scm[ii-1]-u1*ssm[ii-1];
							b[0][0] = (v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj-1];}
					    else{u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = (v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(uq[ii][jj+1]+uq[ii-1][jj+1]);
							u2 = 0.5*(ur[ii-1][jj+1]+ur[ii-1][jj]);
							v2 = u2*scm[ii-1]-u1*ssm[ii-1];
						u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
						
						
						if(bb>s0){ u1 = 0.5*(uq[ii][jj]+uq[ii-1][jj]);
						u2 = 0.5*(ur[ii-1][jj]+ur[ii-1][jj-1]);
						v1 = u2*scm[ii-1]-u1*ssm[ii-1];
							b[0][0] = (v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj];}
					     else{
							u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						  u2 = 0.5*(ur[ii][jj-1]+ur[ii][jj]);
						 v3 = u2*scm[ii]-u1*ssm[ii];
							b[0][0] = (v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];}}
				}
					uz[i][j] = b[0][0]*s0+b[0][1]*r0+b[0][2];

				}
			}
		}
	}
	/* x Vector */
	for(i=0;i<NZ;i++)
	{
		for(j=1;j<NX;j++)
		{
			r0 = sqrt(x[j]*x[j]+(zm[i]-ld)*(zm[i]-ld));
			if(r0<rgm[Ng-3])
			{
				cxflg[i][j] = 1;
				if(r0>rgm[0])
				{
					cxflg[i][j] = 2;
					s0 = acos((zm[i]-ld)/r0);
					flg = 0;
					k = 0;
					do
					{
						if(rgm[k]>r0)
						{
							flg = 1;
							jj = k;
						}
						k+=1;
					}while(flg==0);
					flg = 0;
					k = 0;
					do
					{
						if(((double)k+0.5)*hs>s0)
						{
							flg = 1;
							ii = k;
						}
						k+=1;
					}while(flg==0&&k<NQ);


					if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.75)*hs;

					   

					if(aa>r0)
					{    u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
						v1 =  u2*ssm[ii]+u1*scm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj-1]-ur[ii-1][jj-1]);
						v3 =  u2*ss[NQ]+u1*sc[NQ];
							
						
						if(bb>s0){u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
							v2 =  u2*ssm[ii]+u1*scm[ii];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];}

					    else{u1 = 0.0;
							u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj+1]-ur[ii-1][jj+1]);
							v4 =  u2*ss[NQ]+u1*sc[NQ];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
						u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
						v2 =  u2*ssm[ii]+u1*scm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[ii][jj+1]-ur[ii-1][jj+1]+9.0*ur[ii][jj]-ur[ii-1][jj]);
						v4 =  u2*ss[NQ]+u1*sc[NQ];
												
						if(bb>s0){u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
							v1 =  u2*ssm[ii]+u1*scm[ii];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];}
					     else{u1 = 0.0;
							u2 = 0.0625*(9.0*ur[ii][jj]-ur[ii-1][jj]+9.0*ur[ii][jj-1]-ur[ii-1][jj-1]);
							v3 =  u2*ss[NQ]+u1*sc[NQ];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj];}}
				}
				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.25*hs;
					if(aa>r0)
					{  
						
						u1 = 0.0;
						u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj-1]-ur[1][jj-1]);
						v1 =  u2*ss[ii]+u1*sc[ii];
						u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
						v3 =  u2*ssm[ii]+u1*scm[ii];
												
						if(bb>s0){u1 = 0.0;
							u2 = 0.0625*(9.0*ur[0][jj+1]-ur[1][jj+1]+9.0*ur[0][jj]-ur[1][jj]);
							v2 =  u2*ss[ii]+u1*sc[ii];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*(double)(ii+0.5)*hs - b[0][1]*rgm[jj-1];}
					    else{u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 = u2*ssm[ii]+u1*scm[ii];
							b[0][0] = 2.0*(v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*0.5*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.0;
						u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj+1]-ur[1][jj+1]);
						v2 =  u2*ss[ii]+u1*sc[ii];
						u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
						u2 = 0.5*(ur[ii][jj+1]+ur[ii][jj]);
						v4 = u2*ssm[ii]+u1*scm[ii];
						if(bb>s0){
							u1 = 0.0;
							u2 = 0.0625*(9.0*ur[0][jj]-ur[1][jj]+9.0*ur[0][jj-1]-ur[1][jj-1]);
							v1 = u2*ss[ii]+u1*sc[ii];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*0.5*hs - b[0][1]*rgm[jj];}
					     else{u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj-1]);
							v3 =  u2*ssm[ii]+u1*scm[ii];
							b[0][0] = 2.0*(v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*0.5*hs - b[0][1]*rgm[jj];}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa>r0)
					{   u1 = 0.5*(uq[ii][jj]+uq[ii-1][jj]);
						u2 = 0.5*(ur[ii-1][jj]+ur[ii-1][jj-1]);
						v1 =  u2*ssm[ii-1]+u1*scm[ii-1];
						u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						u2 = 0.5*(ur[ii][jj-1]+ur[ii][jj]);
						v3 =  u2*ssm[ii]+u1*scm[ii];
						
						
						if(bb>s0){u1 = 0.5*(uq[ii][jj+1]+uq[ii-1][jj+1]);
							u2 = 0.5*(ur[ii-1][jj+1]+ur[ii-1][jj]);
							v2 =  u2*ssm[ii-1]+u1*scm[ii-1];
							b[0][0] = (v3-v1)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v1 - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj-1];}
					    else{u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 =  u2*ssm[ii]+u1*scm[ii];
							b[0][0] = (v3-v1)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v3 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(uq[ii][jj+1]+uq[ii-1][jj+1]);
							u2 = 0.5*(ur[ii-1][jj+1]+ur[ii-1][jj]);
							v2 =  u2*ssm[ii-1]+u1*scm[ii-1];
						u1 = 0.5*(uq[ii+1][jj+1]+uq[ii][jj+1]);
							u2 = 0.5*(ur[ii][jj]+ur[ii][jj+1]);
							v4 =  u2*ssm[ii]+u1*scm[ii];
						
						
						if(bb>s0){ u1 = 0.5*(uq[ii][jj]+uq[ii-1][jj]);
						u2 = 0.5*(ur[ii-1][jj]+ur[ii-1][jj-1]);
						v1 =  u2*ssm[ii-1]+u1*scm[ii-1];
							b[0][0] = (v4-v2)*byhs;
							b[0][1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v2 - b[0][0]*((double)ii-0.5)*hs - b[0][1]*rgm[jj];}
					     else{
							u1 = 0.5*(uq[ii+1][jj]+uq[ii][jj]);
						  u2 = 0.5*(ur[ii][jj-1]+ur[ii][jj]);
						 v3 =  u2*ssm[ii]+u1*scm[ii];
							b[0][0] = (v4-v2)*byhs;
							b[0][1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[0][2] = v4 - b[0][0]*((double)ii+0.5)*hs - b[0][1]*rgm[jj];}}
				}


					

					ux[i][j] = b[0][0]*s0+b[0][1]*r0+b[0][2];
				}
			}
		}
	}
	freedmat(b,spmax+1,3);

}

void pint()
{
	int i,j,sp;
	double PI;

	Rs = Rs0;

	PI = 4.0*atan(1.0);

	hrl = 1.0 / (double)Nl;
	hrg = 1.0 / (double)Ng;
	hs = PI / (double)NQ;
	hx = Xmax / (double)NX;
	hz = Lmax / (double)NZ;
	byhrl = 1.0 / hrl;
	byhrg = 1.0 / hrg;
	byhs = 1.0 / hs;
	byhx = 1.0 / hx;
	byhz = 1.0 / hz;

	/* cylindrical coordinate system gas phase */
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			Tgc[i][j] = qTgc[i][j] = Ta;
			Pc[i][j] = dPc[i][j] = 0.0;
			for(sp=0;sp<spmax;sp++)
			{
				Yc[i][j][sp] = qYc[i][j][sp] = Y0[sp];
			}
			propertygc_cal(i,j);
			qrogc[i][j] = rogc[i][j] = rogcal(Mmc[i][j],Tgc[i][j],Pa);
			dPc[i][j] = 0.0;
			ux[i][j] = qux[i][j] = 0.0;
			uz[i][j] = quz[i][j] = 0.0;
		}
	}
	for(j=0;j<=NX;j++)
	{
		x[j] = (double)j*hx;
	}
	for(j=0;j<NX;j++)
	{
		xm[j] = ((double)j+0.5)*hx;
	}
	for(i=0;i<=NZ;i++)
	{
		z[i] = (double)i*hz;
	}
	for(i=0;i<NZ;i++)
	{
		zm[i] = ((double)i+0.5)*hz;
	}

	/* spherical coordinate system gas phase */
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<=Ng;j++)
		{
			Tgs[i][j] = qTgs[i][j] = Ta;
			Ps[i][j] = dPs[i][j] = 0.0;
			for(sp=0;sp<spmax;sp++)
			{
				Ys[i][j][sp] = qYs[i][j][sp] = Y0[sp];
			}
			propertygs_cal(i,j);
			rogs[i][j] = qrogs[i][j] = rogcal(Mms[i][j],Tgs[i][j],Pa);
			ur[i][j] = 0.0;
			qur[i][j] = 0.0;
		}
	}
	for(i=0;i<=NQ;i++)
	{
		for(j=0;j<=Ng;j++)
		{
			uq[i][j] = 0.0;
			quq[i][j] = 0.0;
		}
	}
	rgcal();
	for(i=0;i<=Ng;i++) eta[i] = (double)i*hrg;
	for(i=0;i<Ng;i++) etam[i] = ((double)i+0.5)*hrg;

	/* spherical coordinate system liquid phase */
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl;j++)
		{
			Tls[i][j] = Tl0;
			propertyl_cal(i,j);
		}
		
		Tls[i][Nl] = Tl0;
		propertyl_cal(i,Nl);
	}
	rlcal();
	for(i=0;i<NQ;i++)
	{
		ssm[i] = sin(((double)i+0.5)*hs);
		scm[i] = cos(((double)i+0.5)*hs);
	}
	for(i=0;i<=NQ;i++)
	{
		ss[i] = sin((double)i*hs);
		sc[i] = cos((double)i*hs);
	}
	for(i=0;i<Nl;i++)
	{
		xim[i] = ((double)i+0.5)*hrl;
	}
	for(i=0;i<=Nl;i++)
	{
		xi[i] = (double)i*hrl;
	}

	Search();
}

/* liquid phase temperature calculate */
void Tl_cal() /*Completed*/
{
	int i,j,ii=NQ-1,jj=Nl-1;
	double RsdotbyRs,A,B;

	RsdotbyRs = Rdot / Rs;

	/*i=0&&j=0*/
	qTls[0][0] = (0.5*rol*Cpl[0][0]*RsdotbyRs*xi[1]*(Tls[0][1]-Tls[0][0])*byhrl
				 +0.5*rl[1]*rl[1]*(ll[0][1]+ll[0][0])*(Tls[0][1]-Tls[0][0])*byhrl*byhrl/(rlm[0]*rlm[0]*Rs*Rs)
				 +0.5*ss[1]*(ll[1][0]+ll[0][0])*(Tls[1][0]-Tls[0][0])*byhs*byhs/(rlm[0]*rlm[0]*ssm[0]))/(rol*Cpl[0][0])*dt;
	qTls[0][0]+= Tls[0][0];

	/*i=NQ-1&&j=0*/
	qTls[ii][0] = ((0.5*rol*Cpl[ii][0]*RsdotbyRs*xi[1]*(Tls[ii][1]-Tls[ii][0])*byhrl
				 +0.5*rl[1]*rl[1]*(ll[ii][1]+ll[ii][0])*(Tls[ii][1]-Tls[ii][0])*byhrl*byhrl/(rlm[0]*rlm[0]*Rs*Rs)
				 -0.5*ss[ii]*(ll[ii][0]+ll[ii-1][0])*(Tls[ii][0]-Tls[ii-1][0])*byhs*byhs/(rlm[0]*rlm[0]*ssm[ii]))/(rol*Cpl[ii][0]))*dt;
	qTls[ii][0]+= Tls[ii][0];

	/*i=0&&j=Nl-1*/
	qTls[0][jj] = ((0.5*rol*Cpl[0][jj]*RsdotbyRs*(xi[jj+1]*(8.0*Tls[0][jj+1]-9.0*Tls[0][jj]+Tls[0][jj-1])/3.0
		                                        +xi[jj]*(Tls[0][jj]-Tls[0][jj-1]))*byhrl
			      +0.5*(2.0*rl[jj+1]*rl[jj+1]*ll[0][jj+1]*(8.0*Tls[0][jj+1]-9.0*Tls[0][jj]+Tls[0][jj-1])/3.0
					       -rl[jj]*rl[jj]*(ll[0][jj]+ll[0][jj-1])*(Tls[0][jj]-Tls[0][jj-1]))*byhrl*byhrl/(rlm[jj]*rlm[jj]*Rs*Rs)
				  +0.5*ss[1]*(ll[1][jj]+ll[0][jj])*(Tls[1][jj]-Tls[0][jj])*byhs*byhs/(rlm[jj]*rlm[jj]*ssm[0]))/(rol*Cpl[0][jj]))*dt;
	qTls[0][jj]+= Tls[0][jj];

	/*i=NQ-1&&j=Nl-1*/
	qTls[ii][jj] = ((0.5*rol*Cpl[ii][jj]*RsdotbyRs*(xi[jj+1]*(8.0*Tls[ii][jj+1]-9.0*Tls[ii][jj]+Tls[ii][jj-1])/3.0
		                                        +xi[jj]*(Tls[ii][jj]-Tls[ii][jj-1]))*byhrl
			      +0.5*(2.0*rl[jj+1]*rl[jj+1]*ll[ii][jj+1]*(8.0*Tls[ii][jj+1]-9.0*Tls[ii][jj]+Tls[ii][jj-1])/3.0
					       -rl[jj]*rl[jj]*(ll[ii][jj]+ll[ii][jj-1])*(Tls[ii][jj]-Tls[ii][jj-1]))*byhrl*byhrl/(rlm[jj]*rlm[jj]*Rs*Rs)
				  -0.5*ss[ii]*(ll[ii][jj]+ll[ii-1][jj])*(Tls[ii][jj]-Tls[ii-1][jj])*byhs*byhs/(rlm[jj]*rlm[jj]*ssm[ii]))/(rol*Cpl[ii][jj]))*dt;
	qTls[ii][jj]+= Tls[0][jj];

	for(i=1;i<NQ-1;i++)
	{
		/* spherical center */
		qTls[i][0] = ((0.5*rol*Cpl[i][0]*RsdotbyRs*xi[1]*(Tls[i][1]-Tls[i][0])*byhrl
					 +0.5*rl[1]*rl[1]*(ll[i][1]+ll[i][0])*(Tls[i][1]-Tls[i][0])*byhrl*byhrl/(rlm[0]*rlm[0]*Rs*Rs)
					 +0.5*(ss[i+1]*(ll[i+1][0]+ll[i][0])*(Tls[i+1][0]-Tls[i][0])
				          -ss[i]*(ll[i][0]+ll[i-1][0])*(Tls[i][0]-Tls[i-1][0]))*byhs*byhs/(rlm[0]*rlm[0]*ssm[i]))/(rol*Cpl[i][0]))*dt;
		qTls[i][0]+= Tls[i][0];

		/* vapor-liquid interface */
		qTls[i][jj] = ((0.5*rol*Cpl[i][jj]*RsdotbyRs*(xi[jj+1]*(8.0*Tls[i][jj+1]-9.0*Tls[i][jj]+Tls[i][jj-1])/3.0
			                                        +xi[jj]*(Tls[i][jj]-Tls[i][jj-1]))*byhrl
				      +0.5*(2.0*rl[jj+1]*rl[jj+1]*ll[i][jj+1]*(8.0*Tls[i][jj+1]-9.0*Tls[i][jj]+Tls[i][jj-1])/3.0
						       -rl[jj]*rl[jj]*(ll[i][jj]+ll[i][jj-1])*(Tls[i][jj]-Tls[i][jj-1]))*byhrl*byhrl/(rlm[jj]*rlm[jj]*Rs*Rs)
					  +0.5*(ss[i+1]*(ll[i+1][jj]+ll[i][jj])*(Tls[i+1][jj]-Tls[i][jj])
						   -ss[i]*(ll[i][jj]+ll[i-1][jj])*(Tls[i][jj]-Tls[i-1][jj]))*byhs*byhs/(rlm[jj]*rlm[jj]*ssm[i]))/(rol*Cpl[i][jj]))*dt;
		qTls[i][jj]+= Tls[i][jj];
	}

	for(j=1;j<Nl-1;j++)
	{
		/* i = 0 */
		qTls[0][j] = ((0.5*rol*Cpl[0][j]*RsdotbyRs*(xi[j+1]*(Tls[0][j+1]-Tls[0][j])+xi[j]*(Tls[0][j]-Tls[0][j-1]))*byhrl
					 +0.5*(rl[j+1]*rl[j+1]*(ll[0][j+1]+ll[0][j])*(Tls[0][j+1]-Tls[0][j])
						  -rl[j]*rl[j]*(ll[0][j]+ll[0][j-1])*(Tls[0][j]-Tls[0][j-1]))*byhrl*byhrl/(rlm[j]*rlm[j]*Rs*Rs)
					 +0.5*ss[1]*(ll[1][j]+ll[0][j])*(Tls[1][j]-Tls[0][j])*byhs*byhs/(rlm[j]*rlm[j]*ssm[0]))/(rol*Cpl[0][j]))*dt;
		qTls[0][j]+= Tls[0][j];

		/* i = NQ-1 */
		qTls[ii][j] = ((0.5*rol*Cpl[ii][j]*RsdotbyRs*(xi[j+1]*(Tls[ii][j+1]-Tls[ii][j])+xi[j]*(Tls[ii][j]-Tls[ii][j-1]))*byhrl
					 +0.5*(rl[j+1]*rl[j+1]*(ll[ii][j+1]+ll[ii][j])*(Tls[ii][j+1]-Tls[ii][j])
					      -rl[j]*rl[j]*(ll[ii][j]+ll[ii][j-1])*(Tls[ii][j]-Tls[ii][j-1]))*byhrl*byhrl/(rlm[j]*rlm[j]*Rs*Rs)
				     -0.5*ss[ii]*(ll[ii][j]+ll[ii-1][j])*(Tls[ii][j]-Tls[ii-1][j])*byhs*byhs/(rlm[j]*rlm[j]*ssm[ii]))/(rol*Cpl[ii][j]))*dt;
		qTls[ii][j]+= Tls[ii][j];
	}
	for(i=1;i<NQ-1;i++)
	{
		for(j=1;j<Nl-1;j++)
		{
			qTls[i][j] = ((0.5*rol*Cpl[i][j]*RsdotbyRs*(xi[j+1]*(Tls[i][j+1]-Tls[i][j])+xi[j]*(Tls[i][j]-Tls[i][j-1]))*byhrl
						 +0.5*(rl[j+1]*rl[j+1]*(ll[i][j+1]+ll[i][j])*(Tls[i][j+1]-Tls[i][j])
						      -rl[j]*rl[j]*(ll[i][j]+ll[i][j-1])*(Tls[i][j]-Tls[i][j-1]))*byhrl*byhrl/(rlm[j]*rlm[j]*Rs*Rs)
					     +0.5*(ss[i+1]*(ll[i+1][j]+ll[i][j])*(Tls[i+1][j]-Tls[i][j])
						      -ss[i]*(ll[i][j]+ll[i-1][j])*(Tls[i][j]-Tls[i-1][j]))*byhs*byhs/(rlm[j]*rlm[j]*ssm[i]))/(rol*Cpl[i][j]))*dt;
			qTls[i][j]+= Tls[i][j];
		}
	}

}

/* spherical coordinate region density calculate */
void Ys_cal()
{
	int i,j,ii=NQ-1,sp;
	double AA,RdotbyRs,byRs;

	byRs = 1.0/Rs;
	RdotbyRs = Rdot*byRs;
	AA = 1.0 / log(rmax/Rs);

	for(sp=0;sp<spmax-1;sp++)
	{
		/* i=0&&j=1 */		
		qYs[0][1][sp] = (0.5*AA*RdotbyRs*rogs[0][1]*((1.0-eta[1])*(Ys[0][2][sp]-Ys[0][1][sp])
												   +(-8.0*Ys[0][0][sp]+9.0*Ys[0][1][sp]-Ys[0][2][sp])/3.0)*byhrg
					   //-0.25*AA*rogs[0][1]*(ur[0][1]+ur[0][0])*((Ys[0][2][sp]+Ys[0][1][sp])-2.0*Ys[0][0][sp])*byhrg/rgm[0]
					   //-0.25*rogs[0][1]*(uq[1][1]+uq[0][1])*(Ys[1][1][sp]-Ys[0][1][sp])*byhs/rgm[0]
					   -0.25*AA*rogs[0][1]*byhrg/rgm[0]*
					   ((fabs(ur[0][1]+ur[0][0])+(ur[0][1]+ur[0][0]))*(8.0*(Ys[0][1][sp]-Ys[0][0][sp])+Ys[0][1][sp]-Ys[0][2][sp])/3.0
					   -(fabs(ur[0][1]+ur[0][0])-(ur[0][1]+ur[0][0]))*(Ys[0][2][sp]-Ys[0][1][sp]))
					   -0.25*rogs[0][1]*byhs/rgm[0]*
					   ((fabs(uq[1][1]+uq[0][1])+(uq[1][1]+uq[0][1]))*(Ys[0][1][sp]-Ys[0][1][sp])
					   -(fabs(uq[1][1]+uq[0][1])-(uq[1][1]+uq[0][1]))*(Ys[1][1][sp]-Ys[0][1][sp]))
					   +0.25*AA*AA*(rg[1]*(rogs[0][2]+rogs[0][1])*(Dgs[0][2][sp]+Dgs[0][1][sp])*(Ys[0][2][sp]-Ys[0][1][sp])
						-4.0*rg[0]*rogs[0][0]*Dgs[0][0][sp]*(-8.0*Ys[0][0][sp]+9.0*Ys[0][1][sp]-Ys[0][2][sp])/3.0)
						*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
					   +0.25*ss[1]*(rogs[1][1]+rogs[0][1])*(Dgs[1][1][sp]+Dgs[0][1][sp])*(Ys[1][1][sp]-Ys[0][1][sp])
					   *byhs*byhs/(rgm[0]*ssm[0])+omegas[0][1][sp])/rogs[0][1]*dt;
		qYs[0][1][sp]+= Ys[0][1][sp];
	
	   	/* i=NQ-1&&j=1 */
		qYs[ii][1][sp] = (0.5*AA*RdotbyRs*rogs[ii][1]*((1.0-eta[1])*(Ys[ii][2][sp]-Ys[ii][1][sp])
													 +(-8.0*Ys[ii][0][sp]+9.0*Ys[ii][1][sp]-Ys[ii][2][sp])/3.0)*byhrg
					    //-0.25*AA*rogs[ii][1]*(ur[ii][1]+ur[ii][0])*((Ys[ii][2][sp]+Ys[ii][1][sp])-2.0*Ys[ii][0][sp])*byhrg/rgm[0]
					    //-0.25*rogs[ii][1]*(uq[ii+1][1]+uq[ii][1])*(Ys[ii][1][sp]-Ys[ii-1][1][sp])*byhs/rgm[0]
					    -0.25*AA*rogs[ii][1]*byhrg/rgm[0]*
					    ((fabs(ur[ii][1]+ur[ii][0])+(ur[ii][1]+ur[ii][0]))*(8.0*(Ys[ii][1][sp]-Ys[ii][0][sp])+Ys[ii][1][sp]-Ys[ii][2][sp])/3.0
					    -(fabs(ur[ii][1]+ur[ii][0])-(ur[ii][1]+ur[ii][0]))*(Ys[ii][2][sp]-Ys[ii][1][sp]))
					    -0.25*rogs[ii][1]*byhs/rgm[0]*
					    ((fabs(uq[ii+1][1]+uq[ii][1])+(uq[ii+1][1]+uq[ii][1]))*(Ys[ii][1][sp]-Ys[ii-1][1][sp])
					    -(fabs(uq[ii+1][1]+uq[ii][1])-(uq[ii+1][1]+uq[ii][1]))*(Ys[ii][1][sp]-Ys[ii][1][sp]))
					    +0.25*AA*AA*(rg[1]*(rogs[ii][2]+rogs[ii][1])*(Dgs[ii][2][sp]+Dgs[ii][1][sp])*(Ys[ii][2][sp]-Ys[ii][1][sp])
						-4.0*rg[0]*rogs[ii][0]*Dgs[ii][0][sp]*(-8.0*Ys[ii][0][sp]+9.0*Ys[ii][1][sp]-Ys[ii][2][sp])/3.0)
						*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
					   -0.25*ss[ii]*(rogs[ii][1]+rogs[ii-1][1])*(Dgs[ii][1][sp]+Dgs[ii-1][1][sp])*(Ys[ii][1][sp]-Ys[ii-1][1][sp])
					   *byhs*byhs/(rgm[0]*ssm[ii])+omegas[ii][1][sp])/rogs[ii][1]*dt;
		qYs[ii][1][sp]+= Ys[ii][1][sp];

		/* vapor-liquid interface */
		for(i=1;i<NQ-1;i++)
		{
			qYs[i][1][sp] = (0.5*AA*RdotbyRs*rogs[i][1]*((1.0-eta[1])*(Ys[i][2][sp]-Ys[i][1][sp])
													   +(-8.0*Ys[i][0][sp]+9.0*Ys[i][1][sp]-Ys[i][2][sp])/3.0)*byhrg
						   //-0.25*AA*rogs[i][1]*(ur[i][1]+ur[i][0])*((Ys[i][2][sp]+Ys[i][1][sp])-2.0*Ys[i][0][sp])*byhrg/rgm[0]
						   //-0.25*rogs[i][1]*(uq[i+1][1]+uq[i][1])*(Ys[i+1][1][sp]-Ys[i-1][1][sp])*byhs/rgm[0]
						   -0.25*AA*rogs[i][1]*byhrg/rgm[0]*
						   ((fabs(ur[i][1]+ur[i][0])+(ur[i][1]+ur[i][0]))*(8.0*(Ys[i][1][sp]-Ys[i][0][sp])+Ys[i][1][sp]-Ys[i][2][sp])/3.0
						   -(fabs(ur[i][1]+ur[i][0])-(ur[i][1]+ur[i][0]))*(Ys[i][2][sp]-Ys[i][1][sp]))
						   -0.25*rogs[i][1]*byhs/rgm[0]*
						   ((fabs(uq[i+1][1]+uq[i][1])+(uq[i+1][1]+uq[i][1]))*(Ys[i][1][sp]-Ys[i-1][1][sp])
						   -(fabs(uq[i+1][1]+uq[i][1])-(uq[i+1][1]+uq[i][1]))*(Ys[i+1][1][sp]-Ys[i][1][sp]))
						   +0.25*AA*AA*(rg[1]*(rogs[i][2]+rogs[i][1])*(Dgs[i][2][sp]+Dgs[i][1][sp])*(Ys[i][2][sp]-Ys[i][1][sp])
							-4.0*rg[0]*rogs[i][0]*Dgs[i][0][sp]*(-8.0*Ys[i][0][sp]+9.0*Ys[i][1][sp]-Ys[i][2][sp])/3.0)
							*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
						   +0.25*(ss[i+1]*(rogs[i+1][1]+rogs[i][1])*(Dgs[i+1][1][sp]+Dgs[i][1][sp])*(Ys[i+1][1][sp]-Ys[i][1][sp])
						    -ss[i]*(rogs[i][1]+rogs[i-1][1])*(Dgs[i][1][sp]+Dgs[i-1][1][sp])*(Ys[i][1][sp]-Ys[i-1][1][sp]))
							*byhs*byhs/(rgm[0]*ssm[i])+omegas[i][1][sp])/rogs[i][1]*dt;
			qYs[i][1][sp]+= Ys[i][1][sp];
		}

		for(j=2;j<Ng-1;j++)
		{
			/*i=0*/
			qYs[0][j][sp] = (0.5*AA*RdotbyRs*rogs[0][j]*((1.0-eta[j])*(Ys[0][j+1][sp]-Ys[0][j][sp])
														+(1.0-eta[j-1])*(Ys[0][j][sp]-Ys[0][j-1][sp]))*byhrg
						   //-0.25*AA*rogs[0][j]*(ur[0][j]+ur[0][j-1])*(Ys[0][j+1][sp]-Ys[0][j-1][sp])*byhrg/rgm[j-1]
						   //-0.25*rogs[0][j]*(uq[1][j]+uq[0][j])*(Ys[1][j][sp]-Ys[0][j][sp])*byhs/rgm[j-1]
						   -0.25*AA*rogs[0][j]*byhrg/rgm[j-1]*
						   ((fabs(ur[0][j]+ur[0][j-1])+(ur[0][j]+ur[0][j-1]))*(Ys[0][j][sp]-Ys[0][j-1][sp])
						   -(fabs(ur[0][j]+ur[0][j-1])-(ur[0][j]+ur[0][j-1]))*(Ys[0][j+1][sp]-Ys[0][j][sp]))
						   -0.25*rogs[0][j]*byhs/rgm[j-1]*
						   ((fabs(uq[1][j]+uq[0][j])+(uq[1][j]+uq[0][j]))*(Ys[0][j][sp]-Ys[0][j][sp])
						   -(fabs(uq[1][j]+uq[0][j])-(uq[1][j]+uq[0][j]))*(Ys[1][j][sp]-Ys[0][j][sp]))
						   +0.25*AA*AA*(rg[j]*(rogs[0][j+1]+rogs[0][j])*(Dgs[0][j+1][sp]+Dgs[0][j][sp])*(Ys[0][j+1][sp]-Ys[0][j][sp])
									   -rg[j-1]*(rogs[0][j]+rogs[0][j-1])*(Dgs[0][j][sp]+Dgs[0][j-1][sp])*(Ys[0][j][sp]-Ys[0][j-1][sp]))
									   *byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
						   +0.25*ss[1]*(rogs[1][j]+rogs[0][j])*(Dgs[1][j][sp]+Dgs[0][j][sp])*(Ys[1][j][sp]-Ys[0][j][sp])
								 *byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[0])+omegas[0][j][sp])/rogs[0][j]*dt;
			qYs[0][j][sp]+= Ys[0][j][sp];

			/* i=NQ-1 */
			qYs[ii][j][sp] = (0.5*AA*RdotbyRs*rogs[ii][j]*((1.0-eta[j])*(Ys[ii][j+1][sp]-Ys[ii][j][sp])
														  +(1.0-eta[j-1])*(Ys[ii][j][sp]-Ys[ii][j-1][sp]))*byhrg
						   //-0.25*AA*rogs[ii][j]*(ur[ii][j]+ur[ii][j-1])*(Ys[ii][j+1][sp]-Ys[ii][j-1][sp])*byhrg/rgm[j-1]
						   //-0.25*rogs[ii][j]*(uq[ii+1][j]+uq[ii][j])*(Ys[ii][j][sp]-Ys[ii-1][j][sp])*byhs/rgm[j-1]
						   -0.25*AA*rogs[ii][j]*byhrg/rgm[j-1]*
						   ((fabs(ur[ii][j]+ur[ii][j-1])+(ur[ii][j]+ur[ii][j-1]))*(Ys[ii][j][sp]-Ys[ii][j-1][sp])
						   -(fabs(ur[ii][j]+ur[ii][j-1])-(ur[ii][j]+ur[ii][j-1]))*(Ys[ii][j+1][sp]-Ys[ii][j][sp]))
						   -0.25*rogs[ii][j]*byhs/rgm[j-1]*
						   ((fabs(uq[ii+1][j]+uq[ii][j])+(uq[ii+1][j]+uq[ii][j]))*(Ys[ii][j][sp]-Ys[ii-1][j][sp])
						   -(fabs(uq[ii+1][j]+uq[ii][j])-(uq[ii+1][j]+uq[ii][j]))*(Ys[ii][j][sp]-Ys[ii][j][sp]))
						   +0.25*AA*AA*(rg[j]*(rogs[ii][j+1]+rogs[ii][j])*(Dgs[ii][j+1][sp]+Dgs[ii][j][sp])*(Ys[ii][j+1][sp]-Ys[ii][j][sp])
									   -rg[j-1]*(rogs[ii][j]+rogs[ii][j-1])*(Dgs[ii][j][sp]+Dgs[ii][j-1][sp])*(Ys[ii][j][sp]-Ys[ii][j-1][sp]))
									   *byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
						   -0.25*ss[ii]*(rogs[ii][j]+rogs[ii-1][j])*(Dgs[ii][j][sp]+Dgs[ii-1][j][sp])*(Ys[ii][j][sp]-Ys[ii-1][j][sp])
								 *byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[ii])+omegas[ii][j][sp])/rogs[ii][j]*dt;
			qYs[ii][j][sp]+= Ys[ii][j][sp];
		}
		for(i=1;i<NQ-1;i++)
		{
			for(j=2;j<Ng-1;j++)
			{
				qYs[i][j][sp] = (0.5*AA*RdotbyRs*rogs[i][j]*((1.0-eta[j])*(Ys[i][j+1][sp]-Ys[i][j][sp])
															+(1.0-eta[j-1])*(Ys[i][j][sp]-Ys[i][j-1][sp]))*byhrg
							   //-0.25*AA*rogs[i][j]*(ur[i][j]+ur[i][j-1])*(Ys[i][j+1][sp]-Ys[i][j-1][sp])*byhrg/rgm[j-1]
							   //-0.25*rogs[i][j]*(uq[i+1][j]+uq[i][j])*(Ys[i+1][j][sp]-Ys[i-1][j][sp])*byhs/rgm[j-1]
							   -0.25*AA*rogs[i][j]*byhrg/rgm[j-1]*
							   ((fabs(ur[i][j]+ur[i][j-1])+(ur[i][j]+ur[i][j-1]))*(Ys[i][j][sp]-Ys[i][j-1][sp])
							   -(fabs(ur[i][j]+ur[i][j-1])-(ur[i][j]+ur[i][j-1]))*(Ys[i][j+1][sp]-Ys[i][j][sp]))
							   -0.25*rogs[i][j]*byhs/rgm[j-1]*
							   ((fabs(uq[i+1][j]+uq[i][j])+(uq[i+1][j]+uq[i][j]))*(Ys[i][j][sp]-Ys[i-1][j][sp])
							   -(fabs(uq[i+1][j]+uq[i][j])-(uq[i+1][j]+uq[i][j]))*(Ys[i+1][j][sp]-Ys[i][j][sp]))
							   +0.25*AA*AA*(rg[j]*(rogs[i][j+1]+rogs[i][j])*(Dgs[i][j+1][sp]+Dgs[i][j][sp])*(Ys[i][j+1][sp]-Ys[i][j][sp])
										   -rg[j-1]*(rogs[i][j]+rogs[i][j-1])*(Dgs[i][j][sp]+Dgs[i][j-1][sp])*(Ys[i][j][sp]-Ys[i][j-1][sp]))
										   *byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
							   +0.25*(ss[i+1]*(rogs[i+1][j]+rogs[i][j])*(Dgs[i+1][j][sp]+Dgs[i][j][sp])*(Ys[i+1][j][sp]-Ys[i][j][sp])
							         -ss[i]*(rogs[i][j]+rogs[i-1][j])*(Dgs[i][j][sp]+Dgs[i-1][j][sp])*(Ys[i][j][sp]-Ys[i-1][j][sp]))
									 *byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[i])+omegas[i][j][sp])/rogs[i][j]*dt;
				qYs[i][j][sp]+= Ys[i][j][sp];
			}
		}
	}

}

void Tgs_cal()/*Complete*/
{
	int i,j,ii=NQ-1,sp;
	double AA,RdotbyRs,byRs;
	double Asp,Bsp;
	double tau1,tau2,tau3,tau4,tau5;

	byRs = 1.0/Rs;
	RdotbyRs = Rdot*byRs;
	AA = 1.0 / log(rmax/Rs);

	/*i=0&&j=1*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp += Dgs[0][0][sp]*Cpgsi[0][0][sp]*((Ys[0][2][sp]+Ys[0][1][sp])-2.0*Ys[0][0][sp]);
		Bsp += Dgs[0][0][sp]*Cpgsi[0][0][sp]*(Ys[1][0][sp]-Ys[0][0][sp]);
	}
	Asp *= 0.25*rogs[0][0]*((Tgs[0][2]+Tgs[0][1])-2.0*Tgs[0][0])*AA*AA/(rgm[0]*rgm[0])*byhrg*byhrg;
	Bsp *= 0.25*rogs[0][0]*(Tgs[1][0]-Tgs[0][0])*byhs*byhs/(rgm[0]*rgm[0]);
	tau1 = AA*(ur[0][1]-ur[0][0])*byhrg/rgm[0];
	tau2 = (uq[1][1]-uq[0][1])*byhs/rgm[0]+0.5*(ur[0][1]+ur[0][0])/rgm[0];
	tau3 = 0.5*(ur[0][1]+ur[0][0])/rgm[0]+0.5*(uq[1][1]+uq[0][1])*scm[0]/(ssm[0]*rgm[0]);
	tau4 = 0.25*AA*((uq[1][2]+uq[0][2]+uq[1][1]+uq[0][1])-2.0*(uq[1][0]+uq[0][0]))*byhrg/rgm[0]
		  -0.5*(uq[1][1]+uq[0][1])/rgm[0]
		  +0.25*((ur[1][1]+ur[1][0])-(ur[0][1]+ur[0][0]))*byhs/rgm[0];
	tau5 = AA*(rg[1]*rg[1]*ur[0][1]-rg[0]*rg[0]*ur[0][0])*byhrg/(rgm[0]*rgm[0]*rgm[0])
		  +(ss[1]*uq[1][1]-ss[0]*uq[0][1])*byhs/(rgm[0]*ssm[0]);
	qTgs[0][1] = (0.5*AA*RdotbyRs*Cpgs[0][1]*rogs[0][1]*((1.0-eta[1])*(Tgs[0][2]-Tgs[0][1])
														+(-8.0*Tgs[0][0]+9.0*Tgs[0][1]-Tgs[0][2])/3.0)*byhrg
				//-0.25*AA*Cpgs[0][1]*rogs[0][1]*(ur[0][1]+ur[0][0])*((Tgs[0][2]+Tgs[0][1])-2.0*Tgs[0][0])*byhrg/rgm[0]
				//-0.25*Cpgs[0][1]*rogs[0][1]*(uq[1][1]+uq[0][1])*(Tgs[1][1]-Tgs[0][1])*byhs/rgm[0]
				-0.25*AA*Cpgs[0][1]*rogs[0][1]*byhrg/rgm[0]*
				((fabs(ur[0][1]+ur[0][0])+(ur[0][1]+ur[0][0]))*(8.0*(Tgs[0][1]-Tgs[0][0])+Tgs[0][1]-Tgs[0][2])/3.0
				-(fabs(ur[0][1]+ur[0][0])-(ur[0][1]+ur[0][0]))*(Tgs[0][2]-Tgs[0][1]))
				-0.25*Cpgs[0][1]*rogs[0][1]*byhs/rgm[0]*
				 ((fabs(uq[1][1]+uq[0][1])+(uq[1][1]+uq[0][1]))*(Tgs[0][1]-Tgs[0][1])
				 -(fabs(uq[1][1]+uq[0][1])-(uq[1][1]+uq[0][1]))*(Tgs[1][1]-Tgs[0][1]))
			    +Asp+Bsp-CgTs[0][1]
				//+mus[0][1]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
				+0.5*AA*AA*(rg[1]*(lgs[0][2]+lgs[0][1])*(Tgs[0][2]-Tgs[0][1])
						   -2.0*rg[0]*lgs[0][0]*(-8.0*Tgs[0][0]+9.0*Tgs[0][1]-Tgs[0][2])/3.0)*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
			    +0.5*ss[1]*(lgs[1][1]+lgs[0][1])*(Tgs[1][1]-Tgs[0][1])*byhs*byhs/(rgm[0]*rgm[0]*ssm[0]))/(rogs[0][1]*Cpgs[0][1])*dt;
	qTgs[0][1]+= Tgs[0][1];

	/*i=NQ-1&&j=1*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp += Dgs[ii][0][sp]*Cpgsi[ii][0][sp]*((Ys[ii][2][sp]+Ys[ii][1][sp])-2.0*Ys[ii][0][sp]);
		Bsp += Dgs[ii][0][sp]*Cpgsi[ii][0][sp]*(Ys[ii][0][sp]-Ys[ii-1][0][sp]);
	}
	Asp *= 0.25*rogs[ii][0]*((Tgs[ii][2]+Tgs[ii][1])-2.0*Tgs[ii][0])*AA*AA/(rgm[0]*rgm[0])*byhrg*byhrg;
	Bsp *= 0.25*rogs[ii][0]*(Tgs[ii][0]-Tgs[ii-1][0])*byhs*byhs/(rgm[0]*rgm[0]);
	tau1 = AA*(ur[ii][1]-ur[ii][0])*byhrg/rgm[0];
	tau2 = (uq[ii+1][1]-uq[ii][1])*byhs/rgm[0]+0.5*(ur[ii][1]+ur[ii][0])/rgm[0];
	tau3 = 0.5*(ur[ii][1]+ur[ii][0])/rgm[0]+0.5*(uq[ii+1][1]+uq[ii][1])*scm[ii]/(ssm[ii]*rgm[0]);
	tau4 = 0.25*AA*((uq[ii+1][2]+uq[ii][2]+uq[ii+1][1]+uq[ii][1])-2.0*(uq[ii+1][0]+uq[ii][0]))*byhrg/rgm[0]
		  -0.5*(uq[ii+1][1]+uq[ii][1])/rgm[0]
		  +0.25*((ur[ii][1]+ur[ii][0])-(ur[ii-1][1]+ur[ii-1][0]))*byhs/rgm[0];
	tau5 = AA*(rg[1]*rg[1]*ur[ii][1]-rg[0]*rg[0]*ur[ii][0])*byhrg/(rgm[0]*rgm[0]*rgm[0])
		  +(ss[ii+1]*uq[ii+1][1]-ss[ii]*uq[ii][1])*byhs/(rgm[0]*ssm[ii]);
	qTgs[ii][1] = (0.5*AA*RdotbyRs*Cpgs[ii][1]*rogs[ii][1]*((1.0-eta[1])*(Tgs[ii][2]-Tgs[ii][1])
														   +(-8.0*Tgs[ii][0]+9.0*Tgs[ii][1]-Tgs[ii][2])/3.0)*byhrg
				//-0.25*AA*Cpgs[ii][1]*rogs[ii][1]*(ur[ii][1]+ur[ii][0])*((Tgs[ii][2]+Tgs[ii][1])-2.0*Tgs[ii][0])*byhrg/rgm[0]
				//-0.25*Cpgs[ii][1]*rogs[ii][1]*(uq[ii+1][1]+uq[ii][1])*(Tgs[ii][1]-Tgs[ii-1][1])*byhs/rgm[0]
				-0.25*AA*Cpgs[ii][1]*rogs[ii][1]*byhrg/rgm[0]*
				((fabs(ur[ii][1]+ur[ii][0])+(ur[ii][1]+ur[ii][0]))*(8.0*(Tgs[ii][1]-Tgs[ii][0])+Tgs[ii][1]-Tgs[ii][2])/3.0
				-(fabs(ur[ii][1]+ur[ii][0])-(ur[ii][1]+ur[ii][0]))*(Tgs[ii][2]-Tgs[ii][1]))
				-0.25*Cpgs[ii][1]*rogs[ii][1]*byhs/rgm[0]*
				((fabs(uq[ii+1][1]+uq[ii][1])+(uq[ii+1][1]+uq[ii][1]))*(Tgs[ii][1]-Tgs[ii-1][1])
				-(fabs(uq[ii+1][1]+uq[ii][1])-(uq[ii+1][1]+uq[ii][1]))*(Tgs[ii][1]-Tgs[ii][1]))
			    +Asp+Bsp-CgTs[ii][1]
				//+mus[ii][1]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
				+0.5*AA*AA*(rg[1]*(lgs[ii][2]+lgs[ii][1])*(Tgs[ii][2]-Tgs[ii][1])
						   -2.0*rg[0]*lgs[ii][0]*(-8.0*Tgs[ii][0]+9.0*Tgs[ii][1]-Tgs[ii][2])/3.0)*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
			    -0.5*ss[ii]*(lgs[ii][1]+lgs[ii-1][1])*(Tgs[ii][1]-Tgs[ii-1][1])*byhs*byhs/(rgm[0]*rgm[0]*ssm[ii]))/(rogs[ii][1]*Cpgs[ii][1])*dt;
	qTgs[ii][1]+= Tgs[ii][1];

	/* vapor-liquid interface */
	for(i=1;i<NQ-1;i++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp += Dgs[i][0][sp]*Cpgsi[i][0][sp]*((Ys[i][2][sp]+Ys[i][1][sp])-2.0*Ys[i][0][sp]);
			Bsp += Dgs[i][0][sp]*Cpgsi[i][0][sp]*(Ys[i+1][0][sp]-Ys[i-1][0][sp]);
		}
		Asp *= 0.25*rogs[i][0]*((Tgs[i][2]+Tgs[i][1])-2.0*Tgs[i][0])*AA*AA/(rgm[0]*rgm[0])*byhrg*byhrg;
		Bsp *= 0.25*rogs[i][0]*(Tgs[i+1][0]-Tgs[i-1][0])*byhs*byhs/(rgm[0]*rgm[0]);
		tau1 = AA*(ur[i][1]-ur[i][0])*byhrg/rgm[0];
		tau2 = (uq[i+1][1]-uq[i][1])*byhs/rgm[0]+0.5*(ur[i][1]+ur[i][0])/rgm[0];
		tau3 = 0.5*(ur[i][1]+ur[i][0])/rgm[0]+0.5*(uq[i+1][1]+uq[i][1])*scm[i]/(ssm[i]*rgm[0]);
		tau4 = 0.25*AA*((uq[i+1][2]+uq[i][2]+uq[i+1][1]+uq[i][1])-2.0*(uq[i+1][0]+uq[i][0]))*byhrg/rgm[0]
			  -0.5*(uq[i+1][1]+uq[i][1])/rgm[0]
			  +0.25*((ur[i+1][1]+ur[i+1][0])-(ur[i-1][1]+ur[i-1][0]))*byhs/rgm[0];
		tau5 = AA*(rg[1]*rg[1]*ur[i][1]-rg[0]*rg[0]*ur[i][0])*byhrg/(rgm[0]*rgm[0]*rgm[0])
			  +(ss[i+1]*uq[i+1][1]-ss[i]*uq[i][1])*byhs/(rgm[0]*ssm[i]);
		qTgs[i][1] = (0.5*AA*RdotbyRs*Cpgs[i][1]*rogs[i][1]*((1.0-eta[1])*(Tgs[i][2]-Tgs[i][1])
														   +(-8.0*Tgs[i][0]+9.0*Tgs[i][1]-Tgs[i][2])/3.0)*byhrg
					//-0.25*AA*Cpgs[i][1]*rogs[i][1]*(ur[i][1]+ur[i][0])*((Tgs[i][2]+Tgs[i][1])-2.0*Tgs[i][0])*byhrg/rgm[0]
					//-0.25*Cpgs[i][1]*rogs[i][1]*(uq[i+1][1]+uq[i][1])*(Tgs[i+1][1]-Tgs[i-1][1])*byhs/rgm[0]
					-0.25*AA*Cpgs[i][1]*rogs[i][1]*byhrg/rgm[0]*
					((fabs(ur[i][1]+ur[i][0])+(ur[i][1]+ur[i][0]))*(8.0*(Tgs[i][1]-Tgs[i][0])+Tgs[i][1]-Tgs[i][2])/3.0
					-(fabs(ur[i][1]+ur[i][0])-(ur[i][1]+ur[i][0]))*(Tgs[i][2]-Tgs[i][1]))
					-0.25*Cpgs[i][1]*rogs[i][1]*byhs/rgm[0]*
					((fabs(uq[i+1][1]+uq[i][1])+(uq[i+1][1]+uq[i][1]))*(Tgs[i][1]-Tgs[i-1][1])
					-(fabs(uq[i+1][1]+uq[i][1])-(uq[i+1][1]+uq[i][1]))*(Tgs[i+1][1]-Tgs[i][1]))
				    +Asp+Bsp-CgTs[i][1]
					//+mus[i][1]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
					+0.5*AA*AA*(rg[1]*(lgs[i][2]+lgs[i][1])*(Tgs[i][2]-Tgs[i][1])
							   -2.0*rg[0]*lgs[i][0]*(-8.0*Tgs[i][0]+9.0*Tgs[i][1]-Tgs[i][2])/3.0)*byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
				    +0.5*(ss[i+1]*(lgs[i+1][1]+lgs[i][1])*(Tgs[i+1][1]-Tgs[i][1])
					     -ss[i]*(lgs[i][1]+lgs[i-1][1])*(Tgs[i][1]-Tgs[i-1][1]))*byhs*byhs/(rgm[0]*rgm[0]*ssm[i]))/(rogs[i][1]*Cpgs[i][1])*dt;
		qTgs[i][1]+= Tgs[i][1];
	}

	/*i=0*/
	for(j=2;j<Ng-1;j++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp += Dgs[0][j][sp]*Cpgsi[0][j][sp]*(Ys[0][j+1][sp]-Ys[0][j-1][sp]);
			Bsp += Dgs[0][j][sp]*Cpgsi[0][j][sp]*(Ys[1][j][sp]-Ys[0][j][sp]);
		}
		Asp *= 0.25*rogs[0][j]*(Tgs[0][j+1]-Tgs[0][j-1])*AA*AA/(rgm[j-1]*rgm[j-1])*byhrg*byhrg;
		Bsp *= 0.25*rogs[0][j]*(Tgs[1][j]-Tgs[0][j])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
		tau1 = AA*(ur[0][j]-ur[0][j-1])*byhrg/rgm[j-1];
		tau2 = (uq[1][j]-uq[0][j])*byhs/rgm[j-1]+0.5*(ur[0][j]+ur[0][j-1])/rgm[j-1];
		tau3 = 0.5*(ur[0][j]+ur[0][j-1])/rgm[j-1]+0.5*(uq[1][j]+uq[0][j])*scm[0]/(ssm[0]*rgm[j-1]);
		tau4 = 0.25*AA*((uq[1][j+1]+uq[0][j+1])-(uq[1][j-1]+uq[0][j-1]))*byhrg/rgm[j-1]
			  -0.5*(uq[1][j]+uq[0][j])/rgm[j-1]
			  +0.25*((ur[1][j]+ur[1][j-1])-(ur[0][j]+ur[0][j-1]))*byhs/rgm[j-1];
		tau5 = AA*(rg[j]*rg[j]*ur[0][j]-rg[j-1]*rg[j-1]*ur[0][j-1])*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
			  +(ss[1]*uq[1][j]-ss[0]*uq[0][j])*byhs/(rgm[j-1]*ssm[0]);
		qTgs[0][j] = (0.5*AA*RdotbyRs*Cpgs[0][j]*rogs[0][j]*((1.0-eta[j])*(Tgs[0][j+1]-Tgs[0][j])
															+(1.0-eta[j-1])*(Tgs[0][j]-Tgs[0][j-1]))*byhrg
					//-0.25*AA*Cpgs[0][j]*rogs[0][j]*(ur[0][j]+ur[0][j-1])*(Tgs[0][j+1]-Tgs[0][j-1])*byhrg/rgm[j-1]
					//-0.25*Cpgs[0][j]*rogs[0][j]*(uq[1][j]+uq[0][j])*(Tgs[1][j]-Tgs[0][j])*byhs/rgm[j-1]
					-0.25*AA*Cpgs[0][j]*rogs[0][j]*byhrg/rgm[j-1]*
					((fabs(ur[0][j]+ur[0][j-1])+(ur[0][j]+ur[0][j-1]))*(Tgs[0][j]-Tgs[0][j-1])
					-(fabs(ur[0][j]+ur[0][j-1])-(ur[0][j]+ur[0][j-1]))*(Tgs[0][j+1]-Tgs[0][j]))
					-0.25*Cpgs[0][j]*rogs[0][j]*byhs/rgm[j-1]*
					((fabs(uq[1][j]+uq[0][j])+(uq[1][j]+uq[0][j]))*(Tgs[0][j]-Tgs[0][j])
					-(fabs(uq[1][j]+uq[0][j])-(uq[1][j]+uq[0][j]))*(Tgs[1][j]-Tgs[0][j]))
				    +Asp+Bsp-CgTs[0][j]
					//+mus[0][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
					+0.5*AA*AA*(rg[j]*(lgs[0][j+1]+lgs[0][j])*(Tgs[0][j+1]-Tgs[0][j])
							   -rg[j-1]*(lgs[0][j]+lgs[0][j-1])*(Tgs[0][j]-Tgs[0][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
				    +0.5*ss[1]*(lgs[1][j]+lgs[0][j])*(Tgs[1][j]-Tgs[0][j])*byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[0]))/(rogs[0][j]*Cpgs[0][j])*dt;
		qTgs[0][j]+= Tgs[0][j];
	}
	/*i=NQ-1*/
	for(j=2;j<Ng-1;j++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp += Dgs[ii][j][sp]*Cpgsi[ii][j][sp]*(Ys[ii][j+1][sp]-Ys[ii][j-1][sp]);
			Bsp += Dgs[ii][j][sp]*Cpgsi[ii][j][sp]*(Ys[ii][j][sp]-Ys[ii-1][j][sp]);
		}
		Asp *= 0.25*rogs[ii][j]*(Tgs[ii][j+1]-Tgs[ii][j-1])*AA*AA/(rgm[j-1]*rgm[j-1])*byhrg*byhrg;
		Bsp *= 0.25*rogs[ii][j]*(Tgs[ii][j]-Tgs[ii-1][j])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
		tau1 = AA*(ur[ii][j]-ur[ii][j-1])*byhrg/rgm[j-1];
		tau2 = (uq[ii+1][j]-uq[ii][j])*byhs/rgm[j-1]+0.5*(ur[ii][j]+ur[ii][j-1])/rgm[j-1];
		tau3 = 0.5*(ur[ii][j]+ur[ii][j-1])/rgm[j-1]+0.5*(uq[ii+1][j]+uq[ii][j])*scm[ii]/(ssm[ii]*rgm[j-1]);
		tau4 = 0.25*AA*((uq[ii+1][j+1]+uq[ii][j+1])-(uq[ii+1][j-1]+uq[ii][j-1]))*byhrg/rgm[j-1]
			  -0.5*(uq[ii+1][j]+uq[ii][j])/rgm[j-1]
			  +0.25*((ur[ii][j]+ur[ii][j-1])-(ur[ii-1][j]+ur[ii-1][j-1]))*byhs/rgm[j-1];
		tau5 = AA*(rg[j]*rg[j]*ur[ii][j]-rg[j-1]*rg[j-1]*ur[ii][j-1])*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
			  +(ss[ii+1]*uq[ii+1][j]-ss[ii]*uq[ii][j])*byhs/(rgm[j-1]*ssm[ii]);
		qTgs[ii][j] = (0.5*AA*RdotbyRs*Cpgs[ii][j]*rogs[ii][j]*((1.0-eta[j])*(Tgs[ii][j+1]-Tgs[ii][j])
															   +(1.0-eta[j-1])*(Tgs[ii][j]-Tgs[ii][j-1]))*byhrg
					//-0.25*AA*Cpgs[ii][j]*rogs[ii][j]*(ur[ii][j]+ur[ii][j-1])*(Tgs[ii][j+1]-Tgs[ii][j-1])*byhrg/rgm[j-1]
					//-0.25*Cpgs[ii][j]*rogs[ii][j]*(uq[ii+1][j]+uq[ii][j])*(Tgs[ii][j]-Tgs[ii-1][j])*byhs/rgm[j-1]
					-0.25*AA*Cpgs[ii][j]*rogs[ii][j]*byhrg/rgm[j-1]*
					((fabs(ur[ii][j]+ur[ii][j-1])+(ur[ii][j]+ur[ii][j-1]))*(Tgs[ii][j]-Tgs[ii][j-1])
					-(fabs(ur[ii][j]+ur[ii][j-1])-(ur[ii][j]+ur[ii][j-1]))*(Tgs[ii][j+1]-Tgs[ii][j]))
					-0.25*Cpgs[ii][j]*rogs[ii][j]*byhs/rgm[j-1]*
					((fabs(uq[ii+1][j]+uq[ii][j])+(uq[ii+1][j]+uq[ii][j]))*(Tgs[ii][j]-Tgs[ii-1][j])
					-(fabs(uq[ii+1][j]+uq[ii][j])-(uq[ii+1][j]+uq[ii][j]))*(Tgs[ii][j]-Tgs[ii][j]))
				    +Asp+Bsp-CgTs[ii][j]
					//+mus[ii][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
					+0.5*AA*AA*(rg[j]*(lgs[ii][j+1]+lgs[ii][j])*(Tgs[ii][j+1]-Tgs[ii][j])
							   -rg[j-1]*(lgs[ii][j]+lgs[ii][j-1])*(Tgs[ii][j]-Tgs[ii][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
				    -0.5*ss[ii]*(lgs[ii][j]+lgs[ii-1][j])*(Tgs[ii][j]-Tgs[ii-1][j])*byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[ii]))/(rogs[ii][j]*Cpgs[ii][j])*dt;
		qTgs[ii][j]+= Tgs[ii][j];
	}

	for(i=1;i<NQ-1;i++)
	{
		for(j=2;j<Ng-1;j++)
		{
			Asp=Bsp=0.0;
			for(sp=0;sp<spmax;sp++)
			{
				Asp += Dgs[i][j][sp]*Cpgsi[i][j][sp]*(Ys[i][j+1][sp]-Ys[i][j-1][sp]);
				Bsp += Dgs[i][j][sp]*Cpgsi[i][j][sp]*(Ys[i+1][j][sp]-Ys[i-1][j][sp]);
			}
			Asp *= 0.25*rogs[i][j]*(Tgs[i][j+1]-Tgs[i][j-1])*AA*AA/(rgm[j-1]*rgm[j-1])*byhrg*byhrg;
			Bsp *= 0.25*rogs[i][j]*(Tgs[i+1][j]-Tgs[i-1][j])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
			tau1 = AA*(ur[i][j]-ur[i][j-1])*byhrg/rgm[j-1];
			tau2 = (uq[i+1][j]-uq[i][j])*byhs/rgm[j-1]+0.5*(ur[i][j]+ur[i][j-1])/rgm[j-1];
			tau3 = 0.5*(ur[i][j]+ur[i][j-1])/rgm[j-1]+0.5*(uq[i+1][j]+uq[i][j])*scm[i]/(ssm[i]*rgm[j-1]);
			tau4 = 0.25*AA*((uq[i+1][j+1]+uq[i][j+1])-(uq[i+1][j-1]+uq[i][j-1]))*byhrg/rgm[j-1]
				  -0.5*(uq[i+1][j]+uq[i][j])/rgm[j-1]
				  +0.25*((ur[i+1][j]+ur[i+1][j-1])-(ur[i-1][j]+ur[i-1][j-1]))*byhs/rgm[j-1];
			tau5 = AA*(rg[j]*rg[j]*ur[i][j]-rg[j-1]*rg[j-1]*ur[i][j-1])*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
				  +(ss[i+1]*uq[i+1][j]-ss[i]*uq[i][j])*byhs/(rgm[j-1]*ssm[i]);
			qTgs[i][j] = (0.5*AA*RdotbyRs*Cpgs[i][j]*rogs[i][j]*((1.0-eta[j])*(Tgs[i][j+1]-Tgs[i][j])
																+(1.0-eta[j-1])*(Tgs[i][j]-Tgs[i][j-1]))*byhrg
						//-0.25*AA*Cpgs[i][j]*rogs[i][j]*(ur[i][j]+ur[i][j-1])*(Tgs[i][j+1]-Tgs[i][j-1])*byhrg/rgm[j-1]
						//-0.25*Cpgs[i][j]*rogs[i][j]*(uq[i+1][j]+uq[i][j])*(Tgs[i+1][j]-Tgs[i-1][j])*byhs/rgm[j-1]
						-0.25*AA*Cpgs[i][j]*rogs[i][j]*byhrg/rgm[j-1]*
						((fabs(ur[i][j]+ur[i][j-1])+(ur[i][j]+ur[i][j-1]))*(Tgs[i][j]-Tgs[i][j-1])
						-(fabs(ur[i][j]+ur[i][j-1])-(ur[i][j]+ur[i][j-1]))*(Tgs[i][j+1]-Tgs[i][j]))
						-0.25*Cpgs[i][j]*rogs[i][j]*byhs/rgm[j-1]*
						((fabs(uq[i+1][j]+uq[i][j])+(uq[i+1][j]+uq[i][j]))*(Tgs[i][j]-Tgs[i-1][j])
						-(fabs(uq[i+1][j]+uq[i][j])-(uq[i+1][j]+uq[i][j]))*(Tgs[i+1][j]-Tgs[i][j]))
					    +Asp+Bsp-CgTs[i][j]
						//+mus[i][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0*3.0*tau5*tau5)
						+0.5*AA*AA*(rg[j]*(lgs[i][j+1]+lgs[i][j])*(Tgs[i][j+1]-Tgs[i][j])
								   -rg[j-1]*(lgs[i][j]+lgs[i][j-1])*(Tgs[i][j]-Tgs[i][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
					    +0.5*(ss[i+1]*(lgs[i+1][j]+lgs[i][j])*(Tgs[i+1][j]-Tgs[i][j])
						     -ss[i]*(lgs[i][j]+lgs[i-1][j])*(Tgs[i][j]-Tgs[i-1][j]))*byhs*byhs/(rgm[j-1]*rgm[j-1]*ssm[i]))/(rogs[i][j]*Cpgs[i][j])*dt;
			qTgs[i][j]+= Tgs[i][j];
		}
	}
}

void ur_cal() /* Complete */  /* ÇPéüê∏ìx */
{
	int i,j,ii=NQ-1,sp;
	double AA,RdotbyRs,byRs;
	double Asp,Bsp,Csp,Dsp,Esp;
	double tau1,tau2,tau3,tau4,tau5;

	byRs = 1.0/Rs;
	RdotbyRs = Rdot*byRs;
	AA = 1.0 / log(rmax/Rs);


	/* i=0&&j=1 */
	qur[0][1] = 0.25*AA*RdotbyRs*((1.0-etam[1])*((rogs[0][3]+rogs[0][2])*ur[0][2]-(rogs[0][2]+rogs[0][1])*ur[0][1])
								 +(1.0-etam[0])*((rogs[0][2]+rogs[0][1])*ur[0][1]-2.0*rogs[0][0]*ur[0][0]))*byhrg
			   -0.125*AA*(rgm[1]*rgm[1]*
						 ((fabs(ur[0][2]+ur[0][1])+(ur[0][2]+ur[0][1]))*(rogs[0][2]+rogs[0][1])*ur[0][1]
						 -(fabs(ur[0][2]+ur[0][1])-(ur[0][2]+ur[0][1]))*(rogs[0][1+2]+rogs[0][2])*ur[0][2])
						 -rgm[0]*rgm[0]*
						 (2.0*(fabs(ur[0][1]+ur[0][0])+(ur[0][1]+ur[0][0]))*rogs[0][0]*ur[0][0]
						 -(fabs(ur[0][1]+ur[0][0])-(ur[0][1]+ur[0][0]))*(rogs[0][2]+rogs[0][1])*ur[0][1]))
						 *byhrg/(rg[1]*rg[1]*rg[1])
			   -0.125*ss[1]*
					  ((fabs(uq[1][2]+uq[1][1])+(uq[1][2]+uq[1][1]))*(rogs[0][2]+rogs[0][1])*ur[0][1]
					  -(fabs(uq[1][2]+uq[1][1])-(uq[1][2]+uq[1][1]))*(rogs[1][2]+rogs[1][1])*ur[1][1])
					  *byhs/(rg[1]*ssm[0])
			   +0.03125*(rogs[0][2]+rogs[0][1])*(uq[1][2]+uq[1][1])*(uq[1][2]+uq[1][1])/rg[1]
			   -AA/rg[1]*(Ps[0][2]-Ps[0][1])*byhrg;
    //qur[0][1]+= AA*AA*(rgm[1]*mus[0][2]*(ur[0][2]-ur[0][1])
				//	  -rgm[0]*mus[0][1]*(ur[0][1]-ur[0][0]))*byhrg*byhrg/(rg[1]*rg[1]*rg[1])
			 //  +0.25*(mus[1][2]+mus[1][1]+mus[0][2]+mus[0][1])*ss[1]*(ur[1][1]-ur[0][1])*byhs*byhs/(rg[1]*rg[1]*ssm[0])
			 //  -0.25*(mus[1][2]+mus[1][1]+mus[0][2]+mus[0][1])*ss[1]*(uq[1][2]+uq[1][1])*byhs/(rg[1]*rg[1]*ssm[0])
			 //  -(mus[0][2]+mus[0][1])*ur[0][1]/(rg[1]*rg[1])
			 //  +0.5*AA*(0.25*((mus[1][2]+mus[1][1])-(mus[0][2]+mus[0][1]))*((uq[1][2]+uq[0][2])-(uq[1][1]+uq[0][1]))
				//	   -(mus[0][2]-mus[0][1])*((uq[1][2]+uq[1][1])-(uq[0][2]+uq[0][1])))*byhs*byhrg/(rg[1]*rg[1])
		  //     -0.25*AA*scm[0]/ssm[0]*(uq[1][2]+uq[1][1]+uq[0][2]+uq[0][1])*(mus[0][2]-mus[0][1])*byhrg/(rg[1]*rg[1])
			 //  +0.0625*(uq[1][2]+uq[1][1]+uq[0][2]+uq[0][1])*((mus[1][2]+mus[1][1])-(mus[0][2]+mus[0][1]))*byhs/(rg[1]*rg[1])
			 //  +1.0/3.0*AA*AA*(mus[0][2]*(rg[2]*rg[2]*ur[0][2]-rg[1]*rg[1]*ur[0][1])/(rgm[1]*rgm[1]*rgm[1])
				//			  -mus[0][1]*(rg[1]*rg[1]*ur[0][1]-rg[0]*rg[0]*ur[0][0])/(rgm[0]*rgm[0]*rgm[0]))*byhrg*byhrg/rg[1]
			 //  +1.0/3.0*AA*(mus[0][2]*(ss[1]*uq[1][2]-ss[0]*uq[0][2])
				//		   -mus[0][1]*(ss[1]*uq[1][1]-ss[0]*uq[0][1]))*byhrg*byhs/(rg[1]*rg[1]*ssm[0]);
    qur[0][1]*= dt;
	qur[0][1]+= 0.5*(rogs[0][2]+rogs[0][1])*ur[0][1];

	/* i=NQ-1&&j=1 */
	qur[ii][1] = 0.25*AA*RdotbyRs*((1.0-etam[1])*((rogs[ii][3]+rogs[ii][2])*ur[ii][2]-(rogs[ii][2]+rogs[ii][1])*ur[ii][1])
								  +(1.0-etam[0])*((rogs[ii][2]+rogs[ii][1])*ur[ii][1]-2.0*rogs[ii][0]*ur[ii][0]))*byhrg
			   -0.125*AA*(rgm[1]*rgm[1]*
						 ((fabs(ur[ii][2]+ur[ii][1])+(ur[ii][2]+ur[ii][1]))*(rogs[ii][2]+rogs[ii][1])*ur[ii][1]
						 -(fabs(ur[ii][2]+ur[ii][1])-(ur[ii][2]+ur[ii][1]))*(rogs[ii][1+2]+rogs[ii][2])*ur[ii][2])
						 -rgm[0]*rgm[0]*
						 (2.0*(fabs(ur[ii][1]+ur[ii][0])+(ur[ii][1]+ur[ii][0]))*rogs[ii][0]*ur[ii][0]
						 -(fabs(ur[ii][1]+ur[ii][0])-(ur[ii][1]+ur[ii][0]))*(rogs[ii][2]+rogs[ii][1])*ur[ii][1]))
						 *byhrg/(rg[1]*rg[1]*rg[1])
			   +0.125*ss[ii]*
					  ((fabs(uq[ii][2]+uq[ii][1])+(uq[ii][2]+uq[ii][1]))*(rogs[ii-1][2]+rogs[ii-1][1])*ur[ii-1][1]
					  -(fabs(uq[ii][2]+uq[ii][1])-(uq[ii][2]+uq[ii][1]))*(rogs[ii][2]+rogs[ii][1])*ur[ii][1])
					  *byhs/(rg[1]*ssm[ii])
			   +0.03125*(rogs[ii][2]+rogs[ii][1])*(uq[ii][2]+uq[ii][1])*(uq[ii][2]+uq[ii][1])/rg[1]
			   -AA/rg[1]*(Ps[ii][2]-Ps[ii][1])*byhrg;
    //qur[ii][1]+= AA*AA*(rgm[1]*mus[ii][2]*(ur[ii][2]-ur[ii][1])
				//	  -rgm[0]*mus[ii][1]*(ur[ii][1]-ur[ii][0]))*byhrg*byhrg/(rg[1]*rg[1]*rg[1])
			 //  -0.25*(mus[ii][2]+mus[ii][1]+mus[ii-1][2]+mus[ii-1][1])*ss[ii]*(ur[ii][1]-ur[ii-1][1])*byhs*byhs/(rg[1]*rg[1]*ssm[ii])
			 //  +0.25*(mus[ii][2]+mus[ii][1]+mus[ii-1][2]+mus[ii-1][1])*ss[ii]*(uq[ii][2]+uq[ii][1])*byhs/(rg[1]*rg[1]*ssm[ii])
			 //  -(mus[ii][2]+mus[ii][1])*ur[ii][1]/(rg[1]*rg[1])
			 //  +0.5*AA*(0.25*((mus[ii][2]+mus[ii][1])-(mus[ii-1][2]+mus[ii-1][1]))*((uq[ii+1][2]+uq[ii][2])-(uq[ii+1][1]+uq[ii][1]))
				//	   -(mus[ii][2]-mus[ii][1])*((uq[ii+1][2]+uq[ii+1][1])-(uq[ii][2]+uq[ii][1])))*byhs*byhrg/(rg[1]*rg[1])
		  //     -0.25*AA*scm[ii]/ssm[ii]*(uq[ii+1][2]+uq[ii+1][1]+uq[ii][2]+uq[ii][1])*(mus[ii][2]-mus[ii][1])*byhrg/(rg[1]*rg[1])
			 //  +0.0625*(uq[ii+1][2]+uq[ii+1][1]+uq[ii][2]+uq[ii][1])*((mus[ii][2]+mus[ii][1])-(mus[ii-1][2]+mus[ii-1][1]))*byhs/(rg[1]*rg[1])
			 //  +1.0/3.0*AA*AA*(mus[ii][2]*(rg[2]*rg[2]*ur[ii][2]-rg[1]*rg[1]*ur[ii][1])/(rgm[1]*rgm[1]*rgm[1])
				//			  -mus[ii][1]*(rg[1]*rg[1]*ur[ii][1]-rg[0]*rg[0]*ur[ii][0])/(rgm[0]*rgm[0]*rgm[0]))*byhrg*byhrg/rg[1]
			 //  +1.0/3.0*AA*(mus[ii][2]*(ss[ii+1]*uq[ii+1][2]-ss[ii]*uq[ii][2])
				//		   -mus[ii][1]*(ss[ii+1]*uq[ii+1][1]-ss[ii]*uq[ii][1]))*byhrg*byhs/(rg[1]*rg[1]*ssm[ii]);
    qur[ii][1]*= dt;
	qur[ii][1]+= 0.5*(rogs[ii][2]+rogs[ii][1])*ur[ii][1]; 

	/* vapor-liquid interface */
	for(i=1;i<NQ-1;i++)
	{
		qur[i][1] = 0.25*AA*RdotbyRs*((1.0-etam[1])*((rogs[i][3]+rogs[i][2])*ur[i][2]-(rogs[i][2]+rogs[i][1])*ur[i][1])
									 +(1.0-etam[0])*((rogs[i][2]+rogs[i][1])*ur[i][1]-2.0*rogs[i][0]*ur[i][0]))*byhrg
				   -0.125*AA*(rgm[1]*rgm[1]*
							 ((fabs(ur[i][2]+ur[i][1])+(ur[i][2]+ur[i][1]))*(rogs[i][2]+rogs[i][1])*ur[i][1]
							 -(fabs(ur[i][2]+ur[i][1])-(ur[i][2]+ur[i][1]))*(rogs[i][1+2]+rogs[i][2])*ur[i][2])
							 -rgm[0]*rgm[0]*
							 (2.0*(fabs(ur[i][1]+ur[i][0])+(ur[i][1]+ur[i][0]))*rogs[i][0]*ur[i][0]
							 -(fabs(ur[i][1]+ur[i][0])-(ur[i][1]+ur[i][0]))*(rogs[i][2]+rogs[i][1])*ur[i][1]))
							 *byhrg/(rg[1]*rg[1]*rg[1])
				   -0.125*(ss[i+1]*
						  ((fabs(uq[i+1][2]+uq[i+1][1])+(uq[i+1][2]+uq[i+1][1]))*(rogs[i][2]+rogs[i][1])*ur[i][1]
						  -(fabs(uq[i+1][2]+uq[i+1][1])-(uq[i+1][2]+uq[i+1][1]))*(rogs[i+1][2]+rogs[i+1][1])*ur[i+1][1])
						  -ss[i]*
						  ((fabs(uq[i][2]+uq[i][1])+(uq[i][2]+uq[i][1]))*(rogs[i-1][2]+rogs[i-1][1])*ur[i-1][1]
						  -(fabs(uq[i][2]+uq[i][1])-(uq[i][2]+uq[i][1]))*(rogs[i][2]+rogs[i][1])*ur[i][1]))
						  *byhs/(rg[1]*ssm[i])
				   +0.03125*(rogs[i][2]+rogs[i][1])*(uq[i+1][2]+uq[i+1][1]+uq[i][2]+uq[i][1])*(uq[i+1][2]+uq[i+1][1]+uq[i][2]+uq[i][1])/rg[1]
				   -AA/rg[1]*(Ps[i][2]-Ps[i][1])*byhrg;
	    //qur[i][1]+= AA*AA*(rgm[1]*mus[i][2]*(ur[i][2]-ur[i][1])
					//	  -rgm[0]*mus[i][1]*(ur[i][1]-ur[i][0]))*byhrg*byhrg/(rg[1]*rg[1]*rg[1])
				 //  +0.25*((mus[i+1][2]+mus[i+1][1]+mus[i][2]+mus[i][1])*ss[i+1]*(ur[i+1][1]-ur[i][1])
					//	 -(mus[i][2]+mus[i][1]+mus[i-1][2]+mus[i-1][1])*ss[i]*(ur[i][1]-ur[i-1][1]))*byhs*byhs/(rg[1]*rg[1]*ssm[i])
				 //  -0.25*((mus[i+1][2]+mus[i+1][1]+mus[i][2]+mus[i][1])*ss[i+1]*(uq[i+1][2]+uq[i+1][1])
					//	 -(mus[i][2]+mus[i][1]+mus[i-1][2]+mus[i-1][1])*ss[i]*(uq[i][2]+uq[i][1]))*byhs/(rg[1]*rg[1]*ssm[i])
				 //  -(mus[i][2]+mus[i][1])*ur[i][1]/(rg[1]*rg[1])
				 //  +0.5*AA*(0.25*((mus[i+1][2]+mus[i+1][1])-(mus[i-1][2]+mus[i-1][1]))*((uq[i+1][2]+uq[i][2])-(uq[i+1][1]+uq[i][1]))
					//	   -(mus[i][2]-mus[i][1])*((uq[i+1][2]+uq[i+1][1])-(uq[i][2]+uq[i][1])))*byhs*byhrg/(rg[1]*rg[1])
			  //     -0.25*AA*scm[i]/ssm[i]*(uq[i+1][2]+uq[i+1][1]+uq[i][2]+uq[i][1])*(mus[i][2]-mus[i][1])*byhrg/(rg[1]*rg[1])
				 //  +0.0625*(uq[i+1][2]+uq[i+1][1]+uq[i][2]+uq[i][1])*((mus[i+1][2]+mus[i+1][1])-(mus[i-1][2]+mus[i-1][1]))*byhs/(rg[1]*rg[1])
				 //  +1.0/3.0*AA*AA*(mus[i][2]*(rg[2]*rg[2]*ur[i][2]-rg[1]*rg[1]*ur[i][1])/(rgm[1]*rgm[1]*rgm[1])
					//			  -mus[i][1]*(rg[1]*rg[1]*ur[i][1]-rg[0]*rg[0]*ur[i][0])/(rgm[0]*rgm[0]*rgm[0]))*byhrg*byhrg/rg[1]
				 //  +1.0/3.0*AA*(mus[i][2]*(ss[i+1]*uq[i+1][2]-ss[i]*uq[i][2])
					//		   -mus[i][1]*(ss[i+1]*uq[i+1][1]-ss[i]*uq[i][1]))*byhrg*byhs/(rg[1]*rg[1]*ssm[i]);
	    qur[i][1]*= dt;
		qur[i][1]+= 0.5*(rogs[i][2]+rogs[i][1])*ur[i][1];
	}

	for(j=2;j<Ng-1;j++)
	{
		/*i=0*/
		qur[0][j] = 0.25*AA*RdotbyRs*((1.0-etam[j])*((rogs[0][j+2]+rogs[0][j+1])*ur[0][j+1]-(rogs[0][j+1]+rogs[0][j])*ur[0][j])
									 +(1.0-etam[j-1])*((rogs[0][j+1]+rogs[0][j])*ur[0][j]-(rogs[0][j]+rogs[0][j-1])*ur[0][j-1]))*byhrg
				   -0.125*AA*(rgm[j]*rgm[j]*
							 ((fabs(ur[0][j+1]+ur[0][j])+(ur[0][j+1]+ur[0][j]))*(rogs[0][j+1]+rogs[0][j])*ur[0][j]
							 -(fabs(ur[0][j+1]+ur[0][j])-(ur[0][j+1]+ur[0][j]))*(rogs[0][j+2]+rogs[0][j+1])*ur[0][j+1])
							 -rgm[j-1]*rgm[j-1]*
							 ((fabs(ur[0][j]+ur[0][j-1])+(ur[0][j]+ur[0][j-1]))*(rogs[0][j]+rogs[0][j-1])*ur[0][j-1]
							 -(fabs(ur[0][j]+ur[0][j-1])-(ur[0][j]+ur[0][j-1]))*(rogs[0][j+1]+rogs[0][j])*ur[0][j]))
							 *byhrg/(rg[j]*rg[j]*rg[j])
				   -0.125*ss[1]*
						  ((fabs(uq[1][j+1]+uq[1][j])+(uq[1][j+1]+uq[1][j]))*(rogs[0][j+1]+rogs[0][j])*ur[0][j]
						  -(fabs(uq[1][j+1]+uq[1][j])-(uq[1][j+1]+uq[1][j]))*(rogs[1][j+1]+rogs[1][j])*ur[1][j])
						  *byhs/(rg[j]*ssm[0])
				   +0.03125*(rogs[0][j+1]+rogs[0][j])*(uq[1][j+1]+uq[1][j])*(uq[1][j+1]+uq[1][j])/rg[j]
				   -AA/rg[j]*(Ps[0][j+1]-Ps[0][j])*byhrg;
	    //qur[0][j]+= AA*AA*(rgm[j]*mus[0][j+1]*(ur[0][j+1]-ur[0][j])
					//	  -rgm[j-1]*mus[0][j]*(ur[0][j]-ur[0][j-1]))*byhrg*byhrg/(rg[j]*rg[j]*rg[j])
				 //  +0.25*(mus[1][j+1]+mus[1][j]+mus[0][j+1]+mus[0][j])*ss[1]*(ur[1][j]-ur[0][j])*byhs*byhs/(rg[j]*rg[j]*ssm[0])
				 //  -0.25*(mus[1][j+1]+mus[1][j]+mus[0][j+1]+mus[0][j])*ss[1]*(uq[1][j+1]+uq[1][j])*byhs/(rg[j]*rg[j]*ssm[0])
				 //  -(mus[0][j+1]+mus[0][j])*ur[0][j]/(rg[j]*rg[j])
				 //  +0.5*AA*(0.25*((mus[1][j+1]+mus[1][j])-(mus[0][j+1]+mus[0][j]))*((uq[1][j+1]+uq[0][j+1])-(uq[1][j]+uq[0][j]))
					//	   -(mus[0][j+1]-mus[0][j])*((uq[1][j+1]+uq[1][j])-(uq[0][j+1]+uq[0][j])))*byhs*byhrg/(rg[j]*rg[j])
			  //     -0.25*AA*scm[0]/ssm[0]*(uq[1][j+1]+uq[1][j]+uq[0][j+1]+uq[0][j])*(mus[0][j+1]-mus[0][j])*byhrg/(rg[j]*rg[j])
				 //  +0.0625*(uq[1][j+1]+uq[1][j]+uq[0][j+1]+uq[0][j])*((mus[1][j+1]+mus[1][j])-(mus[0][j+1]+mus[0][j]))*byhs/(rg[j]*rg[j])
				 //  +1.0/3.0*AA*AA*(mus[0][j+1]*(rg[j+1]*rg[j+1]*ur[0][j+1]-rg[j]*rg[j]*ur[0][j])/(rgm[j]*rgm[j]*rgm[j])
					//			  -mus[0][j]*(rg[j]*rg[j]*ur[0][j]-rg[j-1]*rg[j-1]*ur[0][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1]))*byhrg*byhrg/rg[j]
				 //  +1.0/3.0*AA*(mus[0][j+1]*(ss[1]*uq[1][j+1]-ss[0]*uq[0][j+1])
					//		   -mus[0][j]*(ss[1]*uq[1][j]-ss[0]*uq[0][j]))*byhrg*byhs/(rg[j]*rg[j]*ssm[0]);
	    qur[0][j]*= dt;
		qur[0][j]+= 0.5*(rogs[0][j+1]+rogs[0][j])*ur[0][j];
		
		/*i=NQ-1*/
		qur[ii][j] = 0.25*AA*RdotbyRs*((1.0-etam[j])*((rogs[ii][j+2]+rogs[ii][j+1])*ur[ii][j+1]-(rogs[ii][j+1]+rogs[ii][j])*ur[ii][j])
									  +(1.0-etam[j-1])*((rogs[ii][j+1]+rogs[ii][j])*ur[ii][j]-(rogs[ii][j]+rogs[ii][j-1])*ur[ii][j-1]))*byhrg
				   -0.125*AA*(rgm[j]*rgm[j]*
							 ((fabs(ur[ii][j+1]+ur[ii][j])+(ur[ii][j+1]+ur[ii][j]))*(rogs[ii][j+1]+rogs[ii][j])*ur[ii][j]
							 -(fabs(ur[ii][j+1]+ur[ii][j])-(ur[ii][j+1]+ur[ii][j]))*(rogs[ii][j+2]+rogs[ii][j+1])*ur[ii][j+1])
							 -rgm[j-1]*rgm[j-1]*
							 ((fabs(ur[ii][j]+ur[ii][j-1])+(ur[ii][j]+ur[ii][j-1]))*(rogs[ii][j]+rogs[ii][j-1])*ur[ii][j-1]
							 -(fabs(ur[ii][j]+ur[ii][j-1])-(ur[ii][j]+ur[ii][j-1]))*(rogs[ii][j+1]+rogs[ii][j])*ur[ii][j]))
							 *byhrg/(rg[j]*rg[j]*rg[j])
				   +0.125*ss[ii]*
						  ((fabs(uq[ii][j+1]+uq[ii][j])+(uq[ii][j+1]+uq[ii][j]))*(rogs[ii-1][j+1]+rogs[ii-1][j])*ur[ii-1][j]
						  -(fabs(uq[ii][j+1]+uq[ii][j])-(uq[ii][j+1]+uq[ii][j]))*(rogs[ii][j+1]+rogs[ii][j])*ur[ii][j])
						  *byhs/(rg[j]*ssm[ii])
				   +0.03125*(rogs[ii][j+1]+rogs[ii][j])*(uq[ii][j+1]+uq[ii][j])*(uq[ii][j+1]+uq[ii][j])/rg[j]
				   -AA/rg[j]*(Ps[ii][j+1]-Ps[ii][j])*byhrg;
	    //qur[ii][j]+= AA*AA*(rgm[j]*mus[ii][j+1]*(ur[ii][j+1]-ur[ii][j])
					//	  -rgm[j-1]*mus[ii][j]*(ur[ii][j]-ur[ii][j-1]))*byhrg*byhrg/(rg[j]*rg[j]*rg[j])
				 //  -0.25*(mus[ii][j+1]+mus[ii][j]+mus[ii-1][j+1]+mus[ii-1][j])*ss[ii]*(ur[ii][j]-ur[ii-1][j])*byhs*byhs/(rg[j]*rg[j]*ssm[ii])
				 //  +0.25*(mus[ii][j+1]+mus[ii][j]+mus[ii-1][j+1]+mus[ii-1][j])*ss[ii]*(uq[ii][j+1]+uq[ii][j])*byhs/(rg[j]*rg[j]*ssm[ii])
				 //  -(mus[ii][j+1]+mus[ii][j])*ur[ii][j]/(rg[j]*rg[j])
				 //  +0.5*AA*(0.25*((mus[ii][j+1]+mus[ii][j])-(mus[ii-1][j+1]+mus[ii-1][j]))*((uq[ii+1][j+1]+uq[ii][j+1])-(uq[ii+1][j]+uq[ii][j]))
					//	   -(mus[ii][j+1]-mus[ii][j])*((uq[ii+1][j+1]+uq[ii+1][j])-(uq[ii][j+1]+uq[ii][j])))*byhs*byhrg/(rg[j]*rg[j])
			  //     -0.25*AA*scm[ii]/ssm[ii]*(uq[ii+1][j+1]+uq[ii+1][j]+uq[ii][j+1]+uq[ii][j])*(mus[ii][j+1]-mus[ii][j])*byhrg/(rg[j]*rg[j])
				 //  +0.0625*(uq[ii+1][j+1]+uq[ii+1][j]+uq[ii][j+1]+uq[ii][j])*((mus[ii][j+1]+mus[ii][j])-(mus[ii-1][j+1]+mus[ii-1][j]))*byhs/(rg[j]*rg[j])
				 //  +1.0/3.0*AA*AA*(mus[ii][j+1]*(rg[j+1]*rg[j+1]*ur[ii][j+1]-rg[j]*rg[j]*ur[ii][j])/(rgm[j]*rgm[j]*rgm[j])
					//			  -mus[ii][j]*(rg[j]*rg[j]*ur[ii][j]-rg[j-1]*rg[j-1]*ur[ii][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1]))*byhrg*byhrg/rg[j]
				 //  +1.0/3.0*AA*(mus[ii][j+1]*(ss[ii+1]*uq[ii+1][j+1]-ss[ii]*uq[ii][j+1])
					//		   -mus[ii][j]*(ss[ii+1]*uq[ii+1][j]-ss[ii]*uq[ii][j]))*byhrg*byhs/(rg[j]*rg[j]*ssm[ii]);
	    qur[ii][j]*= dt;
		qur[ii][j]+= 0.5*(rogs[ii][j+1]+rogs[ii][j])*ur[ii][j];
	}

	for(i=1;i<NQ-1;i++)
	{
		for(j=2;j<Ng-1;j++)
		{
			qur[i][j] = 0.25*AA*RdotbyRs*((1.0-etam[j])*((rogs[i][j+2]+rogs[i][j+1])*ur[i][j+1]-(rogs[i][j+1]+rogs[i][j])*ur[i][j])
										 +(1.0-etam[j-1])*((rogs[i][j+1]+rogs[i][j])*ur[i][j]-(rogs[i][j]+rogs[i][j-1])*ur[i][j-1]))*byhrg
					   -0.125*AA*(rgm[j]*rgm[j]*
								 ((fabs(ur[i][j+1]+ur[i][j])+(ur[i][j+1]+ur[i][j]))*(rogs[i][j+1]+rogs[i][j])*ur[i][j]
								 -(fabs(ur[i][j+1]+ur[i][j])-(ur[i][j+1]+ur[i][j]))*(rogs[i][j+2]+rogs[i][j+1])*ur[i][j+1])
								 -rgm[j-1]*rgm[j-1]*
								 ((fabs(ur[i][j]+ur[i][j-1])+(ur[i][j]+ur[i][j-1]))*(rogs[i][j]+rogs[i][j-1])*ur[i][j-1]
								 -(fabs(ur[i][j]+ur[i][j-1])-(ur[i][j]+ur[i][j-1]))*(rogs[i][j+1]+rogs[i][j])*ur[i][j]))
								 *byhrg/(rg[j]*rg[j]*rg[j])
					   -0.125*(ss[i+1]*
							  ((fabs(uq[i+1][j+1]+uq[i+1][j])+(uq[i+1][j+1]+uq[i+1][j]))*(rogs[i][j+1]+rogs[i][j])*ur[i][j]
							  -(fabs(uq[i+1][j+1]+uq[i+1][j])-(uq[i+1][j+1]+uq[i+1][j]))*(rogs[i+1][j+1]+rogs[i+1][j])*ur[i+1][j])
							  -ss[i]*
							  ((fabs(uq[i][j+1]+uq[i][j])+(uq[i][j+1]+uq[i][j]))*(rogs[i-1][j+1]+rogs[i-1][j])*ur[i-1][j]
							  -(fabs(uq[i][j+1]+uq[i][j])-(uq[i][j+1]+uq[i][j]))*(rogs[i][j+1]+rogs[i][j])*ur[i][j]))
							  *byhs/(rg[j]*ssm[i])
					   +0.03125*(rogs[i][j+1]+rogs[i][j])*(uq[i+1][j+1]+uq[i+1][j]+uq[i][j+1]+uq[i][j])
							   *(uq[i+1][j+1]+uq[i+1][j]+uq[i][j+1]+uq[i][j])/rg[j]
					   -AA/rg[j]*(Ps[i][j+1]-Ps[i][j])*byhrg;
		    //qur[i][j]+= AA*AA*(rgm[j]*mus[i][j+1]*(ur[i][j+1]-ur[i][j])
						//	  -rgm[j-1]*mus[i][j]*(ur[i][j]-ur[i][j-1]))*byhrg*byhrg/(rg[j]*rg[j]*rg[j])
					 //  +0.25*((mus[i+1][j+1]+mus[i+1][j]+mus[i][j+1]+mus[i][j])*ss[i+1]*(ur[i+1][j]-ur[i][j])
						//	 -(mus[i][j+1]+mus[i][j]+mus[i-1][j+1]+mus[i-1][j])*ss[i]*(ur[i][j]-ur[i-1][j]))*byhs*byhs/(rg[j]*rg[j]*ssm[i])
					 //  -0.25*((mus[i+1][j+1]+mus[i+1][j]+mus[i][j+1]+mus[i][j])*ss[i+1]*(uq[i+1][j+1]+uq[i+1][j])
						//	 -(mus[i][j+1]+mus[i][j]+mus[i-1][j+1]+mus[i-1][j])*ss[i]*(uq[i][j+1]+uq[i][j]))*byhs/(rg[j]*rg[j]*ssm[i])
					 //  -(mus[i][j+1]+mus[i][j])*ur[i][j]/(rg[j]*rg[j])
					 //  +0.5*AA*(0.25*((mus[i+1][j+1]+mus[i+1][j])-(mus[i-1][j+1]+mus[i-1][j]))*((uq[i+1][j+1]+uq[i][j+1])-(uq[i+1][j]+uq[i][j]))
						//	   -(mus[i][j+1]-mus[i][j])*((uq[i+1][j+1]+uq[i+1][j])-(uq[i][j+1]+uq[i][j])))*byhs*byhrg/(rg[j]*rg[j])
				  //     -0.25*AA*scm[i]/ssm[i]*(uq[i+1][j+1]+uq[i+1][j]+uq[i][j+1]+uq[i][j])*(mus[i][j+1]-mus[i][j])*byhrg/(rg[j]*rg[j])
					 //  +0.0625*(uq[i+1][j+1]+uq[i+1][j]+uq[i][j+1]+uq[i][j])*((mus[i+1][j+1]+mus[i+1][j])-(mus[i-1][j+1]+mus[i-1][j]))*byhs/(rg[j]*rg[j])
					 //  +1.0/3.0*AA*AA*(mus[i][j+1]*(rg[j+1]*rg[j+1]*ur[i][j+1]-rg[j]*rg[j]*ur[i][j])/(rgm[j]*rgm[j]*rgm[j])
						//			  -mus[i][j]*(rg[j]*rg[j]*ur[i][j]-rg[j-1]*rg[j-1]*ur[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1]))*byhrg*byhrg/rg[j]
					 //  +1.0/3.0*AA*(mus[i][j+1]*(ss[i+1]*uq[i+1][j+1]-ss[i]*uq[i][j+1])
						//		   -mus[i][j]*(ss[i+1]*uq[i+1][j]-ss[i]*uq[i][j]))*byhrg*byhs/(rg[j]*rg[j]*ssm[i]);
		    qur[i][j]*= dt;
			qur[i][j]+= 0.5*(rogs[i][j+1]+rogs[i][j])*ur[i][j];
		}
	}
}

void uq_cal()
{
	int i,j,ii=NQ-1,sp;
	double AA,RdotbyRs,byRs;
	double Asp,Bsp,Csp,Dsp,Esp;
	double tau1,tau2,tau3,tau4,tau5;

	byRs = 1.0/Rs;
	RdotbyRs = Rdot*byRs;
	AA = 1.0 / log(rmax/Rs);

	/*i==1&&j==1*/
	quq[1][1] = 0.125*AA*RdotbyRs*(1.0-etam[0])*(rogs[1][2]+rogs[1][1]+rogs[0][2]+rogs[0][1])*(uq[1][2]+uq[1][1])*byhrg
			   -0.125*AA*rg[1]*rg[1]*
						 ((fabs(ur[1][1]+ur[0][1])+(ur[1][1]+ur[0][1]))*(rogs[1][1]+rogs[0][1])*uq[1][1]
						 -(fabs(ur[1][1]+ur[0][1])-(ur[1][1]+ur[0][1]))*(rogs[1][2]+rogs[0][2])*uq[1][2])
				         *byhrg/(rgm[0]*rgm[0]*rgm[0])
			   -0.125*(ssm[1]*
					  ((fabs(uq[2][1]+uq[1][1])+(uq[2][1]+uq[1][1]))*(rogs[1][1]+rogs[0][1])*uq[1][1]
					  -(fabs(uq[2][1]+uq[1][1])-(uq[2][1]+uq[1][1]))*(rogs[2][1]+rogs[1][1])*uq[2][1])
					 -ssm[0]*
					 ((fabs(uq[1][1]+uq[0][1])+(uq[1][1]+uq[0][1]))*(rogs[0][1]+rogs[0][1])*uq[0][1]
					 -(fabs(uq[1][1]+uq[0][1])-(uq[1][1]+uq[0][1]))*(rogs[1][1]+rogs[0][1])*uq[1][1]))
					  *byhs/(rgm[0]*ss[1])
			   -0.125*(rogs[1][1]+rogs[0][1])*(ur[1][1]+ur[0][1]+ur[1][0]+ur[0][0])*uq[1][1]/rgm[0]
			   -(Ps[1][1]-Ps[0][1])*byhs/rgm[0];
	//quq[1][1]+= 0.25*AA*AA*(rg[1]*(mus[1][2]+mus[1][1]+mus[0][2]+mus[0][1])*(uq[1][2]-uq[1][1])
	//					   -rg[0]*(mus[1][0]+mus[0][0]+mus[1][0]+mus[0][0])*(8.0*(uq[1][1]-uq[1][0])+(uq[1][1]-uq[1][2]))/3.0)
	//					   *byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
	//	       +(mus[1][1]*ssm[1]*(uq[2][1]-uq[1][1])-mus[1][1]*ssm[0]*(uq[1][1]-uq[0][1]))*byhs*byhs/(rgm[0]*rgm[0]*ss[1])
	//		   -0.5*(mus[1][1]+mus[0][1])*uq[1][1]/(rgm[0]*rgm[0]*ss[1]*ss[1])
	//		   +0.5*(mus[1][1]+mus[0][1])*((ur[1][1]+ur[1][0])-(ur[0][1]+ur[0][0]))*byhs/(rgm[0]*rgm[0])
	//		   -0.25*AA*uq[1][1]*((mus[1][2]+mus[0][2]+mus[1][1]+mus[0][1])-2.0*(mus[1][0]+mus[0][0]))*byhrg/(rgm[0]*rgm[0])
	//		   +0.5*AA*(0.25*((mus[1][2]+mus[0][2]+mus[1][1]+mus[0][1])-2.0*(mus[1][0]+mus[0][0]))*((ur[1][1]+ur[1][0])-(ur[0][1]+ur[0][0]))
	//				   -(mus[1][1]-mus[0][1])*((ur[1][1]+ur[0][1])-(ur[1][0]+ur[0][0])))*byhrg*byhs/(rgm[0]*rgm[0])
	//		   -sc[1]/ss[1]*uq[1][1]*(mus[1][1]-mus[0][1])*byhs/(rgm[0]*rgm[0])
	//		   +1.0/3.0*AA*(mus[1][1]*(rg[1]*rg[1]*ur[1][1]-rg[0]*rg[0]*ur[1][0])
	//					   -mus[0][1]*(rg[1]*rg[1]*ur[0][1]-rg[0]*rg[0]*ur[0][0]))*byhs*byhrg/(rgm[0]*rgm[0]*rgm[0]*rgm[0])
	//		   +1.0/3.0*(mus[1][1]*(ss[2]*uq[2][1]-ss[1]*uq[1][1])/ssm[1]
	//					-mus[0][1]*(ss[1]*uq[1][1]-ss[0]*uq[0][1])/ssm[0])*byhs*byhs/(rgm[0]*rgm[0]);
	quq[1][1]*= dt;
	quq[1][1]+= 0.5*(rogs[1][1]+rogs[0][1])*uq[1][1];

	/*i==NQ-1&&j==1*/
	quq[ii][1] = 0.125*AA*RdotbyRs*(1.0-etam[0])*(rogs[ii][2]+rogs[ii][1]+rogs[ii-1][2]+rogs[ii-1][1])*(uq[ii][2]+uq[ii][1])*byhrg
			   -0.125*AA*rg[1]*rg[1]*
						 ((fabs(ur[ii][1]+ur[ii-1][1])+(ur[ii][1]+ur[ii-1][1]))*(rogs[ii][1]+rogs[ii-1][1])*uq[ii][1]
						 -(fabs(ur[ii][1]+ur[ii-1][1])-(ur[ii][1]+ur[ii-1][1]))*(rogs[ii][2]+rogs[ii-1][2])*uq[ii][2])
				         *byhrg/(rgm[0]*rgm[0]*rgm[0])
			   -0.125*(ssm[ii]*
					  ((fabs(uq[ii+1][1]+uq[ii][1])+(uq[ii+1][1]+uq[ii][1]))*(rogs[ii][1]+rogs[ii-1][1])*uq[ii][1]
					  -(fabs(uq[ii+1][1]+uq[ii][1])-(uq[ii+1][1]+uq[ii][1]))*(rogs[ii][1]+rogs[ii][1])*uq[ii+1][1])
					 -ssm[ii-1]*
					 ((fabs(uq[ii][1]+uq[ii-1][1])+(uq[ii][1]+uq[ii-1][1]))*(rogs[ii-1][1]+rogs[ii-2][1])*uq[ii-1][1]
					 -(fabs(uq[ii][1]+uq[ii-1][1])-(uq[ii][1]+uq[ii-1][1]))*(rogs[ii][1]+rogs[ii-1][1])*uq[ii][1]))
					  *byhs/(rgm[0]*ss[ii])
			   -0.125*(rogs[ii][1]+rogs[ii-1][1])*(ur[ii][1]+ur[ii-1][1]+ur[ii][0]+ur[ii-1][0])*uq[ii][1]/rgm[0]
			   -(Ps[ii][1]-Ps[ii-1][1])*byhs/rgm[0];
	//quq[ii][1]+= 0.25*AA*AA*(rg[1]*(mus[ii][2]+mus[ii][1]+mus[ii-1][2]+mus[ii-1][1])*(uq[ii][2]-uq[ii][1])
	//					   -rg[0]*(mus[ii][0]+mus[ii-1][0]+mus[ii][0]+mus[ii-1][0])*(8.0*(uq[ii][1]-uq[ii][0])+(uq[ii][1]-uq[ii][2]))/3.0)
	//					   *byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
	//	       +(mus[ii][1]*ssm[ii]*(uq[ii+1][1]-uq[ii][1])-mus[ii][1]*ssm[ii-1]*(uq[ii][1]-uq[ii-1][1]))*byhs*byhs/(rgm[0]*rgm[0]*ss[ii])
	//		   -0.5*(mus[ii][1]+mus[ii-1][1])*uq[ii][1]/(rgm[0]*rgm[0]*ss[ii]*ss[ii])
	//		   +0.5*(mus[ii][1]+mus[ii-1][1])*((ur[ii][1]+ur[ii][0])-(ur[ii-1][1]+ur[ii-1][0]))*byhs/(rgm[0]*rgm[0])
	//		   -0.25*AA*uq[ii][1]*((mus[ii][2]+mus[ii-1][2]+mus[ii][1]+mus[ii-1][1])-2.0*(mus[ii][0]+mus[ii-1][0]))*byhrg/(rgm[0]*rgm[0])
	//		   +0.5*AA*(0.25*((mus[ii][2]+mus[ii-1][2]+mus[ii][1]+mus[ii-1][1])-2.0*(mus[ii][0]+mus[ii-1][0]))*((ur[ii][1]+ur[ii][0])-(ur[ii-1][1]+ur[ii-1][0]))
	//				   -(mus[ii][1]-mus[ii-1][1])*((ur[ii][1]+ur[ii-1][1])-(ur[ii][0]+ur[ii-1][0])))*byhrg*byhs/(rgm[0]*rgm[0])
	//		   -sc[ii]/ss[ii]*uq[ii][1]*(mus[ii][1]-mus[ii-1][1])*byhs/(rgm[0]*rgm[0])
	//		   +1.0/3.0*AA*(mus[ii][1]*(rg[1]*rg[1]*ur[ii][1]-rg[0]*rg[0]*ur[ii][0])
	//					   -mus[ii-1][1]*(rg[1]*rg[1]*ur[ii-1][1]-rg[0]*rg[0]*ur[ii-1][0]))*byhs*byhrg/(rgm[0]*rgm[0]*rgm[0]*rgm[0])
	//		   +1.0/3.0*(mus[ii][1]*(ss[ii+1]*uq[ii+1][1]-ss[ii]*uq[ii][1])/ssm[ii]
	//					-mus[ii-1][1]*(ss[ii]*uq[ii][1]-ss[ii-1]*uq[ii-1][1])/ssm[ii-1])*byhs*byhs/(rgm[0]*rgm[0]);
	quq[ii][1]*= dt;
	quq[ii][1]+= 0.5*(rogs[ii][1]+rogs[ii-1][1])*uq[ii][1];



	/* vapor-liquid interface */
	for(i=2;i<NQ-1;i++)
	{
		quq[i][1] = 0.125*AA*RdotbyRs*(1.0-etam[0])*(rogs[i][2]+rogs[i][1]+rogs[i-1][2]+rogs[i-1][1])*(uq[i][2]+uq[i][1])*byhrg
				   -0.125*AA*rg[1]*rg[1]*
							 ((fabs(ur[i][1]+ur[i-1][1])+(ur[i][1]+ur[i-1][1]))*(rogs[i][1]+rogs[i-1][1])*uq[i][1]
							 -(fabs(ur[i][1]+ur[i-1][1])-(ur[i][1]+ur[i-1][1]))*(rogs[i][2]+rogs[i-1][2])*uq[i][2])
					         *byhrg/(rgm[0]*rgm[0]*rgm[0])
				   -0.125*(ssm[i]*
						  ((fabs(uq[i+1][1]+uq[i][1])+(uq[i+1][1]+uq[i][1]))*(rogs[i][1]+rogs[i-1][1])*uq[i][1]
						  -(fabs(uq[i+1][1]+uq[i][1])-(uq[i+1][1]+uq[i][1]))*(rogs[i+1][1]+rogs[i][1])*uq[i+1][1])
						 -ssm[i-1]*
						 ((fabs(uq[i][1]+uq[i-1][1])+(uq[i][1]+uq[i-1][1]))*(rogs[i-1][1]+rogs[i-2][1])*uq[i-1][1]
						 -(fabs(uq[i][1]+uq[i-1][1])-(uq[i][1]+uq[i-1][1]))*(rogs[i][1]+rogs[i-1][1])*uq[i][1]))
						  *byhs/(rgm[0]*ss[i])
				   -0.125*(rogs[i][1]+rogs[i-1][1])*(ur[i][1]+ur[i-1][1]+ur[i][0]+ur[i-1][0])*uq[i][1]/rgm[0]
				   -(Ps[i][1]-Ps[i-1][1])*byhs/rgm[0];
		//quq[i][1]+= 0.25*AA*AA*(rg[1]*(mus[i][2]+mus[i][1]+mus[i-1][2]+mus[i-1][1])*(uq[i][2]-uq[i][1])
		//					   -rg[0]*(mus[i][0]+mus[i-1][0]+mus[i][0]+mus[i-1][0])*(8.0*(uq[i][1]-uq[i][0])+(uq[i][1]-uq[i][2]))/3.0)
		//					   *byhrg*byhrg/(rgm[0]*rgm[0]*rgm[0])
		//	       +(mus[i][1]*ssm[i]*(uq[i+1][1]-uq[i][1])-mus[i][1]*ssm[i-1]*(uq[i][1]-uq[i-1][1]))*byhs*byhs/(rgm[0]*rgm[0]*ss[i])
		//		   -0.5*(mus[i][1]+mus[i-1][1])*uq[i][1]/(rgm[0]*rgm[0]*ss[i]*ss[i])
		//		   +0.5*(mus[i][1]+mus[i-1][1])*((ur[i][1]+ur[i][0])-(ur[i-1][1]+ur[i-1][0]))*byhs/(rgm[0]*rgm[0])
		//		   -0.25*AA*uq[i][1]*((mus[i][2]+mus[i-1][2]+mus[i][1]+mus[i-1][1])-2.0*(mus[i][0]+mus[i-1][0]))*byhrg/(rgm[0]*rgm[0])
		//		   +0.5*AA*(0.25*((mus[i][2]+mus[i-1][2]+mus[i][1]+mus[i-1][1])-2.0*(mus[i][0]+mus[i-1][0]))*((ur[i][1]+ur[i][0])-(ur[i-1][1]+ur[i-1][0]))
		//				   -(mus[i][1]-mus[i-1][1])*((ur[i][1]+ur[i-1][1])-(ur[i][0]+ur[i-1][0])))*byhrg*byhs/(rgm[0]*rgm[0])
		//		   -sc[i]/ss[i]*uq[i][1]*(mus[i][1]-mus[i-1][1])*byhs/(rgm[0]*rgm[0])
		//		   +1.0/3.0*AA*(mus[i][1]*(rg[1]*rg[1]*ur[i][1]-rg[0]*rg[0]*ur[i][0])
		//					   -mus[i-1][1]*(rg[1]*rg[1]*ur[i-1][1]-rg[0]*rg[0]*ur[i-1][0]))*byhs*byhrg/(rgm[0]*rgm[0]*rgm[0]*rgm[0])
		//		   +1.0/3.0*(mus[i][1]*(ss[i+1]*uq[i+1][1]-ss[i]*uq[i][1])/ssm[i]
		//					-mus[i-1][1]*(ss[i]*uq[i][1]-ss[i-1]*uq[i-1][1])/ssm[i-1])*byhs*byhs/(rgm[0]*rgm[0]);
		quq[i][1]*= dt;
		quq[i][1]+= 0.5*(rogs[i][1]+rogs[i-1][1])*uq[i][1];
	}

	for(j=2;j<Ng-1;j++)
	{
		/* i==1 */
		quq[1][j] = 0.25*AA*RdotbyRs*((1.0-eta[j])*((rogs[1][j+1]+rogs[0][j+1])*uq[1][j+1]-(rogs[1][j]+rogs[0][j])*uq[1][j])
									 +(1.0-eta[j-1])*((rogs[1][j]+rogs[0][j])*uq[1][j]-(rogs[1][j-1]+rogs[0][j-1])*uq[1][j-1]))*byhrg
				   -0.125*AA*(rg[j]*rg[j]*
							 ((fabs(ur[1][j]+ur[0][j])+(ur[1][j]+ur[0][j]))*(rogs[1][j]+rogs[0][j])*uq[1][j]
							 -(fabs(ur[1][j]+ur[0][j])-(ur[1][j]+ur[0][j]))*(rogs[1][j+1]+rogs[0][j+1])*uq[1][j+1])
							 -rg[j-1]*rg[j-1]*
							 ((fabs(ur[1][j-1]+ur[0][j-1])+(ur[1][j-1]+ur[0][j-1]))*(rogs[1][j-1]+rogs[0][j-1])*uq[1][j-1]
							 -(fabs(ur[1][j-1]+ur[0][j-1])-(ur[1][j-1]+ur[0][j-1]))*(rogs[1][j]+rogs[0][j])*uq[1][j]))
					         *byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
				   -0.125*(ssm[1]*
						  ((fabs(uq[2][j]+uq[1][j])+(uq[2][j]+uq[1][j]))*(rogs[1][j]+rogs[0][j])*uq[1][j]
						  -(fabs(uq[2][j]+uq[1][j])-(uq[2][j]+uq[1][j]))*(rogs[2][j]+rogs[1][j])*uq[2][j])
						 -ssm[0]*
						 ((fabs(uq[1][j]+uq[0][j])+(uq[1][j]+uq[0][j]))*(rogs[0][j]+rogs[0][j])*uq[0][j]
						 -(fabs(uq[1][j]+uq[0][j])-(uq[1][j]+uq[0][j]))*(rogs[1][j]+rogs[0][j])*uq[1][j]))
						  *byhs/(rgm[j-1]*ss[1])
				   -0.125*(rogs[1][j]+rogs[0][j])*(ur[1][j]+ur[0][j]+ur[1][j-1]+ur[0][j-1])*uq[1][j]/rgm[j-1]
				   -(Ps[1][j]-Ps[0][j])*byhs/rgm[j-1];
		//quq[1][j]+= 0.25*AA*AA*(rg[j]*(mus[1][j+1]+mus[1][j]+mus[0][j+1]+mus[0][j])*(uq[1][j+1]-uq[1][j])
		//					   -rg[j-1]*(mus[1][j]+mus[0][j]+mus[1][j-1]+mus[0][j-1])*(uq[1][j]-uq[1][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
		//	       +(mus[1][j]*ssm[1]*(uq[2][j]-uq[1][j])-mus[1][j]*ssm[0]*(uq[1][j]-uq[0][j]))*byhs*byhs/(rgm[j-1]*rgm[j-1]*ss[1])
		//		   -0.5*(mus[1][j]+mus[0][j])*uq[1][j]/(rgm[j-1]*rgm[j-1]*ss[1]*ss[1])
		//		   +0.5*(mus[1][j]+mus[0][j])*((ur[1][j]+ur[1][j-1])-(ur[0][j]+ur[0][j-1]))*byhs/(rgm[j-1]*rgm[j-1])
		//		   -0.25*AA*uq[1][j]*((mus[1][j+1]+mus[0][j+1])-(mus[1][j-1]+mus[0][j-1]))*byhrg/(rgm[j-1]*rgm[j-1])
		//		   +0.5*AA*(0.25*((mus[1][j+1]+mus[0][j+1])-(mus[1][j-1]+mus[0][j-1]))*((ur[1][j]+ur[1][j-1])-(ur[0][j]+ur[0][j-1]))
		//				   -(mus[1][j]-mus[0][j])*((ur[1][j]+ur[0][j])-(ur[1][j-1]+ur[0][j-1])))*byhrg*byhs/(rgm[j-1]*rgm[j-1])
		//		   -sc[1]/ss[1]*uq[1][j]*(mus[1][j]-mus[0][j])*byhs/(rgm[j-1]*rgm[j-1])
		//		   +1.0/3.0*AA*(mus[1][j]*(rg[j]*rg[j]*ur[1][j]-rg[j-1]*rg[j-1]*ur[1][j-1])
		//					   -mus[0][j]*(rg[j]*rg[j]*ur[0][j]-rg[j-1]*rg[j-1]*ur[0][j-1]))*byhs*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1]*rgm[j-1])
		//		   +1.0/3.0*(mus[1][j]*(ss[2]*uq[2][j]-ss[1]*uq[1][j])/ssm[1]
		//					-mus[0][j]*(ss[1]*uq[1][j]-ss[0]*uq[0][j])/ssm[0])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
		quq[1][j]*= dt;
		quq[1][j]+= 0.5*(rogs[1][j]+rogs[0][j])*uq[1][j];

		/* i==NQ-1 */
		quq[ii][j] = 0.25*AA*RdotbyRs*((1.0-eta[j])*((rogs[ii][j+1]+rogs[ii-1][j+1])*uq[ii][j+1]-(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j])
									 +(1.0-eta[j-1])*((rogs[ii][j]+rogs[ii-1][j])*uq[ii][j]-(rogs[ii][j-1]+rogs[ii-1][j-1])*uq[ii][j-1]))*byhrg
				   -0.125*AA*(rg[j]*rg[j]*
							 ((fabs(ur[ii][j]+ur[ii-1][j])+(ur[ii][j]+ur[ii-1][j]))*(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j]
							 -(fabs(ur[ii][j]+ur[ii-1][j])-(ur[ii][j]+ur[ii-1][j]))*(rogs[ii][j+1]+rogs[ii-1][j+1])*uq[ii][j+1])
							 -rg[j-1]*rg[j-1]*
							 ((fabs(ur[ii][j-1]+ur[ii-1][j-1])+(ur[ii][j-1]+ur[ii-1][j-1]))*(rogs[ii][j-1]+rogs[ii-1][j-1])*uq[ii][j-1]
							 -(fabs(ur[ii][j-1]+ur[ii-1][j-1])-(ur[ii][j-1]+ur[ii-1][j-1]))*(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j]))
					         *byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
				   -0.125*(ssm[ii]*
					  ((fabs(uq[ii+1][j]+uq[ii][j])+(uq[ii+1][j]+uq[ii][j]))*(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j]
					  -(fabs(uq[ii+1][j]+uq[ii][j])-(uq[ii+1][j]+uq[ii][j]))*(rogs[ii][j]+rogs[ii][j])*uq[ii+1][j])
					 -ssm[ii-1]*
						 ((fabs(uq[ii][j]+uq[ii-1][j])+(uq[ii][j]+uq[ii-1][j]))*(rogs[ii-1][j]+rogs[ii-2][j])*uq[ii-1][j]
						 -(fabs(uq[ii][j]+uq[ii-1][j])-(uq[ii][j]+uq[ii-1][j]))*(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j]))
						  *byhs/(rgm[j-1]*ss[ii])
				   -0.125*(rogs[ii][j]+rogs[ii-1][j])*(ur[ii][j]+ur[ii-1][j]+ur[ii][j-1]+ur[ii-1][j-1])*uq[ii][j]/rgm[j-1]
				   -(Ps[ii][j]-Ps[ii-1][j])*byhs/rgm[j-1];
		//quq[ii][j]+= 0.25*AA*AA*(rg[j]*(mus[ii][j+1]+mus[ii][j]+mus[ii-1][j+1]+mus[ii-1][j])*(uq[ii][j+1]-uq[ii][j])
		//					   -rg[j-1]*(mus[ii][j]+mus[ii-1][j]+mus[ii][j-1]+mus[ii-1][j-1])*(uq[ii][j]-uq[ii][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
		//	       +(mus[ii][j]*ssm[ii]*(uq[ii+1][j]-uq[ii][j])-mus[ii][j]*ssm[ii-1]*(uq[ii][j]-uq[ii-1][j]))*byhs*byhs/(rgm[j-1]*rgm[j-1]*ss[ii])
		//		   -0.5*(mus[ii][j]+mus[ii-1][j])*uq[ii][j]/(rgm[j-1]*rgm[j-1]*ss[ii]*ss[ii])
		//		   +0.5*(mus[ii][j]+mus[ii-1][j])*((ur[ii][j]+ur[ii][j-1])-(ur[ii-1][j]+ur[ii-1][j-1]))*byhs/(rgm[j-1]*rgm[j-1])
		//		   -0.25*AA*uq[ii][j]*((mus[ii][j+1]+mus[ii-1][j+1])-(mus[ii][j-1]+mus[ii-1][j-1]))*byhrg/(rgm[j-1]*rgm[j-1])
		//		   +0.5*AA*(0.25*((mus[ii][j+1]+mus[ii-1][j+1])-(mus[ii][j-1]+mus[ii-1][j-1]))*((ur[ii][j]+ur[ii][j-1])-(ur[ii-1][j]+ur[ii-1][j-1]))
		//				   -(mus[ii][j]-mus[ii-1][j])*((ur[ii][j]+ur[ii-1][j])-(ur[ii][j-1]+ur[ii-1][j-1])))*byhrg*byhs/(rgm[j-1]*rgm[j-1])
		//		   -sc[ii]/ss[ii]*uq[ii][j]*(mus[ii][j]-mus[ii-1][j])*byhs/(rgm[j-1]*rgm[j-1])
		//		   +1.0/3.0*AA*(mus[ii][j]*(rg[j]*rg[j]*ur[ii][j]-rg[j-1]*rg[j-1]*ur[ii][j-1])
		//					   -mus[ii-1][j]*(rg[j]*rg[j]*ur[ii-1][j]-rg[j-1]*rg[j-1]*ur[ii-1][j-1]))*byhs*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1]*rgm[j-1])
		//		   +1.0/3.0*(mus[ii][j]*(ss[ii+1]*uq[ii+1][j]-ss[ii]*uq[ii][j])/ssm[ii]
		//					-mus[ii-1][j]*(ss[ii]*uq[ii][j]-ss[ii-1]*uq[ii-1][j])/ssm[ii-1])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
		quq[ii][j]*= dt;
		quq[ii][j]+= 0.5*(rogs[ii][j]+rogs[ii-1][j])*uq[ii][j];
	}


	for(i=2;i<NQ-1;i++)
	{
		for(j=2;j<Ng-1;j++)
		{
			quq[i][j] = 0.25*AA*RdotbyRs*((1.0-eta[j])*((rogs[i][j+1]+rogs[i-1][j+1])*uq[i][j+1]-(rogs[i][j]+rogs[i-1][j])*uq[i][j])
										 +(1.0-eta[j-1])*((rogs[i][j]+rogs[i-1][j])*uq[i][j]-(rogs[i][j-1]+rogs[i-1][j-1])*uq[i][j-1]))*byhrg
					   -0.125*AA*(rg[j]*rg[j]*
								 ((fabs(ur[i][j]+ur[i-1][j])+(ur[i][j]+ur[i-1][j]))*(rogs[i][j]+rogs[i-1][j])*uq[i][j]
								 -(fabs(ur[i][j]+ur[i-1][j])-(ur[i][j]+ur[i-1][j]))*(rogs[i][j+1]+rogs[i-1][j+1])*uq[i][j+1])
								 -rg[j-1]*rg[j-1]*
								 ((fabs(ur[i][j-1]+ur[i-1][j-1])+(ur[i][j-1]+ur[i-1][j-1]))*(rogs[i][j-1]+rogs[i-1][j-1])*uq[i][j-1]
								 -(fabs(ur[i][j-1]+ur[i-1][j-1])-(ur[i][j-1]+ur[i-1][j-1]))*(rogs[i][j]+rogs[i-1][j])*uq[i][j]))
						         *byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
					   -0.125*(ssm[i]*
							  ((fabs(uq[i+1][j]+uq[i][j])+(uq[i+1][j]+uq[i][j]))*(rogs[i][j]+rogs[i-1][j])*uq[i][j]
							  -(fabs(uq[i+1][j]+uq[i][j])-(uq[i+1][j]+uq[i][j]))*(rogs[i+1][j]+rogs[i][j])*uq[i+1][j])
							 -ssm[i-1]*
							 ((fabs(uq[i][j]+uq[i-1][j])+(uq[i][j]+uq[i-1][j]))*(rogs[i-1][j]+rogs[i-2][j])*uq[i-1][j]
							 -(fabs(uq[i][j]+uq[i-1][j])-(uq[i][j]+uq[i-1][j]))*(rogs[i][j]+rogs[i-1][j])*uq[i][j]))
							  *byhs/(rgm[j-1]*ss[i])
					   -0.125*(rogs[i][j]+rogs[i-1][j])*(ur[i][j]+ur[i-1][j]+ur[i][j-1]+ur[i-1][j-1])*uq[i][j]/rgm[j-1]
					   -(Ps[i][j]-Ps[i-1][j])*byhs/rgm[j-1];
			//quq[i][j]+= 0.25*AA*AA*(rg[j]*(mus[i][j+1]+mus[i][j]+mus[i-1][j+1]+mus[i-1][j])*(uq[i][j+1]-uq[i][j])
			//					   -rg[j-1]*(mus[i][j]+mus[i-1][j]+mus[i][j-1]+mus[i-1][j-1])*(uq[i][j]-uq[i][j-1]))*byhrg*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1])
			//	       +(mus[i][j]*ssm[i]*(uq[i+1][j]-uq[i][j])-mus[i][j]*ssm[i-1]*(uq[i][j]-uq[i-1][j]))*byhs*byhs/(rgm[j-1]*rgm[j-1]*ss[i])
			//		   -0.5*(mus[i][j]+mus[i-1][j])*uq[i][j]/(rgm[j-1]*rgm[j-1]*ss[i]*ss[i])
			//		   +0.5*(mus[i][j]+mus[i-1][j])*((ur[i][j]+ur[i][j-1])-(ur[i-1][j]+ur[i-1][j-1]))*byhs/(rgm[j-1]*rgm[j-1])
			//		   -0.25*AA*uq[i][j]*((mus[i][j+1]+mus[i-1][j+1])-(mus[i][j-1]+mus[i-1][j-1]))*byhrg/(rgm[j-1]*rgm[j-1])
			//		   +0.5*AA*(0.25*((mus[i][j+1]+mus[i-1][j+1])-(mus[i][j-1]+mus[i-1][j-1]))*((ur[i][j]+ur[i][j-1])-(ur[i-1][j]+ur[i-1][j-1]))
			//				   -(mus[i][j]-mus[i-1][j])*((ur[i][j]+ur[i-1][j])-(ur[i][j-1]+ur[i-1][j-1])))*byhrg*byhs/(rgm[j-1]*rgm[j-1])
			//		   -sc[i]/ss[i]*uq[i][j]*(mus[i][j]-mus[i-1][j])*byhs/(rgm[j-1]*rgm[j-1])
			//		   +1.0/3.0*AA*(mus[i][j]*(rg[j]*rg[j]*ur[i][j]-rg[j-1]*rg[j-1]*ur[i][j-1])
			//					   -mus[i-1][j]*(rg[j]*rg[j]*ur[i-1][j]-rg[j-1]*rg[j-1]*ur[i-1][j-1]))*byhs*byhrg/(rgm[j-1]*rgm[j-1]*rgm[j-1]*rgm[j-1])
			//		   +1.0/3.0*(mus[i][j]*(ss[i+1]*uq[i+1][j]-ss[i]*uq[i][j])/ssm[i]
			//					-mus[i-1][j]*(ss[i]*uq[i][j]-ss[i-1]*uq[i-1][j])/ssm[i-1])*byhs*byhs/(rgm[j-1]*rgm[j-1]);
			quq[i][j]*= dt;
			quq[i][j]+= 0.5*(rogs[i][j]+rogs[i-1][j])*uq[i][j];
		}

	}
}

/* cylindrical coordinate system */
void Yc_cal()
{
	int i,j,sp,ii=NZ-1,jj=NX-1;

	for(sp=0;sp<spmax-1;sp++)
	{
		/*i==0&&j==0*/
		qYc[0][0][sp] = (-0.25*rogc[0][0]*byhx*
						 ((fabs(ux[0][1]+ux[0][0])+(ux[0][1]+ux[0][0]))*(Yc[0][0][sp]-Yc[0][0][sp])
						 -(fabs(ux[0][1]+ux[0][0])-(ux[0][1]+ux[0][0]))*(Yc[0][1][sp]-Yc[0][0][sp]))
					      -0.25*rogc[0][0]*byhz*
						 ((fabs(uz[0][0]+uz[1][0])+(uz[0][0]+uz[1][0]))*(8.0*(Yc[0][0][sp]-Yc[NZ][0][sp])+(Yc[0][0][sp]-Yc[1][0][sp]))/3.0
						 -(fabs(uz[0][0]+uz[1][0])-(uz[0][0]+uz[1][0]))*(Yc[1][0][sp]-Yc[0][0][sp])
																					)
						 //-0.25*rogc[0][0]*(ux[0][1]+ux[0][0])*(Yc[0][1][sp]-Yc[0][0][sp])*byhx
					  //   -0.25*rogc[0][0]*(uz[0+1][0]+uz[0][0])*(2.0*Yc[0+1][0][sp]-(Yc[0][0][sp]+Yc[0-1][0][sp]))*byhz
					     +0.25*x[1]*(rogc[0][1]+rogc[0][0])*(Dgc[0][1][sp]+Dgc[0][0][sp])*(Yc[0][1][sp]-Yc[0][0][sp])*byhx*byhx/xm[0]
					     +0.25*((rogc[0][0]+rogc[1][0])*(Dgc[0][0][sp]+Dgc[1][0][sp])*(Yc[1][0][sp]-Yc[0][0][sp])
							   -4.0*rogc[NZ][0]*Dgc[NZ][0][sp]*(9.0*Yc[0][0][sp]-8.0*Yc[NZ][0][sp]-Yc[1][0][sp])/3.0)*byhz*byhz
					     +omegac[0][0][sp])/rogc[0][0]*dt;
		qYc[0][0][sp]+= Yc[0][0][sp];

		/*i==0&&j==NX-1*/
		qYc[0][jj][sp] = (-0.25*rogc[0][jj]*byhx*
						  ((fabs(ux[0][jj+1]+ux[0][jj])+(ux[0][jj+1]+ux[0][jj]))*(Yc[0][jj][sp]-Yc[0][jj-1][sp])
						  -(fabs(ux[0][jj+1]+ux[0][jj])-(ux[0][jj+1]+ux[0][jj]))*(8.0*(Yc[0][jj+1][sp]-Yc[0][jj][sp])
																						 +(Yc[0][jj-1][sp]-Yc[0][jj][sp]))/3.0)
					       -0.25*rogc[0][jj]*byhz*
						  ((fabs(uz[1][jj]+uz[0][jj])+(uz[1][jj]+uz[0][jj]))*(8.0*(Yc[0][jj][sp]-Yc[NZ][jj][sp])
																						 +(Yc[0][jj][sp]-Yc[1][jj][sp]))/3.0
						  -(fabs(uz[1][jj]+uz[0][jj])-(uz[1][jj]+uz[0][jj]))*(Yc[1][jj][sp]-Yc[0][jj][sp]))
						 //-0.25*rogc[0][jj]*(ux[0][jj+1]+ux[0][jj])*(2.0*Yc[0][jj+1][sp]-(Yc[0][jj][sp]+Yc[0][jj-1][sp]))*byhx
					  //   -0.25*rogc[0][jj]*(uz[0+1][jj]+uz[0][jj])*(2.0*Yc[0+1][jj][sp]-(Yc[0][jj][sp]+Yc[0-1][jj][sp]))*byhz
					     +0.25*(4.0*x[jj+1]*rogc[0][jj+1]*Dgc[0][jj+1][sp]*(8.0*Yc[0][jj+1][sp]-9.0*Yc[0][jj][sp]+Yc[0][jj-1][sp])/3.0
							   -x[jj]*(rogc[0][jj]+rogc[0][jj-1])*(Dgc[0][jj][sp]+Dgc[0][jj-1][sp])*(Yc[0][jj][sp]-Yc[0][jj-1][sp]))
							   *byhx*byhx/xm[jj]
					     +0.25*((rogc[0][jj]+rogc[1][jj])*(Dgc[0][jj][sp]+Dgc[1][jj][sp])*(Yc[1][jj][sp]-Yc[0][jj][sp])
							   -4.0*rogc[NZ][jj]*Dgc[NZ][jj][sp]*(9.0*Yc[0][jj][sp]-8.0*Yc[NZ][jj][sp]-Yc[1][jj][sp])/3.0)*byhz*byhz
					     +omegac[0][jj][sp])/rogc[0][jj]*dt;
		qYc[0][jj][sp]+= Yc[0][jj][sp];

		/*i==NZ-1&&j==0*/
		qYc[ii][0][sp] = (-0.25*rogc[ii][0]*byhx*
						 ((fabs(ux[ii][1]+ux[ii][0])+(ux[ii][1]+ux[ii][0]))*(Yc[ii][0][sp]-Yc[ii][0][sp])
						 -(fabs(ux[ii][1]+ux[ii][0])-(ux[ii][1]+ux[ii][0]))*(Yc[ii][1][sp]-Yc[ii][0][sp]))
					      -0.25*rogc[ii][0]*byhz*
						 ((fabs(uz[ii+1][0]+uz[ii][0])+(uz[ii+1][0]+uz[ii][0]))*(Yc[ii][0][sp]-Yc[ii-1][0][sp])
						 -(fabs(uz[ii+1][0]+uz[ii][0])-(uz[ii+1][0]+uz[ii][0]))*(8.0*(Yc[ii+1][0][sp]-Yc[ii][0][sp])
																					+(Yc[ii-1][0][sp]-Yc[ii][0][sp]))/3.0)
						 //-0.25*rogc[ii][0]*(ux[ii][1]+ux[ii][0])*(Yc[ii][1][sp]-Yc[ii][0][sp])*byhx
					  //   -0.25*rogc[ii][0]*(uz[ii+1][0]+uz[ii][0])*(2.0*Yc[ii+1][0][sp]-(Yc[ii][0][sp]+Yc[ii-1][0][sp]))*byhz
					     +0.25*x[1]*(rogc[ii][1]+rogc[ii][0])*(Dgc[ii][1][sp]+Dgc[ii][0][sp])*(Yc[ii][1][sp]-Yc[ii][0][sp])*byhx*byhx/xm[0]
					     +0.25*(4.0*rogc[ii+1][0]*Dgc[ii+1][0][sp]*(8.0*Yc[ii+1][0][sp]-9.0*Yc[ii][0][sp]+Yc[ii-1][0][sp])/3.0
							   -(rogc[ii][0]+rogc[ii-1][0])*(Dgc[ii][0][sp]+Dgc[ii-1][0][sp])*(Yc[ii][0][sp]-Yc[ii-1][0][sp]))*byhz*byhz
					     +omegac[ii][0][sp])/rogc[ii][0]*dt;
		qYc[ii][0][sp]+= Yc[ii][0][sp];

		/*i==NZ-1&&j==NX-1*/
		qYc[ii][jj][sp] = (-0.25*rogc[ii][jj]*byhx*
						  ((fabs(ux[ii][jj+1]+ux[ii][jj])+(ux[ii][jj+1]+ux[ii][jj]))*(Yc[ii][jj][sp]-Yc[ii][jj-1][sp])
						  -(fabs(ux[ii][jj+1]+ux[ii][jj])-(ux[ii][jj+1]+ux[ii][jj]))*(8.0*(Yc[ii][jj+1][sp]-Yc[ii][jj][sp])
																						 +(Yc[ii][jj-1][sp]-Yc[ii][jj][sp]))/3.0)
					       -0.25*rogc[ii][jj]*byhz*
						  ((fabs(uz[ii+1][jj]+uz[ii][jj])+(uz[ii+1][jj]+uz[ii][jj]))*(Yc[ii][jj][sp]-Yc[ii-1][jj][sp])
						  -(fabs(uz[ii+1][jj]+uz[ii][jj])-(uz[ii+1][jj]+uz[ii][jj]))*(8.0*(Yc[ii+1][jj][sp]-Yc[ii][jj][sp])
																						 +(Yc[ii-1][jj][sp]-Yc[ii][jj][sp]))/3.0)
						 //-0.25*rogc[ii][jj]*(ux[ii][jj+1]+ux[ii][jj])*(2.0*Yc[ii][jj+1][sp]-(Yc[ii][jj][sp]+Yc[ii][jj-1][sp]))*byhx
					  //   -0.25*rogc[ii][jj]*(uz[ii+1][jj]+uz[ii][jj])*(2.0*Yc[ii+1][jj][sp]-(Yc[ii][jj][sp]+Yc[ii-1][jj][sp]))*byhz
					     +0.25*(4.0*x[jj+1]*rogc[ii][jj+1]*Dgc[ii][jj+1][sp]*(8.0*Yc[ii][jj+1][sp]-9.0*Yc[ii][jj][sp]+Yc[ii][jj-1][sp])/3.0
							   -x[jj]*(rogc[ii][jj]+rogc[ii][jj-1])*(Dgc[ii][jj][sp]+Dgc[ii][jj-1][sp])*(Yc[ii][jj][sp]-Yc[ii][jj-1][sp]))
							   *byhx*byhx/xm[jj]
					     +0.25*(4.0*rogc[ii+1][jj]*Dgc[ii+1][jj][sp]*(8.0*Yc[ii+1][jj][sp]-9.0*Yc[ii][jj][sp]+Yc[ii-1][jj][sp])/3.0
							   -(rogc[ii][jj]+rogc[ii-1][jj])*(Dgc[ii][jj][sp]+Dgc[ii-1][jj][sp])*(Yc[ii][jj][sp]-Yc[ii-1][jj][sp]))*byhz*byhz
					     +omegac[ii][jj][sp])/rogc[ii][jj]*dt;
		qYc[ii][jj][sp]+= Yc[ii][jj][sp];

		for(i=1;i<NZ-1;i++)
		{
			/*j==0*/
			if(cflg[i][0]==0)
			{
				qYc[i][0][sp] = (-0.25*rogc[i][0]*byhx*
								((fabs(ux[i][1]+ux[i][0])+(ux[i][1]+ux[i][0]))*(Yc[i][0][sp]-Yc[i][0][sp])
								-(fabs(ux[i][1]+ux[i][0])-(ux[i][1]+ux[i][0]))*(Yc[i][1][sp]-Yc[i][0][sp]))
							     -0.25*rogc[i][0]*byhz*
								((fabs(uz[i+1][0]+uz[i][0])+(uz[i+1][0]+uz[i][0]))*(Yc[i][0][sp]-Yc[i-1][0][sp])
								-(fabs(uz[i+1][0]+uz[i][0])-(uz[i+1][0]+uz[i][0]))*(Yc[i+1][0][sp]-Yc[i][0][sp]))
								 //-0.25*rogc[i][0]*(ux[i][1]+ux[i][0])*(Yc[i][1][sp]-Yc[i][0][sp])*byhx
							  //   -0.25*rogc[i][0]*(uz[i+1][0]+uz[i][0])*(Yc[i+1][0][sp]-Yc[i-1][0][sp])*byhz
							     +0.25*x[1]*(rogc[i][1]+rogc[i][0])*(Dgc[i][1][sp]+Dgc[i][0][sp])*(Yc[i][1][sp]-Yc[i][0][sp])*byhx*byhx/xm[0]
							     +0.25*((rogc[i+1][0]+rogc[i][0])*(Dgc[i+1][0][sp]+Dgc[i][0][sp])*(Yc[i+1][0][sp]-Yc[i][0][sp])
									   -(rogc[i][0]+rogc[i-1][0])*(Dgc[i][0][sp]+Dgc[i-1][0][sp])*(Yc[i][0][sp]-Yc[i-1][0][sp]))*byhz*byhz
							     +omegac[i][0][sp])/rogc[i][0]*dt;
				qYc[i][0][sp]+= Yc[i][0][sp];
			}
			/*j==NX-1*/
			qYc[i][jj][sp] = (-0.25*rogc[i][jj]*byhx*
					 		 ((fabs(ux[i][jj+1]+ux[i][jj])+(ux[i][jj+1]+ux[i][jj]))*(Yc[i][jj][sp]-Yc[i][jj-1][sp])
							 -(fabs(ux[i][jj+1]+ux[i][jj])-(ux[i][jj+1]+ux[i][jj]))*(8.0*(Yc[i][jj+1][sp]-Yc[i][jj][sp])
																						+(Yc[i][jj-1][sp]-Yc[i][jj][sp]))/3.0)
						      -0.25*rogc[i][jj]*byhz*
							 ((fabs(uz[i+1][jj]+uz[i][jj])+(uz[i+1][jj]+uz[i][jj]))*(Yc[i][jj][sp]-Yc[i-1][jj][sp])
							 -(fabs(uz[i+1][jj]+uz[i][jj])-(uz[i+1][jj]+uz[i][jj]))*(Yc[i+1][jj][sp]-Yc[i][jj][sp]))
							 //-0.25*rogc[i][jj]*(ux[i][jj+1]+ux[i][jj])*(2.0*Yc[i][jj+1][sp]-(Yc[i][jj][sp]+Yc[i][jj-1][sp]))*byhx
						  //   -0.25*rogc[i][jj]*(uz[i+1][jj]+uz[i][jj])*(Yc[i+1][jj][sp]-Yc[i-1][jj][sp])*byhz
						     +0.25*(4.0*x[jj+1]*rogc[i][jj+1]*Dgc[i][jj+1][sp]*(8.0*Yc[i][jj+1][sp]-9.0*Yc[i][jj][sp]+Yc[i][jj-1][sp])/3.0
								   -x[jj]*(rogc[i][jj]+rogc[i][jj-1])*(Dgc[i][jj][sp]+Dgc[i][jj-1][sp])*(Yc[i][jj][sp]-Yc[i][jj-1][sp]))
								   *byhx*byhx/xm[jj]
						     +0.25*((rogc[i+1][jj]+rogc[i][jj])*(Dgc[i+1][jj][sp]+Dgc[i][jj][sp])*(Yc[i+1][jj][sp]-Yc[i][jj][sp])
								   -(rogc[i][jj]+rogc[i-1][jj])*(Dgc[i][jj][sp]+Dgc[i-1][jj][sp])*(Yc[i][jj][sp]-Yc[i-1][jj][sp]))*byhz*byhz
						     +omegac[i][jj][sp])/rogc[i][jj]*dt;
			qYc[i][jj][sp]+= Yc[i][jj][sp];
		}
		for(j=1;j<NX-1;j++)
		{
			/*i==0*/
			qYc[0][j][sp] = (-0.25*rogc[0][j]*byhx*
							 ((fabs(ux[0][j+1]+ux[0][j])+(ux[0][j+1]+ux[0][j]))*(Yc[0][j][sp]-Yc[0][j-1][sp])
							 -(fabs(ux[0][j+1]+ux[0][j])-(ux[0][j+1]+ux[0][j]))*(Yc[0][j+1][sp]-Yc[0][j][sp]))
						      -0.25*rogc[0][j]*byhz*
							 ((fabs(uz[1][j]+uz[0][j])+(uz[1][j]+uz[0][j]))*(8.0*(Yc[0][j][sp]-Yc[NZ][j][sp])
																				  +(Yc[0][j][sp]-Yc[1][j][sp]))/3.0
							 -(fabs(uz[1][j]+uz[0][j])-(uz[1][j]+uz[0][j]))*(Yc[1][j][sp]-Yc[0][j][sp]))
							 //-0.25*rogc[0][j]*(ux[0][j+1]+ux[0][j])*(Yc[0][j+1][sp]-Yc[0][j-1][sp])*byhx
							 //-0.25*rogc[0][j]*(uz[0+1][j]+uz[0][j])*(2.0*Yc[0+1][j][sp]-(Yc[0][j][sp]+Yc[0-1][j][sp]))*byhz
						     +0.25*(x[j+1]*(rogc[0][j+1]+rogc[0][j])*(Dgc[0][j+1][sp]+Dgc[0][j][sp])*(Yc[0][j+1][sp]-Yc[0][j][sp])
							-x[j]*(rogc[0][j]+rogc[0][j-1])*(Dgc[0][j][sp]+Dgc[0][j-1][sp])*(Yc[0][j][sp]-Yc[0][j-1][sp]))*byhx*byhx/xm[j]
						     +0.25*((rogc[0][j]+rogc[1][j])*(Dgc[0][j][sp]+Dgc[1][j][sp])*(Yc[1][j][sp]-Yc[0][j][sp])
								   -4.0*rogc[NZ][j]*Dgc[NZ][j][sp]*(9.0*Yc[0][j][sp]-8.0*Yc[NZ][j][sp]-Yc[1][j][sp])/3.0)*byhz*byhz
						     +omegac[0][j][sp])/rogc[0][j]*dt;
			qYc[0][j][sp]+= Yc[0][j][sp];
		

			/*i==NZ-1*/
			qYc[ii][j][sp] = (-0.25*rogc[ii][j]*byhx*
							 ((fabs(ux[ii][j+1]+ux[ii][j])+(ux[ii][j+1]+ux[ii][j]))*(Yc[ii][j][sp]-Yc[ii][j-1][sp])
							 -(fabs(ux[ii][j+1]+ux[ii][j])-(ux[ii][j+1]+ux[ii][j]))*(Yc[ii][j+1][sp]-Yc[ii][j][sp]))
						      -0.25*rogc[ii][j]*byhz*
							 ((fabs(uz[ii+1][j]+uz[ii][j])+(uz[ii+1][j]+uz[ii][j]))*(Yc[ii][j][sp]-Yc[ii-1][j][sp])
							 -(fabs(uz[ii+1][j]+uz[ii][j])-(uz[ii+1][j]+uz[ii][j]))*(8.0*(Yc[ii+1][j][sp]-Yc[ii][j][sp])
																				  +(Yc[ii-1][j][sp]-Yc[ii][j][sp]))/3.0)
							 //-0.25*rogc[ii][j]*(ux[ii][j+1]+ux[ii][j])*(Yc[ii][j+1][sp]-Yc[ii][j-1][sp])*byhx
							 //-0.25*rogc[ii][j]*(uz[ii+1][j]+uz[ii][j])*(2.0*Yc[ii+1][j][sp]-(Yc[ii][j][sp]+Yc[ii-1][j][sp]))*byhz
						     +0.25*(x[j+1]*(rogc[ii][j+1]+rogc[ii][j])*(Dgc[ii][j+1][sp]+Dgc[ii][j][sp])*(Yc[ii][j+1][sp]-Yc[ii][j][sp])
							-x[j]*(rogc[ii][j]+rogc[ii][j-1])*(Dgc[ii][j][sp]+Dgc[ii][j-1][sp])*(Yc[ii][j][sp]-Yc[ii][j-1][sp]))*byhx*byhx/xm[j]
						     +0.25*(4.0*rogc[ii+1][j]*Dgc[ii+1][j][sp]*(8.0*Yc[ii+1][j][sp]-9.0*Yc[ii][j][sp]+Yc[ii-1][j][sp])/3.0
								   -(rogc[ii][j]+rogc[ii-1][j])*(Dgc[ii][j][sp]+Dgc[ii-1][j][sp])*(Yc[ii][j][sp]-Yc[ii-1][j][sp]))*byhz*byhz
						     +omegac[ii][j][sp])/rogc[ii][j]*dt;
			qYc[ii][j][sp]+= Yc[ii][j][sp];
		}

		for(i=1;i<NZ-1;i++)
		{
			for(j=1;j<NX-1;j++)
			{
				if(cflg[i][j]==0)
				{
					qYc[i][j][sp] = (-0.25*rogc[i][j]*byhx*
									((fabs(ux[i][j+1]+ux[i][j])+(ux[i][j+1]+ux[i][j]))*(Yc[i][j][sp]-Yc[i][j-1][sp])
									-(fabs(ux[i][j+1]+ux[i][j])-(ux[i][j+1]+ux[i][j]))*(Yc[i][j+1][sp]-Yc[i][j][sp]))
								     -0.25*rogc[i][j]*byhz*
									((fabs(uz[i+1][j]+uz[i][j])+(uz[i+1][j]+uz[i][j]))*(Yc[i][j][sp]-Yc[i-1][j][sp])
									-(fabs(uz[i+1][j]+uz[i][j])-(uz[i+1][j]+uz[i][j]))*(Yc[i+1][j][sp]-Yc[i][j][sp]))
				 					 //-0.25*rogc[i][j]*(ux[i][j+1]+ux[i][j])*(Yc[i][j+1][sp]-Yc[i][j-1][sp])*byhx
								   //  -0.25*rogc[i][j]*(uz[i+1][j]+uz[i][j])*(Yc[i+1][j][sp]-Yc[i-1][j][sp])*byhz
								     +0.25*(x[j+1]*(rogc[i][j+1]+rogc[i][j])*(Dgc[i][j+1][sp]+Dgc[i][j][sp])*(Yc[i][j+1][sp]-Yc[i][j][sp])
									-x[j]*(rogc[i][j]+rogc[i][j-1])*(Dgc[i][j][sp]+Dgc[i][j-1][sp])*(Yc[i][j][sp]-Yc[i][j-1][sp]))*byhx*byhx/xm[j]
								     +0.25*((rogc[i+1][j]+rogc[i][j])*(Dgc[i+1][j][sp]+Dgc[i][j][sp])*(Yc[i+1][j][sp]-Yc[i][j][sp])
										   -(rogc[i][j]+rogc[i-1][j])*(Dgc[i][j][sp]+Dgc[i-1][j][sp])*(Yc[i][j][sp]-Yc[i-1][j][sp]))*byhz*byhz
								     +omegac[i][j][sp])/rogc[i][j]*dt;
					qYc[i][j][sp]+= Yc[i][j][sp];
				}

			}
		}
	}
}
void Tgc_cal()
{
	int i,j,ii=NZ-1,jj=NX-1,sp;
	double Asp,Bsp;
	double tau1,tau2,tau3,tau4,tau5;

	/*i==0&&j==0*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp+= Dgc[0][0][sp]*Cpgci[0][0][sp]*(Yc[0][1][sp]-Yc[0][0][sp]);
		Bsp+= Dgc[0][0][sp]*Cpgci[0][0][sp]*((Yc[0][0][sp]+Yc[1][0][sp])-2.0*Yc[NZ][0][sp]);
	}
	Asp *= 0.25*rogc[0][0]*(Tgc[0][1]-Tgc[0][0])*byhx*byhx;
	Bsp *= 0.25*rogc[0][0]*((Tgc[0][0]+Tgc[1][0])-2.0*Tgc[NZ][0])*byhz*byhz;
	tau1 = (ux[0][1]-ux[0][0])*byhx;
	tau2 = 0.5*(ux[0][1]+ux[0][0])/xm[0];
	tau3 = (uz[1][0]-uz[0][0])*byhz;
	tau4 = 0.25*((uz[0][1]+uz[1][1])-(uz[1][0]+uz[0][0]))*byhx
		  +0.25*((ux[0][1]+ux[0][0]+ux[1][1]+ux[1][0])-2.0*(ux[NZ][1]+ux[NZ][0]))*byhz;
	tau5 = (x[1]*ux[0][1]-x[0]*ux[0][0])*byhx/xm[0]+(uz[0+1][0]-uz[0][0])*byhz;
	qTgc[0][0] = (-0.25*Cpgc[0][0]*rogc[0][0]*byhx*
				  ((fabs(ux[0][1]+ux[0][0])+(ux[0][1]+ux[0][0]))*(Tgc[0][0]-Tgc[0][0])
				  -(fabs(ux[0][1]+ux[0][0])-(ux[0][1]+ux[0][0]))*(Tgc[0][1]-Tgc[0][0]))
				   -0.25*Cpgc[0][0]*rogc[0][0]*byhz*
				  ((fabs(uz[1][0]+uz[0][0])+(uz[1][0]+uz[0][0]))*(Tgc[0][0]-8.0*(Tgc[NZ][0]-Tgc[0][0])-Tgc[1][0])/3.0
				  -(fabs(uz[1][0]+uz[0][0])-(uz[1][0]+uz[0][0]))*(Tgc[1][0]-Tgc[0][0]))
				   //-0.25*Cpgc[0][0]*rogc[0][0]*(ux[0][1]+ux[0][0])*(Tgc[0][1]-Tgc[0][0])*byhx
				   //-0.25*Cpgc[0][0]*rogc[0][0]*(uz[0+1][0]+uz[0][0])*(2.0*Tgc[0+1][0]-(Tgc[0][0]+Tgc[0-1][0]))*byhz
			      +Asp+Bsp-CgTc[0][0]
				  //+muc[0][0]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
				  +0.5*x[1]*(lgc[0][1]+lgc[0][0])*(Tgc[0][1]-Tgc[0][0])*byhx*byhx/xm[0]
				  +0.5*((lgc[0][0]+lgc[1][0])*(Tgc[1][0]-Tgc[0][0])
					  -2.0*lgc[NZ][0]*(9.0*Tgc[0][0]-8.0*Tgc[NZ][0]-Tgc[1][0])/3.0)*byhz*byhz)/(rogc[0][0]*Cpgc[0][0])*dt;
	qTgc[0][0]+= Tgc[0][0];
	
	/*i==NZ-1&&j==0*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp+= Dgc[ii][0][sp]*Cpgci[ii][0][sp]*(Yc[ii][1][sp]-Yc[ii][0][sp]);
		Bsp+= Dgc[ii][0][sp]*Cpgci[ii][0][sp]*(2.0*Yc[ii+1][0][sp]-(Yc[ii][0][sp]+Yc[ii-1][0][sp]));
	}
	Asp *= 0.25*rogc[ii][0]*(Tgc[ii][1]-Tgc[ii][0])*byhx*byhx;
	Bsp *= 0.25*rogc[ii][0]*(2.0*Tgc[ii+1][0]-(Tgc[ii][0]+Tgc[ii-1][0]))*byhz*byhz;
	tau1 = (ux[ii][1]-ux[ii][0])*byhx;
	tau2 = 0.5*(ux[ii][1]+ux[ii][0])/xm[0];
	tau3 = (uz[ii+1][0]-uz[ii][0])*byhz;
	tau4 = 0.25*((uz[ii+1][1]+uz[ii][1])-(uz[ii+1][0]+uz[ii][0]))*byhx
		  +0.25*(2.0*(ux[ii+1][1]+ux[ii+1][0])-(ux[ii][1]+ux[ii][0]+ux[ii-1][1]+ux[ii-1][0]))*byhz;
	tau5 = (x[1]*ux[ii][1]-x[0]*ux[ii][0])*byhx/xm[0]+(uz[ii+1][0]-uz[ii][0])*byhz;
	qTgc[ii][0] = (-0.25*Cpgc[ii][0]*rogc[ii][0]*byhx*
				  ((fabs(ux[ii][1]+ux[ii][0])+(ux[ii][1]+ux[ii][0]))*(Tgc[ii][0]-Tgc[ii][0])
				  -(fabs(ux[ii][1]+ux[ii][0])-(ux[ii][1]+ux[ii][0]))*(Tgc[ii][1]-Tgc[ii][0]))
				   -0.25*Cpgc[ii][0]*rogc[ii][0]*byhz*
				  ((fabs(uz[ii+1][0]+uz[ii][0])+(uz[ii+1][0]+uz[ii][0]))*(Tgc[ii][0]-Tgc[ii-1][0])
				  -(fabs(uz[ii+1][0]+uz[ii][0])-(uz[ii+1][0]+uz[ii][0]))*(8.0*(Tgc[ii+1][0]-Tgc[ii][0])+Tgc[ii-1][0]-Tgc[ii][0])/3.0)
				   //-0.25*Cpgc[ii][0]*rogc[ii][0]*(ux[ii][1]+ux[ii][0])*(Tgc[ii][1]-Tgc[ii][0])*byhx
				   //-0.25*Cpgc[ii][0]*rogc[ii][0]*(uz[ii+1][0]+uz[ii][0])*(2.0*Tgc[ii+1][0]-(Tgc[ii][0]+Tgc[ii-1][0]))*byhz
			      +Asp+Bsp-CgTc[ii][0]
				  //+muc[ii][0]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
				  +0.5*x[1]*(lgc[ii][1]+lgc[ii][0])*(Tgc[ii][1]-Tgc[ii][0])*byhx*byhx/xm[0]
				  +0.5*(2.0*lgc[ii+1][0]*(8.0*Tgc[ii+1][0]-9.0*Tgc[ii][0]+Tgc[ii-1][0])/3.0
					  -(lgc[ii][0]+lgc[ii-1][0])*(Tgc[ii][0]-Tgc[ii-1][0]))*byhz*byhz)/(rogc[ii][0]*Cpgc[ii][0])*dt;
	qTgc[ii][0]+= Tgc[ii][0];

	/*i==0&&j==NX-1*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp+= Dgc[0][jj][sp]*Cpgci[0][jj][sp]*(2.0*Yc[0][jj+1][sp]-(Yc[0][jj][sp]+Yc[0][jj-1][sp]));
		Bsp+= Dgc[0][jj][sp]*Cpgci[0][jj][sp]*((Yc[0][jj][sp]+Yc[1][jj][sp])-2.0*Yc[NZ][jj][sp]);
	}
	Asp *= 0.25*rogc[0][jj]*(2.0*Tgc[0][jj+1]-(Tgc[0][jj]+Tgc[0][jj-1]))*byhx*byhx;
	Bsp *= 0.25*rogc[0][jj]*((Tgc[0][jj]+Tgc[1][jj])-2.0*Tgc[NZ][jj])*byhz*byhz;
	tau1 = (ux[0][jj+1]-ux[0][jj])*byhx;
	tau2 = 0.5*(ux[0][jj+1]+ux[0][jj])/xm[jj];
	tau3 = (uz[1][jj]-uz[0][jj])*byhz;
	tau4 = 0.25*(2.0*(uz[1][jj+1]+uz[0][jj+1])-(uz[1][jj]+uz[0][jj]+uz[1][jj-1]+uz[0][jj-1]))*byhx
		  +0.25*((ux[0][jj+1]+ux[0][jj]+ux[1][jj+1]+ux[1][jj])-2.0*(ux[NZ][jj+1]+ux[NZ][jj]))*byhz;
	tau5 = (x[jj+1]*ux[0][jj+1]-x[jj]*ux[0][jj])*byhx/xm[jj]+(uz[1][jj]-uz[0][jj])*byhz;
	qTgc[0][jj] = (-0.25*Cpgc[0][jj]*rogc[0][jj]*byhx*
				   ((fabs(ux[0][jj+1]+ux[0][jj])+(ux[0][jj+1]+ux[0][jj]))*(Tgc[0][jj]-Tgc[0][jj-1])
				   -(fabs(ux[0][jj+1]+ux[0][jj])-(ux[0][jj+1]+ux[0][jj]))*(8.0*(Tgc[0][jj+1]-Tgc[0][jj])+Tgc[0][jj-1]-Tgc[0][jj])/3.0)
				    -0.25*Cpgc[0][jj]*rogc[0][jj]*byhz*
				   ((fabs(uz[1][jj]+uz[0][jj])+(uz[1][jj]+uz[0][jj]))*(Tgc[0][jj]-8.0*(Tgc[1][jj]-Tgc[0][jj])-Tgc[1][jj])/3.0
				   -(fabs(uz[1][jj]+uz[0][jj])-(uz[1][jj]+uz[0][jj]))*(Tgc[1][jj]-Tgc[0][jj]))
					//-0.25*Cpgc[0][jj]*rogc[0][jj]*(ux[0][jj+1]+ux[0][jj])*(2.0*Tgc[0][jj+1]-(Tgc[0][jj]+Tgc[0][jj-1]))*byhx
				 //   -0.25*Cpgc[0][jj]*rogc[0][jj]*(uz[0+1][jj]+uz[0][jj])*(2.0*Tgc[0+1][jj]-(Tgc[0][jj]+Tgc[0-1][jj]))*byhz
			        +Asp+Bsp-CgTc[0][jj]
					//+muc[0][jj]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
				    +0.5*(2.0*x[jj+1]*lgc[0][jj+1]*(8.0*Tgc[0][jj+1]-9.0*Tgc[0][jj]+Tgc[0][jj-1])/3.0
					     -x[jj]*(lgc[0][jj]+lgc[0][jj-1])*(Tgc[0][jj]-Tgc[0][jj-1]))*byhx*byhx/xm[jj]
				    +0.5*((lgc[1][jj]+lgc[0][jj])*(Tgc[1][jj]-Tgc[0][jj])
					    -2.0*lgc[NZ][jj]*(9.0*Tgc[0][jj]-8.0*Tgc[NZ][jj]-Tgc[1][jj])/3.0)*byhz*byhz)/(rogc[0][jj]*Cpgc[0][jj])*dt;
	qTgc[0][jj]+= Tgc[0][jj];

	/*i==NZ-1&&j==NX-1*/
	Asp=Bsp=0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Asp+= Dgc[ii][jj][sp]*Cpgci[ii][jj][sp]*(2.0*Yc[ii][jj+1][sp]-(Yc[ii][jj][sp]+Yc[ii][jj-1][sp]));
		Bsp+= Dgc[ii][jj][sp]*Cpgci[ii][jj][sp]*(2.0*Yc[ii+1][jj][sp]-(Yc[ii][jj][sp]+Yc[ii-1][jj][sp]));
	}
	Asp *= 0.25*rogc[ii][jj]*(2.0*Tgc[ii][jj+1]-(Tgc[ii][jj]+Tgc[ii][jj-1]))*byhx*byhx;
	Bsp *= 0.25*rogc[ii][jj]*(2.0*Tgc[ii+1][jj]-(Tgc[ii][jj]+Tgc[ii-1][jj]))*byhz*byhz;
	tau1 = (ux[ii][jj+1]-ux[ii][jj])*byhx;
	tau2 = 0.5*(ux[ii][jj+1]+ux[ii][jj])/xm[jj];
	tau3 = (uz[ii+1][jj]-uz[ii][jj])*byhz;
	tau4 = 0.25*(2.0*(uz[ii+1][jj+1]+uz[ii][jj+1])-(uz[ii+1][jj]+uz[ii][jj]+uz[ii+1][jj-1]+uz[ii][jj-1]))*byhx
		  +0.25*(2.0*(ux[ii+1][jj+1]+ux[ii+1][jj])-(ux[ii][jj+1]+ux[ii][jj]+ux[ii-1][jj+1]+ux[ii-1][jj]))*byhz;
	tau5 = (x[jj+1]*ux[ii][jj+1]-x[jj]*ux[ii][jj])*byhx/xm[jj]+(uz[ii+1][jj]-uz[ii][jj])*byhz;
	qTgc[ii][jj] = (-0.25*Cpgc[ii][jj]*rogc[ii][jj]*byhx*
				   ((fabs(ux[ii][jj+1]+ux[ii][jj])+(ux[ii][jj+1]+ux[ii][jj]))*(Tgc[ii][jj]-Tgc[ii][jj-1])
				   -(fabs(ux[ii][jj+1]+ux[ii][jj])-(ux[ii][jj+1]+ux[ii][jj]))*(8.0*(Tgc[ii][jj+1]-Tgc[ii][jj])+Tgc[ii][jj-1]-Tgc[ii][jj])/3.0)
				    -0.25*Cpgc[ii][jj]*rogc[ii][jj]*byhz*
				   ((fabs(uz[ii+1][jj]+uz[ii][jj])+(uz[ii+1][jj]+uz[ii][jj]))*(Tgc[ii][jj]-Tgc[ii-1][jj])
				   -(fabs(uz[ii+1][jj]+uz[ii][jj])-(uz[ii+1][jj]+uz[ii][jj]))*(8.0*(Tgc[ii+1][jj]-Tgc[ii][jj])+Tgc[ii-1][jj]-Tgc[ii][jj])/3.0)
					//-0.25*Cpgc[ii][jj]*rogc[ii][jj]*(ux[ii][jj+1]+ux[ii][jj])*(2.0*Tgc[ii][jj+1]-(Tgc[ii][jj]+Tgc[ii][jj-1]))*byhx
				 //   -0.25*Cpgc[ii][jj]*rogc[ii][jj]*(uz[ii+1][jj]+uz[ii][jj])*(2.0*Tgc[ii+1][jj]-(Tgc[ii][jj]+Tgc[ii-1][jj]))*byhz
			        +Asp+Bsp-CgTc[ii][jj]
					//+muc[ii][jj]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
				    +0.5*(2.0*x[jj+1]*lgc[ii][jj+1]*(8.0*Tgc[ii][jj+1]-9.0*Tgc[ii][jj]+Tgc[ii][jj-1])/3.0
					     -x[jj]*(lgc[ii][jj]+lgc[ii][jj-1])*(Tgc[ii][jj]-Tgc[ii][jj-1]))*byhx*byhx/xm[jj]
				    +0.5*(2.0*lgc[ii+1][jj]*(8.0*Tgc[ii+1][jj]-9.0*Tgc[ii][jj]+Tgc[ii-1][jj])/3.0
					    -(lgc[ii][jj]+lgc[ii-1][jj])*(Tgc[ii][jj]-Tgc[ii-1][jj]))*byhz*byhz)/(rogc[ii][jj]*Cpgc[ii][jj])*dt;
	qTgc[ii][jj]+= Tgc[ii][jj];

	/*i==0*/
	for(j=1;j<NX-1;j++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp+= Dgc[0][j][sp]*Cpgci[0][j][sp]*(Yc[0][j+1][sp]-Yc[0][j-1][sp]);
			Bsp+= Dgc[0][j][sp]*Cpgci[0][j][sp]*((Yc[0][j][sp]+Yc[1][j][sp])-2.0*Yc[NZ][j][sp]);
		}
		Asp *= 0.25*rogc[0][j]*(Tgc[0][j+1]-Tgc[0][j-1])*byhx*byhx;
		Bsp *= 0.25*rogc[0][j]*((Tgc[0][j]+Tgc[1][j])-2.0*Tgc[NZ][j])*byhz*byhz;
		tau1 = (ux[0][j+1]-ux[0][j])*byhx;
		tau2 = 0.5*(ux[0][j+1]+ux[0][j])/xm[j];
		tau3 = (uz[1][j]-uz[0][j])*byhz;
		tau4 = 0.25*((uz[1][j+1]+uz[0][j+1])-(uz[1][j-1]+uz[0][j-1]))*byhx
			  +0.25*((ux[0][j+1]+ux[0][j]+ux[1][j+1]+ux[1][j])-2.0*(ux[NZ][j+1]+ux[NZ][j]))*byhz;
		tau5 = (x[j+1]*ux[0][j+1]-x[j]*ux[0][j])*byhx/xm[j]+(uz[1][j]-uz[0][j])*byhz;
		qTgc[0][j] = (-0.25*Cpgc[0][j]*rogc[0][j]*byhx*
					  ((fabs(ux[0][j+1]+ux[0][j])+(ux[0][j+1]+ux[0][j]))*(Tgc[0][j]-Tgc[0][j-1])
					  -(fabs(ux[0][j+1]+ux[0][j])-(ux[0][j+1]+ux[0][j]))*(Tgc[0][j+1]-Tgc[0][j]))
					   -0.25*Cpgc[0][j]*rogc[0][j]*byhz*
					  ((fabs(uz[1][j]+uz[0][j])+(uz[1][j]+uz[0][j]))*(Tgc[0][j]-8.0*(Tgc[NZ][j]-Tgc[0][j])-Tgc[1][j])/3.0
					  -(fabs(uz[1][j]+uz[0][j])-(uz[1][j]+uz[0][j]))*(Tgc[1][j]-Tgc[0][j]))
					   //-0.25*Cpgc[0][j]*rogc[0][j]*(ux[0][j+1]+ux[0][j])*(Tgc[0][j+1]-Tgc[0][j-1])*byhx
					   //-0.25*Cpgc[0][j]*rogc[0][j]*(uz[0+1][j]+uz[0][j])*(2.0*Tgc[0+1][j]-(Tgc[0][j]+Tgc[0-1][j]))*byhz
				      +Asp+Bsp-CgTc[0][j]
					  //+muc[0][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
					  +0.5*(x[j+1]*(lgc[0][j+1]+lgc[0][j])*(Tgc[0][j+1]-Tgc[0][j])
						   -x[j]*(lgc[0][j]+lgc[0][j-1])*(Tgc[0][j]-Tgc[0][j-1]))*byhx*byhx/xm[j]
					  +0.5*((lgc[0][j]+lgc[1][j])*(Tgc[1][j]-Tgc[0][j])
						  -2.0*lgc[NZ][j]*(9.0*Tgc[0][j]-8.0*Tgc[NZ][j]-Tgc[1][j])/3.0)*byhz*byhz)/(rogc[0][j]*Cpgc[0][j])*dt;
		qTgc[0][j]+= Tgc[0][j];
	}


	/*j==0*/
	for(i=1;i<NZ-1;i++)
	{
		if(cflg[i][0]==0)
		{
			Asp=Bsp=0.0;
			for(sp=0;sp<spmax;sp++)
			{
				Asp+= Dgc[i][0][sp]*Cpgci[i][0][sp]*(Yc[i][1][sp]-Yc[i][0][sp]);
				Bsp+= Dgc[i][0][sp]*Cpgci[i][0][sp]*(Yc[i+1][0][sp]-Yc[i-1][0][sp]);
			}
			Asp *= 0.25*rogc[i][0]*(Tgc[i][1]-Tgc[i][0])*byhx*byhx;
			Bsp *= 0.25*rogc[i][0]*(Tgc[i+1][0]-Tgc[i-1][0])*byhz*byhz;
			tau1 = (ux[i][1]-ux[i][0])*byhx;
			tau2 = 0.5*(ux[i][1]+ux[i][0])/xm[0];
			tau3 = (uz[i+1][0]-uz[i][0])*byhz;
			tau4 = 0.25*((uz[i+1][1]+uz[i][1])-(uz[i+1][0]+uz[i][0]))*byhx
				  +0.25*((ux[i+1][1]+ux[i+1][0])-(ux[i-1][1]+ux[i-1][0]))*byhz;
			tau5 = (x[1]*ux[i][1]-x[0]*ux[i][0])*byhx/xm[0]+(uz[i+1][0]-uz[i][0])*byhz;
			qTgc[i][0] = (-0.25*Cpgc[i][0]*rogc[i][0]*byhx*
						 ((fabs(ux[i][1]+ux[i][0])+(ux[i][1]+ux[i][0]))*(Tgc[i][0]-Tgc[i][0])
						 -(fabs(ux[i][1]+ux[i][0])-(ux[i][1]+ux[i][0]))*(Tgc[i][1]-Tgc[i][0]))
						  -0.25*Cpgc[i][0]*rogc[i][0]*byhz*
						 ((fabs(uz[i+1][0]+uz[i][0])+(uz[i+1][0]+uz[i][0]))*(Tgc[i][0]-Tgc[i-1][0])
						 -(fabs(uz[i+1][0]+uz[i][0])-(uz[i+1][0]+uz[i][0]))*(Tgc[i+1][0]-Tgc[i][0]))
						  //-0.25*Cpgc[i][0]*rogc[i][0]*(ux[i][1]+ux[i][0])*(Tgc[i][1]-Tgc[i][0])*byhx
						  //-0.25*Cpgc[i][0]*rogc[i][0]*(uz[i+1][0]+uz[i][0])*(Tgc[i+1][0]-Tgc[i-1][0])*byhz
					      +Asp+Bsp-CgTc[i][0]
						 // +muc[i][0]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
						  +0.5*x[1]*(lgc[i][1]+lgc[i][0])*(Tgc[i][1]-Tgc[i][0])*byhx*byhx/xm[0]
						  +0.5*((lgc[i+1][0]+lgc[i][0])*(Tgc[i+1][0]-Tgc[i][0])
							   -(lgc[i][0]+lgc[i-1][0])*(Tgc[i][0]-Tgc[i-1][0]))*byhz*byhz)/(rogc[i][0]*Cpgc[i][0])*dt;
			qTgc[i][0]+= Tgc[i][0];
		}
	}
	
	/*i==NZ-1*/
	for(j=1;j<NX-1;j++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp+= Dgc[ii][j][sp]*Cpgci[ii][j][sp]*(Yc[ii][j+1][sp]-Yc[ii][j-1][sp]);
			Bsp+= Dgc[ii][j][sp]*Cpgci[ii][j][sp]*(2.0*Yc[ii+1][j][sp]-(Yc[ii][j][sp]+Yc[ii-1][j][sp]));
		}
		Asp *= 0.25*rogc[ii][j]*(Tgc[ii][j+1]-Tgc[ii][j-1])*byhx*byhx;
		Bsp *= 0.25*rogc[ii][j]*(2.0*Tgc[ii+1][j]-(Tgc[ii][j]+Tgc[ii-1][j]))*byhz*byhz;
		tau1 = (ux[ii][j+1]-ux[ii][j])*byhx;
		tau2 = 0.5*(ux[ii][j+1]+ux[ii][j])/xm[j];
		tau3 = (uz[ii+1][j]-uz[ii][j])*byhz;
		tau4 = 0.25*((uz[ii+1][j+1]+uz[ii][j+1])-(uz[ii+1][j-1]+uz[ii][j-1]))*byhx
			  +0.25*(2.0*(ux[ii+1][j+1]+ux[ii+1][j])-(ux[ii][j+1]+ux[ii][j]+ux[ii-1][j+1]+ux[ii-1][j]))*byhz;
		tau5 = (x[j+1]*ux[ii][j+1]-x[j]*ux[ii][j])*byhx/xm[j]+(uz[ii+1][j]-uz[ii][j])*byhz;
		qTgc[ii][j] = (-0.25*Cpgc[ii][j]*rogc[ii][j]*byhx*
					  ((fabs(ux[ii][j+1]+ux[ii][j])+(ux[ii][j+1]+ux[ii][j]))*(Tgc[ii][j]-Tgc[ii][j-1])
					  -(fabs(ux[ii][j+1]+ux[ii][j])-(ux[ii][j+1]+ux[ii][j]))*(Tgc[ii][j+1]-Tgc[ii][j]))
					   -0.25*Cpgc[ii][j]*rogc[ii][j]*byhz*
					  ((fabs(uz[ii+1][j]+uz[ii][j])+(uz[ii+1][j]+uz[ii][j]))*(Tgc[ii][j]-Tgc[ii-1][j])
					  -(fabs(uz[ii+1][j]+uz[ii][j])-(uz[ii+1][j]+uz[ii][j]))*(8.0*(Tgc[ii+1][j]-Tgc[ii][j])+Tgc[ii-1][j]-Tgc[ii][j])/3.0)
					   //-0.25*Cpgc[ii][j]*rogc[ii][j]*(ux[ii][j+1]+ux[ii][j])*(Tgc[ii][j+1]-Tgc[ii][j-1])*byhx
					   //-0.25*Cpgc[ii][j]*rogc[ii][j]*(uz[ii+1][j]+uz[ii][j])*(2.0*Tgc[ii+1][j]-(Tgc[ii][j]+Tgc[ii-1][j]))*byhz
				      +Asp+Bsp-CgTc[ii][j]
					  //+muc[ii][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
					  +0.5*(x[j+1]*(lgc[ii][j+1]+lgc[ii][j])*(Tgc[ii][j+1]-Tgc[ii][j])
						   -x[j]*(lgc[ii][j]+lgc[ii][j-1])*(Tgc[ii][j]-Tgc[ii][j-1]))*byhx*byhx/xm[j]
					  +0.5*(2.0*lgc[ii+1][j]*(8.0*Tgc[ii+1][j]-9.0*Tgc[ii][j]+Tgc[ii-1][j])/3.0
						  -(lgc[ii][j]+lgc[ii-1][j])*(Tgc[ii][j]-Tgc[ii-1][j]))*byhz*byhz)/(rogc[ii][j]*Cpgc[ii][j])*dt;
		qTgc[ii][j]+= Tgc[ii][j];
	}

	/*j==NX-1*/
	for(i=1;i<NZ-1;i++)
	{
		Asp=Bsp=0.0;
		for(sp=0;sp<spmax;sp++)
		{
			Asp+= Dgc[i][jj][sp]*Cpgci[i][jj][sp]*(2.0*Yc[i][jj+1][sp]-(Yc[i][jj][sp]+Yc[i][jj-1][sp]));
			Bsp+= Dgc[i][jj][sp]*Cpgci[i][jj][sp]*(Yc[i+1][jj][sp]-Yc[i-1][jj][sp]);
		}
		Asp *= 0.25*rogc[i][jj]*(2.0*Tgc[i][jj+1]-(Tgc[i][jj]+Tgc[i][jj-1]))*byhx*byhx;
		Bsp *= 0.25*rogc[i][jj]*(Tgc[i+1][jj]-Tgc[i-1][jj])*byhz*byhz;
		tau1 = (ux[i][jj+1]-ux[i][jj])*byhx;
		tau2 = 0.5*(ux[i][jj+1]+ux[i][jj])/xm[jj];
		tau3 = (uz[i+1][jj]-uz[i][jj])*byhz;
		tau4 = 0.25*(2.0*(uz[i+1][jj+1]+uz[i][jj+1])-(uz[i+1][jj]+uz[i][jj]+uz[i+1][jj-1]+uz[i][jj-1]))*byhx
			  +0.25*((ux[i+1][jj+1]+ux[i+1][jj])-(ux[i-1][jj+1]+ux[i-1][jj]))*byhz;
		tau5 = (x[jj+1]*ux[i][jj+1]-x[jj]*ux[i][jj])*byhx/xm[jj]+(uz[i+1][jj]-uz[i][jj])*byhz;
		qTgc[i][jj] = (-0.25*Cpgc[i][jj]*rogc[i][jj]*byhx*
					  ((fabs(ux[i][jj+1]+ux[i][jj])+(ux[i][jj+1]+ux[i][jj]))*(Tgc[i][jj]-Tgc[i][jj-1])
					  -(fabs(ux[i][jj+1]+ux[i][jj])-(ux[i][jj+1]+ux[i][jj]))*(8.0*(Tgc[i][jj+1]-Tgc[i][jj])+Tgc[i][jj-1]-Tgc[i][jj])/3.0)
					   -0.25*Cpgc[i][jj]*rogc[i][jj]*byhz*
					  ((fabs(uz[i+1][jj]+uz[i][jj])+(uz[i+1][jj]+uz[i][jj]))*(Tgc[i][jj]-Tgc[i-1][jj])
					  -(fabs(uz[i+1][jj]+uz[i][jj])-(uz[i+1][jj]+uz[i][jj]))*(Tgc[i+1][jj]-Tgc[i][jj]))
					   //-0.25*Cpgc[i][jj]*rogc[i][jj]*(ux[i][jj+1]+ux[i][jj])*(2.0*Tgc[i][jj+1]-(Tgc[i][jj]+Tgc[i][jj-1]))*byhx
					   //-0.25*Cpgc[i][jj]*rogc[i][jj]*(uz[i+1][jj]+uz[i][jj])*(Tgc[i+1][jj]-Tgc[i-1][jj])*byhz
				       +Asp+Bsp-CgTc[i][jj]
					  // +muc[i][jj]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
					   +0.5*(2.0*x[jj+1]*lgc[i][jj+1]*(8.0*Tgc[i][jj+1]-9.0*Tgc[i][jj]+Tgc[i][jj-1])/3.0
						    -x[jj]*(lgc[i][jj]+lgc[i][jj-1])*(Tgc[i][jj]-Tgc[i][jj-1]))*byhx*byhx/xm[jj]
					   +0.5*((lgc[i+1][jj]+lgc[i][jj])*(Tgc[i+1][jj]-Tgc[i][jj])
						    -(lgc[i][jj]+lgc[i-1][jj])*(Tgc[i][jj]-Tgc[i-1][jj]))*byhz*byhz)/(rogc[i][jj]*Cpgc[i][jj])*dt;
		qTgc[i][jj]+= Tgc[i][jj];
	}

	for(i=1;i<NZ-1;i++)
	{
		for(j=1;j<NX-1;j++)
		{
			if(cflg[i][j]==0)
			{
				Asp=Bsp=0.0;
				for(sp=0;sp<spmax;sp++)
				{
					Asp+= Dgc[i][j][sp]*Cpgci[i][j][sp]*(Yc[i][j+1][sp]-Yc[i][j-1][sp]);
					Bsp+= Dgc[i][j][sp]*Cpgci[i][j][sp]*(Yc[i+1][j][sp]-Yc[i-1][j][sp]);
				}
				Asp *= 0.25*rogc[i][j]*(Tgc[i][j+1]-Tgc[i][j-1])*byhx*byhx;
				Bsp *= 0.25*rogc[i][j]*(Tgc[i+1][j]-Tgc[i-1][j])*byhz*byhz;
				tau1 = (ux[i][j+1]-ux[i][j])*byhx;
				tau2 = 0.5*(ux[i][j+1]+ux[i][j])/xm[j];
				tau3 = (uz[i+1][j]-uz[i][j])*byhz;
				tau4 = 0.25*((uz[i+1][j+1]+uz[i][j+1])-(uz[i+1][j-1]+uz[i][j-1]))*byhx
					  +0.25*((ux[i+1][j+1]+ux[i+1][j])-(ux[i-1][j+1]+ux[i-1][j]))*byhz;
				tau5 = (x[j+1]*ux[i][j+1]-x[j]*ux[i][j])*byhx/xm[j]+(uz[i+1][j]-uz[i][j])*byhz;
				qTgc[i][j] = (-0.25*Cpgc[i][j]*rogc[i][j]*byhx*
							 ((fabs(ux[i][j+1]+ux[i][j])+(ux[i][j+1]+ux[i][j]))*(Tgc[i][j]-Tgc[i][j-1])
							 -(fabs(ux[i][j+1]+ux[i][j])-(ux[i][j+1]+ux[i][j]))*(Tgc[i][j+1]-Tgc[i][j]))
							  -0.25*Cpgc[i][j]*rogc[i][j]*byhz*
							 ((fabs(uz[i+1][j]+uz[i][j])+(uz[i+1][j]+uz[i][j]))*(Tgc[i][j]-Tgc[i-1][j])
							 -(fabs(uz[i+1][j]+uz[i][j])-(uz[i+1][j]+uz[i][j]))*(Tgc[i+1][j]-Tgc[i][j]))
							  //-0.25*Cpgc[i][j]*rogc[i][j]*(ux[i][j+1]+ux[i][j])*(Tgc[i][j+1]-Tgc[i][j-1])*byhx
							  //-0.25*Cpgc[i][j]*rogc[i][j]*(uz[i+1][j]+uz[i][j])*(Tgc[i+1][j]-Tgc[i-1][j])*byhz
						      +Asp+Bsp-CgTc[i][j]
							 // +muc[i][j]*(2.0*(tau1*tau1+tau2*tau2+tau3*tau3)+tau4*tau4-2.0/3.0*tau5*tau5)
							  +0.5*(x[j+1]*(lgc[i][j+1]+lgc[i][j])*(Tgc[i][j+1]-Tgc[i][j])
								   -x[j]*(lgc[i][j]+lgc[i][j-1])*(Tgc[i][j]-Tgc[i][j-1]))*byhx*byhx/xm[j]
							  +0.5*((lgc[i+1][j]+lgc[i][j])*(Tgc[i+1][j]-Tgc[i][j])
								   -(lgc[i][j]+lgc[i-1][j])*(Tgc[i][j]-Tgc[i-1][j]))*byhz*byhz)/(rogc[i][j]*Cpgc[i][j])*dt;
				qTgc[i][j]+= Tgc[i][j];
			}
		}
	}
}

void ux_cal()
{
	int i,j,ii=NZ-1,jj=NX-1;

	/* i==0&&j==1 */
	qux[0][1] = -0.125*(xm[1]*
					    ((fabs(ux[0][2]+ux[0][1])+(ux[0][2]+ux[0][1]))*(rogc[0][1]+rogc[0][0])*ux[0][1]
					    -(fabs(ux[0][2]+ux[0][1])-(ux[0][2]+ux[0][1]))*(rogc[0][2]+rogc[0][1])*ux[0][2])
					    +xm[0]*(fabs(ux[0][1]+ux[0][0])-(ux[0][1]+ux[0][0]))*(rogc[0][1]+rogc[0][0])*ux[0][1])
					    *byhx/x[1]
				  -0.125*((fabs(uz[1][1]+uz[1][0])+(uz[1][1]+uz[1][0]))*(rogc[0][1]+rogc[0][0])*ux[0][1]
						 -(fabs(uz[1][1]+uz[1][0])-(uz[1][1]+uz[1][0]))*(rogc[1][1]+rogc[1][0])*ux[1][1])
						 *byhz
				 -(Pc[0][1]-Pc[0][0])*byhx;
	//qux[0][1]+= (xm[1]*muc[0][1]*(ux[0][2]-ux[0][1])
	//			 -xm[0]*muc[0][0]*(ux[0][1]-ux[0][0]))*byhx*byhx/x[1]
	//		     +0.25*((muc[0+1][1]+muc[0+1][0]+muc[0+1][1]+muc[0+1][0])*(8.0*(ux[0+1][1]-ux[0][1])+(ux[0-1][1]-ux[0][1]))/3.0
	//				   -(muc[0][1]+muc[0][0]+muc[0-1][1]+muc[0-1][0])*(ux[0][1]-ux[0-1][1]))*byhz*byhz
	//		     -0.5*(muc[0][1]+muc[0][0])*ux[0][1]/(x[1]*x[1])
	//		     -ux[0][1]*(muc[0][1]-muc[0][0])*byhx/x[1]
	//		     +0.125*((2.0*(muc[0+1][1]+muc[0+1][0])-(muc[0][1]+muc[0][0]+muc[0-1][1]+muc[0-1][0]))
	//				   *((uz[0+1][1]+uz[0][1])-(uz[0+1][0]+uz[0][0])))*byhz*byhx
	//		     -0.5*(muc[0][1]-muc[0][0])*((uz[0+1][1]+uz[0+1][0])-(uz[0][1]+uz[0][0]))*byhz*byhx
	//		     +1.0/3.0*(muc[0][1]*(x[2]*ux[0][2]-x[1]*ux[0][1])/xm[1]
	//					  -muc[0][0]*(x[1]*ux[0][1]-x[0]*ux[0][0])/xm[0])*byhx*byhx
	//		     +1.0/3.0*(muc[0][1]*(uz[0+1][1]-uz[0][1])-muc[0][0]*(uz[0+1][0]-uz[0][0]))*byhx*byhz;
	qux[0][1]*= dt;
	qux[0][1]+= 0.5*(rogc[0][1]+rogc[0][0])*ux[0][1];

	/* i==0&&j==NX-1 */
   qux[0][jj] = -0.125*(xm[jj]*(fabs(ux[0][jj+1]+ux[0][jj])+(ux[0][jj+1]+ux[0][jj]))*(rogc[0][jj]+rogc[0][jj-1])*ux[0][jj]
					    -xm[jj-1]*
					    ((fabs(ux[0][jj]+ux[0][jj-1])+(ux[0][jj]+ux[0][jj-1]))*(rogc[0][jj-1]+rogc[0][jj-2])*ux[0][jj-1]
					    -(fabs(ux[0][jj]+ux[0][jj-1])-(ux[0][jj]+ux[0][jj-1]))*(rogc[0][jj]+rogc[0][jj-1])*ux[0][jj]))
					    *byhx/x[jj]
				 +0.125*((fabs(uz[1][jj]+uz[1][jj-1])+(uz[1][jj]+uz[1][jj-1]))*(rogc[0][jj]+rogc[0][jj-1])*ux[0][jj]
						-(fabs(uz[1][jj]+uz[1][jj-1])-(uz[1][jj]+uz[1][jj-1]))*(rogc[1][jj]+rogc[1][jj-1])*ux[1][jj])
						*byhz
				 -(Pc[0][jj]-Pc[0][jj-1])*byhx;
	//qux[0][jj]+= (xm[jj]*muc[0][jj]*(ux[0][jj+1]-ux[0][jj])
	//			 -xm[jj-1]*muc[0][jj-1]*(ux[0][jj]-ux[0][jj-1]))*byhx*byhx/x[jj]
	//		     +0.25*((muc[0+1][jj]+muc[0+1][jj-1]+muc[0+1][jj]+muc[0+1][jj-1])*(8.0*(ux[0+1][jj]-ux[0][jj])+(ux[0-1][jj]-ux[0][jj]))/3.0
	//				   -(muc[0][jj]+muc[0][jj-1]+muc[0-1][jj]+muc[0-1][jj-1])*(ux[0][jj]-ux[0-1][jj]))*byhz*byhz
	//		     -0.5*(muc[0][jj]+muc[0][jj-1])*ux[0][jj]/(x[jj]*x[jj])
	//		     -ux[0][jj]*(muc[0][jj]-muc[0][jj-1])*byhx/x[jj]
	//		     +0.125*((2.0*(muc[0+1][jj]+muc[0+1][jj-1])-(muc[0][jj]+muc[0][jj-1]+muc[0-1][jj]+muc[0-1][jj-1]))
	//				   *((uz[0+1][jj]+uz[0][jj])-(uz[0+1][jj-1]+uz[0][jj-1])))*byhz*byhx
	//		     -0.5*(muc[0][jj]-muc[0][jj-1])*((uz[0+1][jj]+uz[0+1][jj-1])-(uz[0][jj]+uz[0][jj-1]))*byhz*byhx
	//		     +1.0/3.0*(muc[0][jj]*(x[jj+1]*ux[0][jj+1]-x[jj]*ux[0][jj])/xm[jj]
	//					  -muc[0][jj-1]*(x[jj]*ux[0][jj]-x[jj-1]*ux[0][jj-1])/xm[jj-1])*byhx*byhx
	//		     +1.0/3.0*(muc[0][jj]*(uz[0+1][jj]-uz[0][jj])-muc[0][jj-1]*(uz[0+1][jj-1]-uz[0][jj-1]))*byhx*byhz;
	qux[0][jj]*= dt;
	qux[0][jj]+= 0.5*(rogc[0][jj]+rogc[0][jj-1])*ux[0][jj];

	/* i==NZ-1&&j==1 */
	qux[ii][1] = -0.125*(xm[1]*
					    ((fabs(ux[ii][2]+ux[ii][1])+(ux[ii][2]+ux[ii][1]))*(rogc[ii][1]+rogc[ii][0])*ux[ii][1]
					    -(fabs(ux[ii][2]+ux[ii][1])-(ux[ii][2]+ux[ii][1]))*(rogc[ii][2]+rogc[ii][1])*ux[ii][2])
					    +xm[0]*(fabs(ux[ii][1]+ux[ii][0])-(ux[ii][1]+ux[ii][0]))*(rogc[ii][1]+rogc[ii][0])*ux[ii][1])
					    *byhx/x[1]
				  +0.125*((fabs(uz[ii][1]+uz[ii][0])+(uz[ii][1]+uz[ii][0]))*(rogc[ii-1][1]+rogc[ii-1][0])*ux[ii-1][1]
						 -(fabs(uz[ii][1]+uz[ii][0])-(uz[ii][1]+uz[ii][0]))*(rogc[ii][1]+rogc[ii][0])*ux[ii][1])
						 *byhz
				 -(Pc[ii][1]-Pc[ii][0])*byhx;
	//qux[ii][1]+= (xm[1]*muc[ii][1]*(ux[ii][2]-ux[ii][1])
	//			 -xm[0]*muc[ii][0]*(ux[ii][1]-ux[ii][0]))*byhx*byhx/x[1]
	//		     +0.25*((muc[ii+1][1]+muc[ii+1][0]+muc[ii+1][1]+muc[ii+1][0])*(8.0*(ux[ii+1][1]-ux[ii][1])+(ux[ii-1][1]-ux[ii][1]))/3.0
	//				   -(muc[ii][1]+muc[ii][0]+muc[ii-1][1]+muc[ii-1][0])*(ux[ii][1]-ux[ii-1][1]))*byhz*byhz
	//		     -0.5*(muc[ii][1]+muc[ii][0])*ux[ii][1]/(x[1]*x[1])
	//		     -ux[ii][1]*(muc[ii][1]-muc[ii][0])*byhx/x[1]
	//		     +0.125*((2.0*(muc[ii+1][1]+muc[ii+1][0])-(muc[ii][1]+muc[ii][0]+muc[ii-1][1]+muc[ii-1][0]))
	//				   *((uz[ii+1][1]+uz[ii][1])-(uz[ii+1][0]+uz[ii][0])))*byhz*byhx
	//		     -0.5*(muc[ii][1]-muc[ii][0])*((uz[ii+1][1]+uz[ii+1][0])-(uz[ii][1]+uz[ii][0]))*byhz*byhx
	//		     +1.0/3.0*(muc[ii][1]*(x[2]*ux[ii][2]-x[1]*ux[ii][1])/xm[1]
	//					  -muc[ii][0]*(x[1]*ux[ii][1]-x[0]*ux[ii][0])/xm[0])*byhx*byhx
	//		     +1.0/3.0*(muc[ii][1]*(uz[ii+1][1]-uz[ii][1])-muc[ii][0]*(uz[ii+1][0]-uz[ii][0]))*byhx*byhz;
	qux[ii][1]*= dt;
	qux[ii][1]+= 0.5*(rogc[ii][1]+rogc[ii][0])*ux[ii][1];

	/* i==NZ-1&&j==NX-1 */
	qux[ii][jj] = -0.125*(xm[jj]*(fabs(ux[ii][jj+1]+ux[ii][jj])+(ux[ii][jj+1]+ux[ii][jj]))*(rogc[ii][jj]+rogc[ii][jj-1])*ux[ii][jj]
					    -xm[jj-1]*
					    ((fabs(ux[ii][jj]+ux[ii][jj-1])+(ux[ii][jj]+ux[ii][jj-1]))*(rogc[ii][jj-1]+rogc[ii][jj-2])*ux[ii][jj-1]
					    -(fabs(ux[ii][jj]+ux[ii][jj-1])-(ux[ii][jj]+ux[ii][jj-1]))*(rogc[ii][jj]+rogc[ii][jj-1])*ux[ii][jj]))
					    *byhx/x[jj]
				 +0.125*((fabs(uz[ii][jj]+uz[ii][jj-1])+(uz[ii][jj]+uz[ii][jj-1]))*(rogc[ii-1][jj]+rogc[ii-1][jj-1])*ux[ii-1][jj]
						-(fabs(uz[ii][jj]+uz[ii][jj-1])-(uz[ii][jj]+uz[ii][jj-1]))*(rogc[ii][jj]+rogc[ii][jj-1])*ux[ii][jj])
						*byhz
				 -(Pc[ii][jj]-Pc[ii][jj-1])*byhx;
	//qux[ii][jj]+= (xm[jj]*muc[ii][jj]*(ux[ii][jj+1]-ux[ii][jj])
	//			 -xm[jj-1]*muc[ii][jj-1]*(ux[ii][jj]-ux[ii][jj-1]))*byhx*byhx/x[jj]
	//		     +0.25*((muc[ii+1][jj]+muc[ii+1][jj-1]+muc[ii+1][jj]+muc[ii+1][jj-1])*(8.0*(ux[ii+1][jj]-ux[ii][jj])+(ux[ii-1][jj]-ux[ii][jj]))/3.0
	//				   -(muc[ii][jj]+muc[ii][jj-1]+muc[ii-1][jj]+muc[ii-1][jj-1])*(ux[ii][jj]-ux[ii-1][jj]))*byhz*byhz
	//		     -0.5*(muc[ii][jj]+muc[ii][jj-1])*ux[ii][jj]/(x[jj]*x[jj])
	//		     -ux[ii][jj]*(muc[ii][jj]-muc[ii][jj-1])*byhx/x[jj]
	//		     +0.125*((2.0*(muc[ii+1][jj]+muc[ii+1][jj-1])-(muc[ii][jj]+muc[ii][jj-1]+muc[ii-1][jj]+muc[ii-1][jj-1]))
	//				   *((uz[ii+1][jj]+uz[ii][jj])-(uz[ii+1][jj-1]+uz[ii][jj-1])))*byhz*byhx
	//		     -0.5*(muc[ii][jj]-muc[ii][jj-1])*((uz[ii+1][jj]+uz[ii+1][jj-1])-(uz[ii][jj]+uz[ii][jj-1]))*byhz*byhx
	//		     +1.0/3.0*(muc[ii][jj]*(x[jj+1]*ux[ii][jj+1]-x[jj]*ux[ii][jj])/xm[jj]
	//					  -muc[ii][jj-1]*(x[jj]*ux[ii][jj]-x[jj-1]*ux[ii][jj-1])/xm[jj-1])*byhx*byhx
	//		     +1.0/3.0*(muc[ii][jj]*(uz[ii+1][jj]-uz[ii][jj])-muc[ii][jj-1]*(uz[ii+1][jj-1]-uz[ii][jj-1]))*byhx*byhz;
	qux[ii][jj]*= dt;
	qux[ii][jj]+= 0.5*(rogc[ii][jj]+rogc[ii][jj-1])*ux[ii][jj];

	for(i=1;i<NZ-1;i++)
	{
		if(cxflg[i][1]==0)
		{
			/* j==1 */
			qux[i][1] = -0.125*(xm[1]*
							   ((fabs(ux[i][2]+ux[i][1])+(ux[i][2]+ux[i][1]))*(rogc[i][1]+rogc[i][0])*ux[i][1]
							   -(fabs(ux[i][2]+ux[i][1])-(ux[i][2]+ux[i][1]))*(rogc[i][2]+rogc[i][1])*ux[i][2])
							   +xm[0]*(fabs(ux[i][1]+ux[i][0])-(ux[i][1]+ux[i][0]))*(rogc[i][1]+rogc[i][0])*ux[i][1])
							   *byhx/x[1]
						-0.125*(((fabs(uz[i+1][1]+uz[i+1][0])+(uz[i+1][1]+uz[i+1][0]))*(rogc[i][1]+rogc[i][0])*ux[i][1]
								-(fabs(uz[i+1][1]+uz[i+1][0])-(uz[i+1][1]+uz[i+1][0]))*(rogc[i+1][1]+rogc[i+1][0])*ux[i+1][1])
							   -((fabs(uz[i][1]+uz[i][0])+(uz[i][1]+uz[i][0]))*(rogc[i-1][1]+rogc[i-1][0])*ux[i-1][1]
								-(fabs(uz[i][1]+uz[i][0])-(uz[i][1]+uz[i][0]))*(rogc[i][1]+rogc[i][0])*ux[i][1]))
							   *byhz
						-(Pc[i][1]-Pc[i][0])*byhx;
			//qux[i][1]+= (xm[1]*muc[i][1]*(ux[i][2]-ux[i][1])
			//			-xm[0]*muc[i][0]*(ux[i][1]-ux[i][0]))*byhx*byhx/x[1]
			//		   +0.25*((muc[i+1][1]+muc[i+1][0]+muc[i][1]+muc[i][0])*(ux[i+1][1]-ux[i][1])
			//				 -(muc[i][1]+muc[i][0]+muc[i-1][1]+muc[i-1][0])*(ux[i][1]-ux[i-1][1]))*byhz*byhz
			//		   -0.5*(muc[i][1]+muc[i][0])*ux[i][1]/(x[1]*x[1])
			//		   -ux[i][1]*(muc[i][1]-muc[i][0])*byhx/x[1]
			//		   +0.125*(((muc[i+1][1]+muc[i+1][0])-(muc[i-1][1]+muc[i-1][0]))*((uz[i+1][1]+uz[i][1])-(uz[i+1][0]+uz[i][0])))*byhz*byhx
			//		   -0.5*(muc[i][1]-muc[i][0])*((uz[i+1][1]+uz[i+1][0])-(uz[i][1]+uz[i][0]))*byhz*byhx
			//		   +1.0/3.0*(muc[i][1]*(x[2]*ux[i][2]-x[1]*ux[i][1])/xm[1]
			//					-muc[i][0]*(x[1]*ux[i][1]-x[0]*ux[i][0])/xm[0])*byhx*byhx
			//		   +1.0/3.0*(muc[i][1]*(uz[i+1][1]-uz[i][1])-muc[i][0]*(uz[i+1][0]-uz[i][0]))*byhx*byhz;
			qux[i][1]*= dt;
			qux[i][1]+= 0.5*(rogc[i][1]+rogc[i][0])*ux[i][1];
		}
		/* j==NX-1 */
		qux[i][jj] = -0.125*(xm[jj]*(fabs(ux[i][jj+1]+ux[i][jj])+(ux[i][jj+1]+ux[i][jj]))*(rogc[i][jj]+rogc[i][jj-1])*ux[i][jj]
						    -xm[jj-1]*
						    ((fabs(ux[i][jj]+ux[i][jj-1])+(ux[i][jj]+ux[i][jj-1]))*(rogc[i][jj-1]+rogc[i][jj-2])*ux[i][jj-1]
						    -(fabs(ux[i][jj]+ux[i][jj-1])-(ux[i][jj]+ux[i][jj-1]))*(rogc[i][jj]+rogc[i][jj-1])*ux[i][jj]))
						    *byhx/x[jj]
					 -0.125*(((fabs(uz[i+1][jj]+uz[i+1][jj-1])+(uz[i+1][jj]+uz[i+1][jj-1]))*(rogc[i][jj]+rogc[i][jj-1])*ux[i][jj]
							 -(fabs(uz[i+1][jj]+uz[i+1][jj-1])-(uz[i+1][jj]+uz[i+1][jj-1]))*(rogc[i+1][jj]+rogc[i+1][jj-1])*ux[i+1][jj])
							-((fabs(uz[i][jj]+uz[i][jj-1])+(uz[i][jj]+uz[i][jj-1]))*(rogc[i-1][jj]+rogc[i-1][jj-1])*ux[i-1][jj]
							 -(fabs(uz[i][jj]+uz[i][jj-1])-(uz[i][jj]+uz[i][jj-1]))*(rogc[i][jj]+rogc[i][jj-1])*ux[i][jj]))
							*byhz
					 -(Pc[i][jj]-Pc[i][jj-1])*byhx;
		//qux[i][jj]+= (xm[jj]*muc[i][jj]*(ux[i][jj+1]-ux[i][jj])
		//			 -xm[jj-1]*muc[i][jj-1]*(ux[i][jj]-ux[i][jj-1]))*byhx*byhx/x[jj]
		//		     +0.25*((muc[i+1][jj]+muc[i+1][jj-1]+muc[i][jj]+muc[i][jj-1])*(ux[i+1][jj]-ux[i][jj])
		//				   -(muc[i][jj]+muc[i][jj-1]+muc[i-1][jj]+muc[i-1][jj-1])*(ux[i][jj]-ux[i-1][jj]))*byhz*byhz
		//		     -0.5*(muc[i][jj]+muc[i][jj-1])*ux[i][jj]/(x[jj]*x[jj])
		//		     -ux[i][jj]*(muc[i][jj]-muc[i][jj-1])*byhx/x[jj]
		//		     +0.125*(((muc[i+1][jj]+muc[i+1][jj-1])-(muc[i-1][jj]+muc[i-1][jj-1]))
		//				   *((uz[i+1][jj]+uz[i][jj])-(uz[i+1][jj-1]+uz[i][jj-1])))*byhz*byhx
		//		     -0.5*(muc[i][jj]-muc[i][jj-1])*((uz[i+1][jj]+uz[i+1][jj-1])-(uz[i][jj]+uz[i][jj-1]))*byhz*byhx
		//		     +1.0/3.0*(muc[i][jj]*(x[jj+1]*ux[i][jj+1]-x[jj]*ux[i][jj])/xm[jj]
		//					  -muc[i][jj-1]*(x[jj]*ux[i][jj]-x[jj-1]*ux[i][jj-1])/xm[jj-1])*byhx*byhx
		//		     +1.0/3.0*(muc[i][jj]*(uz[i+1][jj]-uz[i][jj])-muc[i][jj-1]*(uz[i+1][jj-1]-uz[i][jj-1]))*byhx*byhz;
		qux[i][jj]*= dt;
		qux[i][jj]+= 0.5*(rogc[i][jj]+rogc[i][jj-1])*ux[i][jj];
	}

	for(j=2;j<NX-1;j++)
	{
		/*i==0*/
		qux[0][j] = -0.125*(xm[j]*
						    ((fabs(ux[0][j+1]+ux[0][j])+(ux[0][j+1]+ux[0][j]))*(rogc[0][j]+rogc[0][j-1])*ux[0][j]
						    -(fabs(ux[0][j+1]+ux[0][j])-(ux[0][j+1]+ux[0][j]))*(rogc[0][j+1]+rogc[0][j])*ux[0][j+1])
						    -xm[j-1]*
						    ((fabs(ux[0][j]+ux[0][j-1])+(ux[0][j]+ux[0][j-1]))*(rogc[0][j-1]+rogc[0][j-2])*ux[0][j-1]
						    -(fabs(ux[0][j]+ux[0][j-1])-(ux[0][j]+ux[0][j-1]))*(rogc[0][j]+rogc[0][j-1])*ux[0][j]))
						    *byhx/x[j]
					 -0.125*((fabs(uz[1][j]+uz[1][j-1])+(uz[1][j]+uz[1][j-1]))*(rogc[0][j]+rogc[0][j-1])*ux[0][j]
							-(fabs(uz[1][j]+uz[1][j-1])-(uz[1][j]+uz[1][j-1]))*(rogc[1][j]+rogc[1][j-1])*ux[1][j])
							*byhz
					 -(Pc[0][j]-Pc[0][j-1])*byhx;
		//qux[0][j]+= (xm[j]*muc[0][j]*(ux[0][j+1]-ux[0][j])
		//			-xm[j-1]*muc[0][j-1]*(ux[0][j]-ux[0][j-1]))*byhx*byhx/x[j]
		//		   +0.25*(muc[1][j]+muc[1][j-1]+muc[0][j]+muc[0][j-1])*(ux[1][j]-ux[0][j])*byhz*byhz
		//		   -0.5*(muc[0][j]+muc[0][j-1])*ux[0][j]/(x[j]*x[j])
		//		   -ux[0][j]*(muc[0][j]-muc[0][j-1])*byhx/x[j]
		//		   +0.125*(((muc[1][j]+muc[1][j-1])-(muc[0][j]+muc[0][j-1]))*((uz[1][j]+uz[0][j])-(uz[1][j-1]+uz[0][j-1])))*byhz*byhx
		//		   -0.5*(muc[0][j]-muc[0][j-1])*((uz[1][j]+uz[1][j-1])-(uz[0][j]+uz[0][j-1]))*byhz*byhx
		//		   +1.0/3.0*(muc[0][j]*(x[j+1]*ux[0][j+1]-x[j]*ux[0][j])/xm[j]
		//					-muc[0][j-1]*(x[j]*ux[0][j]-x[j-1]*ux[0][j-1])/xm[j-1])*byhx*byhx
		//		   +1.0/3.0*(muc[0][j]*(uz[1][j]-uz[0][j])-muc[0][j-1]*(uz[1][j-1]-uz[0][j-1]))*byhx*byhz;
		qux[0][j]*= dt;
		qux[0][j]+= 0.5*(rogc[0][j]+rogc[0][j-1])*ux[0][j];

		/*i==NZ-1*/
		qux[ii][j] = -0.125*(xm[j]*
							 ((fabs(ux[ii][j+1]+ux[ii][j])+(ux[ii][j+1]+ux[ii][j]))*(rogc[ii][j]+rogc[ii][j-1])*ux[ii][j]
							 -(fabs(ux[ii][j+1]+ux[ii][j])-(ux[ii][j+1]+ux[ii][j]))*(rogc[ii][j+1]+rogc[ii][j])*ux[ii][j+1])
							 -xm[j-1]*
							 ((fabs(ux[ii][j]+ux[ii][j-1])+(ux[ii][j]+ux[ii][j-1]))*(rogc[ii][j-1]+rogc[ii][j-2])*ux[ii][j-1]
							 -(fabs(ux[ii][j]+ux[ii][j-1])-(ux[ii][j]+ux[ii][j-1]))*(rogc[ii][j]+rogc[ii][j-1])*ux[ii][j]))
							 *byhx/x[j]
					  +0.125*((fabs(uz[ii][j]+uz[ii][j-1])+(uz[ii][j]+uz[ii][j-1]))*(rogc[ii-1][j]+rogc[ii-1][j-1])*ux[ii-1][j]
							  -(fabs(uz[ii][j]+uz[ii][j-1])-(uz[ii][j]+uz[ii][j-1]))*(rogc[ii][j]+rogc[ii][j-1])*ux[ii][j])
							  *byhz
					 -(Pc[ii][j]-Pc[ii][j-1])*byhx;
		//qux[ii][j]+= (xm[j]*muc[ii][j]*(ux[ii][j+1]-ux[ii][j])
		//			 -xm[j-1]*muc[ii][j-1]*(ux[ii][j]-ux[ii][j-1]))*byhx*byhx/x[j]
		//		     +0.25*((muc[ii+1][j]+muc[ii+1][j-1]+muc[ii+1][j]+muc[ii+1][j-1])*(8.0*(ux[ii+1][j]-ux[ii][j])+(ux[ii-1][j]-ux[ii][j]))/3.0
		//				   -(muc[ii][j]+muc[ii][j-1]+muc[ii-1][j]+muc[ii-1][j-1])*(ux[ii][j]-ux[ii-1][j]))*byhz*byhz
		//		     -0.5*(muc[ii][j]+muc[ii][j-1])*ux[ii][j]/(x[j]*x[j])
		//		     -ux[ii][j]*(muc[ii][j]-muc[ii][j-1])*byhx/x[j]
		//		     +0.125*((2.0*(muc[ii+1][j]+muc[ii+1][j-1])-(muc[ii][j]+muc[ii][j-1]+muc[ii-1][j]+muc[ii-1][j-1]))
		//				   *((uz[ii+1][j]+uz[ii][j])-(uz[ii+1][j-1]+uz[ii][j-1])))*byhz*byhx
		//		     -0.5*(muc[ii][j]-muc[ii][j-1])*((uz[ii+1][j]+uz[ii+1][j-1])-(uz[ii][j]+uz[ii][j-1]))*byhz*byhx
		//		     +1.0/3.0*(muc[ii][j]*(x[j+1]*ux[ii][j+1]-x[j]*ux[ii][j])/xm[j]
		//					  -muc[ii][j-1]*(x[j]*ux[ii][j]-x[j-1]*ux[ii][j-1])/xm[j-1])*byhx*byhx
		//		     +1.0/3.0*(muc[ii][j]*(uz[ii+1][j]-uz[ii][j])-muc[ii][j-1]*(uz[ii+1][j-1]-uz[ii][j-1]))*byhx*byhz;
		qux[ii][j]*= dt;
		qux[ii][j]+= 0.5*(rogc[ii][j]+rogc[ii][j-1])*ux[ii][j];
	}

	for(i=1;i<NZ-1;i++)
	{
		for(j=2;j<NX-1;j++)
		{
			if(cxflg[i][j]==0)
			{
				qux[i][j] = -0.125*(xm[j]*
								   ((fabs(ux[i][j+1]+ux[i][j])+(ux[i][j+1]+ux[i][j]))*(rogc[i][j]+rogc[i][j-1])*ux[i][j]
								   -(fabs(ux[i][j+1]+ux[i][j])-(ux[i][j+1]+ux[i][j]))*(rogc[i][j+1]+rogc[i][j])*ux[i][j+1])
								   -xm[j-1]*
								   ((fabs(ux[i][j]+ux[i][j-1])+(ux[i][j]+ux[i][j-1]))*(rogc[i][j-1]+rogc[i][j-2])*ux[i][j-1]
								   -(fabs(ux[i][j]+ux[i][j-1])-(ux[i][j]+ux[i][j-1]))*(rogc[i][j]+rogc[i][j-1])*ux[i][j]))
								   *byhx/x[j]
							-0.125*(((fabs(uz[i+1][j]+uz[i+1][j-1])+(uz[i+1][j]+uz[i+1][j-1]))*(rogc[i][j]+rogc[i][j-1])*ux[i][j]
									-(fabs(uz[i+1][j]+uz[i+1][j-1])-(uz[i+1][j]+uz[i+1][j-1]))*(rogc[i+1][j]+rogc[i+1][j-1])*ux[i+1][j])
								   -((fabs(uz[i][j]+uz[i][j-1])+(uz[i][j]+uz[i][j-1]))*(rogc[i-1][j]+rogc[i-1][j-1])*ux[i-1][j]
									-(fabs(uz[i][j]+uz[i][j-1])-(uz[i][j]+uz[i][j-1]))*(rogc[i][j]+rogc[i][j-1])*ux[i][j]))
								   *byhz
							-(Pc[i][j]-Pc[i][j-1])*byhx;
				//qux[i][j]+= (xm[j]*muc[i][j]*(ux[i][j+1]-ux[i][j])
				//			-xm[j-1]*muc[i][j-1]*(ux[i][j]-ux[i][j-1]))*byhx*byhx/x[j]
				//		   +0.25*((muc[i+1][j]+muc[i+1][j-1]+muc[i][j]+muc[i][j-1])*(ux[i+1][j]-ux[i][j])
				//				 -(muc[i][j]+muc[i][j-1]+muc[i-1][j]+muc[i-1][j-1])*(ux[i][j]-ux[i-1][j]))*byhz*byhz
				//		   -0.5*(muc[i][j]+muc[i][j-1])*ux[i][j]/(x[j]*x[j])
				//		   -ux[i][j]*(muc[i][j]-muc[i][j-1])*byhx/x[j]
				//		   +0.125*(((muc[i+1][j]+muc[i+1][j-1])-(muc[i-1][j]+muc[i-1][j-1]))*((uz[i+1][j]+uz[i][j])-(uz[i+1][j-1]+uz[i][j-1])))*byhz*byhx
				//		   -0.5*(muc[i][j]-muc[i][j-1])*((uz[i+1][j]+uz[i+1][j-1])-(uz[i][j]+uz[i][j-1]))*byhz*byhx
				//		   +1.0/3.0*(muc[i][j]*(x[j+1]*ux[i][j+1]-x[j]*ux[i][j])/xm[j]
				//					-muc[i][j-1]*(x[j]*ux[i][j]-x[j-1]*ux[i][j-1])/xm[j-1])*byhx*byhx
				//		   +1.0/3.0*(muc[i][j]*(uz[i+1][j]-uz[i][j])-muc[i][j-1]*(uz[i+1][j-1]-uz[i][j-1]))*byhx*byhz;
				qux[i][j]*= dt;
				qux[i][j]+= 0.5*(rogc[i][j]+rogc[i][j-1])*ux[i][j];
			}
		}
	}
}

void uz_cal()
{
	int i,j,ii=NZ-1,jj=NX-1;

	/*i==1&&j==0*/
	quz[1][0] =  -0.125*x[1]*
					   ((fabs(ux[1][1]+ux[0][1])+(ux[1][1]+ux[0][1]))*(rogc[1][0]+rogc[0][0])*uz[1][0]
					   -(fabs(ux[1][1]+ux[0][1])-(ux[1][1]+ux[0][1]))*(rogc[1][1]+rogc[0][1])*uz[1][1])
					   *byhx/xm[0]
			     -0.125*(((fabs(uz[2][0]+uz[1][0])+(uz[2][0]+uz[1][0]))*(rogc[1][0]+rogc[0][0])*uz[1][0]
						 -(fabs(uz[2][0]+uz[1][0])-(uz[2][0]+uz[1][0]))*(rogc[2][0]+rogc[1][0])*uz[2][0])
					    +(fabs(uz[1][0]+uz[0][0])-(uz[1][0]+uz[0][0]))*(rogc[1][0]+rogc[0][0])*uz[1][0])
						*byhz
				 -(Pc[1][0]-Pc[0][0])*byhz;
	//quz[1][0]+= 0.25*x[1]*(muc[1][1]+muc[0][1]+muc[1][0]+muc[0][0])*(uz[1][1]-uz[1][0])*byhx*byhx/xm[0]
	//		   +(muc[1][0]*(uz[2][0]-uz[1][0])-muc[0][0]*(uz[1][0]-uz[0][0]))*byhz*byhz
	//		   +0.125*((muc[1][1]+muc[0][1])-(muc[1][0]+muc[0][0]))*((ux[1][1]+ux[1][0])-(ux[0][1]+ux[0][0]))*byhx*byhz
	//		   -0.5*(muc[1][0]-muc[0][0])*((ux[1][1]+ux[0][1])-(ux[1][0]+ux[0][0]))*byhx*byhz
	//		   -0.25*(ux[1][1]+ux[0][1]+ux[1][0]+ux[0][0])*(muc[1][0]-muc[0][0])*byhz/xm[0]
	//		   +1.0/3.0*(muc[1][0]*(x[1]*ux[1][1]-x[0]*ux[1][0])-muc[0][0]*(x[1]*ux[0][1]-x[0]*ux[0][0]))*byhx*byhz/xm[0]
	//		   +1.0/3.0*(muc[1][0]*(uz[2][0]-uz[1][0])-muc[0][0]*(uz[1][0]-uz[0][0]))*byhz*byhz;
	quz[1][0]*= dt;
	quz[1][0]+= 0.5*(rogc[1][0]+rogc[0][0])*uz[1][0];

	/*i==1&&j==NX-1*/
	quz[1][jj] =   0.125*x[jj]*
						((fabs(ux[1][jj]+ux[0][jj])+(ux[1][jj]+ux[0][jj]))*(rogc[1][jj-1]+rogc[0][jj-1])*uz[1][jj-1]
						-(fabs(ux[1][jj]+ux[0][jj])-(ux[1][jj]+ux[0][jj]))*(rogc[1][jj]+rogc[0][jj])*uz[1][jj])
						*byhx/xm[jj]
			      -0.125*(((fabs(uz[2][jj]+uz[1][jj])+(uz[2][jj]+uz[1][jj]))*(rogc[1][jj]+rogc[0][jj])*uz[1][jj]
						  -(fabs(uz[2][jj]+uz[1][jj])-(uz[2][jj]+uz[1][jj]))*(rogc[2][jj]+rogc[1][jj])*uz[2][jj])
					     +(fabs(uz[1][jj]+uz[0][jj])-(uz[1][jj]+uz[0][jj]))*(rogc[1][jj]+rogc[0][jj])*uz[1][jj])
						*byhz
				 -(Pc[1][jj]-Pc[0][jj])*byhz;
	//quz[1][jj]+= 0.25*(2.0*x[jj+1]*(muc[1][jj+1]+muc[0][jj+1])*(8.0*(uz[1][jj+1]-uz[1][jj])+(uz[1][jj-1]-uz[1][jj]))/3.0
	//				 -x[jj]*(muc[1][jj]+muc[0][jj]+muc[1][jj-1]+muc[0][jj-1])*(uz[1][jj]-uz[1][jj-1]))*byhx*byhx/xm[jj]
	//		   +(muc[1][jj]*(uz[2][jj]-uz[1][jj])-muc[0][jj]*(uz[1][jj]-uz[0][jj]))*byhz*byhz
	//		   +0.125*(2.0*(muc[1][jj+1]+muc[0][jj+1])-(muc[1][jj]+muc[0][jj]+muc[1][jj-1]+muc[0][jj-1]))
	//				 *((ux[1][jj+1]+ux[1][jj])-(ux[0][jj+1]+ux[0][jj]))*byhx*byhz
	//		   -0.5*(muc[1][jj]-muc[0][jj])*((ux[1][jj+1]+ux[0][jj+1])-(ux[1][jj]+ux[0][jj]))*byhx*byhz
	//		   -0.25*(ux[1][jj+1]+ux[0][jj+1]+ux[1][jj]+ux[0][jj])*(muc[1][jj]-muc[0][jj])*byhz/xm[jj]
	//		   +1.0/3.0*(muc[1][jj]*(x[jj+1]*ux[1][jj+1]-x[jj]*ux[1][jj])-muc[0][jj]*(x[jj+1]*ux[0][jj+1]-x[jj]*ux[0][jj]))*byhx*byhz/xm[jj]
	//		   +1.0/3.0*(muc[1][jj]*(uz[2][jj]-uz[1][jj])-muc[0][jj]*(uz[1][jj]-uz[0][jj]))*byhz*byhz;
	quz[1][jj]*= dt;
	quz[1][jj]+= 0.5*(rogc[1][jj]+rogc[0][jj])*uz[1][jj];

	/*i==NZ-1&&j==0*/
	quz[ii][0] = -0.125*x[1]*
						((fabs(ux[ii][1]+ux[ii-1][1])+(ux[ii][1]+ux[ii-1][1]))*(rogc[ii][0]+rogc[ii-1][0])*uz[ii][0]
						-(fabs(ux[ii][1]+ux[ii-1][1])-(ux[ii][1]+ux[ii-1][1]))*(rogc[ii][1]+rogc[ii-1][1])*uz[ii][1])
						*byhx/xm[0]
			     -0.125*((fabs(uz[ii+1][0]+uz[ii][0])+(uz[ii+1][0]+uz[ii][0]))*(rogc[ii][0]+rogc[ii-1][0])*uz[ii][0]
					    -((fabs(uz[ii][0]+uz[ii-1][0])+(uz[ii][0]+uz[ii-1][0]))*(rogc[ii-1][0]+rogc[ii-2][0])*uz[ii-1][0]
						 -(fabs(uz[ii][0]+uz[ii-1][0])-(uz[ii][0]+uz[ii-1][0]))*(rogc[ii][0]+rogc[ii-1][0])*uz[ii][0]))
						*byhz
				 -(Pc[ii][0]-Pc[ii-1][0])*byhz;
	//quz[ii][0]+= 0.25*x[1]*(muc[ii][1]+muc[ii-1][1]+muc[ii][0]+muc[ii-1][0])*(uz[ii][1]-uz[ii][0])*byhx*byhx/xm[0]
	//		    +(muc[ii][0]*(uz[ii+1][0]-uz[ii][0])-muc[ii-1][0]*(uz[ii][0]-uz[ii-1][0]))*byhz*byhz
	//		    +0.125*((muc[ii][1]+muc[ii-1][1])-(muc[ii][0]+muc[ii-1][0]))*((ux[ii][1]+ux[ii][0])-(ux[ii-1][1]+ux[ii-1][0]))*byhx*byhz
	//		    -0.5*(muc[ii][0]-muc[ii-1][0])*((ux[ii][1]+ux[ii-1][1])-(ux[ii][0]+ux[ii-1][0]))*byhx*byhz
	//		    -0.25*(ux[ii][1]+ux[ii-1][1]+ux[ii][0]+ux[ii-1][0])*(muc[ii][0]-muc[ii-1][0])*byhz/xm[0]
	//		    +1.0/3.0*(muc[ii][0]*(x[1]*ux[ii][1]-x[0]*ux[ii][0])-muc[ii-1][0]*(x[1]*ux[ii-1][1]-x[0]*ux[ii-1][0]))*byhx*byhz/xm[0]
	//		    +1.0/3.0*(muc[ii][0]*(uz[ii+1][0]-uz[ii][0])-muc[ii-1][0]*(uz[ii][0]-uz[ii-1][0]))*byhz*byhz;
	quz[ii][0]*= dt;
	quz[ii][0]+= 0.5*(rogc[ii][0]+rogc[ii-1][0])*uz[ii][0];

	/*i==NZ-1&&j==NX-1*/
	quz[ii][jj] =   0.125*x[jj]*
						 ((fabs(ux[ii][jj]+ux[ii-1][jj])+(ux[ii][jj]+ux[ii-1][jj]))*(rogc[ii][jj-1]+rogc[ii-1][jj-1])*uz[ii][jj-1]
						 -(fabs(ux[ii][jj]+ux[ii-1][jj])-(ux[ii][jj]+ux[ii-1][jj]))*(rogc[ii][jj]+rogc[ii-1][jj])*uz[ii][jj])
						 *byhx/xm[jj]
				   -0.125*((fabs(uz[ii+1][jj]+uz[ii][jj])+(uz[ii+1][jj]+uz[ii][jj]))*(rogc[ii][jj]+rogc[ii-1][jj])*uz[ii][jj]
					     -((fabs(uz[ii][jj]+uz[ii-1][jj])+(uz[ii][jj]+uz[ii-1][jj]))*(rogc[ii-1][jj]+rogc[ii-2][jj])*uz[ii-1][jj]
						 -(fabs(uz[ii][jj]+uz[ii-1][jj])-(uz[ii][jj]+uz[ii-1][jj]))*(rogc[ii][jj]+rogc[ii-1][jj])*uz[ii][jj]))
						 *byhz
 				 -(Pc[ii][jj]-Pc[ii-1][jj])*byhz;
	//quz[ii][jj]+= 0.25*(2.0*x[jj+1]*(muc[ii][jj+1]+muc[ii-1][jj+1])*(8.0*(uz[ii][jj+1]-uz[ii][jj])+(uz[ii][jj-1]-uz[ii][jj]))/3.0
	//			 -x[jj]*(muc[ii][jj]+muc[ii-1][jj]+muc[ii][jj-1]+muc[ii-1][jj-1])*(uz[ii][jj]-uz[ii][jj-1]))*byhx*byhx/xm[jj]
	//			 +(muc[ii][jj]*(uz[ii+1][jj]-uz[ii][jj])-muc[ii-1][jj]*(uz[ii][jj]-uz[ii-1][jj]))*byhz*byhz
	//		     +0.125*(2.0*(muc[ii][jj+1]+muc[ii-1][jj+1])-(muc[ii][jj]+muc[ii-1][jj]+muc[ii][jj-1]+muc[ii-1][jj-1]))
	//				   *((ux[ii][jj+1]+ux[ii][jj])-(ux[ii-1][jj+1]+ux[ii-1][jj]))*byhx*byhz
	//		     -0.5*(muc[ii][jj]-muc[ii-1][jj])*((ux[ii][jj+1]+ux[ii-1][jj+1])-(ux[ii][jj]+ux[ii-1][jj]))*byhx*byhz
	//		     -0.25*(ux[ii][jj+1]+ux[ii-1][jj+1]+ux[ii][jj]+ux[ii-1][jj])*(muc[ii][jj]-muc[ii-1][jj])*byhz/xm[jj]
	//		     +1.0/3.0*(muc[ii][jj]*(x[jj+1]*ux[ii][jj+1]-x[jj]*ux[ii][jj])-muc[ii-1][jj]*(x[jj+1]*ux[ii-1][jj+1]-x[jj]*ux[ii-1][jj]))*byhx*byhz/xm[jj]
	//		     +1.0/3.0*(muc[ii][jj]*(uz[ii+1][jj]-uz[ii][jj])-muc[ii-1][jj]*(uz[ii][jj]-uz[ii-1][jj]))*byhz*byhz;
	quz[ii][jj]*= dt;
	quz[ii][jj]+= 0.5*(rogc[ii][jj]+rogc[ii-1][jj])*uz[ii][jj];

	for(j=1;j<NX-1;j++)
	{
		/*i==1*/
		quz[1][j] =  -0.125*(x[j+1]*
							((fabs(ux[1][j+1]+ux[0][j+1])+(ux[1][j+1]+ux[0][j+1]))*(rogc[1][j]+rogc[0][j])*uz[1][j]
							-(fabs(ux[1][j+1]+ux[0][j+1])-(ux[1][j+1]+ux[0][j+1]))*(rogc[1][j+1]+rogc[0][j+1])*uz[1][j+1])
							-x[j]*
							((fabs(ux[1][j]+ux[0][j])+(ux[1][j]+ux[0][j]))*(rogc[1][j-1]+rogc[0][j-1])*uz[1][j-1]
							-(fabs(ux[1][j]+ux[0][j])-(ux[1][j]+ux[0][j]))*(rogc[1][j]+rogc[0][j])*uz[1][j]))
							*byhx/xm[j]
				     -0.125*(((fabs(uz[2][j]+uz[1][j])+(uz[2][j]+uz[1][j]))*(rogc[1][j]+rogc[0][j])*uz[1][j]
							 -(fabs(uz[2][j]+uz[1][j])-(uz[2][j]+uz[1][j]))*(rogc[2][j]+rogc[1][j])*uz[2][j])
						    +(fabs(uz[1][j]+uz[0][j])-(uz[1][j]+uz[0][j]))*(rogc[1][j]+rogc[0][j])*uz[1][j])
							*byhz
					 -(Pc[1][j]-Pc[0][j])*byhz;
		//quz[1][j]+= 0.25*(x[j+1]*(muc[1][j+1]+muc[0][j+1]+muc[1][j]+muc[0][j])*(uz[1][j+1]-uz[1][j])
		//				 -x[j]*(muc[1][j]+muc[0][j]+muc[1][j-1]+muc[0][j-1])*(uz[1][j]-uz[1][j-1]))*byhx*byhx/xm[j]
		//		   +(muc[1][j]*(uz[2][j]-uz[1][j])-muc[0][j]*(uz[1][j]-uz[0][j]))*byhz*byhz
		//		   +0.125*((muc[1][j+1]+muc[0][j+1])-(muc[1][j-1]+muc[0][j-1]))*((ux[1][j+1]+ux[1][j])-(ux[0][j+1]+ux[0][j]))*byhx*byhz
		//		   -0.5*(muc[1][j]-muc[0][j])*((ux[1][j+1]+ux[0][j+1])-(ux[1][j]+ux[0][j]))*byhx*byhz
		//		   -0.25*(ux[1][j+1]+ux[0][j+1]+ux[1][j]+ux[0][j])*(muc[1][j]-muc[0][j])*byhz/xm[j]
		//		   +1.0/3.0*(muc[1][j]*(x[j+1]*ux[1][j+1]-x[j]*ux[1][j])-muc[0][j]*(x[j+1]*ux[0][j+1]-x[j]*ux[0][j]))*byhx*byhz/xm[j]
		//		   +1.0/3.0*(muc[1][j]*(uz[2][j]-uz[1][j])-muc[0][j]*(uz[1][j]-uz[0][j]))*byhz*byhz;
		quz[1][j]*= dt;
		quz[1][j]+= 0.5*(rogc[1][j]+rogc[0][j])*uz[1][j];

		/*i==NZ-1*/
		quz[ii][j] = -0.125*(x[j+1]*
							((fabs(ux[ii][j+1]+ux[ii-1][j+1])+(ux[ii][j+1]+ux[ii-1][j+1]))*(rogc[ii][j]+rogc[ii-1][j])*uz[ii][j]
							-(fabs(ux[ii][j+1]+ux[ii-1][j+1])-(ux[ii][j+1]+ux[ii-1][j+1]))*(rogc[ii][j+1]+rogc[ii-1][j+1])*uz[ii][j+1])
							-x[j]*
							((fabs(ux[ii][j]+ux[ii-1][j])+(ux[ii][j]+ux[ii-1][j]))*(rogc[ii][j-1]+rogc[ii-1][j-1])*uz[ii][j-1]
							-(fabs(ux[ii][j]+ux[ii-1][j])-(ux[ii][j]+ux[ii-1][j]))*(rogc[ii][j]+rogc[ii-1][j])*uz[ii][j]))
							*byhx/xm[j]
				     -0.125*((fabs(uz[ii+1][j]+uz[ii][j])+(uz[ii+1][j]+uz[ii][j]))*(rogc[ii][j]+rogc[ii-1][j])*uz[ii][j]
						    -((fabs(uz[ii][j]+uz[ii-1][j])+(uz[ii][j]+uz[ii-1][j]))*(rogc[ii-1][j]+rogc[ii-2][j])*uz[ii-1][j]
							 -(fabs(uz[ii][j]+uz[ii-1][j])-(uz[ii][j]+uz[ii-1][j]))*(rogc[ii][j]+rogc[ii-1][j])*uz[ii][j]))
							*byhz
					 -(Pc[ii][j]-Pc[ii-1][j])*byhz;
		//quz[ii][j]+= 0.25*(x[j+1]*(muc[ii][j+1]+muc[ii-1][j+1]+muc[ii][j]+muc[ii-1][j])*(uz[ii][j+1]-uz[ii][j])
		//				  -x[j]*(muc[ii][j]+muc[ii-1][j]+muc[ii][j-1]+muc[ii-1][j-1])*(uz[ii][j]-uz[ii][j-1]))*byhx*byhx/xm[j]
		//		    +(muc[ii][j]*(uz[ii+1][j]-uz[ii][j])-muc[ii-1][j]*(uz[ii][j]-uz[ii-1][j]))*byhz*byhz
		//		    +0.125*((muc[ii][j+1]+muc[ii-1][j+1])-(muc[ii][j-1]+muc[ii-1][j-1]))*((ux[ii][j+1]+ux[ii][j])-(ux[ii-1][j+1]+ux[ii-1][j]))*byhx*byhz
		//		    -0.5*(muc[ii][j]-muc[ii-1][j])*((ux[ii][j+1]+ux[ii-1][j+1])-(ux[ii][j]+ux[ii-1][j]))*byhx*byhz
		//		    -0.25*(ux[ii][j+1]+ux[ii-1][j+1]+ux[ii][j]+ux[ii-1][j])*(muc[ii][j]-muc[ii-1][j])*byhz/xm[j]
		//		    +1.0/3.0*(muc[ii][j]*(x[j+1]*ux[ii][j+1]-x[j]*ux[ii][j])-muc[ii-1][j]*(x[j+1]*ux[ii-1][j+1]-x[j]*ux[ii-1][j]))*byhx*byhz/xm[j]
		//		    +1.0/3.0*(muc[ii][j]*(uz[ii+1][j]-uz[ii][j])-muc[ii-1][j]*(uz[ii][j]-uz[ii-1][j]))*byhz*byhz;
		quz[ii][j]*= dt;
		quz[ii][j]+= 0.5*(rogc[ii][j]+rogc[ii-1][j])*uz[ii][j];
	}

	for(i=2;i<NZ-1;i++)
	{
		if(czflg[i][0]==0)
		{
			/*j==0*/
			quz[i][0] =  -0.125*x[1]*
							   ((fabs(ux[i][1]+ux[i-1][1])+(ux[i][1]+ux[i-1][1]))*(rogc[i][0]+rogc[i-1][0])*uz[i][0]
							   -(fabs(ux[i][1]+ux[i-1][1])-(ux[i][1]+ux[i-1][1]))*(rogc[i][1]+rogc[i-1][1])*uz[i][1])
							   *byhx/xm[0]
						 -0.125*(((fabs(uz[i+1][0]+uz[i][0])+(uz[i+1][0]+uz[i][0]))*(rogc[i][0]+rogc[i-1][0])*uz[i][0]
								 -(fabs(uz[i+1][0]+uz[i][0])-(uz[i+1][0]+uz[i][0]))*(rogc[i+1][0]+rogc[i][0])*uz[i+1][0])
							    -((fabs(uz[i][0]+uz[i-1][0])+(uz[i][0]+uz[i-1][0]))*(rogc[i-1][0]+rogc[i-2][0])*uz[i-1][0]
								 -(fabs(uz[i][0]+uz[i-1][0])-(uz[i][0]+uz[i-1][0]))*(rogc[i][0]+rogc[i-1][0])*uz[i][0]))
							   *byhz
						 -(Pc[i][0]-Pc[i-1][0])*byhz;
			//quz[i][0]+= 0.25*x[1]*(muc[i][1]+muc[i-1][1]+muc[i][0]+muc[i-1][0])*(uz[i][1]-uz[i][0])*byhx*byhx/xm[0]
			//		   +(muc[i][0]*(uz[i+1][0]-uz[i][0])-muc[i-1][0]*(uz[i][0]-uz[i-1][0]))*byhz*byhz
			//		   +0.125*((muc[i][1]+muc[i-1][1])-(muc[i][0]+muc[i-1][0]))*((ux[i][1]+ux[i][0])-(ux[i-1][1]+ux[i-1][0]))*byhx*byhz
			//		   -0.5*(muc[i][0]-muc[i-1][0])*((ux[i][1]+ux[i-1][1])-(ux[i][0]+ux[i-1][0]))*byhx*byhz
			//		   -0.25*(ux[i][1]+ux[i-1][1]+ux[i][0]+ux[i-1][0])*(muc[i][0]-muc[i-1][0])*byhz/xm[0]
			//		   +1.0/3.0*(muc[i][0]*(x[1]*ux[i][1]-x[0]*ux[i][0])-muc[i-1][0]*(x[1]*ux[i-1][1]-x[0]*ux[i-1][0]))*byhx*byhz/xm[0]
			//		   +1.0/3.0*(muc[i][0]*(uz[i+1][0]-uz[i][0])-muc[i-1][0]*(uz[i][0]-uz[i-1][0]))*byhz*byhz;
			quz[i][0]*= dt;
			quz[i][0]+= 0.5*(rogc[i][0]+rogc[i-1][0])*uz[i][0];
		}
	
		/*j==NX-1*/
		quz[i][jj] =  0.125*x[jj]*
							((fabs(ux[i][jj]+ux[i-1][jj])+(ux[i][jj]+ux[i-1][jj]))*(rogc[i][jj-1]+rogc[i-1][jj-1])*uz[i][jj-1]
							-(fabs(ux[i][jj]+ux[i-1][jj])-(ux[i][jj]+ux[i-1][jj]))*(rogc[i][jj]+rogc[i-1][jj])*uz[i][jj])
							*byhx/xm[jj]
				      -0.125*(((fabs(uz[i+1][jj]+uz[i][jj])+(uz[i+1][jj]+uz[i][jj]))*(rogc[i][jj]+rogc[i-1][jj])*uz[i][jj]
							  -(fabs(uz[i+1][jj]+uz[i][jj])-(uz[i+1][jj]+uz[i][jj]))*(rogc[i+1][jj]+rogc[i][jj])*uz[i+1][jj])
						     -((fabs(uz[i][jj]+uz[i-1][jj])+(uz[i][jj]+uz[i-1][jj]))*(rogc[i-1][jj]+rogc[i-2][jj])*uz[i-1][jj]
							  -(fabs(uz[i][jj]+uz[i-1][jj])-(uz[i][jj]+uz[i-1][jj]))*(rogc[i][jj]+rogc[i-1][jj])*uz[i][jj]))
							 *byhz
					 -(Pc[i][jj]-Pc[i-1][jj])*byhz;
		//quz[i][jj]+= 0.25*(2.0*x[jj+1]*(muc[i][jj+1]+muc[i-1][jj+1])*(8.0*(uz[i][jj+1]-uz[i][jj])+(uz[i][jj-1]-uz[i][jj]))/3.0
		//				 -x[jj]*(muc[i][jj]+muc[i-1][jj]+muc[i][jj-1]+muc[i-1][jj-1])*(uz[i][jj]-uz[i][jj-1]))*byhx*byhx/xm[jj]
		//		   +(muc[i][jj]*(uz[i+1][jj]-uz[i][jj])-muc[i-1][jj]*(uz[i][jj]-uz[i-1][jj]))*byhz*byhz
		//		   +0.125*(2.0*(muc[i][jj+1]+muc[i-1][jj+1])-(muc[i][jj]+muc[i-1][jj]+muc[i][jj-1]+muc[i-1][jj-1]))
		//				 *((ux[i][jj+1]+ux[i][jj])-(ux[i-1][jj+1]+ux[i-1][jj]))*byhx*byhz
		//		   -0.5*(muc[i][jj]-muc[i-1][jj])*((ux[i][jj+1]+ux[i-1][jj+1])-(ux[i][jj]+ux[i-1][jj]))*byhx*byhz
		//		   -0.25*(ux[i][jj+1]+ux[i-1][jj+1]+ux[i][jj]+ux[i-1][jj])*(muc[i][jj]-muc[i-1][jj])*byhz/xm[jj]
		//		   +1.0/3.0*(muc[i][jj]*(x[jj+1]*ux[i][jj+1]-x[jj]*ux[i][jj])-muc[i-1][jj]*(x[jj+1]*ux[i-1][jj+1]-x[jj]*ux[i-1][jj]))*byhx*byhz/xm[jj]
		//		   +1.0/3.0*(muc[i][jj]*(uz[i+1][jj]-uz[i][jj])-muc[i-1][jj]*(uz[i][jj]-uz[i-1][jj]))*byhz*byhz;
		quz[i][jj]*= dt;
		quz[i][jj]+= 0.5*(rogc[i][jj]+rogc[i-1][jj])*uz[i][jj];
	}

	for(i=2;i<NZ-1;i++)
	{
		for(j=1;j<NX-1;j++)
		{
			if(czflg[i][j]==0)
			{
				quz[i][j] =  -0.125*(x[j+1]*
									((fabs(ux[i][j+1]+ux[i-1][j+1])+(ux[i][j+1]+ux[i-1][j+1]))*(rogc[i][j]+rogc[i-1][j])*uz[i][j]
									-(fabs(ux[i][j+1]+ux[i-1][j+1])-(ux[i][j+1]+ux[i-1][j+1]))*(rogc[i][j+1]+rogc[i-1][j+1])*uz[i][j+1])
									-x[j]*
									((fabs(ux[i][j]+ux[i-1][j])+(ux[i][j]+ux[i-1][j]))*(rogc[i][j-1]+rogc[i-1][j-1])*uz[i][j-1]
									-(fabs(ux[i][j]+ux[i-1][j])-(ux[i][j]+ux[i-1][j]))*(rogc[i][j]+rogc[i-1][j])*uz[i][j]))
									*byhx/xm[j]
						     -0.125*(((fabs(uz[i+1][j]+uz[i][j])+(uz[i+1][j]+uz[i][j]))*(rogc[i][j]+rogc[i-1][j])*uz[i][j]
									 -(fabs(uz[i+1][j]+uz[i][j])-(uz[i+1][j]+uz[i][j]))*(rogc[i+1][j]+rogc[i][j])*uz[i+1][j])
								    -((fabs(uz[i][j]+uz[i-1][j])+(uz[i][j]+uz[i-1][j]))*(rogc[i-1][j]+rogc[i-2][j])*uz[i-1][j]
									 -(fabs(uz[i][j]+uz[i-1][j])-(uz[i][j]+uz[i-1][j]))*(rogc[i][j]+rogc[i-1][j])*uz[i][j]))
									*byhz
							 -(Pc[i][j]-Pc[i-1][j])*byhz;
				//quz[i][j]+= 0.25*(x[j+1]*(muc[i][j+1]+muc[i-1][j+1]+muc[i][j]+muc[i-1][j])*(uz[i][j+1]-uz[i][j])
				//				 -x[j]*(muc[i][j]+muc[i-1][j]+muc[i][j-1]+muc[i-1][j-1])*(uz[i][j]-uz[i][j-1]))*byhx*byhx/xm[j]
				//		   +(muc[i][j]*(uz[i+1][j]-uz[i][j])-muc[i-1][j]*(uz[i][j]-uz[i-1][j]))*byhz*byhz
				//		   +0.125*((muc[i][j+1]+muc[i-1][j+1])-(muc[i][j-1]+muc[i-1][j-1]))*((ux[i][j+1]+ux[i][j])-(ux[i-1][j+1]+ux[i-1][j]))*byhx*byhz
				//		   -0.5*(muc[i][j]-muc[i-1][j])*((ux[i][j+1]+ux[i-1][j+1])-(ux[i][j]+ux[i-1][j]))*byhx*byhz
				//		   -0.25*(ux[i][j+1]+ux[i-1][j+1]+ux[i][j]+ux[i-1][j])*(muc[i][j]-muc[i-1][j])*byhz/xm[j]
				//		   +1.0/3.0*(muc[i][j]*(x[j+1]*ux[i][j+1]-x[j]*ux[i][j])-muc[i-1][j]*(x[j+1]*ux[i-1][j+1]-x[j]*ux[i-1][j]))*byhx*byhz/xm[j]
				//		   +1.0/3.0*(muc[i][j]*(uz[i+1][j]-uz[i][j])-muc[i-1][j]*(uz[i][j]-uz[i-1][j]))*byhz*byhz;
				quz[i][j]*= dt;
				quz[i][j]+= 0.5*(rogc[i][j]+rogc[i-1][j])*uz[i][j];
			}
		}
	}
}

void density_cal()
{
	int i,j,sp;

	/* spherical coordinate system */
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl;j++)
		{
			Tls[i][j] = qTls[i][j];
		}
	}
	for(i=0;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			Tgs[i][j] = qTgs[i][j];
			for(sp=0;sp<spmax-1;sp++)
			{
				Ys[i][j][sp] = qYs[i][j][sp];
			}
			Ys[i][j][sN2] = N2sum(Ys[i][j]);
			Mms[i][j] = Mmcal(Ys[i][j]);
			qrogs[i][j] = rogcal(Mms[i][j],Tgs[i][j],Pa);
			propertygs_cal(i,j);
		}
	}
	
	/* cylindrical coordinate system */
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				Tgc[i][j] = qTgc[i][j];
				for(sp=0;sp<spmax-1;sp++)
				{
					Yc[i][j][sp] = qYc[i][j][sp];
				}
				Yc[i][j][sN2] = N2sum(Yc[i][j]);
				Mmc[i][j] = Mmcal(Yc[i][j]);
				qrogc[i][j] = rogcal(Mmc[i][j],Tgc[i][j],Pa);
				propertygc_cal(i,j);
			}
		}
	}
}

void Search_SMAC()
{
	int i,j,ii,jj,k,flg;
	double r0,s0,u1,u2,aa;
	double *b;
	double v1,v2,v3;

	b = dvec(3);

	/* from spherical to cylindrical */
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			/*qux*/
			if(cxflg[i][j]==2)
			{
				r0 = sqrt(x[j]*x[j]+(zm[i]-ld)*(zm[i]-ld));
				s0 = acos((zm[i]-ld)/r0);
				flg = 0;
				k = 0;
				do
				{
					if(rgm[k]>r0)
					{
						flg = 1;
						jj = k;
					}
					k+=1;
				}while(flg==0);
				flg = 0;
				k = 0;
				do
				{
					if(((double)k+0.5)*hs>s0)
					{
						flg = 1;
						ii = k;
					}
					k+=1;
				}while(flg==0&&k<NQ);
				if(flg==0)
				{
					ii = NQ-1;
					u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
					u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
					v1 = u2*ssm[ii]+u1*scm[ii];
					u1 = 0.0;
					u2 = 0.0625*(9.0*qur[ii][jj+1]-qur[ii-1][jj+1]+9.0*qur[ii][jj]-qur[ii-1][jj]);
					v2 = u2*ss[NQ]+u1*sc[NQ];
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii+0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v3 = u2*ssm[ii]+u1*scm[ii];
						b[0] = 2.0*(v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];
					}
					else
					{
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
						v3 = u2*ss[NQ]+u1*sc[NQ];
						b[0] = 2.0*(v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*(double)NQ*hs - b[1]*rgm[jj-1];
					}
				}
				else if(ii==0)
				{
					u1 = 0.0;
					u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
					v1 = u2*ss[0]+u1*sc[0];
					u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
					u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
					v2 = u2*ssm[ii]+u1*scm[ii];
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*s0+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj+1]-qur[1][jj+1]+9.0*qur[0][jj]-qur[1][jj]);
						v3 = u2*ss[0]+u1*sc[0];
						b[0] = 2.0*(v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[1]*rgm[jj];
					}
					else
					{
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 = u2*ssm[ii]+u1*scm[ii];
						b[0] = 2.0*(v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*0.5*hs - b[1]*rgm[jj-1];
					}
				}
				else
				{
					u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
					u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
					v1 = u2*ssm[ii-1]+u1*scm[ii-1];
					u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
					u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
					v2 = u2*ssm[ii]+u1*scm[ii];
					aa = (rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii-0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
						u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
						v3 = u2*ssm[ii-1]+u1*scm[ii-1];
						b[0] = (v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];
					}
					else
					{
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 = u2*ssm[ii]+u1*scm[ii];
						b[0] = (v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];
					}
				}
				qux[i][j] = b[0]*s0+b[1]*r0+b[2];
			}

			/*quz*/
			if(czflg[i][j]==2)
			{
				r0 = sqrt(xm[j]*xm[j]+(z[i]-ld)*(z[i]-ld));
				s0 = acos((z[i]-ld)/r0);
				flg = 0;
				k = 0;
				do
				{
					if(rgm[k]>r0)
					{
						flg = 1;
						jj = k;
					}
					k+=1;
				}while(flg==0);
				flg = 0;
				k = 0;
				do
				{
					if(((double)k+0.5)*hs>s0)
					{
						flg = 1;
						ii = k;
					}
					k+=1;
				}while(flg==0&&k<NQ);
				if(flg==0)
				{
					ii = NQ-1;
					u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
					u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
					v1 = u2*scm[ii]-u1*ssm[ii];
					u1 = 0.0;
					u2 = 0.0625*(9.0*qur[ii][jj+1]-qur[ii-1][jj+1]+9.0*qur[ii][jj]-qur[ii-1][jj]);
					v2 = u2*sc[NQ]-u1*ss[NQ];
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii+0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v3 = u2*scm[ii]-u1*ssm[ii];
						b[0] = 2.0*(v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];
					}
					else
					{	
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
						v3 = u2*sc[NQ]-u1*ss[NQ];
						b[0] = 2.0*(v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*(double)NQ*hs - b[1]*rgm[jj-1];
					}
				}
				else if(ii==0)
				{
					u1 = 0.0;
					u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
					v1 = u2*sc[0]-u1*ss[0];
					u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
					u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
					v2 = u2*scm[ii]-u1*ssm[ii];
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*s0+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj+1]-qur[1][jj+1]+9.0*qur[0][jj]-qur[1][jj]);
						v3 = u2*sc[0]-u1*ss[0];
						b[0] = 2.0*(v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[1]*rgm[jj];
					}
					else
					{
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 = u2*scm[ii]-u1*ssm[ii];
						b[0] = 2.0*(v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*0.5*hs - b[1]*rgm[jj-1];
					}
				}
				else
				{
					u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
					u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
					v1 = u2*scm[ii-1]-u1*ssm[ii-1];
					u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
					u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
					v2 = u2*scm[ii]-u1*ssm[ii];
					aa = (rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii-0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
						u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
						v3 = u2*scm[ii-1]-u1*ssm[ii-1];
						b[0] = (v2-v3)*byhs;
						b[1] = (v3-v1)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];
					}
					else
					{
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 = u2*scm[ii]-u1*ssm[ii];
						b[0] = (v3-v1)*byhs;
						b[1] = (v2-v3)/(rgm[jj]-rgm[jj-1]);
						b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];
					}
				}
				quz[i][j] = b[0]*s0+b[1]*r0+b[2];
			}
		}
	}

	freedvec(b,3);
	
}

void Search_SMAC_Spherical()
{
	int i,j,k,flg,ii,jj;
	double *b;
	double x0,z0,aa;

	b = dvec(3);

	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{
			x0 = rgm[i-1]*ssm[j];
			z0 = rgm[i-1]*scm[j] + ld;
			flg = 0;
			k = 0;
			do
			{
				if(z0<zm[k])
				{
					flg = 1;
					ii = k;
				}
				k += 1;
			}while(flg==0&&k<=NZ);
			flg = 0;
			k = 0;
			do
			{
				if(x0<xm[k])
				{
					flg = 1;
					jj = k;
				}
				k+=1;
			}while(flg==0&&k<=NX);
			if(jj==0)
			{
				aa = 0.5*hx*byhz*(z0-zm[ii-1]);
				if(x0>aa)
				{
					b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = 0.25*(dPc[ii-1][1]-dPc[ii-1][0])*byhx;
					b[2] = dPc[ii][0] - b[0]*zm[ii] - b[1]*xm[0];
				}
				else
				{
					b[0] = 0.125*(9.0*(dPc[ii][0]-dPc[ii-1][0])-(dPc[ii][1]-dPc[ii-1][1]))*byhz;
					b[1] = 0.25*(dPc[ii][1]-dPc[ii][0])*byhx;
					b[2] = dPc[ii][0] - b[0]*zm[ii] - b[1]*xm[0];
				}
			}
			else
			{
				aa = hx*byhz*(z0-zm[ii-1])+xm[jj-1];
				if(x0>aa)
				{
					b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = (dPc[ii-1][jj]-dPc[ii-1][jj-1])*byhx;
					b[2] = dPc[ii-1][jj] - b[0]*zm[ii-1] - b[1]*xm[jj];
				}
				else
				{
					b[0] = (dPc[ii][jj-1]-dPc[ii-1][jj-1])*byhz;
					b[1] = (dPc[ii][jj]-dPc[ii][jj-1])*byhx;
					b[2] = dPc[ii][jj-1] - b[0]*zm[ii] - b[1]*xm[jj-1];
				}
			}
			dPs[j][i] = b[0]*z0+b[1]*x0+b[2];
		}
	}

	freedvec(b,3);
}

void Search_SMAC_Cylindrical()
{
	int i,j,k,flg,ii,jj,g;
	double *b;
	double r0,s0,aa;

	b = dvec(3);

	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==2)
			{
				r0 = sqrt(xm[j]*xm[j]+(zm[i]-ld)*(zm[i]-ld));
				s0 = acos((zm[i]-ld)/r0);
				flg = 0;
				k = 0;
				do
				{
					if(rgm[k]>r0)
					{
						flg = 1;
						jj = k;
					}
					k+=1;
				}while(flg==0);
				flg = 0;
				k = 0;
				do
				{
					if(((double)k+0.5)*hs>s0)
					{
						flg = 1;
						ii = k;
					}
					k+=1;
				}while(flg==0&&k<NQ);
				if(flg==0)
				{
					ii = NQ-1;
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii+0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						b[0] = 0.25*(dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];
					}
					else
					{
						b[0] = 0.25*(dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = 0.125*(9.0*(dPs[ii][jj+1]-dPs[ii][jj])-(dPs[ii-1][jj+1]-dPs[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];
					}
				}
				else if(ii==0)
				{
					aa = 2.0*(rgm[jj]-rgm[jj-1])*byhs*s0+rgm[jj-1];
					if(aa>r0)
					{
						b[0] = 0.25*(dPs[1][jj+1]-dPs[0][jj+1])*byhs;
						b[1] = 0.125*(9.0*(dPs[0][jj+1]-dPs[0][jj])-(dPs[1][jj+1]-dPs[1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];
					}
					else
					{
						b[0] = 0.25*(dPs[1][jj]-dPs[0][jj])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];
					}
				}
				else
				{
					aa = (rgm[jj]-rgm[jj-1])*byhs*(s0-((double)ii-0.5)*hs)+rgm[jj-1];
					if(aa>r0)
					{
						b[0] = (dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = (dPs[ii-1][jj+1]-dPs[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii-1][jj+1] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];
					}
					else
					{
						b[0] = (dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];
					}
				}
				dPc[i][j] = b[0]*s0+b[1]*r0+b[2];
			}
		}
	}
	freedvec(b,3);
}

int SMAC_method()
{
	int i,j,m=0,flg;
	double er;
	double RdotbyRs,AA,AAold,PDIFF;
	double OMG = 0.8;

	AA = 1.0 / log(rmax/Rs);
	AAold = 1.0 / log(rmax/Rsold);
	RdotbyRs = Rdot/Rs;

	Search_SMAC();

	/*spherical*/
	for(i=0;i<NQ;i++)
	{
		Dnp[i][1] = -(qrogs[i][1]-rogs[i][1])/dt+0.5*RdotbyRs*AA*((1.0-eta[1])*(qrogs[i][1+1]-qrogs[i][1])
																 +(1.0-eta[0])*(-8.0*qrogs[i][0]+9.0*qrogs[i][1]-qrogs[i][2])/3.0)*byhrg;
		Dnm[i][1] = AA*(rg[1]*rg[1]*qur[i][1]-rg[0]*rg[0]*qur[i][0])/(rgm[0]*rgm[0]*rgm[0])*byhrg
				   +(ss[i+1]*quq[i+1][1]-ss[i]*quq[i][1])/(rgm[0]*ssm[i])*byhs;
		dPs[i][1] = 0.0;
		for(j=2;j<=Ng;j++)
		{
			Dnp[i][j] = -(qrogs[i][j]-rogs[i][j])/dt+0.5*RdotbyRs*AA*((1.0-eta[j])*(qrogs[i][j+1]-qrogs[i][j])
																	 +(1.0-eta[j-1])*(qrogs[i][j]-qrogs[i][j-1]))*byhrg;
			Dnm[i][j] = AA*(rg[j]*rg[j]*qur[i][j]-rg[j-1]*rg[j-1]*qur[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*byhrg
					   +(ss[i+1]*quq[i+1][j]-ss[i]*quq[i][j])/(rgm[j-1]*ssm[i])*byhs;
			dPs[i][j] = 0.0;
		}
	}
	/*cylindrical*/
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				Gnp[i][j] = -(qrogc[i][j]-rogc[i][j])/dt;
				Gnm[i][j] = (x[j+1]*qux[i][j+1]-x[j]*qux[i][j])*byhx/xm[j]+(quz[i+1][j]-quz[i][j])*byhz;
			}
		}
	}
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			dPc[i][j] = 0.0;
		}
	}


	do
	{
		m+=1;
		er = 0.0;
		flg = 0;
		/*spherical coordinate system*/
		for(i=0;i<NQ;i++)
		{
			for(j=1;j<Ng-1;j++)
			{
				if(i==0&&j==1)
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]*dPs[i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;	
				}
				else if(i==NQ-1&&j==1)
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]*dPs[i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
				}
				else if(i==0)
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]*dPs[i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;					
				}
				else if(i==NQ-1)
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]*dPs[i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
				}
				else if(j==1)
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]*dPs[i+1][j]+ss[i]*dPs[i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
				}
				else
				{
					PDIFF = Dnp[i][j]-Dnm[i][j]
						   +AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]*dPs[i+1][j]+ss[i]*dPs[i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
					PDIFF/= AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs;
				}
				PDIFF-= dPs[i][j];
				if(fabs(PDIFF)>er)er=fabs(PDIFF);
				dPs[i][j]+= OMG*PDIFF;
				if(m!=1&&flg==0)
				{
					if(PDIFF*dPDIFFs[i][j]<0.0)
					{
						flg = 1;
					}
				}
				dPDIFFs[i][j] = PDIFF;
			}
		}

		Search_SMAC_Cylindrical();
		/*cylindrical coordinate system*/
		for(i=0;i<NZ;i++)
		{
			for(j=0;j<NX;j++)
			{
				if(cflg[i][j]==0)
				{
					if(i==0&&j==0)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt;
						PDIFF/= dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else if(i==0&&j==NX-1)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt;
						PDIFF/= dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else if(i==NZ-1&&j==0)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt;
						PDIFF/= dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else if(i==NZ-1&&j==NX-1)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt;
						PDIFF/= dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else if(i==0)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt;
						PDIFF/= 2.0*dt*(byhx*byhx+0.5*byhz*byhz);
					}
					else if(j==0)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt;
						PDIFF/= 2.0*dt*(0.5*x[j+1]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else if(i==NZ-1)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt;
						PDIFF/= 2.0*dt*(byhx*byhx+0.5*byhz*byhz);
					}
					else if(j==NX-1)
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt;
						PDIFF/= 2.0*dt*(0.5*x[j]/xm[j]*byhx*byhx+byhz*byhz);
					}
					else
					{
						PDIFF = Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt;
						PDIFF/= 2.0*dt*(byhx*byhx+byhz*byhz);
					}
					PDIFF-= dPc[i][j];
					if(fabs(PDIFF)>er)er=fabs(PDIFF);
					dPc[i][j]+= OMG*PDIFF;
					if(m!=0&&flg==0)
					{
						if(dPDIFFc[i][j]*PDIFF<0.0)
						{
							flg = 1;
						}
					}
					dPDIFFc[i][j] = PDIFF;
				}
			}
		}
		Search_SMAC_Spherical();
		if(flg==0&&OMG<1.2)
		{
			OMG+= 0.01;
		}
		else if(flg == 1 && OMG>0.5)
		{
			OMG-= 0.01;
		}
	}while(er>erefT&&m<Mmax);

	if(m==Mmax)
	{
		printf("à≥óÕÇ™é˚ë©ÇµÇ‹ÇπÇÒÇ≈ÇµÇΩ\ter=%le\n",er);
		exit(1);
	}

	/*spherical*/
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Ng-1;j++)
		{
			rogs[i][j] = qrogs[i][j];
		}
	}
	AA = 1.0/log(rmax/Rsold);
	for(i=0;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			qur[i][j]+= -AAold*dt*(dPs[i][j+1]-dPs[i][j])*byhrg/rg2[j];
			ur[i][j] = 2.0*qur[i][j]/(rogs[i][j+1]+rogs[i][j]);
		}
	}
	for(i=1;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			quq[i][j]+= -(dPs[i][j]-dPs[i-1][j])*dt*byhs/rgm2[j-1];
			uq[i][j] = 2.0*quq[i][j]/(rogs[i][j]+rogs[i-1][j]);
		}
	}

	for(i=0;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			Ps[i][j]+= dPs[i][j];
			propertygs_cal(i,j);
		}
	}

	/*cylindrical*/
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				rogc[i][j] = qrogc[i][j];
			}
		}
	}
	for(i=1;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(czflg[i][j]==0)
			{
				quz[i][j]+= -(dPc[i][j]-dPc[i-1][j])*dt*byhz;
				uz[i][j] = 2.0*quz[i][j]/(rogc[i][j]+rogc[i-1][j]);
			}
		}
	}
	for(i=0;i<NZ;i++)
	{
		for(j=1;j<NX;j++)
		{
			if(cxflg[i][j]==0)
			{
				qux[i][j]+= -(dPc[i][j]-dPc[i][j-1])*dt*byhx;
				ux[i][j] = 2.0*qux[i][j]/(rogc[i][j]+rogc[i][j-1]);
			}
		}
	}
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				Pc[i][j]+=dPc[i][j];
				propertygc_cal(i,j);
			}
		}
	}

	return m;

}

double Ts_boundary(int i, double Ts)
{
	double Xs,X[spmax];
	double laht;
	int sp;
	double AA,BB;
	double Mm,rog;

	for(sp=0;sp<spmax;sp++)
	{
		X[sp] = Ys[i][0][sp] * Mms[i][0]/M[sp];
	}
	
	Xs = exp(Lf*M[nFUEL]/RR*(1.0/Tb0-1.0/Ts));
	AA = 1.0 - X[nFUEL];
	BB = 1.0 - Xs;
	X[nFUEL] = Xs;
	for(sp=1;sp<spmax-1;sp++)
	{
		X[sp] *= BB/AA;
	}
	AA = 1.0;
	for(sp=0;sp<spmax-1;sp++)
	{
		AA -= X[sp];
	}
	X[nN2] = AA;
	Mm = 0.0;
	for(sp=0;sp<spmax;sp++)
	{
		Mm += X[sp]*M[sp];
	}
	for(sp=0;sp<spmax-1;sp++)
	{
		Ys[i][0][sp] = X[sp]*M[sp]/Mm;
	}
	Ys[i][0][nN2] = N2sum(Ys[i][0]);
	Mms[i][0] = Mm;
	qrogs[i][0] = rogcal(Mms[i][0],Ts,Pa);
	propertygs_cal(i,0);
	laht = dHv_cal(Ts);
	AA = 1.0 / log(rmax/Rs);

	return qrogs[i][0]*Dgs[i][0][nFUEL]/(1.0-Ys[i][0][nFUEL])*AA*(8.0*Ys[i][0][nFUEL]-9.0*Ys[i][1][nFUEL]+Ys[i][2][nFUEL])*byhrg/3.0
			+(lgs[i][0]*AA*(8.0*Ts-9.0*Tgs[i][1]+Tgs[i][2])*byhrg/3.0+ll[i][Nl]*(8.0*Ts-9.0*Tls[i][Nl-1]+Tls[i][Nl-2])*byhrl/3.0)/laht;

}

void boundarys_cal()
{
	int i,sp,flg,k,g=0;
	double Ts1,Ts2,F1,F2,F3,Ts3,dTs,AA,dRs,Rsold2;
	double Ts0=290.0,laht;
	int Tsflg;

	do
	{
		g+=1;
		Rsold2 = Rs;
		for(i=0;i<NQ;i++)
		{
			/* èâä˙ílíTçı */
			Ts1 = Tls[i][Nl]+5.0;
			F1 = Ts_boundary(i,Ts1);
			Ts2 = Tls[i][Nl]-5.0;
			F2 = Ts_boundary(i,Ts2);
			if(F1*F2>0.0)
			{
				Ts1 = Ts0;
				flg = 0;
				k = 0;
				do
				{
					k += 1;
					F1 = Ts_boundary(i,Ts1);
					Ts2 = 0.5*(Ts1+Tb0);
					F2 = Ts_boundary(i,Ts2);
					if(F1*F2>0.0) Ts1 = Ts2;
					else flg = 1;
				}while(flg==0&&k<10);
				if(k==10)
				{
					printf("èâä˙ílÇ™å©Ç¬Ç©ÇËÇ‹ÇπÇÒÇ≈ÇµÇΩ\n");
					exit(1);
				}
			}
			/* 2ï™ñ@ */
			k = 0;
			Tsflg=0;
			do
			{
				k += 1;
				Ts3 = 0.5*(Ts1+Ts2);
				F3 = Ts_boundary(i,Ts3);
				if(fabs(F3)<1.0e-30)
				{
					Tsflg = 1;
				}
				if(F1*F3<0.0)
				{
					Ts2 = Ts3;
					F2 = F3;
				}
				else if(F2*F3<0.0)
				{
					Ts1 = Ts3;
					F1 = F3;
				}
				dTs = fabs(Ts1-Ts2);
			}while(k<Mmax&&dTs>erefZ&&Tsflg==0);
			if(k==Mmax)
			{
				printf("âÇ™å©Ç¬Ç©ÇËÇ‹ÇπÇÒÇ≈ÇµÇΩ\n");
				exit(1);
			}
			Tls[i][Nl] = Tgs[i][0] = 0.5*(Ts1+Ts2);
		}

		/* mdotïΩãœíl */
		AA = 1.0/log(rmax/Rs);
		for(i=0;i<NQ;i++)
		{
			laht = dHv_cal(Tls[i][Nl]);
			mdoti[i] = (lgs[i][0]*AA*(-8.0*Tgs[i][0]+9.0*Tgs[i][1]-Tgs[i][2])/3.0*byhrg
					   -ll[i][Nl]*(8.0*Tls[i][Nl]-9.0*Tls[i][Nl-1]+Tls[i][Nl-2])/3.0*byhrl)/(laht*Rs);
		}
		mdot = 0.0;
		for(i=0;i<NQ;i++)
		{
			mdot += mdoti[i]*(sc[i]-sc[i+1]);
		}
		mdot *= 0.5;
		Rdot = -mdot / rol;

		Rs = Rsold + Rdot*dt;

		dRs = fabs(Rsold2-Rs);

	}while(dRs>erefZ&&g<Mmax);
	if(g==Mmax)
	{
		printf("âtìHåaÇ™é˚ë©ÇµÇ‹ÇπÇÒÇ≈ÇµÇΩ\n");
		exit(1);
	}
	for(i=0;i<NQ;i++)
	{
		ur[i][0] = mdoti[i]/qrogs[i][0] + Rdot;
		qur[i][0] = qrogs[i][0]*ur[i][0];
		for(sp=1;sp<spmax-1;sp++)
		{
			Ys[i][0][sp] = Dgs[i][0][sp]*AA*(9.0*Ys[i][1][sp]-Ys[i][2][sp])/(3.0*hrg*Rs*mdoti[i]+8.0*Dgs[i][0][sp]*AA);
		}
		Ys[i][0][nN2] = N2sum(Ys[i][0]);
	}

}

void output(int k)
{
	char chark[20];
	char filename1[50];
	FILE *fp;
	int i,j,g;

	sprintf(chark,"%d",k);
	strcpy(filename1,filename);
	strcat(filename1,chark);
	strcat(filename1,".dat");

	OUP += 1;

	if((fp=fopen(filename1,"w"))==NULL)
	{
		printf("Can't open datafile\n");
		exit(1);
	}

	fprintf(fp,"TITLE = ""HEAT MASS TRANSFER""\n");
	fprintf(fp,"VARIABLES = ""z"", ""r"", ""T"",""rog"",""rogF"",""P"",""phi"",""UZ"",""UX""\n");
	fprintf(fp,"ZONE N=%d, E=%d, STRANDID=%d, SOLUTIONTIME=%le, DATAPACKING=BLOCK, ZONETYPE=FEQUADRILATERAL\n",NODE,ELM,OUP,tref);
	fprintf(fp,"VARLOCATION=([3-9]=CELLCENTERED)\n");
	fprintf(fp,"# z value\n");
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			fprintf(fp,"%le\t",z[i]);
		}
		fprintf(fp,"\n");
	}
	g = (NX+1)*(NZ+1)+1;
	fprintf(fp,"%le\n",ld);
	for(i=0;i<=NQ;i++)
	{
		for(j=1;j<=Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",rl[j]*sc[i]+ld);
			}
			else
			{
				fprintf(fp,"%le\t",rg[j-Nl]*sc[i]+ld);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# r value\n");
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			fprintf(fp,"%le\t",x[j]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"0.0\n");
	for(i=0;i<=NQ;i++)
	{
		for(j=1;j<=Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",rl[j]*ss[i]);
			}
			else
			{
				fprintf(fp,"%le\t",rg[j-Nl]*ss[i]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# T value (cell-centered)\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",Tgc[i][j]);
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",Tls[i][j]);
			}
			else
			{
				fprintf(fp,"%le\t",Tgs[i][j+1-Nl]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# rog value (cell-centered)\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",rogc[i][j]);
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",rogs[i][j+1-Nl]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# rogF value (cell-centered)\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",rogc[i][j]*Yc[i][j][0]);
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",rogs[i][j+1-Nl]*Ys[i][j+1-Nl][0]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# P value\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",Pc[i][j]);
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",Ps[i][j+1-Nl]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# phi value\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",STN*M[nO2]*Yc[i][j][nFUEL]/(M[nFUEL]*Yc[i][j][nO2]));
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",STN*M[nO2]*Ys[i][j+1-Nl][nFUEL]/(M[nFUEL]*Ys[i][j+1-Nl][nO2]));
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# UZ value\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",0.5*(uz[i+1][j]+uz[i][j]));
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",0.5*(ur[i][j+1-Nl]+ur[i][j-Nl])*scm[i]-0.5*(uq[i+1][j+1-Nl]+uq[i][j+1-Nl])*ssm[i]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# UX value\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%le\t",0.5*(ux[i][j+1]+ux[i][j+1]));
		}
		fprintf(fp,"\n");
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j<Nl)
			{
				fprintf(fp,"%le\t",0.0);
			}
			else
			{
				fprintf(fp,"%le\t",0.5*(ur[i][j+1-Nl]+ur[i][j-Nl])*ssm[i]+0.5*(uq[i+1][j+1-Nl]+uq[i][j+1-Nl])*scm[i]);
			}
		}
		fprintf(fp,"\n");
	}
	fprintf(fp,"# Connectivity list\n");
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			fprintf(fp,"%d\t%d\t%d\t%d\n",i*(NX+1)+j+1,i*(NX+1)+j+2,(i+1)*(NX+1)+j+2,(i+1)*(NX+1)+j+1);
		}
	}
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Nl+Ng;j++)
		{
			if(j==0)
			{
				fprintf(fp,"%d\t%d\t%d\t%d\n",g,g,g+i*(Nl+Ng)+1,g+(i+1)*(Nl+Ng)+1);
			}
			else
			{
				fprintf(fp,"%d\t%d\t%d\t%d\n",g+i*(Nl+Ng)+j,g+i*(Nl+Ng)+j+1,g+(i+1)*(Nl+Ng)+j+1,g+(i+1)*(Nl+Ng)+j);
			}
		}
	}
	fclose(fp);

}

void dt_cal(int k)
{
	int i,j;
	double u1,u2;
	double CRN1,CRN2;
	double CRN;
	double dtold=dt;

	CRN=0.0;
	for(i=0;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			u1 = 0.5*(ur[i][j]+ur[i][j-1]);
			u2 = 0.5*(uq[i+1][j]+uq[i][j]);
			CRN1 = u1*dt/(rg[j]-rg[j-1]);
			CRN2 = u2*dt/(rgm[j-1])*byhs;
			if(fabs(CRN1)>CRN)CRN=fabs(CRN1);
			if(fabs(CRN2)>CRN)CRN=fabs(CRN2);
		}
	}
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				u1 = 0.5*(uz[i+1][j]+uz[i][j]);
				u2 = 0.5*(ux[i][j+1]+ux[i][j]);
				CRN1 = u1*dt*byhz;
				CRN2 = u2*dt*byhx;
				if(fabs(CRN1)>CRN)CRN=fabs(CRN1);
				if(fabs(CRN2)>CRN)CRN=fabs(CRN2);
			}
		}
	}
	
	dt*= 0.01/CRN;
	if(k<dtchange){
	if(dt>dtdefini)
	{
		dt = dtdefini;
	}
	}
	else{
		if(dt>dtdef)
	{
		dt = dtdef;
	}
	}
}

void iteration()
{
	int i,j,k,m;
	FILE *fp;
	double treq;

	if((fp=fopen("Rs.txt","w"))==NULL)
	{
		printf("Can't open datafile\n");
		exit(1);
	}

	k = 0;
	output(k);
	treq = tfreq;
	do
	{
		Rsold = Rs;
		k+=1;
		tref+= dt;

		/* liquid phase */
		Tl_cal();

		/* spherical coordinate system gas phase */
		Ys_cal();
		Tgs_cal();
		ur_cal();
		uq_cal();

		/* cylindrical coordinate system gas phase */
		Yc_cal();
		Tgc_cal();
		ux_cal();
		uz_cal();

		density_cal();

		boundarys_cal();

		rlcal();
		rgcal();

		if(k>kchange){m = SMAC_method();}
		else{
			m = CG_method();
		}


		Search();

		dt_cal(k);

		if(k%1==0)
		{
			printf("t=%le\tm=%d\tRs=%le\tdt=%le\n",tref,m,Rs*1000.0,dt);
		}

		if(tref>treq)
		//if(k%1==0)
		{
			treq += tfreq;
			output(k);
	
			output2(k);

			output3(k);
			output4(k);
		
			fprintf(fp,"%d\t%le\t%le\t%le\n",k,tref,Rs*1000.0,Rdot*1000.0);
			fflush(fp);
		}
			
	if(tref>tmax)
	{ exit(1);}

	}while(k<Kmax&&Rs>Rsmin);
		
	fclose(fp);														

}

void main(int argc,char *argv[])
{
	static int i,sp,g;

	if (argc<4){
		printf("\n example : %s [Pa(bar)] [Ta(K)] [l(mm)]",argv[0]);
		printf("\n  *******   option   *******");
		printf("\n D0 [mm] Dstop[mm]      : d [] []       %1.0lf,%1.0lf",2.0*Rs0,2.0*Rsmin);
		printf("\n Nl Ng NQ rmax          : g [] [] [] [] %d,%d,%d,%4.1lf",Nl,Ng,NQ,rmax);
		printf("\n NX NZ                  : f [] []       %d,%d",NX,NZ);
		printf("\n droplet temperature    : l []          %1.0lf",Tl0);
		printf("\n calcstop (time, sec)   : e []          %2.1lf (off)",timestop);
		printf("\n");
		exit(1);
	}
	Pa=atol(argv[1])*1.0e5;
	Ta=atol(argv[2]);
	ld=atol(argv[3])*0.5e-3;

	for (i=3;i<argc;i++){
		switch (*argv[i]){
			case 'd' : case 'D' :{
				Rs0		=atof(argv[++i])*0.5*1.0e-3;
				Rsmin	=atof(argv[++i])*0.5*1.0e-3;
			} break;
			case 'g' : case 'G' :{
				Nl=atoi(argv[++i]);
				Ng=atoi(argv[++i]);
				NQ=atoi(argv[++i]);
				rmax=atof(argv[++i])*1.0e-3;
			} break;
			case 'f' : case 'F' :{
				NX=atoi(argv[++i]);
				NZ=atoi(argv[++i]);
			} break;
			case 'l' : case 'L' :	Tl0=atof(argv[++i]);		break;
			case 'e' : case 'E' :	timestop=atof(argv[++i]);	break;
		}
	}
	NODE = (NX+1)*(NZ+1) + (NQ+1)*(Nl+Ng)+1;
	ELM = NX*NZ + NQ*(Nl+Ng);


	alm();
	pint();
	filedef();

	printf("start iteration!!\n");
	iteration();

}


void Search_CG()
{
	int i,j,ii,jj,k,flg;
	double r0,s0,u1,u2,aa,bb;
	double *b;
	double v1,v2,v3,v4;

	b = dvec(3);

	/* from spherical to cylindrical */
		/* z Vector */
	for(i=1;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(czflg[i][j]==2)
			{
				r0 = sqrt(xm[j]*xm[j]+(z[i]-ld)*(z[i]-ld));
					s0 = acos((z[i]-ld)/r0);
					flg = 0;
					k = 0;
					do
					{
						if(rgm[k]>r0)
						{
							flg = 1;
							jj = k;
						}
						k+=1;
					}while(flg==0);
					flg = 0;
					k = 0;
					do
					{
						if(((double)k+0.5)*hs>s0)
						{
							flg = 1;
							ii = k;
						}
						k+=1;
					}while(flg==0&&k<NQ);
				if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.5)*hs;

					   

					if(aa>r0)
					{    u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v1 = u2*scm[ii]-u1*ssm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
						v3 = u2*sc[NQ]-u1*ss[NQ];
							
						
						if(bb>s0){u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
							v2 = u2*scm[ii]-u1*ssm[ii];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}

					    else{u1 = 0.0;
							u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj+1]-qur[ii-1][jj+1]);
							v4 = u2*sc[NQ]-u1*ss[NQ];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v2 = u2*scm[ii]-u1*ssm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj+1]-qur[ii-1][jj+1]+9.0*qur[ii][jj]-qur[ii-1][jj]);
						v4 = u2*sc[NQ]-u1*ss[NQ];
												
						if(bb>s0){u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
							v1 = u2*scm[ii]-u1*ssm[ii];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					     else{u1 = 0.0;
							u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
							v3 = u2*sc[NQ]-u1*ss[NQ];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj];}}
				}
				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.5*hs;
					if(aa>r0)
					{  
						
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
						v1 = u2*sc[0]-u1*ss[0];
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 = u2*scm[ii]-u1*ssm[ii];
												
						if(bb>s0){u1 = 0.0;
							u2 = 0.0625*(9.0*qur[0][jj+1]-qur[1][jj+1]+9.0*qur[0][jj]-qur[1][jj]);
							v2 = u2*sc[0]-u1*ss[0];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj-1];}
					    else{u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*0.5*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj+1]-qur[1][jj+1]);
						v2 = u2*sc[0]-u1*ss[0];
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v4 = u2*scm[ii]-u1*ssm[ii];
						if(bb>s0){
							u1 = 0.0;
							u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
							v1 = u2*sc[0]-u1*ss[0];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*0.5*hs - b[1]*rgm[jj];}
					     else{u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
							v3 = u2*scm[ii]-u1*ssm[ii];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*0.5*hs - b[1]*rgm[jj];}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa>r0)
					{   u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
						u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
						v1 = u2*scm[ii-1]-u1*ssm[ii-1];
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj-1]+qur[ii][jj]);
						v3 = u2*scm[ii]-u1*ssm[ii];
						
						
						if(bb>s0){u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
							u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
							v2 = u2*scm[ii-1]-u1*ssm[ii-1];
							b[0] = (v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj-1];}
					    else{u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
							b[0] = (v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
							u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
							v2 = u2*scm[ii-1]-u1*ssm[ii-1];
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 = u2*scm[ii]-u1*ssm[ii];
						
						
						if(bb>s0){ u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
						u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
						v1 = u2*scm[ii-1]-u1*ssm[ii-1];
							b[0] = (v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];}
					     else{
							u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						  u2 = 0.5*(qur[ii][jj-1]+qur[ii][jj]);
						 v3 = u2*scm[ii]-u1*ssm[ii];
							b[0] = (v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}}
				}
					quz[i][j] = b[0]*s0+b[1]*r0+b[2];

				}
			}
		}
		/* x Vector */
	for(i=0;i<NZ;i++)
	{
		for(j=1;j<NX;j++)
		{
			
			if(cxflg[i][j]==2)
			{       r0 = sqrt(x[j]*x[j]+(zm[i]-ld)*(zm[i]-ld));
			 		s0 = acos((zm[i]-ld)/r0);
					flg = 0;
					k = 0;
					do
					{
						if(rgm[k]>r0)
						{
							flg = 1;
							jj = k;
						}
						k+=1;
					}while(flg==0);
					flg = 0;
					k = 0;
					do
					{
						if(((double)k+0.5)*hs>s0)
						{
							flg = 1;
							ii = k;
						}
						k+=1;
					}while(flg==0&&k<NQ);


					if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.5)*hs;

					   

					if(aa>r0)
					{    u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v1 =  u2*ssm[ii]+u1*scm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
						v3 =  u2*ss[NQ]+u1*sc[NQ];
							
						
						if(bb>s0){u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
							v2 =  u2*ssm[ii]+u1*scm[ii];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}

					    else{u1 = 0.0;
							u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj+1]-qur[ii-1][jj+1]);
							v4 =  u2*ss[NQ]+u1*sc[NQ];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v2 =  u2*ssm[ii]+u1*scm[ii];
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[ii][jj+1]-qur[ii-1][jj+1]+9.0*qur[ii][jj]-qur[ii-1][jj]);
						v4 =  u2*ss[NQ]+u1*sc[NQ];
												
						if(bb>s0){u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
							v1 =  u2*ssm[ii]+u1*scm[ii];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					     else{u1 = 0.0;
							u2 = 0.0625*(9.0*qur[ii][jj]-qur[ii-1][jj]+9.0*qur[ii][jj-1]-qur[ii-1][jj-1]);
							v3 =  u2*ss[NQ]+u1*sc[NQ];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj];}}
				}
				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.5*hs;
					if(aa>r0)
					{  
						
						u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
						v1 =  u2*ss[ii]+u1*sc[ii];
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
						v3 =  u2*ssm[ii]+u1*scm[ii];
												
						if(bb>s0){u1 = 0.0;
							u2 = 0.0625*(9.0*qur[0][jj+1]-qur[1][jj+1]+9.0*qur[0][jj]-qur[1][jj]);
							v2 =  u2*ss[ii]+u1*sc[ii];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*(double)(ii+0.5)*hs - b[1]*rgm[jj-1];}
					    else{u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 = u2*ssm[ii]+u1*scm[ii];
							b[0] = 2.0*(v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*0.5*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.0;
						u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj+1]-qur[1][jj+1]);
						v2 =  u2*ss[ii]+u1*sc[ii];
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
						u2 = 0.5*(qur[ii][jj+1]+qur[ii][jj]);
						v4 = u2*ssm[ii]+u1*scm[ii];
						if(bb>s0){
							u1 = 0.0;
							u2 = 0.0625*(9.0*qur[0][jj]-qur[1][jj]+9.0*qur[0][jj-1]-qur[1][jj-1]);
							v1 = u2*ss[ii]+u1*sc[ii];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*0.5*hs - b[1]*rgm[jj];}
					     else{u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj-1]);
							v3 =  u2*ssm[ii]+u1*scm[ii];
							b[0] = 2.0*(v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*0.5*hs - b[1]*rgm[jj];}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa>r0)
					{   u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
						u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
						v1 =  u2*ssm[ii-1]+u1*scm[ii-1];
						u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						u2 = 0.5*(qur[ii][jj-1]+qur[ii][jj]);
						v3 =  u2*ssm[ii]+u1*scm[ii];
						
						
						if(bb>s0){u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
							u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
							v2 =  u2*ssm[ii-1]+u1*scm[ii-1];
							b[0] = (v3-v1)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v1 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj-1];}
					    else{u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 =  u2*ssm[ii]+u1*scm[ii];
							b[0] = (v3-v1)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v3 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}}
					else
					{   u1 = 0.5*(quq[ii][jj+1]+quq[ii-1][jj+1]);
							u2 = 0.5*(qur[ii-1][jj+1]+qur[ii-1][jj]);
							v2 =  u2*ssm[ii-1]+u1*scm[ii-1];
						u1 = 0.5*(quq[ii+1][jj+1]+quq[ii][jj+1]);
							u2 = 0.5*(qur[ii][jj]+qur[ii][jj+1]);
							v4 =  u2*ssm[ii]+u1*scm[ii];
						
						
						if(bb>s0){ u1 = 0.5*(quq[ii][jj]+quq[ii-1][jj]);
						u2 = 0.5*(qur[ii-1][jj]+qur[ii-1][jj-1]);
						v1 =  u2*ssm[ii-1]+u1*scm[ii-1];
							b[0] = (v4-v2)*byhs;
							b[1] = (v2-v1)/(rgm[jj]-rgm[jj-1]);
							b[2] = v2 - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];}
					     else{
							u1 = 0.5*(quq[ii+1][jj]+quq[ii][jj]);
						  u2 = 0.5*(qur[ii][jj-1]+qur[ii][jj]);
						 v3 =  u2*ssm[ii]+u1*scm[ii];
							b[0] = (v4-v2)*byhs;
							b[1] = (v4-v3)/(rgm[jj]-rgm[jj-1]);
							b[2] = v4 - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}}
				}


					

					qux[i][j] = b[0]*s0+b[1]*r0+b[2];
				}
			}
		
	}
	
	

	freedvec(b,3);
	
}

void Search_CG_Spherical()
{
	int i,j,k,flg,ii,jj;
	double *b;
	double x0,z0,aa,bb;

	b = dvec(3);

	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{
			x0 = rgm[i-1]*ssm[j];
			z0 = rgm[i-1]*scm[j] + ld;
			sPflgx[j][i]=x0;
			sPflgz[j][i]=z0;
			flg = 0;
			k = 0;
			do
			{
				if(z0<zm[k])
				{
					flg = 1;
					ii = k;
				}
				k += 1;
			}while(flg==0&&k<=NZ);
			flg = 0;
			k = 0;
			do
			{
				if(x0<xm[k])
				{
					flg = 1;
					jj = k;
				}
				k+=1;
			}while(flg==0&&k<=NX);
			if(jj==0)
			{
				aa = 0.5*xm[jj];
				bb = 0.5*(zm[ii]+zm[ii-1]);
				if(x0>aa)
				{ if(z0>bb){sPflg[j][i]=0;}
				  else{sPflg[j][i]=1;}}
				else
				{  if(z0>bb){sPflg[j][i]=2;}
				else{sPflg[j][i]=3;}}
			}
			else
			{
				aa = (xm[jj-1]+xm[jj])/2;
				bb = (zm[ii-1]+zm[ii])/2;
				//printf("%f,%f\n",bb,z0);
				if(x0>aa)
				{ if(z0>bb){sPflg[j][i]=4;}
				else{//printf("a\n");
				sPflg[j][i]=5;}}
				else
				{  if(z0>bb){sPflg[j][i]=6;}
				else{//printf("a\n");
				sPflg[j][i]=7;}}
			}
			sPflgi[j][i]=ii;
			sPflgj[j][i]=jj;
			
			
		}
	}

	freedvec(b,3);
}

void Search_CG_Cylindrical()
{
	int i,j,k,flg,ii,jj,g;
	
	double r0,s0,aa,bb;

	

	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==2)
			{
				r0 = sqrt(xm[j]*xm[j]+(zm[i]-ld)*(zm[i]-ld));
				s0 = acos((zm[i]-ld)/r0);
				cPflgr[i][j]=r0;
				cPflgs[i][j]=s0;
				flg = 0;
				k = 0;
				do
				{
					if(rgm[k]>r0)
					{
						flg = 1;
						jj = k;
					}
					k+=1;
				}while(flg==0);
				flg = 0;
				k = 0;
				do
				{
					if(((double)k+0.5)*hs>s0)
					{
						flg = 1;
						ii = k;
					}
					k+=1;
				}while(flg==0&&k<NQ);
				
				

				if(flg==0)
				{
					ii = NQ-1;
					aa = (rgm[jj-1]+rgm[jj])/2;
					bb = (ii+0.75)*hs;
					if(aa<r0)
					{   if(bb<s0){cPflg[i][j]=0;}
					    else{cPflg[i][j]=1;}}
					else
					{  if(bb<s0){cPflg[i][j]=2;}
					     else{cPflg[i][j]=3;}}
				}
				else if(ii==0)
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb = 0.25*hs;
					if(aa<r0)
					{   if(bb<s0){cPflg[i][j]=4;}
					    else{cPflg[i][j]=5;}}
					else
					{  if(bb<s0){cPflg[i][j]=6;}
					     else{cPflg[i][j]=7;}}
				}
				else
				{
					aa =  (rgm[jj-1]+rgm[jj])/2;
					bb =ii*hs;
					if(aa<r0)
					{   if(bb<s0){cPflg[i][j]=8;}
					    else{cPflg[i][j]=9;}}
					else
					{  if(bb<s0){cPflg[i][j]=10;}
					     else{cPflg[i][j]=11;}}
				}
				  cPflgi[i][j]=ii;
				  cPflgj[i][j]=jj;
			}
		}
	}
	
}

int CG_method()
{
	int i,j,ii,jj,k,m=0,mm=0,flg;
	double er,er1,er2,er3;
	double RdotbyRs,AA,AAold,PDIFF;
	double *b;

	AA = 1.0 / log(rmax/Rs);
	AAold = 1.0 / log(rmax/Rsold);
	RdotbyRs = Rdot/Rs;
	b = dvec(3);


	Search_CG();
    Search_CG_Cylindrical();
	Search_CG_Spherical();

	/*spherical*/
	for(i=0;i<NQ;i++)
	{
		Dnp[i][1] = -(qrogs[i][1]-rogs[i][1])/dt+0.5*RdotbyRs*AA*((1.0-eta[1])*(qrogs[i][1+1]-qrogs[i][1])
																 +(1.0-eta[0])*(-8.0*qrogs[i][0]+9.0*qrogs[i][1]-qrogs[i][2])/3.0)*byhrg;
		Dnm[i][1] = AA*(rg[1]*rg[1]*qur[i][1]-rg[0]*rg[0]*qur[i][0])/(rgm[0]*rgm[0]*rgm[0])*byhrg
				   +(ss[i+1]*quq[i+1][1]-ss[i]*quq[i][1])/(rgm[0]*ssm[i])*byhs;
		dPs[i][1] = 0.0;

		for(j=2;j<=Ng;j++)
		{
			Dnp[i][j] = -(qrogs[i][j]-rogs[i][j])/dt+0.5*RdotbyRs*AA*((1.0-eta[j])*(qrogs[i][j+1]-qrogs[i][j])
																	 +(1.0-eta[j-1])*(qrogs[i][j]-qrogs[i][j-1]))*byhrg;
			Dnm[i][j] = AA*(rg[j]*rg[j]*qur[i][j]-rg[j-1]*rg[j-1]*qur[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*byhrg
					   +(ss[i+1]*quq[i+1][j]-ss[i]*quq[i][j])/(rgm[j-1]*ssm[i])*byhs;
			dPs[i][j] = 0.0;
		}
	}
	
	/*cylindrical*/
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				Gnp[i][j] = -(qrogc[i][j]-rogc[i][j])/dt;
				Gnm[i][j] = (x[j+1]*qux[i][j+1]-x[j]*qux[i][j])*byhx/xm[j]+(quz[i+1][j]-quz[i][j])*byhz;
				
			}
		}
	}
	for(i=0;i<=NZ;i++)
	{
		for(j=0;j<=NX;j++)
		{
			dPc[i][j] = 0.0;
			Es[i][j]=0;
		}
	}

   for(i=0;i<NQ;i++)
	{   
		for(j=0;j<=Ng;j++)
		{
			Es[i+NZ+1][j]=dPs[i][j];}}
	
	
   
  	for(i=0;i<NZ;i++)
		{
			for(j=0;j<NX;j++)
			{
				if(cflg[i][j]==0)
				{
					if(i==0&&j==0)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j] = Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt
							   -dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
					}
					else if(i==0&&j==NX-1)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt
						-dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
							 
					}
					else if(i==NZ-1&&j==0)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt
						-dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
							  
					}
					else if(i==NZ-1&&j==NX-1)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt
						-dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
					}
					else if(i==0)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +dPc[i+1][j]*byhz*byhz*dt
								-2.0*dt*(byhx*byhx+0.5*byhz*byhz)*dPc[i][j];
					}
					else if(j==0)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]=Gnp[i][j]-Gnm[i][j]
							   +x[j+1]*dPc[i][j+1]*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt
						-2.0*dt*(0.5*x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
					}
					else if(i==NZ-1)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +dPc[i-1][j]*byhz*byhz*dt
								-2.0*dt*(byhx*byhx+0.5*byhz*byhz)*dPc[i][j];
					}
					else if(j==NX-1)
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +x[j]*dPc[i][j-1]*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt
						       -2.0*dt*(0.5*x[j]/xm[j]*byhx*byhx+byhz*byhz)*dPc[i][j];
					}
					else
					{
						Fs[i][j]=Is[i][j]=IIs[i][j]= Gnp[i][j]-Gnm[i][j]
							   +(x[j+1]*dPc[i][j+1]+x[j]*dPc[i][j-1])*byhx*byhx*dt/xm[j]
							   +(dPc[i+1][j]+dPc[i-1][j])*byhz*byhz*dt
								-2.0*dt*(byhx*byhx+byhz*byhz)*dPc[i][j];
					}
					
				}
				else if(cflg[i][j]==2){

					ii=cPflgi[i][j];
					jj=cPflgj[i][j];

					if(cPflg[i][j]==0){b[0] = 0.25*(dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = 0.125*(9.0*(dPs[ii][jj+1]-dPs[ii][jj])-(dPs[ii-1][jj+1]-dPs[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==1){b[0] = 0.25*(dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==2){b[0] = 0.25*(dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = 0.125*(9.0*(dPs[ii][jj+1]-dPs[ii][jj])-(dPs[ii-1][jj+1]-dPs[ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==3){b[0] = 0.25*(dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==4){b[0] = 0.25*(dPs[1][jj+1]-dPs[0][jj+1])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==5){b[0] = 0.25*(dPs[1][jj+1]-dPs[0][jj+1])*byhs;
						b[1] =  0.125*(9.0*(dPs[0][jj+1]-dPs[0][jj])-(dPs[1][jj+1]-dPs[1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==6){b[0] = 0.25*(dPs[1][jj]-dPs[0][jj])*byhs;
						b[1] =(dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==7){b[0] = 0.25*(dPs[1][jj]-dPs[0][jj])*byhs;
						b[1] = 0.125*(9.0*(dPs[0][jj+1]-dPs[0][jj])-(dPs[1][jj+1]-dPs[1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==8){b[0] = (dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==9){b[0] = (dPs[ii][jj+1]-dPs[ii-1][jj+1])*byhs;
						b[1] = (dPs[ii-1][jj+1]-dPs[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii-1][jj+1] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==10){b[0] = (dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = (dPs[ii][jj+1]-dPs[ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else{b[0] = (dPs[ii][jj]-dPs[ii-1][jj])*byhs;
						b[1] = (dPs[ii-1][jj+1]-dPs[ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = dPs[ii-1][jj] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj-1];}

					Fs[i][j]=Is[i][j]=IIs[i][j]= b[0]*cPflgs[i][j]+b[1]*cPflgr[i][j]+b[2]-dPc[i][j];}



				
				
			}
		}
	
			//Fs[i][j]=Dnp[i][j]-Dnm[i][j];
			//Is[i][j]=Dnp[i][j]-Dnm[i][j];
			//if(Fs[i][j]!=0){printf("f%f\n",Fs[i][j]);}
			//if(Is[i][j]!=0){printf("i%f\n",Is[i][j]);}


		
	for(i=0;i<NQ;i++)
		{
			for(j=1;j<Ng-1;j++)
			{
				if(i==0&&j==1)
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j]+AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]*dPs[i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
							-(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						  +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];	
				}
				else if(i==NQ-1&&j==1)
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j]+AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]*dPs[i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					-( AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];
				}
				else if(i==0)
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j]+AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]*dPs[i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
				    -(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];					
				}
				else if(i==NQ-1)
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j] +AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]*dPs[i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					-(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];
				}
				else if(j==1)
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j] +AAold*AA*rg[j]*rg[j]/rg2[j]*dPs[i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]*dPs[i+1][j]+ss[i]*dPs[i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					-(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];
				}
				else
				{
					Fs[i+NZ+1][j]=Is[i+NZ+1][j]=IIs[i+NZ+1][j]=Dnp[i][j]-Dnm[i][j] +AAold*AA*(rg[j]*rg[j]/rg2[j]*dPs[i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*dPs[i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]*dPs[i+1][j]+ss[i]*dPs[i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					        -(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*dPs[i][j];
				}

				
			//printf("%f,",Fs[i][j]);
				//if(Gs[i][j]!=0){printf("%f,",Gs[i][j]);}
			}}
	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{   ii=sPflgi[j][i];
	 	    jj=sPflgj[j][i];			
			
			if(sPflg[j][i]==0){b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = 0.25*(dPc[ii][1]-dPc[ii][0])*byhx;
					b[2] = dPc[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==1){b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = 0.25*(dPc[ii-1][1]-dPc[ii-1][0])*byhx;
					b[2] = dPc[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==2){b[0] = 0.125*(9.0*(dPc[ii][0]-dPc[ii-1][0])-(dPc[ii][1]-dPc[ii-1][1]))*byhz;
					b[1] = 0.25*(dPc[ii][1]-dPc[ii][0])*byhx;
					b[2] = dPc[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==3){b[0] = 0.125*(9.0*(dPc[ii][0]-dPc[ii-1][0])-(dPc[ii][1]-dPc[ii-1][1]))*byhz;
					b[1] = 0.25*(dPc[ii-1][1]-dPc[ii-1][0])*byhx;
					b[2] = dPc[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==4){b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = (dPc[ii][jj]-dPc[ii][jj-1])*byhx;
					b[2] = dPc[ii][jj] - b[0]*zm[ii] - b[1]*xm[jj];}
			else if(sPflg[j][i]==5){b[0] = (dPc[ii][jj]-dPc[ii-1][jj])*byhz;
					b[1] = (dPc[ii-1][jj]-dPc[ii-1][jj-1])*byhx;
					b[2] = dPc[ii-1][jj] - b[0]*zm[ii-1] - b[1]*xm[jj];}
			else if(sPflg[j][i]==6){b[0] = (dPc[ii][jj-1]-dPc[ii-1][jj-1])*byhz;
					b[1] = (dPc[ii][jj]-dPc[ii][jj-1])*byhx;
					b[2] = dPc[ii][jj-1] - b[0]*zm[ii] - b[1]*xm[jj-1];}
			else{b[0] = (dPc[ii][jj-1]-dPc[ii-1][jj-1])*byhz;
					b[1] = (dPc[ii-1][jj]-dPc[ii-1][jj-1])*byhx;
					b[2] = dPc[ii-1][jj-1] - b[0]*zm[ii-1] - b[1]*xm[jj-1];}

			Fs[j+NZ+1][i]=Is[j+NZ+1][i]=IIs[j+NZ+1][i]= b[0]*sPflgz[j][i]+b[1]*sPflgx[j][i]+b[2]-dPs[j][i];}}
	//for(i=0;i<NQ;i++){
//	Fs[NZ+1][0]=Is[NZ+1][0]=IIs[NZ+1][0]=-dPs[0][1];
//}
	
   /*for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			
			Fc[i][j]=Gnp[i][j]-Gnm[i][j];
			Ic[i][j]=Gnp[i][j]-Gnm[i][j];
		}
	}
  */

    Ks=0;
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			Ks+=Is[i][j]*IIs[i][j];
			
		}
	}
for(i=0;i<NQ;i++)
	{
		for(j=1;j<=Ng;j++)
		{
			Ks+=Is[i+NZ+1][j]*IIs[i+NZ+1][j];
			
		}
	}
//Ks+=Is[NZ+1][0]*IIs[NZ+1][0];
er2=Ks;
er=1;

//printf("%f",Ks);
 /*for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			Kc+=Ic[i][j]*Ic[i][j];
		}
	}*/


	do
	{ 
		/*
		if(mm!=0){
			for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Is[i][j]=Ms[i][j]-Ns*Ws[i][j];
					Ks+=Is[i][j]*IIs[i][j];
				}
			}
		for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{   Is[i+NZ+1][j]=Ms[i+NZ+1][j]-Ns*Ws[i+NZ+1][j];
					Ks+=Is[i+NZ+1][j]*IIs[i+NZ+1][j];
				}
			}


				for(i=0;i<NZ;i++)
					{
						for(j=0;j<NX;j++)
						{
							Fs[i][j]=Is[i][j]+Hs/Ns*Ks/Js*(Fs[i][j]-Ns*Gs[i][j]);
						}
					}
				for(i=0;i<NQ;i++)
					{
						for(j=1;j<=Ng;j++)
						{
							Fs[i+NZ+1][j]=Is[i+NZ+1][j]+Hs/Ns*Ks/Js*(Fs[i+NZ+1][j]-Ns*Gs[i+NZ+1][j]);
						}
					}*/

			          //printf("b%f,%f,%f,%f,\n",Is[0][1],Ns,Ks/Js,Fs[0][1]-Ns*Gs[0][1]);
				/*for(i=0;i<NZ;i++)
				{
					for(j=0;j<NX;j++)
					{
						Ic[i][j]-=Hc*Gc[i][j];
						Kc+=Ic[i][j]*Ic[i][j];}}
				for(i=0;i<NZ;i++)
				{
					for(j=0;j<NX;j++)
					{
						Fc[i][j]=Ic[i][j]+Kc*Fc[i][j];}}}*/

        mm++;
		flg=0;
		er3 = 0.0;
		Js=Ks;
		Ls=0;
		Ks=0;
		Us=0;
		Vs=0;
		//Jc=Kc;
		//Lc=0;
		//Kc=0; 
		
		
		
		for(i=0;i<NZ;i++)
		{
			for(j=0;j<NX;j++)
			{
				if(cflg[i][j]==0)
				{
					if(i==0&&j==0)
					{
						Gs[i][j]=-x[j+1]*Fs[i][j+1]*byhx*byhx*dt/xm[j]
							   -Fs[i+1][j]*byhz*byhz*dt
							   +dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else if(i==0&&j==NX-1)
					{
						Gs[i][j]=-x[j]*Fs[i][j-1]*byhx*byhx*dt/xm[j]
							   -Fs[i+1][j]*byhz*byhz*dt
						+dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else if(i==NZ-1&&j==0)
					{
						Gs[i][j]=-x[j+1]*Fs[i][j+1]*byhx*byhx*dt/xm[j]
							   -Fs[i-1][j]*byhz*byhz*dt
						+dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else if(i==NZ-1&&j==NX-1)
					{
						Gs[i][j]= -x[j]*Fs[i][j-1]*byhx*byhx*dt/xm[j]
							   -Fs[i-1][j]*byhz*byhz*dt
						+dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else if(i==0)
					{
						Gs[i][j]= -(x[j+1]*Fs[i][j+1]+x[j]*Fs[i][j-1])*byhx*byhx*dt/xm[j]
							   -Fs[i+1][j]*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+0.5*byhz*byhz)*Fs[i][j];
					}
					else if(j==0)
					{
						Gs[i][j]=-x[j+1]*Fs[i][j+1]*byhx*byhx*dt/xm[j]
							   -(Fs[i+1][j]+Fs[i-1][j])*byhz*byhz*dt
						+2.0*dt*(0.5*x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else if(i==NZ-1)
					{
						Gs[i][j]=-(x[j+1]*Fs[i][j+1]+x[j]*Fs[i][j-1])*byhx*byhx*dt/xm[j]
							   -Fs[i-1][j]*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+0.5*byhz*byhz)*Fs[i][j];
					}
					else if(j==NX-1)
					{
						Gs[i][j]= -x[j]*Fs[i][j-1]*byhx*byhx*dt/xm[j]
							   -(Fs[i+1][j]+Fs[i-1][j])*byhz*byhz*dt
						       +2.0*dt*(0.5*x[j]/xm[j]*byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					else
					{
						Gs[i][j]= -(x[j+1]*Fs[i][j+1]+x[j]*Fs[i][j-1])*byhx*byhx*dt/xm[j]
							   -(Fs[i+1][j]+Fs[i-1][j])*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+byhz*byhz)*Fs[i][j];
					}
					
				}
				else if(cflg[i][j]==2){

					ii=cPflgi[i][j];
					jj=cPflgj[i][j];

					
					if(cPflg[i][j]==0){b[0] = 0.25*(Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii-1][jj+1])*byhs;
						b[1] = 0.125*(9.0*(Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])-(Fs[NZ+1+ii-1][jj+1]-Fs[NZ+1+ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==1){b[0] = 0.25*(Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==2){b[0] = 0.25*(Fs[NZ+1+ii][jj]-Fs[NZ+1+ii-1][jj])*byhs;
						b[1] = 0.125*(9.0*(Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])-(Fs[NZ+1+ii-1][jj+1]-Fs[NZ+1+ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==3){b[0] = 0.25*(Fs[NZ+1+ii][jj]-Fs[NZ+1+ii-1][jj])*byhs;
						b[1] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==4){b[0] = 0.25*(Fs[NZ+1+1][jj+1]-Fs[NZ+1+0][jj+1])*byhs;
						b[1] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==5){b[0] = 0.25*(Fs[NZ+1+1][jj+1]-Fs[NZ+1+0][jj+1])*byhs;
						b[1] =  0.125*(9.0*(Fs[NZ+1+0][jj+1]-Fs[NZ+1+0][jj])-(Fs[NZ+1+1][jj+1]-Fs[NZ+1+1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==6){b[0] = 0.25*(Fs[NZ+1+1][jj]-Fs[NZ+1+0][jj])*byhs;
						b[1] =(Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==7){b[0] = 0.25*(Fs[NZ+1+1][jj]-Fs[NZ+1+0][jj])*byhs;
						b[1] = 0.125*(9.0*(Fs[NZ+1+0][jj+1]-Fs[NZ+1+0][jj])-(Fs[NZ+1+1][jj+1]-Fs[NZ+1+1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==8){b[0] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==9){b[0] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Fs[NZ+1+ii-1][jj+1]-Fs[NZ+1+ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii-1][jj+1] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==10){b[0] = (Fs[NZ+1+ii][jj]-Fs[NZ+1+ii-1][jj])*byhs;
						b[1] = (Fs[NZ+1+ii][jj+1]-Fs[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else{b[0] = (Fs[NZ+1+ii][jj]-Fs[NZ+1+ii-1][jj])*byhs;
						b[1] = (Fs[NZ+1+ii-1][jj+1]-Fs[NZ+1+ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Fs[NZ+1+ii-1][jj] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj-1];}

					Gs[i][j]= -b[0]*cPflgs[i][j]-b[1]*cPflgr[i][j]-b[2]+Fs[i][j];}
				else{Gs[i][j]=0;}				
				
			}
		}
	
			//Fs[i][j]=Dnp[i][j]-Dnm[i][j];
			//Is[i][j]=Dnp[i][j]-Dnm[i][j];
			//if(Fs[i][j]!=0){printf("f%f\n",Fs[i][j]);}
			//if(Is[i][j]!=0){printf("i%f\n",Is[i][j]);}


		
	for(i=0;i<NQ;i++)
		{
			for(j=1;j<Ng-1;j++)
			{
				if(i==0&&j==1)
				{
					Gs[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i+1]*Fs[NZ+1+i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
							+(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						  + ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];	
				}
				else if(i==NQ-1&&j==1)
				{
					Gs[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i]*Fs[NZ+1+i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];

					

				}
				else if(i==0)
				{
					Gs[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Fs[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i+1]*Fs[NZ+1+i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
				+(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];					
				}
				else if(i==NQ-1)
				{
				    Gs[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Fs[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i]*Fs[NZ+1+i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];
				}
				else if(j==1)
				{
					Gs[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -(ss[i+1]*Fs[NZ+1+i+1][j]+ss[i]*Fs[NZ+1+i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];
				}
				else
				{
					Gs[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Fs[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Fs[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -(ss[i+1]*Fs[NZ+1+i+1][j]+ss[i]*Fs[NZ+1+i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					        +(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Fs[NZ+1+i][j];
				}
               
                
			//printf("%f,",Fs[i][j]);
				//if(Gs[i][j]!=0){printf("%f,",Gs[i][j]);}
			}}


	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{   ii=sPflgi[j][i];
	 	    jj=sPflgj[j][i];
			if(sPflg[j][i]==0){b[0] = (Fs[ii][jj]-Fs[ii-1][jj])*byhz;
					b[1] = 0.25*(Fs[ii][1]-Fs[ii][0])*byhx;
					b[2] = Fs[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==1){b[0] = (Fs[ii][jj]-Fs[ii-1][jj])*byhz;
					b[1] = 0.25*(Fs[ii-1][1]-Fs[ii-1][0])*byhx;
					b[2] = Fs[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==2){b[0] = 0.125*(9.0*(Fs[ii][0]-Fs[ii-1][0])-(Fs[ii][1]-Fs[ii-1][1]))*byhz;
					b[1] = 0.25*(Fs[ii][1]-Fs[ii][0])*byhx;
					b[2] = Fs[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==3){b[0] = 0.125*(9.0*(Fs[ii][0]-Fs[ii-1][0])-(Fs[ii][1]-Fs[ii-1][1]))*byhz;
					b[1] = 0.25*(Fs[ii-1][1]-Fs[ii-1][0])*byhx;
					b[2] = Fs[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==4){b[0] = (Fs[ii][jj]-Fs[ii-1][jj])*byhz;
					b[1] = (Fs[ii][jj]-Fs[ii][jj-1])*byhx;
					b[2] = Fs[ii][jj] - b[0]*zm[ii] - b[1]*xm[jj];}
			else if(sPflg[j][i]==5){b[0] = (Fs[ii][jj]-Fs[ii-1][jj])*byhz;
					b[1] = (Fs[ii-1][jj]-Fs[ii-1][jj-1])*byhx;
					b[2] = Fs[ii-1][jj] - b[0]*zm[ii-1] - b[1]*xm[jj];}
			else if(sPflg[j][i]==6){b[0] = (Fs[ii][jj-1]-Fs[ii-1][jj-1])*byhz;
					b[1] = (Fs[ii][jj]-Fs[ii][jj-1])*byhx;
					b[2] = Fs[ii][jj-1] - b[0]*zm[ii] - b[1]*xm[jj-1];}
			else{b[0] = (Fs[ii][jj-1]-Fs[ii-1][jj-1])*byhz;
					b[1] = (Fs[ii-1][jj]-Fs[ii-1][jj-1])*byhx;
					b[2] = Fs[ii-1][jj-1] - b[0]*zm[ii-1] - b[1]*xm[jj-1];}

			
			Gs[j+NZ+1][i]= -b[0]*sPflgz[j][i]-b[1]*sPflgx[j][i]-b[2]+Fs[NZ+1+j][i];}}
	//for(i=0;i<NQ;i++){ 
		//Gs[NZ+1][0]=Fs[NZ+1][1];
		

			for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Ls+=IIs[i][j]*Gs[i][j];
				}
			}
			for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{
					Ls+=IIs[i+NZ+1][j]*Gs[i+NZ+1][j];
				}
			}
          // Ls+=IIs[NZ+1][0]*Gs[NZ+1][0];
	        Hs=Js/Ls;

			
			for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Ms[i][j]=Is[i][j]-Hs*Gs[i][j];
				}
			}
			for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{
					Ms[i+NZ+1][j]=Is[i+NZ+1][j]-Hs*Gs[i+NZ+1][j];
				}
			}
	
		  	//Ms[NZ+1][0]=Is[NZ+1][0]-Hs*Gs[NZ+1][0];


			for(i=0;i<NZ;i++)
		{
			for(j=0;j<NX;j++)
			{
				if(cflg[i][j]==0)
				{
					if(i==0&&j==0)
					{
						Ws[i][j]=-x[j+1]*Ms[i][j+1]*byhx*byhx*dt/xm[j]
							   -Ms[i+1][j]*byhz*byhz*dt
							   +dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else if(i==0&&j==NX-1)
					{
						Ws[i][j]=-x[j]*Ms[i][j-1]*byhx*byhx*dt/xm[j]
							   -Ms[i+1][j]*byhz*byhz*dt
						+dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else if(i==NZ-1&&j==0)
					{
						Ws[i][j]=-x[j+1]*Ms[i][j+1]*byhx*byhx*dt/xm[j]
							   -Ms[i-1][j]*byhz*byhz*dt
						+dt*(x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else if(i==NZ-1&&j==NX-1)
					{
						Ws[i][j]= -x[j]*Ms[i][j-1]*byhx*byhx*dt/xm[j]
							   -Ms[i-1][j]*byhz*byhz*dt
						+dt*(x[j]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else if(i==0)
					{
						Ws[i][j]= -(x[j+1]*Ms[i][j+1]+x[j]*Ms[i][j-1])*byhx*byhx*dt/xm[j]
							   -Ms[i+1][j]*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+0.5*byhz*byhz)*Ms[i][j];
					}
					else if(j==0)
					{
						Ws[i][j]=-x[j+1]*Ms[i][j+1]*byhx*byhx*dt/xm[j]
							   -(Ms[i+1][j]+Ms[i-1][j])*byhz*byhz*dt
						+2.0*dt*(0.5*x[j+1]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else if(i==NZ-1)
					{
						Ws[i][j]=-(x[j+1]*Ms[i][j+1]+x[j]*Ms[i][j-1])*byhx*byhx*dt/xm[j]
							   -Ms[i-1][j]*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+0.5*byhz*byhz)*Ms[i][j];
					}
					else if(j==NX-1)
					{
						Ws[i][j]= -x[j]*Ms[i][j-1]*byhx*byhx*dt/xm[j]
							   -(Ms[i+1][j]+Ms[i-1][j])*byhz*byhz*dt
						       +2.0*dt*(0.5*x[j]/xm[j]*byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					else
					{
						Ws[i][j]= -(x[j+1]*Ms[i][j+1]+x[j]*Ms[i][j-1])*byhx*byhx*dt/xm[j]
							   -(Ms[i+1][j]+Ms[i-1][j])*byhz*byhz*dt
								+2.0*dt*(byhx*byhx+byhz*byhz)*Ms[i][j];
					}
					
				}
				else if(cflg[i][j]==2){

					ii=cPflgi[i][j];
					jj=cPflgj[i][j];

					if(cPflg[i][j]==0){b[0] = 0.25*(Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii-1][jj+1])*byhs;
						b[1] = 0.125*(9.0*(Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])-(Ms[NZ+1+ii-1][jj+1]-Ms[NZ+1+ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==1){b[0] = 0.25*(Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==2){b[0] = 0.25*(Ms[NZ+1+ii][jj]-Ms[NZ+1+ii-1][jj])*byhs;
						b[1] = 0.125*(9.0*(Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])-(Ms[NZ+1+ii-1][jj+1]-Ms[NZ+1+ii-1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==3){b[0] = 0.25*(Ms[NZ+1+ii][jj]-Ms[NZ+1+ii-1][jj])*byhs;
						b[1] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==4){b[0] = 0.25*(Ms[NZ+1+1][jj+1]-Ms[NZ+1+0][jj+1])*byhs;
						b[1] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==5){b[0] = 0.25*(Ms[NZ+1+1][jj+1]-Ms[NZ+1+0][jj+1])*byhs;
						b[1] =  0.125*(9.0*(Ms[NZ+1+0][jj+1]-Ms[NZ+1+0][jj])-(Ms[NZ+1+1][jj+1]-Ms[NZ+1+1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj+1] - b[0]*0.5*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==6){b[0] = 0.25*(Ms[NZ+1+1][jj]-Ms[NZ+1+0][jj])*byhs;
						b[1] =(Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}
					else if(cPflg[i][j]==7){b[0] = 0.25*(Ms[NZ+1+1][jj]-Ms[NZ+1+0][jj])*byhs;
						b[1] = 0.125*(9.0*(Ms[NZ+1+0][jj+1]-Ms[NZ+1+0][jj])-(Ms[NZ+1+1][jj+1]-Ms[NZ+1+1][jj]))/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj] - b[0]*0.5*hs - b[1]*rgm[jj-1];}

					else if(cPflg[i][j]==8){b[0] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj+1] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==9){b[0] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii-1][jj+1])*byhs;
						b[1] = (Ms[NZ+1+ii-1][jj+1]-Ms[NZ+1+ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii-1][jj+1] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj];}
					else if(cPflg[i][j]==10){b[0] = (Ms[NZ+1+ii][jj]-Ms[NZ+1+ii-1][jj])*byhs;
						b[1] = (Ms[NZ+1+ii][jj+1]-Ms[NZ+1+ii][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii][jj] - b[0]*((double)ii+0.5)*hs - b[1]*rgm[jj-1];}
					else{b[0] = (Ms[NZ+1+ii][jj]-Ms[NZ+1+ii-1][jj])*byhs;
						b[1] = (Ms[NZ+1+ii-1][jj+1]-Ms[NZ+1+ii-1][jj])/(rgm[jj]-rgm[jj-1]);
						b[2] = Ms[NZ+1+ii-1][jj] - b[0]*((double)ii-0.5)*hs - b[1]*rgm[jj-1];}
					Ws[i][j]= -b[0]*cPflgs[i][j]-b[1]*cPflgr[i][j]-b[2]+Ms[i][j];}
				else{Ws[i][j]=0;}				
				
			}
		}
	
			//Ms[i][j]=Dnp[i][j]-Dnm[i][j];
			//Is[i][j]=Dnp[i][j]-Dnm[i][j];
			//if(Ms[i][j]!=0){printf("f%f\n",Ms[i][j]);}
			//if(Is[i][j]!=0){printf("i%f\n",Is[i][j]);}


		
	for(i=0;i<NQ;i++)
		{
			for(j=1;j<Ng-1;j++)
			{
				if(i==0&&j==1)
				{
					Ws[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i+1]*Ms[NZ+1+i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
							+(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						  + ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];	
				}
				else if(i==NQ-1&&j==1)
				{
					Ws[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i]*Ms[NZ+1+i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+( AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];
				}
				else if(i==0)
				{
					Ws[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Ms[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i+1]*Ms[NZ+1+i+1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
				+(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i+1]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];					
				}
				else if(i==NQ-1)
				{
				    Ws[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Ms[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -ss[i]*Ms[NZ+1+i-1][j]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +ss[i]/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];
				}
				else if(j==1)
				{
					Ws[i+NZ+1][j]=-AAold*AA*rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -(ss[i+1]*Ms[NZ+1+i+1][j]+ss[i]*Ms[NZ+1+i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					+(AAold*AA*rg[j]*rg[j]/rg2[j]/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];
				}
				else
				{
					Ws[i+NZ+1][j]=-AAold*AA*(rg[j]*rg[j]/rg2[j]*Ms[NZ+1+i][j+1]+rg[j-1]*rg[j-1]/rg2[j-1]*Ms[NZ+1+i][j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   -(ss[i+1]*Ms[NZ+1+i+1][j]+ss[i]*Ms[NZ+1+i-1][j])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs
					        +(AAold*AA*(rg[j]*rg[j]/rg2[j]+rg[j-1]*rg[j-1]/rg2[j-1])/(rgm[j-1]*rgm[j-1]*rgm[j-1])*dt*byhrg*byhrg
						   +(ss[i+1]+ss[i])/(rgm[j-1]*rgm2[j-1]*ssm[i])*dt*byhs*byhs)*Ms[NZ+1+i][j];
				}
                    
               
				
			//printf("%f,",Ms[i][j]);
				//if(Ws[i][j]!=0){printf("%f,",Ws[i][j]);}
			}}
	//for(i=0;i<NQ;i++){ 
	//Ws[NZ+1][0]=Ms[NZ+1][1];
//}
			
	for(i=Ng-1;i<=Ng;i++)
	{
		for(j=0;j<NQ;j++)
		{   ii=sPflgi[j][i];
	 	    jj=sPflgj[j][i];			
			
			if(sPflg[j][i]==0){b[0] = (Ms[ii][jj]-Ms[ii-1][jj])*byhz;
					b[1] = 0.25*(Ms[ii][1]-Ms[ii][0])*byhx;
					b[2] = Ms[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==1){b[0] = (Ms[ii][jj]-Ms[ii-1][jj])*byhz;
					b[1] = 0.25*(Ms[ii-1][1]-Ms[ii-1][0])*byhx;
					b[2] = Ms[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==2){b[0] = 0.125*(9.0*(Ms[ii][0]-Ms[ii-1][0])-(Ms[ii][1]-Ms[ii-1][1]))*byhz;
					b[1] = 0.25*(Ms[ii][1]-Ms[ii][0])*byhx;
					b[2] = Ms[ii][0] - b[0]*zm[ii] - b[1]*xm[0];}
            else if(sPflg[j][i]==3){b[0] = 0.125*(9.0*(Ms[ii][0]-Ms[ii-1][0])-(Ms[ii][1]-Ms[ii-1][1]))*byhz;
					b[1] = 0.25*(Ms[ii-1][1]-Ms[ii-1][0])*byhx;
					b[2] = Ms[ii-1][0] - b[0]*zm[ii-1] - b[1]*xm[0];}
			else if(sPflg[j][i]==4){b[0] = (Ms[ii][jj]-Ms[ii-1][jj])*byhz;
					b[1] = (Ms[ii][jj]-Ms[ii][jj-1])*byhx;
					b[2] = Ms[ii][jj] - b[0]*zm[ii] - b[1]*xm[jj];}
			else if(sPflg[j][i]==5){b[0] = (Ms[ii][jj]-Ms[ii-1][jj])*byhz;
					b[1] = (Ms[ii-1][jj]-Ms[ii-1][jj-1])*byhx;
					b[2] = Ms[ii-1][jj] - b[0]*zm[ii-1] - b[1]*xm[jj];}
			else if(sPflg[j][i]==6){b[0] = (Ms[ii][jj-1]-Ms[ii-1][jj-1])*byhz;
					b[1] = (Ms[ii][jj]-Ms[ii][jj-1])*byhx;
					b[2] = Ms[ii][jj-1] - b[0]*zm[ii] - b[1]*xm[jj-1];}
			else{b[0] = (Ms[ii][jj-1]-Ms[ii-1][jj-1])*byhz;
					b[1] = (Ms[ii-1][jj]-Ms[ii-1][jj-1])*byhx;
					b[2] = Ms[ii-1][jj-1] - b[0]*zm[ii-1] - b[1]*xm[jj-1];}

			Ws[j+NZ+1][i]= -b[0]*sPflgz[j][i]-b[1]*sPflgx[j][i]-b[2]+Ms[NZ+1+j][i];}}




		   for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Us+=Ws[i][j]*Ms[i][j];
					Vs+=Ws[i][j]*Ws[i][j];
				}
			}
			for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{
					Us+=Ws[i+NZ+1][j]*Ms[i+NZ+1][j];
					Vs+=Ws[i+NZ+1][j]*Ws[i+NZ+1][j];
				}
			}
			
					//Us+=Ws[NZ+1][0]*Ms[NZ+1][0];
					//Vs+=Ws[NZ+1][0]*Ws[NZ+1][0];
            Ns=Us/Vs;


			//printf("%f,%f,%f,",Hs,Js,Ls);
			 for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Es[i][j]+=Hs*Fs[i][j]+Ns*Ms[i][j];
					//dPc[i][j]= Es[i][j];
					//if(fabs(Hs*Fs[i][j]+Ns*Ms[i][j])>er){er=fabs(Hs*Fs[i][j]+Ns*Ms[i][j]);}
				}
			}
			for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{   Es[i+NZ+1][j]+=Hs*Fs[i+NZ+1][j]+Ns*Ms[i+NZ+1][j];
				   // dPs[i][j]= Es[i+NZ+1][j];
					//if(fabs(Hs*Fs[i+NZ+1][j]+Ns*Ms[i+NZ+1][j])>er){ii=i;jj=j;er=fabs(Hs*Fs[i+NZ+1][j]+Ns*Ms[i+NZ+1][j]);}
					
				}
			}
			//Es[NZ+1][0]+=Hs*Fs[NZ+1][0]+Ns*Ms[NZ+1][0];
		//output6(mm);

  //if(mm<10){
//		   output2(mm);
	 // }


		for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					Is[i][j]=Ms[i][j]-Ns*Ws[i][j];
					Ks+=Is[i][j]*IIs[i][j];
				}
			}
		for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{   Is[i+NZ+1][j]=Ms[i+NZ+1][j]-Ns*Ws[i+NZ+1][j];
					Ks+=Is[i+NZ+1][j]*IIs[i+NZ+1][j];
				}
			}
		//Is[NZ+1][0]=Ms[NZ+1][0]-Ns*Ws[NZ+1][0];
					//Ks+=Is[NZ+1][0]*IIs[NZ+1][0];


				for(i=0;i<NZ;i++)
					{
						for(j=0;j<NX;j++)
						{
							Fs[i][j]=Is[i][j]+Hs/Ns*Ks/Js*(Fs[i][j]-Ns*Gs[i][j]);
						}
					}
				for(i=0;i<NQ;i++)
					{
						for(j=1;j<=Ng;j++)
						{
							Fs[i+NZ+1][j]=Is[i+NZ+1][j]+Hs/Ns*Ks/Js*(Fs[i+NZ+1][j]-Ns*Gs[i+NZ+1][j]);
							//Fs[i+NZ+1][j]=Fs[NZ+1][j];
						}
					}
				//Fs[NZ+1][0]=Is[NZ+1][0]+Hs/Ns*Ks/Js*(Fs[NZ+1][0]-Ns*Gs[NZ+1][0]);
                er1=0;
				for(i=0;i<NZ;i++)
					{
						for(j=0;j<NX;j++)
						{   er1+=Is[i][j]*Is[i][j];
							
						}
					}
				for(i=0;i<NQ;i++)
					{
						for(j=0;j<=Ng;j++)
						{
							er1+=Is[i+NZ+1][j]*Is[i+NZ+1][j];
						}
					}
	//	 printf("er=%f,er1=%f,er2=%f\n",er,er1,er2);
				er3=sqrt(er1/er2);
     // printf("i=%d,j=%d,er=%f\n",ii,jj,er3);
	
	   if(er>er3){
		   er=er3;
	    for(i=0;i<NZ;i++)
			{
				for(j=0;j<NX;j++)
				{
					
					dPc[i][j]= Es[i][j];
					//if(fabs(Hs*Fs[i][j]+Ns*Ms[i][j])>er){er=fabs(Hs*Fs[i][j]+Ns*Ms[i][j]);}
				}
			}
			for(i=0;i<NQ;i++)
			{
				for(j=1;j<=Ng;j++)
				{  
				    dPs[i][j]= Es[i+NZ+1][j];
					//if(fabs(Hs*Fs[i+NZ+1][j]+Ns*Ms[i+NZ+1][j])>er){ii=i;jj=j;er=fabs(Hs*Fs[i+NZ+1][j]+Ns*Ms[i+NZ+1][j]);}
					
				}
			}}

	  			
		}while(er>erefCGT&&mm<MCGmax);
						

	//if(mm==MCGmax)
	//{
		//printf("à≥óÕÇ™é˚ë©ÇµÇ‹ÇπÇÒÇ≈ÇµÇΩ\ter=%le\n",er);
		printf("er=%le\n",er);
	//	exit(1);
	//}

	/*spherical*/
	for(i=0;i<NQ;i++)
	{
		for(j=0;j<Ng-1;j++)
		{
			rogs[i][j] = qrogs[i][j];
		}
	}
	AA = 1.0/log(rmax/Rsold);
	for(i=0;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			qur[i][j]+= -AAold*dt*(dPs[i][j+1]-dPs[i][j])*byhrg/rg2[j];
			ur[i][j] = 2.0*qur[i][j]/(rogs[i][j+1]+rogs[i][j]);
			//ur[i][j]=ur[0][j];
		}
	}
	for(i=1;i<NQ;i++)
	{
		for(j=1;j<Ng-1;j++)
		{
			quq[i][j]+= -(dPs[i][j]-dPs[i-1][j])*dt*byhs/rgm2[j-1];
			uq[i][j] = 2.0*quq[i][j]/(rogs[i][j]+rogs[i-1][j]);
			//uq[i][j] =0;
		}
	}
	
	
	for(i=0;i<NQ;i++)
	{
		for(j=1;j<=Ng;j++)
		{
			Ps[i][j]+= dPs[i][j];
			propertygs_cal(i,j);
		}
	}

	/*cylindrical*/
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				rogc[i][j] = qrogc[i][j];
			}
		}
	}
	for(i=1;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(czflg[i][j]==0)
			{
				quz[i][j]+= -(dPc[i][j]-dPc[i-1][j])*dt*byhz;
				uz[i][j] = 2.0*quz[i][j]/(rogc[i][j]+rogc[i-1][j]);
			}
		}
	}
	for(i=0;i<NZ;i++)
	{
		for(j=1;j<NX;j++)
		{
			if(cxflg[i][j]==0)
			{
				qux[i][j]+= -(dPc[i][j]-dPc[i][j-1])*dt*byhx;
				ux[i][j] = 2.0*qux[i][j]/(rogc[i][j]+rogc[i][j-1]);
			}
		}
	}
	for(i=0;i<NZ;i++)
	{
		for(j=0;j<NX;j++)
		{
			if(cflg[i][j]==0)
			{
				Pc[i][j]+=dPc[i][j];
				propertygc_cal(i,j);
			}
		}
	}
    freedvec(b,3);
	return mm;

}





void output2(int k)
{
	char chark[20];
	char filename2[50];
	FILE *fp;
	int i,j,g,mm,nn;

	sprintf(chark,"%d",k);
	strcpy(filename2,filenamea);
	strcat(filename2,chark);
	strcat(filename2,".dat");

	OUP += 1;

	if((fp=fopen(filename2,"w"))==NULL)
	{
		printf("Can't open datafile\n");
		exit(1);
	}

	for(j=0;j<NX;j++)
		{  nn=0;
	   for(i=0;i<NZ;i++)
	       {  
			if(cflg[i][j]==1){
				if(nn!=1){
					nn=1;
					mm=0;
									  do{mm++;}while(xm[j]<rg[0]*ssm[NQ/2+mm]&&mm<NQ/2-1);
									  fprintf(fp,"%le\t",STN*M[nO2]*Ys[NQ/2+mm][0][nFUEL]/(M[nFUEL]*Ys[NQ/2+mm][0][nO2]));
									  fprintf(fp,"%le\t",Tgs[NQ/2+mm][0]);
									  
									  fprintf(fp,"%le\t",ur[NQ/2+mm][0]);
									  
									  fprintf(fp,"\n");
									  fprintf(fp,"%le\t",STN*M[nO2]*Ys[NQ/2-mm][0][nFUEL]/(M[nFUEL]*Ys[NQ/2-mm][0][nO2]));
									  fprintf(fp,"%le\t",Tgs[NQ/2-mm][0]);
									
									  fprintf(fp,"%le\t",ur[NQ/2-mm][0]);
									 
									  fprintf(fp,"\n");
				}}
		else{fprintf(fp,"%le\t",STN*M[nO2]*Yc[i][j][nFUEL]/(M[nFUEL]*Yc[i][j][nO2]));
			fprintf(fp,"%le\t",Tgc[i][j]);
			
			fprintf(fp,"%le\t",sqrt(pow(ux[i][j]+ux[i][j+1],2)+pow(uz[i][j]+uz[i+1][j],2))/2);
		
			fprintf(fp,"\n");
		 }
	   }fprintf(fp,"\n");
	}
	fclose(fp);

}
void output3(int k)
{
	char chark[20];
	char filename2[50];
	FILE *fp;
	int i,j,g,mm,nn;

	sprintf(chark,"%d",k);
	strcpy(filename2,filenameb);
	strcat(filename2,chark);
	strcat(filename2,".dat");

	OUP += 1;

	if((fp=fopen(filename2,"w"))==NULL)
	{
		printf("Can't open datafile\n");
		exit(1);
	}
	fprintf(fp,"TITLE = ""HEAT MASS TRANSFER""\n");
	fprintf(fp,"VARIABLES = ""z"", ""r"", ""U/RFC"", ""V/RFC"", ""\n");
	fprintf(fp,"ZONE I=%d,J=%d\n",NZ,NX);
	for(j=0;j<NX;j++)
		{  
		for(i=0;i<NZ;i++)
	       {  
			if(cflg[i][j]==0){
				fprintf(fp,"%le\t",zm[i]);
			    fprintf(fp,"%le\t",xm[j]);
				fprintf(fp,"%le\t",(uz[i][j]+uz[i+1][j])/2);
			    fprintf(fp,"%le\t",(ux[i][j]+ux[i][j+1])/2);
			    fprintf(fp,"\n");
		 }
			else{fprintf(fp,"%le\t",0);
			    fprintf(fp,"%le\t",0);
				fprintf(fp,"%le\t",0);
			    fprintf(fp,"%le\t",0);
				fprintf(fp,"\n");}
	   }
	}
	fclose(fp);

}

void output4(int k)
{
	char chark[20];
	char filename2[50];
	FILE *fp;
	int i,j,g,mm,nn,sp;

	sprintf(chark,"%d",k);
	strcpy(filename2,filenamec);
	strcat(filename2,chark);
	strcat(filename2,".dat");

	OUP += 1;

	if((fp=fopen(filename2,"w"))==NULL)
	{
		printf("Can't open datafile\n");
		exit(1);
	}
fprintf(fp,"%le\t",tref);
fprintf(fp,"%le\t",dt);
fprintf(fp,"%le\t",Rs);

	for(j=0;j<NX;j++)
		{  
		for(i=0;i<NZ;i++)
	       {  
				fprintf(fp,"%le\t",Tgc[i][j]);
				
				for(sp=0;sp<spmax;sp++){
					fprintf(fp,"%le\t",Yc[i][j][sp]);
				    }
				fprintf(fp,"%le\t",dPc[i][j]);
				fprintf(fp,"%le\t",Pc[i][j]);
	   }
	}
	for(j=0;j<NX;j++)
		{  
		for(i=0;i<=NZ;i++)
		{  fprintf(fp,"%le\t",uz[i][j]);}}

    for(j=0;j<=NX;j++)
		{  
		for(i=0;i<NZ;i++)
		{ fprintf(fp,"%le\t",ux[i][j]);}}

	for(j=0;j<=Nl;j++)
		{  
		for(i=0;i<NQ;i++)
		{ fprintf(fp,"%le\t",Tls[i][j]);}}


	for(j=0;j<=Ng;j++)
		{  
		for(i=0;i<NQ;i++)
		{ fprintf(fp,"%le\t",Tgs[i][j]);
	
				for(sp=0;sp<spmax;sp++){
					fprintf(fp,"%le\t",Ys[i][j][sp]);
				    }
				fprintf(fp,"%le\t",dPs[i][j]);
				fprintf(fp,"%le\t",Ps[i][j]);
				fprintf(fp,"%le\t",ur[i][j]);
		}}
	for(j=0;j<Ng;j++)
		{  
		for(i=0;i<=NQ;i++)
		{fprintf(fp,"%le\t",uq[i][j]);}}

	fclose(fp);

}

