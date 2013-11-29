#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>

#include "DynamicalSystem.hxx"

using namespace std;

//________________________________________________________________________
DynamicalSystem::DynamicalSystem()
{
	
	SetPrecision();
	fHestimate=1.; //this value is change in DynamicalSystem::RungeKutta
	fHmin=0.;
	fMaxHdid=-1e30;
	fMinHdid=1e30;
	fIsNegativeValueAllowed=true;
}

//________________________________________________________________________
DynamicalSystem::DynamicalSystem(const DynamicalSystem & DS)
{
	fNVar=DS.fNVar;
	fPrecision=DS.fPrecision;
	fHestimate=DS.fHestimate;
	fHmin=DS.fHmin;
	fMaxHdid=DS.fMaxHdid;
	fMinHdid=DS.fMinHdid;
	fIsNegativeValueAllowed=DS.fIsNegativeValueAllowed;
}

//________________________________________________________________________
DynamicalSystem::~DynamicalSystem()
{
	
}

//________________________________________________________________________
void  DynamicalSystem::RungeKutta(double *YStart,double t1, double t2, int EquationNumber)
{
	//double shortestHalfLife=gMURE->GetShortestHalfLife();
	fNVar = EquationNumber;
	int nstp;
	double t,hnext,hdid,h;
	double *yscal=new double[fNVar];
	double *y=new double[fNVar];
	double *dydt=new double[fNVar];
	
	const double MAXSTP = 10000;
	const double TINY = 1.0e-30;
 	
	
	
	if(fabs(t1-t2)<=TINY)
		cout << "Integration time is 0." << endl;
	
	t=t1;
	fHestimate=(t2-t1)/100;
	//h=(t2 > t1) ? fabs(fHestimate) : -fabs(fHestimate);
	h=fHestimate;
	
	// pragma omp parallel for
	for (int i = 0; i < fNVar; i++) y[i] = YStart[i];
	
	for ( nstp = 1; nstp <= MAXSTP; nstp++)
	{
		BuildEqns(t,y,dydt);
		
		// pragma omp parallel for
		for ( int i = 0; i < fNVar; i++)
			yscal[i] = fabs(y[i]) + fabs(dydt[i]*h)+TINY;
		
		if ( (t+h-t2) * (t+h-t1) > 0.0 ) h=t2-t;
		
		AdaptStepSize(y,dydt,&t,h,fPrecision,yscal,&hdid,&hnext);
		
		if(fMaxHdid<hdid) fMaxHdid=hdid;
		if(fMinHdid>hdid) fMinHdid=hdid;
		
		if ((t-t2)*(t2-t1) >= 0.0)
		{
			// pragma omp parallel for
			for (int i=0;i<fNVar;i++) YStart[i]=y[i];
			
			delete [] dydt;
			delete [] y;
			delete [] yscal;
			//cout << "The maximum step used in RK was "<<fMaxHdid<<" Step NUM in RK was "<<nstp << endl;
			return;
		}
		if (fabs(hnext) <= fHmin)
			cout << "Step size too small in RungeKutta" << endl;
		
		h=hnext;
	}
	cout <<  "Too many steps in routine RungeKutta" << endl;
}

//________________________________________________________________________
void  DynamicalSystem::RK4(double *y, double *dydx, double x, double h, double *yout)
{
	//cout<<"Calling Function RK4"<<endl;
	double xh,hh,h6;
	
	double *dym=new double[fNVar];
	double *dyt=new double[fNVar];
	double *yt=new double[fNVar];
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	
	// pragma omp parallel for
	for (int i=0;i<fNVar;i++) yt[i]=y[i]+hh*dydx[i];
	
	BuildEqns(xh,yt,dyt);
	
	// pragma omp parallel for
	for (int i=0;i<fNVar;i++) yt[i]=y[i]+hh*dyt[i];
	
	BuildEqns(xh,yt,dym);
	
	for (int i=0;i<fNVar;i++)
	{
		yt[i]=y[i]+h*dym[i];
		dym[i] += dyt[i];
	}
	
	BuildEqns(x+h,yt,dyt);
	// pragma omp parallel for
	for (int i=0;i<fNVar;i++)
	{
		yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
		if(!fIsNegativeValueAllowed && yout[i]<0)
		{
			cout << "Material composition is negative "
			<<"i="<<i<<" ("/*<<fEvolvingMaterial->GetComposition()[i]->GetZAI()->PrintName()
					*/<<") old proportion="<<y[i]<<" new="<<yout[i]
			<<". Setting to 0." << endl;
			yout[i]=0.;
		}
	}
	delete [] yt;
	delete [] dyt;
	delete [] dym;
}

//________________________________________________________________________
void DynamicalSystem::AdaptStepSize(double *y, double *dydx, double *x, double htry,
				    double eps, double *yscal, double *hdid, double *hnext)
{
	//cout<<"Calling Function AdaptStepSize()"<<endl;
	double xsav,hh,h,temp,errmax;
	double *dysav=new double[fNVar];
	double *ysav=new double[fNVar];
	double *ytemp=new double[fNVar];
	
	const double PGROW =-0.20;
	const double PSHRNK =-0.25;
	const double FCOR =0.06666666 ;
	const double SAFETY =0.9;
	const double ERRCON =6.0e-4;
	
	xsav=(*x);
	
	// pragma omp parallel for
 	for (int i=0;i<fNVar;i++)
	{
		ysav[i]=y[i];
		dysav[i]=dydx[i];
	}
	h=htry;
	for (;;)
	{
		hh=0.5*h;
		RK4(ysav,dysav,xsav,hh,ytemp);
		*x=xsav+hh;
		
		BuildEqns(*x,ytemp,dydx);
		RK4(ytemp,dydx,*x,hh,y);
		*x=xsav+h;
		if (*x == xsav )
		{
			//cout << "Step size ("<<h<<") too small in routine AdaptStepSize" << endl;
		}
		RK4(ysav,dysav,xsav,h,ytemp);
		errmax=0.0;
		
		for (int i=0;i<fNVar;i++)
		{
			ytemp[i]=y[i]-ytemp[i];
			temp=fabs(ytemp[i]/yscal[i]);
			if (errmax < temp) errmax=temp;
		}
		errmax /= eps;
		if (errmax <= 1.0)
		{
			*hdid=h;
			*hnext=(errmax > ERRCON ? SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
			break;
		}
		h=SAFETY*h*exp(PSHRNK*log(errmax));
	}
	
	for (int i=0;i<fNVar;i++) y[i] += ytemp[i]*FCOR;
        
	delete [] ytemp;
	delete [] dysav;
	delete [] ysav;
}

