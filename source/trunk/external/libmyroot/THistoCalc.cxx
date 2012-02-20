#include "THistoCalc.hxx"
#include <TF1.h>
#include <iostream>
#include <math.h>
#include <TSpline.h>
#include <fstream>
using namespace std;

////////////////////////////////////////////////////////////////////
//                           THistoCalc




THistoCalc::THistoCalc(TH1F& histo) : // Constructeur par recopie
  TH1D(NextName(),"",histo.GetNbinsX(),histo.GetXaxis()->GetXbins()->fArray) {

  for (Int_t i=1;i<=GetNbinsX();i++)
    SetBinValErr(i,ValErr_t(histo.GetBinContent(i),histo.GetBinError(i)));

}


ValErr_t THistoCalc::GetBinValErr(Int_t bin) {

  return ValErr_t(GetBinContent(bin),GetBinError(bin));

}

void THistoCalc::SetBinValErr(Int_t bin, const ValErr_t& ve) {

  SetBinContent(bin,ve.fVal);
  SetBinError(bin,ve.fErr);

}

void THistoCalc::Write(ostream& s) {

  s << "@type xydy\n";
  for(Int_t i=1 ; i<=GetNbinsX(); i++) {
    s << GetBinCenter(i) << " " << GetBinContent(i) << " " <<
         GetBinError(i) <<endl;
  }
  s << "&\n";

}

void THistoCalc::AddBinSizeErr() {

  ValErr_t ve=0;
  Int_t n=GetNbinsX();
  if (n<2)
    return;

  ve.fErr=fabs(GetBinContent(2)-GetBinContent(1))/2;
  AddBinValErr(1,ve);

  for (Int_t i=2;i<n;i++) {
    ve.fErr=fabs(GetBinContent(i)-GetBinContent(i-1))/2;
    if (fabs(GetBinContent(i)-GetBinContent(i+1))/2 > ve.fErr)
      ve.fErr=fabs(GetBinContent(i)-GetBinContent(i+1))/2;
    AddBinValErr(i,ve);
  }

  ve.fErr=fabs(GetBinContent(n)-GetBinContent(n-1))/2;
  AddBinValErr(n,ve);

}


void THistoCalc::CenterToMean() {
  
  if (GetNbinsX()<3)
    return;

  TAxis *axe=GetXaxis();
  ValErr_t a;
  ValErr_t y1,y2;
  Double_t x1,x2;
  Double_t w1,w2,w3;

  w1=axe->GetBinWidth(1);
  w2=axe->GetBinWidth(2);
  w3=axe->GetBinWidth(3);
  x1=(w1+w2)/2;
  x2=(w2+w3)/2;
  y1=GetBinValErr(1)-GetBinValErr(2);
  y2=GetBinValErr(3)-GetBinValErr(2);

  a=(y1/x1+y2/x2)/(x1+x2)/12.0;
  SetBinValErr(1,GetBinValErr(1)+a*w1*w1);

  w3=w2;
  x2=x1;
  y2=-y1;

  for (Int_t i=2;i<GetNbinsX();i++) { // au moins une fois

    w2=w3;
    w3=axe->GetBinWidth(i+1);
    x1=x2;
    x2=(w2+w3)/2;
    y1=-y2;
    y2=GetBinValErr(i+1)-GetBinValErr(i);

    a=(y1/x1+y2/x2)/(x1+x2)/12.0;
    SetBinValErr(i,GetBinValErr(i)+a*w2*w2);
  }
  
  SetBinValErr(GetNbinsX(),GetBinValErr(GetNbinsX())+a*w3*w3);

}


ValErr_t THistoCalc::Integrale(Double_t xmin,Double_t xmax) {

  ValErr_t integrale;
  Double_t invert=1;

  TAxis *axe=GetXaxis();

  if (xmin>xmax) {
    invert=-1;
    Double_t temp=xmin;
    xmin=xmax;
    xmax=temp;
  }
  
  if (xmin<axe->GetBinLowEdge(1))
    xmin=axe->GetBinLowEdge(1);
  if (xmax>axe->GetBinUpEdge(GetNbinsX()))
    xmax=axe->GetBinUpEdge(GetNbinsX());

  Int_t bin=1;
  while(axe->GetBinUpEdge(bin)<xmin)
    bin++;

  integrale=GetBinValErr(bin)*(axe->GetBinUpEdge(bin)-xmin);
  bin++;

  while(axe->GetBinUpEdge(bin)<xmax && bin <= GetNbinsX()) {
    integrale += GetBinValErr(bin)*(axe->GetBinUpEdge(bin)-axe->GetBinLowEdge(bin));
    bin++;
  }

  if (bin <= GetNbinsX())
    integrale += GetBinValErr(bin)*(xmax-axe->GetBinLowEdge(bin));

  return integrale*invert;

}
ValErr_t THistoCalc::LogIntegrale(Double_t xmin,Double_t xmax) {

  ValErr_t integrale;
  Double_t invert=1;

  TAxis *axe=GetXaxis();

  if (xmin>xmax) {
    invert=-1;
    Double_t temp=xmin;
    xmin=xmax;
    xmax=temp;
  }
  
  if (xmin<axe->GetBinLowEdge(1))
    xmin=axe->GetBinLowEdge(1);
  if (xmax>axe->GetBinUpEdge(GetNbinsX()))
    xmax=axe->GetBinUpEdge(GetNbinsX());

  Int_t bin=1;
  while(axe->GetBinUpEdge(bin)<xmin)
    bin++;
  
  integrale=GetBinValErr(bin)*(log(axe->GetBinUpEdge(bin))-log(xmin));
  bin++;

  while(axe->GetBinUpEdge(bin)<xmax && bin <= GetNbinsX()) {
    integrale += GetBinValErr(bin)*(log(axe->GetBinUpEdge(bin))-log(axe->GetBinLowEdge(bin)));
	bin++;
  }

  if (bin <= GetNbinsX())
    integrale += GetBinValErr(bin)*(log(xmax)-log(axe->GetBinLowEdge(bin)));
  return integrale*invert;

}


ValErr_t THistoCalc::First_Moment(Double_t xmin,Double_t xmax)
{ ValErr_t integrale;
  Double_t invert=1;
  TAxis *axe=GetXaxis();
  if (xmin>xmax)
     invert = -1 ;

  if (xmin<axe->GetBinLowEdge(1))
    xmin=axe->GetBinLowEdge(1);
  if (xmax>axe->GetBinUpEdge(GetNbinsX()))
    xmax=axe->GetBinUpEdge(GetNbinsX());

  Int_t bin=1;     
  while(axe->GetBinUpEdge(bin)<xmin)
     bin++;
  integrale=GetBinValErr(bin)*(xmin+(axe->GetBinUpEdge(bin)-xmin)/2);
  bin++;
  while(axe->GetBinUpEdge(bin)<xmax && bin <= GetNbinsX()) 
    {integrale += GetBinValErr(bin)*(axe->GetBinLowEdge(bin)+(axe->GetBinUpEdge(bin)-axe->GetBinLowEdge(bin))/2);
     bin++;
    }
  if (bin <= GetNbinsX())
    integrale +=
    GetBinValErr(bin)*(axe->GetBinLowEdge(bin)+(xmax-axe->GetBinLowEdge(bin))/2);

  return integrale*invert;

}



TGraphErrors *THistoCalc::Omega(Double_t tmin, Double_t tmax, Int_t nbins) {

  TF1 Fexpo(NextName(),"expo");

  Int_t n=GetNbinsX();
  Double_t *x=new Double_t[n/nbins];
  Double_t *dx=new Double_t[n/nbins];
  Double_t *y=new Double_t[n/nbins];
  Double_t *dy=new Double_t[n/nbins];
  Double_t xmin,xmax;

  Int_t j=0;
  Int_t i=0;

  // init i.
  while (i+nbins-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbins)+
	 GetBinWidth(i+nbins) < 2.0*tmin)
    i++;

  // calculs omegas sur intervalles disjoints.
  while (i+nbins-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbins)+
	 GetBinWidth(i+nbins) < 2.0*tmax) {
    xmin=GetBinLowEdge(i+1);
    xmax=GetBinLowEdge(i+nbins)+GetBinWidth(i+nbins);
    x[j]=(xmin+xmax)/2.0;
    dx[j]=0;
    Fit(&Fexpo,"QN",0,xmin,xmax);
    y[j]=-Fexpo.GetParameter(1);
    dy[j]=Fexpo.GetParError(1);
    j++;
    i+=nbins;
  }
  
  TGraphErrors *retour=new TGraphErrors(j,x,y,dx,dy);
  delete x;
  delete dx;
  delete y;
  delete dy;
  return retour;

}

TGraphErrors *THistoCalc::OmegaGlissant(Double_t tmin, Double_t tmax, Int_t nbins) {

  TF1 Fexpo(NextName(),"expo");

  Int_t n=GetNbinsX();
  Double_t *x=new Double_t[n];
  Double_t *dx=new Double_t[n];
  Double_t *y=new Double_t[n];
  Double_t *dy=new Double_t[n];
  Double_t xmin,xmax;

  Int_t j=0;
  Int_t i=0;

  // init i.
  while (i+nbins-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbins)+
	 GetBinWidth(i+nbins) < 2.0*tmin)
    i++;

  // calculs omegas sur intervalles disjoints.
  while (i+nbins-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbins)+
	 GetBinWidth(i+nbins) < 2.0*tmax) {
    xmin=GetBinLowEdge(i+1);
    xmax=GetBinLowEdge(i+nbins)+GetBinWidth(i+nbins);
    x[j]=(xmin+xmax)/2.0;
    dx[j]=0;
    Fit(&Fexpo,"QN",0,xmin,xmax);
    y[j]=-Fexpo.GetParameter(1);
    dy[j]=Fexpo.GetParError(1);
    j++;
    i++;
  }
  
  TGraphErrors *retour=new TGraphErrors(j,x,y,dx,dy);
  delete x;
  delete dx;
  delete y;
  delete dy;
  return retour;

}

TGraphErrors *THistoCalc::AlphaPoly(Double_t tmin,Double_t tmax,Int_t nbinfit) 
{

 TF1 FPol(NextName(),"pol1"); 

  Int_t n=GetNbinsX();
  Double_t *x=new Double_t[n/nbinfit];
  Double_t *y=new Double_t[n/nbinfit];
  Double_t *dx=new Double_t[n/nbinfit];
  Double_t *dy=new Double_t[n/nbinfit];
  Double_t xmin,xmax;

  Int_t j=0;
  Int_t i=0;

  // init i.
  while (i+nbinfit-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbinfit)+
	 GetBinWidth(i+nbinfit) < 2.0*tmin)
    i++;

  // calculs omegas sur intervalles disjoints.
  while (i+nbinfit-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbinfit)+
	 GetBinWidth(i+nbinfit) < 2.0*tmax) {
    xmin=GetBinLowEdge(i+1);
    xmax=GetBinLowEdge(i+nbinfit)+GetBinWidth(i+nbinfit);
    x[j]=(xmin+xmax)/2.0;
    dx[j]= 0 ;
    Fit(&FPol,"QN",0,xmin,xmax);

    y[j]=-FPol.GetParameter(1);
    dy[j]=FPol.GetParError(1);
//    y[j]=-FPol.Derivative(x[j]);
//    dy[j]=sqrt( FPol.GetParError(1)*FPol.GetParError(1)*(xmin+xmax)/2.0 );
 
    j++;
    i+=nbinfit;
  }
  
  TGraphErrors *retour=new TGraphErrors(j,x,y,dx,dy);
  delete x;
  delete y;
  delete dx;
  delete dy;
  return retour;

}




TGraphErrors *THistoCalc::AlphaPolyGlissant(Double_t tmin,Double_t tmax,Int_t nbinfit)

{
 TF1 FPol(NextName(),"pol1");   
 int n= GetNbinsX();
 double *x=new Double_t[n];
 double *dx=new Double_t[n];
 double *y=new Double_t[n];
 double *dy=new Double_t[n];

 double xmin,xmax;
  int j=0;
  int i=0;

  while (i+nbinfit-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbinfit)+
	 GetBinWidth(i+nbinfit) < 2.0*tmin)
    i++;

  while (i+nbinfit-1<n &&
	 GetBinLowEdge(i+1)+GetBinLowEdge(i+nbinfit)+
	 GetBinWidth(i+nbinfit) < 2.0*tmax) {
    xmin=GetBinLowEdge(i+1);
    xmax=GetBinLowEdge(i+nbinfit)+GetBinWidth(i+nbinfit);
    x[j]=(xmin+xmax)/2.0;
    dx[j]= 0 ;
    Fit(&FPol,"QN",0,xmin,xmax);
    y[j]=-FPol.GetParameter(1);
    dy[j]=FPol.GetParError(1);
//    y[j]=-FPol.Derivative(x[j]);
//   dy[j]=sqrt( FPol.GetParError(1)*FPol.GetParError(1)*(xmin+xmax)/2.0 );		
 	
    j++;
    i++;
  }

  TGraphErrors *retour=new TGraphErrors(j,x,y,dx,dy);
  delete x;
  delete y;
  delete dx;
  delete dy;
  return retour ;
}

/* Fonction qui renvoie un histo contenant la diff bin a bin de this */
THistoCalc *THistoCalc::DiffBin()

{THistoCalc *Temp = new THistoCalc(this->GetNbinsX()-1,this->GetXaxis()->GetBinCenter(1),this->GetXaxis()->GetBinCenter(this->GetNbinsX() ));
 for (int i=0;i<Temp->GetNbinsX();i++)
	Temp->SetBinValErr(i+1,this->GetBinValErr(i+1)- this->GetBinValErr(i+2));
 return Temp;
}

// The next 3 functions are only used internally
// with errors bars 


void DaubTransc(int j, ValErr_t *vector, int p, double *coefs) {
  // Forward Daubechies wavelet transform of vector (size 2^j)
  // using the filter coefficients coefs (size 2*p)
  // In place.

  ValErr_t tempc[1<<(j-1)],tempd[1<<(j-1)];
  int sgn;
  // Faire une itération cj => moitié de scaling cj-1, moitié de wavelet dj-1
  for (int l=0;l<(1<<(j-1));l++) {
    tempc[l]=tempd[l]=0;
    sgn=1;
    for (int k=0;k<2*p;k++) {
      tempc[l]+=      coefs[k]      *vector[(2*l+k)%(1<<j)];
      tempd[l]+=sgn * coefs[2*p-1-k]*vector[(2*l+k)%(1<<j)];
      sgn=-sgn;
    }
  }
  for (int l=0;l< (1<<(j-1));l++) {
    vector[l]=tempc[l];
    vector[l+(1<<(j-1))]=tempd[l];
  }

  // S'appeler pour la moitié du vecteur

  if (j>1)
    DaubTransc(j-1,vector,p,coefs);
}

void InvDaubTransc(int j, ValErr_t *vector, int p, double *coefs) {

  // Inverse Daubechies wavelet transform of vector (size 2^j)
  // using the filter coefficients coefs (size 2*p)
  // In place.

  if (j>1)
    InvDaubTransc(j-1,vector,p,coefs);
 
  ValErr_t tempc[1<<j];
  int sgn=1;
  // Faire une itération cj-1 et dj-1 => cj.
  for (int l=0;l<(1<<j);l++) {
    tempc[l]=0;

//      colonne=((1<<(j-1))+1-p)%(1<<(j-1));
//      icoef=2*p+(l%2)-2;
//      while(icoef>=0) {
//        tempc[l]+=    coefs[icoef]      *vector[colonne];
//        tempc[l]+=sgn*coefs[2*p-1-icoef]*vector[colonne+(1<<(j-1))];
//        colonne=(colonne+1)%(1<<(j-1));
//        icoef-=2;
//      }


    int indice;
    for (int n=(l-2*p+1)/2;n<=l/2;n++) {
      indice=n;
      while (indice<0)
	indice+=1<<(j-1);
      tempc[l]+=vector[indice]                  * coefs[l-2*n];
      tempc[l]+=vector[indice+(1<<(j-1))] * sgn * coefs[2*p-1-l+2*n];
    }


    sgn=-sgn;
  }

  for (int l=0;l<(1<<j);l++)
    vector[l]=tempc[l];
}

double *DaubCoefs(int p) {

  double *coefs;

  switch(p) {
  case 2:
    coefs=new double[4];
    coefs[0] = 0.482962913145;
    coefs[1] = 0.836516303738;
    coefs[2] = 0.224143868042;
    coefs[3] =-0.129409522551;
    break;
  case 3:
    coefs=new double[6];
    coefs[0] = 0.332670552950;
    coefs[1] = 0.806891509311;
    coefs[2] = 0.459877502118;
    coefs[3] =-0.135011020010;
    coefs[4] =-0.085441273882;
    coefs[5] = 0.035226291882;
    break;
  case 4:
    coefs=new double[8];
    coefs[0] = 0.230377813309;
    coefs[1] = 0.714846570553;
    coefs[2] = 0.630880767930;
    coefs[3] =-0.027983769417;
    coefs[4] =-0.187034811719;
    coefs[5] = 0.030841381836;
    coefs[6] = 0.032883011667;
    coefs[7] =-0.010597401785;
    break;
  case 5:
    coefs=new double[10];
    coefs[0] = 0.160102397974;
    coefs[1] = 0.603829269797;
    coefs[2] = 0.724308528438;
    coefs[3] = 0.138428145901;
    coefs[4] =-0.242294887066;
    coefs[5] =-0.032244869585;
    coefs[6] = 0.077571493840;
    coefs[7] =-0.006241490213;
    coefs[8] =-0.012580751999;
    coefs[9] = 0.003335725285;
    break;
  case 6:
    coefs=new double[12];
    coefs[0] = 0.111540743350;
    coefs[1] = 0.494623890398;
    coefs[2] = 0.751133908021;
    coefs[3] = 0.315250351709;
    coefs[4] =-0.226264693965;
    coefs[5] =-0.129766867567;
    coefs[6] = 0.097501605587;
    coefs[7] = 0.027522865530;
    coefs[8] =-0.031582039317;
    coefs[9] = 0.000553842201;
    coefs[10]= 0.004777257511;
    coefs[11]=-0.001077301085;
    break;
  case 7:
    coefs=new double[14];
    coefs[0] = 0.077852054085;
    coefs[1] = 0.396539319482;
    coefs[2] = 0.729132090846;
    coefs[3] = 0.469782287405;
    coefs[4] =-0.143906003929;
    coefs[5] =-0.224036184994;
    coefs[6] = 0.071309219267;
    coefs[7] = 0.080612609151;
    coefs[8] =-0.038029936935;
    coefs[9] =-0.016574541631;
    coefs[10]= 0.012550998556;
    coefs[11]= 0.000429577973;
    coefs[12]=-0.001801640704;
    coefs[13]= 0.000353713800;
    break;
  case 8:
    coefs=new double[16];
    coefs[0] = 0.054415842243;
    coefs[1] = 0.312871590914;
    coefs[2] = 0.675630736297;
    coefs[3] = 0.585354683654;
    coefs[4] =-0.015829105256;
    coefs[5] =-0.284015542962;
    coefs[6] = 0.000472484574;
    coefs[7] = 0.128747426620;
    coefs[8] =-0.017369301002;
    coefs[9] =-0.044088253931;
    coefs[10]= 0.013981027917;
    coefs[11]= 0.008746094047;
    coefs[12]=-0.004870352993;
    coefs[13]=-0.000391740373;
    coefs[14]= 0.000675449406;
    coefs[15]=-0.000117476784;
    break;
  case 9:
    coefs=new double[18];
    coefs[0] = 0.038077947364;
    coefs[1] = 0.243834674613;
    coefs[2] = 0.604823123690;
    coefs[3] = 0.657288078051;
    coefs[4] = 0.133197385825;
    coefs[5] =-0.293273783279;
    coefs[6] =-0.096840783223;
    coefs[7] = 0.148540749338;
    coefs[8] = 0.030725681479;
    coefs[9] =-0.067632829061;
    coefs[10]= 0.000250947115;
    coefs[11]= 0.022361662124;
    coefs[12]=-0.004723204758;
    coefs[13]=-0.004281503682;
    coefs[14]= 0.001847646883;
    coefs[15]= 0.000230385764;
    coefs[16]=-0.000251963189;
    coefs[17]= 0.000039347320;
    break;
  case 10:
    coefs=new double[20];
    coefs[0] = 0.026670057901;
    coefs[1] = 0.188176800078;
    coefs[2] = 0.527201188932;
    coefs[3] = 0.688459039454;
    coefs[4] = 0.281172343661;
    coefs[5] =-0.249846424327;
    coefs[6] =-0.195946274377;
    coefs[7] = 0.127369340336;
    coefs[8] = 0.093057364604;
    coefs[9] =-0.071394147166;
    coefs[10]=-0.029457536822;
    coefs[11]= 0.033212674059;
    coefs[12]= 0.003606553567;
    coefs[13]=-0.010733175483;
    coefs[14]= 0.001395351747;
    coefs[15]= 0.001992405295;
    coefs[16]=-0.000685856695;
    coefs[17]=-0.000116466855;
    coefs[18]= 0.000093588670;
    coefs[19]=-0.000013264203;
    break;
  default:
    coefs=0;
  }
  return coefs;
}

  // The next 2 functions are the only ones needed to call.

void DaubTrans(int j, ValErr_t *vector, int p) {
  // Forward Daubechies wavelet transform of vector (size 2^j).
  // p vanishing moments.
  double *coefs=DaubCoefs(p);
  if (coefs)
    DaubTransc(j,vector,p,coefs);
  delete coefs;
}

void InvDaubTrans(int j, ValErr_t *vector, int p) {
  // Inverse Daubechies wavelet transform of vector (size 2^j).
  // p vanishing moments.
  double *coefs=DaubCoefs(p);
  if (coefs)
    InvDaubTransc(j,vector,p,coefs);
  delete coefs;
}
 
void THistoCalc::Smooth(int p, double epsilon,int nfreq) {
  
  int pow2 = 0 ;
  
  // Calcul de la longueur du vecteur en puissance de 2
  while ((1<<pow2) < GetNbinsX())
    pow2 += 1 ;

  // Recopie des valeurs de l'histo dans le vecteur
  ValErr_t *y = new ValErr_t[1<<pow2];  
  for (int i=0; i<GetNbinsX(); i++)
    {y[i].fVal = GetBinContent(i+1);
	 y[i].fErr = GetBinError(i+1);}
  // On complète entre GetNbinsX() et 2^pow2 par des 0.
  for (int  i=GetNbinsX()+1; i<(1<<pow2); i++)
    {y[i].fVal=0;
     y[i].fErr=0;}
  // Forward Daubechies
  DaubTrans(pow2,y,p);
  
  //Suppression de tous les coeffs inferieurs a epsilon
  int k=0;

/*  ofstream ofic("Freq.dat");

  for (int i=0;i<(1<<pow2);i++)	
      ofic << i+1 << " " <<fabs(y[i].fVal)<< endl;
      

  for (int i=pow2;i>pow2-nfreq;i--)
	for (int j= 1<<i ;j>(1<<(i-1)) ;j--)
		{y[j-1].fVal= 0 ;y[j-1].fErr = 0 ;}
	


for (int i=pow2;i>pow2-nfreq;i--)
	for (int j= 1<<i ;j>(1<<(i-1)) ;j--)
	  if (fabs(y[j].fVal) < epsilon) 
	    {y[j].fVal= 0 ;y[j].fErr = 0 ;k++;}


      ofic << " & " << endl;
  for (int i=0;i<(1<<pow2);i++)	
      ofic << i+1 << " " <<fabs(y[i].fVal)<< endl;
*/      

  for (int i=0;i<(1<<pow2);i++)
    if (fabs(y[i].fVal) < epsilon  )
      {y[i].fVal= 0 ;y[i].fErr = 0 ;k++;}

  cout << k << " coeffs supprimés / " << (1<<pow2) << " en tout\n";
  
  // Inverse Daubechies
  InvDaubTrans(pow2,y,p);
  
  for (int i=0; i<GetNbinsX(); i++)
    SetBinValErr(i+1,y[i]) ;

  delete y;
}


ValErr_t THistoCalc::Convole(THistoCalc& h, Axis_t t) {
  
  Int_t i1,i2;
  // low[i1]<=t1<up[i1]
  // low[i2]<=t-t1<up[i2]
  Axis_t t1;
  ValErr_t valeur;
  TAxis *axe1=GetXaxis();
  Int_t n1=GetNbinsX();
  TAxis *axe2=h.GetXaxis();
  Int_t n2=h.GetNbinsX();

  valeur=0;

  //initialisation des indices

  // si i?==n?+1 ou i?=0, on considere que h?[i?]==0

  i1=1;
  t1=axe1->GetBinLowEdge(i1);

  if (axe2->GetBinLowEdge(1)>t-t1)
    return 0;

  i2=1;
  while(i2<=n2 && axe2->GetBinLowEdge(i2)<=t-t1)
    i2++;
  i2--;
  if (axe2->GetBinLowEdge(n2)<=t-t1)
    i2=n2;
  if (axe2->GetBinUpEdge(n2)<=t-t1)
    i2=n2+1;

  // ici i1=1 et i2 entre 1 et n2+1

  // boucle hors de la definition de h : le temps t1 augmente sans contribuer a valeur.
  while(i2==n2+1 && i1<=n1) {
    if (axe1->GetBinUpEdge(i1) <= t-axe2->GetBinLowEdge(i2)) {
      // On fait avancer i1
      i1++;
      t1=axe1->GetBinLowEdge(i1);
    }
    else {
      // On fait reculer i2
      t1=t-axe2->GetBinLowEdge(i2);
      i2--;
    }    
  }

  if (i1==n1+1)
    return 0;

  // boucle dedans : valeur est incrementee.
  while(i1<=n1 && i2>=0) {
    if (axe1->GetBinUpEdge(i1) <= t-axe2->GetBinLowEdge(i2)) {
      // On fait avancer i1
      valeur+=(axe1->GetBinUpEdge(i1)-t1)*GetBinValErr(i1)*h.GetBinValErr(i2);
      i1++;
      t1=axe1->GetBinLowEdge(i1);
    }
    else {
      // On fait reculer i2
      valeur+=(t-axe2->GetBinLowEdge(i2)-t1)*GetBinValErr(i1)*h.GetBinValErr(i2);
      t1=t-axe2->GetBinLowEdge(i2);
      i2--;
    }
  }

  return valeur;

}

void THistoCalc::HeavisideConvole() {

  TAxis *axe=GetXaxis();
  for (Int_t i=1;i<=GetNbinsX();i++)
    SetBinValErr(i,Integrale(GetBinCenter(i),axe->GetBinUpEdge(GetNbinsX())));
  CenterToMean();

}

void THistoCalc::HeavisideConvole(ValErr_t integrale) {

  TAxis *axe=GetXaxis();
  for (Int_t i=GetNbinsX();i>=1;i--)
    SetBinValErr(i,integrale-Integrale(axe->GetBinLowEdge(1),GetBinCenter(i)));
  CenterToMean();

}

THistoCalc THistoCalc::Convole(THistoCalc& h) {
  
  THistoCalc result(*this);

  TAxis *axe=GetXaxis();

  for (Int_t i=1;i<=GetNbinsX();i++)
    result.SetBinValErr(i,Convole(h,axe->GetBinCenter(i)));

  result.CenterToMean();

  return result;

}

void THistoCalc::logreptemps()

{   ValErr_t val ;
    Int_t n=GetNbinsX();
	
 
    Int_t k = 1 ;
    while (k <= n && GetBinContent(k) <= 0 )	
     	k++;
    double min = GetBinContent(k);
      for (Int_t i=1; i<=n;i++)
           if (GetBinContent(i)>0. && GetBinContent(i) < min )
	   min = GetBinContent(i); 
	cout << " Fonction logreptemps valeur min = " <<  min << endl ;
      for (Int_t i=1; i<=n;i++)
	{if (GetBinContent(i)<=0)
		SetBinValErr(i, log(min)  );
	 else
		{val = GetBinValErr(i);
		 SetBinContent(i,log(val.Valeur()));
 		 SetBinError(i,val.Erreur()/val.Valeur());
		}
	}
} 


// Remove the intrinsic source and take the logarithmic value
void THistoCalc::LogRemoveSI(double SI_min,double SI_max)

{	TAxis *axe=GetXaxis();
	int bin_min = axe->FindBin(SI_min);
	int bin_max = axe->FindBin(SI_max);

	Int_t bin= bin_min;
	ValErr_t Integrale;
	Int_t compt = 0;	

	while (bin<bin_max)
		{Integrale += GetBinValErr(bin);
		 bin++;
		 compt++;
		}
	Integrale /= compt;
	cout << "Intrinsic source between " << SI_min << " and " << SI_max << " = " << Integrale << endl ;
	(*this) += -Integrale;
	this->logreptemps();
}	

// Remove the intrinsic source 
void THistoCalc::RemoveSI(double SI_min,double SI_max)

{	TAxis *axe=GetXaxis();
	int bin_min = axe->FindBin(SI_min);
	int bin_max = axe->FindBin(SI_max);

	Int_t bin= bin_min;
	ValErr_t Integrale;
	Int_t compt = 0;	

	while (bin<bin_max)
		{Integrale += GetBinValErr(bin);
		 bin++;
		 compt++;
		}
	Integrale /= compt;
	cout << "Intrinsic source between " << SI_min << " and " << SI_max << " = " << Integrale << endl ;
	(*this) += -Integrale;

}	

void THistoCalc::Oppose() {
  for (int i=1;i<GetNbinsX();i++)
    SetBinValErr(i,-GetBinValErr(i));
}

THistoCalc THistoCalc::RemoveBins(double tCut)

{
	
	TAxis *axe=GetXaxis();
	int ind=0;	
	int bincut = FindBin(tCut);
	Double_t tbin[axe->GetNbins()-axe->FindBin(tCut)+1];
	for (int i=axe->FindBin(tCut);i<=axe->GetNbins();i++)
		{tbin[ind]=axe->GetBinLowEdge(i);
		ind++;
		}			
	tbin[ind]=axe->GetBinUpEdge(axe->GetNbins());

	THistoCalc result(axe->GetNbins()-axe->FindBin(tCut),tbin);
	axe=result.GetXaxis();	
	for (int i=0;i<axe->GetNbins();i++)
		{result.SetBinValErr(i+1,this->GetBinValErr(i+bincut));
		}

	return result;
}

   THistoCalc THistoCalc::RemoveBins(double tCutLow,double tCutHigh)

{
	
	TAxis *axe= this->GetXaxis();
	int ind=0;	
	int bincutLow = axe->FindBin(tCutLow);
//	int bincutHigh = axe->FindBin(tCutHigh);

	Double_t tbin[axe->FindBin(tCutHigh)-axe->FindBin(tCutLow)+1];

	for (int i=axe->FindBin(tCutLow);i<=axe->FindBin(tCutHigh);i++)
		{tbin[ind]=axe->GetBinLowEdge(i);
		ind++;
		}			
	tbin[ind]=axe->GetBinUpEdge(axe->FindBin(tCutLow));

	THistoCalc result(axe->FindBin(tCutHigh)-axe->FindBin(tCutLow),tbin);
	axe=result.GetXaxis();	
	for (int i=0;i<axe->GetNbins();i++)
		{result.SetBinValErr(i+1,this->GetBinValErr(i+bincutLow));
		}

	return result;
}

THistoCalc& THistoCalc::operator += (THistoCalc& h) {

  Int_t i,j;
  Axis_t temps;
  TAxis *axe1=GetXaxis();
  TAxis *axe2=h.GetXaxis();
  Int_t n1=GetNbinsX();
  Int_t n2=h.GetNbinsX();

  i=j=1;

  while(j<=n2 && axe2->GetBinUpEdge(j)<axe1->GetBinLowEdge(1))
    j++;
  if (j>n2)
    return *this;

  while(i<=n1 && axe1->GetBinUpEdge(i)<axe2->GetBinLowEdge(1))
    i++;
  if (i>n1)
    return *this;

  temps=axe1->GetBinLowEdge(i);
  if (axe2->GetBinLowEdge(j)>temps)
    temps=axe2->GetBinLowEdge(j);
  
  while(i<=n1 && j<=n2) {
    if (axe2->GetBinUpEdge(j) <= axe1->GetBinUpEdge(i)) {
      AddBinValErr(i,h.GetBinValErr(j)*(axe2->GetBinUpEdge(j)-temps)
		   /(axe1->GetBinUpEdge(i)-axe1->GetBinLowEdge(i)));
      temps = axe2->GetBinUpEdge(j);
      j++;
    }
    else {
      AddBinValErr(i,h.GetBinValErr(j)*(axe1->GetBinUpEdge(i)-temps)
		   /(axe1->GetBinUpEdge(i)-axe1->GetBinLowEdge(i)));
      temps = axe1->GetBinUpEdge(i);
      i++;
    }
  }
  return *this;
}

THistoCalc& THistoCalc::operator -= (THistoCalc& h) {
  h.Oppose();
  (*this)+=h;
  h.Oppose();
  return (*this);
}

THistoCalc& THistoCalc::operator *= (ValErr_t factor) {
  for (Int_t i=1;i<=GetNbinsX();i++) {
    SetBinValErr(i,GetBinValErr(i)*factor);
  }
  return *this;
}

THistoCalc& THistoCalc::operator += (ValErr_t constant) {
  for (Int_t i=1;i<=GetNbinsX();i++) {
    SetBinValErr(i,GetBinValErr(i)+constant);
  }
  return *this;
}

THistoCalc& THistoCalc::operator *= (THistoCalc& h) {

  // Rebinning.
  THistoCalc r=*this;
  for (Int_t i=1;i<=GetNbinsX();i++)
    r.SetBinValErr(i,0);
  r+=h;

  // Division bin a bin
  for (Int_t i=1;i<=GetNbinsX();i++)
    SetBinValErr(i,GetBinValErr(i)*r.GetBinValErr(i));
  return *this;

}

THistoCalc& THistoCalc::operator *= (const TSpline& Spline) {
  TAxis *Axe=GetXaxis();
  Axis_t xmin,xmax;
  Double_t y1,y2;
  xmax=Axe->GetBinLowEdge(1);
  for (Int_t i=1;i<=Axe->GetNbins();i++) {
    xmin=xmax;
    y1=Spline.Eval(xmin);
    xmax=Axe->GetBinUpEdge(i);
    y2=Spline.Eval(xmax);
    SetBinValErr(i,GetBinValErr(i)*ValErr_t((y1+y2)/2,fabs(y2-y1)/2));
  }
  return *this;
}

THistoCalc& THistoCalc::operator /= (THistoCalc& h) {

  // Rebinning.
  THistoCalc r=*this;
  for (Int_t i=1;i<=GetNbinsX();i++)
    r.SetBinValErr(i,0);
  r+=h;

  ValErr_t ve;
  // Division bin a bin
  for (Int_t i=1;i<=GetNbinsX();i++) {
    ve=r.GetBinValErr(i);
    if (ve.fVal)
      SetBinValErr(i,GetBinValErr(i)/ve);
    else
      SetBinValErr(i,0);
  }
  return *this;

}

THistoCalc operator + (const THistoCalc& x, THistoCalc& y) {THistoCalc r=x;r+=y;return r;}
THistoCalc operator * (const THistoCalc& x,Double_t factor) {THistoCalc r=x;r*=factor;return r;}
THistoCalc operator * (const THistoCalc& x, THistoCalc& y) {THistoCalc r=x;r*=y;return r;}
THistoCalc operator * (const THistoCalc& x,const TSpline& Spline) {THistoCalc r=x;r*=Spline;return r;}
THistoCalc operator * (const THistoCalc& x,TGraph& Graph) {THistoCalc r=x;r*=Graph;return r;}
THistoCalc operator / (const THistoCalc& x, THistoCalc& y) {THistoCalc r=x;r/=y;return r;}
