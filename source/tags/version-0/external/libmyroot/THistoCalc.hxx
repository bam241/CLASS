#ifndef _THistoCalc_HXX_
#define _THistoCalc_HXX_

#include <TROOT.h>
#include <TH1.h>
#include <TGraphErrors.h>
#include "NextName.hxx"
#include <libValErr/ValErr.hxx>
#include <TSpline.h>

// La valeur du bin est supposee etre la valeur moyenne ! (pas integrale)


class THistoCalc : public TH1D {
  
public :
  
  THistoCalc() : TH1D() {} // Constructeur par defaut
  THistoCalc(const THistoCalc& histocalc) : // Constructeur par recopie
    TH1D(histocalc) {SetName(NextName());}
  THistoCalc(TH1F& histo); // Constructeur par recopie
  // Constructeur avec bins de taille constante
  THistoCalc(Int_t nbinsx, Axis_t xlow, Axis_t xup) :
    TH1D(NextName(),"",nbinsx,xlow,xup) {}
  // Constructeur avec bins Float_t de taille variable
  THistoCalc(Int_t nbinsx, Float_t *xbins) :
    TH1D(NextName(),"",nbinsx,xbins) {}
  // Constructeur avec bins Double_t de taille variable
  THistoCalc(Int_t nbinsx, Double_t *xbins) :
    TH1D(NextName(),"",nbinsx,xbins) {}

  ValErr_t GetBinValErr(Int_t bin);
  void SetBinValErr(Int_t bin, const ValErr_t& ve);
  void AddBinValErr(Int_t bin, const ValErr_t& ve) {SetBinValErr(bin,GetBinValErr(bin)+ve);}

  void Write(ostream& s); // Ecrit x y dy dans la stream s.

  // Ajoute une erreur proportionnelle a la difference entre deux bins succesifs :
  // erreur due a la taille des bins. In place.
  void AddBinSizeErr();

  // Suppose le contenu des bins est la valeur au centre.
  // Le remplace par la valeur moyenne sur la largeur du bin (in place).
  // (diff. si derivée second non nulle).
  void CenterToMean();

  // Somme des bins pondérée par la largeur.
  // Ne tient pas compte des bins underflow et overflow.
  ValErr_t Integrale(Double_t xmin,Double_t xmax); // *MENU*
  ValErr_t LogIntegrale(Double_t xmin,Double_t xmax); // *MENU*

  // Somme des bins ponderes par la largeur et la valeur en abscisse au milieu
  // du bin 

  ValErr_t First_Moment(Double_t xmin,Double_t xmax);
  
  // Fit exponentiel sur nbins bins pour t = le centre de chaque bin,
  // et renvoie un TGraphErrors
  TGraphErrors *Omega(Double_t tmin, Double_t tmax, Int_t nbins);
  TGraphErrors *OmegaGlissant(Double_t tmin, Double_t tmax, Int_t nbins);
  

  // fit polynomiale glissant sur un nombre de bin renvoie la valeur negative
  // de la derive au milieu du bin
  TGraphErrors *AlphaPoly(Double_t tmin, Double_t tmax, Int_t nbins); 
  TGraphErrors *AlphaPolyGlissant(Double_t tmin,Double_t tmax,Int_t nbinfit);
  THistoCalc   *DiffBin();
  
  // Smooth avec des ondelettes un histo sans calcul d erreurs !!!!
  // p=nombre de moments annulés par les ondelettes (2 à 10)
  // epsilon=seuil sous lequel les coefs sont supprimés.
  void Smooth(int p, double epsilon,int nfreq=0);
  void Smooth(Int_t i) {TH1::Smooth(i);};
  // Valeur de la convolution avec h en t.
  ValErr_t Convole(THistoCalc& h, Axis_t t);
  
  // Convolue de this et h, avec le binning de this. Pas in place.
  THistoCalc Convole(THistoCalc& h);

  // Convolue de this et une fonction de Heaviside (in place)
  void HeavisideConvole();
  // Idem sans connaitre forcément bien la fin mais avec la
  // valeur de l'intégrale du début du binning à +infty
  void HeavisideConvole(ValErr_t integrale);
  
  // remplace le contenu de chaque bin par
  // le log du contenu si positif
  void logreptemps(); 

  // Remove the intrinsic source for PNS Histo and log of each bin
  void LogRemoveSI(double SI_min,double SI_max);

  // Remove the intrinsic source for PNS Histo
  void RemoveSI(double SI_min,double SI_max);
  
  // Replace each value by its opposite
  void Oppose();

  // Remove bins before tCut
  THistoCalc RemoveBins(double tCut);
  // Remove bins before tCutLow and after tCutHigh
   THistoCalc RemoveBins(double tCutLow,double tCutHigh);
  
  // Addition d'un autre histo a celui-ci
  // ATTENTION, des methodes Add sont definies dans TH1 ! in place
  THistoCalc& operator += (THistoCalc& h);

  // Soustraction d'un autre histo a celui-ci
  THistoCalc& operator -= (THistoCalc& h);

  // Multiplication par un facteur constant
  THistoCalc& operator *= (ValErr_t factor);

  // Addition d'une constante
  THistoCalc& operator += (ValErr_t constant);

  // Multiplication par un spline, erreurs prises en compte :
  // erreurs de this, largeur de bin de this. in place
  THistoCalc& operator *= (const TSpline& Spline);

  // Multiplication par un autre histo.
  // l'autre est rebinne, puis multiplication bin par bin. in place
  THistoCalc& operator *= (THistoCalc& h);

  // Multiplication par un TGraph, erreurs prises en compte :
  // erreurs de this, largeur de bin de this.
  THistoCalc& operator *= (TGraph& Graph)
  {
  	TSpline3 *spleen=new TSpline3("",&Graph);
	*this *= *spleen;
 	 return *this;
  }

  // Division par un autre histo.
  // l'autre est rebinne, puis division bin par bin. in place
  THistoCalc& operator /= (THistoCalc& h);
  

};

// Le binning est celui de x
THistoCalc operator + (THistoCalc& x, const THistoCalc& y);
THistoCalc operator * (const THistoCalc& x, Double_t factor);
THistoCalc operator * (const THistoCalc& x, THistoCalc& y);
THistoCalc operator * (const THistoCalc& x, const TSpline& Spline);
THistoCalc operator * (const THistoCalc& x, TGraph& Graph);
THistoCalc operator / (const THistoCalc& x, THistoCalc& y);

#endif
