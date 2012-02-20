#ifndef _THisto2Calc_HXX_
#define _THisto2Calc_HXX_

#include <TROOT.h>
#include <TH2.h>
#include <TGraphErrors.h>
#include "NextName.hxx"
#include <libValErr/ValErr.hxx>
#include <TSpline.h>

// Histo 2D
// La valeur du bin est supposee etre la valeur au centre ! (pas integrale)


class THisto2Calc : public TH2D {
  
public :
  
  THisto2Calc() : TH2D() {} // Constructeur par defaut
  THisto2Calc(const THisto2Calc& histo2calc) : // Constructeur par recopie
    TH2D(histo2calc) {SetName(NextName());}
  // Constructeur avec bins de taille constante en x et y
  THisto2Calc(Int_t nbinsx, Axis_t xlow, Axis_t xup,
	      Int_t nbinsy, Axis_t ylow, Axis_t yup):
    TH2D(NextName(),"",nbinsx,xlow,xup,nbinsy,ylow,yup) {}
  // Constructeur avec bins Double_t en x, taille constante en y
  THisto2Calc(Int_t nbinsx, const Double_t* xbins,
	      Int_t nbinsy, Axis_t ylow, Axis_t yup):
    TH2D(NextName(),"",nbinsx,xbins,nbinsy,ylow,yup) {}
  // Constructeur avec bins Double_t en y, taille constante en x
  THisto2Calc(Int_t nbinsx, Axis_t xlow, Axis_t xup,
	      Int_t nbinsy, const Double_t* ybins):
    TH2D(NextName(),"",nbinsx,xlow,xup,nbinsy,ybins) {}
  // Constructeur avec bins Double_t de taille variable
  THisto2Calc(Int_t nbinsx, const Double_t* xbins,
	      Int_t nbinsy, const Double_t* ybins):
    TH2D(NextName(),"",nbinsx,xbins,nbinsy,ybins) {}
  // Constructeur avec bins Float_t de taille variable
  THisto2Calc(Int_t nbinsx, const Float_t* xbins,
	      Int_t nbinsy, const Float_t* ybins):
    TH2D(NextName(),"",nbinsx,xbins,nbinsy,ybins) {}

  ValErr_t GetBinValErr(Int_t binx, Int_t biny);
  void SetBinValErr(Int_t binx, Int_t biny, const ValErr_t& ve);
  void AddBinValErr(Int_t binx, Int_t biny, const ValErr_t& ve)
  {SetBinValErr(binx,biny,GetBinValErr(binx,biny)+ve);}

  void Write(ostream& s); // Ecrit x y z dz dans la stream s.

};

#endif
