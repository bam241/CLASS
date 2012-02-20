#include "THisto2Calc.hxx"
#include <TF1.h>
#include <iostream>
#include <math.h>
#include <TSpline.h>

using namespace std;


////////////////////////////////////////////////////////////////////
//                           THisto2Calc






ValErr_t THisto2Calc::GetBinValErr(Int_t binx, Int_t biny) {

  return ValErr_t(GetBinContent(binx,biny),GetBinError(binx,biny));

}

void THisto2Calc::SetBinValErr(Int_t binx, Int_t biny, const ValErr_t& ve) {

  SetBinContent(binx,biny,ve.fVal);
  SetBinError(binx,biny,ve.fErr);

}

void THisto2Calc::Write(ostream& s) {

  s << "@type xyz\n";
  for(Int_t i=1 ; i<=GetNbinsX(); i++)
    for(Int_t j=1 ; j<=GetNbinsY(); j++) {
      s << GetXaxis()->GetBinCenter(i) << " "
	<< GetYaxis()->GetBinCenter(j) << " "
	<< GetBinContent(i,j) //<< " "
	//<< GetBinError(i,j)
	<< endl;
    }
  s << "&\n";

}
