#ifndef _ValErr_HXX_
#define _ValErr_HXX_

#include <iostream>
using namespace std;
class ValErr_t;

// Implementation reelle des calculs
// additions
ValErr_t& __doapl (ValErr_t* ths, const ValErr_t& r); // somme des variances (supp indep)
ValErr_t& __doapl (ValErr_t* ths, const double& d); // erreur inchangee
// soustractions
ValErr_t& __doami (ValErr_t* ths, const ValErr_t& r); // somme des variances (supp indep)
ValErr_t& __doami (ValErr_t* ths, const double& d); // erreur inchangee
// multiplications
ValErr_t& __doaml (ValErr_t* ths, const ValErr_t& r); // 'a la derivee'
ValErr_t& __doaml (ValErr_t* ths, const double& d); // erreur multipliee
// divisions
ValErr_t& __doadi (ValErr_t* ths, const ValErr_t& r); // 'a la derivee'
ValErr_t& __doadi (ValErr_t* ths, const double& d); // erreur divisee
// opposition
ValErr_t& __doaop (ValErr_t* ths); // erreur inchangee

class ValErr_t {
  
public :
  
  double fVal; // Moyenne
  double fErr; // Ecart-type
  
  ValErr_t(const double Val=0, const double Err=0) {fVal=Val;fErr=Err;}
  ValErr_t(const ValErr_t& ve) {fVal=ve.fVal;fErr=ve.fErr;}
  ValErr_t& operator = (const ValErr_t& ve) {fVal=ve.fVal;fErr=ve.fErr; return *this;}
  ValErr_t& operator = (const double& d) {fVal=d;fErr=0;return *this;}
  ValErr_t& operator += (const ValErr_t& ve) {return __doapl(this,ve);}
  ValErr_t& operator += (const double& d) {return __doapl(this,d);}
  ValErr_t& operator -= (const ValErr_t& ve) {return __doami(this,ve);}
  ValErr_t& operator -= (const double& d) {return __doami(this,d);}
  ValErr_t& operator *= (const ValErr_t& ve) {return __doaml(this,ve);}
  ValErr_t& operator *= (const double& d) {return __doaml(this,d);}
  ValErr_t& operator /= (const ValErr_t& ve) {return __doadi(this,ve);}
  ValErr_t& operator /= (const double& d) {return __doadi(this,d);}
  
  double Valeur() {return fVal;}
  double Erreur() {return fErr;}

};

ValErr_t operator + (const ValErr_t& x, const ValErr_t& y);
ValErr_t operator + (const ValErr_t& x, const double& y);
ValErr_t operator + (const double& x, const ValErr_t& y);
ValErr_t operator - (const ValErr_t& x);
ValErr_t operator - (const ValErr_t& x, const ValErr_t& y);
ValErr_t operator - (const ValErr_t& x, const double& y);
ValErr_t operator - (const double& x, const ValErr_t& y);
ValErr_t operator * (const ValErr_t& x, const ValErr_t& y);
ValErr_t operator * (const ValErr_t& x, const double& y);
ValErr_t operator * (const double& x, const ValErr_t& y);
ValErr_t operator / (const ValErr_t& x, const ValErr_t& y);
ValErr_t operator / (const ValErr_t& x, const double& y);
bool operator == (const ValErr_t& x, const ValErr_t& y);
bool operator != (const ValErr_t& x, const ValErr_t& y);
istream& operator >> (istream& is, ValErr_t& x);
ostream& operator << (ostream& os, const ValErr_t& x);

#endif
