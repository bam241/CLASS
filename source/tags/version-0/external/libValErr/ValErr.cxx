#include "ValErr.hxx"
#include <cmath>


////////////////////////////////////////////////////////////////
//                       ValErr_t




ValErr_t& __doapl (ValErr_t* ths, const ValErr_t& r) {

  ths->fVal+=r.fVal;
  ths->fErr=sqrt(ths->fErr*ths->fErr+r.fErr*r.fErr);
  // ths->fErr=ths->fErr+r.fErr;
  return *ths;

}

ValErr_t& __doapl (ValErr_t* ths, const double& d) {

  ths->fVal+=d;
  return *ths;

}

ValErr_t& __doami (ValErr_t* ths, const ValErr_t& r) {

  ths->fVal-=r.fVal;
  ths->fErr=sqrt(ths->fErr*ths->fErr+r.fErr*r.fErr);
  // ths->fErr=ths->fErr+r.fErr;
  return *ths;

}

ValErr_t& __doami (ValErr_t* ths, const double& d) {

  ths->fVal-=d;
  return *ths;

}

ValErr_t& __doaml (ValErr_t* ths, const ValErr_t& r) {
  
  double tmpVal=ths->fVal*r.fVal;
  ths->fErr=ths->fErr*fabs(r.fVal)+r.fErr*fabs(ths->fVal);
  ths->fVal=tmpVal;
  return *ths;
  
}

ValErr_t& __doaml (ValErr_t* ths, const double& d) {
  
  ths->fVal*=d;
  ths->fErr*=fabs(d);
  return *ths;
  
}

ValErr_t& __doadi (ValErr_t* ths, const ValErr_t& r) {

	if(r.fVal>0)
	{
		ths->fErr=ths->fErr/fabs(r.fVal)+r.fErr*fabs(ths->fVal)/r.fVal/r.fVal;
  		ths->fVal/=r.fVal;
	}
  return *ths;

}

ValErr_t& __doadi (ValErr_t* ths, const double& d) {
  ths->fVal/=d;
  ths->fErr/=fabs(d);
  return *ths;
}

ValErr_t& __doaop (ValErr_t* ths) {
  
  ths->fVal=-ths->fVal;
  return *ths;

}

ValErr_t operator + (const ValErr_t& x, const ValErr_t& y) {ValErr_t r=x;r+=y;return r;}
ValErr_t operator + (const ValErr_t& x, const double& y) {ValErr_t r=x;r+=y;return r;}
ValErr_t operator + (const double& x, const ValErr_t& y) {return y+x;}
ValErr_t operator - (const ValErr_t& x) {ValErr_t r=x; return __doaop(&r);}
ValErr_t operator - (const ValErr_t& x, const ValErr_t& y) {ValErr_t r=x;r-=y;return r;}
ValErr_t operator - (const ValErr_t& x, const double& y) {ValErr_t r=x;r-=y;return r;}
ValErr_t operator - (const double& x, const ValErr_t& y) {ValErr_t r=-y;r+=x;return r;}
ValErr_t operator * (const ValErr_t& x, const ValErr_t& y) {ValErr_t r=x;r*=y;return r;}
ValErr_t operator * (const ValErr_t& x, const double& y) {ValErr_t r=x;r*=y;return r;}
ValErr_t operator * (const double& x, const ValErr_t& y) {ValErr_t r=y;r*=x;return r;}
ValErr_t operator / (const ValErr_t& x, const ValErr_t& y) {ValErr_t r=x;r/=y;return r;}
ValErr_t operator / (const ValErr_t& x, const double& y) {ValErr_t r=x;r/=y;return r;}
bool operator == (const ValErr_t& x, const ValErr_t& y) {return x.fVal==y.fVal && x.fErr==y.fErr;}
bool operator != (const ValErr_t& x, const ValErr_t& y) {return x.fVal!=y.fVal || x.fErr!=y.fErr;}

istream& operator >> (istream& is, ValErr_t& x) {

  double val, err = 0;
  char ch = 0;
  
  if (is) {
      is >> val;
      is >> ch;
      if (ch == '+') {
	is >> ch;
	if (ch=='-')
	  is >> err;
      }
  }
  
  if (ch != 0 && ch != '-')
    is.setstate (ios::failbit);
  else if (is.good ())
    x = ValErr_t(val,err);
  
  return is;

}

ostream& operator << (ostream& os, const ValErr_t& x) {
  
  return os << x.fVal << "+-" << x.fErr;

}
