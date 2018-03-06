#include "ZAIMass.hxx"

//#include "CLASSConstante.hxx"

#include <fstream>
#include <string>
#include "stdlib.h"

#include "IsotopicVector.hxx"
//________________________________________________________________________
//
//		ZAI
//
//
//
//
//________________________________________________________________________
//____________________________InClass Operator____________________________
//________________________________________________________________________

ZAIMass::ZAIMass() {
  string CLASSPATH = getenv("CLASS_PATH");
  string MassDataFile = CLASSPATH + "/data/Mass.dat";

  ifstream infile(MassDataFile.c_str());

  if (!infile.good()) {
    cout << " ZAIMass Error.\n can't find/open file " << MassDataFile << endl;
    exit(1);
  }

  int Z, A;
  string Name;
  double MassUnity, MassDec, error;
  while (infile >> Z >> A >> Name >> MassUnity >> MassDec >> error) {
    double Masse = MassUnity + MassDec * 1e-6;
    fZAIMass.insert(pair<ZAI, double>(ZAI(Z, A, 0), Masse));
  }

  infile.close();
}

ZAIMass::~ZAIMass() { fZAIMass.clear(); }

double ZAIMass::GetMass(ZAI zai) const {
  if (zai.A() < 0 || zai.Z() < 0) return 0;  // If isotopes -1, -2, -3 return 0

  map<ZAI, double>::const_iterator MassIT =
      fZAIMass.find(ZAI(zai.Z(), zai.A(), 0));

  if (MassIT == fZAIMass.end())
    return zai.A();
  else
    return MassIT->second;
}

double ZAIMass::GetMass(const IsotopicVector& IV) const {
  double AVOGADRO = 6.02214129e23;
  double TotalMass = 0;

  for (IsotopicVector::const_iterator it = IV.begin(); it != IV.end(); it++) {
    TotalMass += it->second / AVOGADRO * GetMass(it->first);
  }

  return TotalMass * 1e-6;
}
