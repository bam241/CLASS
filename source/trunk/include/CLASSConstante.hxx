#ifndef _CLASSConstante_
#define _CLASSConstante_

//CLASS library
#include "ZAIMass.hxx"
#include "ZAIHeat.hxx"
#include "ZAITox.hxx"

typedef long long int cSecond;

const double AVOGADRO = 6.02214129e23;	// Avogadro Number [1/mol]

const ZAIMass cZAIMass;					// Mass list of all nuclei stored in [g/mol]
const ZAIHeat cZAIHeat;					// Thermal power list of all nuclei  [W/nucleus]
const ZAITox  cZAITox;

const cSecond cYear = 3600*24*365.25;	// Seconds in a year


#endif
