#ifndef _CLASSConstante_
#define _CLASSConstante_

//CLASS library
#include "ZAIMass.hxx"
#include "ZAIHeat.hxx"
#include "ZAITox.hxx"
#include "DecayDataBank.hxx"

typedef long long int cSecond;

static const double AVOGADRO = 6.02214129e23;	// Avogadro Number [1/mol]

static const ZAIMass cZAIMass;					// Mass list of all nuclei stored in [g/mol]
static const ZAIHeat cZAIHeat;					// Thermal power list of all nuclei  [W/nucleus]
static const ZAITox  cZAITox;

static const cSecond cYear = 3600*24*365.25;	// Seconds in a year
static DecayDataBank cDecayData;	// Seconds in a year


#endif
