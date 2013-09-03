#ifndef _ZAI_
#define _ZAI_ 

/*!
 \file
 \brief Header file for ZAI classes. 
 Define a nuclei as : Z A I.

 
 @author BaM
 @version 2.0
 */

#include <string>
#include "TObject.h"
#include <iostream>

using namespace std;




///< A ZAI defined for a Nucleus. 
/*!

*/
class ZAI : public TObject
{
	
	
 public:
	///< Default constructor
	ZAI();
 	///< Normal Constructor. 
 	/*!
 	 Default: No parent
 	 \param Z : number of protons
 	 \param A : number of nucleons (A=0 means natural isotopes) 
 	*/
	ZAI(int Z, int A, int I=0);

	///< Normal Destructor.
	~ZAI(); 

	/*!
	\name ZAI main attributes 
	*/
	//@{
	int  Z()const{ return fZ; } //!< returns the number of protons
	int  A()const{ return fA; } //!< returns the number of nucleons
	int  I()const{ return fI; } //!< returns the Isomeric State
	int  N()const{ return fA-fZ; } //!< returns the number of neutrons
		
	void SetMass(double m)	{ fMass=m; }	///< set the mass of a ZAI
	double GetMass();			///< get the mass of a ZAI


	ZAI& operator=(ZAI& IVa);		//!< ...
	bool operator <(const ZAI& zai) const	{ return (fZ != zai.Z())?  
							(fZ < zai.Z()) : ( (fA != zai.A())?
								 (fA < zai.A()) : (fI < zai.I()) ); }
						//!< ZAI Comparator
	bool operator !=(const ZAI& zai) const	{ return ( fZ != zai.Z() ) || ( fA != zai.A() ) || ( fI < zai.I() ); }
						//!< ZAI Comparator
	void Print() const	{ cout << fZ << " " << fA << " " << fI << endl;}
	protected :
 	
 	string 	fName;		///< Name of the ZAI
	short	fZ;		///< number of protons
	short	fA;		///< number of nucleons (A=0 means natural isotopes) 
	short	fI;		///< Isomeric state
	double	fMass;		///< Mass of a ZAI (from the BaseSummary.dat file
	ClassDef(ZAI,1);
};


#endif
