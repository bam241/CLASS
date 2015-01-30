#ifndef _ZAI_
#define _ZAI_ 

/*!
 \file
 \brief Header file for ZAI classes. 
 */

#include <string>
#include "TObject.h"
#include <iostream>

using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a nuclei as 3 integer Z,A,I:
 \li its charge number : Z
 \li its mass number : A
 \li its Isomeric state : I (0 : fundamental, 1 : first isomeric state ,...)

 The aim of this class is to discribe each ZAI.

 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class ZAI : public TObject
{
public:


//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	ZAI();	///< Default constructor

	//{
 	///< Normal Constructor. 
 	/*!
 	 Default: No parent
 	 \param Z : number of protons
 	 \param A : number of nucleons (A=0 means natural isotopes) 
 	*/
	ZAI(int Z, int A, int I=0);
	//}


	~ZAI();	///< Normal Destructor.


//********* ZAI main attributes Method *********//

	/*!
	\name ZAI main attributes 
	*/
	//@{
	int  Z()	const	{ return fZ; } //!< returns the number of protons
	int  A()	const	{ return fA; } //!< returns the number of nucleons
	int  I()	const	{ return fI; } //!< returns the Isomeric State
	int  N()	const	{ return fA-fZ; } //!< returns the number of neutrons
		
	//@}



	ZAI operator=(ZAI IVa);		//!<
	bool operator <(const ZAI& zai)		const	{ return (fZ != zai.Z())?
							(fZ < zai.Z()) : ( (fA != zai.A())?
								 (fA < zai.A()) : (fI < zai.I()) ); }//!< Compartor operator

	bool operator !=(const ZAI& zai)	const	{ return ( fZ != zai.Z() ) || ( fA != zai.A() ) || ( fI != zai.I() ); }//!< Compartor operator
	bool operator ==(const ZAI& zai)	const	{ return ( fZ == zai.Z()  && fA == zai.A() &&  fI == zai.I()); }//!< Compartor operator
	void Print()				const	{ cout << fZ << " " << fA << " " << fI << endl;}//!< Print in standard output : Z A I


protected :
 	
	short	fZ;		///< number of protons
	short	fA;		///< number of nucleons (A=0 means natural isotopes) 
	short	fI;		///< Isomeric state

	ClassDef(ZAI,1);
};


#endif
