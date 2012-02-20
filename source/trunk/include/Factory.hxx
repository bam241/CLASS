#ifndef __FACTORY_HXX__
#define __FACTORY_HXX__

/*!
 \file 
 \brief Header file for Factory class. 
*/
#include "CLASSHeaders.hxx"
using namespace std;
/*


class Factory
{
public :
	//! Normal Constructor.
 	Factory();
 	Factory(IsotopicVector isotopicvector, EvolutiveProduct* evolutivevie = NULL) 
 	{
 		fFactoryIsotopicVector = isotopicvector; 
 		fEvolutionDataBase = evolutivevie; 
 		
 	}

 	
 	//! Normal Destructor.
 	~Factory();
	
	IsotopicVector GetFactoryIsotopicvector() {return fFactoryIsotopicVector;} //<! Return the IV contain in the factory
	void SetIsotopicVector(IsotopicVector isotopicvector) {fFactoryIsotopicVector = isotopicvector;}
	
	virtual IsotopicVector* GetProductAt(double t); 	//!< Get IsotopicVector composition at the t time
	void Evolve(double t);
	
protected :
	
	IsotopicVector		fFactoryIsotopicVector;
	EvolutiveProduct*	fEvolutionDataBase;
 
};
*/
#endif
