#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file 
 \brief Header file for Reactor class. 
*/
#include "CLASSHeaders.hxx"
#include "CLASS.hxx"
#include "TreatmentFactory.hxx"
using namespace std;


class Reactor
{
public :
	//! Normal Constructor.
 	Reactor();
 	Reactor(EvolutiveProduct* evolutivedb = NULL, 
 		TreatmentFactory* TreatmentFactory = NULL,
 		long int creationtime = 0, long int lifetime = (long int)(3600*24*365.4)*50 );

 	
 	//! Normal Destructor.
 	~Reactor();
	

//********* Get Method *********//

	IsotopicVector 		GetIVReactor()		{return fIVReactor;} 	//<! Return the IV contain in the Reactor
	IsotopicVector		GetIVBeginCycle()	{return fIVBeginCycle;}
	IsotopicVector		GetIVOutCycle()		{return fIVOutCycle;}
	IsotopicVector		GetIVInCycle()		{return fIVInCycle;}
	long int 		GetCycleTime()		{return fCycleTime;} 	//!< Return the cycle time of the Reactor
	long int 		GetCreationTime()	{return fCreationTime;}	//!< Return the creation time of the Reactor
	long int 		GetLifeTime()		{return fLifeTime;}	//!< Return the creation time of the Reactor
	EvolutiveProduct*	GetEvolutionDB()	{return fEvolutionDB;}	//!< Return the database of the Fuel composition evolution
	TreatmentFactory*	GetAssociedTreatmentFactory() {return fAssociedTreatmentFactory;}


//********* Set Method *********//
	void SetParc(CLASS* parc)				{fParc = parc;}
	void SetIVReactor(IsotopicVector isotopicvector)	{fIVReactor = isotopicvector;}
	void SetIVBeginCycle(IsotopicVector isotopicvector)	{fIVBeginCycle = isotopicvector;}
	void SetIVOutCycle(IsotopicVector isotopicvector)	{fIVOutCycle = isotopicvector;}
	void SetIVInCycle(IsotopicVector isotopicvector)	{fIVInCycle = isotopicvector;}
	void SetCycleTime(long int cycletime)			{fCycleTime = cycletime;} 
	void SetEvolutionDB(EvolutiveProduct* evolutionDB)	{fEvolutionDB = evolutionDB;}

//********* Modification Method *********//
//	IsotopicVector* GetProductAt(double t); 	//!< Get IsotopicVector composition at the t time
	void Evolution(long int t);
	
//********* Other Method *********//
	void Write(string Rbasename);
	
	
protected :
	long int			fInternalTime;		//!< Internal Clock
	long int 			fInCycleTime;		
	
	
//********* Internal Parameter *********//
	long int			fCreationTime;		//!< CLASS Universal Time of Creation
	long int			fLifeTime;		//!< LifeTime Of the Reactor
	long int 			fCycleTime;		//!< Cycle Time

	EvolutiveProduct*	fEvolutionDB;	//!< Pointer to the Evolution DataBase

	IsotopicVector		fIVReactor;
	IsotopicVector		fIVBeginCycle;
	IsotopicVector		fIVInCycle;
	IsotopicVector		fIVOutCycle;
	
	
	CLASS*			fParc;
	TreatmentFactory*	fAssociedTreatmentFactory;
 
 	bool 			IsStarted;
};


#endif
