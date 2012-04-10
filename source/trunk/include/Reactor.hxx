#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file 
 \brief Header file for Reactor class. 
*/

#include "IsotopicVector.hxx"
#include "TreatmentFactory.hxx"

using namespace std;


class CLASS;
class TreatmentFactory;
class EvolutiveProduct;

class Reactor : public TObject
{
public :
	//! Normal Constructor.
 	Reactor();
 
 	Reactor(EvolutiveProduct* evolutivedb , 
 		TreatmentFactory* TreatmentFactory ,
 		long int creationtime = 0, long int lifetime = (long int)(3600*24*365.4)*50 );

 	
 	//! Normal Destructor.
 	~Reactor();
	

//********* Get Method *********//

	IsotopicVector 		GetIVReactor()		{return fIVReactor;} 	//<! Return the IV contain in the Reactor
	IsotopicVector		GetIVBeginCycle()	{return fIVBeginCycle;}	//<! Return the Starting Cycle IV (Note : IVBegin != IVIn, only if using charging plan)
	IsotopicVector		GetIVOutCycle()		{return fIVOutCycle;}	//<! Return the Out Cycle IV 
	IsotopicVector		GetIVInCycle()		{return fIVInCycle;}	//<! Return the In Cycle IV (Note : IVIn != IVBegin, only if using charging plan)
										
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
	void Dump();
	
//********* Other Method *********//
	void Write(string Rbasename);
	
	
protected :
	long int			fInternalTime;		//!< Internal Clock
	long int 			fInCycleTime;		//!< Time spend since the beginning of the last Cycle
	
	
//********* Internal Parameter *********//
	long int			fCreationTime;		//!< CLASS Universal Time of Creation
	long int			fLifeTime;		//!< LifeTime Of the Reactor
	long int 			fCycleTime;		//!< Cycle Time

	IsotopicVector		fIVReactor;			//!< Fuel evoluated IV in the reactor
	IsotopicVector		fIVBeginCycle;			//!< Fuel IV at the Beginning of a Cycle
	IsotopicVector		fIVInCycle;			//!< IVBegin add at the Beginning of the Cycle
	IsotopicVector		fIVOutCycle;			//!< IV wich get out at the End of a Cycle
	
	
	EvolutiveProduct*	fEvolutionDB;			//!< Pointer to the Evolution DataBase
	CLASS*			fParc;				//!< Pointer to the main Parc
	TreatmentFactory*	fAssociedTreatmentFactory;	//!< Pointer to the TF which collect the spend fuel
 
 	bool 			fIsStarted;
 	bool 			fShutDown;
 	bool			fEndOfCycle;
 
 
 	ClassDef(Reactor,1);
 };


#endif
