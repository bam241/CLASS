#ifndef __Reactor_HXX__
#define __Reactor_HXX__

/*!
 \file 
 \brief Header file for Reactor class. 
*/

#include "IsotopicVector.hxx"
//#include "TreatmentFactory.hxx"
#include "LogFile.hxx"

using namespace std;


class CLASS;
class TreatmentFactory;
class EvolutiveProduct;
class LogFile;

class Reactor : public TObject
{
public :
	
 	Reactor();	///< Normal Constructor.
 
 	Reactor(EvolutiveProduct* evolutivedb , 
 		TreatmentFactory* TreatmentFactory ,
 		long int creationtime = 0, long int lifetime = (long int)(3600*24*365.4)*50 );	///< Advbanced Constructor.

 	
 	
 	~Reactor();	///< Normal Destructor.
	

//********* Get Method *********//

	IsotopicVector 		GetIVReactor()		{return fIVReactor;} 	//!< Return the IV contain in the Reactor
	IsotopicVector		GetIVBeginCycle()	{return fIVBeginCycle;}	//!< Return the Starting Cycle IV (Note : IVBegin != IVIn, only if using charging plan)
	IsotopicVector		GetIVOutCycle()		{return fIVOutCycle;}	//!< Return the Out Cycle IV 
	IsotopicVector		GetIVInCycle()		{return fIVInCycle;}	//!< Return the In Cycle IV (Note : IVIn != IVBegin, only if using charging plan)
	long int 		GetCycleTime()		{return fCycleTime;} 	//!< Return the cycle time of the Reactor
	long int 		GetCreationTime()	{return fCreationTime;}	//!< Return the creation time of the Reactor
	long int 		GetLifeTime()		{return fLifeTime;}	//!< Return the creation time of the Reactor

	EvolutiveProduct*	GetEvolutionDB()		{return fEvolutionDB;}			//!< Return the database of the Fuel composition evolution
	TreatmentFactory*	GetAssociedTreatmentFactory()	{return fAssociedTreatmentFactory;}	//!< Return the pointer to Associeted TF
	LogFile*		GetLog()			{return fLog;}				//!< Return the pointer to Log

//********* Set Method *********//
	void SetParc(CLASS* parc)				{fParc = parc;}				//!< Set the Pointer to the Parc
	void SetLog(LogFile* LOG)				{fLog = LOG;}				//!< Set the Pointer to the Log
	void SetIVReactor(IsotopicVector isotopicvector)	{fIVReactor = isotopicvector;}		//!< Set the IV inside the Reactor Core
	void SetIVBeginCycle(IsotopicVector isotopicvector)	{fIVBeginCycle = isotopicvector;}	//!< Set the IV at the Beginging of the Reactor Cycle
	void SetIVOutCycle(IsotopicVector isotopicvector)	{fIVOutCycle = isotopicvector;}		//!< Set the IV Going Out at the End of the Cycle
	void SetIVInCycle(IsotopicVector isotopicvector)	{fIVInCycle = isotopicvector;}		//!< Set the IV Coming In at the Beginning of the Cycle 
	void SetCycleTime(long int cycletime)			{fCycleTime = cycletime;}		//!< Set the Cycle time (Cycle of the loading Plan)
	void SetEvolutionDB(EvolutiveProduct* evolutionDB)	{fEvolutionDB = evolutionDB;}		//!< Set the Pointer to the DB Evolution of the Reactor
	
//********* Modification Method *********//
	void Evolution(long int t);									//!< Performe the Evolution until the Time t
	void Dump();											//!< Write Modification (IV In/Out, filling the TF...)
	
//********* Other Method *********//
	
	
protected :
	long int		fInternalTime;		///< Internal Clock
	long int 		fInCycleTime;		///< Time spend since the beginning of the last Cycle
 	bool 			fIsStarted;		///< True if Running, False Otherwise
 	bool 			fShutDown;		///< True if ShutDown
 	bool			fEndOfCycle;		///< True if Reaching the End of a Reactor Cycle
	
	
//********* Internal Parameter *********//
 	LogFile*		fLog;				//!< Pointer to the Log
	CLASS*			fParc;				//!< Pointer to the main Parc
	TreatmentFactory*	fAssociedTreatmentFactory;	//!< Pointer to the TF which collect the spend fuel
	EvolutiveProduct*	fEvolutionDB;			//!< Pointer to the Evolution DataBase

	
	long int		fCreationTime;	///< CLASS Universal Time of Creation
	long int		fLifeTime;	///< LifeTime Of the Reactor
	long int 		fCycleTime;	///< Cycle Time

	IsotopicVector		fIVReactor;	///< Fuel evoluated IV in the reactor
	IsotopicVector		fIVBeginCycle;	///< Fuel IV at the Beginning of a Cycle
	IsotopicVector		fIVInCycle;	///< IVBegin add at the Beginning of the Cycle
	IsotopicVector		fIVOutCycle;	///< IV wich get out at the End of a Cycle
	
 
 	ClassDef(Reactor,3);
 };


#endif
