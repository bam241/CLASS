#ifndef _CLASS_HXX_
#define _CLASS_HXX_

/*!
 \file 
 \brief Header file for CLASS class. 
*/
#include "CLASSHeaders.hxx"
#include "Reactor.hxx"
#include "TreatmentFactory.hxx"

using namespace std;


class CLASS
{
public :
	//! Normal Constructor.
 	CLASS();
 	CLASS(long int abstime);
 	
 	//! Normal Destructor.
 	~CLASS();
	

//********* Get Method *********//
	long int			GetAbsoluteTime()	{return fAbsoluteTime;}		//!< Absolute Clock
 	map<long int, int>		GetTimeStep()		{return fTimeStep;}
 	vector<TreatmentFactory*>	GetTreatmentFactory()	{return fTreatmentFactory;}
 	vector<Reactor*>		GetReactor()		{return fReactor;}
	long int			GetPrintSet()		{return fPrintStep;}

//********* Set Method *********//
//	void	SetAbsoluteTime(double time);		//!< Absolute Clock
	void	SetTimeStep(int timestep) 				{fPrintStep = timestep;}

//********* Modification Method *********//
	IsotopicVector	BuildIsotopicVector(IsotopicVector isotopicvector );
	void	BuildTimeVector(long int t);
	void	Evolution(long int t);
 	void	AddTreatmentFactory(TreatmentFactory* treatmentfactory);
 	void	AddReactor(Reactor* reactor);
	
//********* Other Method *********//
	void	TreatmentEvolution();
	void	ReactorEvolution();
	void	RemoveReactor();
	void	Print();
	void	Write();

	
protected :
 	long int			fPrintStep;
	long int			fAbsoluteTime;		//!< Absolute Clock
 	map<long int, int>		fTimeStep;		//!< 1 to print, 0 otherwise (can et new property...)
 	vector<TreatmentFactory*>	fTreatmentFactory;
 	vector<int>			fTreatmentFactoryIndex;
 	int				fTreatmentFactoryLastIndex;
 	vector<Reactor*>		fReactor;
 	vector<int>			fReactorIndex;
	int				fReactorLastIndex;

};


#endif
