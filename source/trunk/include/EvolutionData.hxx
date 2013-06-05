#ifndef __EvolutionData_HXX__
#define __EvolutionData_HXX__

/*!
 \file
 \brief Header file for EvolutionData classes. 
  The aim of this Class is to manage evolution of a system, such as a radioactiv nuclei or a reactor. It store the evolution of daughter nuclei proportion as the function of time of the 
 
 @author BaM
 @version 0.
 */

#include <string>
#include <map>
#include "IsotopicVector.hxx"

#include "ZAI.hxx"
#include "TObject.h"
#include "TMatrix.h"

class TGraph;
class EvolutionData;
class LogFile;

using namespace std;
typedef long long int cSecond;

///< A ZAIIDataBase defined a database which contain the evolution of faction of all product, subproduct (or sub(sub...sub)product) for a nucleus. 
/*!
 The aim of this class is to handle the evolution all Information of all Nucleus product, subproduct (or sub(sub...sub)product) as the fonction of the time.
 
 @author BaM
 @version 1.0
 */

EvolutionData operator*(EvolutionData const& evol, double F);
EvolutionData operator*(double F, EvolutionData const& evol);
EvolutionData operator/(EvolutionData const& evol, double F);



class EvolutionData : public TObject
{
	
public :

//********* Constructor/Destructor Method *********//
	///< Normal DB Constructor.
	
	EvolutionData();  
	EvolutionData(LogFile* Log); 	///< Make a new Evolutive Product evolution 
	EvolutionData(LogFile* Log, string DB_file, bool oldread = false, ZAI zai = ZAI(0,0,0) ); 	///< Make a new Evolutive Product evolution
	 
	///< Normal Destructor.
	~EvolutionData();

//********* Set Method *********//
	void 	SetReactorType(string reactortype)	{ fReactorType = reactortype; }
	void	SetFuelType(string fueltype)		{ fFuelType = fueltype; }
	void 	SetPower(double power)			{ fPower = power; }
	void 	SetHMMass(double HMMass)		{ fHMMass = HMMass; }
	void	SetFlux(TGraph* flux )			{ fFlux = flux; }
	
//********* Get Method *********//
	map<ZAI ,TGraph* >	GetEvolutionData()	const { return fEvolutionData; }	//!<
	map<ZAI ,TGraph* >	GetFissionXS()		const { return fFissionXS; }		//!<
	map<ZAI ,TGraph* >	GetCaptureXS()		const { return fCaptureXS; }		//!<
	map<ZAI ,TGraph* >	Getn2nXS()		const { return fn2nXS; }		//!<
	TGraph*			GetKeff()		const { return fKeff; }
	TGraph*			GetFlux()		const { return fFlux; }

	double	GetCycleTime()		const { return fCycleTime; }
	double	GetPower()		const { return fPower; }		//!<
	double	GetHMMass()		const { return fHMMass; }
	string	GetDB_file()		const { return fDB_file; }
	string	GetReactorType()	const { return fReactorType; }
	TGraph*	GetEvolutionTGraph(const ZAI& zai); 
								///< Return the A,Z product proportion evolution TGraph
	IsotopicVector	GetIsotopicVectorAt(double t); 		///< Return the Product IsotopicVector at t time
	
	double	GetGetXSForAt(double t, ZAI zai, int ReactionId); 		///< Return the XS for a reactionId on zai at t time
										///< ReactionId : 1 Fission,
										///<		2 Capture,
										///<		3 (n, 2n) ,




	bool NucleiInsert(pair<ZAI, TGraph*> zaitoinsert);			//!<
	bool FissionXSInsert(pair<ZAI, TGraph*> zaitoinsert);
	bool CaptureXSInsert(pair<ZAI, TGraph*> zaitoinsert);
	bool n2nXSInsert(pair<ZAI, TGraph*> zaitoinsert);
	
	
//********* Get Method *********//
	EvolutionData GenerateDBFor(IsotopicVector isotopicvector);	///< Build A DB from a close one


protected :
	
	string	fDB_file;
	map<ZAI ,TGraph* >	fEvolutionData;	//!< 
	map<ZAI ,TGraph* >	fFissionXS;	//!< 
	map<ZAI ,TGraph* >	fCaptureXS;	//!< 
	map<ZAI ,TGraph* >	fn2nXS;	//!< 
	TGraph*	fKeff;
	TGraph*	fFlux;
	
	cSecond	fFinalTime;
	bool	fIsCrossSection;
	
	
	
	string	fReactorType;
	string	fFuelType;
	double	fPower;
	double	fCycleTime;
	double 	fHMMass;
	
    
	void	OldReadDB(string DBfile);
	void	ReadDB(string DBfile, bool oldread = false);
	void	ReadKeff(string line, double* time);
	void	ReadFlux(string line, double* time);
	void	ReadInv(string line, double* time);
	void	ReadXSFis(string line, double* time);
	void	ReadXSCap(string line, double* time);
	void	ReadXSn2n(string line, double* time);
	void	ReadInfo();

	
	double	Interpolate(double t, TGraph& EvolutionGraph);
								///< Interpolating the value of EvolutionGraph at the t time
	void	AddAsStable(ZAI zai);
	LogFile*	fLog;

	ClassDef(EvolutionData,0);
};

#endif
