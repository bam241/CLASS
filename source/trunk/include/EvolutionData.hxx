#ifndef __EvolutionData_HXX__
#define __EvolutionData_HXX__

/*!
 \file
 \brief Header file for EvolutionData classes. 
  The aim of this Class is to manage evolution of a system, such as a radioactiv nuclei or a reactor. It store the evolution of daughter nuclei proportion as the function of time of the 
 

 @author BaM
 @version 2.0
 */

#include <string>
#include <map>

#include "IsotopicVector.hxx"
#include "CLSSObject.hxx"
#include "ZAI.hxx"

#include "TMatrix.h"

class TGraph;
class EvolutionData;
class LogFile;

using namespace std;
typedef long long int cSecond;


EvolutionData operator*(EvolutionData const& evol, double F);
EvolutionData operator*(double F, EvolutionData const& evol);
EvolutionData operator/(EvolutionData const& evol, double F);



class EvolutionData : public CLSSObject
{
	
public :

//********* Constructor/Destructor Method *********//
	///< Normal DB Constructor.
	
	EvolutionData();  
	EvolutionData(LogFile* Log); 	///< Make a new Evolutive Product evolution 
	EvolutionData(LogFile* Log, string DB_file, bool oldread = true, ZAI zai = ZAI(0,0,0) ); 	///< Make a new Evolutive Product evolution
	 
	///< Normal Destructor.
	~EvolutionData();

//********* Set Method *********//
	void 	SetReactorType(string reactortype)	{ fReactorType = reactortype; }
	void	SetFuelType(string fueltype)		{ fFuelType = fueltype; }
	void 	SetPower(double power)			{ fPower = power; }
	void	SetFlux(TGraph* flux )			{ fFlux = flux; }
	void	SetFissionXS(map<ZAI, TGraph*> maptoinsert)	{ fFissionXS = maptoinsert;}
	void	SetCaptureXS(map<ZAI, TGraph*> maptoinsert)	{ fCaptureXS = maptoinsert;}
	void	Setn2nXS(map<ZAI, TGraph*> maptoinsert)	{ fn2nXS = maptoinsert;}
	void	SetCycleTime(cSecond cycletime)		{ fCycleTime = cycletime; }

	
	
//********* Get Method *********//
#ifndef __CINT__
	map<ZAI ,TGraph* >	GetEvolutionData()	const { return fEvolutionData; }	//!<
	map<ZAI ,TGraph* >	GetFissionXS()		const { return fFissionXS; }		//!<
	map<ZAI ,TGraph* >	GetCaptureXS()		const { return fCaptureXS; }		//!<
	map<ZAI ,TGraph* >	Getn2nXS()		const { return fn2nXS; }		//!<
	TGraph*			GetKeff()		const { return fKeff; }
	TGraph*			GetFlux()		const { return fFlux; }
#endif

	double	GetFinalTime()		const { return fFinalTime; }
	double	GetCycleTime()		const { return fCycleTime; }
	double	GetPower()		const { return fPower; }		//!<
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


protected :
	
	string	fDB_file;
	
	
#ifndef __CINT__
	map<ZAI ,TGraph* >	fEvolutionData;	//!<
	map<ZAI ,TGraph* >	fFissionXS;	//!< 
	map<ZAI ,TGraph* >	fCaptureXS;	//!< 
	map<ZAI ,TGraph* >	fn2nXS;		//!< 
	TGraph*	fKeff;
	TGraph*	fFlux;
#endif
	
	cSecond	fFinalTime;
	bool	fIsCrossSection;
	
	
	
	string	fReactorType;
	string	fFuelType;
	double	fPower;
	double	fCycleTime;
	double	fNormFactor;

    
	void	OldReadDB(string DBfile);
	void	ReadDB(string DBfile, bool oldread = false);
	void	ReadKeff(string line, double* time, int NTimeStep);
	void	ReadFlux(string line, double* time, int NTimeStep);
	void	ReadInv(string line, double* time, int NTimeStep);
	void	ReadXSFis(string line, double* time, int NTimeStep);
	void	ReadXSCap(string line, double* time, int NTimeStep);
	void	ReadXSn2n(string line, double* time, int NTimeStep);
	void	ReadInfo();

	
	double	Interpolate(double t, TGraph& EvolutionGraph);
								///< Interpolating the value of EvolutionGraph at the t time
	void	AddAsStable(ZAI zai);

	ClassDef(EvolutionData,0);
};

#endif
