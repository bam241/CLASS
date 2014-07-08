#ifndef __EvolutionData_HXX__
#define __EvolutionData_HXX__

/*!
 \file
 \brief Header file for EvolutionData class.
 @version 2.0
 */

#include <string>
#include <map>

#include "IsotopicVector.hxx"
#include "CLASSObject.hxx"
#include "ZAI.hxx"

#include "TMatrix.h"

class TGraph;
class EvolutionData;
class CLASSLogger;

using namespace std;
typedef long long int cSecond;


EvolutionData operator*(EvolutionData const& evol, double F);
EvolutionData operator*(double F, EvolutionData const& evol);
EvolutionData operator/(EvolutionData const& evol, double F);
EvolutionData Sum(EvolutionData const& evol1, EvolutionData const& evol2);
EvolutionData Multiply(EvolutionData const& evol, double F);
EvolutionData Multiply(double F, EvolutionData const& evol);

double 	Distance(IsotopicVector IV1, EvolutionData Evd1 );
double 	Distance(EvolutionData Evd1, IsotopicVector IV1 );

//-----------------------------------------------------------------------------//
/*!
 Define a EvolutionData.
 The aim of these class is describe the evolution of a single evoluting system in CLASS.
 The system can either be a fuel evolution trough irradiation or a nuclei which produce, trough his decay, a large nuclei tree.
 
 The nuclei tree resulting of the evolution are stored in a map of ZAI and TGraph, each TGraph correspond to the evolution of the quantity of the associeted ZAI.

 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class EvolutionData : public CLASSObject
{
	
public :

//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	EvolutionData(); 	///< Normal DB Constructor.


	//{
	/// CLASSLogger Constructor.
	/*!
	 Use create an empty EvolutionData loading a CLASSLogger
	 \param CLASSLogger CLASSLogger used for the log...
	 */
	EvolutionData(CLASSLogger* Log); 	///< Make a new Evolutive Product evolution
	//}

	//{
	/// Special Constructor.
	/*!
	 Make a new EvolutionData
	 \param Log CLASSLogger used for the log...
	 \param DB_file path to the DataBase file
	 \param oldread true if the oldmethod should be use to read the DatBase File
	 \param zai set the ZAI if you want to add a stable nuclei.
	 */
	EvolutionData(CLASSLogger* Log, string DB_file, bool oldread = true, ZAI zai = ZAI(0,0,0) );
	//}


	//{
	/// Normal Destructor.
	/*!
	 Only remove the map without deleting the pointer to TGraph...
	 One need to call the DeleteEvolutionData() method to fully delete the EvolutionData, and then avoiding memory leak...
	 */
	~EvolutionData();
	//}

	//{
	/// Delete the EvolutionData.
	/*!
	 Use to fully delete the EvolutionData and all associeted TGraph.
	 In some case needed to be called to avoid memory leak.
	 */
	void DeleteEvolutionData();
	//}
	
	//@}


	

//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{
	void	SetHeavyMetalMass(double Mass)			{fHeavyMetalMass = Mass;}	//!< Set the HeavyMetal Mass 

	void 	SetReactorType(string reactortype)	{ fReactorType = reactortype; }		///< Set the reactor Type (string)
	void	SetFuelType(string fueltype)		{ fFuelType = fueltype; }		///< Set the fuel Type (string)
	void 	SetPower(double power)			{ fPower = power; }			///< Set the power of the EvolutionData (double)
	void	SetFlux(TGraph* flux )			{ fFlux = flux; }			///< Set the neutron flux of the EvolutionData (double)
	void	SetCycleTime(cSecond cycletime)		{ fCycleTime = cycletime; }		///< Set cycletime of the EvolutionData (double)


	void	SetEvolutionData(map<ZAI, TGraph*> maptoinsert)	{ fEvolutionData = maptoinsert;}///< Set EvolutionData map
	void	SetFissionXS(map<ZAI, TGraph*> maptoinsert)	{ fFissionXS = maptoinsert;}	///< Set fission cross section map
	void	SetCaptureXS(map<ZAI, TGraph*> maptoinsert)	{ fCaptureXS = maptoinsert;}	///< Set capture cross section map
	void	Setn2nXS(map<ZAI, TGraph*> maptoinsert)		{ fn2nXS = maptoinsert;}	///< Set (n,2n) cross section map

	//@}

	


//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{

#ifndef __CINT__
	map<ZAI ,TGraph* >	GetEvolutionData()	const { return fEvolutionData; }	//!< return the EvolutionData map
	map<ZAI ,TGraph* >	GetFissionXS()		const { return fFissionXS; }		//!< return the fission cross section map
	map<ZAI ,TGraph* >	GetCaptureXS()		const { return fCaptureXS; }		//!< return the capture cross section map
	map<ZAI ,TGraph* >	Getn2nXS()		const { return fn2nXS; }		//!< return the (n,2n) cross section map
	TGraph*			GetKeff()		const { return fKeff; }			//!< return the evolution of the keff (TGraph*)
	TGraph*			GetFlux()		const { return fFlux; }			//!< return the evolution of the neutron flux (TGraph*)
#endif

	double	GetFinalTime()		const { return fFinalTime; }			//!< return the final time - last point (double)
	double	GetCycleTime()		const { return fCycleTime; }			//!< return the cycletime (double)
	double	GetPower()		const { return fPower; }			//!< return the power (double)
	string	GetDB_file()		const { return fDB_file; }			//!< return the name of the Database file (string)
	string	GetReactorType()	const { return fReactorType; }			//!< return the type of reactor (string)
	TGraph*	GetEvolutionTGraph(const ZAI& zai);					//!< return the evolution of the ZAI quantity (TGraph*)

	IsotopicVector	GetIsotopicVectorAt(double t); 		///< Return the Product IsotopicVector at t time

	double	GetHeavyMetalMass()	const	{ return fHeavyMetalMass; }	//!< Return the HeavyMetal Mass in the Core at the begining of the cycle


	//{
	/// Return the XS for a reactionId on zai at t time
	/*!
	 // This method cross section of a reaction for a ZAI at a time
	 // \param t time
	 // \param ZAI ZAI for which the cross section if asked
	 // \param ReactionId ID of the reaction asked

	 // The different reaction ID are :
		\li 1 fission,
		\li 2 capture,
		\li 3 (n,2n).
	 */
	double	GetXSForAt(double t, ZAI zai, int ReactionId); 		
	//}

	//@}




//********* Insertion Method *********//

	//@}
	/*!
	 \name Insertion  Method
	 */
	//@{

	bool NucleiInsert(pair<ZAI, TGraph*> zaitoinsert);		//!< Add a nuclei evolution to the evolution map
	bool FissionXSInsert(pair<ZAI, TGraph*> zaitoinsert);		//!< Add a nuclei to the fission cross section map
	bool CaptureXSInsert(pair<ZAI, TGraph*> zaitoinsert);		//!< Add a nuclei to the capture cross section map
	bool n2nXSInsert(pair<ZAI, TGraph*> zaitoinsert);		//!< Add a nuclei to the (n,2n) cross section map
	

	//@}
	


protected :
	
	string	fDB_file;			///!< path to the DataBase file
	
	
#ifndef __CINT__
	map<ZAI ,TGraph* >	fEvolutionData;	//!< evolution map
	map<ZAI ,TGraph* >	fFissionXS;	//!< fission cross section map
	map<ZAI ,TGraph* >	fCaptureXS;	//!< capture cross section map
	map<ZAI ,TGraph* >	fn2nXS;		//!< (n,2n) cross section map
	TGraph*	fKeff;				//!< Keff evolution
	TGraph*	fFlux;				//!< Flux evolution
#endif
	
	cSecond	fFinalTime;			///< time of the last point
	bool	fIsCrossSection;		///< true if some cross section are present in the database
	
	
	
	string	fReactorType;			///< Type of reactor
	string	fFuelType;			///< Type of fuel
	double	fPower;				///< Power in W
	double	fCycleTime;			///< Cycle time of the DataBase
	double	fHeavyMetalMass;		///< Cycle time of the DataBase

    
	void	OldReadDB(string DBfile);				//!< Read old format database
	void	ReadDB(string DBfile, bool oldread = false);		//!< Main function to read database
	void	ReadKeff(string line, double* time, int NTimeStep);	//!< Read the Keff in the database
	void	ReadFlux(string line, double* time, int NTimeStep);	//!< Read the Flux in the database
	void	ReadInv(string line, double* time, int NTimeStep);	//!< Read the Inventory evolution in the database
	void	ReadXSFis(string line, double* time, int NTimeStep);	//!< Read the fission cross section evolution in the database
	void	ReadXSCap(string line, double* time, int NTimeStep);	//!< Read the capture cross evolution in the database
	void	ReadXSn2n(string line, double* time, int NTimeStep);	//!< Read the (n,2n) cross evolution in the database
	void	ReadInfo();						//!< Read the info file of the database

	
	double	Interpolate(double t, TGraph& EvolutionGraph);		///< Interpolating the value of EvolutionGraph at the t time

	void	AddAsStable(ZAI zai);					///< USe when adding a EvolutionData of a stable nuclei (for "non" decay)

	ClassDef(EvolutionData,0);
};




#endif
