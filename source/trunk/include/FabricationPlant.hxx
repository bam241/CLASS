#ifndef __FabricationPlant_HXX__
#define __FabricationPlant_HXX__

/*!
 \file
 \brief Header file for FabricationPlant class.
 @version 2.0
 */



#include <vector>
#include <map>

#include "CLASSConstante.hxx"
#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"
#include "Scenario.hxx"
#include "Storage.hxx"
#include "Reactor.hxx"
#include "CLASSLogger.hxx"
#include "ZAI.hxx"

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
//! CLASS object to build the fresh fuel (do chemical separation) & store it until core loading

/*!
 Define a FabricationPLant.
 The aim of these class is to manage the manufacturing of reprocessed fuel.
 It includes the fabrication of the fuel from a stock of material, using the appropriate
 algrorithm, and the storage of the fresh fuel until reactor loading.
 The parameters used for the fuel fabrication are recover from a PhysicsModels.
 The PhysicsModels MUST include an EquivalenceModel to build the fuel.
Some EquivalenceModel are available in the CLASS package, but an user can make his own.

Once the fuel is built, the FabricationPlant store the corresponding EvolutionData 
 generated using the PhysicsModels.
 @see PhysicsModels.hxx
 @see EquivalenceModel.hxx
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________

class DecayDataBank;
class FuelDataBank;


class FabricationPlant : public CLASSFacility
{

public :

//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	FabricationPlant();	///< Normal constructor


	//{
	/// Special Constructor.
	/*!
	 Make a new FabricationPlant
	 \param log : used for the log.
	 \param fabricationtime duration of the fabrication process (default : 2 years) in [s].
	 */
	FabricationPlant(CLASSLogger* log, double fabricationtime = cYear*2);
	//}

	~FabricationPlant(); 	///< Normal Destructor.

	//@}



	
//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetDecayDataBank(DecayDataBank* decayDB) {fDecayDataBase = decayDB;}	//! Set the Decay DataBank

	void SetFiFo(bool bval = true)	{ fFiFo = bval;}				//!< Set the chronological priority (true for chronological, false instead)
	
	void SetSubstitutionFuel(EvolutionData fuel);					//!< To use a substitution fuel if the fabrication fail (not enough material in stock)
	
	void AddReactor(int reactorid, double creationtime)
			{ fReactorNextStep.insert( pair<int,cSecond> (reactorid, (cSecond)creationtime-GetCycleTime() ) ); }	//!< Add a new reactor to be filled with the fresh fuel build by the FabricationPlant

#ifndef __CINT__
	void SetReUsableStorage(Storage* store) { fReUsable = store; fIsReusable = true;} //!< Set the Storage where all the separated matetial not used in the fabrication process will be sent. (if not present it goes to WASTE)
#endif

	using CLASSFacility::SetName;

	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
#ifndef __CINT__
	vector<Storage*>	GetFissileStorage()		{ return fFissileStorage; }		//!< Return the Pointer to the fissile Storage
	vector<Storage*>	GetFertileStorage()		{ return fFertileStorage; }		//!< Return the Pointer to the fertile Storage

	EvolutionData GetReactorEvolutionDB(int ReactorId);			//!< Return the EvolutionData of Reactor ReactorId
#endif
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at time t

	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}	//!< Return the list of the futur fuel IV


	//@}



#ifndef __CINT__
	void AddFissileStorage(Storage* stock) { fFissileStorage.push_back(stock); } //!< Add a new Storage to the list of fissile material provider.
	void AddFertileStorage(Storage* stock) { fFertileStorage.push_back(stock); } //!< Add a new Storage to the list of fertile material provider.
#endif

//********* Fabrication & Evolution Method *********//

	/*!
	 \name Fabrication & Evolution Method
	 */
	//@{

	void SetSeparartionEfficiencyIV(ZAI zai, double factor);	//!< Set the extraction efficiency of zai to factor (0<=factor<=1)
	void Evolution(cSecond t);					//!< Perform the FabricationPlant evolution
	
	void DumpStock(vector<double> lambdaArray);			//!< Update the Stock status after building process

	void TakeReactorFuel(int ReactorId) ;				//!< Remove the fuel of reactor ReactorId from stock
	void UpdateInsideIV();

	IsotopicVector BuildFuelFromEqModel(vector<double> LambdaArray); //!<Build the fresh fuel for the reactor according the results of the EquivalenceModel (@see  EquivalenceModel)
	void BuildFissileArray();					//!< virtualy extract fissile nuclei from Storage according EquivalenceModel fFissileList and make it virtually decay FabricationTime
	void BuildFertileArray();					//!< virtualy extract fertile nuclei from Storage according EquivalenceModel fFertileList and make it virtually decay FabricationTime

#ifndef __CINT__
	void BuildFuelForReactor(int ReactorId, cSecond t);			//!< Build a fuel for the reactor ReactorId
#endif

	void SortArray(int i); //!< Sort IsotopicVector array according priority preferences (e.g first in first out)

	//@}




protected :



//********* Internal Parameter *********//
	IsotopicVector	 fSeparationLostFraction;	///< The lost fraction table during separation (1- efficiency)
	map<int, cSecond >	fReactorNextStep;	///< Next time step to build a new fuel

#ifndef __CINT__
	map< int,EvolutionData >	fReactorFuturDB; ///< List of the futur EvolutionData use in the reactor
#endif
	map< int,IsotopicVector >	fReactorFuturIV; ///< List of the futur fuel IsotopicVector used in the reactor




	bool	fFiFo;					//!< First In First Out flag

	bool	fSubstitutionFuel;			//!< true if a substitution fuel as been set

	void	FabricationPlantEvolution(cSecond t);	//!< Deal the FabricationPlant evolution
	void 	ResetArrays();				//!< empty the fFertileArray and fFissileArray


#ifndef __CINT__
	
	vector<Storage*>	fFissileStorage;	//!< Pointer to the Storage used to get the fissile part of the fuel
	vector<IsotopicVector>  fFissileArray;		//!< The vector of isotopicVector use as fissile material
	vector<cSecond>		fFissileArrayTime;	//!< Time when a IsotopicVector arrives in its storage
	vector< pair<int,int> > fFissileArrayAdress;
	IsotopicVector		fFissileList;		//!< The list of fissile ZAI to consider

	vector<Storage*>	fFertileStorage;	//!< Pointer to the Storage used to get the fertile part of the fuel
	vector<IsotopicVector>  fFertileArray;		//!< The vector of isotopicVector used as fissile material
	vector<cSecond>		fFertileArrayTime;	//!< Time when a IsotopicVector arrives in its storage

	vector< pair<int,int> > fFertileArrayAdress;
	IsotopicVector		fFertileList;		//!< The List of fertile ZAI to consider

	Storage*	fReUsable;			//!< Pointer to the Storage using for storing unused material
	bool		fIsReusable;

	EvolutionData	fSubstitutionEvolutionData;	//!< EvolutionData of the subtitution fuel

	DecayDataBank*	fDecayDataBase;			//!< Pointer to the DecayDataBank


	//{
	/// Separation Method
	/*!
	 Make the Separation
	 \li IV[0] -> To Keep
	 \li IV[1] -> To Waste
	 */
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector, IsotopicVector ExtractedList);
	//}

#endif

	
	ClassDef(FabricationPlant,3);

};

#endif
