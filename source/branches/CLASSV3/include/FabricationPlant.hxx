#ifndef __FabricationPlant_HXX__
#define __FabricationPlant_HXX__

/*!
 \file
 \brief Header file for FabricationPlant class.
 @version 2.0
 */



#include <vector>
#include <map>

#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"
#include "CLASS.hxx"
#include "Storage.hxx"
#include "Reactor.hxx"
#include "LogFile.hxx"
#include "ZAI.hxx"

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
/*!
 Define a FabricationPLant.
 The aim of these class is describe the deal all the reprocessed fuel.
 It includes the fabrication of the fuel from a stock of used fuel, using the aproprieted algrorythm, and the storage of this fuel before putting it into a reactor.
 The parameter used for the fuel fabrication are recover from the DataBank.
 The Databank MUST include an equivalence model to build the fuel. This model is not necessary provided, each user need to put his own. By default a equivalence model is provided for PWR MOX fuel.

 The FabricationPlant once the fuel is builded, also store the corresponding EvolutionData generated using the DataBank.

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
	 Make a new FabricationPlant evolution
	 \param LogFile LogFile used for the log...
	 \param storage storage used to build the reprocessed fuel
	 \param reusable storage used to store all separated material not used in the fabrication process
	 \param fabricationtime duration of the fabrication process (2 years by default).
	 */
	FabricationPlant(LogFile* log, double fabricationtime = 365.25*24*3600*2);
	//}

	~FabricationPlant(); 	///< Normal Destructor.

	//@}



	
//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetDecayDataBank(DecayDataBank* decayDB) {fDecayDataBase = decayDB;}	//! Set the Decay DataBank

	void	SetUpdateReferenceDBatEachStep(bool val){ fUpdateReferenceDBatEachStep = val;}	//!< Set fUpdateReferenceDBatEachStep variable

	void	SetFiFo(bool bval = true)	{ fFiFo = bval;}	//!< Set the chronological priority (true for chronological, false unstead)
	
	void	SetSubstitutionFuel(EvolutionData fuel);				//!< To use a subtition fuel if the fabrication fail (not enough material in stock)
	
	void	AddReactor(int reactorid, double creationtime)
			{ fReactorNextStep.insert( pair<int,cSecond> (reactorid, (cSecond)creationtime-GetCycleTime() ) ); }	//!< Add a new reactor
	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
	vector<Storage*>	GetFissileStorage()		{ return fFissileStorage; }		//!< Return the Pointer to the Storage
	vector<Storage*>	GetFertileStorage()		{ return fFertileStorage; }		//!< Return the Pointer to the Storage

	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}	//!< Return the List of the Futur Fuel IV

	EvolutionData GetReactorEvolutionDB(int ReactorId);			//!< Return the EvolutionData of Reactor ReactorId
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time

	//@}



	void AddFissileStorage(Storage* stock) { fFissileStorage.push_back(stock); } //!< Add a new Storage to the list of Fissile material provider...
	void AddFertileStorage(Storage* stock) { fFertileStorage.push_back(stock); } //!< Add a new Storage to the list of Fertile material provider...

//********* Fabrication & Evolution Method *********//

	/*!
	 \name Fabrication & Evolution Method
	 */
	//@{

	void SetSeparartionEfficiencyIV(ZAI zai, double factor);	///< Add Valorisable Element
	void Evolution(cSecond t);					//!< Perform the Evolution
	
	void BuildFuelForReactor(int ReactorId);			//!< Build a Fuel for the reactor ReactorId
	void DumpStock(vector<double> lambdaArray);			//!< Update the Stock status after building process

	void TakeReactorFuel(int ReactorId) ;				//!< Remove the Fuel of reactor ReactorId


	void BuildFissileArray();
	void BuildFertileArray();
	IsotopicVector BuildFuelFromEqModel(vector<double> LambdaArray);

	void SortArray(int i);

	//@}




protected :
	bool		fUpdateReferenceDBatEachStep;	///< Set to true if the Reference Evolution Product must be updated at each calculation step (in the DataBank calculation)

//********* Internal Parameter *********//
	IsotopicVector	 fSeparationLostFraction;	///< The speration efficiency Table
	map<int, cSecond >	fReactorNextStep;	///< Next Time Step to Build a New Fuel

	map< int,EvolutionData >	fReactorFuturDB; ///< List of the Futur EvolutionData use in the reactor
	map< int,IsotopicVector >	fReactorFuturIV; ///< List of the Futur Fuel Isotopic Vector used in the reactor



	vector<Storage*>	fFissileStorage;		//!< Pointer to the Storage to recycle used to get the fissile part of the fuel
	vector<IsotopicVector>  fFissileArray;
	vector<cSecond>		fFissileArrayTime;
	vector< pair<int,int> > fFissileArrayAdress;
	IsotopicVector		fFissileList;

	vector<Storage*>	fFertileStorage;		//!< Pointer to the Storage to recycle used to get the fertile part of the fuel
	vector<IsotopicVector>  fFertileArray;
	vector<cSecond>		fFertileArrayTime;
	vector< pair<int,int> > fFertileArrayAdress;
	IsotopicVector		fFertileList;

	Storage*		fReUsable;			//!< Pointer to the Storage using for recycling unused Product

	bool		fFiFo;	//!< Set the First In First Out

	bool		fSubstitutionFuel;		//!< true if a subtitution fuel as been set
	EvolutionData	fSubstitutionEvolutionData;	//!< EvolutionData of the subtitution fuel
	
	DecayDataBank*	fDecayDataBase;			//!< Pointer to the Decay DataBase

//********* Private Method *********//
	void	FabricationPlantEvolution(cSecond t);	//!< Deal the FabricationPlant Evolution

	//{
	/// Separation Method
	/*!
	 Make the Separation
		\li IV[0] -> To Keep
		\li IV[1] -> To Waste
	 */
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector, IsotopicVector ExtractedList);
	//}

	
	ClassDef(FabricationPlant,3);

};

#endif
