#ifndef _FabricationPlant_
#define _FabricationPlant_

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

//_________________________________________________________________________________________________________
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
 @author BLG
 @author FaC
 
 @version 2.0
 */
//________________________________________________________________________

	//! Define the storage management for fuel fresh construction.
	/*!
	// Posible  Storage Management are
		using : YourFabPlant->SetStorageManagement(key);
		kpFiFo : First In First Out (i.e the older storage first)
		kpLiFo : Last In First Out  (i.e the youger storage first)
		kpMix : IVs are sorted that way : First: The younger , 
		Second: The older, Third: The second younger ,4th : the second older ....
		kpRand : IVs order in storage is random

	*/
	enum StorageManagement{ kpFiFo, kpLiFo ,kpMix, kpRand};

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

    void SetFiFo(bool bval = true)	{ if(bval) fStorageManagement=kpFiFo; else fStorageManagement=kpLiFo ;}	//!< Set the chronological priority (true for chronological, false instead).Equivalent to SetStorageManagement(kpFiFo) or SetStorageManagement(kpLiFo)
    void SetStorageManagement(StorageManagement SM){fStorageManagement = SM ;} //!<  The storage management : either kpFiFo, kpLiFo , kpMix or kpRand

	void SetSubstitutionMaterialFromIV(string keyword, IsotopicVector SubstitutionIV) 						//!< If the construction fails : it creates a substitution material according to the IV defined by the user
		{fSubstitutionMaterialFromIV[keyword] = true; fSubstitutionIV[keyword]= SubstitutionIV;} 	
	
	void SetSubstitutionFuel(EvolutionData fuel);											//!< To use a substitution fuel if the fabrication fail (not enough material in stock)

    void SetSeparationManagement(bool bval = true)	{ fIsSeparationManagement = bval;}				//!< Set the separation managmeent for the fabrication plant

	void AddReactor(int reactorid, double creationtime)
		{ fReactorNextStep.insert( pair<int,cSecond> (reactorid, (cSecond)creationtime-GetCycleTime() ) ); }		//!< Add a new reactor to be filled with the fresh fuel build by the FabricationPlant

#ifndef __CINT__
	void SetReUsableStorage(Storage* store) { fReUsable = store; fIsReusable = true;} 					//!< Set the Storage where all the separated matetial not used in the fabrication process will be sent. (if not present it goes to WASTE)
#endif

	using CLASSFacility::SetName;

	//@}

#ifndef __CINT__
		
	void AddStorage(string keyword, Storage* Stock, double MassFractionMin = 0, double MassFractionMax = 1., int Priority = 0) ;	//!< Fill the storage vector for a list
	void AddInfiniteStorage(string keyword, double MassFractionMin = 0, double MassFractionMax = 1., int Priority = 0);		//!< Creates an infinite stock of this material according to the list defined in the EqM
	void AddFuelBuffer(string keyword);													//!< Tell the buffer for this fuel. Creates an infinite stock of this material according to the list defined in the EqM
	void AddFuelBuffer(string keyword, Storage* Stock);											//!< Tell the buffer for this fuel taken from the storage


#endif
	
//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
#ifndef __CINT__
	map < string , vector <Storage*> > GetAllStorage() {return fStorage;}			//!< Return the map containing all the storage vectors (useful in CLASS Reactor to check list consistency)
	
	vector<Storage*> GetStorage(string keyword) { return fStorage[keyword]; }		//!< Return the Pointer to Storage associated to a StreamList 

	EvolutionData GetReactorEvolutionDB(int ReactorId);						//!< Return the EvolutionData of Reactor ReactorId
#endif
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);				//!< Get IsotopicVector Decay at time t

	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}				//!< Return the list of the futur fuel IV


	//@}

//********* Fabrication & Evolution Method *********//

	/*!
	 \name Fabrication & Evolution Method
	 */
	//@{
	void SetSeparationEfficiency(IsotopicVector IV,  cSecond TimeOfSep);	//!< Set the extraction efficiency of IsotopicVector IV.This separation efficiency is effectove at time TimeOfSep
	IsotopicVector GetSeparationEfficiencyAt(cSecond time);


	void Evolution(cSecond t);									//!< Perform the FabricationPlant evolution
	void DumpStock(map <string , vector<double> > LambdaArray);				//!< Update the Stock status after building process
	void TakeReactorFuel(int ReactorId) ;								//!< Remove the fuel of reactor ReactorId from stock
	void UpdateInsideIV();										

	IsotopicVector BuildFuelFromEqModel(map <string , vector<double> > LambdaArray); 	//!<Build the fresh fuel for the reactor according the results of the EquivalenceModel (@see  EquivalenceModel)
	void BuildArray(int ReactorId, cSecond ReactorLoadingTime);									//!< virtualy extract fissile nuclei from Storage according EquivalenceModel fStreamList and make it virtually decay FabricationTime

#ifndef __CINT__
	void BuildFuelForReactor(int ReactorId, cSecond t);			//!< Build a fuel for the reactor ReactorId
#endif

	void SortArray(); //!< Sort IsotopicVector array according priority preferences (given by key in YourFabPlant->SetStorageManagement(key);)

	void SortFiFo(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray);
	void SortLiFo(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray);
	void SortMix(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray);
	void SortRandom(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray);

	//@}




protected:



//********* Internal Parameter *********//
	IsotopicVector	 fSeparationLostFraction;	///< The lost fraction table during separation (1- efficiency)

	map<cSecond, IsotopicVector> fSeparationStrategy;

	map<int, cSecond >	fReactorNextStep;	///< Next time step to build a new fuel

#ifndef __CINT__
	map< int,EvolutionData >	fReactorFuturDB; ///< List of the futur EvolutionData use in the reactor
#endif
	map< int,IsotopicVector >	fReactorFuturIV; ///< List of the futur fuel IsotopicVector used in the reactor



    StorageManagement fStorageManagement;	//!< The storage management : either kpFiFo, kpLiFo , kpMix or kpRand

    bool	fIsSeparationManagement; //!< Separation managment

    bool	fSubstitutionFuel;										//!< True if a substitution fuel as been set

	void	FabricationPlantEvolution(cSecond t);							//!< Deal the FabricationPlant evolution
	void 	ResetArrays();										//!< Empty Arrays


#ifndef __CINT__

	map < string , IsotopicVector>	fStreamList;						//!< Map that contains lists of stream according to the EqModel with corresponding isotopes list
	map < string , double>		fStreamListFPMassFractionMax;			//!< Map that contains lists of stream according to the EqModel with mass maximum fraction
	map < string , double>		fStreamListFPMassFractionMin;			//!< Map that contains lists of stream according to the EqModel with mass minimum fraction
	map < string , int>			fStreamListFPPriority;					//!< Map that contains lists of stream according to the EqModel with priority (1 = first, 2 = second, etc...)
	map < string , bool>			fStreamListFPIsBuffer;					//!< Map that contains lists of stream according to the EqModel saying if fuel buffer

	map < string , vector <Storage*> >		fStorage;					//!< Pointer to the Storages defined for each list
	map < string , vector <IsotopicVector> >  	fStreamArray;					//!< The vector of isotopicVector of each material and each stock
	map < string , vector <cSecond> >	  	fStreamArrayTime;				//!< Time when a IsotopicVector arrives in its storage
	map < string , vector < pair<int,int> > >  	fStreamArrayAdress;
	map < string , IsotopicVector>		fSubstitutionIV;					//!< contains the susbstitution IV defined by the user 

	map < string , bool > fSubstitutionMaterialFromIV;						//!< True = a substitution IV is set for this material in case of failure in fuel building 
	map < string , bool > fInfiniteMaterialFromList;						//!< True = an infinite stock of this material is created according to the list defined in the EqM

	map < string , bool > fErrorOnLambda;							//!< True = lambdas haven't been well calculated for this material (not enough material in stock....)

	EvolutionData	fSubstitutionEvolutionData;							//!< EvolutionData of the subtitution fuel

	Storage*	fReUsable;									//!< Pointer to the Storage used to storing unused material
	bool		fIsReusable;									//!< Sets a storage used to storing unused material
	bool		fFuelCanBeBuilt;								//!< Default fuel fabrication process has failed
	DecayDataBank*	fDecayDataBase;							//!< Pointer to the DecayDataBank


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
