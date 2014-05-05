#ifndef __FabricationPlant_HXX__
#define __FabricationPlant_HXX__

/*!
 \file
 \brief Header file for FabricationPlant class.
 @version 2.0
 */



#include <vector>
#include <map>

#include "CLSSFacility.hxx"
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


class FabricationPlant : public CLSSFacility
{

public :

//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	FabricationPlant();	///< Normal constructor


	//{
	/// LogFile Constructor.
	/*!
	 Use create an empty FabricationPlant loading a LogFile
	 \param LogFile LogFile used for the log...
	 */
	FabricationPlant(LogFile* log);
	//}


	//{
	/// Special Constructor.
	/*!
	 Make a new FabricationPlant evolution
	 \param LogFile LogFile used for the log...
	 \param storage storage used to build the reprocessed fuel
	 \param reusable storage used to store all separated material not used in the fabrication process
	 \param fabricationtime duration of the fabrication process (2 years by default).
	 */
	FabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabricationtime = 365.25*24*3600*2);
	//}

	~FabricationPlant(); 	///< Normal Destructor.

	//@}



	
//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void	SetUpdateReferenceDBatEachStep(bool val){ fUpdateReferenceDBatEachStep = val;}	//!< Set fUpdateReferenceDBatEachStep variable
	void	SetStorage(Storage* storage)		{ fStorage = storage; }			//!< Set the Pointer to the Storage
	
	void	SetChronologicalTimePriority(bool bval)	{ fChronologicalTimePriority = bval;}	//!< Set the chronological priority (true for chronological, false unstead)
	
	void	SetSubstitutionFuel(EvolutionData fuel);				//!< To use a subtition fuel if the fabrication fail (not enough material in stock)
	
	void	AddReactor(int reactorid, double creationtime)
			{ fReactorNextStep.insert( pair<int,cSecond> (reactorid, (cSecond)creationtime-GetCycleTime() ) ); }	//!< Add a new reactor
	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
	Storage*	GetStorage()		{ return fStorage; }		//!< Return the Pointer to the Storage
	
	map<int, IsotopicVector >	GetReactorFuturIncome() const
						{ return fReactorFuturIV;}	//!< Return the List of the Futur Fuel IV

	IsotopicVector GetFullFabrication();					//!< Return the Sum of all Fuel waiting to be put in a reator

	EvolutionData GetReactorEvolutionDB(int ReactorId);			//!< Return the EvolutionData of Reactor ReactorId

	//@}





//********* Fabrication & Evolution Method *********//

	/*!
	 \name Fabrication & Evolution Method
	 */
	//@{

	void	AddValorisableIV(ZAI zai, double factor);						///< Add Valorisable Element
	void	Evolution(cSecond t);									//!< Perform the Evolution
	virtual void	BuildFuelForReactor(int ReactorId);							//!< Build a Fuel for the reactor ReactorId
	void	RecycleStock(double fraction);								//!< Take a franction of the current stock
	IsotopicVector	GetStockToRecycle();								//!< Get the next stock to recycle
	void 	DumpStock();										//!< Update the storage
	EvolutionData	BuildEvolutiveDB(int ReactorId, IsotopicVector isotopicvector);		//!< Build the Evolution Database for the Reactir ReactorId Fuel
	void	TakeReactorFuel(int ReactorId) ;								//!< Remove the Fuel of reactor ReactorId
	
	//@}




protected :
	bool		fUpdateReferenceDBatEachStep;	///< Set to true if the Reference Evolution Product must be updated at each calculation step (in the DataBank calculation)

//********* Internal Parameter *********//
	map<ZAI ,double>	fValorisableIV;		///< The Valorisable Table
	map<int, cSecond >	fReactorNextStep;	///< Next Time Step to Build a New Fuel

	map< int,EvolutionData >	fReactorFuturDB; ///< List of the Futur EvolutionData use in the reactor
	map< int,IsotopicVector >	fReactorFuturIV; ///< List of the Futur Fuel Isotopic Vector used in the reactor

	Storage*	fStorage;			//!< Pointer to the Storage to recycle
	Storage*	fReUsable;			//!< Pointer to the Storage using for recycling unused Product

	vector< pair<int, double> >	fFractionToTake;	//!< The Temporary Storage IsotopicVector

		//	double		fFabricationTime;		///< Fabrication Duration Time
	bool		fChronologicalTimePriority;	//!< Set the Chronological Priotity (for the Stock Management) or the anti-chronological one

	bool		fSubstitutionFuel;		//!< true if a subtitution fuel as been set
	EvolutionData	fSubstitutionEvolutionData;	//!< EvolutionData of the subtitution fuel
	

//********* Private Method *********//
	void	FabricationPlantEvolution(cSecond t);				//!< Deal the FabricationPlant Evolution

	//{
	/// Separation Method
	/*!
	 Make the Separation
		\li IV[0] -> To Keep
		\li IV[1] -> To Waste
	 */
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector);
	//}
	ClassDef(FabricationPlant,2);

};

#endif
