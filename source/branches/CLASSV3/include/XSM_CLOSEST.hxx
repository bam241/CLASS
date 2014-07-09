
#ifndef _XSM_CLOSEST_HXX
#define _XSM_CLOSEST_HXX


/*!
 \file
 \brief Header file for XS_CLOSEST class.
 
 
 @authors BaM,BLG
 @version 1.0
 */
#include "XSModel.hxx"
#include <string>
#include <fstream>
#include <iostream> 
#include <map>
#include <vector>
typedef long long int cSecond;
using namespace std;

//-----------------------------------------------------------------------------//
/*!
 Define a XSM_CLOSEST.
CLASS to get cross sections from a set of pre-calculation 
(with MURE,or other depletion code) 
get cross sections of the closest (in composition) calculation 

 @authors BaM,BLG
 @version 1.0
 */
//________________________________________________________________________


class XSM_CLOSEST : public XSModel
{

public : 

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

	XSM_CLOSEST(LogFile* Log, string DB_index_file, bool oldreadmethod = true );	
	~XSM_CLOSEST();
	//{


	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
 	EvolutionData GetCrossSections(IsotopicVector isotopicvector,double t=0) ; //!< Reason to live of this CLASS Return the closest Evolutiondata
	map<IsotopicVector ,EvolutionData >	GetFuelDataBank()	const	{ return fFuelDataBank; }	//!< Return the FuelDataBank
	string 					GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index Name
	string					GetFuelType()		const	{ return fFuelType; }		//!< Return the fuel type of the DB
	vector<double>			GetFuelParameter()	const	{ return fFuelParameter; }	//!< Return the Fuel parameter of the DB
	pair<double,double>		GetBurnUpRange()	const	{ return fBurnUpRange;}		//!< Return the BurnUp range of the DB
	bool 					IsDefine(IsotopicVector IV)	const;					//!< True the key is define, false unstead

	map<double, EvolutionData>	GetDistancesTo(IsotopicVector isotopicvector, double t = 0) const;	//! Return a map containing the distance of each EvolutionData in the DataBase to the set IV at the t time
	//@}

	//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetFuelDataBank(map< IsotopicVector ,EvolutionData > mymap)	{ fFuelDataBank = mymap; }	//!< Set the FuelDataBank map

	void SetDataBaseIndex(string database) { fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the Name of the database index

	void SetOldReadMethod(bool val)			{ fOldReadMethod = val; ReadDataBase();}			///< use the old reading method



	void SetWeightedDistanceCalculation(bool val = true) { fWeightedDistance = val;}		///< Set weighted Distance calculation
	void SetEvolutionDataInterpolation(bool val = true) { fEvolutionDataInterpolation = val;}		///< Set weighted Distance calculation
	void SetDistanceParameter(IsotopicVector DistanceParameter);		///< Define mannually the weight for each ZAI in the distance calculation


	//{
	/// Define the way to decide if two isotopic vectors are close.
	/*!
	 // The different algorythm are:
	 // \li 0 is for the standard norme,
	 // \li 1 for each ZAI weighted with its XS,
	 // \li 2 for each ZAI weighted with coefficient given by the user.
	 */
	void SetDistanceType(int DistanceType);
	//}

	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void	ReadDataBase();				///< read the index file and fill the evolutionData map
	void	CalculateDistanceParameter();		///< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the FuelDataBank

	//@}

 private :

	map<IsotopicVector, EvolutionData>	fFuelDataBank;		///< DataBanck map
	map<IsotopicVector, EvolutionData>	fFuelDataBankCalculated;	///< Map of the already calculated EvolutionData (to avoid recalculation...)

 	string			fDataBaseIndex;			///< Name of the index

	bool			fOldReadMethod;		///< use old DB format
	bool			fWeightedDistance;	///< USe XS weighted distance calculation
	bool			fEvolutionDataInterpolation;	///< USe XS weighted distance calculation


 	string 				fFuelType;		///< Type of fuel of the FuelDataBank
 	pair<double,double>	fBurnUpRange;		///< Range of the Burn-up range of the FuelDataBank
 	vector<double>		fFuelParameter;		///< Parameter needed by the equivalence model



	int 	fDistanceType;		///< Set the distance calculation algorytm
	/// \li 0 is for the standard norm (Default = 0),
	/// \li 1 for each ZAI weighted with its XS,
	/// \li 2 for each ZAI weighted with coefficient given by the user.

	IsotopicVector		fDistanceParameter;	///< weight for each ZAI in the distance calculation

};

#endif

