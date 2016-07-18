
#ifndef _XSM_CLOSEST_HXX
#define _XSM_CLOSEST_HXX


/*!
 \file
 \brief Header file for XSM_MLP_PWR_MOX class.


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
//! Defines a XSModel getting mean cross sections from the closest EvolutionData 

/*!
 Define a XSM_CLOSEST.
 Class to get cross sections from a set of pre-calculation
 (with MURE,or other depletion code).
With this class, cross sections needed to solves Bateman equation come from an
 already performed depletion calculation contained in a set of many depletion calculations. 
 The way to pick up these or these depletion calculation is related to a distance :
 The depletion calculation (or EvolutionData) the closest from the new fresh fuel composition
 is selected. Different distances are defined : 
 
 \li Standard euclidean distance (weights : 1) IS THE DEFAULT
 \li Euclidean distance with weights assigned according mean cross section values
 \li Euclidean distance with weights assigned by user

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
	
	//{
	///  normal constructor
	//
	/*!
	 Make a new XSM_CLOSEST
	 \param DB_index_file : File listing the path of all depletion calculations (in EvolutionData format (.dat file) )
	 \param oldreadmethod :
	 */
	XSM_CLOSEST(string DB_index_file, bool oldreadmethod = false );
	//}
	
	//{
	///  Logger constructor
	//
	/*!
	 Make a new XSM_CLOSEST
	 \param log : Use for the log
	 \param DB_index_file : File listing the path of all depletion calculations (in EvolutionData format (.dat file) )
	 \param oldreadmethod :
	 */
	XSM_CLOSEST(CLASSLogger* Log, string DB_index_file, bool oldreadmethod = false );
	//}
	
	~XSM_CLOSEST();
	//@}


	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
 	virtual EvolutionData GetCrossSections( const IsotopicVector & isotopicvector,double t = 0) ;	//!< Reason to live of this CLASS : Return the closest Evolutiondata
	vector< EvolutionData >	GetFuelDataBank()	const	{ return fFuelDataBank; }	//!< Return the FuelDataBank
	string 			GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index Name
	string			GetFuelType()		const	{ return fFuelType; }		//!< Return the fuel type of the DB
	vector<double>		GetFuelParameter()	const	{ return fFuelParameter; }	//!< Return the Fuel parameters of the DB
	pair<double,double>	GetBurnUpRange()	const	{ return fBurnUpRange;}		//!< Return the Burnup range of the DB
	bool 			IsDefine(IsotopicVector IV)	const;				//!< True if the key is define, false instead

	map<double, int>	GetDistancesTo( const IsotopicVector & isotopicvector, double t = 0);	//! Return a map containing the distance of each EvolutionData in the DataBase to the set IV at the time t
	//@}

	//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetFuelDataBank(vector< EvolutionData > mymap)	{ fFuelDataBank = mymap; }	//!< Set the FuelDataBank map

	void SetDataBaseIndex(string database)	{ fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the Name of the database index
	void SetOldReadMethod(bool val)		{ fOldReadMethod = val; ReadDataBase();}	//!< use the old reading method



	void SetWeightedDistanceCalculation(bool val = true) { fWeightedDistance = val;}		//!< Set weighted distance calculation
	void SetInventoryEvolutionInterpolation(bool val = true) { fEvolutionDataInterpolation = val;}	//!< \deprecated The 64 closest EvolutionData (ED) are used to build a new ED. The cross section from each ED are weighted according to their distance
	void SetDistanceParameter( const IsotopicVector & DistanceParameter);					//!< Define mannually the weight for each ZAI in the distance calculation
	//@}

	//{
	/// Choose the way to calculate the distance between two isotopic vectors.
	/*!
	 // The different algorythms are:
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
	void	ReadDataBase();				//!< read the index file and fill the evolutionData map
	void	CalculateDistanceParameter();		//!< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the FuelDataBank

	//@}

	private :

	vector< EvolutionData > fFuelDataBank;			//!< DataBank map

 	string			fDataBaseIndex;			//!< Name of the index

	bool			fOldReadMethod;			//!< use old DB format
	bool			fWeightedDistance;		//!< USe XS weighted distance calculation
	bool			fEvolutionDataInterpolation;	//!< USe XS weighted distance calculation


 	string 				fFuelType;	//!< Type of fuel of the FuelDataBank
 	pair<double,double>	fBurnUpRange;		//!< Range of the Burn-up range of the FuelDataBank
 	vector<double>		fFuelParameter;		//!< Parameter needed by the equivalence model



	int 	fDistanceType;		//!< Set the distance calculation algorythm
	/// \li 0 is for the standard norm (Default = 0),
	/// \li 1 for each ZAI weighted with its XS,
	/// \li 2 for each ZAI weighted with coefficient given by the user.

	IsotopicVector		fDistanceParameter;	//!< weight for each ZAI in the distance calculation
	
};

#endif

