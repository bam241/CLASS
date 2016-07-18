#ifndef _ISOTOPICVECTOR_
#define _ISOTOPICVECTOR_

/*!
 \file
 \brief Header file for IsotopicVector class. 
 @version 2.0
 */

#include "ZAI.hxx"

#include <TObject.h>

#include <iostream>
#include <map>
#include <algorithm>
#include <iterator>
#include <string>
#include <sstream>
#include <cmath>
#include <fstream>

typedef long long int cSecond;

//----------------------------------------------------------------------------
//! Allows to store & operate on radioactive sample

/*!
 Defines an Isotopicvector.
 An isotopicVector is a map of ZAI and double (e.g number of atoms).
 Its aim is to define a radioactive sample.
 
 @author BaM
 @author BLG
 @author Marc
 @version 2.1
*/
//____________________________________________________________________________
class IsotopicVector : public TObject
{
	public :
		typedef std::map<ZAI,double>::const_iterator	const_iterator;
		typedef std::map<ZAI,double>::iterator			iterator;

	/*!
		\name Constructor/Destructor
	!*/
	//@{
// CONSTRUCTOR ---------------------------------------------------------------
		IsotopicVector ();                         ///< Default constructor
		IsotopicVector ( const IsotopicVector & ); ///< Copy constructor

// DESTRUCTOR ----------------------------------------------------------------
		~IsotopicVector (); ///< Destructor
 
	//@}

// OPERATOR ------------------------------------------------------------------
	/*!
	 \name Internal Operator Method
	 */
	//@{
		IsotopicVector & operator = ( const IsotopicVector & ); ///< Assignment operator
		
		double   operator [] ( const ZAI & ) const; ///< Access operator
		double & operator [] ( const ZAI & );       ///< Access operator

		IsotopicVector & operator += ( const IsotopicVector & ); //!< Operator += definition
		IsotopicVector & operator -= ( const IsotopicVector & ); //!< Operator -= definition
		IsotopicVector & operator *= ( const IsotopicVector & ); //!< Operator *= definition
		IsotopicVector & operator *= ( double );                 //!< Operator *= definition (scalar)
	//@}

// GETTER --------------------------------------------------------------------
	/*!
	 \name Get Method
	 */
	//@{
		std::map<ZAI ,double> GetIsotopicQuantity       () const; //!< Return the IsotopicVector as a map
		std::map<ZAI ,double> GetIsotopicQuantityNeeded () const; //!< Return the needed IsotopicVector as a map

		IsotopicVector GetSpeciesComposition ( const int ) const; //!< Return the Species composition of the "z"
		IsotopicVector GetThisComposition ( const IsotopicVector & ) const; //!< Return the composition according the IV list...
		std::vector< ZAI > GetZAIList () const; //!< Return the list of ZAI present in the IV
		IsotopicVector GetActinidesComposition () const; //!< Return the Actinides composition (89 <= Z <= 103)

		double GetZAIIsotopicQuantity ( const short int , const short int , const short int ) const; ///< Return the quantity of the ZAI
		double GetZAIIsotopicQuantity ( const ZAI & ) const; ///< Return the quantity of the ZAI
		double GetQuantity ( const short int , const short int , const short int ) const;
		double GetQuantity ( const ZAI & ) const;

		double GetTotalMass     () const; //!< Return the mass (in tons) of the isotopic vector
		double GetMeanMolarMass () const; //<! Return the mean molar mass of the isotopic vector

		std::vector< int > GetChemicalSpecies () const; //!< Return the list of chemichal species contained

		double GetSumOfAll () const; //!< Return the Sum of nuclei in the IsotopicVector

		std::size_t size () const; ///< Get size of the list
		
		const_iterator begin () const; ///< Return an iterator at the begining of IsotopicVector
		iterator       begin (); ///< Return an iterator at the begining of IsotopicVector
		const_iterator end   () const; ///< Return an iterator at the end of IsotopicVector
		iterator       end   (); ///< Return an iterator at the end of IsotopicVector

		const_iterator find ( const ZAI & ) const; ///< Emulate map comportement with find
		iterator       find ( const ZAI & ); ///< Emulate map comportement with find
	//@}


// SETTER --------------------------------------------------------------------
	/*!
	 \name Set Method
	 */
	//@{
		void Initialize ( double );
		void Clear ();
		void ClearNeed();

		void Add ( const short int , const short int , const short int , double );
		void Add ( const ZAI & , double );
		void Add ( const std::pair<ZAI,double> & );
		void Add ( const IsotopicVector & );
		void Add ( const std::map<ZAI,double> & );

		void Remove ( const short int , const short int , const short int , double );
		void Remove ( const ZAI & , double );
		void Remove ( const std::pair<ZAI,double> & );
		void Remove ( const IsotopicVector & );

		void Multiply ( double );

		void ApplyZAIThreshold ( int );
	//@}

// METHOD --------------------------------------------------------------------
	/*!
	 \name  In/Out Method
	 */
	//@{
		void Write ( std::string filename , cSecond time = -1 ) const;
		void Print ( std::string o  = " ") const;
		std::string sPrint() const;
		void PrintList( std::string o  = " " ) const;
	//@}

	ClassDef(IsotopicVector,1);
// ATTRIBUT ------------------------------------------------------------------
	private :
		std::map< ZAI , double > fdata; // data of the IsotopicVector, this is isotopic quantities
		std::map< ZAI , double > fIsotopicQuantityNeeded; // quantity needed, never use, here just for retrocompatibility

};

/*!
 \name  Boolean IsotopicVector comparator
 */
//@{
bool operator == ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
bool operator != ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
bool operator <  ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
bool operator >  ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
bool operator <= ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
bool operator >= ( const IsotopicVector & , const IsotopicVector & ); //!< IsotopicVector Comparator
//@}


/*!
 \name  Aritmetic IsotopicVector comparator
 */
//@{
IsotopicVector operator / ( const IsotopicVector & , double );
IsotopicVector operator / ( const ZAI & , double ); ///< Build an IsotopicVector from a ZAI with an initial quantity
IsotopicVector operator * ( const IsotopicVector & , double );
IsotopicVector operator * ( const ZAI & , double ); ///< Build an IsotopicVector from a ZAI with an initial quantity
IsotopicVector operator * ( double , const IsotopicVector & );
IsotopicVector operator * ( double , const ZAI & ); ///< Build an IsotopicVector from a ZAI with an initial quantity

IsotopicVector operator + ( const IsotopicVector & , const IsotopicVector & );
IsotopicVector operator - ( const IsotopicVector & , const IsotopicVector & );
IsotopicVector operator * ( const IsotopicVector & , const IsotopicVector & );
//@}

std::ostream & operator << ( std::ostream & , const IsotopicVector & );

double 	RelativDistance ( const IsotopicVector & , const IsotopicVector & ); //!< return the euclidean distance between two IV. The two IV are normalize to unity
double 	Distance ( const IsotopicVector & , const IsotopicVector & , int DistanceType = 0 , const IsotopicVector & DistanceParameter = IsotopicVector() ); //!< return weighted euclidean distance between two IV
double	DistanceStandard( const IsotopicVector & , const IsotopicVector & ); //!< return the euclidean distance between two IV
double	DistanceAdjusted( const IsotopicVector & , const IsotopicVector & , const IsotopicVector & ); //!< return the weighted euclidean distance between two IV
double 	Norme( const IsotopicVector & , int DistanceType = 0 , const IsotopicVector & DistanceParameter = IsotopicVector() ); //!< return the norm of an IV

#endif
