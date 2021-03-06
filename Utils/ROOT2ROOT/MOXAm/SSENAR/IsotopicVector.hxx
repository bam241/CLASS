#ifndef _ISOTOPICVECTOR_
#define _ISOTOPICVECTOR_


/*!
 \file
 \brief Header file for IsotopicVector class. 
 @version 2.0
 */
#include "ZAI.hxx"

#include "TObject.h"
#include <string>
#include <vector>
#include <map>

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
//! Allows to store & operate on radioactive sample

/*!
 Defines an Isotopicvector.
 An isotopicVector is a map of ZAI and double (e.g number of atoms).
 Its aim is to define a radioactive sample.
 
 @author BaM
 @author BLG
 @author Marc
 @version 2.0
 */
//________________________________________________________________________



class IsotopicVector : public TObject
{
	public :
		typedef std::map<ZAI,double>::const_iterator	const_iterator;
		typedef std::map<ZAI,double>::iterator			iterator;

	//********* Constructor/Destructor Method *********//
	/*!
	 \name Constructor/Desctructor
	 */
	//@{

		IsotopicVector();	///< Normal Constructor.
		IsotopicVector ( IsotopicVector const& IVa ); ///< Copy Constructor.


		~IsotopicVector();	///< Normal Destructor.

	//@}

	//********* Operator Method *********//
	/*!
	 \name Internal Operator
	 */
	//@{

		IsotopicVector& operator += ( IsotopicVector const& IVb );	//!< Operator += definition
		IsotopicVector& operator -= ( IsotopicVector const& IVb );	//!< Operator -=  definition
		IsotopicVector& operator *= ( IsotopicVector const& IVb );	//!< Operator *=  definition
		IsotopicVector& operator *= ( double const& factor );	//!< Operator *=  definition (scalar)
		bool operator < (const IsotopicVector& isotopicvector) const;	//!< IsotopicVector Comparator
	
	//@}

	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{

		map<ZAI,double> GetIsotopicQuantity () const
							{ return fIsotopicQuantity; }		//!< Return the IsotopicVector as a map
		map<ZAI,double> GetIsotopicQuantityNeeded () const
							{ return fIsotopicQuantityNeeded; }	//!< Return the needed IsotopicVector as a map 
		IsotopicVector GetSpeciesComposition ( int z ) const;   //!< Return the Species composition of the "z" atom
		IsotopicVector GetThisComposition    ( IsotopicVector IV ) const; //!< Return the composition according the IV list...
		vector<ZAI>    GetZAIList () const; //!< Return the list of ZAI present in the IV
		IsotopicVector GetActinidesComposition() const; //!< Return the Actinides composition (Z >= 89)

		double GetZAIIsotopicQuantity ( const ZAI& zai ) const; ///< Return the quantity of the ZAI
		double GetZAIIsotopicQuantity ( const int z, const int a, const int i) const; ///< Return the quantity of the ZAI
		double GetQuantity ( const int z, const int a, const int i ) const { return GetZAIIsotopicQuantity(z,a,i); }
		double GetQuantity ( const ZAI& zai ) const { return GetZAIIsotopicQuantity(zai); }
		
		void Initiatlize ( double val );
		
		
		double	GetTotalMass     () const; //!< Return the mass (in tons) of the isotopic vector
		double	GetMeanMolarMass () const; //<! Return the mean molar mass of the isotopic vector
		
		vector<int> GetChemicalSpecies () const; //!< Return the list of chemichal species contained
		int         GetZAIQuantity     () const
							{return  fIsotopicQuantity.size(); } //!< Return the number of different ZAI in the IsotopicVector

		double      GetSumOfAll () const; //!< Return the Sum of nuclei in the IsotopicVector

		iterator       begin ()
			{ return fIsotopicQuantity.begin(); } //!< Return an iterator on the begin of fIsotopicQuantity
		const_iterator begin () const
			{ return fIsotopicQuantity.begin(); } //!< Return a constant iterator on the begin of fIsotopicQuantity
		iterator       end   ()
			{ return fIsotopicQuantity.end();   } //!< Return an iterator on the end of fIsotopicQuantity
		const_iterator end   () const
			{ return fIsotopicQuantity.end();   } //!< Return a constant iterator on the end of fIsotopicQuantity

		iterator       find ( const ZAI & zai )       { return fIsotopicQuantity.find(zai); }
		const_iterator find ( const ZAI & zai ) const { return fIsotopicQuantity.find(zai); }
	//@}




//*********  Internal Operation Method *********//

	/*!
	 \name Internal Operation Method
	 */
	//@{


	void 	Clear();					//!< Empty all the IV
	void 	ClearNeed();					//!< Empty Need componant of the IV 

	void	Add(const ZAI& zai, double quantity); 		//!< Add Quantity gramme of the ZAI Element
	void	Add(const IsotopicVector& isotopicvector); 	//!< Add IsotopicVector to the existing IsotopicVector
	void	Add(const map<ZAI ,double>& quantity);		//!< Add IsotopicVector to the existing IsotopicVector
	void	Add(int Z, int A, int I, double quantity)
		{ (*this).Add(ZAI(Z,A,I), quantity); } 		//!< Add Quantity gramme of the ZAI Element
	void  Add ( const pair<ZAI,double> & zaiQ ) { Add( zaiQ.first , zaiQ.second ); }


	void	Need(const ZAI& zai, double quantity);		//!< Fill the fIsotopicQuantityNeeded
	void	Need(const IsotopicVector& isotopicvector);	//!< Fill the fIsotopicQuantityNeeded
	void	Need(const map<ZAI ,double>& quantityneeded) { fIsotopicQuantityNeeded = quantityneeded; }	
								//!< Fill the fIsotopicQuantityNeeded

	void	Remove(const ZAI& zai, double quantity); 	//!< Remove Quantity gramme of the ZAI Element
	void	Remove(const IsotopicVector& isotopicvector); 	//!< Remove IsotopicVector to the existing IsotopicVector
	void Remove ( const pair<ZAI,double> & zaiQ ) { Remove( zaiQ.first , zaiQ.second ); }

	void 	Multiply(double factor);			//!< Multiply the IV by a Factor

	void	ApplyZAIThreshold(int z = 90);			//!< Put all nuclei below the threshold in -2 -2 -2 ZAI...
	//@}



//********* In/Out related Method *********//

	/*!
	 \name  In/Out Method
	 */
	//@{

	void	Write(string filename, cSecond time = -1 ) const;	//!< Write the Content of the IV in the filename file

	void	Print(string o  = " ") const ;				//!< Print the composition of the IV in terminal
	string	sPrint() const ;					//!< Print the composition of the IV in a string

	void	PrintList(string o  = " ") const ;			//!< Print the composition of the IV

	//@}
	

//***************************************************///< 

	
	protected :

	map<ZAI ,double>	fIsotopicQuantity;		///< Isotopic vector composition in atomes Number
	map<ZAI ,double>	fIsotopicQuantityNeeded;	///< map where negative value are saved

	ClassDef(IsotopicVector,1);
};

IsotopicVector operator/(IsotopicVector const& IVA, double F);
IsotopicVector operator/(IsotopicVector const& IVA, IsotopicVector const& IVB);
IsotopicVector operator/(ZAI const& zai, double F);
IsotopicVector operator*(IsotopicVector const& IVA, double F);
IsotopicVector operator*(ZAI const& zai, double F);
IsotopicVector operator*(double F, IsotopicVector const& IVA);
IsotopicVector operator*(double F, ZAI const& zai);
IsotopicVector operator+(IsotopicVector const& IVa, IsotopicVector const& IVb);
IsotopicVector operator-(IsotopicVector const& IVa, IsotopicVector const& IVb);

IsotopicVector operator*(IsotopicVector const& IVa, IsotopicVector const& IVb);


double 	RelativDistance ( const IsotopicVector & , const IsotopicVector & ); //!< return the euclidean distance between two IV. The two IV are normalize to unity
double 	Distance ( const IsotopicVector & , const IsotopicVector & , int DistanceType = 0 , const IsotopicVector & DistanceParameter = IsotopicVector() ); //!< return weighted euclidean distance between two IV
double	DistanceStandard( const IsotopicVector & , const IsotopicVector & ); //!< return the euclidean distance between two IV
double	DistanceAdjusted( const IsotopicVector & , const IsotopicVector & , const IsotopicVector & ); //!< return the weighted euclidean distance between two IV
double 	Norme( const IsotopicVector & , int DistanceType = 0 , const IsotopicVector & DistanceParameter = IsotopicVector() ); //!< return the norm of an IV


#endif
