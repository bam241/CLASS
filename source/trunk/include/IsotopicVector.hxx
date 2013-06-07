#ifndef __ISOTOPICVECTOR_HXX__
#define __ISOTOPICVECTOR_HXX__


/*!
 \file
 \brief Header file for IsotopicVector class. 
  The aim of this Class is to manage any kind of IsotopicVector, and any operation between them : sum, substraction....
  
 @author BaM, Marc
 @version 2.0
 */
#include "ZAI.hxx"

#include "TObject.h"
#include <string>
#include <map>

using namespace std;
typedef long long int cSecond;

class IsotopicVector;





class IsotopicVector : public TObject
{
public :

//********* Constructor/Destructor Method *********//
	///< Normal Constructor.
 	IsotopicVector();
 	
 	///< Normal Destructor.
 	~IsotopicVector();	 

	
//********* Get Method *********//

	map<ZAI ,double>	GetIsotopicQuantity() const { return fIsotopicQuantity; }			//!< Return the IVQuantity map 
	map<ZAI ,double>	GetIsotopicQuantityNeeded() const { return fIsotopicQuantityNeeded; }	//!< Return the IVQuantityNeeded map
	IsotopicVector		GetSpeciesComposition(int z) const;					//!< Return the Species composition of the "z" atom
	IsotopicVector		GetActinidesComposition() const;					//!< Return the Actinides composition of the "z" atom
	double	GetZAIIsotopicQuantity(const ZAI& zai) const;						///< Return the composition of the IsotopicVector
	double	GetZAIIsotopicQuantity(const int z, const int a, const int i) const;			///< Return the composition of the IsotopicVector

	vector<int>		GetChemicalSpecies() const;						//!< Return the Species Species contained

	
//********* Modification Method *********//
	void 	Clear();					//!< Empty all the IV 
	void 	ClearNeed();					//!< Empty Need componant of the IV 

	void	Add(const ZAI& zai, double quantity); 		//!< Add Quantity gramme of the ZAI Element
	void	Add(const IsotopicVector& isotopicvector); 	//!< Add IsotopicVector to the existing IsotopicVector
	void	Add(const map<ZAI ,double>& quantity);		//!< Add IsotopicVector to the existing IsotopicVector

	void	Need(const ZAI& zai, double quantity);		//!< Fill the fIsotopicQuantityNeeded
	void	Need(const IsotopicVector& isotopicvector);	//!< Fill the fIsotopicQuantityNeeded
	void	Need(const map<ZAI ,double>& quantityneeded) { fIsotopicQuantityNeeded = quantityneeded; }	
								//!< Fill the fIsotopicQuantityNeeded

	void	Remove(const ZAI& zai, double quantity); 	//!< Remove Quantity gramme of the ZAI Element
	void	Remove(const IsotopicVector& isotopicvector); 	//!< Remove IsotopicVector to the existing IsotopicVector

	void 	Multiply(double factor);			//!< Multiply the IV by a Factor
	
	void	Write(string filename, cSecond time = -1 ) const;	///< Write the Content of the IV in the filename file

	void	Print(string o =" ") const ;			///< Print the composition of the IV
	

	
	
//******* Set Operator between IsotopicVector *******//

	IsotopicVector& operator+=(IsotopicVector const& IVb);	//!<....
	IsotopicVector& operator-=(IsotopicVector const& IVb);	//!<....
	bool operator <(const IsotopicVector& isotopicvector) const;	//!< IsotopicVector Comparator

//***************************************************///< 

	
	protected :

	map<ZAI ,double>	fIsotopicQuantity;		///< Isotopic vector composition in Atome Number
	map<ZAI ,double>	fIsotopicQuantityNeeded;	///< Isotopic vector request and not present

	ClassDef(IsotopicVector,2);
};

IsotopicVector operator/(IsotopicVector const& IVA, double F);
IsotopicVector operator/(ZAI const& zai, double F);
IsotopicVector operator*(IsotopicVector const& IVA, double F);
IsotopicVector operator*(ZAI const& zai, double F);
IsotopicVector operator*(double F, IsotopicVector const& IVA);
IsotopicVector operator*(double F, ZAI const& zai);
IsotopicVector operator+(IsotopicVector const& IVa, IsotopicVector const& IVb);
IsotopicVector operator-(IsotopicVector const& IVa, IsotopicVector const& IVb);


double 	RelativDistance(IsotopicVector IV1, IsotopicVector IV2 );
double 	Distance(IsotopicVector IV1, IsotopicVector IV2 ,int DistanceType=0, IsotopicVector DistanceParameter=IsotopicVector());
double	DistanceStandard(IsotopicVector IV1, IsotopicVector IV2);
double	DistanceAdjusted(IsotopicVector IV1, IsotopicVector IV2, IsotopicVector DistanceParameter);
double 	Norme(IsotopicVector IV1,int DistanceType=0, IsotopicVector DistanceParameter=IsotopicVector());


#endif
