#include "EQM_FBR_BakerRoss_MOX.hxx"
#include "CLASSLogger.hxx"

#include <vector>

//________________________________________________________________________
//
//		EQM_FBR_BakerRoss_MOX
//
//	Equivalenve Model based on Plutonium 239 equivalent  
//	using the formula of Baker&Ross
//
//________________________________________________________________________


EQM_FBR_BakerRoss_MOX::EQM_FBR_BakerRoss_MOX(double Weight_U_235, double Weight_Pu_238, double Weight_Pu_240, double Weight_Pu_241, double Weight_Pu_242, double Weight_Am_241, double EquivalentFissile):EquivalenceModel(new CLASSLogger("EQM_FBR_BakerRoss_MOX.log"))
{
	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich = 0.0025;

	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	FertileList = U5*U5_enrich + U8*(1-U5_enrich); //Default Fertile composition (if no Fertile Storage is set for the FabricationPlant  )	
	FissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1; //ZAI to be extracted by the FabricationPlant in its affiliated Fissile Storage 
	fStreamList["Fissile"] = FissileList;
	fStreamList["Fertile"] = FertileList;

	fReferenceFissilContent = EquivalentFissile;//Equivalent fissile content

	fFirstGuessContent["Fissile"] = fReferenceFissilContent;
	fFirstGuessContent["Fertile"] = 1 - fFirstGuessContent["Fissile"];

	fWeight_U_235  = Weight_U_235 ;
	fWeight_Pu_238 = Weight_Pu_238;
	fWeight_Pu_240 = Weight_Pu_240;
	fWeight_Pu_241 = Weight_Pu_241;
	fWeight_Pu_242 = Weight_Pu_242;
	fWeight_Am_241 = Weight_Am_241;


	INFO << "__An equivalence model of FBR MOX has been defined__" << endl;
	INFO << "\tThis model is based on Plutonium 239 equivalent" << endl;
	INFO << "\t\t Weigt values : " << endl;
	INFO << "\t\t Fertile :  " << endl;
	INFO << "\t\t\tU_235 " <<   Weight_U_235 << endl;
	INFO << "\t\t\tU_238 0 (by definition)" << endl;
	INFO << "\t\t Fissile :  " << endl;
	INFO << "\t\t\tPu_238 " << Weight_Pu_238 << endl;
	INFO << "\t\t\tPu_239 1 (by definition)" << endl;
	INFO << "\t\t\tPu_240 " << Weight_Pu_240 << endl;
	INFO << "\t\t\tPu_241 " << Weight_Pu_241 << endl;
	INFO << "\t\t\tPu_242 " << Weight_Pu_242 << endl;
	INFO << "\t\t\tAm_241 " << Weight_Am_241 << endl;

}
//________________________________________________________________________
EQM_FBR_BakerRoss_MOX::EQM_FBR_BakerRoss_MOX(CLASSLogger* log, double Weight_U_235, double Weight_Pu_238, double Weight_Pu_240, double Weight_Pu_241, double Weight_Pu_242, double Weight_Am_241, double EquivalentFissile ):EquivalenceModel(log)
{
	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich = 0.0025;
	
	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);

	FertileList = U5*U5_enrich + U8*(1-U5_enrich); 	//Default Fertile composition (if no Fertile Storage is set for the FabricationPlant  )	
	FissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1; 	//ZAI to be extracted by the FabricationPlant in its affiliated Fissile Storage 
	fStreamList["Fissile"] = FissileList;
	fStreamList["Fertile"] = FertileList;

	fReferenceFissilContent = EquivalentFissile;		//Equivalent fissile content

	fFirstGuessContent["Fissile"] = fReferenceFissilContent;
	fFirstGuessContent["Fertile"] = 1 - fFirstGuessContent["Fissile"];

	fWeight_U_235  = Weight_U_235 ;
	fWeight_Pu_238 = Weight_Pu_238;
	fWeight_Pu_240 = Weight_Pu_240;
	fWeight_Pu_241 = Weight_Pu_241;
	fWeight_Pu_242 = Weight_Pu_242;
	fWeight_Am_241 = Weight_Am_241;


	INFO << "__An equivalence model of FBR MOX has been defined__" << endl;
	INFO << "\tThis model is based on Plutonium 239 equivalent" << endl;
	INFO << "\t\t Weigt values : " << endl;
	INFO << "\t\t Fertile : " << endl;
	INFO << "\t\t\tU_235 " <<   Weight_U_235 << endl;
	INFO << "\t\t\tU_238 0 (by definition)" << endl;
	INFO << "\t\t Fissile : " << endl;
	INFO << "\t\t\tPu_238 " << Weight_Pu_238 << endl;
	INFO << "\t\t\tPu_239 1 (by definition)" << endl;
	INFO << "\t\t\tPu_240 " << Weight_Pu_240 << endl;
	INFO << "\t\t\tPu_241 " << Weight_Pu_241 << endl;
	INFO << "\t\t\tPu_242 " << Weight_Pu_242 << endl;
	INFO << "\t\t\tAm_241 " << Weight_Am_241 << endl;


}
//________________________________________________________________________
map < string , double> EQM_FBR_BakerRoss_MOX::GetMolarFraction(map < string , IsotopicVector> IVStream, double BurnUp)
{

	if( BurnUp != 0 )
		WARNING << "Burn up (third argument) has no effect " << endl;
	
	IsotopicVector Fissile = IVStream["Fissile"];
	IsotopicVector Fertile = IVStream["Fertile"];

	double FissileContent = 0.;       

	IsotopicVector FissileListPlusDecay;
	FissileListPlusDecay.Add(94,238,0,1);
	FissileListPlusDecay.Add(94,239,0,1);
	FissileListPlusDecay.Add(94,240,0,1);
	FissileListPlusDecay.Add(94,241,0,1);
	FissileListPlusDecay.Add(94,242,0,1);
	FissileListPlusDecay.Add(95,241,0,1);

	IsotopicVector FertileList;
	FertileList.Add(92,238,0,1);
	FertileList.Add(92,239,0,1);

	//Getting the fissile from the Fissile input & normalize it
	IsotopicVector FissileFromInput = Fissile.GetThisComposition(FissileListPlusDecay);
	FissileFromInput = FissileFromInput/FissileFromInput.GetSumOfAll();

	//Getting the fissile from the Fissile input & normalize it
	IsotopicVector FertileFromInput = Fertile.GetThisComposition(FertileList);
	FertileFromInput = FertileFromInput/FertileFromInput.GetSumOfAll();

	double  SumWeightNFissile = fWeight_Pu_238 * FissileFromInput.GetZAIIsotopicQuantity(94,238,0) 
								+ 1		 * FissileFromInput.GetZAIIsotopicQuantity(94,239,0)
								+ fWeight_Pu_240 * FissileFromInput.GetZAIIsotopicQuantity(94,240,0) 
								+ fWeight_Pu_241 * FissileFromInput.GetZAIIsotopicQuantity(94,241,0) 
								+ fWeight_Pu_242 * FissileFromInput.GetZAIIsotopicQuantity(94,242,0) 
								+ fWeight_Am_241 * FissileFromInput.GetZAIIsotopicQuantity(95,241,0) ;

	double  SumWeightNFertile = fWeight_U_235 * FertileFromInput.GetZAIIsotopicQuantity(92,235,0);

	FissileContent = (fReferenceFissilContent - SumWeightNFertile )/(SumWeightNFissile-SumWeightNFertile); //Baker & Ross formula

	map < string , double> MolarFraction;
	
	MolarFraction["Fissile"] = FissileContent;
	MolarFraction["Fertile"] = 1.- FissileContent;

	return MolarFraction; //return Molar content of each component in the fuel
}
