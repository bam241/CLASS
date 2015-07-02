#include "EquivalenceModel.hxx"
#include "EQM_MLP_Keff_FBR.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

//________________________________________________________________________
//
//		EQM_MLP_Keff_FBR
//
//	Equivalenve Model based on multi layer perceptron from TMVA (root cern)
//	For FBR 
//
//________________________________________________________________________



//________________________________________________________________________
//
//	Objects & Methods for content prediction with calculation of Keff averaged
//	over batches.
//
//________________________________________________________________________

//________________________________________________________________________
EQM_MLP_Keff_FBR::EQM_MLP_Keff_FBR(string TMVAWeightPath,  int NumOfBatch, double LowerKeffective, double UpperKeffective, string InformationFile):EquivalenceModel(new CLASSLogger("EQM_MLP_Keff_FBR.log"))
{
	fIsAverageKeff = true;

	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;
	
	/*INFORMATION FILE HANDLING*/
	if(InformationFile=="")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMLPInformationFile = InformationFile;

	GetModelInformation();//Getting information from fMLPInformationFile

	/*OTHER MODEL PARAMETERS*/
	fNumberOfBatch = NumOfBatch;
	fKmin = LowerKeffective ;
	fKmax = UpperKeffective ;

	/*MODEL PARAMETERS INITIALIZATION */
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO<<"__An equivalence model has been define__"<<endl;
	INFO<<"\tThis model is based on the prediction of keff averaged over the number of batch"<<endl;
	INFO<<"\tThe TMVA (weight | information) files are :"<<endl;
	INFO<<"\t"<<"( "<<fTMVAWeightPath[0]<<" | "<<fMLPInformationFile<<" )"<<endl;
}
//________________________________________________________________________
EQM_MLP_Keff_FBR::EQM_MLP_Keff_FBR(CLASSLogger* log, string TMVAWeightPath,  int NumOfBatch, double LowerKeffective, double UpperKeffective, string InformationFile):EquivalenceModel(log)
{
	fIsAverageKeff = true;
	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;

	/*INFORMATION FILE HANDLING*/
	if(InformationFile=="")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMLPInformationFile = InformationFile;
	GetModelInformation();//Getting information from fMLPInformationFile

	/*OTHER MODEL PARAMETERS*/
	fNumberOfBatch = NumOfBatch;
	fKmin = LowerKeffective ;
	fKmax = UpperKeffective ;

	/*MODEL PARAMETERS INITIALIZATION */
	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO<<"__An equivalence model has been define__"<<endl;
	INFO<<"\tThis model is based on the prediction of keff averaged over the number of batch"<<endl;
	INFO<<"\tThe TMVA (weight | information) files are :"<<endl;
	INFO<<"\t"<<"( "<<fTMVAWeightPath[0]<<" | "<<fMLPInformationFile<<" )"<<endl;
}
//________________________________________________________________________
TGraph* EQM_MLP_Keff_FBR::BuildKeffGraph(IsotopicVector FreshFuel)
{
	if(!fIsAverageKeff)
		{ERROR<<" Can't be used with EQM_MLP_Keff_FBR(string TMVAWeightPath, double keff_target, string InformationFile) constructor";exit(1);}

	TGraph * keffGraph = new TGraph();
	for(int i = 0 ; i < (int) fMLP_Time.size() ; i++)
	{
		double keff_t = ExecuteTMVA( CreateTMVAInputTree(FreshFuel,(float) fMLP_Time[i]), true );
		keffGraph->SetPoint(i, (double)fMLP_Time[i], keff_t );
	}

	return keffGraph;
}
//________________________________________________________________________
TGraph* EQM_MLP_Keff_FBR::BuildAverageKeffGraph(TGraph* GRAPH_KEFF)
{
	if(!fIsAverageKeff)
		{ERROR<<" Can't be used with EQM_MLP_Keff_FBR(string TMVAWeightPath, double keff_target, string InformationFile) constructor";exit(1);}

	TGraph * AveragekeffGraph = new TGraph();
	int NumberOfPoint = 50;

	int NumOfInputGraphPoint = GRAPH_KEFF->GetN()-1 ;
	double TimeFinal = 0;
	double KFinal = 0;

	GRAPH_KEFF->GetPoint(NumOfInputGraphPoint, TimeFinal, KFinal);

	double step = TimeFinal/NumberOfPoint;
	int p=0;
 	for(int n = 0 ;n<NumberOfPoint;n++)
 	{	
 		double k_av = 0;
 		for(int b=0;b<fNumberOfBatch;b++)
 		{
 			if((step*p)> TimeFinal/(double)fNumberOfBatch)
 				p=0;
 			k_av +=	GRAPH_KEFF->Eval( (step*p + b*TimeFinal/(double)fNumberOfBatch) , 0 , "S" );

 		}
 		p++;
 		k_av/=(double)fNumberOfBatch;

 		//cout<< step*n<<" "<<k_av<<endl;
 		AveragekeffGraph->SetPoint(n, step*n, k_av);
 	}	

 	return AveragekeffGraph;
}
//________________________________________________________________________
double EQM_MLP_Keff_FBR::GetKeffAt(TGraph* GRAPH_KEFF, int Step)
{	
	if(!fIsAverageKeff)
		{ERROR<<" Can't be used with EQM_MLP_Keff_FBR(string TMVAWeightPath, double keff_target, string InformationFile) constructor";exit(1);}

	double Time = 0;
	double Keff=0;
	GRAPH_KEFF->GetPoint(Step, Time,  Keff);

	return Keff;
}
//________________________________________________________________________
double EQM_MLP_Keff_FBR::GetFissileContent_keffAveragePredictor(IsotopicVector Fissile, IsotopicVector Fertile)
{
	/**Algorithm  not so clever ...**/
	/**need improvements to make it faster*/

	if(!fIsAverageKeff)
		{ERROR<<" Can't be used with EQM_MLP_Keff_FBR(string TMVAWeightPath, double keff_target, string InformationFile) constructor";exit(1);}

	double FissileContent = 0.01;
	double test_Keff_beg = fKmin - 0.01;
	double test_keff_end = fKmax + 0.01;

	double speedstep = 1;

	while(fKmin >= test_Keff_beg || fKmax <=  test_keff_end)
	{	
		IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());

		TGraph* KEFF 	 = BuildKeffGraph(FreshFuel);
		TGraph* KEFF_avg = BuildAverageKeffGraph(KEFF);

		test_Keff_beg = GetKeffAt(KEFF_avg , 0);
		test_keff_end = GetKeffAt(KEFF_avg , fMLP_Time.size()-1);
		delete KEFF;
		delete KEFF_avg;

		if(test_Keff_beg < 0.9) //why 0.9 ? exactly
			speedstep = 10;
		else
			speedstep = 1;

		FissileContent+=0.001*speedstep;

		if( test_Keff_beg > 1.30 )
		{
			cout<<"This plutonium can not satisfy the criticality condition imposed"<<endl;
			FissileContent = -1;
			break;
		}

	}

	return FissileContent;
}

//________________________________________________________________________
//
//	Objects & Methods for content prediction with calculation of Keff at
//	one given time (often BOC or EOC)
//
//________________________________________________________________________

//________________________________________________________________________
EQM_MLP_Keff_FBR::EQM_MLP_Keff_FBR(string TMVAWeightPath, double keff_target, string InformationFile):EquivalenceModel(new CLASSLogger("EQM_MLP_Keff_FBR.log"))
{
	fIsAverageKeff = false;

	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;

	fTargetKeff = keff_target;

	if(InformationFile=="")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMLPInformationFile = InformationFile;
	GetModelInformation();//Getting information from fMLPInformationFile

	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO<<"__An equivalence model has been define__"<<endl;
	INFO<<"\tThis model is based on the prediction of keff at a specific time"<<endl;
	INFO<<"\tThe TMVA (weight | information) files are :"<<endl;
	INFO<<"\t"<<"( "<<fTMVAWeightPath[0]<<" | "<<fMLPInformationFile<<" )"<<endl;

}
//________________________________________________________________________
EQM_MLP_Keff_FBR::EQM_MLP_Keff_FBR(CLASSLogger* log, string TMVAWeightPath, double keff_target, string InformationFile):EquivalenceModel(log)
{
	fIsAverageKeff = false;

	/**The tmva weight*/
	fTMVAWeightPath = TMVAWeightPath;

	fTargetKeff = keff_target;

	if(InformationFile=="")
		InformationFile = StringLine::ReplaceAll(TMVAWeightPath,".xml",".nfo");

	fMLPInformationFile = InformationFile;
	GetModelInformation();//Getting information from fMLPInformationFile

	SetPCMprecision(10);
	SetBuildFuelFirstGuess(0.15);//First fissile content guess for the EquivalenceModel::BuildFuel algorithm
	fActualFissileContent = fFirstGuessFissilContent ;

	INFO<<"__An equivalence model has been define__"<<endl;
	INFO<<"\tThis model is based on the prediction of keff at a specific time"<<endl;
	INFO<<"\tThe TMVA (weight | information) files are :"<<endl;
	INFO<<"\t"<<"( "<<fTMVAWeightPath[0]<<" | "<<fMLPInformationFile<<" )"<<endl;

}
//________________________________________________________________________
double EQM_MLP_Keff_FBR::GetFissileMolarFraction_keffPredictor(IsotopicVector Fissile,IsotopicVector Fertile)
{
	if(fIsAverageKeff)
		{ERROR<<" Can't be used with EQM_MLP_Keff_FBR(string TMVAWeightPath,  int NumOfBatch, double LowerKeffective, double UpperKeffective, string InformationFile) constructor";exit(1);}

	//initialization
	double FissileContent = GetActualFissileContent(); 
	double OldFissileContentMinus	= 0;
	double OldFissileContentPlus	= fMaximalContent;
	double PredictedKeff = 0 ;
	IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
	double OldPredictedKeff =  GetKeffAtFixedTime(FreshFuel);

	double Precision = fPCMprecision/1e5*fTargetKeff; //pcm to 1 

	int count = 0;
	int MaximumLoopCount = 100;
	do	
	{
		if(count > MaximumLoopCount )
		{
			ERROR<<"CRITICAL ! Can't manage to predict fissile content\nHint : Try to decrease the precision on keff using :\nYourEQM_MLP_Keff_FBR->SetPCMPrecision(prop); with prop the precision  (default 0.5percent :  0.005) INCREASE IT\nIf this message still appear mail to leniau@subatech.in2p3.fr\nor nicolas.thiolliere@subatech.in2p3.fr "<<endl;
			exit(1);
		}

		if( (OldPredictedKeff - fTargetKeff) < 0 ) //The Content can be increased
		{
			OldFissileContentMinus = FissileContent;
			FissileContent = FissileContent + fabs(OldFissileContentPlus-FissileContent)/2.;
		}
		else if( (OldPredictedKeff - fTargetKeff) > 0) //The Content is too high
		{
			OldFissileContentPlus = FissileContent;
			FissileContent = FissileContent - fabs(OldFissileContentMinus-FissileContent)/2.;
		}

		IsotopicVector FreshFuel = (1-FissileContent)*(Fertile/Fertile.GetSumOfAll()) + FissileContent*(Fissile/Fissile.GetSumOfAll());
	//	FreshFuel.Print();
		PredictedKeff = GetKeffAtFixedTime(FreshFuel);
	//	cout<<PredictedKeff<<"  "<<FissileContent<<endl;;

		OldPredictedKeff = PredictedKeff;
		count ++;

	}while(fabs(fTargetKeff-PredictedKeff)>Precision);

	//cout<<"Predicted keff "<<PredictedKeff<<" FissileContent "<<FissileContent<<endl;;
return FissileContent;
}

//________________________________________________________________________
//
//	COMMON METHODS
//	
//
//________________________________________________________________________

//________________________________________________________________________
TTree* EQM_MLP_Keff_FBR::CreateTMVAInputTree(IsotopicVector TheFreshfuel, double ThisTime)
{
	/******Create Input data tree to be interpreted by TMVA::Reader***/
	TTree*   InputTree = new TTree("InTMPKef", "InTMPKef");

	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);

	float Time=0;

	IsotopicVector IVInputTMVA;
	map<ZAI ,string >::iterator it;
	int j=0;

	for( it = fMapOfTMVAVariableNames.begin()  ; it!=fMapOfTMVAVariableNames.end() ; it++)
	{
		InputTree->Branch( ((*it).second).c_str()	,&InputTMVA[j], ((*it).second + "/F").c_str());
		IVInputTMVA+= ((*it).first)*1;
		j++;
	}

	if(ThisTime!=-1)
		InputTree->Branch(	"Time"	,&Time	,"Time/F"	);

	IsotopicVector IVAccordingToUserInfoFile = TheFreshfuel.GetThisComposition(IVInputTMVA);

	double Ntot = IVAccordingToUserInfoFile.GetSumOfAll();

	IVAccordingToUserInfoFile = IVAccordingToUserInfoFile/Ntot;

	j=0;
	map<ZAI ,string >::iterator it2;

	for( it2 = fMapOfTMVAVariableNames.begin() ; it2!=fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTMVA[j] = IVAccordingToUserInfoFile.GetZAIIsotopicQuantity( (*it2).first ) ;
		//DBGV((*it2).first.Z()<<" "<<(*it2).first.A()<<" "<<InputTMVA[j]);
		j++;
	}

	Time=ThisTime;

	InputTree->Fill();

	return InputTree;

}
//________________________________________________________________________
double EQM_MLP_Keff_FBR::ExecuteTMVA(TTree* InputTree, bool IsTimeDependent)
{

	// --- Create the Reader object
	TMVA::Reader *reader = new TMVA::Reader( "Silent" );

	// Create a set of variables and declare them to the reader
	// - the variable names MUST corresponds in name and type to those given in the weight file(s) used
	vector<float> 	InputTMVA;
	for(int i = 0 ; i< (int)fMapOfTMVAVariableNames.size() ; i++)
		InputTMVA.push_back(0);
	Float_t Time;

	map<ZAI ,string >::iterator it;
	int j=0;
	for( it = fMapOfTMVAVariableNames.begin()  ; it!=fMapOfTMVAVariableNames.end() ; it++)
	{	reader->AddVariable( ( (*it).second ).c_str(),&InputTMVA[j]);
		j++;
	}

	if(IsTimeDependent)
		reader->AddVariable( "Time" ,&Time);

	// --- Book the MVA methods

	// Book method MLP
	TString methodName = "MLP method";
	reader->BookMVA( methodName, fTMVAWeightPath );

	map<ZAI ,string >::iterator it2;
	j=0;
	for( it2 = fMapOfTMVAVariableNames.begin()  ; it2!=fMapOfTMVAVariableNames.end() ; it2++)
	{
		InputTree->SetBranchAddress(( (*it2).second ).c_str(),&InputTMVA[j]);
		j++;
	}

	if(IsTimeDependent)
		InputTree->SetBranchAddress( "Time" ,&Time  );

	InputTree->GetEntry(0);
	Float_t val = (reader->EvaluateRegression( methodName ))[0];

	delete reader;

	return (double)val;//return k_{eff}(t=Time)
}
//________________________________________________________________________
void EQM_MLP_Keff_FBR::GetModelInformation()
{
DBGL
	ifstream FILE(fMLPInformationFile.c_str());
	if(FILE.good())
	{
		while(!FILE.eof())
		{
			string line;
			getline(FILE, line);  
			if(line=="") //ie if eof
				break;

			size_t foundSpecPow	= line.find("Specific Power (W/gHM) :");
			size_t foundContent	= line.find("Maximal fissile content (molar proportion) :");
			size_t foundTime    = line.find("Time (s) :");
			size_t foundZAI	  	= line.find("Z A I Name (input MLP) :");
			size_t foundFissile = line.find("Fissile Liste (Z A I) :");
			size_t foundFertile = line.find("Fertile Liste (Z A I Default Proportion) :");

			int pos=0;
			if(foundSpecPow !=std::string::npos)
			{	StringLine::NextWord(line,pos,':');
				fSpecificPower = atof( (StringLine::NextWord(line,pos,':') ).c_str() );
			}
			pos=0;
			if(foundContent!=std::string::npos)
			{
				StringLine::NextWord(line,pos,':');
				fMaximalContent = atof( (StringLine::NextWord(line,pos,':') ).c_str() );
			}
			pos=0;
			if(foundTime!=std::string::npos)
			{
				StringLine::NextWord(line,pos,':');
				while( pos< (int)line.size() )
					fMLP_Time.push_back( atof( (StringLine::NextWord(line,pos,' ')).c_str() ));
			}
 			pos=0;
			if(foundZAI != std::string::npos)
			{	string Z;
				string A;
				string I;
				string Name;
				int posoflinebeforbadline=0;
				do
				{	posoflinebeforbadline = FILE.tellg();
					getline(FILE, line);

					stringstream ssline;
					ssline<<line;
					ssline>>Z>>A>>I>>Name;
					if(StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I) )
					{	
						fMapOfTMVAVariableNames.insert( pair<ZAI,string>(ZAI(atoi(Z.c_str()),atoi(A.c_str()),atoi(I.c_str())),Name) );
					}

				}while((StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I)) && !FILE.eof());

				FILE.seekg(posoflinebeforbadline); //return one line before

			}
 			pos=0;			
			if(foundFissile != std::string::npos)
			{	string Z;
				string A;
				string I;

				int posoflinebeforbadline=0;
				do
				{	posoflinebeforbadline = FILE.tellg();
					getline(FILE, line);

					stringstream ssline;
					ssline<<line;
					ssline>>Z>>A>>I;
					if(StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I) )
						fFissileList.Add(atoi(Z.c_str()),atoi(A.c_str()),atoi(I.c_str()),1.0);
					

				}while((StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I)) && !FILE.eof());

				FILE.seekg(posoflinebeforbadline); //return one line before

			}
			pos=0;
			if(foundFertile != std::string::npos)
			{	string Z;
				string A;
				string I;
				string prop;
				int posoflinebeforbadline=0;
				do
				{	posoflinebeforbadline = FILE.tellg();
					getline(FILE, line);

					stringstream ssline;
					ssline<<line;
					ssline>>Z>>A>>I>>prop;
					if(StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I) )
						fFertileList.Add(atoi(Z.c_str()),atoi(A.c_str()),atoi(I.c_str()),atof(prop.c_str()) );
					
				}while((StringLine::IsDouble(Z) && StringLine::IsDouble(A) && StringLine::IsDouble(I)) && !FILE.eof());

				FILE.seekg(posoflinebeforbadline); //return one line before
			}

		}
	}
	else
	{
		ERROR << "Can't find/open file " << fMLPInformationFile << endl;
		exit(1);
	}

	/********DEBUG*************************************/
	INFO<<"\tMLP kinf Information : "<<endl;
	INFO<<"\t\tSpecific Power (W/gHM) :"<<fSpecificPower<<endl;	
	if(fIsAverageKeff)
	{	if((int)fMLP_Time.size()==0)
		{	ERROR<<" Time vector is empty"<<endl;
			exit(1);
		}	
		INFO<<"\t\tTime (s) :"<<endl;
		for (int i = 0; i < (int)fMLP_Time.size(); ++i)
		INFO<<"\t\t\t"<<fMLP_Time[i]<<endl;
	}
	else 
		INFO<<"\t\tMaximal fissile content (molar proportion) :"<<fMaximalContent<<endl;

	INFO<<"\t\tZ A I Name (input MLP) :"<<endl;	
	map<ZAI ,string >::iterator it;
	for (it= fMapOfTMVAVariableNames.begin();it!=fMapOfTMVAVariableNames.end();it++)
		INFO<<"\t\t\t"<< it->first.Z()<<" "<<it->first.A()<<" "<<it->second<<endl;
	INFO<<"Fissile Liste (Z A I) :"<<endl;
	INFO<<fFissileList.sPrint()<<endl;
	INFO<<"Fertile Liste (Z A I Default Proportion) :"<<endl;
	INFO<<fFertileList.sPrint()<<endl;
DBGL
}
//________________________________________________________________________
double EQM_MLP_Keff_FBR::GetFissileMolarFraction(IsotopicVector Fissile,IsotopicVector Fertile,double TargetBU)
{
	if(TargetBU != 0)
		WARNING<<"The third arguement : Burnup has no effect here.";

	double FissileContent = 0;

	if(fIsAverageKeff)
		FissileContent = GetFissileContent_keffAveragePredictor(Fissile, Fertile);

	else
		FissileContent = GetFissileMolarFraction_keffPredictor(Fissile, Fertile);


return FissileContent;
}
//________________________________________________________________________
