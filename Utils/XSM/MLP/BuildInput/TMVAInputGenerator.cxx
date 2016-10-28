/**********************************************************/
//	Make the input file for the MLPs training
//
//	This programs reads a set of .dat files which are the
//	results of a depletion calculation (see manual and 
//  looks for XS_CLOSEST). From the reading it fills a 
//  TTree (Data) and write it in a file named
//	TrainingInput.root . 
//	The file TrainingInput.cxx is the list of MLP outputs 
// (cross sections)
//
//@author BaM, BaL
/**********************************************************/
#include "TMVAInputGenerator.hxx"
#include <TH1F.h>
#include <TH2D.h>
#include <TFile.h>
#include <TTree.h>
#include "../../../../source/external/StringLine.hxx"
#include <TString.h>
#include <string>
#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <sstream>
#include <numeric>
#include <functional>
#include <algorithm>    
#include <ctime>

using namespace std;

ZAIMass cZAIMass; //atomic masses

string ElNames[110]={"  ","H","He","Li","Be",
						"B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
						"K","Ca","Sc","Ti","V ","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge",
						"As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd",
						"Ag","Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
						"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta","W",
						"Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn","Fr","Ra",
						"Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",
						"Lr","Rf","Db","Sg","Bh","Hs","Mt"};

//--------------------------------------------------------------------------------------------------
/*************************
 		MAIN
*************************/
int main(int argc, char ** argv){
	
	if(argc!=2)
	{
		cout << "Usage : TMVAInputGenerator Path" << endl;
		cout << " Where Path is the path to the folder containing the Evolution Datas" << endl;
		cout << " i.e the (.dat) files" << endl;
		exit(0);
	}	

	fEvolutionDataFolder = argv[1];

	CheckJob();	// looks fot the .dat files in the fEvolutionDataFolder

	cout << "Load your EvolutionDatas to R.A.M" << endl;
	cout<<endl;

	for(int i = 0; i < (int)JobName.size(); i++)
	{
		ReadAndFill(JobName[i]);
		ProgressBar(i,JobName.size());
	}

	FillMapName();
	DumpInputNeuron("TrainingInput.root");
	CreateInfoFile();

	cout << "╭───────────────────────────────────────────────╮" << endl; 
	cout << "│                GENERATED FILES:               │" << endl; 
	cout << "├───────────────────────────────────────────────┤" << endl; 
	cout << "│#1 Input for TMVA training: \033[36mTrainingInput.root\033[0m │" << endl; 
	cout << "│#2 Target names for TMVA:   \033[36mTrainingInput.cxx\033[0m  │" << endl; 
	cout << "│#3 Model Information for CLASS:                │" << endl; 
	cout << "│                           \033[36mData_Base_Info.nfo\033[0m  │" << endl; 
	cout << "╰───────────────────────────────────────────────╯" << endl; 
	cout << endl;
	cout << "╭───────────────────────────────────────────────╮" << endl; 
	cout << "│                NEXT STEPS:                    │" << endl; 
	cout << "├───────────────────────────────────────────────┤" << endl; 
	cout << "│1. Train your MLPs with ../Train/Train_XS.cxx  │" << endl; 
	cout << "│2. Test MLPs performance using information in: │" << endl; 
	cout << "│	     ../Train/EvaluateTrainingCommands.dat   │" << endl; 
	cout << "│3. Put the file #3 in ../Train/weights then    │" << endl; 
    cout << "│ mouve this folder to $CLASS_PATH/DATA_BASES   │" << endl; 
	cout << "╰───────────────────────────────────────────────╯" << endl; 

	cout << "=> Do you want me to train the MLPs on your local machine on one cpu ? It can take a (huge) while.  [y/n]" << endl; 

	if(UserSayYes())
	{

	}
	else
	{
		cout<<" Ok , I suggest you run on a grid " <<endl;
		cout << "You can proceed like so and run this script on  as many  NumberOfProcessor that you want:"<<endl;
		cout << "#bin/bash"<<endl;
		cout<< " for ((i=0 ; NumberOfReaction/NumberOfProcessor  ; i++)) "<<endl;
		cout<<"    do Train_XS($i)  "<<endl;
		cout<<"  done"

	}


}
//--------------------------------------------------------------------------------------------------
// Function definitions
//--------------------------------------------------------------------------------------------------
//Convert int to string
string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os << num;
	return os.str();
}
//--------------------------------------------------------------------------------------------------
//Give a name to a cross section 
//!!!!!!!!!!!!!!!!!!!!!!!!!!//
// YOU MUST KEEP THE FORMAT :
// XS_Z_A_I_fis	(fission cross section)
// XS_Z_A_I_cap (n,gamma cross section)
// XS_Z_A_I_n2n	(n,2n cross section)
//!!!!!!!!!!!!!!!!!!!!!!!!!!//
string NameXS(ZAI act,string xs)
{
	stringstream Name;
	Name<<"XS_"<<act.Z<<"_"<<act.A<<"_"<<act.I<<"_"<<xs;

	return Name.str();
}

//--------------------------------------------------------------------------------------------------
void  FillMapName()
{
	cout<<endl;
	cout<< "╭───────────────────────────────────────────────╮"<<endl;
	cout<< "│  Looking for fresh fuel composition @ t=0  to │"<<endl;
	cout<< "│  set up TMVA MLPs inputs                      │"<<endl;
	cout<< "╰───────────────────────────────────────────────╯"<<endl;
	cout<<endl;

	vector<ZAI> ZAI_T0 =  fActinideCompoInit[0].GetNonZeroZAIList();

	for(int zai = 0 ; zai < ZAI_T0.size() ; zai++ )
	{	stringstream ssname;
		int Z = ZAI_T0[zai].Z;
		int A = ZAI_T0[zai].A;
		int I = ZAI_T0[zai].I;
		string sI="";
		if(I == 1)
			sI = "m";
		else if(I == 2) 
			sI = "2m";
		else if(I == 3) 
			sI = "3m";

		ssname<<A<<ElNames[Z]<<sI;
		if(zai == 0)
			cout<< "Add this nuclei with this name [y/n] ?  (if you don't know say yes to all)  "<<endl;
		else
			cout<< "Add [y/n] ?  "<<endl;

		cout<< "[Z\tA\tI\tName]"<<endl;
		cout<<"\033[36m"<<Z<<"\t"<<A<<"\t"<<I<<"\t"<<ssname.str()<<"\033[0m"<<endl;

		if (UserSayYes())
			fMapName.insert(pair<ZAI,string> ( ZAI(Z,A,I) , ssname.str()) );

	}

	bool UserWantToAdd = true;
	while(UserWantToAdd)
	{
		cout<< "Do you want to add additional nuclei  [y/n] ? (if you don't know say no) "<<endl;
		if(UserSayYes())
		{	
			int Z = 0;
			int A = 0;
			int I = 0;
			string Name;
			cout << "Z -> "; cin >> Z ; cout <<" A -> "; cin >> A ; cout <<" I -> "; cin >> I ;cout <<" Name -> "; cin >> Name;
			fMapName.insert(pair<ZAI,string> ( ZAI(Z,A,I) ,Name) );
		}

		else
			UserWantToAdd = false;
	}	

}
//--------------------------------------------------------------------------------------------------
bool UserSayYes()
{
		bool AnswerIsNotGiven = true;
		bool isYES = false;
		while(AnswerIsNotGiven)
		{
			string answer;
			std::getline(std::cin, answer);

	
			if(answer == "y" || answer == "yes" || answer == "Yes" || answer == "Y")
			{
				isYES = true;
				AnswerIsNotGiven = false;
			}
			else if (answer == "n" || answer == "no" || answer == "No" || answer == "N")
			{
				isYES = false;
				AnswerIsNotGiven = false;
			}

			else{
				cout << "Yes OR No ???!"<<endl;
			}
		}
	return isYES;	
}

//--------------------------------------------------------------------------------------------------
void CreateInfoFile()
{

	/**************************************************/
	// GETTING USER INFO
	/**************************************************/
	double sum = 0 ;
	for (int i=0 ; i < fHMMass.size() ; i++)
			sum+=fHMMass[i];

    double MeanHMMass = sum / (double)fHMMass.size();
    string ReactorType,FuelType,Author,Mail,XSBase,HLCut,EnergyDisc,FPYBase,SABase,Geom,AddInfo,DepCode ;
    double Power = 0;


	int start=0;
	string AnInfoFile = JobName[0];
	int pos=AnInfoFile.find(".dat",start);
	AnInfoFile.replace(pos,4,".info");

    ReadInfo(AnInfoFile,ReactorType,FuelType,Power);

	cout<<endl;
	cout<<endl;
	cout<<"╭───────────────────────────────────────────────╮"<<endl;
	cout<<"│        XSM_MLP .NFO FILE GENERATOR            │"<<endl;
	cout<<"╰───────────────────────────────────────────────╯"<<endl;
	cout<<" Answer following questions "<<endl;
	cout<<"-> Author(s) name(s) : "<<endl;
	std::getline(std::cin, Author);
	cout<<endl;

	cout<<"-> email adress(es) : "<<endl;
	std::getline(std::cin, Mail);	
	cout<<endl;

	cout<<"-> Depletion code used : "<<endl;
	std::getline(std::cin, DepCode);
	cout<<endl;

	cout<<"-> Cross section data base (e.g ENSDF7.1) : "<<endl;
	std::getline(std::cin, XSBase);
	cout<<endl;

	cout<<"-> Fission yield data base (e.g ENSDF7.1) : "<<endl;
	std::getline(std::cin, FPYBase);
	cout<<endl;

	cout<<"-> S(alpha,beta) data base (e.g ENSDF7.1) : "<<endl;
	std::getline(std::cin, SABase);
	cout<<endl;

	cout<<"-> Geometry simulated  (e.g Cubic Assembly with mirror boundary) : "<<endl;
	std::getline(std::cin, Geom);
	cout<<endl;

	cout<<"-> Half life cut [s] (if any) : "<<endl;
	std::getline(std::cin, HLCut);	
	cout<<endl;

	cout<<"-> Multi group treatment (yes/no if yes give the number of groups) : "<<endl;
	std::getline(std::cin, EnergyDisc);
	cout<<endl;	

	cout<<"-> Additional informations : "<<endl;
	std::getline(std::cin, AddInfo);	
	cout<<endl;	


	cout<<"-> Reactor type (e.g PWR, FBR,...) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<ReactorType<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "<<endl;
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> ReactorType;
	}	
	cout<<endl;

	cout<<"-> Fuel type (e.g UOX, MOX,...) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<FuelType<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "<<endl;
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> FuelType;
	}	
	cout<<endl;

	cout<<"-> Simulated heavy metal mass (tons) : "<<endl;
	cout<< "\t According your evolution datas the AVERAGE heavy metal mass is : "<<endl;
	cout<<"\033[36m"<<MeanHMMass<<"\033[0m"<<" tons"<<endl;
	cout<< "Is that corect ? [y/n] "<<endl;
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> MeanHMMass;
	}
	cout<<endl;

	cout<<"-> Simulated thermal power (W) : "<<endl;
	cout<<"Found in a .info file :"<<endl;
	cout<<"\033[36m"<<Power<<"\033[0m"<<endl;
	cout<< "Is that corect ? [y/n] "<<endl;
	if(!UserSayYes())
	{	cout<<"\t So what it is ?"<<endl;
		cin>> Power;
	}	
	cout<<endl;

	/**************************************************/
	// BUILDING FILE
	/**************************************************/
	ofstream InfoFile("Data_Base_Info.nfo");

  	InfoFile <<"============================================" << endl;
  	InfoFile <<"    Informations needed by XSM_MLP model    " <<endl;
  	InfoFile <<"============================================" << endl;
	InfoFile << endl;
	InfoFile << "Reactor Type :"<<endl;
	InfoFile << "K_REACTOR "<< ReactorType <<endl;
	InfoFile << endl;
	InfoFile << "Fuel Type :"<<endl;
	InfoFile << "K_FUEL "<< FuelType <<endl;
	InfoFile << endl;
	InfoFile << "Heavy Metal [t] :"<<endl;
	InfoFile << "K_MASS  "<< MeanHMMass <<endl;
	InfoFile << endl;
	InfoFile << "Thermal Power [W] :"<<endl;
	InfoFile << "K_POWER  "<< Power <<endl;
	InfoFile << endl;
	InfoFile << "Irradiation time steps [s] :"<<endl;
	InfoFile << "K_TIMESTEP";
		for( int t = 0 ; t < fNOfTimeStep ; t++ )
			InfoFile <<" "<<fTime[t];
	InfoFile << endl<<endl;
	InfoFile << "Z A I Name (input MLP) :"<<endl;
	for(map<ZAI,string>::iterator it = fMapName.begin() ; it != fMapName.end() ; it++ )
			InfoFile <<"K_ZAINAME "<<it->first.Z <<" " <<it->first.A <<" " <<it->first.I<<" " << it->second <<endl ;
	InfoFile <<endl;
	InfoFile << "Fuel range (Z A I min max) :"<<endl;
	for(map<ZAI,string>::iterator it = fMapName.begin() ; it != fMapName.end() ; it++ )
	{	vector <double> AllCompoOfZAI = GetAllCompoOf(it->first);
		vector <double>::iterator Min = std::min_element(AllCompoOfZAI.begin(),AllCompoOfZAI.end()); 
		vector <double>::iterator Max = std::max_element(AllCompoOfZAI.begin(),AllCompoOfZAI.end());
		InfoFile <<"K_ZAIL "<<it->first.Z <<" " <<it->first.A <<" " <<it->first.I<<" " << *Min << " "  << *Max <<endl ;
  	}
    InfoFile << endl;
  	InfoFile <<"============================================" << endl;
  	InfoFile <<"     Data base generation informations      " <<endl;
  	InfoFile <<"============================================" << endl;
  	InfoFile << endl;
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t ); 
    InfoFile <<" Date: "<< now->tm_mday<< '/'<<  (now->tm_mon + 1) << '/'  <<(now->tm_year + 1900)<< endl;
	InfoFile <<" Author(s): "<< Author <<endl;
	InfoFile <<" Author(s) contact: "<< Mail <<endl;
	InfoFile <<" Depletion code: "<< DepCode <<endl;
	InfoFile <<" Simulated geometry: "<<  Geom <<endl;
	InfoFile <<" Nuclear data used in "<< DepCode <<endl;
	InfoFile <<"\tCross section library: "<< XSBase <<endl;
	InfoFile <<"\tFission yield library: "<< FPYBase <<endl;
	InfoFile <<"\tS(alpha,beta) library: "<< SABase <<endl;
	InfoFile <<" Half life cut [s] : "<<  HLCut <<endl;
	InfoFile <<" Multi-group treatment: "<<  EnergyDisc <<endl;
	InfoFile <<" Additional informations: "<< endl << AddInfo <<endl;


}
//--------------------------------------------------------------------------------------------------
vector<double> GetAllCompoOf(ZAI zai)
{
	vector<double> AllCompoOfZAI;

		for( int b = 0 ; b < fActinideCompoInit.size() ; b++ )
			AllCompoOfZAI.push_back(fActinideCompoInit[b].GetZAIIsotopicQuantity(zai)); 

return AllCompoOfZAI;
}
//--------------------------------------------------------------------------------------------------
void ProgressBar(double loopindex, double totalindex)
{
	// Reset the line
	for(int i = 0; i < 10; i++)
		cout << "  ";
	cout << flush ;

	cout << "\r[\033[42m";
	for(int i = 0; i < (int)(loopindex/totalindex*20.0); i++)
		cout << " ";
	cout<<"\033[0m";
	for(int i = 20; i >= (int)(loopindex/totalindex*20.0); i--)
		cout << " ";
	cout << "] ";

	cout << (int)(loopindex/totalindex*100) << "%\r";
	cout << flush;
	//cout << endl;

}
//--------------------------------------------------------------------------------------------------
void ReadInfo(string InfoDBFile,string &ReactorType,string &FuelType,double &Power)
{	
	ifstream InfoDB(InfoDBFile.c_str());				// Open the File
	if(!InfoDB)
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
	}
	else
	{
		int start = 0;
		string line;
		getline(InfoDB, line);
		string Next = StringLine::NextWord(line, start, ' ');
		StringLine::ToLower(Next) ;
		if ( Next ==  "reactor")
			ReactorType = StringLine::NextWord(line, start, ' ');
		
		start = 0;
		getline(InfoDB, line);
		Next = StringLine::NextWord(line, start, ' ');
		StringLine::ToLower(Next) ;
		if ( Next ==  "fueltype")
			FuelType = StringLine::NextWord(line, start, ' ');
	
		start = 0;
		getline(InfoDB, line);
		Next = StringLine::NextWord(line, start, ' ');
		StringLine::ToLower(Next) ;
		if ( Next  ==  "cycletime")
			double cycletime = atof(StringLine::NextWord(line, start, ' ').c_str());
	
		getline(InfoDB, line); // Assembly HM Mass DONT TRUST THIS ONE CALCULATED WITH A instead of real atomic mass
	
		start = 0;
		getline(InfoDB, line);
		Next = StringLine::NextWord(line, start, ' ');
		StringLine::ToLower(Next) ;
		if ( Next ==  "constantpower")
			Power = atof(StringLine::NextWord(line, start, ' ').c_str());
		InfoDB.close();
	}

}


//--------------------------------------------------------------------------------------------------
void DumpInputNeuron(string filename)
{
	TFile*   fOutFile = new TFile(filename.c_str(),"RECREATE");
	TTree*   fOutT = new TTree("Data", "Data");


/**********************INITIALISATIONNN********************/

	////////////////////////////////////////////////////////
	// INIT FRESH FUEL COMPOSITION and TIME
	////////////////////////////////////////////////////////

	double *FreshCompo = new double[fMapName.size()]; 

	for(int i = 0 ; i < fMapName.size() ; i++ )
		FreshCompo[i] = 0;

	////////////////////////////////////////////////////////
	double Time 		     = 0;	

/**********************init map********************/
	map < ZAI,vector<double> > mAllXS;
	map < ZAI, vector<double>  > mAllInventories;
	for(int act=0;act<int(fAllNuclei.size());act++ )	
	{	 
		vector<double>  InitVect;
		for(int Tstep=0 ;Tstep<fNOfTimeStep;Tstep++)
		{
			InitVect.push_back(0);
		}
		mAllInventories.insert( pair<ZAI,vector<double> >(fAllNuclei[act],InitVect) );
	}

	for(int act=0;act<int(fAllNuclei.size());act++ )	
	{	
		vector< double>  InitVect;
		for(int xs=0;xs<3;xs++)
		{
			
			InitVect.push_back(0);
				
		}

		mAllXS.insert(pair<ZAI,vector< double> >(fAllNuclei[act], InitVect) );
	}


/**********************BRANCHING**************************************************/
/**********************Fresh fuel**************************************************/

	map<ZAI,string>::iterator it;
	int index = 0;
	for(it = fMapName.begin() ; it != fMapName.end() ; it++ )
	{	string Name = it->second + "/D"; 
		fOutT->Branch( it->second.c_str() , &FreshCompo[index] , Name.c_str());
		index++;
	}

	////////////////////////////////////////////////////////
	fOutT->Branch(	"Time"			,&Time			,"Time/D"			);


/**********************cross section**************************************************/

	string XSType[3]={"fis","cap","n2n"}; //!!!!!!DO NOT TOUCH THIS

	for(int act=0;act<int(fAllNuclei.size());act++ )
	{	
		for(int xs=0;xs<3;xs++)
		{
				string NamedXS= NameXS(fAllNuclei[act],XSType[xs]);
				string NameXSBis=NamedXS+"/D";

		if(	xs==0 && fXSFis[0][fAllNuclei[act]].size()!=0)	//Check if the reaction we want to branch exists
			 if(fXSFis[0][fAllNuclei[act]][0]!=0) 
			 		fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);			

		if(xs==1 && fXSCap[0][fAllNuclei[act]].size()!=0)
			 if(fXSCap[0][fAllNuclei[act]][0]!=0) 
			 		fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);				

		if(xs==2 && fXSN2N[0][fAllNuclei[act]].size()!=0)
			 if(fXSN2N[0][fAllNuclei[act]][0]!=0) 
					fOutT->Branch(	NamedXS.c_str() ,  &mAllXS[fAllNuclei[act]][xs],   NameXSBis.c_str()	);				


		}
	}
/**********************FILLING THE TTREE**************************************************/
		cout<<endl;
		cout<<endl;
		cout<<"╭───────────────────────────────────────────────╮"<<endl;
		cout<<"│                 FILLING TTREE                 │"<<endl;
		cout<<"│         (building TrainingInput.root)         │"<<endl;		         
		cout<<"╰───────────────────────────────────────────────╯"<<endl;
		cout<<endl;
	//File containing all the output of the networks to train
	 ofstream  InputNetwork("TrainingInput.cxx");

	 int NumOfBase=fActinideCompoInit.size();
	for(int b=0;b<NumOfBase;b++) 
	{ 

		ProgressBar(b,NumOfBase);

		///////////////////////////////////////////////////////
		int index =0;
		for(it = fMapName.begin() ; it != fMapName.end() ; it++ )
		{	
			int Z = it->first.Z;
			int A = it->first.A;
			int I = it->first.I;
			FreshCompo[index] = fActinideCompoInit[b].GetZAIIsotopicQuantity(ZAI(Z,A,I));

			index++;
		}	
		///////////////////////////////////////////////////////

			for(int Tstep=0 ;Tstep<fNOfTimeStep;Tstep++ )	
			{	
	 			Time=fTime[Tstep];
				for(int act=0;act<int(fAllNuclei.size());act++)
				{
	
					if(fXSFis[b][fAllNuclei[act]].size()!=0) // reaction may not be present
					{	
						if(fXSFis[b][fAllNuclei[act]][Tstep] !=0 )
						{	mAllXS[fAllNuclei[act]][0]   = 	fXSFis[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)
								InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[0])<<"\");"<<endl;
						}
					}		
					if(fXSCap[b][fAllNuclei[act]].size()!=0)
				 	{	
				 		if(fXSCap[b][fAllNuclei[act]][Tstep] !=0) 
				 		{	mAllXS[fAllNuclei[act]][1]   = 	fXSCap[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)		 			
				 				InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[1])<<"\");"<<endl;
						}
				 	}	
				 	if( fXSN2N[b][fAllNuclei[act]].size()!=0)
				 	{	
				 		if(fXSN2N[b][fAllNuclei[act]][Tstep] !=0)
				 		{	mAllXS[fAllNuclei[act]][2]   = 	fXSN2N[b][fAllNuclei[act]][Tstep];
							if(b==0 && Tstep==0)
				 				InputNetwork<<"OUTPUT.push_back(\""<<NameXS(fAllNuclei[act],XSType[2])<<"\");"<<endl;
						}
				 	}
				 	
				 } 
				 fOutT->Fill();	
			}		
	}


	fOutFile->Write();
	delete fOutT;
	fOutFile-> Close();
	delete fOutFile;

}

//--------------------------------------------------------------------------------------------------
void CheckJob()
{	//LOAD THE LIST OF EvolutionData

		cout<<endl;
		cout<<"╭───────────────────────────────────────────────╮"<<endl;
		cout<<"    Scanning :                                   "<<endl;
		cout<< "  " << fEvolutionDataFolder                     <<endl;	
		cout<<"          for EvolutionData (.dat files)         "<<endl;
		cout<<"╰───────────────────────────────────────────────╯"<<endl;

	string Command = "find "+ fEvolutionDataFolder + " -name \"*.dat\" > JOB.tmp";
	system(Command.c_str());
	
	ifstream JOB("JOB.tmp");
	if (JOB.is_open())
	{
		while (!JOB.eof())
		{
			string tmp;
			getline(JOB,tmp);
			JobName.push_back(tmp);
		}
	} JOB.close(); JobName.pop_back();
	
	// Remove temporary files...
	Command = "\\rm -f JOB.tmp";
	system(Command.c_str());
	random_shuffle(JobName.begin(), JobName.end());
	cout << "Scan complete" <<endl; 
}

//--------------------------------------------------------------------------------------------------
void ReadAndFill(string jobname)
{	//Read a .dat file and fill XS maps and the fuel initial composition

	vector<double>	vT;
	vector<int>	Z;
	vector<int>	A;
	vector<int>	I;
	vector<double>	Q;
	
	ifstream DecayDB(jobname.c_str());							// Open the File
	if(!DecayDB)
	{
		cout << "!!Warning!! !!!EvolutiveProduct!!! \n Can't open \"" << jobname << "\"\n" << endl;
	}
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	
	/******Getting Time vecotr ....******/
	if( StringLine::NextWord(line, start, ' ') != "time")
	{
		cout << "!!Bad Trouble!! !!!EvolutiveProduct!!! Bad Database file : " <<  jobname << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vT.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	fNOfTimeStep=int(vT.size());

	/****Getting Inventories***/	
	getline(DecayDB, line);
	do
	{	
		start = 0;
		int z;
		string tmp2 = StringLine::NextWord(line, start, ' ');
		if (tmp2 == "Inv") {
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
		}
		else z = atoi(tmp2.c_str());
		int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
		int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
		if(a!=0 && z!=0)
		{
			
			ZAI zaitmp(z, a, i);
			Z.push_back(z);
			A.push_back(a);
			I.push_back(i);
			if(!fIsAllNucleiAlreadyFill)
			{	
				fAllNuclei.push_back(zaitmp);
				fTime=vT;	
			}	

			long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
			Q.push_back(q);	

		}

		getline(DecayDB, line);
		start = 0;
		tmp2 = StringLine::NextWord(line, start, ' ');
		
		if(line == "" || line == "CrossSection" || tmp2 == "XSFis" || tmp2 == "XSCap" || tmp2 == "XSn2n") break;
	}while (!DecayDB.eof() );

	if(fAllNuclei.size()!=0)
	 fIsAllNucleiAlreadyFill=true;

	//XS  
	 map<ZAI, vector <double> > mapFistmp;
	 map<ZAI, vector <double> > mapCaptmp;
	 map<ZAI, vector <double> > mapN2Ntmp;
	do
	{
		
		start = 0;
		int z;
		string tmp2 = StringLine::NextWord(line, start, ' ');
		if (tmp2 == "XSFis") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{
				
				ZAI zaitmp(z, a, i);

				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapFistmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));		
			}
		}
		else if (tmp2 == "XSCap") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{
				
				ZAI zaitmp(z, a, i);
				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapCaptmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));	
				
			}
		}	
		else if (tmp2 == "XSn2n") 
		{
			z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int a = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int i = atoi(StringLine::NextWord(line, start, ' ').c_str());
		
			if(a!=0 && z!=0)
			{		
				ZAI zaitmp(z, a, i);
				vector<double> XSTime;
				for(int i = 0; i < (int)vT.size(); i++)
				{	long double q = atof(StringLine::NextWord(line, start, ' ').c_str());
					XSTime.push_back(q);
				}
				mapN2Ntmp.insert( pair<ZAI,vector<double> >(zaitmp,XSTime));					
			}
		}
		getline(DecayDB, line);
		start = 0;
		tmp2 = StringLine::NextWord(line, start, ' ');
		
		if(line == "") break;
	}while (!DecayDB.eof() );
	
	fXSFis.push_back(mapFistmp);
	fXSCap.push_back(mapCaptmp);
	fXSN2N.push_back(mapN2Ntmp);
	
	DecayDB.close();
	
	double N = 0;
	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(  Z[i]>89 ) 
			N += Q[i];
	}
	

	IsotopicVector CompoBasei; 
	IsotopicVector CompoBaseiUnormalize;

	for(int i=0; i < (int)Z.size()-2; i++)
	{
		if(Z[i]>89)
		{
			ZAI zai = ZAI(Z[i], A[i], I[i]);
			CompoBasei.IVquantity.insert(pair<ZAI,double>(zai,Q[i]/N));
			CompoBaseiUnormalize.IVquantity.insert(pair<ZAI,double>(zai,Q[i]));
		}	
	}
	
	fActinideCompoInit.push_back(CompoBasei);
	double MassOfThisFuel = cZAIMass.GetMass(CompoBaseiUnormalize);
	//cout<<MassOfThisFuel<<endl;
	fHMMass.push_back(MassOfThisFuel);

GoodJobName.push_back(jobname);

}
/*-------------------------------------------------------------------------------------------------
COMPILATION :

g++ -o TMVAInputGenerator TMVAInputGenerator.cxx `root-config --cflags` `root-config --libs`

*/