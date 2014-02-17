#include <stdio.h>
#include <stdlib.h>
#include <istream>
#include <math.h>
#include <vector.h>
#include <string>
#include <TMath.h>
#include <TGraph.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TPaveLabel.h>
#include <TF2.h>
#include <TMinuit.h>
#include <TVirtualFitter.h>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <TColor.h>
#include <TGraph2D.h>
#include <TColor.h>
#include <TCutG.h>
#include <TColor.h>
#include <map>

using namespace std;

map<string,int> ReactorIndex;
map<string,int> PoolIndex;
map<string,int> StorageIndex;
map<string,int> FabricationPlantIndex;

string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os<<setprecision(3)<<num;
	return os.str();
}


vector<ZAI> PrintZAIList(TTree *T)
{
	IsotopicVector IVTot;

	IsotopicVector *IVIn=0;
	T->SetBranchStatus("TOTAL.",1);
	T->SetBranchAddress("TOTAL.", &IVIn);

	Long64_t nentries = T->GetEntries();

	T->GetEntry(nentries-1);
	IVTot = (*IVIn);


	vector<ZAI> ZAIvector = IVTot.GetZAIList();
	return ZAIvector;
}


void ReadPower(TTree *T, char* opt = "L*")
{

	T->Draw("ParcPower*1e-9:AbsTime/3600./24./365.25","AbsTime/3600./24./365.25 < 1e3",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle("Power [GWth]");
	htemp->SetTitle("Electronuclear fleet power");
	T_g->Draw(opt);


}

void PrintNames(TTree *T, string option, bool PrintResult=true)
{
	if(option==""){
		if(PrintResult){
		cout<<"Avaliable options :"<<endl;
		cout<<"\t reactor : print the names of the reactors"<<endl;
		cout<<"\t pool : print the names of the cooling pools"<<endl;
		cout<<"\t storage : print the names of the storage facilities"<<endl;
		cout<<"\t fabrication : print the names of the fabrication plants"<<endl;
		cout<<"\t all : print all of the above"<<endl;}
	}
	else{
	if(option=="all"||option=="reactor"){
		if(PrintResult)cout<<"Reactors :"<<endl;
	
		Int_t nBranches =  T->GetNbranches();
		int place=0;
		for(Int_t i=0;i<nBranches;i++){
			string Name;
			Name=T->GetListOfBranches()->At(i)->GetName();
			Name=Name.substr(0,7);
			if(Name=="Reactor"){
				Reactor *reactor;
				TString Rname;
				Rname.Form("Reactor%d.",place);
		
				T->SetBranchStatus(Rname,1);
				T->SetBranchAddress(Rname, &reactor);
	
				T->GetEntry(i);
				ReactorIndex[reactor->GetName()]=place;
				if(PrintResult)cout<<"\t "<<Rname<<"\t"<<reactor->GetName()<<endl;
				place++;
			}
		}
		if(PrintResult)cout<<endl;
	}
	if(option=="all"||option=="pool"){
		if(PrintResult)cout<<"Cooling Pools :"<<endl;
	
		Int_t nBranches =  T->GetNbranches();
		int place=0;
		for(Int_t i=0;i<nBranches;i++){
			string Name;
			Name=T->GetListOfBranches()->At(i)->GetName();
			Name=Name.substr(0,4);
			if(Name=="Pool"){
				Pool *pool;
				TString Rname;
				Rname.Form("Pool%d.",place);
		
				T->SetBranchStatus(Rname,1);
				T->SetBranchAddress(Rname, &pool);
	
				T->GetEntry(i);
				PoolIndex[pool->GetName()]=place;
				if(PrintResult)cout<<"\t "<<Rname<<"\t\t"<<pool->GetName()<<endl;
				place++;
			}
		}
		if(PrintResult)cout<<endl;
	}
	if(option=="all"||option=="storage"){
		if(PrintResult)cout<<"Storages :"<<endl;
	
		Int_t nBranches =  T->GetNbranches();
		int place=0;
		for(Int_t i=0;i<nBranches;i++){
			string Name;
			Name=T->GetListOfBranches()->At(i)->GetName();
			Name=Name.substr(0,7);
			if(Name=="Storage"){
				Storage *storage;
				TString Rname;
				Rname.Form("Storage%d.",place);
		
				T->SetBranchStatus(Rname,1);
				T->SetBranchAddress(Rname, &storage);
	
				T->GetEntry(i);
				StorageIndex[storage->GetName()]=place;
				if(PrintResult)cout<<"\t "<<Rname<<"\t"<<storage->GetName()<<endl;
				place++;
			}
		}
		if(PrintResult)cout<<endl;
	}
	if(option=="all"||option=="fabrication"){
		if(PrintResult)cout<<"Fabrication plants :"<<endl;
	
		Int_t nBranches =  T->GetNbranches();
		int place=0;
		for(Int_t i=0;i<nBranches;i++){
			string Name;
			Name=T->GetListOfBranches()->At(i)->GetName();
			Name=Name.substr(0,16);
			if(Name=="FabricationPlant"){
				FabricationPlant *fabrication;
				TString Rname;
				Rname.Form("FabricationPlant%d.",place);
		
				T->SetBranchStatus(Rname,1);
				T->SetBranchAddress(Rname, &fabrication);
	
				T->GetEntry(i);
				FabricationPlantIndex[fabrication->GetName()]=place;
				if(PrintResult)cout<<"\t "<<Rname<<"\t"<<fabrication->GetName()<<endl;
				place++;
			}
		}
		if(PrintResult)cout<<endl;
	}
	}
	T->ResetBranchAddresses();

};


void Read(TTree *T, TString IV, int Z, int A, int I, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	IsotopicVector *IVIn;
	TString IVname = IV +".";
	TString IVnameStatus = IV +".*";
	
	T->SetBranchStatus(IVname,1);
	T->SetBranchAddress(IVnameStatus, &IVIn);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	BranchName.Form("%d_%d_%d",Z,A,I);
	BranchName = IV + "_"+ BranchName;

	TString TitleName;
	TitleName.Form("%d %d %d",Z,A,I);
	TitleName = IV + " " + TitleName + " [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = IVIn->GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void Read(TTree *T, TString IV, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	IsotopicVector *IVIn;
	TString IVname = IV +".";
	TString IVnameStatus = IV +"*";

	T->SetBranchStatus(IVnameStatus,1);
	T->SetBranchAddress(IVname, &IVIn);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	BranchName = IV + "_Total";

	TString TitleName;
	TitleName = IV + " Total [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 1; i < nentries-1; i++)
	{

		T->GetEntry(i);
		ZAIQ = 0;
		for(int z = 96; z < 97; z++)
		{
			for(int a = 230; a < 260; a++)
			{
				for (int j=0; j < 2; j++)
				{
					ZAIQ += IVIn->GetZAIIsotopicQuantity(z,a,j)*a/6.02e23*1e-3;
				}
			}

		}

		cout << Time/365.25/24./3600. << " " << ZAIQ*1e-3 << endl;;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadReactor(TTree *T, string ReactorName, int Z, int A, int I, char* opt = "L*")
{
	if(ReactorIndex.size()==0)PrintNames(T,"reactor",false);
	int ID=ReactorIndex[ReactorName];
	ReadReactor(T,ID,Z,A,I,opt);
};

void ReadReactor(TTree *T, int ReactorId, int Z, int A, int I, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Reactor *reactor;

	string Rname = "Reactor" + itoa(ReactorId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &reactor);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	Name.Form("%d_%d_%d",Z,A,I);
	BranchName.Form("Reactor_%d",ReactorId);
	BranchName += "_"+ Name;

	TString TitleName;
	Name.Form("%d %d %d",Z,A,I);
	TitleName = "Reactor ";
	TitleName.Form("Reactor %d",ReactorId);
	TitleName += " " + Name + " [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = reactor->GetIVReactor().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadReactorList(TTree *T, int* liste, int Z, int A, int I)
{
	if( sizeof(liste) / sizeof(int) >0 )
	{
		ReadReactor(T, liste[0], Z, A, I, "L*");
		for (int i = 1; i < (int)sizeof(liste) / sizeof(int); i++)
		{
			ReadReactor(T, liste[i], Z, A, I, "same L*");
		}
	}
};

void ReadReactorList(TTree *T, int Z, int A, int I,char *opt = "L*", int ParamNumber, ...)
{
	va_list ap;
	int rId = -1;
	va_start(ap, ParamNumber);
	for (int i = 0; i < ParamNumber; i++) {
		if(i == 0)
		{

			rId = va_arg(ap, int);
			ReadReactor(T, rId, Z, A, I, opt);
		}
		else
		{
			rId = va_arg(ap, int);
			ReadReactor(T, rId, Z, A, I, "same L*");
		}
	}
	va_end(ap);
};

void ReadStorage(TTree *T, string StorageName, int Z, int A, int I, char* opt = "L*")
{
	if(StorageIndex.size()==0)PrintNames(T,"storage",false);
	int ID=StorageIndex[StorageName];
	ReadStorage(T,ID,Z,A,I,opt);
};

void ReadStorage(TTree *T, int StorageId, int Z, int A, int I, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Storage *storage;
	
	string Rname = "Storage" + itoa(StorageId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &storage);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	Name.Form("%d_%d_%d",Z,A,I);
	BranchName.Form("Storage_%d",StorageId);
	BranchName += "_"+ Name;

	TString TitleName;
	Name.Form("%d %d %d",Z,A,I);
	TitleName = "Storage ";
	TitleName.Form("Storage %d",StorageId);
	TitleName += " " + Name + " [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = storage->GetFullStock().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadStorage(TTree *T, string StorageName, char* opt = "L*")
{
	if(StorageIndex.size()==0)PrintNames(T,"storage",false);
	int ID=StorageIndex[StorageName];
	ReadStorage(T,ID,opt);
};

void ReadStorage(TTree *T, int StorageId, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Storage *storage;
	
	string Rname = "Storage" + itoa(StorageId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &storage);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	BranchName.Form("Storage_%d",StorageId);
	BranchName += "_total";

	TString TitleName;
	TitleName = "Storage ";
	TitleName.Form("Storage %d",StorageId);
	TitleName += " Total [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{
		T->GetEntry(i);
		ZAIQ = 0;
		//if(Time/365.25/24./3600.==(int)(Time/365.25/24./3600.))
		{
		for(int z = 10; z < 100; z++)
		{
			for(int a = 0; a < 260; a++)
			{
				for (int j=0; j < 2; j++)
				{
					ZAIQ += storage->GetFullStock().GetZAIIsotopicQuantity(z,a,j)*a/6.02e23*1e-3;
				}
			}
		}

		ZAIQ += storage->GetFullStock().GetZAIIsotopicQuantity(-2,-2,-2)*239/2./6.02e23*1e-3;

		cout << Time/365.25/24./3600. << " " << ZAIQ*1e-3 << endl;;
		}
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadCooling(TTree *T, string PoolName, int Z, int A, int I, char* opt = "L*")
{
	if(PoolIndex.size()==0)PrintNames(T,"pool",false);
	int ID=PoolIndex[PoolName];
	ReadCooling(T,ID,Z,A,I,opt);
};

void ReadCooling(TTree *T, int CoolingId, int Z, int A, int I, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Pool *cooling;
	
	string Rname = "Pool" + itoa(CoolingId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &cooling);
	
	
	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	Name.Form("%d_%d_%d",Z,A,I);
	BranchName.Form("Cooling_%d",CoolingId);
	BranchName += "_"+ Name;

	TString TitleName;
	Name.Form("%d %d %d",Z,A,I);
	TitleName = "Cooling ";
	TitleName.Form("Cooling %d",CoolingId);
	TitleName += " " + Name + " [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = cooling->GetFullCooling().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};


void ReadCooling(TTree *T, string PoolName, char* opt = "L*")
{
	if(PoolIndex.size()==0)PrintNames(T,"pool",false);
	int ID=PoolIndex[PoolName];
	ReadCooling(T,ID,opt);
};


void ReadCooling(TTree *T, int CoolingId, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Pool *cooling;
	
	string Rname = "Pool" + itoa(CoolingId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &cooling);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	BranchName.Form("Cooling_%d",CoolingId);
	BranchName += "_total";

	TString TitleName;
	TitleName = "Cooling ";
	TitleName.Form("Cooling %d",CoolingId);
	TitleName += " Total [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = 0;
		for(int z = 90; z < 98; z++)
		{
			for(int a = 230; a < 260; a++)
			{
				for (int j=0; j < 2; j++)
				{
					ZAIQ += cooling->GetFullCooling().GetZAIIsotopicQuantity(z,a,j)*a/6.02e23*1e-3;
				}
			}
		}


		ZAIQ += cooling->GetFullCooling().GetZAIIsotopicQuantity(-2,-2,-2)*239/2./6.02e23*1e-3;



			cout << Time/365.25/24./3600. << " " << ZAIQ*1e-3 << endl;;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadFabrication(TTree *T, string FabricationName, int Z, int A, int I, char* opt = "L*")
{
	if(FabricationIndex.size()==0)PrintNames(T,"fabrication",false);
	int ID=FabricationIndex[FabricationName];
	ReadFabrication(T,ID,Z,A,I,opt);
};

void ReadFabrication(TTree *T, int FabricationId, int Z, int A, int I, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	FabricationPlant *fabrication;
	
	string Rname = "FabricationPlant" + itoa(FabricationId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &fabrication);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	Name.Form("%d_%d_%d",Z,A,I);
	BranchName.Form("Fabrication_%d",FabricationId);
	BranchName += "_"+ Name;

	TString TitleName;
	Name.Form("%d %d %d",Z,A,I);
	TitleName = "Fabrication ";
	TitleName.Form("Fabrication %d",FabricationId);
	TitleName += " " + Name + " [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{
		T->GetEntry(i);
		ZAIQ = fabrication->GetFullFabrication().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void ReadFabrication(TTree *T, string FabricationName, char* opt = "L*")
{
	if(FabricationIndex.size()==0)PrintNames(T,"fabrication",false);
	int ID=FabricationIndex[FabricationName];
	ReadFabrication(T,ID,opt);
};


void ReadFabrication(TTree *T, int FabricationId, char* opt = "L*")
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	FabricationPlant *fabrication;
	
	string Rname = "FabricationPlant" + itoa(FabricationId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &fabrication);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	BranchName.Form("Fabrication_%d",FabricationId);
	BranchName += "_total";

	TString TitleName;
	TitleName = "Fabrication ";
	TitleName.Form("Fabrication %d",FabricationId);
	TitleName += " Total [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	for (Long64_t i = 0; i < nentries; i++)
	{

		T->GetEntry(i);
		ZAIQ = 0;
		for(int z = 90; z < 98; z++)
		{
			for(int a = 230; a < 260; a++)
			{
				for (int j=0; j < 2; j++)
				{
					ZAIQ += fabrication->GetFullFabrication().GetZAIIsotopicQuantity(z,a,j)*a/6.02e23*1e-3;
				}
			}
		}
		newBranch->Fill();
	}

	BranchName = BranchName + ":AbsTime/3600/24/365.25";
	T->Draw(BranchName,"",opt);
	TGraph *T_g = (TGraph*)gPad->GetPrimitive("Graph");
	TH2F   *htemp = (TH2F*)gPad->GetPrimitive("htemp"); // empty, but has axes

	htemp->GetXaxis()->SetTitle("Time [y]");
	htemp->GetYaxis()->SetTitle(TitleName);
	T_g->Draw(opt);

	T->ResetBranchAddresses();

};

void GetStockAt(TTree *T, int StorageId, int date)
{

	T->SetBranchStatus("*", 0);
	T->SetBranchStatus("AbsTime", 1);

	Storage *storage;
	
	string Rname = "Storage" + itoa(StorageId) + ".";
	string RnameStatus = Rname+"*";
	
	T->SetBranchStatus(RnameStatus.c_str(),1);
	T->SetBranchAddress(Rname.c_str(), &storage);

	Long64_t Time;
	T->SetBranchAddress("AbsTime", &Time);
	Long64_t nentries = T->GetEntries();
	Double_t ZAIQ;
	TString BranchName;

	TString Name;
	BranchName.Form("Storage_%d",StorageId);
	BranchName += "_total";

	TString TitleName;
	TitleName = "Storage ";
	TitleName.Form("Storage %d",StorageId);
	TitleName += " Total [kg]";


	TString Branchdescription;
	Branchdescription = BranchName +"/D";

  	TBranch *newBranch = T->Branch(BranchName, &ZAIQ, Branchdescription);
	T->SetTitle(TitleName);
	int CorrectEntrie = 0;
	long long int TimeDiff = abs(date*365.25*3600*24);
	for (Long64_t i = 1; i < nentries-1; i++)
	{
		T->GetEntry(i);
		if(TimeDiff > abs(date*365.25*3600*24 - Time) )
		{
			TimeDiff = abs(date*365.25*3600*24 - Time);
			CorrectEntrie = i;
		}

	}
	T->GetEntry(CorrectEntrie);
	string filename = "stock";
	storage->Write(filename.c_str(),date*365.25*3600*24.);


};
