#include "CLASSRead.hxx"

#include <TApplication.h>
#include <TGTableLayout.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <cmath>
#include <map>


#include "TFile.h"
#include "TBranch.h"
#include "TTree.h"
#include <TCanvas.h>
#include "TH2F.h"
#include "TGraph.h"
#include "TString.h"


using namespace std;
// for gcc3.2.3 only
char* operator+( std::streampos&, char* );
//end of gcc3.2.3



//________________________________________________________________________
string ReadNucleusName[] = {
	"n","H","He","Li","Be","B","C","N","O","F","Ne",			// List of nucleus used
	"Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
	"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
	"Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
	"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
	"Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
	"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
	"Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
	"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
	"Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm", //100
	"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","X",
	"-","-","-","-","-","-","-","-","-","-",				// add nucleus for DRAGON lattice code
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-",
	"-","-","-","-","-","-","-","-","-","-", //200
	"H_nat","He_nat","Li_nat","Be_nat","B_nat","C_nat","N_nat","O_nat","F_nat","Ne_nat",
	"Na_nat","Mg_nat","Al_nat","Si_nat","P_nat","S_nat","Cl_nat","Ar_nat","K_nat","Ca_nat",
	"Sc_nat","Ti_nat","V_nat","Cr_nat","Mn_nat","Fe_nat","Co_nat","Ni_nat","Cu_nat","Zn_nat",
	"Ga_nat","Ge_nat","As_nat","Se_nat","Br_nat","Kr_nat","Rb_nat","Sr_nat","Y_nat","Zr_nat",
	"Nb_nat","Mo_nat","Tc_nat","Ru_nat","Rh_nat","Pd_nat","Ag_nat","Cd_nat","In_nat","Sn_nat",
	"Sb_nat","Te_nat","I_nat","Xe_nat","Cs_nat","Ba_nat","La_nat","Ce_nat","Pr_nat","Nd_nat",
	"Pm_nat","Sm_nat","Eu_nat","Gd_nat","Tb_nat","Dy_nat","Ho_nat","Er_nat","Tm_nat","Yb_nat",
	"Lu_nat","Hf_nat","Ta_nat","W_nat","Re_nat","Os_nat","Ir_nat","Pt_nat","Au_nat","Hg_nat",
	"Tl_nat","Pb_nat","Bi_nat","Po_nat","At_nat","Rn_nat","Fr_nat","Ra_nat","Ac_nat","Th_nat",
	"Pa_nat","U_nat","Np_nat","Pu_nat","Am_nat","Cm_nat","Bk_nat","Cf_nat","Es_nat","Fm_nat"
};

int CurveColor(int graph_num)
{
	int ColorTable[]={1,kBlue,kRed,kGreen+1,kMagenta,kCyan+2,kOrange-3,kRed+2,kBlue-2,
		kSpring+9,kGreen+3,kAzure+8,kMagenta+2,kYellow+2,kBlue-9,kOrange+2};
	return ColorTable[graph_num%16];

}

string itoa(int num)
{
	ostringstream os(ostringstream::out);
	os<<setprecision(3)<<num;
	return os.str();
}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
CLASSRead::CLASSRead(TString filename)
{
	TFile *FileIn;
	FileIn = TFile::Open(filename);
	fFileIn.push_back(FileIn);

	for( int i =0; i < fFileIn.back()->GetNkeys(); i++)
	{
		//cout<<"KeyNum "<<i<<endl;
		fData.push_back( (TTree*)gDirectory->Get(fFileIn.back()->GetListOfKeys()->At(fFileIn.back()->GetNkeys()-1)->GetName() ) );
	}

	fCNucleiInv = 0;
	fGraphInv = 0;
	fLegendInv = 0;
	fNumberGraphInvIterator = 0;
	fGraphInvSumOfSelected = 0;
	fLegendInvSumOfSelected = 0;

	fCNucleiTox = 0;
	fGraphTox = 0;
	fLegendTox = 0;
	fNumberGraphToxIterator = 0;
	fGraphToxSumOfSelected = 0;
	fLegendToxSumOfSelected = 0;

	fCNucleiHeat = 0;
	fGraphHeat = 0;
	fLegendHeat = 0;
	fNumberGraphHeatIterator = 0;
	fGraphHeatSumOfSelected = 0;
	fLegendHeatSumOfSelected = 0;

	
	fCPower = 0;
	fGraphPower = 0;
	fLegendPower = 0;

	fNumberGraphPowerIterator = 0;
}

//________________________________________________________________________
CLASSRead::~CLASSRead()
{
	for(int i=fData.size()-1; i !=0; i--)
		delete fData[i];
	for(int i=fFileIn.size()-1; i !=0; i--)
		delete fFileIn[i];
}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
void CLASSRead::AddFile(TString filename)
{

	TFile *FileIn;
	FileIn = TFile::Open(filename);
	fFileIn.push_back(FileIn);
	for( int i =0; i < fFileIn.back()->GetNkeys(); i++)
	{
		fData.push_back( (TTree*)gDirectory->Get(fFileIn.back()->GetListOfKeys()->At(i)->GetName() ) );
	}
}

//________________________________________________________________________
void CLASSRead::ReadName()
{

	for (int j = 0; j < (int)fData.size() ; j++)
	{
		vector<TString>	ReactorName;
		vector<TString>	PoolName;
		vector<TString>	FabricationName;
		vector<TString>	StockName;

		Int_t nBranches = fData[j]->GetNbranches();
		int place_reactor=0;
		int place_pool=0;
		int place_fabrication=0;
		int place_stock=0;
		for(Int_t i=0;i<nBranches;i++)
		{
			string Name;
			Name=fData[j]->GetListOfBranches()->At(i)->GetName();
			string TypeName = Name;
			TypeName = TypeName.substr(0,2);

			if(TypeName=="R_")
			{
				ReactorName.push_back(Name.substr(2,Name.size()-1));

			}
			else if(TypeName=="P_")
			{
				PoolName.push_back(Name.substr(2,Name.size()-1));

			}
			else if(TypeName=="S_")
			{
				StockName.push_back(Name.substr(2,Name.size()-1));

			}
			else if(TypeName=="F_")
			{
				FabricationName.push_back(Name.substr(2,Name.size()-1));
				
			}
		}
		fReactorName.push_back(ReactorName);
		fPoolName.push_back(PoolName);
		fFabricationName.push_back(FabricationName);
		fStockName.push_back(StockName);
	}


}

//________________________________________________________________________
void CLASSRead::ReadZAI()
{

	IsotopicVector IVTot;
	for (int i = 0; i < (int)fData.size() ; i++)
	{
		IsotopicVector *IVIn=0;
		fData[i]->SetBranchStatus("TOTAL.",1);
		fData[i]->SetBranchAddress("TOTAL.", &IVIn);

		Long64_t nentries = fData[i]->GetEntries();

		for(int j = 0; j < nentries; j++)
		{
			fData[i]->GetEntry(j);
			IVTot += (*IVIn);
		}
		fData[i]->ResetBranchAddresses();
	}
	fZAIvector = IVTot.GetZAIList();
	fDecay = new IM_Matrix_Decay(IVTot);

}


//________________________________________________________________________
void CLASSRead::ReadTime()
{
	
	vector< vector<cSecond> > FullTimeVector;

	for (int i = 0; i < (int)fData.size() ; i++)
	{
		vector<cSecond> TimeVector;
		cSecond timeStep;
		fData[i]->SetBranchStatus("AbsTime",1);
		fData[i]->SetBranchAddress("AbsTime", &timeStep);
		
		Long64_t nentries = fData[i]->GetEntries();
		
		for(int j = 0; j < nentries; j++)
		{
			fData[i]->GetEntry(j);
			TimeVector.push_back(timeStep);
		}
		FullTimeVector.push_back(TimeVector);
		fData[i]->ResetBranchAddresses();
	}
	fTimeVector = FullTimeVector;
	
}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
void CLASSRead::PlotInv(vector<CLASSPlotElement> toplot, bool DecayChain, int StartingStep, cSecond FinalTime, int StepNUmber, bool LinBin, string opt)
{

	if(fGraphInv)
	{
		for(int i=0; i < fNumberGraphInvIterator;i++) delete fGraphInv[i];
		delete [] fGraphInv;
	}
	if(fLegendInvSumOfSelected)
		delete fLegendInvSumOfSelected;
	if(fGraphInvSumOfSelected)
		delete fGraphInvSumOfSelected;

	if(fLegendInv)
	{
		for(int i=0; i < fNumberGraphInvIterator;i++) delete fLegendInv[i];
		delete [] fLegendInv;
	}
	if(fCNucleiInv && gROOT->FindObject("c_NucleiInv"))
	{	delete fCNucleiInv;
		fCNucleiInv=0;
	}
	fCNucleiInv = new TCanvas("c_NucleiInv","NucleiInv",50,110,400,300);


	fGraphInv = new TGraph*[toplot.size()];
	fLegendInv = new TLatex*[toplot.size()];
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		fGraphInv[i] = 0;
		fLegendInv[i] = 0;
	}

	vector<CLASSPlotElement> toplotTTree[fData.size()];



	Xmin = +1.e36;
	Xmax =  -1.e36;
	Ymin = 1.e36;
	Ymax = -1.e36;

	fNumberGraphInvIterator = 0;
	bool SumOfSelected = false;

	if(toplot[0].fTreeId == -1 )
	{
		SumOfSelected = true;
		toplot.erase(toplot.begin());
	}

	for (int i = 0; i < (int)toplot.size(); i++)
	{
		toplotTTree[toplot[i].fTreeId].push_back(toplot[i]);
	}

	string out = opt;
	for (int i = 0; i < (int)fData.size(); i++)
	{
		if(i == 1) out += " same";
		if(toplotTTree[i].size() !=0)
			if(!DecayChain)
				BuildTGraph(toplotTTree[i], 0, out);

	}
	fCNucleiInv->cd();

	double X_Sum[fGraphInv[0]->GetN()];
	double Y_Sum[fGraphInv[0]->GetN()];

	if(SumOfSelected)
	{
		for (int i = 0; i < (int)fNumberGraphInvIterator; i++)
		{

			for(int j = 0; j < fGraphInv[i]->GetN(); j++)
			{
				double x;
				double y;
				fGraphInv[i]->GetPoint(j, x, y);
				if(i == 0)
					X_Sum[j] = x;
				if(i == 0)
					Y_Sum[j] = y;
				else
					Y_Sum[j] += y;
			}
		}


		for (int i =0; i < fGraphInv[0]->GetN(); i++)
		{
			if(X_Sum[i] > Xmax) Xmax = X_Sum[i];
			if(X_Sum[i] < Xmin) Xmin = X_Sum[i];
			if(Y_Sum[i] > Ymax) Ymax = Y_Sum[i];
			if(Y_Sum[i] < Ymin) Ymin = Y_Sum[i];
		}
	}



	TH1F*	  fhr = fCNucleiInv->DrawFrame(Xmin,Ymin*0.95,Xmax,Ymax*1.05);
	string Xtitle="Time [year]";
	string Ytitle="Mass [kg]";
	fhr->SetXTitle(Xtitle.c_str());
	fhr->SetYTitle(Ytitle.c_str());
	fhr->GetXaxis()->CenterTitle();
	fhr->GetYaxis()->CenterTitle();
	fhr->GetYaxis()->SetTitleOffset(1.25);


	if(SumOfSelected)
	{
		fGraphInvSumOfSelected = new TGraph(fGraphInv[0]->GetN(), X_Sum, Y_Sum );
		fGraphInvSumOfSelected->SetName("Sum_Of_Selected");
		fGraphInvSumOfSelected->SetTitle("Sum Of Selected");
		fGraphInvSumOfSelected->SetLineColor(CurveColor(fNumberGraphInvIterator));
		fGraphInvSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphInvIterator));
		fGraphInvSumOfSelected->SetMarkerStyle(10);
		fGraphInvSumOfSelected->Draw(out.c_str());
		fGraphInvSumOfSelected->SetLineColor(CurveColor(fNumberGraphInvIterator));
		fGraphInvSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphInvIterator));

		double x;
		double y;
		double x_0;
		double y_0;
		fGraphInvSumOfSelected->GetPoint(1, x_0, y_0);
		fGraphInvSumOfSelected->GetPoint(fGraphInvSumOfSelected->GetN()-1, x, y);

		fLegendInvSumOfSelected = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),"Sum_Of_Selected");
		fLegendInvSumOfSelected->SetTextSize(0.05);
		fLegendInvSumOfSelected->SetTextFont(132);
		fLegendInvSumOfSelected->SetTextColor(CurveColor(fNumberGraphInvIterator));
		fLegendInvSumOfSelected->Draw();


	}


	for (int i = 0; i < (int)fNumberGraphInvIterator; i++)
	{
		if( i !=0 || SumOfSelected) out += " same";

		fGraphInv[i]->SetName(GetTittleOutName(toplot[i]).c_str());
		fGraphInv[i]->SetTitle(GetTittleOutName(toplot[i]).c_str());
		fGraphInv[i]->SetLineColor(CurveColor(i));
		fGraphInv[i]->SetMarkerColor(CurveColor(i));
		fGraphInv[i]->SetMarkerStyle(10);
		fGraphInv[i]->Draw(out.c_str());
		fGraphInv[i]->SetLineColor(CurveColor(i));
		fGraphInv[i]->SetMarkerColor(CurveColor(i));

		double x;
		double y;
		double x_0;
		double y_0;
		fGraphInv[i]->GetPoint(1, x_0, y_0);
		
		fGraphInv[i]->GetPoint(fGraphInv[i]->GetN()-1, x, y);

		fLegendInv[i] = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),GetLegendOutName(toplot[i]).c_str());
		fLegendInv[i]->SetTextSize(0.05);
		fLegendInv[i]->SetTextFont(132);
		fLegendInv[i]->SetTextColor(CurveColor(i));
		fLegendInv[i]->Draw();
	}

	fCNucleiInv->Update();




}

//________________________________________________________________________
void CLASSRead::PlotTox(vector<CLASSPlotElement> toplot, bool DecayChain, int StartingStep, cSecond FinalTime, int StepNUmber, bool LinBin, string opt)
{
	if(fGraphTox)
	{
		for(int i=0; i < fNumberGraphToxIterator;i++) delete fGraphTox[i];
		delete [] fGraphTox;
	}
	if(fLegendToxSumOfSelected)
		delete fLegendToxSumOfSelected;
	if(fGraphToxSumOfSelected)
		delete fGraphToxSumOfSelected;
	
	if(fLegendTox)
	{
		for(int i=0; i < fNumberGraphToxIterator;i++) delete fLegendTox[i];
		delete [] fLegendTox;
	}
	if(fCNucleiTox && gROOT->FindObject("c_NucleiTox"))
	{	delete fCNucleiTox;
		fCNucleiTox=0;
	}
	fCNucleiTox = new TCanvas("c_NucleiTox","NucleiTox",50,110,400,300);
	
	
	fGraphTox = new TGraph*[toplot.size()];
	fLegendTox = new TLatex*[toplot.size()];
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		fGraphTox[i] = 0;
		fLegendTox[i] = 0;
	}
	
	vector<CLASSPlotElement> toplotTTree[fData.size()];

	
	Xmin = +1.e36;
	Xmax =  -1.e36;
	Ymin = 1.e36;
	Ymax = -1.e36;
	
	fNumberGraphToxIterator = 0;
	bool SumOfSelected = false;
	
	if(toplot[0].fTreeId == -1 )
	{
		SumOfSelected = true;
		toplot.erase(toplot.begin());
	}
	
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		toplotTTree[toplot[i].fTreeId].push_back(toplot[i]);
	}
	
	string out = opt;
	for (int i = 0; i < (int)fData.size(); i++)
	{
		if(i == 1) out += " same";
		if(toplotTTree[i].size() !=0)
			if(!DecayChain)
				BuildTGraph(toplotTTree[i], 1, out);
		
	}
	fCNucleiTox->cd();
	
	double X_Sum[fGraphTox[0]->GetN()];
	double Y_Sum[fGraphTox[0]->GetN()];
	
	if(SumOfSelected)
	{
		for (int i = 0; i < (int)fNumberGraphToxIterator; i++)
		{
			
			for(int j = 0; j < fGraphTox[i]->GetN(); j++)
			{
				double x;
				double y;
				fGraphTox[i]->GetPoint(j, x, y);
				if(i == 0)
					X_Sum[j] = x;
				if(i == 0)
					Y_Sum[j] = y;
				else
					Y_Sum[j] += y;
			}
		}
		
		
		for (int i =0; i < fGraphTox[0]->GetN(); i++)
		{
			if(X_Sum[i] > Xmax) Xmax = X_Sum[i];
			if(X_Sum[i] < Xmin) Xmin = X_Sum[i];
			if(Y_Sum[i] > Ymax) Ymax = Y_Sum[i];
			if(Y_Sum[i] < Ymin) Ymin = Y_Sum[i];
		}
	}
	
	
	
	TH1F*	  fhr = fCNucleiTox->DrawFrame(Xmin,Ymin*0.95,Xmax,Ymax*1.05);
	string Xtitle="Time [year]";
	string Ytitle="Radio-Toxicity [Sv]";
	fhr->SetXTitle(Xtitle.c_str());
	fhr->SetYTitle(Ytitle.c_str());
	fhr->GetXaxis()->CenterTitle();
	fhr->GetYaxis()->CenterTitle();
	fhr->GetYaxis()->SetTitleOffset(1.25);
	
	
	if(SumOfSelected)
	{
		fGraphToxSumOfSelected = new TGraph(fGraphTox[0]->GetN(), X_Sum, Y_Sum );
		fGraphToxSumOfSelected->SetName("Sum_Of_Selected");
		fGraphToxSumOfSelected->SetTitle("Sum Of Selected");
		fGraphToxSumOfSelected->SetLineColor(CurveColor(fNumberGraphToxIterator));
		fGraphToxSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphToxIterator));
		fGraphToxSumOfSelected->SetMarkerStyle(10);
		fGraphToxSumOfSelected->Draw(out.c_str());
		fGraphToxSumOfSelected->SetLineColor(CurveColor(fNumberGraphToxIterator));
		fGraphToxSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphToxIterator));
		
		double x;
		double y;
		double x_0;
		double y_0;
		fGraphToxSumOfSelected->GetPoint(1, x_0, y_0);
		fGraphToxSumOfSelected->GetPoint(fGraphToxSumOfSelected->GetN()-1, x, y);
		
		fLegendToxSumOfSelected = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),"Sum_Of_Selected");
		fLegendToxSumOfSelected->SetTextSize(0.05);
		fLegendToxSumOfSelected->SetTextFont(132);
		fLegendToxSumOfSelected->SetTextColor(CurveColor(fNumberGraphToxIterator));
		fLegendToxSumOfSelected->Draw();
		
		
	}
	
	
	for (int i = 0; i < (int)fNumberGraphToxIterator; i++)
	{
		if( i !=0 || SumOfSelected) out += " same";
		
		fGraphTox[i]->SetName(GetTittleOutName(toplot[i]).c_str());
		fGraphTox[i]->SetTitle(GetTittleOutName(toplot[i]).c_str());
		fGraphTox[i]->SetLineColor(CurveColor(i));
		fGraphTox[i]->SetMarkerColor(CurveColor(i));
		fGraphTox[i]->SetMarkerStyle(10);
		fGraphTox[i]->Draw(out.c_str());
		fGraphTox[i]->SetLineColor(CurveColor(i));
		fGraphTox[i]->SetMarkerColor(CurveColor(i));
		
		double x;
		double y;
		double x_0;
		double y_0;
		fGraphTox[i]->GetPoint(1, x_0, y_0);
		
		fGraphTox[i]->GetPoint(fGraphTox[i]->GetN()-1, x, y);
		
		fLegendTox[i] = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),GetLegendOutName(toplot[i]).c_str());
		fLegendTox[i]->SetTextSize(0.05);
		fLegendTox[i]->SetTextFont(132);
		fLegendTox[i]->SetTextColor(CurveColor(i));
		fLegendTox[i]->Draw();
	}
	
	fCNucleiTox->Update();
	
	
	
	
}

//________________________________________________________________________
void CLASSRead::PlotHeat(vector<CLASSPlotElement> toplot, bool DecayChain, int StartingStep, cSecond FinalTime, int StepNUmber, bool LinBin, string opt)
{
	
	if(fGraphHeat)
	{
		for(int i=0; i < fNumberGraphHeatIterator;i++) delete fGraphHeat[i];
		delete [] fGraphHeat;
	}
	if(fLegendHeatSumOfSelected)
		delete fLegendHeatSumOfSelected;
	if(fGraphHeatSumOfSelected)
		delete fGraphHeatSumOfSelected;
	
	if(fLegendHeat)
	{
		for(int i=0; i < fNumberGraphHeatIterator;i++) delete fLegendHeat[i];
		delete [] fLegendHeat;
	}
	if(fCNucleiHeat && gROOT->FindObject("c_NucleiHeat"))
	{	delete fCNucleiHeat;
		fCNucleiHeat=0;
	}
	fCNucleiHeat = new TCanvas("c_NucleiHeat","NucleiHeat",50,110,400,300);
	
	
	fGraphHeat = new TGraph*[toplot.size()];
	fLegendHeat = new TLatex*[toplot.size()];
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		fGraphHeat[i] = 0;
		fLegendHeat[i] = 0;
	}
	
	vector<CLASSPlotElement> toplotTTree[fData.size()];
	
	
	
	Xmin = +1.e36;
	Xmax =  -1.e36;
	Ymin = 1.e36;
	Ymax = -1.e36;
	
	fNumberGraphHeatIterator = 0;
	bool SumOfSelected = false;
	
	if(toplot[0].fTreeId == -1 )
	{
		SumOfSelected = true;
		toplot.erase(toplot.begin());
	}
	
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		toplotTTree[toplot[i].fTreeId].push_back(toplot[i]);
	}
	
	string out = opt;
	for (int i = 0; i < (int)fData.size(); i++)
	{
		if(i == 1) out += " same";
		if(toplotTTree[i].size() !=0)
			if(!DecayChain)
				BuildTGraph(toplotTTree[i], 2, out);
		
	}
	fCNucleiHeat->cd();
	
	double X_Sum[fGraphHeat[0]->GetN()];
	double Y_Sum[fGraphHeat[0]->GetN()];
	
	if(SumOfSelected)
	{
		for (int i = 0; i < (int)fNumberGraphHeatIterator; i++)
		{
			
			for(int j = 0; j < fGraphHeat[i]->GetN(); j++)
			{
				double x;
				double y;
				fGraphHeat[i]->GetPoint(j, x, y);
				if(i == 0)
					X_Sum[j] = x;
				if(i == 0)
					Y_Sum[j] = y;
				else
					Y_Sum[j] += y;
			}
		}
		
		
		for (int i =0; i < fGraphHeat[0]->GetN(); i++)
		{
			if(X_Sum[i] > Xmax) Xmax = X_Sum[i];
			if(X_Sum[i] < Xmin) Xmin = X_Sum[i];
			if(Y_Sum[i] > Ymax) Ymax = Y_Sum[i];
			if(Y_Sum[i] < Ymin) Ymin = Y_Sum[i];
		}
	}
	
	
	
	TH1F*	  fhr = fCNucleiHeat->DrawFrame(Xmin,Ymin*0.95,Xmax,Ymax*1.05);
	string Xtitle="Time [year]";
	string Ytitle="Decay Heat [W]";
	fhr->SetXTitle(Xtitle.c_str());
	fhr->SetYTitle(Ytitle.c_str());
	fhr->GetXaxis()->CenterTitle();
	fhr->GetYaxis()->CenterTitle();
	fhr->GetYaxis()->SetTitleOffset(1.25);
	
	
	if(SumOfSelected)
	{
		fGraphHeatSumOfSelected = new TGraph(fGraphHeat[0]->GetN(), X_Sum, Y_Sum );
		fGraphHeatSumOfSelected->SetName("Sum_Of_Selected");
		fGraphHeatSumOfSelected->SetTitle("Sum Of Selected");
		fGraphHeatSumOfSelected->SetLineColor(CurveColor(fNumberGraphHeatIterator));
		fGraphHeatSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphHeatIterator));
		fGraphHeatSumOfSelected->SetMarkerStyle(10);
		fGraphHeatSumOfSelected->Draw(out.c_str());
		fGraphHeatSumOfSelected->SetLineColor(CurveColor(fNumberGraphHeatIterator));
		fGraphHeatSumOfSelected->SetMarkerColor(CurveColor(fNumberGraphHeatIterator));
		
		double x;
		double y;
		double x_0;
		double y_0;
		fGraphHeatSumOfSelected->GetPoint(1, x_0, y_0);
		fGraphHeatSumOfSelected->GetPoint(fGraphHeatSumOfSelected->GetN()-1, x, y);
		
		fLegendHeatSumOfSelected = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),"Sum_Of_Selected");
		fLegendHeatSumOfSelected->SetTextSize(0.05);
		fLegendHeatSumOfSelected->SetTextFont(132);
		fLegendHeatSumOfSelected->SetTextColor(CurveColor(fNumberGraphHeatIterator));
		fLegendHeatSumOfSelected->Draw();
		
		
	}
	
	
	for (int i = 0; i < (int)fNumberGraphHeatIterator; i++)
	{
		if( i !=0 || SumOfSelected) out += " same";
		
		fGraphHeat[i]->SetName(GetTittleOutName(toplot[i]).c_str());
		fGraphHeat[i]->SetTitle(GetTittleOutName(toplot[i]).c_str());
		fGraphHeat[i]->SetLineColor(CurveColor(i));
		fGraphHeat[i]->SetMarkerColor(CurveColor(i));
		fGraphHeat[i]->SetMarkerStyle(10);
		fGraphHeat[i]->Draw(out.c_str());
		fGraphHeat[i]->SetLineColor(CurveColor(i));
		fGraphHeat[i]->SetMarkerColor(CurveColor(i));
		
		double x;
		double y;
		double x_0;
		double y_0;
		fGraphHeat[i]->GetPoint(1, x_0, y_0);
		
		fGraphHeat[i]->GetPoint(fGraphHeat[i]->GetN()-1, x, y);
		
		fLegendHeat[i] = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),GetLegendOutName(toplot[i]).c_str());
		fLegendHeat[i]->SetTextSize(0.05);
		fLegendHeat[i]->SetTextFont(132);
		fLegendHeat[i]->SetTextColor(CurveColor(i));
		fLegendHeat[i]->Draw();
	}
	
	fCNucleiHeat->Update();
	
	
	
	
}


//________________________________________________________________________
void CLASSRead::BuildTGraph(vector<CLASSPlotElement> toplot, int PlotId, string opt)
{

	TGraph** Graph;
	int NumberGraphIterator = 0;

	if(PlotId == 0)
		Graph = fGraphInv;
	else if(PlotId == 1)
		Graph = fGraphTox;
	else if(PlotId == 2)
		Graph = fGraphHeat;
	else
	{
		cout << "Bad PlotId" << endl;
		return;
	}
	
	fData[toplot[0].fTreeId]->SetBranchStatus("*", 0);
	fData[toplot[0].fTreeId]->SetBranchStatus("AbsTime", 1);



	string out = opt;
	Long64_t nentries = fData[toplot[0].fTreeId]->GetEntries();

	Long64_t Time = 0;
	fData[toplot[0].fTreeId]->SetBranchAddress("AbsTime", &Time);

	Reactor* reactor[fReactorName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fReactorName[toplot[0].fTreeId].size(); i++ )
		reactor[i] = 0;

	Pool* pool[fPoolName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fPoolName[toplot[0].fTreeId].size(); i++ )
		pool[i] = 0;

	FabricationPlant* fabricationplant[fFabricationName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fFabricationName[toplot[0].fTreeId].size(); i++ )
		fabricationplant[i] = 0;

	Storage* stock[fStockName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fStockName[toplot[0].fTreeId].size(); i++ )
		stock[i] = 0;

	IsotopicVector* IV[8];
	for(int i = 0; i < 8; i++ )
		IV[i] = 0;

	vector< double > vTime;
	vector< double > vQuantity[toplot.size()];

	for (int i=0; i < (int)toplot.size(); i++)
	{

		string InBranchName = GetBranchInName(toplot[i]);

		string ActiveInBranchName = InBranchName + "*";
		fData[toplot[0].fTreeId]->SetBranchStatus(ActiveInBranchName.c_str(),1);

		if(toplot[i].fFacilityId == 0)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &IV[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 1)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &reactor[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 2)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &stock[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 3)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &pool[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 4)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &fabricationplant[toplot[i].fFacylityNumber]);
	}


	for (Long64_t  j = 0; j < nentries; j++)
	{
		fData[toplot[0].fTreeId]->GetEntry(j);

		vTime.push_back(Time/3600./24./365.25);

		if(Xmin>vTime.back()) Xmin = vTime.back();
		if(Xmax<vTime.back()) Xmax = vTime.back();

		for (int i=0; i < (int)toplot.size(); i++)
		{

			if(toplot[i].fFacilityId == 0)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();
				double ZAIQuantity = IV[toplot[i].fFacylityNumber]->GetZAIIsotopicQuantity(Z,A,I);
			
				if(PlotId == 0)
					ZAIQuantity *= cZAIMass.GetMass(Z,A)/AVOGADRO*1e-3;
				else if(PlotId == 1)
					ZAIQuantity *= cZAITox.GetRadioTox(Z,A,I);
				else if(PlotId == 2)
					ZAIQuantity *= cZAIHeat.GetHeat(Z,A,I);
				else
				{
					cout << "Bad PlotId" << endl;
					return;
				}

				
				vQuantity[i].push_back(ZAIQuantity);



				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 1)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I);
				else
				{
					cout << "Bad IVNumber" << endl;
					return;
				}

				if(PlotId == 0)
					ZAIQuantity *= cZAIMass.GetMass(Z,A)/AVOGADRO*1e-3;
				else if(PlotId == 1)
					ZAIQuantity *= cZAITox.GetRadioTox(Z,A,I);
				else if(PlotId == 2)
					ZAIQuantity *= cZAIHeat.GetHeat(Z,A,I);
				else
				{
					cout << "Bad PlotId" << endl;
					return;
				}

				vQuantity[i].push_back(ZAIQuantity);

				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 2)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I);
				else
				{
					cout << "Bad IVNumber" << endl;
					return;
				}
				
				if(PlotId == 0)
					ZAIQuantity *= cZAIMass.GetMass(Z,A)/AVOGADRO*1e-3;
				else if(PlotId == 1)
					ZAIQuantity *= cZAITox.GetRadioTox(Z,A,I);
				else if(PlotId == 2)
					ZAIQuantity *= cZAIHeat.GetHeat(Z,A,I);
				else
				{
					cout << "Bad PlotId" << endl;
					return;
				}

				vQuantity[i].push_back(ZAIQuantity);


				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 3)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I);
				else
				{
					cout << "Bad IVNumber" << endl;
					return;
				}
				
				if(PlotId == 0)
					ZAIQuantity *= cZAIMass.GetMass(Z,A)/AVOGADRO*1e-3;
				else if(PlotId == 1)
					ZAIQuantity *= cZAITox.GetRadioTox(Z,A,I);
				else if(PlotId == 2)
					ZAIQuantity *= cZAIHeat.GetHeat(Z,A,I);
				else
				{
					cout << "Bad PlotId" << endl;
					return;
				}

				vQuantity[i].push_back(ZAIQuantity);


				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 4)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I);
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I);
				else
				{
					cout << "Bad IVNumber" << endl;
					return;
				}
				
				if(PlotId == 0)
					ZAIQuantity *= cZAIMass.GetMass(Z,A)/AVOGADRO*1e-3;
				else if(PlotId == 1)
					ZAIQuantity *= cZAITox.GetRadioTox(Z,A,I);
				else if(PlotId == 2)
					ZAIQuantity *= cZAIHeat.GetHeat(Z,A,I);
				else
				{
					cout << "Bad PlotId" << endl;
					return;
				}

				vQuantity[i].push_back(ZAIQuantity);

				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
		}


	}



	for (int i = 0; i < (int)toplot.size(); i++)
	{
		Graph[NumberGraphIterator] = new TGraph(vTime.size(), &vTime[0], &(vQuantity[i])[0]);
		NumberGraphIterator++;
	}

	fData[toplot[0].fTreeId]->ResetBranchAddresses();
	{
		for(int i=0; i< (int)fReactorName[toplot[0].fTreeId].size(); i++) delete reactor[i];

		for(int i=0; i< (int)fPoolName[toplot[0].fTreeId].size(); i++) delete pool[i];

		for(int i=0; i< (int)fFabricationName[toplot[0].fTreeId].size(); i++) delete fabricationplant[i];

		for(int i=0; i< (int)fStockName[toplot[0].fTreeId].size(); i++) delete stock[i];
		for(int i=0; i< 8; i++) delete IV[i];
	}

	
	if(PlotId == 0)
		fNumberGraphInvIterator += NumberGraphIterator;
	else if(PlotId == 1)
		fNumberGraphToxIterator += NumberGraphIterator;
	else if(PlotId == 2)
		fNumberGraphHeatIterator += NumberGraphIterator;
	else
	{
		cout << "Bad PlotId" << endl;
		return;
	}
	

}


//________________________________________________________________________
//________________________________________________________________________
void CLASSRead::PlotPower(vector<CLASSPlotElement> toplot, string opt)
{
	if(fGraphPower)
	{
		for(int i=0; i < fNumberGraphPowerIterator;i++) delete fGraphPower[i];
		delete [] fGraphPower;
	}


	if(fLegendPower)
	{
		for(int i=0; i < fNumberGraphPowerIterator;i++) delete fLegendPower[i];
		delete [] fLegendPower;
	}

	if(fCPower && gROOT->FindObject("fCPower"))
	{	delete fCPower;
		fCPower=0;
	}


	fCPower = new TCanvas("fCPower","Power",50,110,400,300);
	fGraphPower = new TGraph*[fData.size()];
	fLegendPower = new TLatex*[fData.size()];


	for (int i = 0; i < (int)fData.size(); i++)
	{
		fGraphPower[i] = 0;
		fLegendPower[i] = 0;
	}

	vector<CLASSPlotElement> toplotTTree[fData.size()];

	fNumberGraphPowerIterator = 0;
	Xmin = +1.e36;
	Xmax =  -1.e36;
	Ymin = 1.e36;
	Ymax = -1.e36;




	for (int i = 0; i < (int)toplot.size(); i++)
	{
		toplotTTree[toplot[i].fTreeId].push_back(toplot[i]);
	}

	string out = opt;
	for (int i = 0; i < (int)fData.size(); i++)
	{

		if(i != 0) out += " same";
		if(toplotTTree[i].size() !=0)
			PlotTTreePower(toplotTTree[i], out);


	}

	fCPower->cd();

	TH1F*	  fhrPower= fCPower->DrawFrame(Xmin,Ymin*0.95,Xmax,Ymax*1.05);
	string Xtitle="Time [year]";
	string Ytitle="Total Thermal Power [GW]";
	fhrPower->SetXTitle(Xtitle.c_str());
	fhrPower->SetYTitle(Ytitle.c_str());
	fhrPower->GetXaxis()->CenterTitle();
	fhrPower->GetYaxis()->CenterTitle();
	fhrPower->GetYaxis()->SetTitleOffset(1.25);

	for (int i = 0; i < (int)fNumberGraphPowerIterator; i++)
	{
		fCPower->cd();
		if( i !=0 ) out += " same";

		fGraphPower[i]->SetName(GetTittleOutName(toplot[i]).c_str());
		fGraphPower[i]->SetTitle(GetTittleOutName(toplot[i]).c_str());
		fGraphPower[i]->SetLineColor(CurveColor(i));
		fGraphPower[i]->SetMarkerColor(CurveColor(i));
		fGraphPower[i]->SetMarkerStyle(10);
		fGraphPower[i]->Draw(out.c_str());

		fCPower->cd();

		double x;
		double y;
		double x_0;
		double y_0;
		fGraphPower[i]->GetPoint(1, x_0, y_0);
		
		fGraphPower[i]->GetPoint(fGraphPower[i]->GetN()-1, x, y);
				
		fLegendPower[i] = new TLatex(x_0+0.6*(x-x_0),y_0+1.05*(y-y_0),GetLegendOutName(toplot[i]).c_str());
		fLegendPower[i]->SetTextSize(0.05);
		fLegendPower[i]->SetTextFont(132);
		fLegendPower[i]->SetTextColor(CurveColor(i));
		fLegendPower[i]->Draw();
	}
	fCPower->Update();
	
	
	
}

//________________________________________________________________________
void CLASSRead::PlotTTreePower(vector<CLASSPlotElement> toplot, string opt)
{

	fData[toplot[0].fTreeId]->SetBranchStatus("*", 0);
	fData[toplot[0].fTreeId]->SetBranchStatus("AbsTime", 1);
	fData[toplot[0].fTreeId]->SetBranchStatus("ParcPower", 1);


	string out = opt;
	Long64_t nentries = fData[toplot[0].fTreeId]->GetEntries();


	Long64_t Time = 0;
	fData[toplot[0].fTreeId]->SetBranchAddress("AbsTime", &Time);
	
	double Power = 0;
	fData[toplot[0].fTreeId]->SetBranchAddress("ParcPower", &Power);


	double  vTime[nentries];
	double  vPower[nentries];



	for (Long64_t  j = 0; j < nentries; j++)
	{
		fData[toplot[0].fTreeId]->GetEntry(j);

		vTime[j] = Time/3600./24./365.25;
		vPower[j] = Power/1.e9;

		if(Xmin>vTime[j]) Xmin = vTime[j];
		if(Xmax<vTime[j]) Xmax = vTime[j];

		if(Ymin>vPower[j]) Ymin = vPower[j];
		if(Ymax<vPower[j]) Ymax = vPower[j];
	}


	fGraphPower[fNumberGraphPowerIterator] = new TGraph(nentries, vTime, vPower);

	fData[toplot[0].fTreeId]->ResetBranchAddresses();


	fNumberGraphPowerIterator++;

}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
void CLASSRead::Write(string filename, string fileformat)
{
	if(fileformat == "ASCII")
		ASCIIWrite(filename);

}
void CLASSRead::ASCIIWrite(string filename)
{

	ofstream outfile;
	outfile.open(filename.c_str());
	if(!outfile)
	{
		cout << "Could not open : " << filename << " !" << endl;
		exit(-1);
	}

	cout << "WARNING!! not working if using many CLASS.root file with diffenret timestep!!!"<<endl;

	if (fGraphInv)
	{
		double* X = fGraphInv[0]->GetX();

		outfile << "time";
		for(int i= 0; i < fGraphInv[0]->GetN(); i++)
			outfile << "\t" << X[i];

		outfile << endl;
	}


	if (fGraphInvSumOfSelected) {

		outfile << fGraphInvSumOfSelected->GetTitle();
		double* Y = fGraphInvSumOfSelected->GetY();
		for(int j= 0; j < fGraphInvSumOfSelected->GetN(); j++)
			outfile << "\t" << Y[j];
		outfile << endl;
	}


	for(int i = 0; i < fNumberGraphInvIterator; i++)
	{
		outfile << fGraphInv[i]->GetTitle();
		double* Y = fGraphInv[i]->GetY();
		for(int j= 0; j < fGraphInv[i]->GetN(); j++)
			outfile << "\t" << Y[j];
		outfile << endl;
	}


}



//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
void CLASSRead::ConvertxmlTTreeMass(vector<CLASSPlotElement> toplot, string filename)
{



	fData[toplot[0].fTreeId]->SetBranchStatus("*", 0);
	fData[toplot[0].fTreeId]->SetBranchStatus("AbsTime", 1);



	Long64_t nentries = fData[toplot[0].fTreeId]->GetEntries();


	Long64_t Time = 0;
	fData[toplot[0].fTreeId]->SetBranchAddress("AbsTime", &Time);

	Reactor* reactor[fReactorName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fReactorName[toplot[0].fTreeId].size(); i++ )
		reactor[i] = 0;

	Pool* pool[fPoolName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fPoolName[toplot[0].fTreeId].size(); i++ )
		pool[i] = 0;

	FabricationPlant* fabricationplant[fFabricationName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fFabricationName[toplot[0].fTreeId].size(); i++ )
		fabricationplant[i] = 0;

	Storage* stock[fStockName[toplot[0].fTreeId].size()];
	for(int i = 0; i < (int)fStockName[toplot[0].fTreeId].size(); i++ )
		stock[i] = 0;

	IsotopicVector* IV[8];
	for(int i = 0; i < 8; i++ )
		IV[i] = 0;

	vector< double > vTime;
	vector< double > vQuantity[toplot.size()];

	for (int i=0; i < (int)toplot.size(); i++)
	{

		string InBranchName = GetBranchInName(toplot[i]);

		string ActiveInBranchName = InBranchName + "*";
		fData[toplot[0].fTreeId]->SetBranchStatus(ActiveInBranchName.c_str(),1);

		if(toplot[i].fFacilityId == 0)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &IV[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 1)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &reactor[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 2)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &stock[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 3)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &pool[toplot[i].fFacylityNumber]);
		else if(toplot[i].fFacilityId == 4)
			fData[toplot[i].fTreeId]->SetBranchAddress(InBranchName.c_str(), &fabricationplant[toplot[i].fFacylityNumber]);
	}


	for (Long64_t  j = 0; j < nentries; j++)
	{
		fData[toplot[0].fTreeId]->GetEntry(j);

		vTime.push_back(Time/3600./24./365.25);

		if(Xmin>vTime.back()) Xmin = vTime.back();
		if(Xmax<vTime.back()) Xmax = vTime.back();

		for (int i=0; i < (int)toplot.size(); i++)
		{

			if(toplot[i].fFacilityId == 0)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();
				double ZAIQuantity = IV[toplot[i].fFacylityNumber]->GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;

				vQuantity[i].push_back(ZAIQuantity);



				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 1)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = reactor[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else
				{
					cout << "Bad IVNumber" << endl;
					break;
				}

				vQuantity[i].push_back(ZAIQuantity);

				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 2)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = stock[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else
				{
					cout << "Bad IVNumber" << endl;
					break;
				}

				vQuantity[i].push_back(ZAIQuantity);


				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 3)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = pool[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else
				{
					cout << "Bad IVNumber" << endl;
					break;
				}

				vQuantity[i].push_back(ZAIQuantity);


				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
			else if(toplot[i].fFacilityId == 4)
			{
				int Z = toplot[i].fZAI.Z();
				int A = toplot[i].fZAI.A();
				int I = toplot[i].fZAI.I();

				double ZAIQuantity = 0;

				if( toplot[i].fIVNumber == 0 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetInsideIV().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 1 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetCumulativeIVIn().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else if( toplot[i].fIVNumber == 2 )
					ZAIQuantity = fabricationplant[toplot[i].fFacylityNumber]->GetCumulativeIVOut().GetZAIIsotopicQuantity(Z,A,I)*A/6.02e23*1e-3;
				else
				{
					cout << "Bad IVNumber" << endl;
					break;
				}

				vQuantity[i].push_back(ZAIQuantity);

				if(Ymin>ZAIQuantity) Ymin = ZAIQuantity;
				if(Ymax<ZAIQuantity) Ymax = ZAIQuantity;

			}
		}


	}







	// Beginning of the document XML
	ofstream f (filename.c_str());
	cout << f.is_open();
	if (!f.is_open())
		cout << "Impossible d'ouvrir le fichier en ecriture !" << endl;
	else
	{
		f << "<?xml version=\"1.0\" encoding=\"UTF-8\" ?>" << endl;
		f << "<file filename = \""<< filename << "\" created=\" \""<<endl;
		f << "<unit time=\"years\"  Masse=\"kg\" power=\"MW\" />"<<endl;
                f << "<material matid=\" \">"<<endl;
		f << "\t<time>" << endl;


		//Print all information about time(always the same for 0 to 20 yeaur) and the module involved


		for (int i =0; i < (int)vTime.size(); i++) {

			f << "\t\t<timestamp tid=\"" << i <<"\"  time=\""<< vTime[i] << "\" />" << endl;
		}
		f << "\t</time>" << endl;


		f<<"\t\t<Massehistory>"<<endl;




		//Print information about the masse(kg) of isotope in facility chossen



		for (int i =0; i < (int)vTime.size(); i++)
		{

			f << "\t\t\t<compositiondata tid=\""<< i <<"\">"<<endl;
			f << "\t\t\t\t<composition>"<<endl;



			for (int j =0; j< (int)toplot.size(); j++)
                        {

				string name;
				switch (toplot[j].fFacilityId)
				{
					case 0:
						switch (toplot[j].fFacylityNumber)
					{
						case 0:
							name = "TOTAL.";


							break;

						case 1:
							name = "INCYCLE.";


							break;

						case 2:
							name = "WASTE.";


							break;

						case 3:
							name = "OUTINCOME.";



							break;

						case 4:
							name = "REACTOR.";



							break;

						case 5:
							name = "COOLING.";



							break;

						case 6:
							name = "STOCK.";



							break;

						case 7:
							name = "FUELFABRICATION.";



							break;

						default:
							break;
					}
						break;

					case 1:
						name = "Reactor" + itoa(toplot[j].fFacylityNumber) + ".";



						break;

					case 2:
						name = "Storage" + itoa(toplot[j].fFacylityNumber) + ".";



						break;

					case 3:
						name = "Pool" + itoa(toplot[j].fFacylityNumber) + ".";




						break;

					case 4:
						name = "FabricationPlant" + itoa(toplot[j].fFacylityNumber) + ".";



						break;


					default:
						break;
				}


				f<<"\t\t\t\t\t<isotope zamid=\""<< itoa(toplot[j].fZAI.Z())<<itoa(toplot[j].fZAI.A())<<"\">"<<endl;
				f<<"\t\t\t\t\t\t<facility "<< name <<">"<<"<masse>"<<vQuantity[j][i]<<"</masse>"<<"</facility>"<<endl;


				f<<"\t\t\t\t\t</isotope>"<<endl;

			}

			f << "\t\t\t\t</composition>"<<endl;
			f << "\t\t\t</compositiondata>"<<endl;

		}

		f <<"\t\t</Massehistory>" << endl;
		f <<"</material>"<<endl;
		//------------------------------------------------------------------End of document XML

		f.close();
		cout << "  Your conversion is successful" << endl;
	};


	fData[toplot[0].fTreeId]->ResetBranchAddresses();
	{
		for(int i=0; i< (int)fReactorName[toplot[0].fTreeId].size(); i++) delete reactor[i];

		for(int i=0; i< (int)fPoolName[toplot[0].fTreeId].size(); i++) delete pool[i];

		for(int i=0; i< (int)fFabricationName[toplot[0].fTreeId].size(); i++) delete fabricationplant[i];

		for(int i=0; i< (int)fStockName[toplot[0].fTreeId].size(); i++) delete stock[i];
		for(int i=0; i< 8; i++) delete IV[i];
	}
	
	
	
	
}

//________________________________________________________________________
void CLASSRead::ConvertXmlMass(vector<CLASSPlotElement> toplot, string filename)
{
	
	vector<CLASSPlotElement> toplotTTree[fData.size()];
	
	
	
	fNumberGraphInvIterator = 0;
	for (int i = 0; i < (int)toplot.size(); i++)
	{
		toplotTTree[toplot[i].fTreeId].push_back(toplot[i]);
	}
	
	for (int i = 0; i < (int)fData.size(); i++)
	{
		if(toplotTTree[i].size() !=0)
			ConvertxmlTTreeMass(toplotTTree[i], filename);
		
	}
}
//________________________________________________________________________

//________________________________________________________________________

//________________________________________________________________________
string CLASSRead::GetBranchInName(CLASSPlotElement toplot)
{
	string name;
	switch (toplot.fFacilityId)
	{
		case 0:
			switch (toplot.fFacylityNumber)
		{
			case 0:
				name = "TOTAL.";
				return name;
				break;
				break;

			case 1:
				name = "INCYCLE.";
				return name;
				break;

			case 2:
				name = "WASTE.";
				return name;
				break;

			case 3:
				name = "OUTINCOME.";
				return name;
				break;

			case 4:
				name = "REACTOR.";
				return name;
				break;

			case 5:
				name = "COOLING.";
				return name;
				break;

			case 6:
				name = "STOCK.";
				return name;
				break;

			case 7:
				name = "FUELFABRICATION.";
				return name;
				break;

			default:
				break;
		}
			break;

		case 1:
			name = "R_" + fReactorName[toplot.fTreeId][toplot.fFacylityNumber];
			return name;
			break;

		case 2:
			name = "S_" + fStockName[toplot.fTreeId][toplot.fFacylityNumber];
			return name;
			break;

		case 3:
			name = "P_" + fPoolName[toplot.fTreeId][toplot.fFacylityNumber];
			return name;

			break;

		case 4:
			name = "F_" + fFabricationName[toplot.fTreeId][toplot.fFacylityNumber];
			return name;
			break;


		default:
			break;
	}

	return name;
}


string CLASSRead::GetLegendOutName(CLASSPlotElement toplot)
{
	string name;

	switch (toplot.fFacilityId)
	{
		case -2:
			name = "P_{" + itoa(toplot.fTreeId) + "} POWER ";
			return name;
			break;

		case 0:
			switch (toplot.fFacylityNumber)
		{



			case 0:
				name = "P_{" + itoa(toplot.fTreeId) + "} TOT ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				return name;
				break;

			case 1:
				name = "P_{" + itoa(toplot.fTreeId) + "} INcl ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				return name;
				break;

			case 2:
				name = "P_{" + itoa(toplot.fTreeId) + "} Wst ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				return name;
				break;

			case 3:
				name = "P_{" + itoa(toplot.fTreeId) + "} OUTINCOME ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				break;

			case 4:
				name = "P_{" + itoa(toplot.fTreeId) + "} R_{tot} ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				break;

			case 5:
				name = "P_{" + itoa(toplot.fTreeId) + "} Pl_{tot} ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				return name;
				break;

			case 6:
				name = "P_{" + itoa(toplot.fTreeId) + "} Stk_{tot} ^{" + itoa(toplot.fZAI.A()) + "}" + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				break;

			case 7:
				name = "P_{" + itoa(toplot.fTreeId) + "} FP_{tot} ^{" + itoa(toplot.fZAI.A()) + "}"  + ReadNucleusName[toplot.fZAI.Z()];
				for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
				return name;
				break;

			default:
				break;
		}
			break;

		case 1:
			name = "P_{" + itoa(toplot.fTreeId) + "} " + fReactorName[toplot.fTreeId][toplot.fFacylityNumber] + " ^{" + itoa(toplot.fZAI.A()) + "}" + ReadNucleusName[toplot.fZAI.Z()];
			for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
			break;

		case 2:
			name = "P_{" + itoa(toplot.fTreeId) + "} " + fStockName[toplot.fTreeId][toplot.fFacylityNumber] + " ^{" + itoa(toplot.fZAI.A()) + "}" + ReadNucleusName[toplot.fZAI.Z()];
			for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
			break;

		case 3:
			name = "P_{" + itoa(toplot.fTreeId) + "} " + fPoolName[toplot.fTreeId][toplot.fFacylityNumber] + " ^{" + itoa(toplot.fZAI.A()) + "}" + ReadNucleusName[toplot.fZAI.Z()];
			for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";

			break;

		case 4:
			name = "P_{" + itoa(toplot.fTreeId) + "} " + fFabricationName[toplot.fTreeId][toplot.fFacylityNumber] + " ^{" + itoa(toplot.fZAI.A()) + "}" + ReadNucleusName[toplot.fZAI.Z()];
			for (int i = 0; i < toplot.fZAI.I(); i++) name+= "*";
			break;


		default:
			break;
	}

	switch (toplot.fIVNumber)
	{

		case 0:
			name+= " Inside";
			return name;
			break;

		case 1:
			name+= " CumuIN";
			return name;
			break;

		case 2:
			name+= " CumuOUT";
			return name;
			break;

		default:
			break;

	}



	return name;

}

string CLASSRead::GetTittleOutName(CLASSPlotElement toplot)
{
	string name;

	switch (toplot.fFacilityId)
	{
		case 0:
			switch (toplot.fFacylityNumber)
		{
			case 0:
				name = "PARC "+ itoa(toplot.fTreeId) + " TOTAL " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;
				break;

			case 1:
				name = "PARC "+ itoa(toplot.fTreeId) +  " INCYCLE " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 2:
				name = "PARC "+ itoa(toplot.fTreeId) +  " WASTE " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 3:
				name = "PARC "+ itoa(toplot.fTreeId) +  " OUTINCOME " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 4:
				name = "PARC "+ itoa(toplot.fTreeId) +  " REACTOR " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 5:
				name = "PARC "+ itoa(toplot.fTreeId) +  " COOLING " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 6:
				name = "PARC "+ itoa(toplot.fTreeId) +  " STOCK " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			case 7:
				name = "PARC "+ itoa(toplot.fTreeId) +  " FUELFABRICATION " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
				return name;
				break;

			default:
				break;
		}
			break;

		case 1:
			name = "PARC "+ itoa(toplot.fTreeId) +  " " + fReactorName[toplot.fTreeId][toplot.fFacylityNumber] + " " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
			return name;
			break;

		case 2:
			name = "PARC "+ itoa(toplot.fTreeId) +  " " + fStockName[toplot.fTreeId][toplot.fFacylityNumber] + " " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
			return name;
			break;

		case 3:
			name = "PARC "+ itoa(toplot.fTreeId) +  " " + fPoolName[toplot.fTreeId][toplot.fFacylityNumber] + " " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
			return name;

			break;

		case 4:
			name = "PARC "+ itoa(toplot.fTreeId) +  " " + fFabricationName[toplot.fTreeId][toplot.fFacylityNumber] + " " + itoa(toplot.fZAI.Z()) + " " + itoa(toplot.fZAI.A()) + " " + itoa(toplot.fZAI.I());
			return name;
			break;
			
			
		default:
			break;
	}
	return name;
	
}

















