#ifndef _CLASSRead
#define _CLASSRead_

#include <TROOT.h>
#include <TTree.h>
#include <TString.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLatex.h>

#include <vector>
#include <iostream>

#include "ZAI.hxx"
#include "CLASSPlotElement.hxx"
#include "CLASSHeaders.hxx"
#include "Irradiation/IM_Matrix_Decay.hxx"

using namespace std;

	//________________________________________________________________________
	//
	//		CLASSRead
	//@{
	//	CLASSRead is a class that allows to retrieve results in Binary or Ascii Data
	// 	files written by CLASS. By default binary files are read.
	//
	// Thank to Francisco, Guerin, Louard, Shi for their participation to the XML ouput implementation
	// @author BLG
	// @author BAM
	// @version 1.0
	// @version 0.1
	//
	//@}
	//________________________________________________________________________
class CLASSRead
{
public :
	
	CLASSRead(TString filemname);
	~CLASSRead();		//@- destructor
	
	vector< vector< TString > >	GetReactorName()	{return fReactorName;}
	vector< vector< TString > >	GetPoolName()		{return fPoolName;}
	vector< vector< TString > >	GetFabricationName()	{return fFabricationName;}
	vector< vector< TString > >	GetStockName()		{return fStockName;}
	vector< ZAI >		GetZAIvector()		{return fZAIvector;}
	vector< TTree* >	GetData()		{return fData;}
	
	string GetBranchInName(CLASSPlotElement toplot);
	string GetLegendOutName(CLASSPlotElement toplot);
	string GetTittleOutName(CLASSPlotElement toplot);
	
	
	vector< vector<cSecond> > GetTimeVector() {return fTimeVector;}

	
	void ReadName();
	void ReadZAI();
	void ReadTime();

	void BuildTGraph(vector<CLASSPlotElement> toplot, int PlotId = 0, string opt = "L");
	void BuildTGraphUsingDecayChain(vector<CLASSPlotElement> toplot, int PlotId, int StartingStep, cSecond FinalTime, int StepNUmber, bool LinBin, string opt);

	
	void PlotInv(vector<CLASSPlotElement> toplot, bool DecayChain = false, int StartingStep = 0, cSecond FinalTime = 0, int StepNUmber = 0, bool LinBin = true , string opt = "L");
	
	void PlotTox(vector<CLASSPlotElement> toplot, bool DecayChain = false, int StartingStep = 0, cSecond FinalTime = 0, int StepNUmber = 0, bool LinBin = true , string opt = "L");
	void PlotHeat(vector<CLASSPlotElement> toplot, bool DecayChain = false, int StartingStep = 0, cSecond FinalTime = 0, int StepNUmber = 0, bool LinBin = true , string opt = "L");
	
	void PlotPower(vector<CLASSPlotElement> toplot, string opt = "L");
	void PlotTTreePower(vector<CLASSPlotElement> toplot, string opt = "L");

	void ConvertxmlTTreeMass(vector<CLASSPlotElement> toplot, string filename);
	void ConvertXmlMass(vector<CLASSPlotElement> toplot, string opt);           
	void ConvertxmlTTreePower(vector<CLASSPlotElement> toplotPower, string opt);
	void ConvertXmlPower(vector<CLASSPlotElement> toplotPower, string opt );

	void Write(string filename, string fileformat, string PadName);
	void ASCIIWrite(string filename , string PadName);


	void AddFile(TString filemname);
	TCanvas* fCNucleiInv;
	TCanvas* fCNucleiTox;
	TCanvas* fCNucleiHeat;
	TCanvas* fCPower;
	
	
private :
		
	
	vector< vector< TString > > fReactorName;
	vector< vector< TString > > fPoolName;
	vector< vector< TString > > fFabricationName;
	vector< vector< TString > > fStockName;
	vector< ZAI >	fZAIvector;
	
	
	vector< TTree* > fData;
	vector<TFile* >	 fFileIn;

	vector<TGraph*> fGraphInv;
	TGraph* fGraphInvSumOfSelected;
	TLatex** fLegendInv;
	TLatex* fLegendInvSumOfSelected;

	vector<TGraph*> fGraphTox;
	TGraph* fGraphToxSumOfSelected;
	TLatex** fLegendTox;
	TLatex* fLegendToxSumOfSelected;
	
	vector<TGraph*> fGraphHeat;
	TGraph* fGraphHeatSumOfSelected;
	TLatex** fLegendHeat;
	TLatex* fLegendHeatSumOfSelected;
	
	
	TGraph**  fGraphPower;
	TLatex** fLegendPower;
	int fNumberGraphPowerIterator;

	vector< vector<cSecond> > fTimeVector;

	
	double Xmin;
	double Xmax;
	double Ymin;
	double Ymax;

	
};



#endif
