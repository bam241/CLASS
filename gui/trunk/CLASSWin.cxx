#include "CLASSWin.hxx"

//#include "Convert.hxx"

#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <cstdlib>
#include <unistd.h>
#include <cmath>


#include <TApplication.h>
#include <TGTableLayout.h>
#include <TGResourcePool.h>
#include <TLine.h>
#include "StringLine.hxx"

// For all class not defined in file.hxx see http://root.cern.ch/root/html/
// For example : for TGraphErrors see http://root.cern.ch/root/html/TGraphErrors.html


using namespace std;
// for gcc3.2.3 only
char* operator+( std::streampos&, char* );
//end of gcc3.2.3

string NucleusName[] = {
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

const double densalmost0 = 1.0e-20;
const double Xlogmin = 1.0/365.25/24./3600.; // minimum time = 1s after core exit
const double Ylogmin = 1.0e-16;
const int WIDE = 0;

MainWin::MainWin(CLASSRead * DATA,vector<string> VFileName)
{
	fDATA=DATA;
	fSaveFileFormat = "ASCII";
	Start( VFileName);
}


MainWin::~MainWin()
{
	delete fDATA;
}


//_____________________________________________________________________________________________
void MainWin::Start(vector<string> VFileName)
{
	
	fToxNstep=50;		// default number of time step for spent fuel radiotoxicity calculations
	fToxTimeFirst=1.0;	// default first time for spent fuel radiotoxicity calculations
	fToxTimeLast=1.0E08;	// default last time for spent fuel radiotoxicity calculations
	fMotherIsVisible=false;
	fIsByMother=false;
	fIsLinear=false;
	//
	// Set the line width for all graph: for PPT presentation a value of 2 is probably better
	//
	fGraphLineWidth = 1;  //default=1
	//
	// Set the marker size for all graph: for PPT presentation a value of 1 is probably better
	//
	fGraphMarkerSize = 0.5;   //default=0.5
	
	//const TGFont *fontS = gClient->GetFont("-*-helvetica-medium-r-*-*-8-*-*-*-*-*-*");
	const TGFont *fontS = gClient->GetFont("-*-helvetica-medium-r-*-*-10-*-*-*-*-*-*");
	if (!fontS)
		fontS = gClient->GetResourcePool()->GetDefaultFont();
	fLabelFontS = fontS->GetFontStruct();
	
	//	const TGFont *fontB = gClient->GetFont("-*-helvetica-medium-r-*-*-10*-*-*-*-*-*-*");
	const TGFont *fontB = gClient->GetFont("-*-helvetica-medium-r-*-*-12*-*-*-*-*-*-*");
	if (!fontB)
		fontB = gClient->GetResourcePool()->GetDefaultFont();
	fLabelFontB = fontB->GetFontStruct();
	
	
	/*****************************/
	//Compteur d'objet selectionné
	/*****************************/
	fNselectedNucleus = 0;	// number of nucleus selected
	
	fNumberOfParc = fDATA->GetData().size();			//@@@nombre de parcs
	fNumberOfTOT = 9; //Reacteur + Stock + Pooling + Fabrication + Power
	fNumberOfReactor = new int[fNumberOfParc];
	fNumberOfStock = new int[fNumberOfParc];
	fNumberOfPool = new int[fNumberOfParc];
	fNumberOfFab = new int[fNumberOfParc];
	
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fNumberOfReactor[i] = fDATA->GetReactorName()[i].size();
		fNumberOfStock[i] = fDATA->GetStockName()[i].size();
		fNumberOfPool[i] = fDATA->GetPoolName()[i].size();
		fNumberOfFab[i] = fDATA->GetFabricationName()[i].size();
	}
	
	// Canvas Style
	gStyle->SetOptStat(000);
	gStyle->SetFrameFillColor(10);
	gStyle->SetFrameFillStyle(6);
	gStyle->SetFrameBorderMode(0);
	gStyle->SetFrameBorderSize(0);
	gStyle->SetCanvasBorderMode(0);
	gStyle->SetCanvasBorderSize(0);
	gStyle->SetCanvasColor(10);
	gStyle->SetPadColor(10);
	gStyle->SetPadBorderMode(0);
	gStyle->SetPadBorderSize(0);
	gStyle->SetTitleFont(22);
	gStyle->SetLabelFont(22,"xyz");
	gStyle->SetHistLineWidth(1);
	
	
	//first define Layaout (left, rigth, top, bottom margins)
	TGLayoutHints* fL5555 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX , 5, 5, 5, 5);
	TGLayoutHints* fL2222 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2, 2, 2, 2);
	TGLayoutHints* fL2200 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2, 2, 0, 0);
	TGLayoutHints* fL55100 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,5, 5, 0, 0);
	//
	// the Main frame where all will be insert
	//
	
	fGeneF0 = new TGVerticalFrame(this, 100, 150);
	AddFrame(fGeneF0, fL5555);
	SetHeight(150);
	SetWidth(60);
	int MAXPATHLEN=256;
	char cDir[MAXPATHLEN];
	getcwd(cDir, MAXPATHLEN);
	stringstream tmp;
	tmp.str("");
	tmp << "CLASSGui " ;
	for(int i=0;i<int(VFileName.size());i++)
		tmp <<VFileName[i]<<" " ;
	
	this->SetWindowName(tmp.str().c_str());
	
	
	/***************** LES Differents Parcs  ***********************/
	
	fParcTab = new TGTab(fGeneF0,100,150);
	fGeneF0->AddFrame(fParcTab,fL2222);
	fParcTab->Associate(this);
	
	fParcTabFoil = new TGCompositeFrame*[fNumberOfParc];
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fParcTabFoil[i] = new TGCompositeFrame;
		stringstream ParcName;
		ParcName<<"Park "<<i;
		fParcTabFoil[i] = fParcTab->AddTab(ParcName.str().c_str());
	}
	
	
	/*************Les differentes Facilities Tab****/
	fFacilitiesTab = new TGTab*[fNumberOfParc];
	
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fFacilitiesTab[i] = new TGTab(fParcTabFoil[i],100,150);
		fParcTabFoil[i]->AddFrame(fFacilitiesTab[i],fL2222);
		fFacilitiesTab[i]->Associate(this);
		
	}
	
	fFacilitiesTabFoil = new TGCompositeFrame**[fNumberOfParc];
	//1 jeu de facility tab par Parc
	for(int i = 0; i < fNumberOfParc; i++)
	{
		
		fFacilitiesTabFoil[i] = new TGCompositeFrame*[5];
		
		fFacilitiesTabFoil[i][0] = new TGCompositeFrame;
		fFacilitiesTabFoil[i][0] = fFacilitiesTab[i]->AddTab("Total");
		fFacilitiesTabFoil[i][1] = new TGCompositeFrame;
		fFacilitiesTabFoil[i][1] = fFacilitiesTab[i]->AddTab("Reactor(s)");
		fFacilitiesTabFoil[i][2] = new TGCompositeFrame;
		fFacilitiesTabFoil[i][2] = fFacilitiesTab[i]->AddTab("Stock(s)");
		fFacilitiesTabFoil[i][3] = new TGCompositeFrame;
		fFacilitiesTabFoil[i][3] = fFacilitiesTab[i]->AddTab("Pool(s)");
		fFacilitiesTabFoil[i][4] = new TGCompositeFrame;
		fFacilitiesTabFoil[i][4] = fFacilitiesTab[i]->AddTab("Fabrication Plant(s)");
	}
	
	//1 jeu d'ItemTab par facility tab
	
	fItemTab = new TGTab**[fNumberOfParc];
	
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fItemTab[i] = new TGTab*[5];
		for(int j = 0; j < 5; j++)
		{
			fItemTab[i][j] = new TGTab(fFacilitiesTabFoil[i][j],100,150);//@@ hxx non défini
			fFacilitiesTabFoil[i][j]->AddFrame(fItemTab[i][j]);
			fItemTab[i][j]->Associate(this);
			//construire les foils en fonction du nombre de item possible par foil
			
		}
		
	}
	
	for(int j = 0; j < 5; j++)
		FillItemTab(j);
	
	/***************** What you wana plot ***********************/
	
	//Configuration Table
	fPlotConfigTab = new TGTab(fGeneF0,100,150);
	fGeneF0->AddFrame(fPlotConfigTab,fL2222);
	fPlotConfigTab->Associate(this);
	
	fPlotConfigFoil = new TGCompositeFrame*[2];
	
	fPlotConfigFoil[0] = new TGCompositeFrame;
	fPlotConfigFoil[0] = fPlotConfigTab->AddTab("Inventory");
	
	fPlotConfigFoil[1] = new TGCompositeFrame;
	fPlotConfigFoil[1] = fPlotConfigTab->AddTab("DecayHeat / Radiotox.");
	
	
	//In inventory
	fInventoryFrame = new TGGroupFrame(fPlotConfigFoil[0],"", kHorizontalFrame );
	fPlotConfigFoil[0]->AddFrame(fInventoryFrame,fL5555);
	
	fCheckIVPlot = new TGCheckButton*[3];
	fCheckIVPlot[0] = new TGCheckButton(fInventoryFrame,"Inside",M_CHECK_INSIDE);
	fCheckIVPlot[0]->SetFont(fLabelFontS);
	fInventoryFrame->AddFrame(fCheckIVPlot[0],fL2222);
	fCheckIVPlot[0]->Associate(this);
	fCheckIVPlot[0]->SetState(kButtonDown);
	
	fCheckIVPlot[1] = new TGCheckButton(fInventoryFrame,"Cumul In",M_CHECK_CUMIN);
	fCheckIVPlot[1]->SetFont(fLabelFontS);
	fInventoryFrame->AddFrame(fCheckIVPlot[1],fL2222);
	fCheckIVPlot[1]->Associate(this);
	
	fCheckIVPlot[2] = new TGCheckButton(fInventoryFrame,"Cumul out",M_CHECK_CUMOUT);
	fCheckIVPlot[2]->SetFont(fLabelFontS);
	fInventoryFrame->AddFrame(fCheckIVPlot[2],fL2222);
	fCheckIVPlot[2]->Associate(this);
	
	//Radio/Heat frame
	//Radio Or Decay subframe
	fDecayOrRadioFrame = new  TGGroupFrame(fPlotConfigFoil[1],"Radio-toxicity / decay heat" );
	fPlotConfigFoil[1]->AddFrame(fDecayOrRadioFrame,fL5555);
	
	fButtonHeat=new TGRadioButton(fDecayOrRadioFrame,"Decay Heat",M_RADIO_DECAY_HEAT);
	fButtonRadiotox=new TGRadioButton(fDecayOrRadioFrame,"Radio-toxicity",M_RADIO_RADIOTOX);
	
	fButtonHeat->SetState(kButtonDown);
	fButtonRadiotox->SetState(kButtonUp);
	
	fButtonHeat->Associate(this);
	fButtonRadiotox->Associate(this);
	
	fDecayOrRadioFrame->AddFrame(fButtonHeat);
	fDecayOrRadioFrame->AddFrame(fButtonRadiotox);
	
	
	//By mother sub frame
	fByMotherFrame = new  TGGroupFrame(fGeneF0,"Decay chain (by mother)" );
	fGeneF0->AddFrame(fByMotherFrame,fL2222);
	
	fByMotherMore = new TGPictureButton(fGeneF0 ,gClient->GetPicture("arrow_down.xpm"),M_BUTTON_MOTHER_MORE);
	fByMotherMore->Resize(350,25);
	fByMotherMore->Associate(this);
	fGeneF0->AddFrame(fByMotherMore,fL2222);
	
	//the button
	fByMotherButton=new TGCheckButton(fByMotherFrame,"By Mother",M_CHECK_BY_MOTHER);
	fButtonHeat->SetState(kButtonUp);
	fByMotherButton->Associate(this);
	fByMotherFrame->AddFrame(fByMotherButton);
	
	
	//the time choosed for the end of the scenario
	fScenarTimeFrame = new TGGroupFrame(fByMotherFrame,"Final scenario time (year)" );
	((TGGroupFrame*)fScenarTimeFrame)->SetTextFont(fLabelFontS);
	fByMotherFrame->AddFrame(fScenarTimeFrame,fL2222);
	fScenarTimeCBox = new TGComboBox(fScenarTimeFrame,M_CB_SCENAR_Time);
	fScenarTimeFrame->AddFrame(fScenarTimeCBox,fL2222);
	fScenarTimeCBox->Resize(80,17);
	
	int NumOfTimeStep = fDATA->GetTimeVector()[0].size();
	for(int step=0; step<NumOfTimeStep; step++)
	{
		stringstream tmp;
		tmp.str("");
		tmp<<"Dir#"<< step <<" @ "<<fDATA->GetTimeVector()[0][step]/cYear;
		fScenarTimeCBox->AddEntry(tmp.str().c_str(),step);
	}
	
	fScenarTimeCBox->Associate(this);
	fTimeStep = NumOfTimeStep-1;
	fScenarTimeCBox->Select(fTimeStep);
	
	// Parameters for the geological time
	fTimeParametersFrame = new TGGroupFrame(fByMotherFrame,"Evol. Period [year]: first, last, n_step)", kHorizontalFrame );
	((TGGroupFrame*)fTimeParametersFrame)->SetTextFont(fLabelFontS);
	fByMotherFrame->AddFrame(fTimeParametersFrame,fL2222);
	
	TGTextBuffer *TBtoxfirst = new TGTextBuffer(100);
	TGTextBuffer *TBtoxlast = new TGTextBuffer(100);
	TGTextBuffer *TBtoxnstep = new TGTextBuffer(100);
	
	TBtoxfirst->AddText(0, StringLine::convert<string>(fToxTimeFirst).c_str());
	TBtoxlast->AddText(0, StringLine::convert<string>(fToxTimeLast).c_str());
	TBtoxnstep->AddText(0, StringLine::convert<string>(fToxNstep).c_str());
	
	TEtoxfirst = new TGTextEntry(fTimeParametersFrame,TBtoxfirst , M_TE_toxfirst);// box to enter the first time
	TEtoxlast = new TGTextEntry(fTimeParametersFrame,TBtoxlast , M_TE_toxlast);// box to enter the last time
	TEtoxnstep = new TGTextEntry(fTimeParametersFrame,TBtoxnstep , M_TE_toxnstep);// box to enter the number of step time
	
	fCheckLinear = new TGCheckButton(fTimeParametersFrame,"Linear Binning",M_CHECK_LINEAR_Tox);
	
	TEtoxfirst->Resize(TEtoxfirst->GetDefaultWidth(),17);
	TEtoxlast->Resize(TEtoxlast->GetDefaultWidth(),17);
	TEtoxnstep->Resize(TEtoxnstep->GetDefaultWidth(),17);
	
	fTimeParametersFrame->AddFrame(TEtoxfirst, fL2222);
	fTimeParametersFrame->AddFrame(TEtoxlast, fL2222);
	fTimeParametersFrame->AddFrame(TEtoxnstep, fL2222);
	fTimeParametersFrame->AddFrame(fCheckLinear,fL2222);
	
	fCheckLinear->SetFont(fLabelFontS);
	fCheckLinear->SetState(kButtonUp);
	
	TEtoxfirst->Associate(this);
	TEtoxlast->Associate(this);
	TEtoxnstep->Associate(this);
	fCheckLinear->Associate(this);
	
	//
	//	The Miscellaneous part (AM check box, ...)
	//
	fMiscFrame = new TGGroupFrame(fGeneF0,"Miscellaneous", kVerticalFrame );
	((TGGroupFrame*)fMiscFrame)->SetTextFont(fLabelFontB);
	fGeneF0->AddFrame(fMiscFrame,fL2222);
	
	fMiscHzFrame=new TGHorizontalFrame(fMiscFrame,400, 50 );
	fMiscFrame->AddFrame(fMiscHzFrame,fL2222);
	
	//misc check buttons
	fCheckAMNuc=new TGCheckButton(fMiscHzFrame,"MA",M_CHECK_AM_NUC);
	fCheckFPNuc=new TGCheckButton(fMiscHzFrame,"FP",M_CHECK_FP_NUC);
	fCheckSumOfSelected=new TGCheckButton(fMiscHzFrame,"Sum of Selected");
	
	fCheckFPNuc->SetFont(fLabelFontS);
	fCheckAMNuc->SetFont(fLabelFontS);
	fCheckSumOfSelected->SetFont(fLabelFontS);
	
	fMiscHzFrame->AddFrame(fCheckAMNuc,fL2222);
	fMiscHzFrame->AddFrame(fCheckFPNuc,fL2222);
	fMiscHzFrame->AddFrame(fCheckSumOfSelected,fL2222);
	
	fCheckFPNuc->Associate(this);
	fCheckAMNuc->Associate(this);
	fCheckSumOfSelected->Associate(this);
	
	//1 jeu de NucleusTab pour tous
	
	//
	//	The Nucleus Tab
	//
	fTabNuc = new TGTab(fGeneF0,200,150);
	fGeneF0->AddFrame(fTabNuc,fL2222);
	fTabNuc->Associate(this);
	FillNucTab();
	
	// The Plot, Save, Macro and Quit buttons
	//
	fPlotSaveQuitFrame = new TGHorizontalFrame(fGeneF0, 100, 50 );
	fGeneF0->AddFrame(fPlotSaveQuitFrame,fL5555);
	fButtonPlot = new TGTextButton(fPlotSaveQuitFrame,"&Plot (All)",M_BUTTON_PLOT);
	fButtonSave = new TGTextButton(fPlotSaveQuitFrame,"&Save Data",M_BUTTON_SAVE);
	fButtonQuit = new TGTextButton(fPlotSaveQuitFrame,"&Quit",M_BUTTON_QUIT);
	fPlotSaveQuitFrame->AddFrame(fButtonPlot, fL2222);
	fPlotSaveQuitFrame->AddFrame(fButtonSave, fL2222);
	fPlotSaveQuitFrame->AddFrame(fButtonQuit, fL2222);
	fPlotSaveQuitFrame->Resize(150, fButtonPlot->GetDefaultHeight());
	fButtonPlot->Associate(this);
	fButtonSave->Associate(this);
	fButtonQuit->Associate(this);
	
	MapSubwindows();						// Map all sub windows that are part of the composite frame
	MapWindow();
	Resize(GetDefaultSize()); 					// fit to the exact size
	Move(410,30);
	fMainWidth=fGeneF0->GetWidth();
	fGeneF0->HideFrame(fByMotherFrame);
	//Resize(550,670); 					// fit to the exact size
	Resize(GetDefaultSize());
	fMainWidth=fGeneF0->GetWidth();
	
}
//_____________________________________________________________________________________________
bool MainWin::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
	SubWin *nsw = 0;
	string EnterTextP;
	switch (GET_MSG(msg))
	{
		case kC_COMMAND:
			switch (GET_SUBMSG(msg))
		{
			case kCM_BUTTON:
				switch (parm1)
			{
				case M_BUTTON_QUIT:
					CloseWindow();
					break;
				case M_BUTTON_PLOT:
					Plot();
					break;
				case M_BUTTON_SAVE:
					if(!gPad)
					{
						cout << "Please Plot something to save before" << endl;
						break;
					}
					
					nsw = new SubWin("SaveAs", fClient->GetRoot(), this, 200, 200);
					if (fSaveFileFormat == "XML")
						Conversionxml();
					else
						if(fSaveFileName != "" /*&& fSaveFileFormat !=""*/)
						{	//si plusieur canvas ouvert sauver l'actif exemple :
							string PadName=gPad->GetName();
							if(PadName=="c_Nuclei")
								fDATA->Write(fSaveFileName, fSaveFileFormat);
						}
					break;
					
				case M_BUTTON_MOTHER_MORE:
					if(!fMotherIsVisible)
					{
						fByMotherMore->SetPicture(gClient->GetPicture("arrow_up.xpm"));
						fByMotherMore->Resize(350,10);
						MapSubwindows();
						MapWindow();
						Resize(fMainWidth,0);
						fGeneF0->ShowFrame(fByMotherFrame);
						fMotherIsVisible=true;
					}
					else
					{
						fByMotherMore->SetPicture(gClient->GetPicture("arrow_down.xpm"));
						fByMotherMore->Resize(350,25);
						MapSubwindows();
						MapWindow();
						fGeneF0->HideFrame(fByMotherFrame);
						Resize(fMainWidth,0);
						fMotherIsVisible=false;
					}
					
					
			}
				
			case kCM_CHECKBUTTON:
				switch (parm1)
			{
				case M_CHECK_INSIDE:
					fButtonHeat->SetState(kButtonUp);
					fButtonRadiotox->SetState(kButtonUp);
					break;
					
				case M_CHECK_CUMIN:
					fButtonHeat->SetState(kButtonUp);
					fButtonRadiotox->SetState(kButtonUp);
					break;
					
				case M_CHECK_CUMOUT:
					fButtonHeat->SetState(kButtonUp);
					fButtonRadiotox->SetState(kButtonUp);
					break;
					
				case M_CHECK_BY_MOTHER:
					
					if(fByMotherButton->GetState()==kButtonDown )
						fIsByMother=true;
					else
						fIsByMother=false;
					break;
					
				case M_CHECK_LINEAR_Tox:
					if(fCheckLinear->GetState()==kButtonDown)
						fIsLinear=true;
					else
						fIsLinear=false;
					
					break;
			}
			case kCM_RADIOBUTTON:
				switch (parm1)
			{
				case M_RADIO_DECAY_HEAT:
					fButtonHeat->SetState(kButtonDown);
					fButtonRadiotox->SetState(kButtonUp);
					for(int i = 0 ; i<3 ; i++ )
						fCheckIVPlot[i]->SetState(kButtonUp);
					break;
					
				case M_RADIO_RADIOTOX:
					fButtonHeat->SetState(kButtonUp);
					fButtonRadiotox->SetState(kButtonDown);
					for(int i = 0 ; i<3 ; i++ )
						fCheckIVPlot[i]->SetState(kButtonUp);
					break;
					
			}
		}
		case kC_TEXTENTRY:
			switch (GET_SUBMSG(msg))
		{
			case kTE_TEXTCHANGED:
				switch (parm1)
			{
				case M_TE_toxfirst:   	// get first time
					EnterTextP=TEtoxfirst->GetBuffer()->GetString();
					fToxTimeFirst=StringLine::convert<double>(EnterTextP);		// change string in double
					if(fToxTimeFirst<=0)
					{
						fToxTimeFirst=Xlogmin;
						TEtoxfirst->SetText(StringLine::convert<string>(fToxTimeFirst).c_str());
					}
					break;
					
				case M_TE_toxlast:   	// get last time
					
					EnterTextP=TEtoxlast->GetBuffer()->GetString();
					fToxTimeLast=StringLine::convert<double>(EnterTextP);
					break;
					
				case M_TE_toxnstep:   	// get number of time steps
					
					EnterTextP=TEtoxnstep->GetBuffer()->GetString();
					fToxNstep=StringLine::convert<int>(EnterTextP);
					break;
					
			}
			case kCM_COMBOBOX:
				switch (parm1)
			{
				case M_CB_SCENAR_Time:
					fTimeStep = fScenarTimeCBox->GetSelected();
					break;
			}
				break;
		}
			break;
	}
	
	return kTRUE;
}

//_____________________________________________________________________________________________
void MainWin::CloseWindow()
{
	gApplication->Terminate(0);
}
//_____________________________________________________________________________________________
void MainWin::Plot()
{
	int Nnucleus = fDATA->GetZAIvector().size();
	
	vector<CLASSPlotElement> toplot;
	vector<CLASSPlotElement> toplotPower;
	
	if(fCheckSumOfSelected->GetState()==kButtonDown)
	{
		toplot.push_back( CLASSPlotElement(-1, -1, -1, -1,-1,-1,-1) );
	}
	//Power
	for(int i=0; i < fNumberOfParc; i++)
	{
		if(fCheckArrayTotal[i][fNumberOfTOT-1]->GetState()==kButtonDown)
			toplotPower.push_back( CLASSPlotElement(i, -2, -2, -2,-2,-2,-2) );
		
	}
	
	if(toplotPower.size() != 0)
		fDATA->PlotPower(toplotPower);
	
	
	for(int i=0; i < fNumberOfParc; i++)
	{
		
		for(int j=0; j < fNumberOfTOT-1; j++) //fNumberOfTOT -1 ?? ->All except power
		{
			if(fCheckArrayTotal[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 0, j,0, fDATA->GetZAIvector()[k]));
				}
		}
		
		for(int j=0; j < fNumberOfReactor[i]; j++)
		{
			if(fCheckArrayReactor[i][j]->GetState()==kButtonDown)
			{
				for(int l=0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i,1,j,l, fDATA->GetZAIvector()[k]));
						}
			}
		}
		
		for(int j=0; j < fNumberOfStock[i]; j++)
		{
			if(fCheckArrayStock[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 2, j,l, fDATA->GetZAIvector()[k]));
						}
		}
		
		for(int j=0; j < fNumberOfPool[i]; j++)
		{
			if(fCheckArrayPool[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 3, j,l, fDATA->GetZAIvector()[k]));
						}
		}
		
		for(int j=0; j < fNumberOfFab[i]; j++)
		{
			if(fCheckArrayFab[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 4, j,l, fDATA->GetZAIvector()[k]));
						}
		}
	}
	
	
	int StartingStep = fTimeStep;
	cSecond FinalTime =(cSecond) fToxTimeLast * cYear;
	int NStep = fToxNstep;
	bool IsLinear = fIsLinear;
	if(!fIsByMother)
	{
		StartingStep  = 0;
		FinalTime = 0;
		NStep = 0;
		IsLinear = true;
	}
	
	
	if(toplot.size() != 0)
	{
		for(int i = 0 ; i<3 ; i++ )
		{
			if(fCheckIVPlot[i]->GetState()==kButtonDown)
			{	fDATA->PlotInv(toplot,fIsByMother,StartingStep,FinalTime,NStep,IsLinear);
				break;
			}
		}
		
		if(fButtonRadiotox->GetState()==kButtonDown)
			fDATA->PlotTox(toplot,fIsByMother,StartingStep,FinalTime,NStep,IsLinear);
		
		if(fButtonHeat->GetState()==kButtonDown)
			fDATA->PlotHeat(toplot,fIsByMother,StartingStep,FinalTime,NStep,IsLinear);
		
	}
	
}
void MainWin::Conversionxml()
{
	int Nnucleus = fDATA->GetZAIvector().size();
	
	vector<CLASSPlotElement> toplot;
	vector<CLASSPlotElement> toplotPower;
	
	for(int i=0; i < fNumberOfParc; i++)
	{
		
		for(int j=0; j < fNumberOfTOT-1; j++) //fNumberOfTOT -1 ?? ->All except power
		{
			if(fCheckArrayTotal[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 0, j,0, fDATA->GetZAIvector()[k]));
				}
		}
		
		for(int j=0; j < fNumberOfReactor[i]; j++)
		{
			if(fCheckArrayReactor[i][j]->GetState()==kButtonDown)
			{
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 1, j,l, fDATA->GetZAIvector()[k]));
						}
			}
		}
		
		for(int j=0; j < fNumberOfStock[i]; j++)
		{
			if(fCheckArrayStock[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 2, j,l, fDATA->GetZAIvector()[k]));
						}
		}
		
		for(int j=0; j < fNumberOfPool[i]; j++)
		{
			if(fCheckArrayPool[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 3, j,l, fDATA->GetZAIvector()[k]));
						}
		}
		
		for(int j=0; j < fNumberOfFab[i]; j++)
		{
			if(fCheckArrayFab[i][j]->GetState()==kButtonDown)
				for(int l =0;l<3;l++)
					if (fCheckIVPlot[l]->GetState()==kButtonDown)
						for(int k=0; k < Nnucleus; k++)
						{
							if(fCheckArrayNuc[k]->GetState()==kButtonDown)
								toplot.push_back( CLASSPlotElement(i, 4, j,l, fDATA->GetZAIvector()[k]));
						}
		}
	}
	
	if(toplot.size() != 0)
		fDATA->ConvertXmlMass(toplot, fSaveFileName);
	
}

//_____________________________________________________________________________________________
void MainWin::FillNucTab()		// fill the Inventory Tab foil
{
	static bool first=true;
	
	int NCheck = fDATA->GetZAIvector().size();
	//number of lines and columns per foil
	int Nline = 14-WIDE;
	int Ncol = 8+WIDE;
	//number of FP tab
	static int NFP_tab = 0;
	static int NTab = 0;
	static int NACT_tab = 0;
	static int NGAZ_tab = 0;
	
	if(first)
	{
		NTab = NCheck/(Nline*Ncol)+1;
	}
	//Are we resorting nuclei ?
	if(!first)
	{
		for(int l = 0; l < NTab ; l++)
			fTabFoilNuc[l]->Cleanup();
		delete [] fCheckArrayNuc;
	}
	// the foils
	if(first) fTabFoilNuc = new TGCompositeFrame*[NTab];
	string *TabName = new string[NTab];
	for(int i = 0; i < NTab; i++)
	{
		stringstream tmp;
		tmp.str("");
		tmp << "Nuclei " << i;
		TabName[i] = tmp.str();
	}
	
	// the Check buttons
	fCheckArrayNuc = new TGCheckButton*[NCheck];
	
	if( first )
		for(int l = 0; l < NTab ; l++)
		{
			fTabFoilNuc[l] = fTabNuc->AddTab(TabName[l].c_str());
			fTabFoilNuc[l]->SetLayoutManager(new TGMatrixLayout(fTabFoilNuc[l], Nline, 0, 5));
		}
	
	//
	// Fill the Miscellaneous TAB
	//
	FillNucFoil( NCheck, Ncol, Nline );
	//
	// Remap after resorting
	//
	for(int l = 0; l < NTab; l++)
	{
		fTabFoilNuc[l]->MapSubwindows();
		fTabFoilNuc[l]->Layout();
	}
	//select the second foil then the first foil to avoid "save char in back ground storage"
	first=false;
	
}

//_____________________________________________________________________________________________
void MainWin::FillNucFoil(int n_item, int Ncol,int Nline)
{
	int current_foil = 0;
	int current_item = 0;
	for( int l = n_item-1; l >= 0; l-- )
	{
		
		stringstream name;
		int Atmp = fDATA->GetZAIvector()[l].A();
		if(Atmp !=0)
		{
			
			if(fDATA->GetZAIvector()[l].Z() == -3)
			{
				name << Atmp << "TMP"; //@@@@
				
			}
			else if(fDATA->GetZAIvector()[l].Z() == -2)
			{
				name << Atmp << "PF"; //@@@@
				
			}
			else if(fDATA->GetZAIvector()[l].Z() == -1)
			{
				name << Atmp << "ERR"; //@@@@
				
			}
			else
				name << Atmp << NucleusName[fDATA->GetZAIvector()[l].Z()]; //@@@@
			
			if( fDATA->GetZAIvector()[0].I() > 0 )
				name << "*";
			
			
			fCheckArrayNuc[l] = new TGCheckButton( fTabFoilNuc[current_foil], name.str().c_str(), l + M_CHECK_PLOTALL+1 );
			fCheckArrayNuc[l]->SetFont( fLabelFontS );
			
			fTabFoilNuc[current_foil]->AddFrame( fCheckArrayNuc[l] );
			fCheckArrayNuc[l]->Associate( this );
			
			current_item++;
			
			if( current_item >= Ncol * Nline )
			{
				current_item = 0;
				current_foil++;
			}
		}
	}
	
	current_foil++;
	
}
//_____________________________________________________________________________________________
void MainWin::FillItemTab(int current)
{
	if (current==0)
		FillTotalTab();
	if (current==1)
		FillReactorTab();
	if (current==2)
		FillStockTab();
	if (current==3)
		FillPoolTab();
	if (current==4)
		FillFabricationTab();
	
	if(current>5 || current<0)
	{cout<<"BUG"<<endl; exit(0);}
	
	
}

//_____________________________________________________________________________________________
void MainWin::FillTotalTab()
{
	fCheckArrayTotal = new TGCheckButton**[fNumberOfParc];
	for(int i=0;i<	fNumberOfParc;i++)
		fCheckArrayTotal[i]= new TGCheckButton*[fNumberOfTOT];
	double count=0;
	for(int p=0; p<fNumberOfParc; p++)
	{
		for(int i=0; i < fNumberOfTOT; i++)
		{
			string name;
			if(i==0) name="TOTAL";
			if(i==1) name="INCYCLE";
			if(i==2) name="WASTE";
			if(i==3) name="OUTINCOME";
			if(i==4) name="REACTOR";
			if(i==5) name="COOLING";
			if(i==6) name="STOCK";
			if(i==7) name="FUELFABRICATION";
			if(i==8) name="POWER";
			
			fCheckArrayTotal[p][i]=new TGCheckButton(fFacilitiesTabFoil[p][0],name.c_str());
			fCheckArrayTotal[p][i]->SetFont(fLabelFontS);
			fFacilitiesTabFoil[p][0]->AddFrame(fCheckArrayTotal[p][i]);
			fCheckArrayTotal[p][i]->Associate(this);
			count++;
		}
		
		fFacilitiesTabFoil[p][0]->MapSubwindows();
		fFacilitiesTabFoil[p][0]->Layout();
	}
}
//_____________________________________________________________________________________________
void MainWin::FillReactorTab()
{
	int Nline=7-WIDE;
	int Ncol=6+WIDE;
	
	fTabFoilReactor= new TGCompositeFrame**[fNumberOfParc];
	fCheckArrayReactor= new TGCheckButton**[fNumberOfParc];
	int *NTab=new int[fNumberOfParc];
	
	for (int i=0 ; i<fNumberOfParc;i++)
	{
		
		NTab[i]=fNumberOfReactor[i]/(Nline*Ncol)+1;
		fTabFoilReactor[i]= new TGCompositeFrame*[NTab[i]];
		fCheckArrayReactor[i]= new TGCheckButton*[fNumberOfReactor[i]];
		
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		string TabName[NTab[p]];
		
		for(int l=0; l<NTab[p] ; l++)
		{
			stringstream tmp;
			tmp<<"Reactor "<<l;
			TabName[l]=tmp.str();
			//cout<<l<<" "<<TabName[l]<<endl;
			
			fTabFoilReactor[p][l] = fItemTab[p][1]->AddTab(TabName[l].c_str());
			fTabFoilReactor[p][l]->SetLayoutManager(new TGMatrixLayout(fTabFoilReactor[p][l], Nline, 0, 5));
		}
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		int current_foil=0;
		int current_item=0;
		int current_item_in_the_foil=0;
		for(int n=0;n<fNumberOfReactor[p];n++)
		{
			
			fCheckArrayReactor[p][current_item]= new TGCheckButton(fTabFoilReactor[p][current_foil],fDATA->GetReactorName()[p][n]);
			fCheckArrayReactor[p][current_item]->SetFont(fLabelFontS);
			
			fTabFoilReactor[p][current_foil]->AddFrame(fCheckArrayReactor[p][current_item]);
			fCheckArrayReactor[p][current_item]->Associate(this);
			
			current_item_in_the_foil++;
			current_item++;
			if(current_item_in_the_foil>=Ncol*Nline)
			{
				current_item_in_the_foil=0;
				current_foil++;
			}
			
			
		}
		
	}
}
//_____________________________________________________________________________________________
void MainWin::FillStockTab()
{
	int Nline=7-WIDE;
	int Ncol=6+WIDE;
	
	fTabFoilStock= new TGCompositeFrame**[fNumberOfParc];
	fCheckArrayStock= new TGCheckButton**[fNumberOfParc];
	int *NTab=new int[fNumberOfParc];
	
	for (int i=0 ; i<fNumberOfParc;i++)
	{
		NTab[i]=fNumberOfStock[i]/(Nline*Ncol)+1;
		fTabFoilStock[i]= new TGCompositeFrame*[NTab[i]];
		fCheckArrayStock[i]= new TGCheckButton*[fNumberOfStock[i]];
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		string TabName[NTab[p]];
		for(int l=0; l<NTab[p] ; l++)
		{
			stringstream tmp;
			tmp.str("");
			tmp<<"Stock "<<l;
			TabName[l]=tmp.str();
			fTabFoilStock[p][l]=fItemTab[p][2]->AddTab(TabName[l].c_str());
			fTabFoilStock[p][l]->SetLayoutManager(new TGMatrixLayout(fTabFoilStock[p][l], Nline, 0, 5));
		}
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		int current_foil=0;
		int current_item=0;
		int current_item_in_the_foil=0;
		for(int n=0;n<fNumberOfStock[p];n++)
		{
			
			fCheckArrayStock[p][current_item]= new TGCheckButton(fTabFoilStock[p][current_foil],fDATA->GetStockName()[p][n]);
			fCheckArrayStock[p][current_item]->SetFont(fLabelFontS);
			
			fTabFoilStock[p][current_foil]->AddFrame(fCheckArrayStock[p][current_item]);
			fCheckArrayStock[p][current_item]->Associate(this);
			
			current_item_in_the_foil++;
			current_item++;
			if(current_item_in_the_foil>=Ncol*Nline)
			{
				current_item_in_the_foil=0;
				current_foil++;
			}
			
			
		}
	}
	
}
//_____________________________________________________________________________________________
void MainWin::FillPoolTab()
{
	int Nline=7-WIDE;
	int Ncol=6+WIDE;
	
	fTabFoilPool= new TGCompositeFrame**[fNumberOfParc];
	fCheckArrayPool= new TGCheckButton**[fNumberOfParc];
	int *NTab=new int[fNumberOfParc];
	
	for (int i=0 ; i<fNumberOfParc;i++)
	{
		NTab[i]=fNumberOfPool[i]/(Nline*Ncol)+1;
		fTabFoilPool[i]= new TGCompositeFrame*[NTab[i]];
		fCheckArrayPool[i]= new TGCheckButton*[fNumberOfPool[i]];
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		string TabName[NTab[p]];
		for(int l=0; l<NTab[p] ; l++)
		{
			stringstream tmp;
			tmp.str("");
			tmp<<"Pool "<<l;
			TabName[l]=tmp.str();
			fTabFoilPool[p][l]=fItemTab[p][3]->AddTab(TabName[l].c_str());
			fTabFoilPool[p][l]->SetLayoutManager(new TGMatrixLayout(fTabFoilPool[p][l], Nline, 0, 5));
		}
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		int current_foil=0;
		int current_item=0;
		int current_item_in_the_foil=0;
		for(int n=0;n<fNumberOfPool[p];n++)
		{
			
			fCheckArrayPool[p][current_item]= new TGCheckButton(fTabFoilPool[p][current_foil],fDATA->GetPoolName()[p][n]);
			fCheckArrayPool[p][current_item]->SetFont(fLabelFontS);
			
			fTabFoilPool[p][current_foil]->AddFrame(fCheckArrayPool[p][current_item]);
			fCheckArrayPool[p][current_item]->Associate(this);
			
			current_item_in_the_foil++;
			current_item++;
			if(current_item_in_the_foil>=Ncol*Nline)
			{
				current_item_in_the_foil=0;
				current_foil++;
			}
			
			
		}
	}
	
	
	
}
//_____________________________________________________________________________________________
void MainWin::FillFabricationTab()
{
	int Nline=7-WIDE;
	int Ncol=6+WIDE;
	
	fTabFoilFab= new TGCompositeFrame**[fNumberOfParc];
	fCheckArrayFab= new TGCheckButton**[fNumberOfParc];
	int *NTab=new int[fNumberOfParc];
	
	for (int i=0 ; i<fNumberOfParc;i++)
	{
		NTab[i]=fNumberOfFab[i]/(Nline*Ncol)+1;
		fTabFoilFab[i]= new TGCompositeFrame*[NTab[i]];
		fCheckArrayFab[i]= new TGCheckButton*[fNumberOfFab[i]];
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		string TabName[NTab[p]];
		for(int l=0; l<NTab[p] ; l++)
		{
			stringstream tmp;
			tmp.str("");
			tmp<<"Fab "<<l;
			TabName[l]=tmp.str();
			fTabFoilFab[p][l]=fItemTab[p][4]->AddTab(TabName[l].c_str());
			fTabFoilFab[p][l]->SetLayoutManager(new TGMatrixLayout(fTabFoilFab[p][l], Nline, 0, 5));
		}
		
	}
	
	for (int p=0;p<fNumberOfParc;p++)
	{
		
		int current_foil=0;
		int current_item=0;
		int current_item_in_the_foil=0;
		for(int n=0;n<fNumberOfFab[p];n++)
		{
			
			fCheckArrayFab[p][current_item]= new TGCheckButton(fTabFoilFab[p][current_foil],fDATA->GetFabricationName()[p][n]);
			fCheckArrayFab[p][current_item]->SetFont(fLabelFontS);
			
			fTabFoilFab[p][current_foil]->AddFrame(fCheckArrayFab[p][current_item]);
			fCheckArrayFab[p][current_item]->Associate(this);
			
			current_item_in_the_foil++;
			current_item++;
			if(current_item_in_the_foil>=Ncol*Nline)
			{
				current_item_in_the_foil=0;
				current_foil++;
			}
			
			
		}
	}
	
}
//_____________________________________________________________________________________________
SubWin::SubWin(string event, const TGWindow *p, const TGWindow *main, UInt_t w,UInt_t h,UInt_t options ):
TGTransientFrame(p, main, w, h, options)
{
	fParent=(MainWin*)main;
	
	fS0 = new TGCompositeFrame(this, 50, 20, kVerticalFrame);
	fL0= new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX  , 2, 2, 2, 2);
	AddFrame(fS0, fL0);
	
	fSH1=new TGHorizontalFrame(fS0, 200, 20);
	fS0->AddFrame(fSH1,fL0);
	fLmsg=0;
	TEName=0;
	
	if(event=="SaveAs") SaveAs();
	// Frame for Ok/Cancel button
	//
	fSH2=new TGHorizontalFrame(this, 100, 20,kFixedWidth);
	AddFrame(fSH2,new TGLayoutHints(kLHintsBottom | kLHintsRight,0,0,10,0));
	
	fButtonOK = new TGTextButton(fSH2, "&Ok",M_BUT_OK);
	fButtonCan= new TGTextButton(fSH2,"&Cancel",M_BUT_CANCEL);
	fSH2->AddFrame(fButtonOK, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2, 2, 2, 2));
	fSH2->AddFrame(fButtonCan, new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,2, 2, 2, 2));
	fButtonOK->Associate(this);
	fButtonCan->Associate(this);
	
	
	//***************************************************************************
	//	The REMAP Section
	//***************************************************************************
	MapSubwindows();
	
	Resize(GetDefaultSize()); //fit to the exact size
	// position relative to the parent's window
	Window_t wdum;
	int ax, ay;
	gVirtualX->TranslateCoordinates(main->GetId(), GetParent()->GetId(),
					((TGFrame *) main)->GetWidth()-(int)(1.1*fWidth),
					(((TGFrame *) main)->GetHeight() - (int)(1.3*fHeight)) ,ax, ay, wdum);
	Move(ax, ay);
	MapSubwindows();
	MapWindow();
	fClient->WaitFor(this);
	
}
//_____________________________________________________________________________________________
void SubWin::SaveAs()
{
 fLmsg=new TGLabel(fSH1, new TGString("File Name:"));
 fSH1->AddFrame(fLmsg, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
 TGTextBuffer *TBName = new TGTextBuffer(100);
 TBName->AddText(0, "");
 TEName = new TGTextEntry(fSH1, TBName);
 TEName->Resize(200, TEName->GetDefaultHeight());
 fSH1->AddFrame(TEName, fL0);
 TEName->Associate(this);
	
 fRadioASCIISave=new TGRadioButton(fSH1,"ASCII",M_RADIO_ASCII_SAVE);
 fRadioXMLSave=new TGRadioButton(fSH1,"XML",M_RADIO_XML_SAVE);
 fRadioASCIISave->SetState(kButtonDown);
 fRadioASCIISave->Associate(this);
 fRadioXMLSave->Associate(this);
 fSH1->AddFrame(fRadioASCIISave);
 fSH1->AddFrame(fRadioXMLSave);
 
 SetWindowName("Saving Plotted Data");
 
}
//_____________________________________________________________________________________________
SubWin::~SubWin()
{
	delete fButtonOK ;
	delete fButtonCan;
	delete fSH2;
	if(TEName) delete TEName;
	if(fLmsg)
		delete fLmsg;
	delete fSH1;
	delete fS0 ;
	delete fL0 ;
}

//_____________________________________________________________________________________________
Bool_t SubWin::ProcessMessage(Long_t msg, Long_t parm1, Long_t parm2)
{
	switch (GET_MSG(msg))
	{
		case kC_TEXTENTRY:
			switch (GET_SUBMSG(msg))
		{
			case kTE_TEXTCHANGED:
				
				break;
		}
			break;
		case kC_COMMAND:
			switch (GET_SUBMSG(msg))
		{
			case kCM_RADIOBUTTON:
				switch (parm1)
			{
				case M_RADIO_ASCII_SAVE:
					fRadioASCIISave->SetState(kButtonDown);
					fRadioXMLSave->SetState(kButtonUp);
					fParent->fSaveFileFormat="ASCII";
					break;
				case M_RADIO_XML_SAVE:
					fRadioASCIISave->SetState(kButtonUp);
					fRadioXMLSave->SetState(kButtonDown);
					fParent->fSaveFileFormat="XML";
					break;
			}
			case kCM_BUTTON:
				switch (parm1)
			{
				case M_BUT_CANCEL:
					CloseWindow();
					break;
				case M_BUT_OK:
					if(TEName)
						fParent->fSaveFileName=TEName->GetBuffer()->GetString();
					CloseWindow();
					break;
					
			}
				break;
			default	:break;
		}
	}
	return kTRUE;
}


