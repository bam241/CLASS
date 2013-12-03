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

MainWin::MainWin(CLASSRead * DATA)
{
	fDATA=DATA;
	Start();
}


MainWin::~MainWin()
{
	delete fDATA;
}


//_____________________________________________________________________________________________
void MainWin::Start()
{
	
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
	//Compteur d'objet selectionn�
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
	
	/*Compteurs d'objets selectionn� dans chaque parc****/
	//exemple fNselectedReac[i][j]=1 -> alors le reacteur numero j du parc i a �t� sel�ctionn�
	// les valeurs de ces matrices sont soit 0 soit 1
	
	/*	fNselectedTOT=new int*[NumberOfParc];
	 fNselectedReac=new int*[NumberOfParc];
	 fNselectedStock=new int*[NumberOfParc];
	 fNselectedPool=new int*[NumberOfParc];
	 fNselectedFab=new int*[NumberOfParc];
	 
	 for (int i=0;i<NumberOfParc;i++)
	 {	for (int j=0;j<NumberOfReactor[i];j++)
	 fNselectedReac[i][j]=0;
	 for (int j=0;j<NumberOfStock[i];j++)
	 fNselectedStock[i][j]=0;
	 for (int j=0;j<NumberOfReactor[i];j++)
	 NumberOfPool[i][j]=0;
	 for (int j=0;j<NumberOfFab[i];j++)
	 fNselectedFab[i][j]=0;
	 for (int j=0;j<NumberOfTOT;j++)
	 fNselectedTOT[i][j]=0;
	 }
	 
	 */
	//....
	//
	//  The canvas (plotting window)
	//
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
	TGLayoutHints* fL55100 = new TGLayoutHints(kLHintsTop | kLHintsLeft | kLHintsExpandX,5, 5, 10, 0);
	//
	// the Main frame where all will be insert
	//
	
	fGeneF0 = new TGVerticalFrame(this, 200, 150);
	AddFrame(fGeneF0, fL5555);
	SetHeight(220);
	SetWidth(160);
	int MAXPATHLEN=256;
	char cDir[MAXPATHLEN];
	getcwd(cDir, MAXPATHLEN);
	stringstream tmp;
	tmp.str("");
	tmp << "CLASSGui " << "@@@NOM_Du_SCENAR";
	this->SetWindowName(tmp.str().c_str());
	
	
	/***************** LES Differents Parcs  ***********************/
	
	fParcTab = new TGTab(fGeneF0,200,150);
	fGeneF0->AddFrame(fParcTab,fL2222);
	fParcTab->Associate(this);
	
	fParcTabFoil = new TGCompositeFrame*[fNumberOfParc];
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fParcTabFoil[i] = new TGCompositeFrame;
		string ParcName = "Park ";
		fParcTabFoil[i] = fParcTab->AddTab(ParcName.c_str());
	}
	
	
	/*************Les differentes Facilities Tab****/
	fFacilitiesTab = new TGTab*[fNumberOfParc];
	
	for(int i = 0; i < fNumberOfParc; i++)
	{
		fFacilitiesTab[i] = new TGTab(fParcTabFoil[i],200,150);
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
			fItemTab[i][j] = new TGTab(fFacilitiesTabFoil[i][j],200,150);//@@ hxx non d�fini
			fFacilitiesTabFoil[i][j]->AddFrame(fItemTab[i][j]);
			fItemTab[i][j]->Associate(this);
			//construire les foils en fonction du nombre de item possible par foil
			FillItemTab(j);
		}
		
	}
	
	
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
	fPlotSaveQuitFrame = new TGHorizontalFrame(fGeneF0, 400, 50 );
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
	//------------------------------------------------------------------------
	//The Nucleus Foil
	//------------------------------------------------------------------------
	
	//
	//	The Miscellaneous part (AM check box, ...)
	//
	/*	fMiscGrpFNuc = new TGGroupFrame(fGeneF0,"Miscellaneous", kVerticalFrame );
	 ((TGGroupFrame*)fMiscGrpFNuc)->SetTextFont(fLabelFontB);
	 fGeneF0->AddFrame(fMiscGrpFNuc,fL2222);
	 
	 fMiscHzFrame1=new TGHorizontalFrame(fMiscGrpFNuc,400, 50 );
	 fMiscGrpFNuc->AddFrame(fMiscHzFrame1,fL2222);
	 //misc check buttons
	 fCheckSumNuc=new TGCheckButton(fMiscHzFrame1,"Sum of Selected",M_CHECK_SUM_NUC);
	 fCheckSumNuc->SetFont(fLabelFontS);
	 fMiscHzFrame1->AddFrame(fCheckSumNuc,fL2222);
	 fCheckSumNuc->Associate(this);
	 */
	/***En mol ou en masse?**/
	/*	fRadioMol=new TGRadioButton(fMiscHzFrame1,"Mol",M_RADIO_MOL);
	 fRadioMol->SetFont(fLabelFontS);
	 fMiscHzFrame1->AddFrame(fRadioMol,fL2222);
	 fRadioMol->Associate(this);
	 
	 fRadioMass=new TGRadioButton(fMiscHzFrame1, "Mass",M_RADIO_MASS);
	 fRadioMass->SetFont(fLabelFontS);
	 fMiscHzFrame1->AddFrame(fRadioMass,fL2222);
	 fRadioMass->Associate(this);
	 
	 fRadioMol->SetState(kButtonDown);
	 */
	/*********Choix du vecteur Temps************/
	
	/*	fTimeScale=new TGComboBox(fMiscHzFrame1,M_CB_TIMESCALE_NUC);
	 fMiscHzFrame1->AddFrame(fTimeScale,fL2222);
	 fTimeScale->Resize(80,17);
	 fTimeScale->AddEntry("years",1);
	 fTimeScale->AddEntry("months",2);
	 fTimeScale->AddEntry("days",3);
	 fTimeScale->AddEntry("seconds",4);
	 fTimeScale->Select(1);
	 fTimeScale->Associate(this);
	 fAxisHzFrame->Resize(350, 40);
	 
	 
	 //
	 // apply the correct unit to all X axis combo box
	 //
	 fOldTimeUnits=4;
	 */
	//
	// Final part: this allows the effective display of widgets
	//
	
	MapSubwindows();						// Map all sub windows that are part of the composite frame
	MapWindow();
	Resize(GetDefaultSize()); 					// fit to the exact size
	Move(410,30);
	fMainWidth=fGeneF0->GetWidth();
	
	Resize(550,750  ); 					// fit to the exact size

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
					if(fSaveFileName != "")
					{	//si plusieur canvas ouvert sauver l'actif exemple :
						//	string PadName=gPad->GetName();
						//	if(PadName=="c_Gam")
						//	{
						//		Save(0);
						//	}
						
					}
					break;
			}
				
				break;
				/*	case kCM_COMBOBOX:
				 switch (parm1)
				 {
				 case M_CB_TIMESCALE:
				 int code=fTimeScale->GetSelected();
				 cout<<"user want time scale code :"<<code<<endl;
				 break;
				 }
				 break;
				 */
				/*				case kCM_RADIOBUTTON:
				 switch (parm1)
				 {
				 case M_RADIO_MOL:
				 fRadioMol->SetState(kButtonDown);
				 fRadioMass->SetState(kButtonUp);
				 cout<<"user want results in mol"<<endl;
				 break;
				 case M_RADIO_MASS:
				 fRadioMol->SetState(kButtonUp);
				 fRadioMass->SetState(kButtonDown);
				 cout<<"user want results in kg"<<endl;
				 break;
				 }
				 break;
				 case kCM_CHECKBUTTON:
				 switch (parm1)
				 {
				 default:
				 //Noyaux
				 if(parm1>M_CHECK_PLOTALL)
				 {	int l=parm1-M_CHECK_PLOTALL-1;
				 if(fNAddedCheck[l]>1)
				 {
				 //bool shit: a "same" button can be in more than one category
				 //in fact it is not the same button thus find all the correponding button
				 //and change all of them
				 EButtonState WantedState;
				 bool IsMainBut=true;
				 if(parm2==Long_t(fCheckArrayNuc[l]->GetUserData()))
				 {
				 WantedState=fCheckArrayNuc[l]->GetState();
				 }
				 else
				 {
				 IsMainBut=false;
				 if(fCheckArrayNuc[l]->GetState()==kButtonDown)
				 WantedState=kButtonUp;
				 else
				 WantedState=kButtonDown;
				 }
				 for(int i=0; i<int(fSaveCheckButton.size());i++)
				 {
				 if(int(parm1)==fSaveCheckButton[i]->WidgetId())
				 fSaveCheckButton[i]->SetState(WantedState);
				 }
				 fCheckArrayNuc[l]->SetState(WantedState);
				 }
				 if(fCheckArrayNuc[l]->GetState()==kButtonDown)
				 fNselectedNucleus++;
				 if(fCheckArrayNuc[l]->GetState()==kButtonUp)
				 fNselectedNucleus--;
				 }
				 //total
				 if(parm1>M_CHECK_PLOTALL && fFacilitiesTab->GetCurrent()==0)
				 {
				 int l=parm1-M_CHECK_PLOTALL-1-1000;
				 if(fCheckArrayTOT[l]->GetState()==kButtonDown)
				 fNselectedTOT[fParcTab->GetCurrent()]++;
				 if(fCheckArrayTOT[l]->GetState()==kButtonUp)
				 fNselectedTOT[fParcTab->GetCurrent()]--;
				 }
				 //Reacteurs
				 if(parm1>M_CHECK_PLOTALL && fFacilitiesTab->GetCurrent()==1)
				 {
				 int l=parm1-M_CHECK_PLOTALL-1-5000;
				 if(fCheckArrayReac[l]->GetState()==kButtonDown)
				 fNselectedReac[fParcTab->GetCurrent()]++;
				 if(fCheckArrayReac[l]->GetState()==kButtonUp)
				 fNselectedReac[fParcTab->GetCurrent()]--;
				 }
				 //Stock
				 if(parm1>M_CHECK_PLOTALL && fFacilitiesTab->GetCurrent()==2)
				 {
				 int l=parm1-M_CHECK_PLOTALL-1-9000;
				 if(fCheckArrayStock[l]->GetState()==kButtonDown)
				 fNselectedStock[fParcTab->GetCurrent()]++;
				 if(fCheckArrayStock[l]->GetState()==kButtonUp)
				 fNselectedStock[fParcTab->GetCurrent()]--;
				 }
				 //Pooling
				 if(parm1>M_CHECK_PLOTALL && fFacilitiesTab->GetCurrent()==3)
				 {
				 int l=parm1-M_CHECK_PLOTALL-1-13000;
				 if(fCheckArrayPool[l]->GetState()==kButtonDown)
				 fNselectedPool[fParcTab->GetCurrent()]++;
				 if(fCheckArrayPool[l]->GetState()==kButtonUp)
				 fNselectedPool[fParcTab->GetCurrent()]--;
				 }
				 //Fabrication Plant
				 if(parm1>M_CHECK_PLOTALL && fFacilitiesTab->GetCurrent()==4)
				 {
				 int l=parm1-M_CHECK_PLOTALL-1-17000;
				 if(fCheckArrayFab[l]->GetState()==kButtonDown)
				 fNselectedFab[fParcTab->GetCurrent()]++;
				 if(fCheckArrayFab[l]->GetState()==kButtonUp)
				 fNselectedFab[fParcTab->GetCurrent()]--;
				 }
				 
				 }*/
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

	//Power
	for(int i=0; i < fNumberOfParc; i++)
	{
		if(fCheckArrayTotal[i][fNumberOfTOT-1]->GetState()==kButtonDown)
			toplotPower.push_back( CLASSPlotElement(i, -2, -2, -2,-2,-2) );

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
						toplot.push_back( CLASSPlotElement(i, 0, j, fDATA->GetZAIvector()[k]));
				}
		}
		
		for(int j=0; j < fNumberOfReactor[i]; j++)
		{
			if(fCheckArrayReactor[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 1, j, fDATA->GetZAIvector()[k]));
				}

		}
		
		for(int j=0; j < fNumberOfStock[i]; j++)
		{
			if(fCheckArrayStock[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 2, j, fDATA->GetZAIvector()[k]));
				}
		}
		
		for(int j=0; j < fNumberOfPool[i]; j++)
		{
			if(fCheckArrayPool[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 3, j, fDATA->GetZAIvector()[k]));
				}
		}
		
		for(int j=0; j < fNumberOfFab[i]; j++)
		{
			if(fCheckArrayFab[i][j]->GetState()==kButtonDown)
				for(int k=0; k < Nnucleus; k++)
				{
					if(fCheckArrayNuc[k]->GetState()==kButtonDown)
						toplot.push_back( CLASSPlotElement(i, 4, j, fDATA->GetZAIvector()[k]));
				}
		}
	}
	if(toplot.size() != 0)
		fDATA->Plot(toplot);
	
}
//_____________________________________________________________________________________________
void MainWin::FillNucTab()		// fill the Inventory Tab foil
{
	static bool first=true;
	
	int NCheck=fDATA->GetZAIvector().size();
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
	//	fTabNuc->SetTab(TabName[1].c_str()) ;
	//	fTabNuc->SetTab(TabName[0].c_str()) ;
	first=false;
	
}

//_____________________________________________________________________________________________
void MainWin::FillNucFoil(int n_item, int Ncol,int Nline)
{
	int current_foil = 0;
	int current_item = 0;
	for( int l = 0; l < n_item; l++ )
	{
		
		stringstream name;
		int Atmp = fDATA->GetZAIvector()[l].A();
		
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
			
		} else
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
			if(i==3) name="GOD";
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
	int Nline=15-WIDE;
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
			tmp.str("");
			//tmp<<"Reac "<<l;
			TabName[l]=tmp.str();
			//cout<<l<<" "<<TabName[l]<<endl;
			fTabFoilReactor[p][l]=fItemTab[p][1]->AddTab(TabName[l].c_str());
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
	int Nline=15-WIDE;
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
	int Nline=15-WIDE;
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
	int Nline=15-WIDE;
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
	fSH2=new TGHorizontalFrame(this, 200, 20,kFixedWidth);
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
/*void SubWin::SaveAs()
 {
 fLmsg=new TGLabel(fSH1, new TGString("File Name:"));
 fSH1->AddFrame(fLmsg, new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2));
 TGTextBuffer *TBName = new TGTextBuffer(100);
 TBName->AddText(0, "");
 TEName = new TGTextEntry(fSH1, TBName);
 TEName->Resize(200, TEName->GetDefaultHeight());
 fSH1->AddFrame(TEName, fL0);
 TEName->Associate(this);
 
 SetWindowName("Saving Plotted Data");
 
 }*/
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

