#include <fstream>
#include <sstream>
#include <cmath>
#include "fenetre.hxx"
using namespace std;


////////////////////////////////////////////////////////////////
//                       LitmVue




LitmVue::LitmVue() {

//    fMin=fMax=0;
  fNext=0;
  fFile=0;
  fNoFile=0;
  fTally=0;
  fHisto=0;
  fHisto2=0;
  fNomFichier="NO FILE";
  fNoTally=0;
  fType=0;
  fnCbins=fnEbins=fnTbins=0;
  fCbins=fEbins=fTbins=0;
  zeroE=zeroT=0;
  
}

LitmVue::LitmVue(LitmVue& vue) {
  
//    fMin=vue.fMin;
//    fMax=vue.fMax;
  fNext=0;
  fFile=vue.fFile;
  fNoFile=vue.fNoFile;
  fNomFichier=vue.fNomFichier;
  fTally=vue.fTally;
  fNoTally=vue.fNoTally;
  fType=vue.fType;
  fL=vue.fL;
  fD=vue.fD;
  fU=vue.fU;
  fS=vue.fS;
  fM=vue.fM;
  fC=vue.fC;
  fnCbins=vue.fnCbins;
  if (fnCbins>0) {
    fCbins=new Double_t[fnCbins+1];
    for (unsigned long i=0;i<=fnCbins;i++)
      fCbins[i]=vue.fCbins[i];
  }
  else
    fCbins=0;

  fE=vue.fE;
  fnEbins=vue.fnEbins;
  if (fnEbins>0) {
    fEbins=new Double_t[fnEbins+1];
    for (unsigned long i=0;i<=fnEbins;i++)
      fEbins[i]=vue.fEbins[i];
  }
  else
    fEbins=0;
  zeroE=vue.zeroE;
  fT=vue.fT;
  fnTbins=vue.fnTbins;
  if (fnTbins>0) {
    fTbins=new Double_t[fnTbins+1];
    for (unsigned long i=0;i<=fnTbins;i++)
      fTbins[i]=vue.fTbins[i];
  }
  else
    fTbins=0;
  zeroT=vue.zeroT;
  
  if (!vue.fHisto)
    fHisto=0;
  else
    fHisto=new THistoCalc(*(vue.fHisto));

  if (!vue.fHisto2)
    fHisto2=0;
  else
    fHisto2=new THisto2Calc(*(vue.fHisto2));

}
LitmVue::LitmVue(LitmVue& vue, double MyConst) {
  
//    fMin=vue.fMin;
//    fMax=vue.fMax;
  fNext=0;
  fFile=vue.fFile;
  fNoFile=vue.fNoFile;
  fNomFichier=vue.fNomFichier;
  fTally=vue.fTally;
  fNoTally=vue.fNoTally;
  fType=vue.fType;
  fL=vue.fL;
  fD=vue.fD;
  fU=vue.fU;
  fS=vue.fS;
  fM=vue.fM;
  fC=vue.fC;
  fnCbins=vue.fnCbins;
  if (fnCbins>0) {
    fCbins=new Double_t[fnCbins+1];
    for (unsigned long i=0;i<=fnCbins;i++)
      fCbins[i]=vue.fCbins[i];
  }
  else
    fCbins=0;

  fE=vue.fE;
  fnEbins=vue.fnEbins;
  if (fnEbins>0) {
    fEbins=new Double_t[fnEbins+1];
    for (unsigned long i=0;i<=fnEbins;i++)
      fEbins[i]=vue.fEbins[i];
  }
  else
    fEbins=0;
  zeroE=vue.zeroE;
  fT=vue.fT;
  fnTbins=vue.fnTbins;
  if (fnTbins>0) {
    fTbins=new Double_t[fnTbins+1];
    for (unsigned long i=0;i<=fnTbins;i++)
      fTbins[i]=vue.fTbins[i];
  }
  else
    fTbins=0;
  zeroT=vue.zeroT;
  
  if (!vue.fHisto)
    fHisto=0;
  else
  {
    fHisto=new THistoCalc(*(vue.fHisto));
	ValErr_t MyValErr(MyConst,0.);
	(*fHisto)*=MyValErr;
  }
  if (!vue.fHisto2)
    fHisto2=0;
  else
    fHisto2=new THisto2Calc(*(vue.fHisto2));

}

LitmVue::~LitmVue() {

  delete fHisto;
  delete fHisto2;
  delete fNext;
  delete fCbins;
  delete fEbins;
  delete fTbins;


}

int LitmVue::HVarie() {
  if (!fFile)
    return 0;
  return (fFile->GetNbTal()>1);
}

int LitmVue::LVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fLN>1 && (fType>=10));
}

int LitmVue::DVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fDN>1 && (fType<10 || fType>=20));
}

int LitmVue::UVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fUN>1 && (fType<20 || fType>=30));
}

int LitmVue::SVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fSN>1 && (fType<30 || fType>=40));
}

int LitmVue::MVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fMN>1 && (fType<40 || fType>=50));
}

int LitmVue::CVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fCN>1 && (fType<50 || fType>=60));
}

int LitmVue::EVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fEN>1 && (fType<60 || fType>=70) && fType<100);
}

int LitmVue::TVarie() {
  if (!fFile || !fTally)
    return 0;
  return (fTally->fTN>1 && (fType<70 || fType>=80) && fType<100);
}

// ecrit sur une ligne tout ce qui definit l'histo !
// utilise cout
void LitmVue::Print(int verbose,ostream& out) {
  
  if (!fFile) {
    out << "Robert";
    if (verbose)
      out << endl;
    return;
  }
  out << fNomFichier;
  if (HVarie())
  {
    if (verbose)
      out << ", tally " << fTally->fNo;
    else
      out << "_f" << fTally->fNo;
  }
  if (LVarie())
  {
    if (verbose)
      out << ", l=" << fL << '(' << fTally->fLV[fL] << ")/" << fTally->fLN;
    else
      out << "_l" << fL;
  }
  if (DVarie())
  {
    if (verbose)
      out << ", d=" << fD << "/" << fTally->fDN;
    else
      out << "_d" << fD;
  }
  if (UVarie()) {
    if (fTally->fUT && fU==fTally->fUN-1)
    {
	  if (verbose)
		out << ", u=t";
      else
		out << "_ut";
    }
	else
	{
	  if (verbose)
		out << ", u=" << fU << "/" << fTally->fUN;
      else
		out << "_u" << fU;
	}
  }
  if (SVarie()) {
    if (fTally->fST && fS==fTally->fSN-1)
	{
      if (verbose)
		out << ", s=t";
      else
		out << "_st";
	}
    else
	{
      if (verbose)
		out << ", s=" << fS << "/" << fTally->fSN;
      else
		out << "_s" << fS;
	}
  }
  if (MVarie()) {
    if (fTally->fMT && fM==fTally->fMN-1)
      if (verbose)
	out << ", m=t";
      else
	out << "_mt";
    else
      if (verbose)
	out << ", m=" << fM << "/" << fTally->fMN;
      else
	out << "_m" << fM;
  }
  /**/ // toujours pas de traitement des bins centres !!
  if (CVarie()) {
    if (fTally->fCT && fC==fTally->fCN-1)
	{
      if (verbose)
		out << ", c=t";
      else
		out << "_ct";
    }
	else if (fCbins)
    {
	  if (verbose)
		out << ", c=[" << fCbins[fC] << ',' << fCbins[fC+1] << "]/["
	    << fCbins[0] << ',' << fCbins[fnCbins] << ']';
      else
		out << "_c" << fC;
	}
  }
  if (EVarie()) {
    if (fTally->fET && fE+zeroE==fTally->fEN-1)
	{
      if (verbose)
		out << ", e=t";
      else
		out <<  "_et";
	}
    else if (fEbins)
	{
      if (verbose)
		out << ", e=[" << fEbins[fE] << ',' << fEbins[fE+1] << "]/["
	    << fEbins[0] << ',' << fEbins[fnEbins] << ']';
      else
		out << "_e" << fE;
	}
  }
  if (TVarie()) {
    if (fTally->fTT && fT+zeroT==fTally->fTN-1)
	{
      if (verbose)
		out << ", t=t";
      else
		out << "_tt";
	}
    else if (fTbins)
	{
      if (verbose)
		out << ", t=[" << fTbins[fT] << ',' << fTbins[fT+1] << "]/["
	    << fTbins[0] << ',' << fTbins[fnTbins] << ']';
      else
		out << "_t" << fT;
	}
  }
  if(fHisto)
  	out<<" (histo name: "<<fHisto->GetName()<<")";
  if (verbose)
    out << endl;
  
}

void LitmVue::PrintAll(int n,unsigned type) {

  if(fNext)
    fNext->PrintAll(n+1,type);
  else
    cout << endl;

  if (fType==type)
    cout << '*';
  else
    cout << ' ';

  //cout.form("%3d : ",n);
  cout<<" : "<<n;
  Print();
}

void LitmVue::Draw() {

  if(fType<100 && fHisto)
    fHisto->Draw();
  if(fType>=100 && fHisto2)
    fHisto2->Draw("COLZ");

}

void LitmVue::DrawSame() {

  if(fType<100 && fHisto)
    fHisto->Draw("SAME");  
  if(fType>=100 && fHisto2)
    fHisto2->Draw("SAME");

}

void LitmVue::DrawAll(unsigned n,unsigned type) {
  
  if (fType==type) {
    if (n==1)
      Draw();
    else
      DrawSame();
  }

  if(fNext)
    fNext->DrawAll(n+1,type);
  
}

void LitmVue::ChangeFile(TMctal *fic, char *nom, unsigned long numero) {

  fFile=fic;
  fNoFile=numero;
  fNomFichier=nom;
  fNoTally=0;
  fTally=fic->GetTally(0);
  if (!fTally) {
    fic=0;
    fNomFichier="NO FILE";
  }
  else {
    fType=70;
    ChangeHisto();
  }
}

unsigned LitmVue::SetType(unsigned type) {

  if (!fFile || !fTally)
    return 0;
  if (fType==type)
    return 1;

  fType=type;
  ChangeHisto();
  return 1;

}

void LitmVue::ChangeHisto() {
  
  zeroE=0;
  zeroT=0;

  if (!fFile)
    return;
  
  delete fHisto;
  fHisto=0;
  delete fHisto2;
  fHisto2=0;
  delete fCbins;
  delete [] fEbins;
  delete [] fTbins;
  fL=fD=fU=fS=fM=fC=fE=fT=0;

  // initialisation des binnings a partir du tally.

  fnCbins=fTally->fCN-fTally->fCT;
  if (fnCbins>0) {
    fCbins=new Double_t[fnCbins+1];
    fCbins[0]=0;
    for (unsigned long i=0;i<fnCbins;i++)
      fCbins[i+1]=fTally->fCV[i];
  }
  else
    fCbins=0;

  fnEbins=fTally->fEN-fTally->fET;
  if (fnEbins>0 && fTally->fEV[0]==0) {
    zeroE=1;
    fnEbins--;
  }
  if (fnEbins>0) {
    fEbins=new Double_t[fnEbins+1];
    fEbins[0]=0;
    for (unsigned long i=0;i<fnEbins;i++)
      fEbins[i+1]=fTally->fEV[i+zeroE];
  }
  else
    fEbins=0;

  fnTbins=fTally->fTN-fTally->fTT;
  if (fnTbins>0 && fTally->fTV[0]==0) {
    zeroT=1;
    fnTbins--;
  }
  if (fnTbins>0) {
    fTbins=new Double_t[fnTbins+1];
    fTbins[0]=0;
    for (unsigned long i=0;i<fnTbins;i++)
      fTbins[i+1]=fTally->fTV[i+zeroT]/100.0; // Passage en microsecondes !
  }
  else
    fTbins=0;
  

  // creation histo
  switch (fType) {

  case 0:
    fHisto=new THistoCalc(fTally->fLN,0,fTally->fLN);
    fHisto->SetStats(0);
    fHisto->SetXTitle("Location");
    break;

  case 70:
    fHisto=new THistoCalc(fnTbins,fTbins);
    fHisto->SetStats(0);
    fHisto->SetXTitle("t (#mu s)");
    break;
    
  case 60:
    fHisto=new THistoCalc(fnEbins,fEbins);
    fHisto->SetStats(0);
    fHisto->SetXTitle("E (MeV)");
    break;

  case 100:
    fHisto2=new THisto2Calc(fnTbins,fTbins,fnEbins,fEbins);
    fHisto2->SetStats(0);
    fHisto2->SetXTitle("t (#mu s)");
    fHisto2->SetYTitle("E (MeV)");
    break;
  }

  FillHisto();

}

void LitmVue::FillHisto() {

  if (!fFile || !fTally)
    return;

  ValErr_t val;
  switch (fType) {
  case 0:
    if (!fHisto)
      return;
    for (unsigned long l=0;l<fTally->fLN;l++) {
      val=fTally->fVal[l][fD][fU][fS][fM][fC][fE+zeroE][fT+zeroT];
      if (fCbins && fC<fnCbins)
	val/=fCbins[fC+1]-fCbins[fC];
      if (fTbins && fT<fnTbins)
	val/=fTbins[fT+1]-fTbins[fT];
      if (fEbins && fE<fnEbins)
	val/=log(fEbins[fE+1])-log(fEbins[fE]);
      fHisto->SetBinValErr(l+1,val);
    }
    break;

  case 70:
    if (!fHisto)
      return;
    for (unsigned long t=0;t<fnTbins;t++) {
      val=fTally->fVal[fL][fD][fU][fS][fM][fC][fE+zeroE][t+zeroT];
      if (fCbins && fC<fnCbins)
	val/=fCbins[fC+1]-fCbins[fC];
      if (fTbins && t<fnTbins)
	val/=fTbins[t+1]-fTbins[t];
      if (fEbins && fE<fnEbins)
	val/=log(fEbins[fE+1])-log(fEbins[fE]);
      fHisto->SetBinValErr(t+1,val);
    }
    break;
    
  case 60:
    if (!fHisto)
      return;
    for (unsigned long e=0;e<fnEbins;e++) {
      val=fTally->fVal[fL][fD][fU][fS][fM][fC][e+zeroE][fT+zeroT];
      if (fCbins && fC<fnCbins)
	val/=fCbins[fC+1]-fCbins[fC];
      if (fTbins && fT<fnTbins)
	val/=fTbins[fT+1]-fTbins[fT];
      if (fEbins && e<fnEbins)
	val/=log(fEbins[e+1])-log(fEbins[e]);
      fHisto->SetBinValErr(e+1,val);
    }
	if (fEbins)
	{
		cout<<"Integrale of "<<fHisto->GetName()<<"="<<fHisto->Integrale(fEbins[1],fEbins[fnEbins-1]).fVal<<endl;
		cout<<"Log Integrale of "<<fHisto->GetName()<<"="<<fHisto->LogIntegrale(fEbins[1],fEbins[fnEbins-1]).fVal<<endl;
    }
	break;
    
  case 100:
    if (!fHisto2)
      return;
    for (unsigned long t=0;t<fnTbins;t++)
      for (unsigned long e=0;e<fnEbins;e++) {
	val=fTally->fVal[fL][fD][fU][fS][fM][fC][e+zeroE][t+zeroT];
	if (fCbins && fC<fnCbins)
	  val/=fCbins[fC+1]-fCbins[fC];
	if (fTbins && t<fnTbins)
	  val/=log(fTbins[t+1])-log(fTbins[t]);
	if (fEbins && e<fnEbins)
	  val/=log(fEbins[e+1])-log(fEbins[e]);
	fHisto2->SetBinValErr(t+1,e+1,val);
      }
  }

}

void LitmVue::HandleInput(char c) {

  switch(c) {

    // DEPLACEMENTS

  case 'l':
    if (!LVarie())
      break;
    if (fL>0)
      fL--;
    else
      fL=fTally->fLN-1;
    FillHisto();
    break;
  case 'L':
    if (!LVarie())
      break;
    fL++;
    if (fL>=fTally->fLN)
      fL=0;
    FillHisto();
    break;

  case 'd':
    if (!DVarie())
      break;
    if (fD>0)
      fD--;
    else
      fD=fTally->fDN-1;
    FillHisto();
    break;
  case 'D':
    if (!DVarie())
      break;
    fD++;
    if (fD>=fTally->fDN)
      fD=0;
    FillHisto();
    break;

  case 'u':
    if (!UVarie())
      break;
    if (fU>0)
      fU--;
    else
      fU=fTally->fUN-1;
    FillHisto();
    break;
  case 'U':
    if (!UVarie())
      break;
    fU++;
    if (fU>=fTally->fUN)
      fU=0;
    FillHisto();
    break;

  case 's':
    if (!SVarie())
      break;
    if (fS>0)
      fS--;
    else
      fS=fTally->fSN-1;
    FillHisto();
    break;
  case 'S':
    if (!SVarie())
      break;
    fS++;
    if (fS>=fTally->fSN)
      fS=0;
    FillHisto();
    break;

  case 'm':
    if (!MVarie())
      break;
    if (fM>0)
      fM--;
    else
      fM=fTally->fMN-1;
    FillHisto();
    break;
  case 'M':
    if (!MVarie())
      break;
    fM++;
    if (fM>=fTally->fMN)
      fM=0;
    FillHisto();
    break;

  case 'c':
    if (!CVarie())
      break;
    if (fC>0)
      fC--;
    else
      fC=fTally->fCN-1;
    FillHisto();
    break;
  case 'C':
    if (!CVarie())
      break;
    fC++;
    if (fC>=fTally->fCN)
      fC=0;
    FillHisto();
    break;

  case 'e':
    if (!EVarie())
      break;
    if (fE>0)
      fE--;
    else
      fE=fTally->fEN-1-zeroE;
    FillHisto();
    break;
  case 'E':
    if (!EVarie())
      break;
    fE++;
    if (fE>=fTally->fEN-zeroE)
      fE=0;
    FillHisto();
    break;

  case 't':
    if (!TVarie())
      break;
    if (fT>0)
      fT--;
    else
      fT=fTally->fTN-1-zeroT;
    FillHisto();
    break;
  case 'T':
    if (!TVarie())
      break;
    fT++;
    if (fT>=fTally->fTN-zeroT)
      fT=0;
    FillHisto();
    break;

  case 'h':
    if (!HVarie())
      break;
    if (fNoTally>0)
      fNoTally--;
    else
      fNoTally=fFile->GetNbTal()-1;
    fTally=fFile->GetTally(fNoTally);
    ChangeHisto();
    break;
  case 'H':
    if (!HVarie())
      break;
    fNoTally++;
    if (fNoTally>=fFile->GetNbTal())
      fNoTally=0;
    fTally=fFile->GetTally(fNoTally);
    ChangeHisto();
    break;

  case 'o':
    
    if (!fFile || !fTally)
      break;
    if (LVarie())
      fL=0;
    if (DVarie())
      fD=0;
    if (UVarie())
      fU=0;
    if (SVarie())
      fS=0;
    if (MVarie())
      fM=0;
    if (CVarie())
      fC=0;
    if (EVarie())
      fE=0;
    if (TVarie())
      fT=0;
    FillHisto();
    break;
 

    //  CHOIX DE MODE

    // Mode 0;
  case ';':
    SetType(0);
    break;
    // Mode 70;
  case 'y':
    SetType(70);
    break;
    // Mode 60;
  case 'r':
    SetType(60);
    break;
    // Mode 100;
  case 'b':
    SetType(100);
    break;
    
    // ECRITURE

  case 'W':
    {
      stringstream filename;
      Print(0,filename);
      filename << ".dat" << ends;
      ofstream Out(filename.str().c_str());
      if (fType<100 && fHisto)
	fHisto->Write(Out);
      if (fType>=100 && fHisto2)
	fHisto2->Write(Out);
      cout << endl << "Wrote " << filename.str() << endl;
	Out.close();
    }
  break;    
  
   
  }

}

LitmVue& LitmVue::operator *= (LitmVue& v) {
  
  if (v.fType!=fType)
    return *this;

  fNoFile=0;
  fFile=0;
  fNomFichier="";
  fTally=0;
  fNoTally=0;
  fType=0;
  fL=fD=fU=fS=fM=fC=fE=fT=0;
  fnCbins=fnEbins=fnTbins=0;
  delete fCbins;
  delete fEbins;
  delete fTbins;
  fCbins=fEbins=fTbins=0;

  if (fType<100 && fHisto!=0 && v.GetHisto()!=0)
    *fHisto*=*(v.GetHisto());

  return *this;

}

LitmVue& LitmVue::operator /= (LitmVue& v) {
  
  if (v.fType!=fType)
    return *this;

  fNoFile=0;
  fFile=0;
  fNomFichier="";
  fTally=0;
  fNoTally=0;
  fType=0;
  fL=fD=fU=fS=fM=fC=fE=fT=0;
  fnCbins=fnEbins=fnTbins=0;
  delete fCbins;
  delete fEbins;
  delete fTbins;
  fCbins=fEbins=fTbins=0;

  if (fType<100 && fHisto!=0 && v.GetHisto()!=0)
    *fHisto/=*(v.GetHisto());
  
  return *this;

}

LitmVue& LitmVue::operator += (LitmVue& v) {
  
  if (v.fType!=fType)
    return *this;
  
  fNoFile=0;
  fFile=0;
  fNomFichier="";
  fTally=0;
  fNoTally=0;
  fType=0;
  fL=fD=fU=fS=fM=fC=fE=fT=0;
  fnCbins=fnEbins=fnTbins=0;
  delete fCbins;
  delete fEbins;
  delete fTbins;
  fCbins=fEbins=fTbins=0;

  if (fType<100 && fHisto!=0 && v.GetHisto()!=0)
    *fHisto+=*(v.GetHisto());

  return *this;
}

LitmVue& LitmVue::operator -= (LitmVue& v) {
  
  if (v.fType!=fType)
    return *this;
  
  fNoFile=0;
  fFile=0;
  fNomFichier="";
  fTally=0;
  fNoTally=0;
  fType=0;
  fL=fD=fU=fS=fM=fC=fE=fT=0;
  fnCbins=fnEbins=fnTbins=0;
  delete fCbins;
  delete [] fEbins;
  delete [] fTbins;
  fCbins=fEbins=fTbins=0;

  if (fType<100 && fHisto!=0 && v.GetHisto()!=0)
    *fHisto-=*(v.GetHisto());

  return *this;
}







//////////////////////////////////////////////////////////////////
//                            LitmFenetre





void LitmFenetre::Init() {

  fStack=0;
//    fMin=fMax=0;
  SetLogx(0);
  SetLogy(1);
  SetFillColor(10);
}

void LitmFenetre::Free() {

  for (unsigned i=0;i<fFichiers.size();i++)
    delete fFichiers[i];
  fFichiers.clear();
  fNoms.clear();

  delete fStack;

}

unsigned LitmFenetre::LoadMctal(char *nomfichier) {

  TMctal *fic=new TMctal;
  ifstream ific(nomfichier);
  if (!fic->Read(ific)) {
    delete fic;
    return 0;
  }
  fFichiers.push_back(fic);
  fNoms.push_back(nomfichier);

  return 1;

}

void LitmFenetre::HandleInput(EEventType event, Int_t px, Int_t py) {

  LitmVue *pvue;
  unsigned long no;

  TCanvas::HandleInput(event,px,py);

  if (event==kKeyPress) {
    switch(px) {

//  GESTION DE PILE

      // Nouvelle vue.
    case 'n':
      if (fFichiers.size()==0 || fNoms.size()==0)
	break;
      pvue=new LitmVue();
      pvue->ChangeFile(fFichiers[0],fNoms[0],0);
      pvue->fNext=fStack;
      fStack=pvue;
      break;
      // Duplication de vue.
    case '\r':
      if (!fStack)
	break;
      pvue=new LitmVue(*fStack);
      pvue->fNext=fStack;
      fStack=pvue;
      break;
      // Suppression de vue.
    case '\b':
      if (!fStack)
	break;
      pvue=fStack->fNext;
      fStack->fNext=0;
      delete fStack;
      fStack=pvue;
      break;
      // Swappe les deux dernieres vues.
    case ' ':
      if (!fStack || !(fStack->fNext))
	break;
      pvue=fStack->fNext->fNext;
      fStack->fNext->fNext=fStack;
      fStack=fStack->fNext;
      fStack->fNext->fNext=pvue;
      break;

      // roule dans un sens
    case '>':
      if (!fStack)
	break;
      pvue=fStack;
      while(pvue->fNext!=0)
	pvue=pvue->fNext;
      pvue->fNext=fStack;
      fStack=fStack->fNext;
      pvue->fNext->fNext=0;
      break;

//  CHANGEMENT DE FICHIER

    case 'F':
      if (!fStack || fFichiers.size()<=1)
	break;
      no=fStack->fNoFile;
      no++;
      if (no >= fFichiers.size())
	no=0;
      fStack->ChangeFile(fFichiers[no],fNoms[no],no);
      break;

    case 'f':
      if (!fStack || fFichiers.size()<=1)
	break;
      no=fStack->fNoFile;
      if (no==0)
	no=fFichiers.size()-1;
      else no--;
      fStack->ChangeFile(fFichiers[no],fNoms[no],no);
      break;


//  OPERATIONS SUR LES HISTOS

     case 'p':
	 {
      if (!fStack)
	break;
      double cte;
	  cout<<"Enter the Multiplicative Constant: "<<flush; 
	  cin>>cte;
	  pvue=new LitmVue(*fStack,cte);
      pvue->fNext=fStack;
      fStack=pvue;
      break;
 	  }
    case '+':
      if (!fStack || !fStack->fNext)
	break;
      *(fStack->fNext)+=*fStack;
      pvue=fStack->fNext;
      fStack->fNext=0;
      delete fStack;
      fStack=pvue;
      break;

    case '-':
      if (!fStack || !fStack->fNext)
	break;
      *(fStack->fNext)-=*fStack;
      pvue=fStack->fNext;
      fStack->fNext=0;
      delete fStack;
      fStack=pvue;
      break;

    case '*':
      if (!fStack || !fStack->fNext)
	break;
      *(fStack->fNext)*=*fStack;
      pvue=fStack->fNext;
      fStack->fNext=0;
      delete fStack;
      fStack=pvue;
      break;

    case '/':
      if (!fStack || !fStack->fNext)
	break;
      *(fStack->fNext)/=*fStack;
      pvue=fStack->fNext;
      fStack->fNext=0;
      delete fStack;
      fStack=pvue;
      break;

//  AUTRE CHOSE : ON PASSE A LA VUE DU DESSUS DE LA PILE
    default :
      if (fStack)
	fStack->HandleInput(px);
      
    }
    Affiche();
  }
}


void LitmFenetre::Affiche() {

//    LitmVue *vue;
  cd();
  
  if (fStack) {
    switch(fStack->GetType()) {
    case 0:
      SetLogx(0);
      SetLogy(0);
      break;
    case 60:
      SetLogx(1);
      SetLogy(1);
      break;
    case 70:
      SetLogx(0);
      SetLogy(1);
      break;
    case 100:
      SetLogx(1);
      SetLogy(1);
      SetLogz(1);
      break;
    default:
      SetLogx(0);
      SetLogy(0);
      SetLogz(0);      
    }
    fStack->DrawAll(1,fStack->GetType());
    fStack->PrintAll(1,fStack->GetType());
  }
  
  Update();

}
