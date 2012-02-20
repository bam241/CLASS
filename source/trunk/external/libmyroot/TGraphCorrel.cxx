#include "TGraphCorrel.hxx"
#include <tnt/lapack.h>
#include "NextName.hxx"
#include <fstream>
using namespace std;

bool TGraphCorrel::InitJackKnife(int nGraphs, TGraph **Graphs) {

  // Verification que les graphes sont la.
  if (!nGraphs)
    return false;
  for (Int_t i=0;i<nGraphs;i++)
    if (!Graphs[i])
      return false;

  // pointeurs sur les x et les y des graphes.
  Double_t **x =new (Double_t *)[nGraphs];
  Double_t **y =new (Double_t *)[nGraphs];
  Double_t *xthis;
  Double_t *ythis;

  for (Int_t i=0;i<nGraphs;i++) {
    x[i]=Graphs[i]->GetX();
    y[i]=Graphs[i]->GetY();
  }

  // Verification que les graphes ont le meme nombre de points.
  Int_t nx;
  nx=Graphs[0]->GetN();
  for (Int_t i=1;i<nGraphs;i++)
    if (Graphs[i]->GetN() != nx)
      return false;
  Set(nx);
  fCorrelations.newsize(nx,nx);

  xthis=GetX();
  ythis=GetY();

  // Verification qu'ils ont les memes x et remplissage du graphe
  // avec la moyenne des entrees.
  for (Int_t n=0;n<nx;n++) {
    xthis[n]=x[0][n];
    ythis[n]=y[0][n];
    for (Int_t i=1;i<nGraphs;i++) {
      if (x[i][n]!=xthis[n])
	return false;
      ythis[n]+=y[i][n];
    }
    ythis[n]/=nGraphs;
  }

  // Remplissage des variances et correlations.
  // Les correlations dans la moyenne sont (n-1)^2/n fois
  // les correlations des entrees.

  for (Int_t i=0;i<nx;i++)
    for (Int_t j=i; j<nx;j++) {
      Double_t correl=0;
      
      for (Int_t n=0;n<nGraphs;n++)
	correl+=(y[n][i]-ythis[i])*(y[n][j]-ythis[j]);

      fCorrelations(i+1,j+1)=fCorrelations(j+1,i+1)=correl*(nGraphs-1)/(nGraphs);
      if (i==j)
	SetPointError(i,0,sqrt(fCorrelations(i+1,i+1)));
    }
  int info=Upper_symmetric_eigenvaluevectors_solve(fCorrelations,fValeursPropres, fVecteursPropres);

  if (info)
    return false;

  MinVP=0;

  return true;

}

bool TGraphCorrel::Init(TGraphErrors *Graph) {

  if (!Graph)
    return false;
  
  int nx=Graph->GetN();
  Set(nx);
  fCorrelations.newsize(nx,nx);
  Double_t *xG,*yG,*xT,*yT;
  xG=Graph->GetX();
  yG=Graph->GetY();
  xT=GetX();
  yT=GetY();

  // Verification qu'ils ont les memes x et remplissage du graphe
  // avec la moyenne des entrees.
  for (Int_t n=0;n<nx;n++) {
    xT[n]=xG[n];
    yT[n]=yG[n];
  }

  // Remplissage des variances et correlations.

  for (Int_t i=0;i<nx;i++)
    for (Int_t j=0;i<nx;i++)
      fCorrelations(i+1,i+1)=0;
  
  for (Int_t i=0;i<nx;i++) {
    fCorrelations(i+1,i+1)=Graph->GetErrorY(i)*Graph->GetErrorY(i);
    SetPointError(i,0,Graph->GetErrorY(i));
  }
  int info=Upper_symmetric_eigenvaluevectors_solve(fCorrelations,fValeursPropres, fVecteursPropres);
  
  if (info)
    return false;
  
  MinVP=0;
  
  return true;
  
}

Double_t TGraphCorrel::Chi2(TF1 *F1) {

  Int_t N=GetN();
  Double_t *Xtab=GetX();
  Double_t *Ytab=GetY();
  Vector<double> vec;
  Double_t chi2=0;
  vec.newsize(N);

  for (int i=0;i<N;i++)
    vec(i+1)=F1->Eval(Xtab[i])-Ytab[i];

  chi2=0;

  Double_t s;
  for (int vp=0;vp<N;vp++)
    if (fValeursPropres(vp+1)>=MinVP) {
      s=0;
      for (int i=0;i<N;i++)
	s+=fVecteursPropres(i+1,vp+1)*vec(i+1);
      chi2+=s*s/fValeursPropres(vp+1);
    }
  
  return chi2;

}

Double_t TGraphCorrel::Chi2Indep(TF1 *F1) {

  Double_t temp;
  Double_t *Xtab=GetX();
  Double_t *Ytab=GetY();
  Double_t chi2=0;

  for (int i=0;i<GetN();i++) {
    temp=F1->Eval(Xtab[i])-Ytab[i];
    chi2+=temp*temp/fCorrelations(i+1,i+1);
  }

  return chi2;

}

Int_t TGraphCorrel::NbVp() {

  Int_t retour=0;

  for (int vp=0;vp<GetN();vp++)
    if (fValeursPropres(vp+1)>=MinVP)
      retour ++;
  
  return retour;

}

TH2D *TGraphCorrel::FillTH2withCorrel() {

  int N=GetN();

  TH2D *histo=new TH2D(NextName(),"Correlation matrix",N,0,N,N,0,N);
  
  Double_t s;
  
  for (Int_t i=0;i<N;i++)
    for (Int_t j=0; j<N;j++) {
      s=0;
      for (int vp=0;vp<N;vp++)
	if (fValeursPropres(vp+1)>=MinVP)
 	  s+=fValeursPropres(vp+1)*fVecteursPropres(i+1,vp+1)*fVecteursPropres(j+1,vp+1);
      histo->SetBinContent(i+1,j+1,s);
    }

  histo->SetStats(0);

  return histo;

}

TH2D *TGraphCorrel::FillTH2withInvCorrel() {

  int N=GetN();

  TH2D *histo=new TH2D(NextName(),"Inverse Correlation matrix",N,0,N,N,0,N);
  
  Double_t s;
  
  for (Int_t i=0;i<N;i++)
    for (Int_t j=0; j<N;j++) {
      s=0;
      for (int vp=0;vp<N;vp++)
	if (fValeursPropres(vp+1)>=MinVP)
	  s+=fVecteursPropres(i+1,vp+1)*fVecteursPropres(j+1,vp+1)/fValeursPropres(vp+1);
      histo->SetBinContent(i+1,j+1,s+histo->GetBinContent(i+1,j+1));
    }
  
  histo->SetStats(0);

  return histo;

}

TH2D *TGraphCorrel::FillTH2withVecteursPropres() {

  int N=GetN();

  TH2D *histo=new TH2D(NextName(),"Eigen Vectors",N,0,N,N,0,N);
  
  for (Int_t i=0;i<N;i++)
    for (Int_t j=0; j<N;j++)
      histo->SetBinContent(i+1,j+1,fVecteursPropres(j+1,N-i));
  
  histo->SetStats(0);

  return histo;

}

TH1D* TGraphCorrel::FillTH1WithVP() {

  Int_t N=GetN();
  Int_t nbins=N;
  if (nbins<1)
    nbins=1;
  Double_t min,max;

  min=max=fValeursPropres(1);
  for (int i=1;i<N;i++) {
    if ((fValeursPropres(i+1)<min || min<0) && fValeursPropres(i+1)>0)
      min=fValeursPropres(i+1);
    if (fValeursPropres(i+1)>max)
      max=fValeursPropres(i+1);
  }
  
  min/=100;
  max*=10;
  min=log10(min);
  max=log10(max);
  
  TH1D *VPH1=new TH1D(NextName(),"Eigenvalues repartition",nbins,min,max);
  for (int i=0;i<N;i++)
    if (fValeursPropres(i+1)>0)
      VPH1->Fill(log10(fValeursPropres(i+1)));

  VPH1->SetStats(0);

  return VPH1;

}
