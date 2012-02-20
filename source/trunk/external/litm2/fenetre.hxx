#ifndef _LitmFenetre_HXX_
#define _LitmFenetre_HXX_

#include <libmctal/TMctal.hxx>
#include <TROOT.h>
#include <TCanvas.h>
#include <libmyroot/THistoCalc.hxx>
#include <libmyroot/THisto2Calc.hxx>
#include <vector>

// Correspond a un histo dans la fenetre, plus toutes ses coordonnees
// de tally,...

class LitmVue {

public :
  
  LitmVue();
  LitmVue(LitmVue& vue);
  LitmVue(LitmVue& vue,double MyConst);
  ~LitmVue();

  void HandleInput(char c);
// DEPLACEMENTS
  // 'l' d�cr�mente, 'L' incr�mente le lieu.
  // 'd' d�cr�mente, 'D' incr�mente total/direct ou flagged/unflagged.
  // 'u' d�cr�mente, 'U' incr�mente le bin utilisateur.
  // 's' d�cr�mente, 'S' incr�mente le segment.
  // 'm' d�cr�mente, 'M' incr�mente le multiplicateur.
  // 'c' d�cr�mente, 'C' incr�mente le cosinus.
  // 'e' d�cr�mente, 'E' incr�mente l'�nergie.
  // 't' d�cr�mente, 'T' incr�mente le temps.
  // 'h' d�cr�mente, 'H' incr�mente le tally.
  // 'o' ("origine") remet le lieu, le temps, l'�nergie et le cosinus � 0.
//  CHOIX DE MODE
  // ';' Mode 0
  // 'r' Mode 60
  // 'y' Mode 70
  // 'b' Mode 100
//  ECRITURE
  // 'W' cr�e un .dat

  void ChangeFile(TMctal *fic, char *nom="", unsigned long numero=0); // s'initialise avec un nouveau fichier !

  void Print(int verbose=1,ostream& out=cout);
  // appel r�cursif des Print() sur la pile avec indexage.
  // Une * indique ceux qui sont affich�s (compatibles avec le type du
  // premier de la pile).
  void PrintAll(int n,unsigned type=70);

  void Draw();
  void DrawSame();
  // appel r�cursif des Draw() sur la pile
  // Le dessin ne se fait que si le type courant est type (m�me(s) axe(s)).
  void DrawAll(unsigned n,unsigned type);

  unsigned SetType(unsigned type);
  unsigned GetType() {return fType;}

  THistoCalc *GetHisto() {return fHisto;}

  //  void SetHistoMin(double min) {if(fHisto) fHisto->SetMinimum(min);}
  //  void SetHistoMax(double max) {if(fHisto) fHisto->SetMaximum(max);}

  //  float fMin,fMax;  // Les min et max sur tout le tally (depend de fType)
  unsigned long fNoFile; // Numero du fichier;
  LitmVue* fNext;

  // les op�rations entre vues ne se font que sur des vues de m�me type.
  LitmVue& operator += (LitmVue& v);
  LitmVue& operator -= (LitmVue& v);
  LitmVue& operator *= (LitmVue& v);
  LitmVue& operator /= (LitmVue& v);

protected :

  TMctal *fFile; // fichier m affiche. 0 FIGE LES DEPLACEMENTS DANS LE
 // FICHIER ET LES CHANGEMENTS DE TYPE. n'appartient pas.
  string fNomFichier; // nom du fichier m. Ne lui appartient pas.
  TMTally *fTally; // tally affiche. Ne lui appartient pas.
  unsigned long fNoTally; // Num�ro du tally dans le fichier

  unsigned fType;   // Type de plot (* pour ceux qui sont implementes)
  // Imp�rativement les ou les axes doivent �tre en accord avec fType,
  // m�me pour les histos produits par des op�rations.

  // 0 a 9 : lieu en abscisse,
  // ...
  // 70 a 79 : temps en abscisse
  // >=100 : 2D, energie versus temps

// *0   : Lieu (num�ro) en abscisse, valeur du tally en ordonn�e (n�-lin)
//                                          /dt/dlne/dc

// *60  : Energie (MeV) en abscisse, valeur du tally en ordonn�e (log-log)
//                                          /dt/dlne/dc

// *70  : Temps (microsecondes) en abscisse, valeur du tally en ordonn�e (lin-log) [DEFAUT]
//                                          /dt/dlne/dc

// 100 : bipas x=temps(microsecondes)-y=energie(MeV), valeur du tally en ordonn�e (log-log-log)
//                                                    /dlnt/dlne/dc

  unsigned long fL; // Les 8 coordonnees : lieu
  unsigned long fD; // Total/direct ou flagged/unflagged
  unsigned long fU; // Bin utilisateur
  unsigned long fS; // Segment
  unsigned long fM; // Multiplicateur
  unsigned long fC; // Cosinus
  unsigned long fnCbins; // Nombre de bins de cosinus
  Double_t *fCbins; // bins de cosinus (appartiennent)
  unsigned long fE; // Energie
  unsigned long fnEbins; // Nombre de bins d'energie
  Double_t *fEbins; // bins d'energie (appartiennent)
  bool zeroE; // Bin initial en energie de taille nulle ?
  unsigned long fT; // Temps
  unsigned long fnTbins; // Nombre de bins de temps
  Double_t *fTbins; // bins de temps (appartiennent)
  bool zeroT; // Bin initial en temps de taille nulle ?

  THistoCalc *fHisto;  // l'histo resultat. Appartient.
  THisto2Calc *fHisto2; // l'histo 2D resultat. Appartient.

  // Determinent les vrais degres de libert�s
  // Depend du mode et des dimensions du fichier et du tally
  // 1 si oui, 0 si non.
  int HVarie();
  int LVarie();
  int DVarie();
  int UVarie();
  int SVarie();
  int MVarie();
  int CVarie();
  int EVarie();
  int TVarie();

  void ChangeHisto();
  void FillHisto();

};


// Correspond a une fenetre avec evt plusieurs fichiers,plusieurs
// tallys,...

// Tous les histos ne sont pas forc�ment de m�me type.
// Seuls s'affichent ceux qui sont du type de l'histo en haut de la pile.
class LitmFenetre : public TCanvas {

public :

  LitmFenetre() : TCanvas("Litm2","Litm2") {Init();}
  ~LitmFenetre() {Free();}

  unsigned LoadMctal(char *nomfichier);
  void HandleInput(EEventType event, Int_t px, Int_t py); // appele par ROOT
//  GESTION DE PILE
  // 'n' Nouvelle vue
  // '\r' Duplication de vue
  // '\b' Suppression de vue
  // ' ' Swappe les deux dernieres vues
  // '>' Roule dans un sens
//  CHANGEMENT DE FICHIER
  // 'F' Incremente, 'f' decremente le fichier de la vue du dessus de la pile
//  OPERATIONS SUR LES HISTOS
  // '+' Addition
  // '-' Soustraction
  // '*' Multiplication
  // '/' Divise le deuxieme par le premier

protected :

  void Affiche();  // affiche la pile dans l'etat courant.

  void Init();
  void Free();
  double fMin,fMax; // min et max des histos sur toute la pile.
  vector<TMctal *> fFichiers; // les TMctal lui appartiennent.
  vector<char *> fNoms; // les noms ne lui appartiennent pas.
  LitmVue *fStack;

};

#endif
