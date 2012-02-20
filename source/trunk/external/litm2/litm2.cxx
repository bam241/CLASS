#include "fenetre.hxx"
#include <TROOT.h>
#include <assert.h>
#include <TApplication.h>
#include <libmyroot/NextName.hxx>


int main(int argc,char **argv) {

  assert(argc>1);

  TROOT MyRoot(NextName(),"");  
  TApplication theApp(NextName(), &argc, argv);
  argc=theApp.Argc();
  argv=theApp.Argv();

  cout << "Appel : litm2 fichier1m fichier2m ...\n";
  cout << "\n";
  cout << "Un fichier contient plusieurs tallies.\n";
  cout << "Un tally contient des valeurs index�es par :\n";
  cout << "- lieu (cellule ou surface)\n";
  cout << "- flag total/direct ou flagged/unflagged\n";
  cout << "       (tallies de d�tecteur, cartes SF, CF)\n";
  cout << "- bin utilisateur (carte FU)\n";
  cout << "- segment (carte FS) \n";
  cout << "- multiplicateur (carte FM)\n";
  cout << "- bin de cosinus (carte C) *** non test� ***\n";
  cout << "- bin d'�nergie (carte E)\n";
  cout << "- bin de temps (carte T)\n";
  cout << "\n";
  cout << "Litm2 permet d'afficher � l'aide de ROOT un ou plusieurs\n";
  cout << "histogrammes 1D ou 2D qui sont des coupes de tallies.\n";
  cout << "Les histogrammes sont g�r�s comme une pile.\n";
  cout << "L'�tat de la pile est affich� dans la sortie standard (terminal).\n";
  cout << "Certaines op�rations sont possibles sur les histogrammes.\n";
  cout << "\n";
  cout << "Commandes interactives (la fen�tre ROOT doit avoir le focus) :\n";
  cout << "(toutes les commandes s'appliquent aux histogrammes\n";
  cout << "       en haut de la pile)\n";
  cout << "n : ajout d'un nouvel histogramme\n";
  cout << "<enter> duplication d'histogramme\n";
  cout << "<backspace> suppression d'histogramme\n";
  cout << "<space> : �change de deux histogrammes\n";
  cout << "> : rotation de la pile\n";
  cout << "+ : addition de deux histogrammes\n";
  cout << "- : soustraction de deux histogrammes\n";
  cout << "* : multiplication de deux histogrammes\n";
  cout << "/ : division de deux histogrammes\n";
  cout << "p : (p comme produit) Multiplie un histo par une constante entr�e au clavier\n";
  cout << "W : �criture d'un histogramme dans un fichier .dat\n";
  cout << "             lisible par xmgrace\n";
  cout << "F,f : passage au fichier suivant, pr�c�dent\n";
  cout << "L,l : passage au lieu suivant, pr�c�dent\n";
  cout << "D,d : passage direct <-> total ou flagged <-> unflagged\n";
  cout << "U,u : passage au bin utilisateur suivant, pr�c�dent\n";
  cout << "S,s : passage au segment suivant, pr�c�dent\n";
  cout << "M,m : passage au multiplicateur suivant, pr�c�dent\n";
  cout << "C,c : passage au bin de cosinus suivant, pr�c�dent\n";
  cout << "E,e : passage au bin d'�nergie suivant, pr�c�dent\n";
  cout << "T,t : passage au bin de temps suivant, pr�c�dent\n";
  cout << "H,h : passage au tally suivant, pr�c�dent\n";
  cout << "o : passage au bin de temps 0, bin d'�nergie 0, bin de cosinus 0\n";
  cout << "; : Lieu (num�ro) en abscisse (n�-lin) (/dt/dlne/dc)\n";
  cout << "r : Energie (MeV) en abscisse (log-log) (/dt/dlne/dc)\n";
  cout << "y : Temps (microsecondes) en abscisse (lin-log) (/dt/dlne/dc) [DEFAUT]\n";
  cout << "b : bipas x=temps(microsecondes)-y=energie(MeV) (log-log-log) (/dlnt/dlne/dc)\n";
  cout << "\n";
  cout << "\n";

  LitmFenetre MaFenetre;
  for (int i=1;i<argc;i++)
    MaFenetre.LoadMctal(argv[i]);
  
  theApp.Run();
  return 0;

}
