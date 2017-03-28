class Scenar_t
==============

Les objets de type **Scenar_t** sont utilisés pour décrire un *scénario* CLASS. En pratique les attributs de **Scenar_t** sont les paramètres qui identifient le scénario. Par exemple, dans une étude multiparamétrique ces attributs sont (au minimum) les paramètres avec lesquels on exécutes CLASS_Exec :

* Burn-Up
* Nbre de rechargements
* Facteur de charge
* Temps de refroidissement
* Temps de fabrication
* Vecteurs Isotopiques (Dans les différentes étapes du cycle combustible)
* ...

L'utilisation des objets **Scenar_t** implique que l'ensemble des observables d'intérêt, par exemple l'inventaire total des noyaux, ne soient accessibles que par des methodes internes de la class. L'objectif est de pouvoir calculer ce que l'on veut à la demande.
