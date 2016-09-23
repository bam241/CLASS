# `cgui`

`cgui` (pour *CLASS GUI*) est une page Web permettant la génération automatique de code pour l'input de *CLASS*.

## Lancement du serveur

`cgui` fonctionne avec un serveur Web, il est nécessaire dans un premier temps de lancer celui-ci. Le serveur Web est compatible avec Python 2.7 et 3.*, pour le lancer il suffit d'exécuter le script `code/server.py` :

```bash
cd code
./server.py
```

La sortie devrait alors ressembler à ceci :

```bash
$ ./server.py 
Server HTTP on localhost -- port: 8080
	pid: 16174
```

La valeur du `pid` est stocké dans un fichier (`code/server.pid`), cela sert à fermer le serveur proprement.

Voilà, vous avez lancer le serveur, il ne reste plus qu'à vous connecter !

### Problèmes potentiels

Certains ports de votre machine peuvent être fermé par défaut, vous pouvez alors essayer un autre port que le `8080` en modifiant la valeur de la variable `port` dans le fichier `code/server.py`, un autre port souvent utilisé pour ce type d'utilisation est le `8000`.

L'utilisation des ports inférieurs à `1024` requière des droits *root*, c'est pour contourner ce problème que j'ai décidé de lancer le serveur sur le port `8080`.

## Utilisation de la GUI

### Connexion

La *GUI* s'utilise à l'aide d'un navigateur Internet, pour le moment elle a été testée sur les suivants :

* Opera 37, 39
* Firefox 46, 48
* Chromium 50, 52
* Midori 0.4.3

> À aucun moment il n'a été envisagé de tester la *GUI* sur InternetExplorer ou son successeur

Elle est écrite en HTML/JavaScript/CSS, il est donc nécessaire d'activer JavaScript (fonctionnement par défaut d'un navigateur), donc la *GUI* ne fonctionnera pas avec un navigateur en ligne de commande comme *lynx*.

Pour l'utiliser il suffit d'aller sur la page : [http://localhost:8080/cgui.py](http://localhost:8080/cgui.py).

### How to

La page est un long formulaire donc le résultat apparaît dynamiquement en bas dans la section **Code** (la section **Diagram** est en cours de développement mais devrait permettre à terme d'afficher un petit diagramme représentant le scénario).

Le formulaire est en trois grandes parties :

* **Simulation configuration** il s'agit de la configuration globale de la simulation.
* **Stock, pool and SP creation** il s'agit de la création des *backend facilites*. Le bouton vert à droite **Add unity** permet d'ajouter des unités. Vous pouvez supprimez les unités indésirables en cliquant sur la petite croix jaune. Pour chaque unité choisissez son type et suivez les instruction. Pensez à lui donner un nom, celui-ci sera le nom de la variable *C++*.
* **Models and reactors creation** il s'agit de la création des modèles. La GUI est centrée sur le modèle, c'est à partir de lui que l'on va créer une liste de réacteurs utilisant ce modèle. L'ajout de modèle s'effectue en cliquant sur le petit bouton vert **Add model**. L'ajout de réacteur propre à ce modèle s'effectue en cliquant sur le petit bouton jaune **Add reactor**.

Les formulaires sont normalement suffisamment clairs (en tout cas autant que la doc), mais n’hésitez pas à venir poser des questions, signaler des bugs ou effectuer des demandes d'ajouts de fonctionnalités.

### « J'ai fini ! »

Une fois le code généré vous pouvez le télécharger (ou faire un copier/coller ça marche aussi).

> Pensez à éteindre le serveur, pour cela un joli bouton rouge est présent tout en bas de la page ;( (s'il n'a pas été retiré), ou allez sur la page [http://localhost:8080/kill.py](http://localhost:8080/kill.py) ou un simple `ctrl`+`C` dans le terminal de lancement du serveur fonctionne aussi.


## Utilisation avancée

Il est possible d'avoir le serveur Web sur une machine distante, par exemple à Subatech je l'ai testé sur `nansl3`. Il faut vérifier si le port `8080` est ouvert sur la machine d'accueil du serveur (ou s'il n'est pas déjà utilisé par une autre application). La connexion à l'interface ne se fera plus depuis `localhost` mais depuis `http://adresse.ip.du.serveur:8080/cgui.py` (cela fonctionne aussi avec le DNS, donc à Subatech sur : `http://nansl3.in2p3.fr:8080/cgui.py`).

-----

# Modifier la GUI

En cas de modification de *CLASS* il est nécessaire de modifier certains éléments de la GUI.

## Ajout d'une *backend facility*

Dans notre exemple on ajoutera une *backend facility* qui s'appelle **Mine**.

> Je sais la mine est une *frontend facility* mais imaginons que dans un futur très lointain on oublie l'emplacement de certains stocks et que l'on se mette à creuser à cet endroit, l'ancien stock deviendra une mine. Par conséquent la mine est bien la *backend facility* d'un réacteur de notre époque. CQFD.
> Plus sérieusement, le travail à faire du des hypothétiques futurs *frontend facilities* serait similaire à celui sur les *backend facilities* (sauf la *fabrication plant* qui est intégré au modèle).

### 1. Ajouter dans la liste des *backend facilites*

Dans le fichier `html/backend.py` ajouter un élément au tableau `backendType` (l6).

```python
backendType = ['Pool','SeparationPlant','Stock','Mine']
```

Dans le fichier `js/app-cgui.js` ajouter un élément au tableau `$scope.backEndType` (l65).

```js
$scope.backendType=['Pool','SeparationPlant','Stock','Mine']; // l'ordre ici sera l'ordre d'apparition dans la liste déroulante
```

### 2. Ajouter le formulaire de construction

Dans le dossier `html/backend/` ajouter un fichier HTML contenant le formulaire de construction de la nouvelle *backend facility* : `Mine.html` dans notre exemple (**/!\\** le nom doit être le même que celui renseigné dans la liste `backendType` du fichier `html/backend.py`). Un fichier `example.html` est déjà présent comme squelette pour une nouvelle *backend facility*. Cette exemple est plutôt complet et propose les différentes possibilités d'un formulaire.

L'application utilise le *framework* **AngularJS** pour gérer les formulaires et donc les objets et variables qui contiennent les informations permettant de générer le code C++. De plus quelques connaissances de base en HTML sont préférables (juste savoir que c'est un langage de balises et comment s'écrit une balise).

Pour des bases d'Angular je vous conseille la playlist de [Grafikart](https://www.youtube.com/watch?v=aBE0St5yI7U&list=PLjwdMgw5TTLUDlJyx4yIPQjoI-w-7Zs1r) (en français), les 2 premières vidéos devraient suffire à comprendre et modifier un formulaire AngularJS de base.

Pour notre exemple notre fichier `Mine.html` ressemblera à ceci :

```xml
<div class="row" >
	<div class="col-sm-12" >
		<div class="input-group " >
			<label class="input-group-addon" >Ville : </label>
			<input class="form-control" type="text" placeholder="description de ce qu'il faut mettre ici..." ng-model="unity.town" />
		</div>
	</div>
</div>
<div class="row" >
	<div class="col-sm-6" >
		<div class="input-group " >
			<label class="input-group-addon" >Nombre de travailleurs</label>
			<input class="form-control" type="number" ng-model="unity.workers" />
		</div>
	</div>
	<div class="col-sm-6" >
		<div class="input-group " >
			<label class="input-group-addon" >Bonheur des travailleurs</label>
			<input class="form-control" type="range" min="-1" max="1" step="0.421" ng-model="unity.happy" />
			<span class="input-group-addon" data-toogle="En labor happiness" >{{ unity.happy | number:2 }} $L \cdot H$</span>
		</div>
	</div>
</div>
```

### 3. Modifier la génération de code

L'étape sans doute la plus compliquée est de modifier le fichier `html/code.html`. Tout d'abord repérez la zone où sont définis les *Facilites* (l52) autre que les stocks (l59). (Un filtre Angular permet de gérer dans un premier temps seulement les stocks et dans un second temps tout le reste).

Ajouter à la suite du `ng-switch` le code correspondant à la construction de votre nouvelle unité. Pour notre mine cela pourrait ressembler à ceci.

```xml
<span ng-switch-when="Mine" >
	Mine * {{ unity.id }} = new Mine( gCLASS->GetLog() , {{ unity.workers }} ); // constructeur de ma mine
	{{ unity.id }}->SetTownName( "{{ unity.town }}" );                          // appel de setters en fonction des paramètres
	{{ unity.id }}->SetHappinessOfWorkers( 1.0*{{ unity.happy }} );
	<span ng-show="unity.happy > 0.75" ><!-- affichage conditionnel des valeurs -->
	{{ unity.id }}->SetStrike( false );
	</span>
</span>
```

> Pour connecter une unité à plusieurs unités de sorties ou d'entrée, le formulaire de la `SeparationPlant` pourrait sans doute aider.


## Ajouter une *Evolution Data*

La GUI étant liée à CLASS, elle récupère directement les fichiers `*.dat` dans tous les dossiers `FixedFuel` dans `DATA_BASES` (quelque soit la profondeur du dossier `FixedFuel`). Donc si vous avez bien fait votre boulot c'est fini.

## Ajouter un *Equivalence Model*

La GUI étant liée à CLASS, elle récupère directement la liste des classes qui héritent publiquement d'`EquivalenceModel` et dont le fichier d'entête (extension en `.hxx`) est situé dans `$CLASS_PATH/source/Model/Equivalence/`. Il est maintenant nécessaire d'ajouter le constructeur de votre nouveau modèle à la *GUI*.

### 1. Ajouter le formulaire de construction

Comme pour les *backend facilites* il faut ajouter un fichier HTML avec le formulaire du constructeur, celui-ci se trouvera dans `html/module/EQM/`. Un fichier `example.html` présente un squelette pour un nouveau modèle d'équivalence (le principe est similaire aux unités).

### 2. Modifier la génération de code

L'étape sans doute la plus compliquée est de modifier le fichier `html/code.html`. Il faut repérer la zone où sont définis les *Equivalence Models* (l79).

Puis il faut ajouter la directive `ng-switch` correspondant à vote nouveau modèle (actuellement le dernier modèle se fini à la ligne 133).

```xml
<span ng-switch-when="EQM_Mon_NouveauModel" >
	EQM_Mon_NouveauModel * eqm_{{ m.id }} = new EQM_Mon_NouveauModel( gCLASS->GetLog() , {{ m.tmvaWeightPath }} , ... );
</span>
```

## Ajouter un *XS Model*

La GUI étant liée à CLASS elle récupère directement la liste des dossiers contenu dans un dossier `XSModel` contenant un fichier `*.nfo`. Attention la GUI est sensible à la casse !

## Ajouter des poids pour les reseaux de neuronnes

La GUI étant liée à CLASS elle récupère directement la liste des fichiers `*.xml` contenu dans un dossier `EQModel`. Attention la GUI est sensible à la casse !

## Ajouter un *Irradiation model*

La GUI étant liée à CLASS elle récupère la liste des fichiers d'entête (`*.hxx`) contenu dans `$CLASS_PATH/source/Model/Irradiation/`, il est donc nécessaire que le nom du fichier soit le même que celui de la classe.

## Modifier un formulaire

* Le formulaire de configuration de la simulation est : `html/config.html`.
* Le formulaire mère des *backend facilities* est : `html/backend/backend.html`. L'utilisation du `{0}` indique à Python (script `html/backend.py`) de remplacer cette section par la liste des formulaires des *backend facilities*.
* Les formulaires des *backend facilites* sont dans `html/backend/backend/`.
* Le formulaire mère des *models* est : `html/model/model.html`. L'utilisation des `{0}`, `{1}` et `{2}` indique à Python (script `html/model.py`) de remplacer ces sections par les formulaires des données d'évolution, des modèles d'équivalence et des réacteurs.
* Le formulaire des *evolution data* est : `html/model/evolutionData.html`.
* Le formulaire des *equivalence models* est : `html/model/equivalenceModel.html`.
* Le formulaire des réacteurs est : `html/model/reactor.html`.

## La génération de code

Pour le moment toute la génération de code s'effectue dans `html/code.html`.


-----

# Organisation du projet

## Choix des technologies

Pour le *CSS* il s'agit d'un choix de rapidité de développement :

* [Bootstrap](http://getbootstrap.com/)

Pour le *JS* il y a plusieurs solutions pour des problèmes différents :
	
* [AngularJS](https://docs.angularjs.org/api) pour la gestion des formulaires et génération du code
* [Highcharts](http://www.highcharts.com/docs) pour la génération de graphiques
* [JointJS](http://www.jointjs.com/api) pour la génération de diagrammes (celle-ci nécessite *Backbone*, *JQuery* et *Lodash*)
* [KaTeX](https://khan.github.io/KaTeX/) pour l'affichage des formules (plus rapide et surtout moins lourd que le plus connu [MathJax](https://www.mathjax.org/) mais moins joli)

## Arborescence

Les fichiers `__init__.py` sont des fichiers de configuration pour Python, les fichiers `*.pyc` sont les fichiers générés par Python (Python opère une phase de compilation avant de s'exécuter, cette phase est automatique mais laisse des fichiers `*.pyc`). Ces deux types de fichiers ne seront pas présentés ici.

En plus de ces fichiers, le serveur génère un fichier `server.pid` contenant uniquement le *PID* du serveur Web pour le killer via la page [http://localhost:8080/kill.py](http://localhost:8080/kill.py).

### `.`

Il s'agit ici de la racine du serveur Web, on y trouve donc le script de lancement du serveur : `server.py` ainsi que la page principale de l'application : `cgui.py`.

* `css/` dossier contenant les feuilles de style, différentes polices et le logo de *CLASS*
* `html/` dossier contenant la génération de la page HTML de l'application
* `js` dossier contenant le JavaScript et les différentes libraires JavaScript utilisées (aucun [CDN](https://en.wikipedia.org/wiki/Content_delivery_network#Free_CDNs) n'est utilisé)
* `cgui.py` page principale (unique ?) de l'application
* `kill.py` page permettant de killer le serveur (non indispensable)
* `run.py` page qui permettrait d'envoyer les données du scénario de l'utilisateur à un serveur de calcul (donc s'occuperait de la génération du code, de vérifications de l'input puis de l'envoie du job), laisser à titre indicatif pour des développements futurs
* `server.py` script lançant le serveur Web Python


### `css/`

Le CSS est contenu dans le dossier `css/`, et il est composé de :

	* Le fichier `css/cgui.css` qui est le fichier d'adaptation de *Bootstrap* à la *GUI* ainsi que quelques ajouts divers.
	* Le fichier `css/fonts.css` qui va charger les polices contenus dans le dossier `css/fonts` (essentiellement esthétique mais aussi icones)
	* Le dossier `css/lib` qui contient une version minimisée de *Bootstrap* (suppression des espaces et commentaires, tout en une seule ligne)
	* L'image `css/logo.png` qui est le logo de *CLASS*

### `html/`

Le HTML et une part importante du Python est contenu dans le dossier `html/`. Les fichiers `*.html` ne sont que des extraits de formulaires qui seront assembler par les différentes fonctions Python. Chaque fichier Python représente une fonctionnalité :

* `html/htmlHelper.py` contient un ensemble de fonction pour simplifier l'écriture de code redondant, ce fichier est importé dans tous les autres fichiers Python en temps que `htmlH`
* `html/backend.py` contient la génération du formulaire traitant des *backend facilities*, pour cela il va utiliser les fichiers HTML contenus dans `html/backend/`
* `html/model.py` contient la génération du formulaire traitant des *models*, pour cela il va utiliser les fichiers HTML contenus dans `html/model/`
* `html/code.html` est la génération du code, pour le moment tout ce fait de façon brute mais il est envisagé d'éclater le fichier pour n'avoir que des formulaires indépendants par constructeur
* `html/config.html` fichier HTML contenant la gestion de la configuration du scénario
* `html/diagram.html` fichier HTML contenant le morceau de HTML où sera généré le diagramme
* `html/footer.html` fichier HTML content le pied de la page HTML (appels aux librairies JS et fermeture de la page)
* `html/graph.html` fichier HTML contenant le morceau de HTML où sera généré le graphique de la puissance électrique
* `html/header.html` fichier HTML contenant l'entête de la page HTML (titre de la page, appel aux fichiers CSS et ouverture de la page)

### `js/`

Tout le JavaScript de la *GUI* est contenu dans le dossier `js/` qui est composé de :

* `js/lib` qui contient les différentes librairies minimisées utilisées dans la *GUI*
* `js/app-cgui.js` qui est l'application AngularJS de la *GUI*

Tout le reste c'est des trucs qui ne devraient pas être là.
