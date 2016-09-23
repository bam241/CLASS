#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import os
import fnmatch, glob, re

import html.htmlHelper as htmlH
import html.model
import html.backend

titlePage = "Plop"

def pattern ( title , options , content ) :
	return htmlH.sectionTitle( title , options ) + "\n" + content

data = { \
'contentType': "Content-type: text/html;" ,\
'header'     : htmlH.header(titlePage) ,\
'config'     : htmlH.strHTMLfile('html/config.html') ,\
'backend'    : pattern( "Stock, pool and SP creation" , htmlH.labelAddUnity("backend",html.backend.emptyBackendFacility, "unity" ) , html.backend.form() ) ,\
'model'      : pattern( "Models and reactors creation" , htmlH.labelAddUnity("models",html.model.emptyModel,"model") , html.model.form() ) ,\
'diagram'    : pattern( "Diagram" , "<button class=\"label label-info btn btn-right\" ng-click='makeDiagram()' >Make diagram</button>" , htmlH.strHTMLfile('html/diagram.html') ) ,\
'graph'      : pattern("Graphic" , "<button class=\"label label-default btn btn-right\" ng-click=\"exportData()\" >Export</button> <button ng-click=\"electricPowerGraph()\" class=\"label label-info btn btn-right\" >Make graph</button>" , htmlH.strHTMLfile("html/graph.html") ) ,\
'code'       : pattern( "Code" , "<button class=\"label label-default btn btn-right\" ng-click=\"downloadCode()\" >Download</button>" , htmlH.strHTMLfile('html/code.html') ) ,\
'footer'     : htmlH.footer() \
}

print( """%(contentType)s

%(header)s

<!-- ===== CONFIGURATION DE LA SIMULATION ================================ -->
%(config)s

<!-- ===== AJOUT DES STOCK/POOL/SEPARATION-PLANT ========================= -->
%(backend)s

<!-- ===== AJOUT DES MODELS/REACTEURS ==================================== -->
%(model)s

<!-- ===== DIAGRAM ======================================================= -->
%(diagram)s

<!-- ===== GRAPHICS ====================================================== -->
%(graph)s

<!-- ===== CODE GENERATOR ================================================ -->
%(code)s

<!--Bouton pour l'envoie des données sur le serveur, pour le moment cette fonctionnalité est désactivée. Merci de votre compréhension.-->
<label class="btn btn-lg btn-block" ng-class="{ 'btn-warning':send=='yes' , 'btn-info':send=='no' , 'btn-success':send=='ok' , 'btn-danger':send=='fail' }" style="margin-bottom: 1em;"  >
	<input type="radio" ng-model="send" value="yes" style="display:none;" />
	<i class="glyphicon glyphicon-send" ></i> Send to server
</label>


<label class="btn btn-lg btn-block btn-danger kill" style="margin-bottom: 1em;" >
	<input type="checkbox" data-ng-model="kill" style="display: none;" />
	Kill server ;(
	<br />
	<div data-ng-show="kill" style="width:25%%; margin: auto;" >
		<div class="row" ><strong class="col-sm-12" >Are you sure ?</strong></div>
		<div class="row" >
			<a href="kill.py" class="col-sm-6 text-muted" >Yes</a>
			<a ng-click="kill = false" class="col-sm-6 text-muted" >No</a>
		</div>
	</div>
</label>
%(footer)s
""" % data )

