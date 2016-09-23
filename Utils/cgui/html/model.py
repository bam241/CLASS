#!/usr/bin/env python 
# -*- coding: utf-8 -*-

import htmlHelper as htmlH

emptyModel = """{
	type:'',
	weightU235:0.791135 , weightPu238:0.686385 , weightPu240:0.13553 , weightPu241:1.54572 ,weightPu242:0.0829001 , weightAm241:-0.336945 , equivalentFissile:10.3213 ,
	tmvaWeightPath:'' , keffTarget:1.00 ,
	nbBatches:3 , kThreshold:1.034 ,
	slist:[] ,
	reaclist:[],
	display:'btn-info',
	fp:{
		time:2,
		storageManagement:'LiFo',
		storageList: [
			{ key:'Fertile' , backendFacility: [] , infinit:false },
			{ key:'Fissile' , backendFacility: [] , infinit:false },
		],
	},
	im:{
		model: 'IM_RK4',
		spectrumType:'thermal',
	}
}"""

def importList () :
	"""
		generate list in ng-init directive
	"""
	data = { 'evolutionData':str(htmlH.listEvoData()) , 'equivalencemodel':str(htmlH.listEQModels()) , 'MLPWeight':str(htmlH.listMLPWeight()) , 'XSModel':str(htmlH.listXSModels()) , 'SMlist':str(htmlH.listStorageManagement()) , 'IrModel':str(htmlH.listIModels()) }

	return """EDM=['Evolution_Data','Equivalence_Model'] ; EvolutionDataList=%(evolutionData)s ; EquivalenceModelList=%(equivalencemodel)s ; MLPWeight=%(MLPWeight)s ; XSModelList=%(XSModel)s ; StorageManagementList=%(SMlist)s ; IrradiationModelList=%(IrModel)s"""%data

def strModel ( eqm ) :
	"""
		return the formulaire of an equivalence model as a string
	"""
	return htmlH.strHTMLfile( "html/model/EQM/" + eqm + ".html" )

def modelList () :
	"""
		browse models and stringify the switch and formulaire for each one
	"""
	r = ""
	for model in htmlH.listEQModels() :
		r += """<!-- %(model)s -->
			<li class="list-group-item" ng-switch-when="%(model)s" >
				%(form)s
			</li>
			""" % {'model':model,'form':strModel(model)}

	return r

def equivalenceModel () :
	"""
		display the all form for Equivalence Model
	"""
	eqmForm = htmlH.strHTMLfile("html/model/equivalenceModel.html")
	return eqmForm%{ 'eqmList':modelList() }


##############################################################################

def form () :
	model = htmlH.strHTMLfile("html/model/model.html").format(htmlH.strHTMLfile("html/model/evolutionData.html"),equivalenceModel(),htmlH.strHTMLfile("html/model/reactor.html"))

	return """
	<div class="col-sm-12" ng-repeat="m in models" ng-init="%(lists)s" >
		%(model)s
	</div>
	""" % { 'lists':importList() , 'model':model }
