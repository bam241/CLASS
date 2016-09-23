var app = angular.module('classGUI', ["highcharts-ng"]);

//____________________________________________________________________________
// filter to display : ZAI(z,a,i)*quantity (usefull in the code generator)
app.filter('ZAIQuantity', function() {
	// angularJS filter to display the code to add a initial IV in stock
	return function(string){
		var tab = JSON.parse("["+string+"]");
		var r = [];
		for ( var i=0 ; i<tab.length ; i++ )
		{
			r.push("ZAI("+tab[i][0]+","+tab[i][1]+","+tab[i][2]+") * "+tab[i][3]);
		}
		return r;
   }
});

//____________________________________________________________________________
// filter to display efficacity in a SP : z,a,i,efficacity.
app.filter('ZAIEfficacity', function() {
	// angularJS filter to display the code to add sepration IV in SP
	return function(string){
		var tab = JSON.parse("["+string+"]");
		var r = [];
		for ( var i=0 ; i<tab.length ; i++ )
		{
			r.push(""+tab[i][0]+","+tab[i][1]+","+tab[i][2]+" , "+tab[i][3]);
		}
		return r;
   }
});

//____________________________________________________________________________
// add a directive to AngularJS, like ng-click but for file selection in input type="file"
app.directive("ngFileSelect",function(){
	return {
		link: function($scope,el){
			el.bind("change", function(e){
				$scope.file = (e.srcElement || e.target).files[0];

				var reader = new FileReader();

				reader.onloadend = function(evt) {
					if (evt.target.readyState == FileReader.DONE) {
						$scope.$apply(function () {
							$scope.config  = JSON.parse( evt.target.result ).config;
							$scope.backend = JSON.parse( evt.target.result ).backend;
							$scope.models  = JSON.parse( evt.target.result ).models;
							setTimeout(reloadMathJax,10)
						});
					}
				};

				var blob = $scope.file.slice(0, $scope.file.size);
				reader.readAsBinaryString(blob);
			})
		}
	}
});

//____________________________________________________________________________
// Global, unique and lovely controller
app.controller('classCtrl', ['$scope', '$timeout', '$http' , function ($scope, $timeout,$http) {

	$scope.backendType=['Pool','SeparationPlant','Stock'];
	
	$scope.required = true;

	// initialize default empty config
	$scope.config = {
		timeStart:0,
		timeEnd:100,
		timeStep:0.01,
		outputFile:'',
		log: {termLevel:0,fileLevel:0,file:''},
		SetStockManagement: true,
		SetZAIThreshold: 89,
	};
	// initialize empty list for backend facilities ands models
	$scope.backend = [];
	$scope.models  = [];

	// init array of graphic of power
	$scope.chartSeries = [];

	// graphic configuration
	$scope.chartConfig = {
		options: {
			chart: {
				type: 'areaspline',
				zoomType: 'x',
			},
			plotOptions: {
				area: {
					stacking: 'normal',
				},
				series: {
					stacking: 'normal',
				}
			},
			navigator: {
				enabled: true
			},
		},
		title: {
			text: 'Electric power in the scenario'
		},

		series: $scope.chartSeries,
		
		loading: false,
		size: {},
		useHighStocks: true,
	};

	//________________________________________________________________________
	// function to download configuration of a scenario as a JSON file
	$scope.downloadScenar = function () {
		var filename = $scope.config.inputFile;

		// build data
		var mimeType = 'application/json';
		var raw_data = {
			config  : $scope.config,
			backend : $scope.backend,
			models  : $scope.models,
		};
		var data = mimeType  +  ';charset=utf-8,' + JSON.stringify(raw_data)

		// name of file
		if ( filename == undefined )
			{ filename = "myAwsomeParc.json"; }
		else { filename += ".json"; }

		// make a pseudo link to download data
		var link = document.createElement('a');
		
		link.setAttribute( 'download' , filename );
		link.setAttribute( 'href'     , 'data:' + data );
		document.body.appendChild(link);
		link.click();
	};

	//________________________________________________________________________
	// function to download input CLASS C++ code
	$scope.downloadCode = function () {
		var filename = $scope.config.inputFile;

		// build data
		var mimeType = 'text/plain';
		var code = document.getElementById("code").textContent;

		var data = mimeType + ';charset=utf-8,' + encodeURIComponent(code);

		// name of file
		if ( filename == undefined )
			{ filename = "myAwsomeParc.cxx"; }

		// make a pseudo link to download data
		var link = document.createElement('a');

		link.setAttribute('download', filename);
		link.setAttribute('href', 'data:' + data );
		document.body.appendChild(link);
		link.click();
	};

	//________________________________________________________________________
	// function to export data of the electric power as a csv (raw data of the graph)
	$scope.exportData = function () {
		var filename = "data.txt";

		data = "";
		for (var obj=1 ; obj<$scope.chartSeries.length ; obj++ ) // 0 is sum and we don't display sum
		{
			data += $scope.chartSeries[obj].name + " ";
		}
		data += "\n";
		for ( var i=0 ; i<$scope.chartSeries[0].data.length ; i++)
		{
			for (var obj=1 ; obj<$scope.chartSeries.length ; obj++ )
			{
				data += $scope.chartSeries[obj].data[i] + " ";
			}
			data += "\n";
		}

		var link = document.createElement('a');
		var mimeType = 'text/plain';

		link.setAttribute('download', filename);
		link.setAttribute('href', 'data:' + mimeType  +  ';charset=utf-8,' + encodeURIComponent(data) );
		document.body.appendChild(link);
		link.click();
	};

	//________________________________________________________________________
	// function which make the graphic of electric power
	$scope.electricPowerGraph = function () {
		$scope.chartSeries.length = 0;

		var size = ($scope.config.timeEnd - $scope.config.timeStart)/$scope.config.timeStep;
		var tr_start, tr_end;
		
		var sum = Array.apply(null, Array(size)).map(Number.prototype.valueOf,0);

		for ( var i=0 ; i<$scope.models.length ; i++ )
		{
			for ( var j=0 ; j<$scope.models[i].reaclist.length ; j++ )
			{
				var data = Array.apply(null, Array(size)).map(Number.prototype.valueOf,0);
				tr_start = ($scope.models[i].reaclist[j].tstart - $scope.config.timeStart)/$scope.config.timeStep;
				tr_end   =  Math.min(tr_start + $scope.models[i].reaclist[j].tlife/$scope.config.timeStep , size);
				
				for ( var t=tr_start ; t<tr_end ; t++)
				{
					data[t] = 1.0*$scope.models[i].reaclist[j].power*$scope.models[i].reaclist[j].loadFactor*$scope.models[i].reaclist[j].efficency;
					sum[t] += data[t];
				}
				
				$scope.chartSeries.push( {'name':$scope.models[i].reaclist[j].id , 'data': data } );
			}
		}
		$scope.chartSeries.unshift( {'name':'Sum', 'data': sum , 'type':'line'} );
	};

	$scope.diagram = { cells:[] };
	//________________________________________________________________________
	// function to make the JSON interpreted by joint.js to make the diagram
	$scope.jsonDiag = function () {
		// first loop is on backend facilities, the second on models
		// in a first time build backend facilities and some connection (pool->stock, sp->stock ...)
		// in a second time build a fabrication plant and the list of reactor and their links (fp->reactor->pool)
		$scope.diagram = { cells:[] };

		for ( var i=0 ; i<$scope.backend.length ; i++ )
		{
			var idUnity = $scope.backend[i].id;
			var idBackend = "";

			if ( $scope.backend[i].type == "Stock" ) {
				// Storage
				$scope.diagram.cells.push({
					"type":"erd.Normal","size":{"width":100,"height":50},
					"position":{"x":800,"y":10*i},"angle":0,"embeds":"",
					"z":2,
					"id": idUnity,
					"attrs":{
						"text":{"text": idUnity ,"fill":"#ffffff","letter-spacing":0,"style":{"text-shadow":"1px 0 1px #333333"}},
						".outer":{"stroke":"#fe854f","fill":"#fe8550","filter":{"name":"dropShadow","args":{"dx":0,"dy":2,"blur":2,"color":"#222138"}}}
					}
				});
			}
			if ( $scope.backend[i].type == "Pool" ) {
				idBackend = $scope.backend[i].backendFacility;

				// Pool
				$scope.diagram.cells.push({
					"type":"erd.Derived","size":{"width":100,"height":50},
					"position":{"x":400,"y":10*i},"angle":0,"embeds":"",
					"z":3,
					"id" : idUnity,
					"attrs":{
						"text":{"text": idUnity ,"fill":"#ffffff","letter-spacing":0,"style":{"text-shadow":"1px 0px 1px #333333"}},".inner":{"stroke":"none","fill":"#fca079","display":"block"},
						".outer":{"stroke-dasharray":"3,1","stroke":"#fe854f","fill":"transparent","filter":{"name":"dropShadow","args":{"dx":0,"dy":2,"blur":2,"color":"#222138"}}}
					}
				});
				// Pool -> Stock
				$scope.diagram.cells.push({
					"type":"erd.Line",
					"source":{ "id":idUnity },
					"target":{ "id":idBackend },
					"id":idUnity + "-" + idBackend,
					"router":{"name":"manhattan"},"connector":{"name":"rounded","args":{"radius":10}},"z":1,"attrs":{},

				});
			}
			if ( $scope.backend[i].type == "SeparationPlant" ) {
				// SeparationPlant
				$scope.diagram.cells.push({
					"type":"erd.ISA","size":{"width":100,"height":50},
					"position":{"x":600,"y":10*i},"angle":0,"embeds":"",
					"z":4,
					"id" : idUnity ,
					"attrs":{
						"text":{"text": idUnity ,"fill":"#ffffff","letter-spacing":0,"style":{"text-align":"center","font-size":"0.85em","font-family":"sans-serif","text-shadow":"1px 0px 1px #333333"}},
						"polygon":{"fill":"#fdb664","stroke":"none","filter":{"name":"dropShadow","args":{"dx":0,"dy":2,"blur":1,"color":"#333333"}}}
					}
				});

				// loop on all links
				for ( var j=0 ; j<$scope.backend[i].ivlist.length ; j++ )
				{	
					idBackend = $scope.backend[i].ivlist[j].out;

					// SP -> Stock
					$scope.diagram.cells.push({
						"type":"erd.Line",
						"source":{ "id":idUnity },
						"target":{ "id":idBackend },
						"id":idUnity + "-" + idBackend + "-" + j,
						"router":{"name":"manhattan"},"connector":{"name":"rounded","args":{"radius":10}},"z":1,"attrs":{},
					});
				}
			}
		}

		for ( var i=0 ; i<$scope.models.length ; i++ )
		{
			var idFP = "FP-" + $scope.models[i].id;

			if ( $scope.models[i].type == "Equivalence_Model" )
			{
				// add the fabrication plant of the current model
				$scope.diagram.cells.push({
					"type":"erd.Relationship","size":{"width":80,"height":80},
					"position":{"x":10,"y":10*i},"angle":0,"embeds":"",
					"z":5,
					"id": idFP ,
					"attrs":{
						"text":{"text": idFP , "fill":"#ffffff","letter-spacing":0,"style":{"text-shadow":"1px 0 1px #333333"}},
						".outer":{"fill":"#797d9a","stroke":"none","filter":{"name":"dropShadow","args":{"dx":0,"dy":2,"blur":1,"color":"#333333"}}}
					}
				});

				for ( var j=0 ; j<$scope.models[i].fp.storageList.length ; j++ )
				{
					if ( !$scope.models[i].fp.storageList[j].infinit )
					{
						var key = $scope.models[i].fp.storageList[j].key;
						for ( var k=0 ; k<$scope.models[i].fp.storageList[j].backendFacility.length ; k++ )
						{
							// stock -> fp
							var idStock = $scope.models[i].fp.storageList[j].backendFacility[k].id;
							$scope.diagram.cells.push({
								"type":"erd.Line",
								"source":{ "id":idStock },
								"target":{ "id":idFP },
								"id":idStock + "-" + idFP + "-" + k,
								"router":{"name":"manhattan"},"connector":{"name":"rounded","args":{"radius":10}},"z":1,
								"attrs":{},
								"labels": [{
									"position": 0.2,
									"attrs": {
										"text": { "dy": -8, "text": key , "fill": '#424242' },
										"rect": { "fill": 'none' }
									}
								}]
							});
						}
					}
				}

			}

			// loop on reactor, each itreration add the reactor and links (fp->reactor->pool)
			for ( var j=0 ; j<$scope.models[i].reaclist.length ; j++ )
			{
				var idReactor = $scope.models[i].reaclist[j].id;
				// reactor
				$scope.diagram.cells.push({
					"type":"erd.Entity","size":{"width":150,"height":60},
					"position":{"x":150,"y":50*i+10*j},"angle":0,"embeds":"",
					"z":6,
					"id": idReactor,
					"attrs":{
						"text":{"text": idReactor , "fill":"#ffffff","letter-spacing":0,"style":{"text-shadow":"1px 0 1px #333333"}},
						".outer, .inner":{"fill":"#31d0c6","stroke":"none","filter":{"name":"dropShadow","args":{"dx":0.5,"dy":2,"blur":2,"color":"#333333"}}}
					}
				});

				if ( $scope.models[i].type == "Equivalence_Model" )
				{
					// fp -> reactor
					$scope.diagram.cells.push({
						"type":"erd.Line",
						"source":{ "id":idFP },
						"target":{ "id":idReactor },
						"id":idFP + "-" + idReactor,
						"router":{"name":"manhattan"},"connector":{"name":"rounded","args":{"radius":10}},"z":1,"attrs":{},
					});
				}
				// reactor -> pool
				$scope.diagram.cells.push({
					"type":"erd.Line",
					"source":{ "id":idReactor },
					"target":{ "id":$scope.models[i].reaclist[j].backendFacility },
					"id":idReactor + "-" + $scope.models[i].reaclist[j].backendFacility,
					"router":{"name":"manhattan"},"connector":{"name":"rounded","args":{"radius":10}},"z":1,"attrs":{},
					//"attrs":{ '.marker-target': { d: 'M 10 0 L 0 5 L 10 10 z' } } // for add arrowhead but... voilà la doc est pourrie...
				});
			}

		}
	};

	// diagram make once (not at each click on "Make diagram")
	$scope.erd = joint.shapes.erd;
	$scope.graph = new joint.dia.Graph;
	$scope.paper = new joint.dia.Paper({
		el: $('#diagram'),
		width: document.getElementById('diagram').offsetWidth,
		height: 0,
		model: $scope.graph,
		gridSize: 1,
		linkPinning: false,
		linkConnectionPoint: joint.util.shapePerimeterConnectionPoint,
	});

	//________________________________________________________________________
	// function to call convertion scenrio into JSON for Joint.js and update diagram
	$scope.makeDiagram = function () {
		$scope.paper.setDimensions( document.getElementById('diagram').offsetWidth , 400 );
		$scope.jsonDiag();
		$scope.graph.fromJSON($scope.diagram);
	};

	// code pour l'envoie de données sur le serveur. Il s'agit des mêmes
	// données que lors de la sauvegarde au format JSON qui sont envoyées puis
	// récupérées par la page run.py
	$scope.send = 'no';
	//________________________________________________________________________
	// watch change of value of send variable
	$scope.$watch('send', function () {
		if ( $scope.send == 'yes' )
		{
			$scope.sendToSrv();
		}
	});

	//________________________________________________________________________
	// function to send data to server
	$scope.sendToSrv = function () {
		// data to send to the server
		var raw_data = {
			config  : $scope.config,
			backend : $scope.backend,
			models  : $scope.models,
		};

		// send raw_data as a string with method POST to the page run.py
		$http.post(
			'/run.py',
			{data:JSON.stringify(raw_data)},
			{
				headers: { 'Content-Type': 'application/x-www-form-urlencoded; charset=UTF-8'},
				transformRequest: transform // <- I don't understand that but it works...
			}
		).then(
			function(){ $scope.send = 'ok';   }, // btn become green
			function(){ $scope.send = 'fail'; }  // btn become red
		);
	};

}]);

var transform = function(data){
	return $.param(data);
}
