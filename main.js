import * as THREE from './CMapJS/Libs/three.module.js';
import * as Meshes from './meshes.js';

import { loadCMap2 } from './CMapJS/IO/SurfaceFormats/CMap2IO.js';
import { loadCmap3 } from './CMapJS/IO/Volumes_Formats/CMap3IO.js';

import { OrbitControls } from './CMapJS/Libs/OrbitsControls.js';
import MeshHandler from './MeshHandler.js';

import { GUI } from './CMapJS/Libs/dat.gui.module.js';

import catmullClark from './CMapJS/Modeling/Subdivision/Surface/CatmullClark.js';
import loop from './CMapJS/Modeling/Subdivision/Surface/Loop.js';
import dooSabin from './CMapJS/Modeling/Subdivision/Surface/DooSabin.js';
import sqrt2 from './CMapJS/Modeling/Subdivision/Surface/Sqrt2.js';
import sqrt3 from './CMapJS/Modeling/Subdivision/Surface/Sqrt3.js';
import butterfly from './CMapJS/Modeling/Subdivision/Surface/Butterfly.js';

const scene = new THREE.Scene();
scene.background = new THREE.Color(0xAAAAAA);

const camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000.0);
camera.position.set(0, 0, 2);

const renderer = new THREE.WebGLRenderer();
renderer.setSize( window.innerWidth, window.innerHeight );
document.body.appendChild( renderer.domElement );

let controls = new OrbitControls(camera, renderer.domElement)


window.addEventListener('resize', function() {
    const width = window.innerWidth;
    const height = window.innerHeight;
    renderer.setSize(width, height);
    camera.aspect = width / height;
    camera.updateProjectionMatrix();
});

let ambientLight = new THREE.AmbientLight(0xFFFFFF, 0.5);
scene.add(ambientLight);
let pointLight = new THREE.PointLight(0xFFFFFF, 1);
pointLight.position.set(10,8,5);
scene.add(pointLight);

let meshHandler;
let cmap;

const settings = new (function() {
	this.showVertex = true;
	this.vertexSize = 0.01;
	this.vertexColor = 0x4EA6BA;
	this.updateVertexColor = function (color) {meshHandler.setVertexColor(color)};
	this.updateVertexVisibility = function (visible) {meshHandler.vertexVisibility(visible)};
	this.showEdge = true;
	this.edgeSize = 1.5;
	this.edgeColor = 0x0A0A20;
	this.updateEdgeColor = function (color) {meshHandler.setEdgeColor(color)};
	this.updateEdgeVisibility = function (visible) {meshHandler.edgeVisibility(visible)};
	this.showFace = false;
	this.faceColor = 0x66AABB;
	this.updateFaceColor = function (color) {meshHandler.setFaceColor(color)};
	this.updateFaceVisibility = function (visible) {meshHandler.faceVisibility(visible)};

	this.vertexResize = function (size) {meshHandler.resizeVertices(size)};
	this.edgeResize = function (size) {meshHandler.resizeEdges(size)};


	this.mesh = 'cube';
	
	this.volume = 1;
	this.volumeTet = 1;
	this.volumeSubdiv = 1;

	this.simulate = false;
	
});

function loadMesh (mesh) {
	// cmap = loadCMap2('off', Meshes[mesh + '_off']);
	cmap = loadCmap3("mesh", Meshes.gridMesh)
	if(meshHandler)
		meshHandler.delete();
	meshHandler = new MeshHandler(cmap, {
		vertexColor: settings.vertexColor,
		edgeColor: settings.edgeColor,
		faceColor: settings.faceColor,
		vertexSize: settings.vertexSize,
		edgeSize: settings.edgeSize,

	});
	meshHandler.initialize({
		vertices: settings.showVertex,
		edges: settings.showEdge,
		faces: settings.showFace,
	});
	meshHandler.addMeshesTo(scene);

	// meshHandler.createSubdivision(scene);
	// meshHandler.updateSubdivision();
	settings.volumeTet = meshHandler.computeVolumeTetDecompAvg();
	settings.volume = meshHandler.computeVolumeAvg();
	// settings.volumeSubdiv = meshHandler.computeSubdivisionVolume();
	meshHandler.initializeSimulation();
}


window.move = function(id, x, y, z) {
	meshHandler.moveVertex(id, x, y, z);
	meshHandler.updateMeshes();
	// meshHandler.updateSubdivision();
	settings.volume = meshHandler.computeVolumeAvg();
	settings.volumeTet = meshHandler.computeVolumeTetDecompAvg() / settings.volume;
	// settings.volumeSubdiv = meshHandler.computeSubdivisionVolume()/ settings.volume;
	// meshHandler.computeVolume();
}

const gui = new GUI({autoPlace: true, hideable: true});
const settingsFolder = gui.addFolder("Settings");
settingsFolder.add(settings, 'showVertex').onChange(settings.updateVertexVisibility);
settingsFolder.add(settings, 'showEdge').onChange(settings.updateEdgeVisibility);
settingsFolder.add(settings, 'showFace').onChange(settings.updateFaceVisibility);
settingsFolder.addColor(settings, 'vertexColor').onChange(settings.updateVertexColor);
settingsFolder.addColor(settings, 'edgeColor').onChange(settings.updateEdgeColor);
settingsFolder.addColor(settings, 'faceColor').onChange(settings.updateFaceColor);
settingsFolder.add(settings, 'vertexSize').min(0.001).max(0.1).step(0.001).onChange(settings.vertexResize);
settingsFolder.add(settings, 'edgeSize').min(0.2).max(5).step(0.05).onChange(settings.edgeResize);

gui.add(settings, 'volume').step(0.0001).listen();
gui.add(settings, 'volumeTet').step(0.0001).listen();
gui.add(settings, 'volumeSubdiv').step(0.0001).listen();
gui.add(settings, 'simulate');

gui.add(settings, 'mesh', ['tetrahedron', 'cube', 'octahedron', 'dodecahedron', 'icosahedron']).onChange(loadMesh);
loadMesh(settings.mesh);

function update()
{
	if(settings.simulate) {
		meshHandler.simulate(0.02, 1, 1);
		meshHandler.updateMeshes();
		settings.volume = meshHandler.computeVolumeAvg();
		settings.volumeTet = meshHandler.computeVolumeTetDecompAvg() / settings.volume;
	}
}

function render()
{
	renderer.render(scene, camera);
}

function mainloop()
{
	update();
    render();
    requestAnimationFrame(mainloop);
}

mainloop();