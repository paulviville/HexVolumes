import CMap0 from './CMapJS/CMap/CMap0.js';
import Graph from './CMapJS/CMap/Graph.js';
import IncidenceGraph from './CMapJS/CMap/IncidenceGraph.js';
import * as THREE from './CMapJS/Libs/three.module.js';
import Renderer from './CMapJS/Rendering/Renderer.js';

export default function MeshHandler (mesh, params = {}) {
	
	// const tetrahedra0 = new Graph()
	// const tetrahedra1 = new Graph()
	// const tetRender0 = new Renderer(tetrahedra0)
	// const tetRender1 = new Renderer(tetrahedra1)

	// this.initializeTetrahedra = function() {
	// 	const pos = mesh.getAttribute(mesh.vertex, "position");

	// 	const pos0 = tetrahedra0.addAttribute(tetrahedra0.vertex, "position");
	// 	const pos1 = tetrahedra1.addAttribute(tetrahedra1.vertex, "position");

	// 	for(let i = 0; i < 8; ++i) {
	// 		tetrahedra0.addVertex();
	// 		tetrahedra1.addVertex();

	// 		pos0[i] = pos[i].clone();
	// 		pos1[i] = pos[i].clone();
	// 	}

	// 	tetrahedra0.connectVertices(1, 3);
	// 	tetrahedra0.connectVertices(4, 3);
	// 	tetrahedra0.connectVertices(6, 3);
	// 	tetrahedra0.connectVertices(6, 1);
	// 	tetrahedra0.connectVertices(4, 1);
	// 	tetrahedra0.connectVertices(4, 3);

	// 	tetrahedra1.connectVertices(2, 7);
	// 	tetrahedra1.connectVertices(2, 5);
	// 	tetrahedra1.connectVertices(2, 0);
	// 	tetrahedra1.connectVertices(7, 5);
	// 	tetrahedra1.connectVertices(0, 5);
	// 	tetrahedra1.connectVertices(0, 7);

	// }



	const renderer = new Renderer(mesh);
	



	const vertexColor = new THREE.Color(params.vertexColor || 0x4EA6BA);
	const edgeColor = new THREE.Color(params.edgeColor || 0x0A0A20);
	const faceColor = new THREE.Color(params.faceColor || 0x66AABB);
	let vertexSize = params.vertexSize || 0.01; 
	let edgeSize = params.edgeSize || 1.5; 

	let parentObject;
	let verticesMesh, edgesMesh, facesMesh;
	this.initialize = function (params = {}) {
		console.log(params)
		if(params.vertices) {
			renderer.vertices.create({size: vertexSize, color: vertexColor}); 
			verticesMesh = renderer.vertices.mesh;
		}
		if(params.edges) {
			renderer.edges.create({size: edgeSize, color: edgeColor}); 
			edgesMesh = renderer.edges.mesh;
		}
		if(params.faces) {
			renderer.faces.create({color: faceColor, side: THREE.DoubleSide}); 
			facesMesh = renderer.faces.mesh;
		}
	};

	this.addMeshesTo = function (parent) {
		parentObject = parentObject || parent;
		if(verticesMesh) renderer.vertices.addTo(parent);
		if(edgesMesh) renderer.edges.addTo(parent);
		if(facesMesh) renderer.faces.addTo(parent);
	};

	this.setVertexColor = function (color) {
		if(verticesMesh) {
			vertexColor.setHex(color);
			verticesMesh.material.color.setHex(color);
			verticesMesh.material.needsUpdate = true;
		}
	};

	this.setEdgeColor = function (color) {
		if(edgesMesh) {
			edgeColor.setHex(color);
			edgesMesh.material.color.setHex(color);
			edgesMesh.material.needsUpdate = true;
		}
	};

	this.setFaceColor = function (color) {
		if(facesMesh) {
			faceColor.setHex(color);
			facesMesh.material.color.setHex(color);
			facesMesh.material.needsUpdate = true;
		}
	};

	this.resizeVertices = function(size) {
		vertexSize = size;
		renderer.vertices.resize(size);
		this.updateVertices();
	}

	this.resizeEdges = function(size) {
		edgeSize = size;
		renderer.edges.resize(size);
		this.updateEdges();
	}

	this.updateVertices = function() {
		const visible = verticesMesh.visible
		renderer.vertices.update();
		verticesMesh = renderer.vertices.mesh;
		verticesMesh.visible = visible
	};

	this.updateEdges = function() {
		const visible = edgesMesh.visible;
		renderer.edges.update();
		edgesMesh = renderer.edges.mesh;
		edgesMesh.visible = visible;
	};

	this.updateFaces = function() {
		const visible = facesMesh.visible;
		renderer.faces.update();
		facesMesh = renderer.faces.mesh;
		facesMesh.visible = visible;
	}

	this.updateMeshes = function () {
		if(verticesMesh) {
			this.updateVertices();
		}
		if(edgesMesh) {
			this.updateEdges();
		}
		if(facesMesh) {
			this.updateFaces();
		}
	}

	this.vertexVisibility = function (visible) {
		if(verticesMesh)
			verticesMesh.visible = visible
		else if(visible) {
			renderer.vertices.create({size: vertexSize, color: vertexColor}); 
			verticesMesh = renderer.vertices.mesh;
			renderer.vertices.addTo(parentObject);
		}
	}

	this.edgeVisibility = function (visible) {
		if(edgesMesh)
			edgesMesh.visible = visible
		else if(visible) {
			renderer.edges.create({size: edgeSize, color: edgeColor}); 
			edgesMesh = renderer.edges.mesh;
			renderer.edges.addTo(parentObject);
		}
	}

	this.faceVisibility = function (visible) {
		if(facesMesh)
			facesMesh.visible = visible
		else if(visible) {
			renderer.faces.create({color: faceColor}); 
			facesMesh = renderer.faces.mesh;
			renderer.faces.addTo(parentObject);
		}
	}


	this.delete = function () {
		if(verticesMesh) renderer.vertices.delete();
		if(edgesMesh) renderer.edges.delete();
		if(facesMesh) renderer.faces.delete();
	}

	this.moveVertex = function (id = 0, x = 0, y = 0, z = 0) {
		const pos = mesh.getAttribute(mesh.vertex, "position");
		pos[id].add(new THREE.Vector3(x, y, z));
	}

	const subGrid = new IncidenceGraph();
	const subGridPos = subGrid.addAttribute(subGrid.vertex, "position");
	subGrid.createEmbedding(subGrid.vertex);
	const subGridRenderer = new Renderer(subGrid);
	const X = 8; const Y = 8; const Z = 8;
	const XY = X * Y;
	const subGridSize = X * Y * Z;
	this.createSubdivision = function (parent) {

		for(let i = 0; i < subGridSize; ++i) {
			subGrid.addVertex();
			subGridPos[i] = new THREE.Vector3();
		}

		this.updateSubdivision()
		this.exportSubdivision()

		// subGridRenderer.vertices.create({size: 0.0025});
		// subGridRenderer.vertices.addTo(parent);
	}

	this.updateSubdivision = function () {
		const pos = mesh.getAttribute(mesh.vertex, "position");
		const p0 = pos[0];
		const p1 = pos[1];
		const p2 = pos[2];
		const p3 = pos[3];
		const p4 = pos[4];
		const p5 = pos[5];
		const p6 = pos[6];
		const p7 = pos[7];

		// const p06 = new THREE.Vector3;
		// const p35 = new THREE.Vector3;
		// const p = new THREE.Vector3;
		// const p = new THREE.Vector3;

		// console.log()

		const p01 = new THREE.Vector3;
		const p32 = new THREE.Vector3;
		const p45 = new THREE.Vector3;
		const p76 = new THREE.Vector3;

		for(let i = 0; i < subGridSize; ++i) {
			// subGridPos[i].copy(pos[i]);
			const x = (i % X)/(X-1);
			const y = (Math.floor((i % XY) / X))/(Y-1);
			const z = (Math.floor(i / XY))/(Z-1);
			// console.log(x / (X-1), y/ (Y-1), z/(Z-1));

			p01.lerpVectors(p0, p1, x);
			p32.lerpVectors(p3, p2, x);
			p76.lerpVectors(p7, p6, x);
			p45.lerpVectors(p4, p5, x);
			
			p01.lerp(p32, y);
			p45.lerp(p76, y);
		
			p01.lerp(p45, z);

			subGridPos[i].copy(p01);

		}
		// subGridRenderer.vertices.update();
	
	}

	function toXYZ (i) {
		return new THREE.Vector3(
			(i % X),
			Math.floor((i % XY) / X),
			Math.floor(i / XY)
		);
	}

	function toI(x, y, z) {
		return z * XY + y * X + x; 
	}

	this.exportSubdivision = function() {
		let str = `Dimension`
		str += `3`
		str += `Vertices`
		str += `${subGrid.nbCells(subGrid.vertex)}\n`;

		subGrid.foreach(subGrid.vertex, v => {
			let vpos = subGridPos[v];
			str += `${vpos.x} ${vpos.y} ${vpos.z} 0\n`;
		});

		str += `Quads`
		str +=`0\n`;

		str += `Hexahedra`
		str += `${(X-1)*(Y-1)*(Z-1)}`;
		for(let x = 0; x < X - 1; ++x) {
			for(let y = 0; y < Y - 1; ++y) {
				for(let z = 0; z < X - 1; ++z) {
					let i = toI(x, y, z);
					str += `${toI(x, y, z)} `;
					str += `${toI(x+1, y, z)} `;
					str += `${toI(x+1, y+1, z)} `;
					str += `${toI(x, y+1, z)} `;
					str += `${toI(x, y, z+1)} `;
					str += `${toI(x+1, y, z+1)} `;
					str += `${toI(x+1, y+1, z+1)} `;
					str += `${toI(x, y+1, z+1)} `;
					str += `0\n`
					// volume += avgVolume(
					// 	subGridPos[toI(x, y, z)], 
					// 	subGridPos[toI(x+1, y, z)], 
					// 	subGridPos[toI(x+1, y+1, z)], 
					// 	subGridPos[toI(x, y+1, z)], 
					// 	subGridPos[toI(x, y, z+1)], 
					// 	subGridPos[toI(x+1, y, z+1)], 
					// 	subGridPos[toI(x+1, y+1, z+1)], 
					// 	subGridPos[toI(x, y+1, z+1)])
				}
			}
		}

		str += `End`;
		console.log(str);
	}


	function tetDecompVolume(p0, p1, p2, p3, p4, p5, p6, p7) {
		let volume = 0;
		const t0 = new THREE.Vector3;
		const t1 = new THREE.Vector3;
		const t2 = new THREE.Vector3;

		t0.subVectors(p3, p0);
		t1.subVectors(p1, p0);
		t2.subVectors(p4, p0);
		volume += t0.cross(t1).dot(t2);

		t0.subVectors(p1, p2);
		t1.subVectors(p3, p2);
		t2.subVectors(p6, p2);
		volume += t0.cross(t1).dot(t2);

		t0.subVectors(p1, p5);
		t1.subVectors(p6, p5);
		t2.subVectors(p4, p5);
		volume += t0.cross(t1).dot(t2);

		t0.subVectors(p3, p7);
		t1.subVectors(p4, p7);
		t2.subVectors(p6, p7);
		volume += t0.cross(t1).dot(t2);

		t0.subVectors(p3, p1);
		t1.subVectors(p6, p1);
		t2.subVectors(p4, p1);
		volume += t0.cross(t1).dot(t2);


		return volume / 6;
	}

	this.computeVolumeTetDecomp = function() {
		const pos = mesh.getAttribute(mesh.vertex, "position");
		return tetDecompVolume(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], pos[6], pos[7]);
	}

	this.computeVolumeTetDecompAvg = function() {
		const pos = mesh.getAttribute(mesh.vertex, "position");
		return (tetDecompVolume(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], pos[6], pos[7]) + tetDecompVolume(pos[1], pos[2], pos[3], pos[0], pos[5], pos[6], pos[7], pos[4]) )/2;
	}

	function avgVolume(p0, p1, p2, p3, p4, p5, p6, p7) {
		let v01 = p1.clone().sub(p0)
		let v03 = p3.clone().sub(p0)
		let v04 = p4.clone().sub(p0)
		let v12 = p2.clone().sub(p1)
		let v15 = p5.clone().sub(p1)
		let v23 = p3.clone().sub(p2)
		let v26 = p6.clone().sub(p2)
		let v37 = p7.clone().sub(p3)
		let v45 = p5.clone().sub(p4)
		let v47 = p7.clone().sub(p4)
		let v56 = p6.clone().sub(p5)
		let v67 = p7.clone().sub(p6)

		let volume = 0;
		const t = new THREE.Vector3;

		volume += t.crossVectors(v03, v01).dot(v04);
		volume += t.crossVectors(v12, v01).dot(v15);
		volume += t.crossVectors(v23, v12).dot(v26);
		volume += t.crossVectors(v23, v03).dot(v37);
		volume += t.crossVectors(v45, v04).dot(v47);
		volume += t.crossVectors(v15, v56).dot(v45);
		volume += t.crossVectors(v26, v67).dot(v56);
		volume += t.crossVectors(v47, v37).dot(v67);

		return volume / 8;
	}

	this.computeVolumeAvg = function() {
		const pos = mesh.getAttribute(mesh.vertex, "position");

		return avgVolume(pos[0], pos[1], pos[2], pos[3], pos[4], pos[5], pos[6], pos[7]);


	}

	this.computeSubdivisionVolume = function() {
		let volume = 0;
		for(let x = 0; x < X - 1; ++x) {
			for(let y = 0; y < Y - 1; ++y) {
				for(let z = 0; z < X - 1; ++z) {
					let i = toI(x, y, z);
					volume += avgVolume(
						subGridPos[toI(x, y, z)], 
						subGridPos[toI(x+1, y, z)], 
						subGridPos[toI(x+1, y+1, z)], 
						subGridPos[toI(x, y+1, z)], 
						subGridPos[toI(x, y, z+1)], 
						subGridPos[toI(x+1, y, z+1)], 
						subGridPos[toI(x+1, y+1, z+1)], 
						subGridPos[toI(x, y+1, z+1)])
				}
			}
		}

		return volume;
	}


	let velocity;
	let position;
	let prevPosition;
	let invMass; 
	let hexVolumes;
	let edgeLength;
	this.initializeSimulation = function() {
		mesh.createEmbedding(mesh.edge)
		mesh.createEmbedding(mesh.volume)

		mesh.setEmbeddings(mesh.edge);
		mesh.setEmbeddings(mesh.volume);

		velocity = mesh.addAttribute(mesh.vertex, "velocity");
		prevPosition = mesh.addAttribute(mesh.vertex, "prevPosition");
		position = mesh.getAttribute(mesh.vertex, "position");
		hexVolumes = mesh.addAttribute(mesh.volume, "hexVolumes");
		edgeLength = mesh.addAttribute(mesh.edge, "edgeLength");
		invMass = mesh.addAttribute(mesh.vertex, "invMass");

		mesh.foreach(mesh.vertex, vd => {
			let m = 0;
			mesh.foreachIncident(mesh.volume, mesh.vertex, vd, wd => {
				if(mesh.isBoundary(wd))
					return;

				m += 1/8;
			});
			// const m = mesh.degree(mesh.vertex, vd);
			const vid = mesh.cell(mesh.vertex, vd);
			invMass[vid] = 1/m;
			velocity[vid] = new THREE.Vector3;
			prevPosition[vid] = new THREE.Vector3;
		});

		const edgeVector = new THREE.Vector3;
		mesh.foreach(mesh.edge, ed => {
			edgeVector.subVectors(
				position[mesh.cell(mesh.vertex, ed)],
				position[mesh.cell(mesh.vertex, mesh.phi2[ed])]
			)

			edgeLength[mesh.cell(mesh.edge, ed)] = edgeVector.length();
			// edge
		});


		const vids = new Array(8);
		const vpos = [
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
		];

		mesh.foreach(mesh.volume, wd => {
			if(mesh.isBoundary(wd))
				return;
			console.log(wd)
			// edge
			let d = wd;
			vids[0] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[1] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[2] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[3] = mesh.cell(mesh.vertex, d);

			d = mesh.phi1[mesh.phi1[mesh.phi2[wd]]];
			vids[4] = mesh.cell(mesh.vertex, d);
			d = mesh.phi2[d];			
			vids[5] = mesh.cell(mesh.vertex, d);
			d = mesh.phi_1[d];			
			vids[6] = mesh.cell(mesh.vertex, d);
			d = mesh.phi_1[d];			
			vids[7] = mesh.cell(mesh.vertex, d);
			for(let i = 0; i < 8; ++i) {
				vpos[i].copy(position[vids[i]]);
			} 
			const vol = -avgVolume(vpos[0], vpos[1], vpos[2], vpos[3], vpos[4], vpos[5], vpos[6], vpos[7]);
			hexVolumes[mesh.cell(mesh.volume, wd)] = vol;
		});

		console.log(invMass)
		console.log(edgeLength)
		console.log(hexVolumes)
	}

	function preSolve(dt) {
		mesh.foreach(mesh.vertex, vd => {
			const vid = mesh.cell(mesh.vertex, vd);
			prevPosition[vid].copy(position[vid]);
			position[vid].addScaledVector(velocity[vid], dt);
		});
	}


	function solveEdges(dt, compliance) {
		const alpha = compliance / (dt * dt);
		const edges = mesh.cache(mesh.edge);
		edges.sort((a, b) => Math.random() - 0.5);

		const grad = new THREE.Vector3;
		edges.forEach(ed => {
			const vid0 = mesh.cell(mesh.vertex, ed);
			const vid1 = mesh.cell(mesh.vertex, mesh.phi2[ed]);

			const w0 = invMass[vid0];
			const w1 = invMass[vid1];
			const w = w0 + w1;

			grad.subVectors(position[vid0], position[vid1]);
			const length = grad.length();

			if(length == 0.0)
				return false;

			grad.multiplyScalar(1/length);
			const c = length - edgeLength[mesh.cell(mesh.edge, ed)];
			const s = -c / (w + alpha);
			position[vid0].addScaledVector(grad, s * w0);
			position[vid1].addScaledVector(grad, -s * w1);
		});
	}

	function solveVolumes(dt, compliance) {
		const alpha = compliance / (dt * dt);
		
		// const volumes = mesh.cache(mesh.volume, d => {!mesh.isBoundary(d)});
		const volumes = mesh.cache(mesh.volume);
		const vids = new Array(8);
		const vpos = [
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
		];
		const grads = [
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
			new THREE.Vector3,
		];
		const temp0 = new THREE.Vector3
		const temp1 = new THREE.Vector3
		const center = new THREE.Vector3;
		const dir = new THREE.Vector3;

		volumes.forEach(wd => {
			if(mesh.isBoundary(wd))
				return;


			let d = wd;
			vids[0] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[1] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[2] = mesh.cell(mesh.vertex, d);
			d = mesh.phi1[d];			
			vids[3] = mesh.cell(mesh.vertex, d);

			d = mesh.phi1[mesh.phi1[mesh.phi2[wd]]];
			vids[4] = mesh.cell(mesh.vertex, d);
			d = mesh.phi2[d];			
			vids[5] = mesh.cell(mesh.vertex, d);
			d = mesh.phi_1[d];			
			vids[6] = mesh.cell(mesh.vertex, d);
			d = mesh.phi_1[d];			
			vids[7] = mesh.cell(mesh.vertex, d);

			center.set(0, 0, 0);
			for(let i = 0; i < 8; ++i) {
				vpos[i].copy(position[vids[i]]);
				center.add(vpos[i]);
			} 
			center.multiplyScalar(1/8);



			temp0.copy(vpos[3]).sub(vpos[1]);
			temp1.copy(vpos[4]).sub(vpos[1]);
			grads[0].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[5]).sub(vpos[0]);
			temp1.copy(vpos[2]).sub(vpos[0]);
			grads[1].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[6]).sub(vpos[1]);
			temp1.copy(vpos[3]).sub(vpos[1]);
			grads[2].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[7]).sub(vpos[2]);
			temp1.copy(vpos[0]).sub(vpos[2]);
			grads[3].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[7]).sub(vpos[0]);
			temp1.copy(vpos[5]).sub(vpos[0]);
			grads[4].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[4]).sub(vpos[1]);
			temp1.copy(vpos[6]).sub(vpos[1]);
			grads[5].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[5]).sub(vpos[2]);
			temp1.copy(vpos[7]).sub(vpos[2]);
			grads[6].crossVectors(temp0, temp1).multiplyScalar(1/8);

			temp0.copy(vpos[6]).sub(vpos[3]);
			temp1.copy(vpos[4]).sub(vpos[3]);
			grads[7].crossVectors(temp0, temp1).multiplyScalar(1/8);

			let w = 0;
			for(let i = 0; i < 8; ++i) {
			// 	vpos[i].sub(center).multiplyScalar(1.0/8.0);
				w += invMass[vids[i]] * grads[i].lengthSq();
			}

			if(w == 0)
				return

			const vol = -avgVolume(vpos[0], vpos[1], vpos[2], vpos[3], vpos[4], vpos[5], vpos[6], vpos[7]);
			const restVol = hexVolumes[mesh.cell(mesh.volume, wd)];

			const c = vol - restVol;
			const s = -c / (w + alpha);
			// console.log(c, s, w)
			for(let i = 0; i < 8; ++i) {
				position[vids[i]].addScaledVector(grads[i], s * invMass[vids[i]]);
			} 

		});

	}

	function postSolve(dt) {
		mesh.foreach(mesh.vertex, vd => {
			const vid = mesh.cell(mesh.vertex, vd);
			velocity[vid].subVectors(position[vid], prevPosition[vid]).multiplyScalar(1/dt).clampLength(0, 0.1).multiplyScalar(0.999);
		});
	}


	this.simulate = function (dt, complianceEdge, complianceVolume) {
		const iterations = 10;
		for(let i = 0; i < iterations; ++i) { 
			preSolve(dt/iterations);
			solveEdges(dt/iterations, complianceEdge);
			solveVolumes(dt/iterations, complianceVolume);
			postSolve(dt/iterations);
		}

		// console.log(edgeLength)
	}
}