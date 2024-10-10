function rnd(n) {
	return Math.floor(Math.random() * n);
}
rn = [[7, 11, 20, 15, 5, 13, 20, 13, 19, 6, 10, 16, 4, 18, 4, 19, 4, 5, 4, 0], [4, 17, 8, 11, 20, 12, 16, 4, 8, 1, 10, 10, 9, 13, 8, 14, 10, 5, 8, 5], [13, 14, 5, 6, 15, 2, 2, 20, 2, 11, 19, 11, 8, 17, 16, 15, 9, 0, 20, 12]];
// Use predefined colors or continue generating random colors
const predefinedColors = [
	[255, 87, 51],   // Red-Orange
	[51, 255, 87],   // Green
	[51, 87, 255],   // Blue
	[255, 51, 168],  // Pink
	[255, 189, 51],  // Yellow-Orange
	[51, 255, 209],  // Aqua
	[168, 51, 255],  // Purple
	[51, 255, 138],  // Mint Green
	[255, 51, 51],   // Bright Red
	[51, 138, 255],  // Sky Blue
	[255, 209, 51],  // Mustard Yellow
	[138, 255, 51],  // Lime Green
	[255, 51, 255],  // Magenta
	[51, 255, 215],  // Light Turquoise
	[168, 255, 51],  // Yellow-Green
	[209, 51, 255],  // Violet
	[255, 51, 136],  // Rose
	[51, 255, 161],  // Pale Green
	[255, 138, 51],  // Burnt Orange
	[51, 209, 255]   // Cyan
];
var d; // dimension of embedding space
var dd = 0; // dimension *as currently displayed on the table*. may be different.
var vec = Array(Array(d),Array(d));
var point = Array(dd);

function dot(v1, v2) // d
{
	var s = 0, i;
	for (i = 0; i < d; i++)
		s += v1[i] * v2[i];

	return s;
}

function mdot(m, v1) // d-2
{
	var v2 = Array(d - 2);
	var i, j;

	for (i = 0; i < d - 2; i++) {
		v2[i] = 0;
		for (j = 0; j < d - 2; j++)
			v2[i] += m[i][j] * v1[j];
	}
	return v2;
}


var hash;
var pstack;

function insertpoint(p) // test (and possibly put a point in the stack) if not already tested
{
	if (typeof (hash[p]) == 'undefined') {
		if (convex(p)) {
			hash[p] = 1; // 1 means on the surface but neighbors not tested yet
			pstack.push(p.slice()); // slice() cause we need a copy of the Array
		}
		else
			hash[p] = 0; // 0 means not on the surface
	}

}


var p0, vproj, w, ww;

var orthogproj = 1;

function proj(p) // projection onto plane
{
	var dp = Array(d);
	var i;

	for (i = 0; i < d; i++)
		dp[i] = -p0[i] + p[i];

	return Array(dot(vproj[0], dp), dot(vproj[1], dp));
}

function convex(p) // test if a point is close enough to the surface to be "in"
{
	var dp = Array(d);
	var pr = Array(d - 2), q;
	var i, j, k, l;

	for (i = 0; i < d; i++)
		dp[i] = -p0[i] + p[i];

	for (i = 0; i < d - 2; i++) {
		pr[i] = 0.;
		for (j = 0; j < d; j++)
			pr[i] += f[i][j] * dp[j]; // projection onto (d-2)-space
	}

	for (i = 0; i < d - 1; i++)
		for (j = i + 1; j < d; j++) {
			q = mdot(w[i][j], pr);
			for (k = 0; k < 4; k++) // for at least one (d-2)-face of the lattice
			{
				l = 0; while ((l < d - 2) && (Math.abs(q[l] + ww[i][j][k][l]) <= 0.5)) l++; // we are inside
				if (l == d - 2) return true;
			}
		}

	return false;
}

var A, B, scale;

function check(p) {
	var pr = proj(p);

	return (pr[0] >= -A - 1) && (pr[0] <= A + 1) && (pr[1] >= -B - 1) && (pr[1] <= B + 1); // overkill but at least we don't get glitches! (we don't?)
}

var canvas, ctx;
var cx, cy; // canvas sizes

function pivot(m) {
	// primitive gauss pivot (size d-2) -> produces inverse of m
	var n = d - 2;
	var i, j, k;
	var mm;
	var coeff;

	mm = Array(n);

	for (i = 0; i < n; i++) {
		mm[i] = Array(n);
		for (j = 0; j < n; j++)
			if (i == j) mm[i][j] = 1.; else mm[i][j] = 0.;
	}

	for (i = 0; i < n - 1; i++) {
		for (j = i + 1; j < n; j++) {
			coeff = m[j][i] / m[i][i];  // better choice???
			for (k = i; k < n; k++)
				m[j][k] -= coeff * m[i][k];
			for (k = 0; k < n; k++)
				mm[j][k] -= coeff * mm[i][k];
		}
	}

	for (i = n - 1; i >= 0; i--) {
		for (k = 0; k < n; k++) {
			for (j = i + 1; j < n; j++)
				mm[i][k] -= m[i][j] * mm[j][k];
			mm[i][k] /= m[i][i];
		}
	}
	return mm;
}


var precomputed = false, computed = false; // booleans

function precompute() {

	d = dd; // update dimensionality
	
	var v = Array(Array(d), Array(d));
	vproj = Array(Array(d), Array(d));
	p0 = Array(d);
	w = Array(d); for (i = 0; i < d; i++) w[i] = Array(d);
	ww = Array(d); for (i = 0; i < d; i++) ww[i] = Array(d);

	v = vec;
	p0 = point;


	// do some initial computations...
	// orthonormalize v[0], v[1] to simplify calculations
	var a = 1 / Math.sqrt(dot(v[0], v[0]));
	for (k = 0; k < d; k++)
		v[0][k] *= a;
	a = -dot(v[0], v[1]);
	for (k = 0; k < d; k++)
		v[1][k] += a * v[0][k];
	a = 1 / Math.sqrt(dot(v[1], v[1]));
	for (k = 0; k < d; k++)
		v[1][k] *= a;

	if (orthogproj) {
		vproj = v; // that's easy. though we should probably copy to table too...
	}
	

	// find covectors defining the v[0], v[1] plane (the 2 steps could be combined into Gram-Schmidt...)
	f = Array(d - 2);
	for (i = 0; i < d - 2; i++) {
		f[i] = Array(d);
		for (j = 0; j < d; j++)
			if (i == j) f[i][j] = 1.; else f[i][j] = -0.01 * Math.random(); // arbitrary choice...
		for (k = 0; k < 2; k++) {
			a = -dot(f[i], v[k]);
			for (j = 0; j < d; j++)
				f[i][j] += a * v[k][j];
		}
	}


	// precompute inverses of sub (d-2)x(d-2) matrices
	var k, kk, l;
	var vec1 = Array(d - 2), vec2 = Array(d - 2);
	var mat = Array(d - 2); for (i = 0; i < d - 2; i++) mat[i] = Array(d - 2);
	var tmp1, tmp2
	for (i = 0; i < d - 1; i++)
		for (j = i + 1; j < d; j++) {
			for (l = 0; l < d - 2; l++) {
				kk = -1;
				for (k = 0; k < d; k++) {
					if (k == i)
						vec1[l] = 0.5 * f[l][k];
					else if (k == j)
						vec2[l] = 0.5 * f[l][k];
					else
						mat[l][++kk] = f[l][k]; // the second index is renumbered from 0 to d-3
				}
			}

			w[i][j] = pivot(mat); // mat is f in the memo, w is g
			// matrix that sends the projected (d-2)-faces of the hypercube onto a regular (d-2)-hypercube (projection onto orthogonal of v1,v2)

			tmp1 = mdot(w[i][j], vec1);
			tmp2 = mdot(w[i][j], vec2);

			ww[i][j] = Array(Array(d - 2), Array(d - 2), Array(d - 2), Array(d - 2));
			for (k = 0; k < d - 2; k++) {
				ww[i][j][0][k] = tmp1[k] + tmp2[k];
				ww[i][j][1][k] = tmp1[k] - tmp2[k];
				ww[i][j][2][k] = -tmp1[k] + tmp2[k];
				ww[i][j][3][k] = -tmp1[k] - tmp2[k];
			}
		}

	//
	hash = new Object();
	//hash.length=0; // for testing purposes
	pstack = Array();

	var p = Array(d);
	for (i = 0; i < d; i++)
		p[i] = Math.round(p0[i]);
	insertpoint(p);

	//    debugarea=document.getElementById("debug");

	precomputed = true;
	compute();
	toggleContent();
	redraw();
	
}

function updateSampleSizes() {
	// update sample sizes, in case canvas' been resized
	A = cx / scale;
	B = cy / scale;
}

function compute() {
	if (!precomputed) return;

	// now the main loop
	var p, i;

	scale = 1000 / +document.getElementById("size").value; // only update scale when recomputing
	updateSampleSizes();

	var istk = 0;
	while (istk < pstack.length) {
		//	p=pstack.pop();
		p = pstack[istk]; istk++;
		if ((hash[p] == 1) && check(p)) // if on surface and on the screen (and neighbors not tested yet, in case of multiple runs)
		{
			//	    debugarea.innerHTML+=p+"<br>";
			for (i = 0; i < d; i++) {
				p[i]++;
				insertpoint(p);
				p[i] -= 2;
				insertpoint(p);
				p[i]++;
			}
			hash[p] = 2; // for future reference (in case of multiple runs), neighbors have been inserted
		}
	}
	//    alert(pstack.length);
	if (pstack.length > 1) computed = true; // to avoid bugs, like invalid input...
}



function initcanvas() {
    canvas = document.getElementById('mycanvas');
    if (canvas.getContext) {
        ctx = canvas.getContext('2d');

        // Get the device pixel ratio, falling back to 1 for unsupported browsers
        var pixelRatio = window.devicePixelRatio || 1;

        // Set the display size of the canvas (CSS size)
        var displayWidth = canvas.clientWidth;
        var displayHeight = canvas.clientHeight;

        // Set the canvas width and height based on the pixel ratio
        canvas.width = displayWidth * pixelRatio;
        canvas.height = displayHeight * pixelRatio;

        // Scale the context to account for the pixel ratio
        ctx.scale(pixelRatio, pixelRatio);

        cx = displayWidth;
        cy = displayHeight;

        // Clear the canvas
        ctx.fillRect(0, 0, canvas.width, canvas.height);
    } else {
        alert("Could not find canvas");
    }
}


var col;
var colph1, colph2; // phases

function initcolors() {
	var max = d * (d - 1) / 2;

	

	col = Array(max);
	for (i = 0; i < max; i++) {
		// Cycle through predefined colors or generate random ones if needed
		col[i] = predefinedColors[i % predefinedColors.length];
	}

	ctx.strokeStyle = "black"; // Line color set to black
}


function fillquad(p, c) {
	/*
	debugarea=document.getElementById("debug");
	debugarea.innerHTML+=p[0][0]/A+","+p[0][1]/A+" / "+p[1][0]/A+","+p[1][1]/A+" / "+p[2][0]/A+","+p[2][1]/A+" / "+p[3][0]/A+","+p[3][1]/A+"<br>";
	*/
	ctx.fillStyle = "rgb(" + col[c] + ")";
	ctx.beginPath();
	ctx.moveTo((p[0][0] + A) * scale / 2, (p[0][1] + B) * scale / 2);
	ctx.lineTo((p[1][0] + A) * scale / 2, (p[1][1] + B) * scale / 2);
	ctx.lineTo((p[2][0] + A) * scale / 2, (p[2][1] + B) * scale / 2);
	ctx.lineTo((p[3][0] + A) * scale / 2, (p[3][1] + B) * scale / 2);
	ctx.closePath();
	ctx.fill();
}

function drawline(p1, p2) {
	ctx.beginPath();
	ctx.moveTo((p1[0] + A) * scale / 2, (p1[1] + B) * scale / 2);
	ctx.lineTo((p2[0] + A) * scale / 2, (p2[1] + B) * scale / 2);
	ctx.stroke();
}


function redraw() {
	if (!computed) { ctx.fillRect(0, 0, canvas.width, canvas.height); return; }
	// nothing to do if computation not already done

	updateSampleSizes();

	var s, pr = Array(4);

	var istk, i, j;

	initcolors();

	for (istk = 0; istk < pstack.length; istk++) {
		p = pstack[istk];
		pr[0] = proj(p);
		for (i = 0; i < d - 1; i++) {
			p[i]++;
			if (hash[p]) {
				pr[1] = proj(p);
				for (j = i + 1; j < d; j++) {
					p[j]++;
					if (hash[p]) {
						pr[2] = proj(p);
						p[i]--;
						if (hash[p]) {
							pr[3] = proj(p);
							fillquad(pr, j * (j - 1) / 2 + i); // draw lozenge
						}
						p[i]++;
					}
					p[j]--;
				}
			}
			p[i]--;
		}
	}

	for (istk = 0; istk < pstack.length; istk++) {
		p = pstack[istk];
		pr[0] = proj(p);
		for (i = 0; i < d; i++) {
			p[i]++;
			if (hash[p]) {
				pr[1] = proj(p);
				drawline(pr[0], pr[1]);
			}
			p[i]--;
		}
	}

	//    fillquad([[-0.1,-0.1],[-0.1,0.1],[0.1,0.1],[0.1,-0.1]],0);

}


function changedim(dim) {
	//document.getElementById("coords").style.display= "none";
	if ((dim < 3) || (dim > 20)) { document.getElementById("dim").value = dd; return; }
	var vv = document.getElementById("tab").rows; // list of rows of table of coord
	var i;

	//  rewriting the html will *not* work here so we have to be more subtle
	while (dd < dim) {
		dd++;
		for (i = 0; i < vv.length; i++)
			vv[i].insertCell(dd).innerHTML = "<input />"; // though inside the cell we still rewrite the html
			
	}
	while (dd > dim) {
		for (i = 0; i < vv.length; i++)
			vv[i].deleteCell(dd);
		dd--;
	}
	document.getElementById("dim").value = dd;

	changevalues(); changeoffset();
}


var mouseX, mouseY, mouseDown, resizeX, resizeY;

function handleMouseDown(event) {
	if (event.which == 1) // left click
	{
		mouseDown = true;
		mouseX = event.pageX - canvas.offsetLeft;
		mouseY = event.pageY - canvas.offsetTop;
		return false; // prevent default behavior -- not that we care so much since cursor style moved to onmousemove
	}
}

function handleMouseUp(event) {
	if (event.which == 1) // left declick
	{
		mouseDown = false;
		document.body.style.cursor = "default";
		//	return false; // nasty bug if uncommented with input ``range''
	}
}

function handleMouseMove(event) {
	var newmouseX = event.pageX - canvas.offsetLeft;
	var newmouseY = event.pageY - canvas.offsetTop;
	if (mouseDown) {
		if (resizeX) {
			cx += newmouseX - mouseX;
			if (cx < 50) cx = 50;
			canvas.width = cx;
			if (newmouseX > mouseX) compute(); // get some new points
		}
		if (resizeY) {
			cy += newmouseY - mouseY;
			if (cy < 50) cy = 50;
			canvas.height = cy;
			if (newmouseY > mouseY) compute(); // get some new points
		}
		if (resizeX || resizeY) redraw();
		mouseX = newmouseX;
		mouseY = newmouseY;
	}
	else {
		resizeX = (newmouseX > cx - 20) && (newmouseX < cx) && (newmouseY >= 0) && (newmouseY < cy);
		resizeY = (newmouseY > cy - 20) && (newmouseY < cy) && (newmouseX >= 0) && (newmouseX < cx);
		// change pointer
		if (resizeX && (!resizeY))
			document.body.style.cursor = "e-resize";
		else if (resizeY && (!resizeX))
			document.body.style.cursor = "s-resize";
		else if (resizeX && resizeY)
			document.body.style.cursor = "se-resize";
		else
			document.body.style.cursor = "default";

	}

}

function validatesize(val) {
	if (val < 1) document.getElementById("size").value = val = 1;
	if (val > 100) document.getElementById("size").value = val = 100;
}

function init() {
	// we can't trust the value="" in the html. some browsers don't respect it when reloading...
	initcanvas();
	document.getElementById("dim").value = 3; changedim(3);
	document.getElementById("size").value = 10;
	document.getElementById("values").value = "random"; changevalues();
	document.getElementById("offset").value = "random"; changeoffset();
	
	


	canvas.onmousedown = handleMouseDown;
	document.onmouseup = handleMouseUp;
	document.onmousemove = handleMouseMove;

	readTableContent();

}

function changevalues() {
	d = dd;
	var i, j;
	var opt = document.getElementById("values").value;
	var vv = document.getElementById("tab").rows; 

	if (opt == "random") {
		for (j = 1; j <= d; j++) {
			vv[0].cells[j].firstChild.value = Math.random() - 0.5;
			
		}
	}
	else if (opt == "empty") {
		for (j = 1; j <= d; j++) {
			vv[0].cells[j].firstChild.value = "";
			
		}
	}
}


function changeoffset(bool) {
	var i, j;
	var opt = document.getElementById("offset").value;
	var vv = document.getElementById("tab").rows; 

	if (opt == "random") {
		for (j = 1; j <= dd; j++) {
			vv[1].cells[j].firstChild.value = Math.random() - 0.5;
		}
	}
	else if (opt == "empty") {
		for (j = 1; j <= dd; j++) {
			vv[1].cells[j].firstChild.value = "";
		}
	}
}



function save() {
	//    window.open(canvas.toDataURL());
	x = window.open();
	iframe = x.document.createElement('iframe')
	iframe.src = canvas.toDataURL();
	x.document.body.appendChild(iframe);
}



function readTableContent() {
	d = dd;
	var vv = document.getElementById("tab");
	var v = Array(d);
	var i;

	//this loop reads the values of equation from table and store it to v.
	for(i = 0; i < d; i++){
		v[i] = +vv.rows[0].cells[i + 1].firstChild.value;
	}

	//this loop reads the point from table and store it to point.
	for(i = 0; i < d; i++){
		point[i] = +vv.rows[1].cells[i+1].firstChild.value;
	}

	// Generate random points and calculate vectors from them
    var pp = [];
    for (let j = 0; j < 3; j++) {
        var sum = 0;
        var col = [];
        for (i = 0; i < d - 1; i++) {
           //var a = Math.floor(Math.random() * 21); // Random integer between 0 and 20
		    var a = rn[j][i];
            sum += v[i] * (a - point[i]); // Sum for the equation
            col.push(a); // Store the random value
        }
        col.push((-sum / v[d - 1]) + point[d - 1]); // Calculate the third coordinate
        pp.push(col); // Push the generated point
    }

    // Calculate vectors
    for(i = 0; i < d; i++) {
        vec[0][i] = (pp[2][i] - pp[0][i])*-1; // v1 = pp[2] - pp[0]
        vec[1][i] = pp[1][i] - pp[0][i]; // v2 = pp[1] - pp[0]
    }
	/*console.log(d);
	console.log(vec);
	console.log(point);*/

	precompute();
	
}

function toggleContent() {
	let content = document.getElementById('content');
	content.style.display = (content.style.display === 'none' || content.style.display === '') ? 'block' : 'none';
}

function showInfo() {
	let content = document.getElementById('pdf-container');
	content.style.display = (content.style.display === 'none' || content.style.display === '') ? 'block' : 'none';
}