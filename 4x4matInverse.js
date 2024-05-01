

/*
 * author: Nadeem Elahi
 * nadeem.elahi@gmail.com
 * nad@3deem.com
 * license: gpl v3
 */

//var mat = new Float32Array(16);

function print4x4mat ( mat ) {
	var idx , lim = 4;

	for ( idx = 0 ; idx < lim ; idx ++ ){

		console.log( mat[idx][0]
			+ " " + mat[idx][1]
			+ " " + mat[idx][2]
			+ " " + mat[idx][3]
		);
	}
}

function make4x4array() {
	return [
		[1,0,0,0],
		[0,1,0,0],
		[0,0,1,0], 
		[0,0,0,1]  
	]
}

/*
 *  https://semath.info/src/cofactor-matrix.html
 *
 *  mij
 *  
 *  m11 = 
 *
 *    a22 $ a33 * a44  +  a23 * a34 * a42  +  a24 * a32 * a43
 *  - a24 $ a33 * a42  -  a23 * a32 * a44  -  a22 * a34 * a43
 *
 * i: 1 , big md sm: [2,3,4]
 *
 *     sm    md    bg      sm   md     bg      sm   md     bg
 *    a22 $ a33 * a44  +  a23 * a34 * a42  +  a24 * a32 * a43
 *   
 *     sm    md    bg      sm    md    bg      sm    md    bg
 *  - a24 $ a33 * a42  -  a23 * a32 * a44  -  a22 * a34 * a43
 *
 * j: 1 , big md sm: [2,3,4]
 *
 *     sm    md    bg      md    bg    sm      bg    sm    md
 *    a22 $ a33 * a44  +  a23 * a34 * a42  +  a24 * a32 * a43
 *   
 *     bg    md    sm      md    sm    bg      sm    bg    md
 *  - a24 $ a33 * a42  -  a23 * a32 * a44  -  a22 * a34 * a43
 *
 */

function cofactor4x4( eh , em , cf ) {

	var idexs , jdexs , 
		idexBig , idexMd , idexSm ,
		jdexBig , jdexMd , jdexSm ,
		cnt ;

	var idx , jdx , lim = 5 ;

	for ( idx = 1 ; idx < lim ; idx ++ ) {
		for ( jdx = 1 ; jdx < lim ; jdx ++ ) {

			idexs = [];
			jdexs = [];
			
			for ( cnt = 1 ; cnt < lim ; cnt ++ ){ 
				if ( cnt != idx )
					idexs.push(cnt);

				if ( cnt != jdx ) 
					jdexs.push(cnt);
			}

			idexBig = idexs.pop(); 
			jdexBig = jdexs.pop();
			idexSm = idexs.shift(); 
			jdexSm = jdexs.shift();
			idexMd = idexs.shift(); 
			jdexMd = jdexs.shift();

			em[idx-1][jdx-1] = 
				eh[idexSm-1][jdexSm-1]
				*eh[idexMd-1][jdexMd-1]
				*eh[idexBig-1][jdexBig-1]

				+ eh[idexSm-1][jdexMd-1]
				*eh[idexMd-1][jdexBig-1]
				*eh[idexBig-1][jdexSm-1]

				+ eh[idexSm-1][jdexBig-1]
				*eh[idexMd-1][jdexSm-1]
				*eh[idexBig-1][jdexMd-1] 

				- eh[idexSm-1][jdexBig-1]
				*eh[idexMd-1][jdexMd-1]
				*eh[idexBig-1][jdexSm-1]

				- eh[idexSm-1][jdexMd-1]
				*eh[idexMd-1][jdexSm-1]
				*eh[idexBig-1][jdexBig-1]

				- eh[idexSm-1][jdexSm-1]
				*eh[idexMd-1][jdexBig-1]
				*eh[idexBig-1][jdexMd-1] 
			;


		}
	}

	function calcSign ( idx , jdx ) {
		if ( (idx + jdx ) % 2 ) return -1; // odd
		else return 1; // even
	}

	var sign ;

	for ( idx = 1 ; idx < lim ; idx ++ ) {
		for ( jdx = 1 ; jdx < lim ; jdx ++ ) {

			sign = calcSign ( idx , jdx ) ;

			cf[idx-1][jdx-1] = sign * em[jdx-1][idx-1];
		}
	}
			
	/*
	 * https://semath.info/src/determinant-four-by-four.html
	 *
	 * det = a11*m11 - a21*m21 + a31*m31 - a41*m41
	 *
	 * offset index starting at 0 instead of so -1
	 * det = a00*m00 - a10*m10 + a20*m20 - a30*m30
	 *
	 */
	return eh[0][0]*em[0][0] - eh[1][0]*em[1][0] + eh[2][0]*em[2][0] - eh[3][0]*em[3][0];
}

// a random 4x4 matria 'a' or 'eh'
var eh = [];
eh[0] = [1,1,1,-1];
eh[1] = [1,1,-1,1];
eh[2] = [1,-1,1,1];
eh[3] = [-1,1,1,1];
console.log("---");
print4x4mat( eh );

var em = make4x4array(); // determinants matrix 'm' or 'em'
var cf = make4x4array(); // cofactor matrix 'cf'

var det = cofactor4x4( eh , em , cf );


//console.log("---");
//print4x4mat ( em )

//console.log("---");
//print4x4mat ( cf )


console.log("---");
console.log("det: " + det );

/* 
 * Inverse:
 * https://semath.info/src/inverse-cofactor-ex4.html
 * 
 * says the cofactor matrix is called the adjucate
 */
 
var idx , jdx , lim = 4;
for ( idx = 0 ; idx < lim ; idx ++ ) {
	for ( jdx = 0 ; jdx < lim ; jdx ++ ) {
		cf[idx][jdx] = cf[idx][jdx] / det;	
	}
}

console.log("---");
print4x4mat ( cf )


