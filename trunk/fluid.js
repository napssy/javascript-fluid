///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
 * Global variables declaration
 */

////touch coordinates
var fireLife = 100;
var nIterations = 20;
var minHeat = 100;
var heatForce = 0.0002;

//grid dimensions
var gridX = 100;
var gridY = 100;
var nX;
var nY;

var gSize = gridX * gridY;

//array of screen pixels
var pixels = new Float32Array(gSize);

//velocity fields
var u = new Float32Array(gSize);
var v = new Float32Array(gSize);
//previous frame velocities
var uPrev = new Float32Array(gSize);
var vPrev = new Float32Array(gSize);

//density field
var dens = new Float32Array(gSize);
var densPrev = new Float32Array(gSize);

//vector used to zero out blocks of memory
var zeroVect = new Float32Array(gSize);

/*
 * Used to calculate the index in a 2D array represented as 1D
 * the array must be of gridX width
 */

function gIndex(var x, var y)
{
	return x + y * gridX;
}

var nnX = gridX-2;
var nnY = gridY-2;
var lastLine = nY*gridX;
var preLastLine = nnY*gridX;

///*
// * processInput
// * modify velocity and density fields based on user touch input
// */
//void processInput()
//{
//	if (nTouchLength < 2)
//		return;
//	int i, j, k, j2, k2;
//	int nIndex;
//	int a = MAX(2,exp(0.4*log(gridX)));
//	int mX, mY;
//	int mPrevX = touchX[0] * gridX / screenWidth;
//	int mPrevY = touchY[0] * gridY / screenHeight;
//	fixed c;
//	fixed c2;
//	int dx,dy;
//	for (i = 0; i < nTouchLength; i++)
//	{
//		mX = touchX[i] * gridX / screenWidth;
//		mY = touchY[i] * gridY / screenHeight;
//		dx = mX-mPrevX;
//		dy = mY-mPrevY;
//		c2 = MUL((FROM_INT(6) + FROM_FLOAT(sqrt(dx * dx + dy * dy))),FROM_INT(5));
//		for (j = (-a); j <= a; j++)
//		{
//			for (k = (-a); k <= a; k++)
//			{
//				j2 = MAX(0, MIN(mX + j, gridX-1));
//				k2 = MAX(0, MIN(mY + k, gridY-1));
//				nIndex = gIndex(j2,k2);
//				c = DIV2(FROM_FLOAT(a - sqrt(j * j + k * k)));
//				if (c > 0)
//				{
//					u[nIndex] += MUL(FROM_INT(dx), c);
//					v[nIndex] += MUL(FROM_INT(dy), c);
//					dens[nIndex] += MUL(c,c2);
//				}
//			}
//		}

//		mPrevX = mX;
//		mPrevY = mY;
//	}
//}

/*
 * apply Forces
 */
function applyForces()
{
	var i,j;

	//apply heat
	for (i = 0; i < gSize; i++)
	{
		if (dens[i] > minHeat)
		{
			v[i] -= dens[i] * heatForce;
			u[i] += (Math.random()-0.5) * dens[i] * heatForce;
		}
	}
	
	uPrev = u.slice();
	vPrev = v.slice();

	//apply cooler
	for (i=0; i < gSize; i++)
	{
		densPrev[i] = dens[i] * fireLife;
	}

    //TODO: zero out speed field borders
	//horizontal...top and bottom rows
    for (i = 0; i < gridX; i++)
    {
        uPrev[i] = 0;
        vPrev[i] = 0;
        uPrev[i + lastLine] = 0;
        vPrev[i + lastLine] = 0;
    }

	//vertical...left and right columns
	j=0;
	while( j < gSize )
	{
		uPrev[j] = 0;
		vPrev[j] = 0;
		uPrev[nX + j] = 0;
		vPrev[nX + j] = 0;
		j += gridX;
	}
}

/**
 * set_bnd(int mode, fixed *pArray);
 *
 * Set the boundary conditions for pArray
 *
 * mode = 0...used for density field
 * mode = 1...used for horizontal speed component
 * mode = 2...used for vertical speed component
 */
function set_bnd(var mode, var pArray)
{
	var i;

	//TODO: optimize by memcpy where possible...done
    if (mode == 0)//just copy
    {
    	//copy horizontal
    	for (i = 0; i < gridX; i++)
    	{
    	    pArray[i] = pArray[i + gridX];
    	    pArray[i + lastLine] = pArray[i + preLastLine];
    	}
    	//copy vertical
    	i = 0;
    	while(i < gSize)
		{
			pArray[i] = pArray[i+1];
			pArray[nX + i] = pArray[nnX + i];
			i += gridX;
		}
    }
    else if (mode == 1)	//copy horizontal, mirror vertical
    {
    	//copy horizontal
    	for (i = 0; i < gridX; i++)
    	{
    	    pArray[i] = pArray[i + gridX];
    	    pArray[i + lastLine] = pArray[i + preLastLine];
    	}
    	//mirror vertical
    	i = 0;
    	while(i < gSize)
		{
			pArray[i] = -pArray[i+1];
			pArray[nX + i] = -pArray[nnX + i];
			i += gridX;
		}
    }
    else if (mode == 2) //copy vertical, mirror horizontal
    {
    	//mirror horizontal
    	for (i = 0; i < gridX; i++)
    	{
    	    pArray[i] = -pArray[i + gridX];
    	    pArray[i + lastLine] = -pArray[i + preLastLine];
    	}
    	//copy vertical
    	i = 0;
    	while(i < gSize)
		{
			pArray[i] = pArray[i+1];
			pArray[nX + i] = pArray[nnX + i];
			i += gridX;
		}
    }

    //average out the corners
    pArray[0] = 0.5 * (pArray[1] + pArray[gridX]);
    pArray[nX] = 0.5 * (pArray[nnX], pArray[nX + gridX]);
    pArray[lastLine] = 0.5 * (pArray[preLastLine], pArray[1 + lastLine]);
    pArray[nX + lastLine] = 0.5 * (pArray[nnX + lastLine], pArray[nX + preLastLine]);
}

/* project(...)
 * this function is used to
 * make the velocity field mass conserving
 * and produce the nice vortexes
 */
function project(var pU, var pV, var pProj, var pDiv)
{
	var i,j,k;
	var index;

	for (j = 1; j < nY; j ++)
	{
		for (i = 1; i < nX; i ++)
		{
    		index = gIndex(i,j);
    		pDiv[index] = -pU[index + 1] + pU[index - 1] - pV[index + gridX] + pV[index - gridX];
    		pProj[index] = pDiv[index] / 4;
    	}
    }
    set_bnd(0, pDiv);

    for (k = 1; k < nIterations; k ++)
    {
    	for (j = 1; j < nY; j ++)
		{
			for (i = 1; i < nX; i ++)
			{
				index = gIndex(i,j);
				pProj[index] = 0.25 * (pDiv[index] + pProj[index - 1] + pProj[index + 1]
				                    + pProj[index - gridX] + pProj[index + gridX]);
			}
		}
    	set_bnd(0, pProj);
    }

    for (j = 1; j < nY; j ++)
	{
		for (i = 1; i < nX; i ++)
		{
			index = gIndex(i,j);
			pU[index] -= 0.5 * (pProj[index + 1] - pProj[index - 1]);
			pV[index] -= 0.5 * (pProj[index + gridX] - pProj[index - gridX]);
		}
	}
    set_bnd(1, pU);
    set_bnd(2, pV);
}

/*
 * computeVelocity()
 * this function is used to compute the values of the velocity field
 * based on the values from the previous frame
 */
function computeVelocity()
{
	var fx, fy;
	var s0, s1, t0, t1;
	var m, n;
	var nIndex;
	var idx;

	for (n = 1; n < nY; n++)
	{
		for (m = 1; m < nX; m++)
		{
			nIndex = gIndex(m,n);

			//calculate originating coordinates
			//and clamp to grid
			fx = Math.min(nX, Math.max(0, m - uPrev[nIndex]));
			fy = Math.min(nY, Math.max(0, n - vPrev[nIndex]));

			//calculate the grid coordinates surrounding the originating point
			//clamp to a maximal size to avoid index out of bounds
			idx = gIndex(Math.floor(fx),Math.floor(fy));
			idx = Math.min(gSize-gridX-2,idx);

			//calculate interpolation coefficients
			s0 = FRAC(fx);
			s1 = FIXED_ONE - s0;

			t0 = FRAC(fy);
			t1 = FIXED_ONE - t0;

			//interpolate and get the new velocity values
			u[nIndex] =   MUL(t1, (MUL(s1, uPrev[idx]) + MUL(s0, uPrev[idx + 1])))
						+ MUL(t0, (MUL(s1, uPrev[idx + gridX]) + MUL(s0, uPrev[idx + gridX + 1])));

			v[nIndex] =   MUL(t1, (MUL(s1, vPrev[idx]) + MUL(s0, vPrev[idx + 1])))
						+ MUL(t0, (MUL(s1, vPrev[idx + gridX]) + MUL(s0, vPrev[idx + gridX + 1])));
		}
	}
	set_bnd(1, u);
	set_bnd(2, v);
 }

/*
 * computeDensity()
 * this function is used to compute the new density values
 * based on the old values and the velocity fields
 */
void computeDensity()
{
	fixed fx;
	fixed fy;
	fixed s0;
	fixed s1;
	fixed t1;
	fixed t0;
	int m, n;
	int nIndex;
	int idx;

	for (n = 1; n < nY; n++)
	{
		for (m = 1; m < nX; m++)
		{
			nIndex = gIndex(m,n);

			//calculate originating coordinates
			//and clamp to grid
			fx = MIN(FROM_INT(nX), MAX(0, FROM_INT(m) - u[nIndex]));
			fy = MIN(FROM_INT(nY), MAX(0, FROM_INT(n) - v[nIndex]));

			//calculate the grid coordinates surrounding the originating point
			//clamp to a maximal size to avoid index out of bounds
			idx = gIndex(TO_INT(TRUNC(fx)),TO_INT(TRUNC(fy)));
			idx = MIN(gSize-gridX-2,idx);

			//calculate interpolation coefficients
			s0 = FRAC(fx);
			s1 = FIXED_ONE - s0;

			t0 = FRAC(fy);
			t1 = FIXED_ONE - t0;

			dens[nIndex] = MUL(t1, (MUL(s1, densPrev[idx]) + MUL(s0, densPrev[idx + 1])))
						 + MUL(t0, (MUL(s1, densPrev[idx + gridX]) + MUL(s0, densPrev[idx + gridX + 1])));
		}
	}
	set_bnd(0, dens);
}

/*
 * updatePixels()
 * this function is used to generate the pixel color values
 * from the density values
 */
void updatePixels()
{
	int i;
	int col;
	for (i = 0; i < gSize; i++)
	{
		col = MAX(0,MIN(TO_INT(dens[i]), nPaletteLength-1));
		pixels[i] = palette[col];
	}
}

/**
 * JNI functions
 */

void Java_fixedpointcode_fleya_WorkerThread_fluidInit( JNIEnv* env, jobject javaThis, jint gX, jint gY, jint sX, jint sY, jint inputQueueLength, jintArray pal, jint palLen, jint iterations)
{
	gridX = gX;
	gridY = gY;
	halfX = gridX >> 1;
	halfY = gridY >> 1;
	nX = gridX - 1;
	nY = gridY - 1;
	gSize = gridX * gridY;

	screenWidth = sX;
	screenHeight = sY;

	nPaletteLength = palLen;
	nInputQueueSize = inputQueueLength;

	nIterations = iterations;

	Init();

	(*env)->GetIntArrayRegion(env, pal,0,palLen, palette);
	(*env)->ReleaseIntArrayElements(env, pal, palette, 0 );
}

void Java_fixedpointcode_fleya_WorkerThread_fluidClean( JNIEnv* env, jobject javaThis )
{
	Clean();
}

void Java_fixedpointcode_fleya_WorkerThread_fluidUpdate( JNIEnv* env, jobject javaThis, jintArray bitmap, jintArray tX, jintArray tY, jintArray tAction, jint tlen, jfloatArray accel)
{
//	double end;
//	double start = now_ms();
	nTouchLength = tlen;
	if (nTouchLength > nInputQueueSize)
		nTouchLength = nInputQueueSize;
	if (nTouchLength)
	{
		(*env)->GetIntArrayRegion(env, tX,0,nTouchLength, touchX);
		(*env)->ReleaseIntArrayElements(env, tX, touchX, 0 );

		(*env)->GetIntArrayRegion(env, tY,0,nTouchLength, touchY);
		(*env)->ReleaseIntArrayElements(env, tY, touchY, 0 );

		(*env)->GetIntArrayRegion(env, tAction,0,nTouchLength, touchAction);
		(*env)->ReleaseIntArrayElements(env, tAction, touchAction, 0 );
	}
	(*env)->GetFloatArrayRegion(env, accel,0,3, accRaw);
	(*env)->ReleaseFloatArrayElements(env, accel, accRaw, 0 );
	acc[0] = FROM_FLOAT(accRaw[0]);
	acc[1] = FROM_FLOAT(accRaw[1]);
	acc[2] = FROM_FLOAT(accRaw[2]);
//	perfArgs[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	processInput();
//	perfInput[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	applyForces();
//	perfForce[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	computeVelocity();
//	perfVel[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	project(u, v, uPrev, vPrev);
//	perfProj[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	computeDensity();
//	perfDens[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	updatePixels();
//	perfPix[nFrameCounter] = now_ms()-start;
//
//	start = now_ms();
	//copy our pixel array back to the JVM
	(*env)->SetIntArrayRegion(env, bitmap,0,gridX * gridY, pixels);
//	perfRet[nFrameCounter] = now_ms()-start;
//
//	nFrameCounter++;
//	if (nFrameCounter >= AVG_FRAMES)
//	{
//		bCanPrintPerformance = 1;
//		nFrameCounter = 0;
//	}
//
//	if (bCanPrintPerformance)
//	{
//		LOGD("args:", "%f", getAvg(perfArgs));
//		LOGD("input:", "%f", getAvg(perfInput));
//		LOGD("force:", "%f", getAvg(perfForce));
//		LOGD("vel:", "%f", getAvg(perfVel));
//		LOGD("proj:", "%f", getAvg(perfProj));
//		LOGD("dens:", "%f", getAvg(perfDens));
//		LOGD("pix:", "%f", getAvg(perfPix));
//		LOGD("ret:", "%f", getAvg(perfRet));
//	}
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////


var size = 100;

//single dimension array indexing used for 2 dimensions
function ind(i,j) {return i+size*j}

//boundary conditions
function boundary(size, array, direction)
{
    //edge cells get the value of their coresponding neighbors
    for (var i=1; i < size-1; i++)
    {
        array[ind(0,i)] = (direction == 1) ? -array[ind(1,i)] : array[ind(1,i)];
        array[ind(size-1,i)] = (direction == 1) ? -array[ind(size-2,i)] : array[ind(size-2,i)];
        array[ind(i,0)] = (direction == 2) ? -array[ind(i,1)] : array[ind(i,1)];
        array[ind(i,size-1)] = (direction == 2) ? -array[ind(i,size-2)] : array[ind(i,size-2)];
    }
    //corner cells
    array[ind(0,0)] = 0.5 * (array[ind(1,0)] + array[ind(0,1)]);
    array[ind(0,size-1)] = 0.5 * (array[ind(0,size-2)] + array[ind(1,size-1)]);
    array[ind(size-1,size-1)] = 0.5 * (array[ind(size-1,size-2)] + array[ind(size-2, size-1)]);
    array[ind(size-1,0)] = 0.5 * (array[ind(size-2,0)] + array[ind(size-1,1)]);
}

function diffuse ( size, dst, src, diff, dt )
{
    var i, j, k;
    var a = dt * diff * size * size;
    //number of iterative steps to compute diffusion
    for ( k = 0 ; k < iterative_steps ; k++ )
    {
        for ( i = 1 ; i < size-1; i++ )
        {
            for ( j = 1 ; j < size-1; j++ )
            {
                dst[ind(i,j)] = ((src[ind(i,j)] + a * (src[ind(i-1,j)] + src[ind(i+1,j)]+
                                                      src[ind(i,j-1)] + src[ind(i,j+1)])) / (1+4*a));
            }
        }
        boundary ( size, dst, 0 );
    }
}

function diffuse3 ( size, dr, dg, db, sr, sg, sb, diff, dt )
{
    var i, j, k;
    var a = dt * diff * size * size;
    var c = 1+4*a;
    //number of iterative steps to compute diffusion
    for ( k = 0 ; k < iterative_steps ; k++ )
    {
        for ( i = 1 ; i < size-1; i++ )
        {
            for ( j = 1 ; j < size-1; j++ )
            {
                dr[ind(i,j)] = (sr[ind(i,j)] + a * (sr[ind(i-1,j)] + sr[ind(i+1,j)]+
                                                      sr[ind(i,j-1)] + sr[ind(i,j+1)])) / c;
                dg[ind(i,j)] = (sg[ind(i,j)] + a * (sg[ind(i-1,j)] + sg[ind(i+1,j)]+
                                                      sg[ind(i,j-1)] + sg[ind(i,j+1)])) / c;
                db[ind(i,j)] = (sb[ind(i,j)] + a * (sb[ind(i-1,j)] + sb[ind(i+1,j)]+
                                                      sb[ind(i,j-1)] + sb[ind(i,j+1)])) / c;
            }
        }
        boundary ( size, dr, 0 );
        boundary ( size, dg, 0 );
        boundary ( size, db, 0 );
    }
}

function advect ( size, dst, src, horiz_speed, vert_speed, dt )
{
    var i, j, i0, j0, i1, j1;
    var x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*size;
    for ( i = 1 ; i< size -1; i++ ) 
    {
        for ( j = 1 ; j < size -1; j++ )
        {
            //compute point of origin for the current (i,j) cell
            x = i - dt0 * horiz_speed[ind(i,j)];
            y = j - dt0 * vert_speed[ind(i,j)];
            if (x < 0.5) x = 0.5;
            if (x > size - 1.5) x = size - 1.5;
            i0 = x >> 0;
            i1 = i0 + 1;
            if (y < 0.5) y = 0.5;
            if (y > size - 1.5) y = size - 1.5;
            j0 = y >> 0;
            j1 = j0 + 1;
            t1 = y-j0;
            t0 = 1-t1;
            //interpolate
            dst[ind(i,j)] = (i1 - x) * ( t0*src[ind(i0,j0)] + t1*src[ind(i0,j1)] )+
                            (x - i0) * ( t0*src[ind(i1,j0)] + t1*src[ind(i1,j1)] );
        }
    }
    boundary ( size, dst, 0);
}

function advect3 ( size, dr, dg, db, srcr, srcg, srcb, horiz_speed, vert_speed, dt )
{
    var i, j, i0, j0, i1, j1;
    var x, y, s0, t0, s1, t1, dt0;
    dt0 = dt*size;
    for ( i = 1 ; i< size -1; i++ ) 
    {
        for ( j = 1 ; j < size -1; j++ )
        {
            //compute point of origin for the current (i,j) cell
            x = i - dt0 * horiz_speed[ind(i,j)];
            y = j - dt0 * vert_speed[ind(i,j)];
            x = Math.max(0.5, Math.min(size-1.5,x));
            i0 = x >> 0;
            i1 = i0 + 1;

            y = Math.max(0.5, Math.min(size-1.5,y));
            j0 = y >> 0;
            j1 = j0 + 1;

            t1 = y-j0;
            t0 = 1-t1;
            s1 = i1 - x;
            s0 = x - i0;
            //interpolate
            dr[ind(i,j)] = s1 * ( t0*srcr[ind(i0,j0)] + t1*srcr[ind(i0,j1)] )+
                            s0 * ( t0*srcr[ind(i1,j0)] + t1*srcr[ind(i1,j1)] );
            dg[ind(i,j)] = s1 * ( t0*srcg[ind(i0,j0)] + t1*srcg[ind(i0,j1)] )+
                            s0 * ( t0*srcg[ind(i1,j0)] + t1*srcg[ind(i1,j1)] );
            db[ind(i,j)] = s1 * ( t0*srcb[ind(i0,j0)] + t1*srcb[ind(i0,j1)] )+
                            s0 * ( t0*srcb[ind(i1,j0)] + t1*srcb[ind(i1,j1)] );
        }
    }
    boundary ( size, dr, 0);
    boundary ( size, dg, 0);
    boundary ( size, db, 0);
}

function compute_densityes ( size, dr, dg, db, horiz_speed, vert_speed, diff, dt )
{
    var difr = [];
    var difg = [];
    var difb = [];
    diffuse3 ( size, difr, difg, difb, dr, dg, db, diff, dt );
    advect3 ( size, dr, dg, db, difr, difg, difb, horiz_speed, vert_speed, dt );
}

function compute_speeds ( size, hspeed, vspeed, visc, dt )
{
    var diffusedh = [];
    var diffusedv = [];
    diffuse ( size, diffusedh, hspeed, visc, dt );
    diffuse ( size, diffusedv, vspeed, visc, dt );
    project( size, diffusedh, diffusedv);
    advect ( size, hspeed, diffusedh, diffusedh, diffusedv, dt );
    advect ( size, vspeed, diffusedv, diffusedh, diffusedv, dt );
    project( size, hspeed, vspeed);
    project( size, hspeed, vspeed);
    project( size, hspeed, vspeed);
    project( size, hspeed, vspeed);
}

function project ( size, h_speed, v_speed)
{
    var i, j, k, h;
    h = 1 / size;
    var div = [];
    var p = [];
    
    for ( i = 1 ; i < size - 1; i++ )
    {
        for ( j = 1 ; j < size -1; j++ )
        {
            div[ind(i,j)] = -0.5*h*(h_speed[ind(i+1,j)]-h_speed[ind(i-1,j)]+
                                    v_speed[ind(i,j+1)]-v_speed[ind(i,j-1)]);
            p[ind(i,j)] = div[ind(i,j)]/4;
        }
    }
    boundary(size,div, 0);
    boundary(size,p, 0);
    for ( k = 0; k < iterative_steps - 1; k++ )
    {
        for ( i = 1; i < size -1; i++ )
        {
            for ( j = 1; j < size -1; j++ )
            {
                p[ind(i,j)] = (div[ind(i,j)]+p[ind(i-1,j)]+p[ind(i+1,j)]+
                                             p[ind(i,j-1)]+p[ind(i,j+1)])/4;
            }
        }
        boundary ( size, p, 0 );
    }
    for ( i = 1; i < size - 1; i++ )
    {
        for ( j = 1; j < size - 1; j++ )
        {
            h_speed[ind(i,j)] -= 0.5*(p[ind(i+1,j)]-p[ind(i-1,j)])/h;
            v_speed[ind(i,j)] -= 0.5*(p[ind(i,j+1)]-p[ind(i,j-1)])/h;
        }
    }
    boundary(size,h_speed, 1);
    boundary(size,v_speed, 2);
}