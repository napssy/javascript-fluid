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