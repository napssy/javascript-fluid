var size = 100;

//single dimension array indexing used for 2 dimensions
function ind(i,j) {return i+size*j}

//boundary conditions
function boundary(size, array)
{
    //edge cells get the value of their coresponding neighbors
    for (var i=1; i < size-1; i++)
    {
        array[ind(i,0)] = array[ind(i,1)];
        array[ind(i,size-1)] = array[ind(i,size-2)];
        array[ind(0,i)] = array[ind(1,i)];
        array[ind(size-1,i)] = array[ind(size-2,i)];
    }
    //corner cells
    array[ind(0,0)] = array[ind(1,1)];
    array[ind(0,size-1)] = array[ind(1,size-2)];
    array[ind(size-1,size-1)] = array[ind(size-2,size-2)];
    array[ind(size-1,0)] = array[ind(size-2,1)];
}

function diffuse ( size, dst, src, diff, dt )
{
    var i, j, k;
    var a = dt * diff * size * size;
    //number of iterative steps to compute diffusion
    for ( k = 0 ; k < 1 ; k++ )
    {
        for ( i = 1 ; i < size-1; i++ )
        {
            for ( j = 1 ; j < size-1; j++ )
            {
                dst[ind(i,j)] = ((src[ind(i,j)] + a * (src[ind(i-1,j)] + src[ind(i+1,j)]+
                                                      src[ind(i,j-1)] + src[ind(i,j+1)])) / (1+4*a))>>0;
            }
        }
        boundary ( size, dst );
    }
}