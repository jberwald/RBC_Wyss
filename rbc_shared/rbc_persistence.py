import numpy as np
from matplotlib import pyplot as plt

# To compute differences between 4-connected neighbors in 3 x 3 patch
Dnorm = np.array( [[ 2, -1, 0,-1, 0,0,0,0,0 ],
                   [-1, 3, -1, 0, -1, 0,0,0,0],
                   [0, -1, 2, 0, 0, -1, 0,0,0],
                   [-1, 0,0,3,-1,0,-1,0,0],
                   [0,-1,0,-1,4,-1,0,-1,0],
                   [0,0,-1,0,-1,3,0,0,-1],
                   [0,0,0,-1,0,0,2,-1,0],
                   [0,0,0,0,-1,0,-1,3,-1],
                   [0,0,0,0,0,-1,0,-1,2]], dtype=np.float )                 
                                              
class Patch( object ):
    """
    A single N x N patch from an image. 
    """
    def __init__( self, pixels ):
        """
        A collection of methods to deal with the individual patches in an image.

        Input:
        -----

        pixels : Square of pixels of size NxN from data image in
        RBCFrame. The pixels are raveled in order to place them into
        an N^2- dimensional space. Note: order='F' implies
        column-major Fortran ordering. This is the raveling used in
        [1].

        self.patch is stored as a column vector and in matrix form for
        proper multiplication below.
        """
        #self.patch = np.asmatrix( np.log( pixels.ravel( order='F' ) ) ).T
        self.patch = np.asmatrix( pixels.ravel( order='F' ) ).T
        self.dnorm = 0

    ## def dist( self, v ):
    ##     """
    ##     Angular distance between two vectors/patches
    ##     """
    ##     np.arccos( s

    def __repr__( self ):
        return repr( self.patch )

    def _subtract_mean( self ):
        return self.patch -self.patch.mean()

    def _D_norm( self, D ):
        """
        This is the 'contrast-normailization' procedure from [1]. In
        order to play safe with memory, D is stored only once and
        passed in for each patch.
        """
        v = self._subtract_mean()
        self.dnorm = np.sqrt( v.T * D * v )

    def ellipsoid_transform( self, D=Dnorm ):
        """
        Subtract the mean and normalize by the 'D-norm'. (See [1]).
        """
        xmean = self._subtract_mean()
        y= self._D_norm( xmean, D )
        self.patch = xnum / xdenom
        

class RBCFrame( Patch ):
    """
    Holds a single image of a red blood cell.

    Input:
    -----

    Either a filename pointing to an image array (matrix, not jpg) on
    disk, or the actual array itself.

    RBCFrame( fname='cell_frame' )

    or

    RBCFrame( data=<np.ndarray> )

    We take the absoute value since eventually we will take logs.

    Methods:
    -------

    mask_image : mask the background so we only 'see' the circular cell 
    
    extract_patches : expects a masked array. Extracts patches (3x3 by default). See [1]

    transform2ellipsoid : Using Dnorm, transform patches to 'contrast-normalized' ellipsoid

    sort_patches : Sort patches by 

    [1] : Lee, et al., 'The Nonlinear Statistics of High-Contrast Patches in Natural Images', Int. Jour. Comp. Vis., 2003
    [2] : Ghrist, 'Barcodes: The Persistence Topology of Data', Bull. AMS, 2008 
    """
    def __init__( self, fname=None, data=None, delimiter=',' ):

        # truncate to 65x70, vestige larger array with multiple modes. FIXME later.
        if fname is not None:
            try:
                self.image = np.absolute( np.loadtxt( fname )[:,:70] )
            except ValueError:
                self.image = np.absolute( np.loadtxt( fname, delimiter=delimiter )[:,:70] )
            except:
                print "Problem loading file. Check the delimiter." 
        # reshape data, just in case we have the original 65x280 version
        elif data is not None:
            self.image = np.absolute( data[:,:70] )
        else:
            raise ValueError, "Must pass in either file name or data array!"

        # Other stuff
        self.debug = True
        self.masked = False
        self.patches = []

    def mask_image( self ):
        self.masked_image = mask_data( self.image, level='bg' )
        self.masked = True

    def extract_patches( self, size=3 ):
        """
        From the masked array, which matches the circular portion,
        construct a collection of nx x ny disjoint patches. See [1] and [2].
        """
        # create stride arrays
        nx, ny = self.image.shape
        sx = np.arange( 0, nx, size )
        sy = np.arange( 0, ny, size )

        # Consider size x size blocks in image. If all entries
        # non-zero, then we assume it is interior to the boundary, so
        # we append it the list of patches.
        z = 0
        for i in range( len(sx)-1 ):
            for j in range( len(sy)-1 ):
                p = self.image[sx[i]:sx[i+1], sy[i]:sy[i+1]]
                z += 1
                if np.all( p != 0 ):
                    self.patches.append( Patch( np.log(p) ) )
        if self.debug:
            print "Found", len(self.patches), "interior patches out of a total of", z, "possible."

    def transform( self ):
        """
        Use D-norm to transfom data to a 7-dimensional ellipsoid
        embedded in R^9.
        """
        self.ellipsoid = [ p.ellipsoid_transform() for p in self.patches ]

    def show_mat( self ):
        if self.masked:
            plt.imshow( self.masked_image )
            plt.colorbar()
        else:
            plt.imshow( self.image )
            plt.colorbar()


def mask_data( data, level='bg' ):
    
    # mask background values
    if level == 'mean':
        md = np.ma.masked_less_equal( data, data.mean() )
    elif level == 'bg':
        md = np.ma.masked_less_equal( data, 0 )
    else:
        print "eek!"
    return md


if __name__ == "__main__":

    patch = np.arange( 9 )
    patch.resize( (3,3) )

    P = Patch ( patch )

    fname = '/Users/jesseberwald/data/rbc/C_4/frames/matC1_100'
    R = RBCFrame( fname=fname )
    R.extract_patches()
    
    
