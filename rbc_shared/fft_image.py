from pylab import *
import numpy as np
from scipy import fftpack

def fft_image( img ):
    """
    Perform 2D fft on a cell frame and return the shifted (centered)
    version.
    """
    transform = fftpack.fft2( img )
    # shifts spatial frequencies to center of image
    return fftpack.fftshift( transform )

def ifft_image( X ):
    """
    Return to the spatial domain.
    """
    itransform = fftpack.ifftshift( X )
    return fftpack.ifft2( itransform )

def ideal_low_pass( X, r=0.2 ):
    """
    Determine all frequency components that a distance < r*C, where C = min(
    from center (k_0, l_0) = (M/2,N/2) and X.shape =
    (M,N). Return a binary array with only Freq < threshold
    unmasked. (Default returns the DC component.)
    """
    nx,ny = X.shape
    max_idx = min( nx, ny )
    # row,col center
    k0,l0 = nx/2, ny/2
    H = zeros_like( X )

    # This code yields all freq. in L1 dist of the center
    # if r == 0:
    #     H[k0,l0] = 1
    #     return H
    # # max distance from center
    # r = int( r*min( k0, l0 ) )
    #  # "extract" region around center
    # for i in xrange(k0-r,k0+r): 
    #     for j in xrange(l0-r,l0+r):
    #         H[i,j] = 1.0
    # return H
    
    d = lambda x,y : x*x + y*y
    # probably a better way, but this'll work for now. We need a set of
    # indices describing a circle around the center of the image
    # (k0,l0)
    # percentage of the minimum distance from center to edge
    r = int( r*min( k0, l0 ) )

    print "r", r
    
    kmax = k0
    lmax = l0
    H = np.zeros_like( X )
    keep_going = True
    for i in xrange( k0, nx ):
        # toggle to True below if we find another column
        keep_going = False
        for j in xrange( l0, ny ):
            # set l0 +/- j  index
            if d( k0-i, l0-j ) <= r:
                # still more columns to try
                keep_going = True
                H[i,j] = X[i,j]
                # l0 - (j-l0) = 2l0 -j, etc.
                H[i, l0+l0-j] = X[i, l0+l0-j]
                H[k0+k0-i, j] = X[k0+k0-i, j] 
                H[k0+k0-i, l0+l0-j] = X[k0+k0-i, l0+l0-j] 
            # otherwise, move to next row
            else:
                break
        # the last row does not contain frequencies within threshold
        if not keep_going:
            break
    return H
            
image_name = '/data/jberwald/rbc/New/frames/new_140125/new_140125-concatenated-ASCII_1000'
image = np.loadtxt( image_name )
X = fft_image( image )
# just keep components within a certain distance of the DC component
Y = ideal_low_pass( X, r=1 )
Yinv = ifft_image( Y )
PX = np.abs( X )**2
PY = np.abs( Y )**2

# #inverse transform. 
# Yinv = fftpack.ifft2( Y ) 
# PYinv = np.abs( Yinv )**2

fig = figure()
# original image
fig.add_subplot( 231 )
imshow( image )
# FT of image
fig.add_subplot( 234 )
imshow( log1p( PX ) )
# filter
fig.add_subplot( 232 )
imshow( log1p( np.abs( Y )**2 ) )
# filtered image
fig.add_subplot( 235 )
imshow( log1p( np.abs( Yinv )**2 ) )

fig.show()

