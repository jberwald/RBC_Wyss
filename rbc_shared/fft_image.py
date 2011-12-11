from pylab import *
import numpy 
from scipy import fftpack
import os, errno
import cPickle as pkl
import chomp_betti
import time, shutil

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
    # row,col center NOTE: l0 == "el zero"
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
    
    # Using the square of the distance
    d = lambda x,y : x*x + y*y
    # probably a better way, but this'll work for now. We need a set of
    # indices describing a circle around the center of the image
    # (k0,l0)
    # percentage of the minimum distance from center to edge
    r = int( r*min( k0, l0 ) )
    # Only keep the DC component
    if r == 0:
        H[k0,l0] = X[k0,l0]
        return H
    
    #  print "r", r
    
    kmax = k0
    lmax = l0
    H = numpy.zeros_like( X )
    keep_going = True
    for i in xrange( k0, nx ):
        # toggle to True below if we find another column
        keep_going = False
        for j in xrange( l0, ny ):
            # set l0 +/- j  index
            if d( k0-i, l0-j ) <= r*r:
                # still more columns to try
                keep_going = True
                H[i,j] = X[i,j]
                # To get three other quadrants: l0 - (j-l0) = 2l0 -j,
                # etc.
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

def plot_3d( img ):


    # import matplotlib as mpl
    # from mpl_toolkits.mplot3d import Axes3D
    # import numpy as numpy
    # import matplotlib.pyplot as plt
    
    # mpl.rcParams['legend.fontsize'] = 10
    
    # fig = plt.figure()
    # ax = fig.gca(projection='3d')
    # # theta = numpy.linspace(-4 * numpy.pi, 4 * numpy.pi, 100)
    # z = numpy.linspace(-2, 2, 100)
    # r = z**2 + 1
    # x = r * numpy.sin(theta)
    # y = r * numpy.cos(theta)
    # ax.plot(x, y, z, label='parametric curve')
    # ax.legend()

    # plt.show()

    from mpl_toolkits.mplot3d import Axes3D
    from matplotlib import cm
    from matplotlib.ticker import LinearLocator, FormatStrFormatter
    import matplotlib.pyplot as plt
    import numpy as numpy
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    kmax, lmax = img.shape
    X = numpy.arange( 0, kmax+1 )
    Y = numpy.arange( 0, lmax+1 )

    print "X"
    print X
    print "Y"
    print Y
    
    X, Y = numpy.meshgrid(X, Y)
    # R = numpy.sqrt(X**2 + Y**2)
    Z = img #numpy.sin(R)
    print Z.shape
    surf = ax.plot_surface(X, Y, Z, rstride=1, cstride=1, cmap=cm.jet,
                           linewidth=0, antialiased=False)
    ax.set_zlim(-1.01, 1.01)
    
    ax.zaxis.set_major_locator(LinearLocator(10))
    ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
    
    fig.colorbar(surf, shrink=0.5, aspect=5)

    plt.show()

def make_dir( fdir ):
    """
    Try to make directory. 
    """
    try:
        os.makedirs( fdir )
    except OSError, e:
        if e.errno != errno.EEXIST:
            raise

def compute_fft( img_name ):
    image = numpy.loadtxt( image_name )
    X = fft_image( image )        

def save_fft( data, fname ):
    """
    Pickle data to file fname.
    """
    with open( fname, 'w' ) as fh:
        # pickl with most efficient protocol available
        pkl.dump( data, fh, protocol=-1 )

def fft2cub( fname ):
    """
    Convert thresholded low pass cell image to cubicle file.
    """
    prefix = fname.rpartition( '/' )[0] + '/'
    stripname = fname.rpartition( '/' )[-1].rstrip( '.pkl' )
    savename = prefix + stripname + '.cub'
    with open( fname ) as fh:
        fft_data = pkl.load( fh )
    m = fft_data['mean']
    nz = fft_data['ifft_nonzero']
    # threshold mask doesn't work with DC component
    if 'r00.' in fname:
        w = numpy.where( nz != 0 )
    else:
        nz_mask = numpy.ma.masked_less( nz, m ).mask
        # False == 0
        w = numpy.where( nz_mask==0 )
    # zip locations of thresholded values to get coords (2D)
    z = zip( w[0], w[1] )
    coords = [ str(x)+'\n' for x in z ]
    # now save coords to disk
    with open( savename, 'w' ) as fh:
        fh.writelines( coords )
        

def run_fft_filter( files, bnd_file, modes, make_cubs=True, run_chomp=True ):
    """
    For each frame in fdir, compute fft( frame ), then take the top
    <modes> percent of the frequencies. Threshold these filtered
    images at the mean (of pixels inside the cell boundary). Save the
    thresholded image to fft_frames/ directory.

    fdir : directory to frames

    modes : percentage of modes (low->high) to use for low-pass filter
    """
    # grab all frames, skip the directories
    if os.path.isdir( files ):
        fdir = files + '/'
        dlist = os.listdir( fdir )
        if os.uname()[0] == 'Linux':
            frames = []
            for f in dlist:
                if f.endswith('npy') and not os.path.isdir( fdir+f ):
                    frames.append( f )
        else:
            frames = [ f for f in dlist if not os.path.isdir( fdir+f ) ]
    else:
        frames = [ files ]
        fdir = files.rpartition( '/' )[0].rpartition( '/' )[0] +'/'
    if not fdir.endswith( '/' ): fdir += '/'
    savedir = fdir + 'fft_frames/'
    # create the directory if necessary
    make_dir( savedir )
    
    for frame in frames:
        if frame.endswith( 'npy' ):
            savename = frame.rstrip( '.npy' )
        elif frame.endswith( 'txt' ):
            savename = frame.rstrip( '.txt' )
        else:
            savename = frame
        fft_data = {}
        # store the modes
        fft_data[ 'modes' ] = modes
        try:
            image = numpy.loadtxt( fdir+frame )
        except ValueError:
            image = numpy.load( fdir+frame )
        except IOError:
            raise
        bnd = numpy.loadtxt( bnd_file )
        X = fft_image( image )        
        Y = ideal_low_pass( X, r=modes )
        Yinv = ifft_image( Y ) 
        Ypow = log1p( numpy.abs( Yinv )**2 ) 
        # crop the log of the low pass to get rid of the noise outside
        # of the cell (zero out external points using boundary array)
        crop_zeros = bnd * Ypow
        crop_zeros.resize( Ypow.shape )
        # append the cell with zeros cropped [-4]
        #Y.append( crop_zeros ) 
        m = crop_zeros[numpy.where( crop_zeros!=0 )].mean()
        # everything above the mean, [-3]
        #Ymean = numpy.ma.masked_less_equal( crop_zeros, m ) 
        # just the boolean mask, [-2]
        #Ymask = numpy.ma.masked_less_equal( crop_zeros, m ).mask 
        # now populate fft_data with the essentials and save to disk
        #fft_data['fft'] = X
        #fft_data['lowpass'] = Y
        #fft_data['ifft'] = Yinv
        #fft_data['fft_mag'] = Ypow
        fft_data['ifft_nonzero'] = crop_zeros
        fft_data['mean'] = m
        # masking is easy, so no need to save.
        #fft_data['ifft_mask_mean'] = Ymean
        #fft_data['ifft_mask'] = Ymask

        # create filename with chomp-readable name (remove the 'dot'
        # from mode)
        mode_part = str( modes ).partition( '.' )
        mode_str = mode_part[0] + mode_part[2]
        dataname = savedir+savename+'_r'+mode_str+'.pkl'
        save_fft( fft_data,  dataname )

        # additional work
        if make_cubs:
            fft2cub( dataname )
        if run_chomp:
            if not make_cubs:
                pass
            else:
                cub_name = dataname.rstrip( '.pkl') + '.cub'
                save_prefix = cub_name.rstrip( '.cub' )
                #save_name = save_prefix + 'cbetti' 
                chomp_betti.run_chomp( cub_name, save_prefix+'.cbetti' )
                chomp_betti.extract_betti( save_prefix )

def concatenate_fft_modes( cell_dir, modes ):
    """
    Return a list of all <cell_name>_rNNN.betti files in cell_dir.
    """
    dlist = os.listdir( cell_dir )
    m = str( modes )
    return [ x for x in dlist if x.endswith( 'r'+m+'.betti' ) ]

def fft_analyze( cell_dir, modes, dim=0, max_dim=2 ):
    """
    Compute statistics on Betti numbers for given cell at given modes.

    cell_dir : path to directory containing .betti files

    modes : Percentage of modes used in low pass filter of original
    FFT of cell image.

    dim : Dimension of homology to compute statistics for.

    max_dim : Set to one greater that the ambient space (i.e., there
    are no 2-spheres in a flat 2D image)

    Returns mean, std, and median for 
    """
    # remove the 'dot'
    if type( modes ) != str:
        mode_part = str( modes ).partition( '.' )
        modes = mode_part[0] + mode_part[2]
    print "Finding all .betti files with mode", str( modes )
    cell_dir += 'r'+modes +'_betti/'
    print "dir", cell_dir
    dlist = os.listdir( cell_dir ) 
    # concatenate_fft_modes( cell_dir, modes )
    loadtxt = numpy.loadtxt
    betti_list = [ loadtxt( cell_dir+f )[:,1] for f in dlist ]
    betti_arr = numpy.array( betti_list )
    return betti_arr[:,dim].mean(), betti_arr[:,dim].var() #, numpy.median( betti_arr[:,dim] )

def run_chomp( cub_name ):
    """
    """
    save_prefix = cub_name.rstrip( '.cub' )
    chomp_betti.run_chomp( cub_name, save_prefix+'.cbetti' )
    chomp_betti.extract_betti( save_prefix )


def fix_zero_mode( chomp=True ):
    """
    A typo messed up the thresholding on DC component (r=0.0).
    """
    new_fdir = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/'
    old_fdir = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/'
    new_cells = [ 'new_110125/', 'new_130125/', 'new_140125/', 'new_40125/', 'new_50125/' ]
    old_cells = [ 'old_100125/', 'old_120125/', 'old_50125/', 'old_90125/' ]

    for cell in new_cells + old_cells:
        if cell.startswith( 'new' ):
            cell_dir = new_fdir
        else:
            cell_dir = old_fdir
        cell_path = cell_dir + cell + 'fft_frames/'
        #fdir = cell_dir + cell + 'fft_frames/'
        dlist = os.listdir( cell_path )
        if chomp:
            suffix = 'cub'
        else:
            suffix = 'pkl'
        dlist = [ f for f in dlist if f.endswith( 'r00.'+suffix ) ]
        if chomp:
            # loop opt.
            C = run_chomp
            for f in dlist:
                C( cell_path + f )
        else:
            for f in dlist:
                fft2cub( cell_path + f )

                
def analyze_all():
    """
    """
    # cells and prefix-directories
    new_fdir = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/'
    old_fdir = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/'
    # new_cells = [ 'new_110125/' ]
    # old_cells = [ 'old_100125/' ]
    new_cells = [ 'new_110125/', 'new_130125/', 'new_140125/', 'new_40125/', 'new_50125/' ]
    old_cells = [ 'old_100125/', 'old_120125/', 'old_50125/', 'old_90125/' ]
    
    all_modes = numpy.linspace( 0, 1, 21 )
    all_dims = [ 0, 1 ]
    all_cells = dict.fromkeys( new_cells + old_cells )
    #dims = dict.fromkeys( all_dims )
    # for d in dims:
    #     dims[ d ] = dict.fromkeys( all_modes )
    t0 = time.time()
    
    for cell in all_cells:
        if cell.startswith( 'new' ):
            cell_dir = new_fdir
        else:
            cell_dir = old_fdir
        cell_path = cell_dir + cell + 'fft_frames/'
        print "cell path", cell_path
        # make a blank copy of the dim
        dim_modes = {}
        
        for dim in all_dims:
            betti_stats = {}
            #dims[ dim ] = dict.fromkeys( all_modes )
            # for dimension dim, run through all modes.
            for r in all_modes:
                betti_stats[r] = fft_analyze( cell_path, r, dim )
            dim_modes[ int(dim) ] = betti_stats
        print "dim modes", dim_modes
        savename = '/data/jberwald/wyss/data/Cells_Jesse/fft_betti_means_'+cell[:-1]+'.pkl'
        with open( savename, 'w' ) as fh:
            pkl.dump( dim_modes, fh )
        print "saved means for", cell
        print ""
    print "Done analyzing betti number for fft!"
    print "Total time:", time.time() - t0


def mode2str( mode ):
    mode_part = str( mode ).partition( '.' )
    mode_str = mode_part[0] + mode_part[2]
    return mode_str

def partition2subdirectories():
    """
    Split up a directory with many files into a number of
    subdirectories. Hopefully this speeds up file IO operations.
    """
    new_fdir = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/'
    old_fdir = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/'
    new_cells = [ 'new_110125/', 'new_130125/', 'new_140125/', 'new_40125/', 'new_50125/' ]
    old_cells = [ 'old_100125/', 'old_120125/', 'old_50125/', 'old_90125/' ]
    
    # new_cells = [ 'new_110125/' ]
    # old_cells = []

    modes = numpy.linspace( 0, 1, 21 )
    # loop opt.
    move = shutil.move
    for cell in new_cells + old_cells:
        print "cell", cell
        if 'old' in cell:
            pre = old_fdir
        else:
            pre = new_fdir
        path = pre + cell + 'fft_frames/'
        dlist = os.listdir( path )
        for mode in modes:
            print "mode", mode
            mode_str = mode2str( mode )
            r_mode = 'r'+str( mode_str )
            subdir = path + r_mode +'_cub/'
            os.mkdir( subdir )
            # now copy files to mode folder
            r_mode_suf = r_mode + '.cub'
            mlist = [ f for f in dlist if f.endswith( r_mode_suf ) ]
            print "moving files..."
            print ""
            for f in mlist:
                move( path + f, subdir + f )

def concat_fft( fdir, dim ):
    """
    Read fft_betti_mean_* dicts from disk, either dim 0 or dim 1, and
    return list of ( mode x mean) arrays.
    """
    new_cells = [ 'new_110125', 'new_130125', 'new_140125', 'new_40125', 'new_50125' ]
    old_cells = [ 'old_100125', 'old_120125', 'old_50125', 'old_90125' ]
    old = {}
    new = {}
    for cell in new_cells + old_cells:
        fname = 'fft_betti_means_' + cell + '.pkl'
        with open( fdir + fname ) as fh:
            a = pkl.load( fh )
        data = a[ dim ]
        nx = data.keys()
        nx.sort()
        # form array of prober size
        mean_arr = numpy.zeros( (2, len(nx)) )
        for i, v in enumerate( nx ):
            mean_arr[0,i] = v
            mean_arr[1,i] = data[v][0]
        if 'new' in cell:
            new[cell] = mean_arr
        else:
            old[cell] = mean_arr
    return new, old
        
def plot_fft_means( data, **args ): #dim=, fig=None, leg=False ):
    """
    data must be a dict keyed by cell name, paired with arrays (mode x
    mean homology generators).
    """
    from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, \
        AnnotationBbox
    from matplotlib._png import read_png
    
    fargs = { 'dim' : 0,
              'fig' : None,
              'leg' : False,
              'zoom' : -1
              }
    fargs.update( args )
    if fargs['fig'] is None:
        fig = figure( frameon=False)
    else:
        fig = fargs['fig']
        #ax = fig.gca()
    ax = fig.add_subplot(111)#, frameon=False, xticks=[], yticks=[])

    for cell in data:
        datum = data[ cell ]
        if 'new' in cell:
            color = 'r'
        else:
            color = 'b'
        
        ax.plot( datum[0][:fargs['zoom']], datum[1][:fargs['zoom']], '-o', c=color, lw=2, ms=3, label=cell )
    ax.hlines( 1, 0, datum[0][fargs['zoom']], linestyle='dashed', linewidth=2, alpha=0.7 )
    # ax.set_xlabel( r'Percent of Fourier modes ($\times 100$)' )
    # ax.set_ylabel( r'Mean number of generators for $H_{'+str(fargs['dim'])+'}$' )
    if fargs['leg']:
        ax.legend( loc='upper left' )

    # offset image. read from previously-saved zoom image.
    #fn = get_sample_data( fargs['zoom_img'] + '.png', asfileobj=False )
    # zoom_img = read_png( fargs['zoom_img'] + '.png' )
    # zoombox = OffsetImage( zoom_img, zoom=0.2 )
    # xy = ( 0, 0 )
    # ab = AnnotationBbox(zoombox, xy,
    #             xybox=( 200., 200.),
    #             xycoords='data',
    #             boxcoords="offset points",
    #             pad=0.5
    #             # arrowprops=dict(arrowstyle="->",
    #             #                 connectionstyle="angle,angleA=0,angleB=90,rad=3")
    #                     )


    # ax.add_artist(ab)
    # draw()

    show()

    return fig

    
    

if __name__ == "__main__":

    analyze_all()
    #    concat_means()

    # sage_session = False
    # try:
    #     from sage.all import *
    #     # from sage.plot.plot3d.list_plot3d import list_plot3d
    #     # from sage.plot.matrix_plot import matrix_plot
    #     sage_session = True
    # except ImportError:
    #     print "Not running sage..."


    # # image_name = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/new_140125/new_140125-concatenated-ASCII_1000.npy'
    # # bnd_name =  '/data/jberwald/wyss/data/Cells_Jesse/New/boundary_Nov_new140125'
    # image_name = '/data/jberwald/rbc/New/frames/new_140125/new_140125-concatenated-ASCII_1000'
    # bnd_name =  '/data/jberwald/rbc/New/boundary_Nov_new140125'
    # # image_name = '/data/jberwald/rbc/Old/frames/old_120125/old_120125-concatenated-ASCII_1000'
    # # bnd_name =  '/data/jberwald/rbc/Old/boundary_Nov_old120125'
    # try:
    #     image = numpy.loadtxt( image_name )
    # except ValueError:
    #     image = numpy.load( image_name )
    # bnd = numpy.loadtxt( bnd_name )
    # X = fft_image( image )
    # # just keep components within a certain distance of the DC component
    # lowpass = {}
    # step = 0.05
    # # for i in numpy.arange(0,2.5,0.5):#+[20]:
    # #     r = i*step
    # #     Y = [ideal_low_pass( X, r=r )]
    # #     Y.append( ifft_image( Y[-1] ) )
    # #     Y.append( log1p( np.abs( Y[-1] )**2 ) )
    # #     # crop the log of the low pass to get rid of the noise outside
    # #     # of the cell (zero out external points using boundary array)
    # #     crop_zeros = bnd * Y[-1]
    # #     crop_zeros.resize( Y[-1].shape )
    # #     # append the cell with zeros cropped [-4]
    # #     Y.append( crop_zeros ) 
    # #     m = crop_zeros[np.where( crop_zeros!=0 )].mean()
    # #     # everything above m [-3]
    # #     Y.append( np.ma.masked_less_equal( crop_zeros, m ) )
    # #     # just the boolean mask [-2]
    # #     Y.append( np.ma.masked_less_equal( crop_zeros, m ).mask )
    # #     Y.append( m )
    # #     lowpass[r] = Y

    # if not sage_session:
    #     freq = lowpass.keys()
    #     freq.sort()
    #     fig = figure()
    #     #title( r'Lowpass filter of RBC' )
    #     i = 1
    #     for r in freq[2:]:
    #         T = lowpass[r][-2]
    #         fig.add_subplot(1,3,i)
    #         ax = fig.gca()
    #         ax.imshow( T )
    #         ax.set_title( r'Lowpass filter' )
    #         ax.set_xlabel( r'First '+str( r*100 )+'% of modes' )
    #         i += 1
    #     fig.show()
        
    #  #Ylist= [ ideal_low_pass( X, r=i*0.05 ) for i in range(1,5) ]
    # # Yinv = [ ifft_image( Y ) for Y in Ylist ]
    # # # fft of original
    # # PX = np.abs( X )**2
    # # # low pass filtered images
    # # logPY = [ log1p( np.abs( Y )**2 ) for Y in Yinv ]
    # # cropPY = []
    # # for Y in logPY:
    # #     crop = bndy * Y
    # #     crop.resize( Y.shape )
    # #     cropPY.append( crop )

    # # #inverse transform. 
    # # Yinv = fftpack.ifft2( Y ) 
    # # PYinv = [ np.abs( Y )**2 for Y in Yinv ]

    # if sage_session:
    #     print "plotting with sage..."
    #     #original
    #     P=list_plot3d( image, texture='jet', frame_aspect_ratio=(1,1,1) )
    #     P.show()
    #     # for y in PY:
    #     #     list_plot3d( y, texture='jet', frame_aspect_ratio=(1,1,1/5) )
    
    # # fig = figure()
    # # # original image
    # # fig.add_subplot( 231 )
    # # imshow( image )
    # # # FT of image
    # # fig.add_subplot( 234 )
    # # imshow( log1p( PX ) )
    # # # filter
    # # fig.add_subplot( 232 )
    # # imshow( log1p( np.abs( Y )**2 ) )
    # # # filtered image
    # # fig.add_subplot( 235 )
    # # imshow( log1p( np.abs( Yinv )**2 ) )
    
    # # fig.show()
