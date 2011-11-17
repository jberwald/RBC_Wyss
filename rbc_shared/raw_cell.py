import numpy
from matplotlib import pyplot as plt
import os, shutil
#from jjb.chomp import chomp_betti as cb

def plot_raw_cell( fname, **args ):
    """
    Read in an ASCII file containing raw data for a red blood cell.

    Original files contain 5000 frames. Limit the frames loaded using
    max_frames arg. Number of frames = 5000 - <max_frames>. Therefore,
    number of rows = 5 implies that the last 5 rows of fname were read
    in (eg., 4995-4999).

    NOTE: bounday matrix is hardcoded for now.
    """
    fargs = { 'savename': None,
             'skiprows': 4997,
             'bnd' : 'boundary_Nov_new110125',
             'plot' : True
             }
    fargs.update( args )
    # set max frames
    skiprows = fargs['skiprows']

    # load boundary matrix for dimensions of cell array
    bnd_arr = numpy.loadtxt( fargs['bnd'] )
    nx,ny = bnd_arr.shape
    frames = load_rbc( fname, skiprows, nx, ny )
    num_cols = len( frames ) # for subplot 

    # compute mean over all extracted frames
    avg = find_mean( frames )
    print "avg", avg

    masked = [ numpy.ma.masked_less( f, avg ) for f in frames ]
    
    # save to disk
    if fargs['savename']:
        for i in range( num_cols ):
            sname = savename + '_' + str(i)
            numpy.save( sname, frames[i] )
    # plot the frames
    if fargs['plot']:
        fig = plt.figure()
        for i in range( num_cols ):
            sax = fig.add_subplot( 2, num_cols, i+1 )
            ma = numpy.ma.masked_less_equal( frames[i], 0 )
            sax.imshow( ma )
            #sax.imshow( frames[i] )
        for i in range( num_cols ):
            sax = fig.add_subplot( 2, num_cols, 4+i )
            sax.imshow( masked[i] )
        fig.savefig( fname+'_frames'+str(num_cols)+'.png' )
        fig.show()

def load_rbc( fname, skiprows, nx, ny ):
    """
    Returns a block frames from <fname> cell data. Reshapes array
    according to boundary data determined from <bnd> arg in
    plot_raw_data().
    """
    C = numpy.loadtxt( fname, skiprows=skiprows )    
    cell_frames = [ C[i].reshape(( nx,ny )) for i in range( 5000-skiprows ) ]
    return cell_frames
    
def extract_frames( fname, bnd ):
    """
    Each line of fname contains a raveled matrix. Read each line,
    reshape it based on array stored in <bnd> file. Save to file.
    """
    fh = open( fname, 'r' )
    bnd_arr = numpy.loadtxt( bnd )
    nx, ny = bnd_arr.shape
    # save files in a special folder
    part = fname.rpartition( '/' )
    cell_name = part[-1]
    savedir = part[0] + '/' + 'frames/'

    print 'savedir', savedir
     
    # loop over lines in <fname> and save each a nx x ny array
    #k = 0
    for k, line in enumerate( fh.readlines() ):
        arr = numpy.fromstring( line, sep='\t' )
        arr.resize( (nx,ny) )
        numpy.save( savedir + cell_name + '_' + str(k), arr )

def frames2png( fdir ):
    """
    Convert matrices in directory <fdir> to images.
    """
    dlist = os.listdir( fdir )
    dlist = [ f for f in dlist if f.endswith('npy') ]
    fig = plt.figure()
    ax = fig.gca()
    ax.hold( False )
    for f in dlist:
        arr = numpy.load( fdir + f )
        ma = numpy.ma.masked_less_equal( arr, 0 )
        ax.imshow( ma )
        # strip the .npy off
        fig.savefig( fdir + f[:-3] + 'png' )
        ax.clear()

def thresh2png( fdir, value ):
    """
    Convert thresholded matrices in directory <fdir> to images.
    """
    dlist = os.listdir( fdir )
    thresh = 'thresh'+str(value)
    dlist = [ f for f in dlist if f.endswith('npy') and thresh in f ]
    fig = plt.figure()
    ax = fig.gca()
    ax.hold( False )
    for f in dlist:
        arr = numpy.load( fdir + f )
        #ma = numpy.ma.masked_less_equal( arr, 0 )
        ax.imshow( arr )
        # strip the .npy off
        fig.savefig( fdir + f[:-3] + 'png' )
        ax.clear()

def threshold_frame( arr, pct ):
    """
    Threshold an cell array at <pct>*mean. <pct> should be entered as
    a decimals (eg., 100% = 1.0).

    1 if > <pct>*mean
    0 if <= <pct>*mean
    """
    # mask the array as a percent of the mean (excluding al zeros first)
    ma = numpy.ma.masked_less_equal( arr, 0 )
    return numpy.ma.masked_less_equal( arr, pct*ma.mean() )

def threshold_all( fdir, pct=1.0 ):
    """
    Threshold all RBC images in directory <fdir>. The originals are
    stored in NPY files. Save the thresholded arrays at them as
    boolean arrays after masking.
    """
    dlist = os.listdir( fdir )
    dlist = [ f for f in dlist if f.endswith('npy') and 'thresh' not in f ]
    for f in dlist:
        arr = numpy.load( fdir + f )
        thresh = threshold_frame( arr, pct )
        numpy.save( fdir + f[:-4] +'_thresh'+str(pct)+'.npy', thresh.mask )

def thresh2cub( fdir, value ):
    """
    Convert thresholded image matrices to Chomp CUB format.
    
    Calls bool2cub since thresholded matrices are stored as bools.
    """
    dlist = os.listdir( fdir )
    thresh = 'thresh'+str(value)
    dlist = [ f for f in dlist if f.endswith('npy') and thresh in f ]
    for f in dlist:
        savename = fdir + f[:-3] + 'cub'
        arr = numpy.load( fdir + f )
        bool2cub( arr, savename )
        
def bool2cub( arr, savename ):
    """
    Convert a thresholded aray, stored as a boolean, to a
    Chomp-readable array. Write the new array to disk.
    """
    c = arr.astype( 'uint' )
    w = numpy.where( c==0 )
    # zip locations of thresholded values to get coords (2D)
    z = zip( w[0], w[1] )
    coords = [ str(x)+'\n' for x in z ]
    with open( savename, 'w' ) as fh:
        fh.writelines( coords )

def find_mean( frames ):
    return numpy.mean( [ f.mean() for f in frames] )

def run_chomp( fdir, value ):
    """
    """
    # grab CUB files
    dlist = os.listdir( fdir )
    thresh = 'thresh'+str( value )
    dlist = [ f for f in dlist if f.endswith('cub') and thresh in f ]
    for f in dlist:
        savename = fdir + 'chomp/' + f[:-4] 
        cb.run_chomp( fdir + f, savename+'.cbetti' )
        cb.extract_betti( savename )

def rename_cub_files( fdir ):
    """
    chomp cannot deal with two "dots" in a file name. 
    """
    dlist = os.listdir( fdir )
    dlist = [ f for f in dlist if f.endswith('cub') ]
    newlist = []
    for f in dlist:
        part = f.partition( '.' )
        newf = part[0] + part[-1]
        shutil.move( fdir + f, fdir + newf )
    

def plot_betti( barr, cell=1, savedir=None, dim=0, fig=None,
               total_cells=2, color='b', showplot=False ):
    """
    Plot betti numbers for each frame for a cell. Obtain a time series
    (time=frame number)
    """
    if fig is None:
        fig = plt.figure()
    ax = fig.gca()
    #ax = fig.add_subplot(total_cells, 1, cell_num+1)
    data = barr[:,dim,:]
    ax.plot( data[1], 'o-', color=color, lw=1.5, ms=2 )
    # record title and some stats
    ax.set_title(  'Betti numbers for cell '+str(cell)+\
                 ' (mean='+str( round(data[1].mean()) )+')' )
    ax.set_xlabel( 'Frame' )
    ax.set_ylabel( r'$H_{'+str(dim)+'}$' )
    if savedir == None:
        fname = './figures_raw/betti_frames_H'+str(dim)+'_cell'+str(cell)+'.png'
    else:
        fname = savedir + '/betti_frames_H'+str(dim)+'_cell'+str(cell)+'.png'
    fig.savefig( fname )
    if showplot:
        fig.show()

def read_betti_dir( fdir ):
    """
    Read all .betti files in a directory and organize them for analysis.
    """
    dlist = os.listdir( fdir )
    betti_list = [ f for f in dlist if f.endswith( '.betti' ) ]

    betti_dict = dir_hash( betti_list )
    frames = betti_dict.keys()
    frames.sort()
    # keep the frame numbers organized in a dict ?
    #betti = {}
    # nah, just list them, but append them in frame order
    betti_arr = []
    for i in frames:
        b = betti_dict[i]
        bnums = numpy.loadtxt( fdir+b, dtype=numpy.uint8 )
        betti_arr.append( bnums )
    betti_arr = numpy.asarray( betti_arr )
    return betti_arr.T

def dir_hash( dlist ):
    """
    """
    files= {}
    for filename in dlist: #os.listdir("C:\\TEMP\\py"):
        basename, extension = filename.split('.')
        # filename: new/old, prefix==id, frame number, threshold value
        celltype, prefix, frame, thresh = basename.split('_')
        files[ int( frame ) ] = filename
    return files


if __name__ == "__main__":

    #fdir = '/data/jberwald/wyss/data/Cells_Jesse/New/frames/'
    fdir = '/data/jberwald/wyss/data/Cells_Jesse/Old/frames/'

    # str_values = [ '09', '10', '12'] #[ '03' , '05', '07' ]
    # for sval in str_values:
    #     thresh2cub( fdir, sval )
        
    values =  [ '09', '10', '12', '03' , '05', '07' ]# [ 0.3, 0.5, 0.7, 0.9, 1.0, 1.2 ] #[ 0.3, 0.5, 0.7,
    # #values = [ '09', '10', '12']
    for val in values:
        #threshold_all( fdir, pct=val )
        #thresh2cub( fdir, val )
        run_chomp( fdir, val )

#    rename_cub_files( fdir )
