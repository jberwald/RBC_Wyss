import numpy


def load_rbc( fname, skiprows, nx, ny ):
    """
    Returns a block of frames from <fname> cell data. Reshapes array
    according to boundary data determined from nx, ny (determined from
    corresponding boundary file).
    """
    C = numpy.loadtxt( fname, skiprows=skiprows )    
    cell_frames = [ C[i].reshape(( nx,ny )) for i in range( 5000-skiprows ) ]
    return cell_frames

def extract_frames( fname, bnd ):
    """
    Each line of the cell in <fname> contains a raveled matrix. Read
    each line of <fname>, reshape it based on array stored in <bnd>
    file. Then save each frame to file.
    """
    if os.uname()[0] == 'Linux':
        savefunc = numpy.save
    elif os.uname()[0] == 'Darwin':
        # why the !@$# doesn't Mac save readable .npy files??
        savefunc = numpy.savetxt
    fh = open( fname, 'r' )
    bnd_arr = numpy.loadtxt( bnd )
    nx, ny = bnd_arr.shape
    bnd_arr = bnd_arr.ravel()
    # save files in a special folder
    part = fname.rpartition( '/' )
    cell_name = part[-1]
    # find the subfolder name from the beginning of the cell
    # name. subfolder must have been created already.
    subfolder = cell_name.partition( '-' )[0]
    savedir = part[0] + '/' + 'frames/' + subfolder +'/'

    make_dir( savedir )
    # print 'savedir', savedir
     
    # loop over lines in <fname> and save each as nx x ny array
    #k = 0. 
    fromstring = numpy.fromstring
    for k, line in enumerate( fh.readlines() ):
        arr = fromstring( line, sep='\t' )
        # remove boundary
        arr = bnd_arr * arr
        arr.resize( (nx,ny) )
        #numpy.save( savedir + cell_name + '_' + str(k), arr )
        savefunc( savedir + cell_name + '_' + str(k), arr )

def thresh2png( fdir, value ):
    """
    Convert thresholded matrices in directory <fdir> to PNG images.
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
