import subprocess, os
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import numpy# as np
import pp
import time
from collections import deque
from scipy.stats import linregress
# Python image libary
try:
    from PIL import Image  
except ImportError:
    print "Python Image Library not installed"
    raise


slash = '/'

def run_chomp( fname, savename ):
    """
    Call chomp to compute the betti numbers of the image in file fname.

    See http://chomp.rutgers.edu
    """
    cmd = ['chomp', fname]
    try:
        with open( savename, 'w' ) as fh:
            p = subprocess.call( cmd, stdout=fh )
    except:
        print "subprocess returned with code", p

def run_mse( fname, dim=1, **args ):
    """
    Input:
    -----

    fname : path to 1D data file containing time series to
    analyze. (This is slightly ugly since it reads in the 3D numpy
    array, extracts the proper betti time series, save it to disk,
    then run MSE (sigh).)

    **args:
    -----

    Eventually will correspond to MSE args (type 'mse --help' for more
    info).
    """
    # read betti time series file
    arr = numpy.load( fname )
    arr = arr[:,dim,:]
    # write 2D array to tmpfile
    fs = fname.split('/')[:-1]
    fs.append( 'tmpfile' )
    tmpfile= slash.join( fs )
    numpy.savetxt( tmpfile, arr[1] )

    # form the command to pass to subprocess
    outfile = fname.strip('.npy') + '_H'+str(dim)+'.mse'
    cmd = 'mse <'+ tmpfile + '>'+outfile
    try:
        p = subprocess.Popen( cmd, shell=True )
    except:
        print "Problems calling MSE"
        raise

def mse_converter( fname ):
    """
    Convert data in .mse file to a more friendly format. 
    """
    lines = []
    with open( fname ) as fh:
        for line in fh.readlines():
            if len(line) > 1:  # avoid empty lines
                if line.startswith('m'):
                    continue
                # strip off \n and split on tabs
                line = line.strip().split( '\t' )
                lines.append( ( float(line[0]), float(line[1]) ) )
    return numpy.array( lines )

def png2chomp( fname ):
    """
    Convert a numpy array to a text file with lines ( , , ) format for
    chomp. Note: suffix for chomp-readable file must be 'cub' (for
    cubicle complex).
    """
    # open PNG with python image library
    im = Image.open( fname )
    arr = numpy.asarray( im )
    # Find where pixels are black. 255 == white. 
    w = numpy.where( arr != 255 )    
    del arr
    # filter the rgb format
    w2 = numpy.where( w[2] == 0 )[0]
    newarr = numpy.vstack( ( w[0][w2], w[1][w2] ) ).T
    chfile = fname.strip('png') + 'cub'
    array2chomp( newarr, chfile )
            
def array2chomp( arr, savename ):
    """
    Convert an array to chomp format, ( , , ). Write the resulting
    column of numbers to disk.
    """
    rows = map( lambda x: str(x)+'\n', map( tuple, iter( arr ) ) ) 
    with open( savename, 'w' ) as fh:
        fh.writelines( rows )

def pix2array( fname, dim=2 ):
    """
    Convert a PIX file with entries ( , , ) to an array.
    """
    rows = []
    with open( fname ) as fh:
        if dim == 2:
            for line in fh.readlines():
                x = line.strip().split( ',' )
                rows.append( [int( x[0][1:] ), int( x[1][:-1] )] )
        elif dim == 3:
            for line in fh.readlines():
                x = line.strip().split( ',' )
                rows.append( [int( x[0][1:] ), int( x[1] ), int( x[2][:-1] )] )
    return numpy.array( rows, dtype=numpy.int )

def stack_images( path, num_frames=50 ):
    """

    THIS IS SLOW. MAYBE WE SHOULD CREATE A SPECIAL DATA TYPE THAT HOLDS
    THE CHOMP TUPLE(?)
    
    Takea sequence of frames and concatenate the PIX files
    'vertically'. The bottom image has 1 appended to all coordinates,
    the next has 2 appened, etc.

    For example, path='/Users/jesseberwald/data/rbc/chomp_betti/mean50/contours/rbcC1_'

    Returns numpy array of stacked images.
    """
    # read in the frames
    frames = []
    for i in xrange( 100, 100+num_frames ):
        f = path + str(i) + ".pix"
        frames.append( pix2array(f) )
    # append the z values to each frame in the sequence
    for i,x in enumerate( frames ):
        tmp = numpy.empty( (x.shape[0], 1), dtype=numpy.int )
        tmp.fill( i )
        frames[i] = numpy.hstack( (x, tmp) )
    return numpy.vstack( frames )
                    
def extract_betti( fname ):
    """
    Read the betti numbers from the file containing the output from chomp.
          .cbetti files hold the complete output from the chomp program
          .betti files hold the truncated output as produced in 
    """
    # open and read chomp-produced data file
    with open( fname + '.cbetti', 'r' ) as fh:
        lines = fh.readlines()

    # grab the line with the Betti numbers
    for line in lines:
        if line.startswith( 'Betti' ):
            # keep only the numbers
            betti_numbers = line.strip().split()[2:] 

    max_dim = len( betti_numbers )
    # open new
    with open( fname + '.betti', 'w' ) as fh:
        for i in range( max_dim ):
            line = str(i) + ' ' + betti_numbers[i] +'\n'
            fh.write( line )

def read_betti_dir( fdir ):
    """
    Read all .betti files in a directory and organize them for analysis.
    """
    dlist = os.listdir( fdir )
    betti_list = [ f for f in dlist if f.endswith( '.betti' ) ]

    # keep the frame numbers organized in a dict ?
    #betti = {}
    # nah, just list them
    betti_arr = []
    for b in betti_list:
        bnums = numpy.loadtxt( fdir+b, dtype=numpy.uint8 )
        betti_arr.append( bnums )
    betti_arr = numpy.asarray( betti_arr )
    return betti_arr.T
    
def plot_betti( barr, cell=1, savedir=None, dim=0, fig=None,
               total_cells=2, color='b' ):
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

def plot_hist( data, cell_num=1 ):

    fig = plt.figure()
    ax = fig.add_subplot(111)
    dmean = round( data.mean(), 1 )
    dvar = round( data.var(), 1 )
    n, bins, patches = ax.hist( data, 10, normed=1, facecolor='green',
                                alpha=0.75,  label=r"Mean="+str(dmean)+"\nVar="+str(dvar) )
    
    ## # add a 'best fit' line
    ## y = mlab.normpdf( bins, mu, sigma)
    ## l = plt.plot(bins, y, 'r--', linewidth=1)

    xmin = data.min()
    xmax = data.max()
    ax.set_xlabel(r'$H_1$ generators')
    ax.set_ylabel('Probability')
    ax.set_title(r"Distribution of $H_1$ Generators, Cell "+str(cell_num) )
    ax.axis([xmin-5, xmax+5, 0, 0.1])
    ax.grid(True)
    ax.legend()
    return fig

def plot_spectrum( data ):
    """
    Plot the power spectrum of the data.
    """
    d = data[1]
    # rfft gives positive frequecies. Square to get power spectrum.
    fp = numpy.absolute( numpy.fft.rfft( d ) )**2
    freq = numpy.fft.fftfreq( d.shape[-1] )
    n = len(fp)

    # reshape stuff a bit. keep only positive freqs.
    fp = fp[1:-1]
    freq = freq[1:n-1]
    lrslope = linregress( numpy.log(freq[30:]), numpy.log(fp[30:]) )[0]
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.loglog( freq, fp, label="Lin. reg.="+str(round( lrslope,1 )) )
    ax.legend( loc='lower left' )
    return fig
 
if __name__ == "__main__":

    import optparse, sys, time

    plot_help = "Plot betti numbers previously computed."
    compute_help = "Compute betti numbers for each frame."
    threshold_help = "Which threshold to analyze. Factors chosen from "\
                     "{50,75,100,150} * mean. [100]"
    mode_help = "Mode to analyze. [1]"
    debug_help = "Toggle debug mode. [False]"

    parser = optparse.OptionParser()
    parser.usage = "python chomp_bett.py [options]"

    parser.add_option("--plot", "-p",
                      help=plot_help,
                      action="store_true",
                      dest="plot",
                      default=False)
    parser.add_option("--compute", "-c",
                      help=compute_help,
                      action="store_true",
                      dest="compute",
                      default=False)
    parser.add_option("--threshold", "-t",
                      help=threshold_help,
                      type="int",
                      action="store",
                      dest="threshold",
                      default=100)
    parser.add_option("--mode", "-m",
                      help=mode_help,
                      type="int",
                      action="store",
                      dest="mode",
                      default=1)
    parser.add_option("--debug", "-d",
                      help=debug_help,
                      action="store_true",
                      dest="debug",
                      default=False)

    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()

    threshold = "mean" + str( options.threshold ) + "/"
    cell_dir = "new_5/"
    print "Working on cell in directory", cell_dir
 
    #prefix = '/Users/jesseberwald/data/rbc/chomp_betti/'+threshold+'bitmaps/'
    prefix = '/data/jberwald/wyss/data/rbc/chomp_betti/mean100/contours/C1/'
    #prefix = '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'contours/'
    sprefix =  '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti/contours/'
    cprefix = sprefix + threshold+'contours/'
    fname_prefix = 'rbcC'

    if options.debug:
        print "pefix", prefix
        print "sprefix", sprefix
        print "cprefix", cprefix 
    # this was used in the C_4 directory when two cells were
    # present. Now every dir. has only 1 cell
    #cell_nums = [1,2]

    # # meanN path
    # if not os.path.exists( sprefix ):
    #     # maybe betti/ level does not exist
    #     if not os.path.exists( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti' ):
    #         os.mkdir( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti' )
    #         if not os.path.exists( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti/contours/' ):
    #              os.mkdir( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti/contours/' )
    #     else:
    #         # maybe contours/ level doesn't exist
    #         if not os.path.exists( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti/contours/' ):
    #             os.mkdir( '/Users/jesseberwald/data/rbc/cells/'+cell_dir+threshold+'betti/contours/' )
            

    if options.compute:
        #ncpus = 12
        pool = pp.Server( )#ncpus=ncpus)
        jobs = []
        ncpus = pool.get_ncpus()
        print "Created job server with", ncpus, "workers."
        
  #      dims = [0,1]
        num_frames = 5000
        # compute betti number with Chomp.xs Run the PNG extraction in
        # parallel chunks of size <ncpus> to speed stuff up
        print "Computing Betti numbers..."
        # C1 and C2, this is not the larger study
        for cell in [1,2]:
            pix_prefix = prefix + fname_prefix + str( cell ) + "_"
            for i in xrange(100, num_frames+1, ncpus):
                for k in range( i, i+ncpus ):
                    if k < num_frames+1:
                        # convert PNG to data readable by chomp (see function)
                        #fname = fname_prefix + '_'
                        #pix_name = prefix + fname_prefix + str(k) + "_M"+str( options.mode )+".png"
                        pix_name = pix_prefix +  str(k) +".png"

                        if options.debug:
                            print "pix name ", pix_name
                        
                        #bitmap_name = prefix + fname + str(i) + '.txt'
                        #chomp_name = png2chomp( pix_name )
                        jobs.append( pool.submit( png2chomp,
                                                  args=( pix_name, ),
                                                  depfuncs=( array2chomp, ),
                                                  modules=( "numpy", "Image" ) ) )
                pool.wait()
        pool.print_stats()             
      #  pool.destroy()
        
        # Now compute the betti numbers with Chomp
        # for j in cell_nums:     # only if using cells in C_4
        cjobs = []
        for cell in [1,2]:
            pix_prefix = prefix + fname_prefix + str( cell ) + "_"
            for i in xrange(100, num_frames+1):
                for k in range( i, i+ncpus ):
                    if k < num_frames+1:
                        #fname = fname_prefix + str(j)+'_'
                        #chomp_name = prefix + fname_prefix + str(i) + "_M"+str( options.mode )+".pix"
                        chomp_name = prefix + fname_prefix + str( cell ) + "_" + str(k) + ".cub"
                        # run chomp on PIX file
                        # savename = sprefix + "betti/C" +str(j) +"/"+ fname + str(i) + ".cbetti"
                        #savename = sprefix + fname_prefix + str(i) + "_M" + str(options.mode) + ".cbetti"
                        savename = prefix + fname_prefix + str( cell ) + "_" + str(k) + ".cbetti"
                        # run_chomp( bitmap_name, savename )
                        cjobs.append( pool.submit( run_chomp,
                                                   args=( chomp_name, savename ),
                                                   modules=( "subprocess", ) ) )
                pool.wait()
        pool.print_stats()
        # we're done with pp, kill the server
        pool.destroy()
        # Extract and save just the betti number from Chomp's output
        print "Extracting Betti numbers..."
        # For C1 and C2
        for cell in [1,2]:
            betti_prefix = prefix + fname_prefix + str( cell ) + "_"
            for i in xrange(100, num_frames):
                # betti_file = sprefix + "betti/C" +str(j) +"/" + fname_prefix + str(j) + '_' + str(i)
                #betti_file = sprefix + fname_prefix + str(i) + '_M' + str( options.mode )
                betti_file = betti_prefix + str(i) 
                extract_betti( betti_file )

            
        # Concatenate the betti numbers and save to file (use .npy since
        # 3d array is problematic with savetxt).
        # for j in cell_nums:
        #     #cells = []
        #     # fdir = sprefix + "betti/C"+str(j) +"/"
        #     fdir = sprefix + "C"+str(j) +"/"
        #     ba = read_betti_dir( fdir )
        #     #cells.append( ba )
        #     arrname = sprefix + 'betti_arrC' + str(j)
        #     numpy.save( arrname, ba )
        #     for k in dims:
        #         # compute MSE at the same time
        #         run_mse( arrname+'.npy', dim=k )

    if options.plot:
        # load the array of betti numbers for each cell
        cells = [numpy.load( sprefix + 'betti_arrC'+str(i)+'.npy' )
                 for i in cell_nums]
        # create a plot for each dimension (0 and 1)
        dims = [0,1]
        for dim in dims:
            fig = plt.figure()
            ax = fig.add_subplot(111)
            ax.set_title( "$H_"+str(dim)+"$ Over Time" )
            for i,cell in enumerate( cells ):
                data = cell[:,dim,:]
                cm = data[1].mean()
                ax.plot( data[1], 'o-', lw=1.5, ms=2, label='Cell '+str(i) )
                ax.plot( len(data[1])*[cm], '--', label='Mean, cell '+str(i) )
                ax.set_xlabel( "Time", fontsize=18 )
                ax.set_ylabel( "$H_"+str(dim)+"$", fontsize=18 )
            ax.legend()
            fig.savefig( sprefix + "rbcH"+str(dim)+".png" )
            fig.show()
        # plot a histogram of the *Betti_1* time series for each
        # cell. Include mean and variance.
        for i,cell in enumerate( cells ):
            data = cell[:,1,:]
            fig2 = plot_hist( data[1], cell_num=i+1 )
            fig2.savefig( sprefix + "rbcH"+str(dim)+"_C"+str(i)+"_hist.png" )
            # plot power spectrum for H_1 time series
            fig2 = plot_spectrum( data )
            fig2.savefig( sprefix + "rbcH"+str(dim)+"_C"+str(i)+"_spectrum.png" )
            fig2.show()

        # loop over cells
        fig3 = plt.figure()
        for j in dims:
            fig3.clf()
            ax = fig3.add_subplot(111)
            colors = ['b','r']
            for i in [1,2]:
                msefile = sprefix + 'betti_arrC'+str(i)+'_H'+str(j)+'.mse'
                arr = mse_converter( msefile )
                ax.plot( arr[:,0], arr[:,1], c=colors[i-1], marker='o', lw=2,
                         label=r'Cell '+str(i) )
                ax.set_title( "Cell "+str(i) )
                ax.set_xlabel( "scale" )
                ax.set_ylabel( "Entropy" )
            ax.legend()
            fig3.savefig( sprefix + "betti_mse_H"+str(j)+".png")
            fig3.show()
