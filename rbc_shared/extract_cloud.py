import numpy
from matplotlib import pyplot
import pp
import os

def mask_data( data, level='bg' ):
    # mask background values
    if level == 'mean':
        md = numpy.ma.masked_less_equal( data, data.mean() )
    elif level == 'bg':
        md = numpy.ma.masked_less_equal( data, 0 )
    else:
        print "eek!"
    return md

def map2euclidean( mdata ):
    """
    Convert masked image array to a Euclidean representation.
    """
    nx,ny = numpy.ma.where( mdata != True )
    return numpy.array( zip( nx,ny ) )


def save_frame_as_png( matdata, savepath, fig=None ):
    """
    Save the contour plot of the extracted frame data to a png image.
    """
    if fig is None:
        fig = pyplot.figure()
    fig.clf()
    ax = fig.add_subplot(111, frame_on=False)
    ax.set_axis_off()
    ax.contour( matdata, colors='k' )
    fig.savefig( savepath )
    

def run_extraction( factor, prefix, fprefix, num_frames,
                    plot_bitmaps=False, do_contour=True,
                    old=True, mode=1 ):
    """
    Important options:
    ---------------

    mode_cols : Specify the width of the matrix to extract. Eg., if
    each frame is 65x70, and we want to extract Mode 2, then the image
    variable must contain columns 70 to 2*70. 
    """
    # keep a single figures around so we don't create a new one each time
    fig = pyplot.figure()
    contour_fig = pyplot.figure()
    
    f = factor/100. # multiplication by percentage of mean
    threshold_path = "mean"+str(factor)+"/"  # path addition
    mat_prefix = "matC_"
    rbc_prefix = "rbc_"
    # These prefixes are only used with the cells in C_4 directory.
    # if old:
    #     mat_prefix = "matC_"
    #     rbc_prefix = "rbcC1_"
    # else:
    #     mat_prefix = "matC2_"
    #     rbc_prefix = "rbcC2_"

    # check paths
    if not os.path.exists( fprefix + threshold_path + "contours/" ):
        if not os.path.exists( fprefix + threshold_path ):
            if not os.path.exists( fprefix ):
                os.mkdir( fprefix )
            else:
                os.mkdir( fprefix + threshold_path )
        else:
            os.mkdir( fprefix + threshold_path + "contours/" )

    mode_suffix = "_M"+str( mode )
    fname_prefix = prefix+mat_prefix
    m = mode-1 # extract the mode
    for i in xrange(100, num_frames+1):
        # path to frame
        fname = fname_prefix+str(i)+mode_suffix
        # load image matrices; each frame in a separate file
        image= numpy.loadtxt( fname, delimiter=',' )  ##[:,min_col:max_col]
       # c2 = numpy.loadtxt( fname2, delimiter=',' )[:,:70]

        # mask the data in two ways. First, get rid of
        # background. Next, only keep what's above the mean
        masked_image = mask_data( image, 'bg' )    
       # mc2 = mask_data( c2, 'bg' )    

        # Save contour plots
        if do_contour:
            # Make a contour graph and save it to disk in "Chomp format"
            chprefix = fprefix + threshold_path+'contours/'
            sname = chprefix + rbc_prefix +str(i) + mode_suffix + '.png'
            save_frame_as_png( masked_image, sname, fig=contour_fig )
           # save_frame_as_png( mc2, chprefix + 'rbcC2_'+str(i)+'.png', fig=contour_fig )

        ## mc1 = mask_data( c1, 'mean' )
        ## mc2 = mask_data( c2, 'mean' )

        # Do the original way of saving just matrix coordinates (creates much smaller matrices)
        else:
            threshold1 = masked_image.mean()
  #          threshold2 = mc2.mean()
            w1 = numpy.where( masked_image > f * threshold1 )
 #           w2 = numpy.where( mc2 > f * threshold2 )

            bmap1 = zip( w1[0], w1[1] )
  #          bmap2 = zip( w2[0], w2[1] )

            # chomp prefix with cell numbers (TODO: get rid of hard coding)
            # for j in [1,2]:
            #     cprefix = fprefix + threshold_path+"bitmaps/"
            #     if j == 1:
            #         bm = bmap1
            #     else:
            #         bm = bmap2
                # format bitmap file for chomp
            with open( chprefix+rbc_prefix+'_'+str(i)+'.txt', 'w' ) as fh:
                for x in bm:
                    line = str(x) + "\n"
                    fh.write( str(line) )     
        
        # plot the bitmaps images
        if plot_bitmaps:
            bm1 = numpy.array( bmap1 )
            bm2 = numpy.array( bmap2 )
            fig.clf()
            ax = fig.add_subplot(121)
            ax.set_title( "Cell 1, threshold "+str( factor ) + " * mean" )
            ax.scatter( bm1[:,0],bm1[:,1], s=4 )
            ax2 = fig.add_subplot(122)
            ax2.set_title( "Cell 2, threshold "+str( factor ) + " * mean" )
            ax2.scatter( bm2[:,0], bm2[:,1], s=4 )
            fig.savefig( fprefix + threshold_path +"figures/rbc_frame"+str(i)+".png" )
            #fig.show()

if __name__ == "__main__":

    import optparse, sys, time

    plot_help = "Plot the bitmaps, if that is what we compute. [False]"
    show_help = "Show the plots? [False]"
    mode_help = "Which Mode to analyze. [1]"
    contour_help = "Analyze the contours instead of the bitmaps. [False]"
    frames_help = "Number of frames to extract data from.\nNote: we start at frame 150 "\
                  "to avoid blank frames at the beginning of the sequence. [150]"
    
    parser = optparse.OptionParser()
    parser.usage = "python extract_cloud.py [options]"
    
    parser.add_option("--plot", "-p",
                      help=plot_help,
                      action="store_true",
                      dest="plot",
                      default=False)
    parser.add_option("--show", "-s",
                      help=show_help,
                      action="store_true",
                      dest="show",
                      default=False)
    parser.add_option("--mode", "-m",
                      help=mode_help,
                      type="int",
                      action="store",
                      dest="mode",
                      default=1)
    parser.add_option("--contour", "-c",
                      help=contour_help,
                      action="store_true",
                      dest="contour",
                      default=False)
    parser.add_option("--frames", "-f",
                      help=frames_help,
                      type="int",
                      action="store",
                      dest="frames",
                      default=150)
    
    (options, args) = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    old_dir = None
    new_dir = None

    # Extract and threshold frames for cells 
    cells = [6]
    cell_dirs = []
    for c in cells:
        old_dir = "old_" + str(c) + "/"
       # new_dir = "new_" + str(c) + "/"
        cell_dirs.append( old_dir )
        #cell_dirs.append( new_dir )
    for cell in cell_dirs:
        print "Extracting data in folder", cell, "..."
        print "Using Mode", options.mode

        # show_plot = False
        # plot_bitmaps = False
        num_frames = options.frames
        #prefix = "/Users/jesseberwald/data/rbc/C_4/frames/"
        prefix = "/Users/jesseberwald/data/rbc/cells/"+cell+"frames/"
        fprefix = "/Users/jesseberwald/data/rbc/cells/"+cell
        factors = [100] #[50,75,100]

        # Do some preprocessing to make sure we have created the proper
        # directories to store stuff.
        # check chomp paths
        if not os.path.exists( fprefix + "chomp_betti/" ):
            if not os.path.exists( fprefix  ):
                os.mkdir( fprefix )
            else:
                os.mkdir( fprefix + "chomp_betti/" )
        # meanN path
        for N in factors:
            if not os.path.exists( fprefix + "mean"+str(N)+"/contours" ):
                if not os.path.exists( fprefix+"mean"+str(N)  ):
                    os.mkdir( fprefix+"mean"+str(N) )
                else:
                    os.mkdir( fprefix + "mean"+str(N) + "/contours")

        print "Extracting and thresholding images..."

        # # keep a single figure around so we don't create a new one each time
        # fig = pyplot.figure()

        for factor in factors:
            print "Threshold factor", factor
            run_extraction( factor, prefix, fprefix, num_frames,
                            plot_bitmaps=options.plot, do_contour=options.contour,
                            mode=options.mode )
  
            # map masked array to Euclidean grid
            ## euc1 = map2euclidean( mc1 )
            ## euc2 = map2euclidean( mc2 )

            # fname1 = "/Users/jesseberwald/data/rbc/image_data_threshold0/cloudC1_"+str(i)
            # fname2 = "/Users/jesseberwald/data/rbc/image_data_threshold0/cloudC2_"+str(i)
            # numpy.savetxt( fname1, euc1, fmt='%8d' )
            # numpy.savetxt( fname2, euc2, fmt='%8d' )

            # if show_plot:
            #     fig = pyplot.figure()
            #     ax = fig.add_subplot(121)
            #     ax.imshow( mc1 )
            #     #pyplot.colorbar()
            #     ax = fig.add_subplot(122)
            #     ax.imshow( mc2 )
            #     #pyplot.colorbar()
            #     fig.savefig("/Users/jesseberwald/data/rbc/figures/mean/rbc_frame"+str(i)+".png")

            #     fig2 = pyplot.figure()
            #     ax = fig2.add_subplot(121)
            #     ax.plot( euc1[:,0], euc1[:,1], 'bo', ms=2 )
            #     ax = fig2.add_subplot(122)
            #     ax.plot( euc2[:,0], euc2[:,1], 'bo', ms=2 )
            #     fig2.savefig("/Users/jesseberwald/data/rbc/figures/mean_eucl/rbc_frame"+str(i)+".png")

            # #pyplot.show()

