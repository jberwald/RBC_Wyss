import scipy.io
import numpy as np
from matplotlib import pyplot as plt

under = "_"

def analyze_births( cell, death=False ):
    """
    Input:
    -----

    Sequence of betti numbers obtain from images collected over time
    for a single cell.

    Returns array of average birth time for the cell

    """
    # compute average birth times over all images in the sequence
    endpts = {}
    max_len = 0
    # use righthand endpoints of intervals
    if death:
        m = 1
    else:
        m = 0
    for i,img in enumerate(cell):
        endpts[i] = []
        # find unique left or right endpoints of intervals
        uniq = np.unique( img[m] )
        # nmber of generators at each nique endpoint
        for j,u in enumerate( uniq ):
            endpts[i].append( len( np.where( img[m]== u )[0] ) )
        # make sure we have maximum number of unique endpts
        if max_len < len( uniq ):
            max_len = len( uniq )
            endpt_times = uniq
    # If no endpts at higher radii, fill with zeros
    for i in endpts:
        if len( endpts[i] ) < max_len:
            diff = max_len - len( endpts[i] )
            endpts[i].extend( diff*[0] )
    return endpts, endpt_times

def avg_generators( cell ):
    """
    Input:
    -----

    Dict. with lists of tuples recording time of generator birth and
    number of generators. 

    Return average over the number of generators at each birth time    
    """
    arr = np.array( cell.values() )
    return arr.mean( axis=0 )

def std_generators( cell ):
    """
    Input:
    -----

    Dict. with lists of tuples recording time of generator birth and
    number of generators. 

    Return std over the number of generators at each birth time    
    """
    arr = np.array( cell.values() )
    return arr.std( axis=0 )

    
def plot_dim_avg( avg_gens, dim=0, fig=None ):
    """
    For each cell plot the average birth time of generators.
    """
    if not fig:
        fig = plt.figure()

    # birth times are stored in second tuple entry. They should all be
    # the same, so we just grab one array.
    ax = fig.add_subplot(1,1,1)
    ax.set_title( "Dimension " + str( dim ) + " ($H_"+str( dim )+"$)")
    # Average umber of generators stored in first entry of tuple. 
    for cell in avg_gens :
        ax.set_xlabel( 'time' )
        ax.set_ylabel( 'Average # of generators' )
        nx = avg_gens[cell][0]
        ny = avg_gens[cell][1]
        err = avg_gens[cell][2]
        #ax.plot( nx, ny, marker='o', lw=2, label="Cell "+str(cell) )
        ax.errorbar( nx, ny, yerr=err, fmt='-o', lw=2, ecolor='r', label="Cell "+str(cell) )
    ax.legend( loc='upper right' )
    
    
if __name__ == "__main__":

    prefix = "/Users/jesseberwald/data/rbc/betti_r1.5/rbc_betti"
    fig_prefix = "/Users/jesseberwald/data/rbc/figures/maxr1.5/"
    death = True
    dims = [0,1]

    img_size = {}

    # loop over the cells. collect all [birth,death) intervals for each image in list
    dim_avg = {}
    for k in dims:
        betti_ints = {}
        for i in range(1,3):
            # the betti numbers are stored in .MAT arrays. We extract
            # these using scipy.io.loadmat. This gives a dict, with data
            # stored in the 'interval_mat' key. Kinda ugly, but whatever...
            betti_ints[i] = []
            img_size[i] = []
            for j in range( 100, 201 ):
                fname = under.join( [prefix, "C"+str(i), str(j), "dim"+str(k)+".mat"] )
                betti_ints[i].append( scipy.io.loadmat( fname )['interval_mat'] )
                img_size[i].append( betti_ints[i][-1].shape[1] )

        # find the birth or death time for generators for each cell
        cell_times= {}
        times = {}
        for cell in betti_ints:
            cell_times[cell], times[cell] = analyze_births( betti_ints[cell], death=death )
            
        # compute average number of generators and stdev at each step in time
        avg_gens = {}
        for cell in cell_times:
            avg_gens[cell] = ( times[cell],
                               avg_generators( cell_times[cell] ),
                               std_generators( cell_times[cell] ) )      
        dim_avg[k] = avg_gens

    for d in dim_avg:
        fig = plt.figure( d+1 )
        plot_dim_avg( dim_avg[d], dim=d, fig=fig )
        fig.savefig( fig_prefix + "avg_betti_dim"+str(d)+".png" ) 
    plt.show()
