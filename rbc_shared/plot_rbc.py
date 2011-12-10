from pylab import *
#import raw_cell as R
import cPickle as pkl

def read_cell_betti_numbers( fname ):
    """
    Assumes that betti numbers in chomp/ directory have been collated
    into one pickle file (per threshold value). This is a dictionary
    keyed by cell name.
    """
    with open( fname, 'r' ) as fh:
        betti = pkl.load( fh )
    return betti

def write_betti_pkl( arr, fname ):
    """
    Just in case we want to write the dictionary of betti numbers to disk.
    """
    with open( fname, 'w' ) as fh:
        pkl.dump( arr, fname )

def find_means( betti_arrs, dim=0, r=None ):
    """
    """
    means = {}
    for k in betti_arrs:
        a = betti_arrs[k]
        a_dim = a[:,dim,:]
        means[k] = a_dim[1].mean()
    return means

if __name__ == "__main__":

    # dimension to analyze
    dim = 1

    new_cells = 'new_betti_thresh1.25.pkl'
    old_cells = 'old_betti_thresh1.25.pkl'

    new_betti = read_cell_betti_numbers( new_cells )
    old_betti = read_cell_betti_numbers( old_cells )

    new_means = find_means( new_betti, dim=dim )
    newavg = average( new_means.values() )
    old_means = find_means( old_betti, dim=dim )
    oldavg = average( old_means.values() )

    fig = figure()
    ax = fig.gca()

    for k in new_betti:
        ax.plot( new_betti[k][:,dim,:][1], '-',  lw=2, alpha=0.2 )
    for k in old_betti:
        ax.plot( old_betti[k][:,dim,:][1], '-', lw=2, alpha=0.2 )
        ax.hlines( newavg, 0, 5000, linestyle='dashed', linewidth=2, label=r'New cell avg $H_'+str(dim)+'$')
        ax.hlines( oldavg, 0, 5000, linestyle='dashdot', linewidth=2, label=r'Old cell avg $H_'+str(dim)+'$' )
        ax.legend( loc='upper right' )
        fig.show()
