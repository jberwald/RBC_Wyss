#!/usr/bin/python
import sys, os, subprocess
"""
Make a movie from PNG files. Assumed that this script resides in the
folder containing a sequence of frames (PNG files).
"""

slash = "/"

def sort_files( dlist ):
    """ Sort files in dlist according to the order
    of creation.  """
    return sorted(dlist, key=order)

def order( line ):
    """split line and return the index of the partial cover"""
    return int( line.rstrip().split('.')[-2] )

def prefix_files( dlist, prefix ):
    """Return only files beginning with <prefix>"""
    plist = [x for x in dlist if x.startswith( prefix )]
    return plist

def suffix_files( dlist, suffix ):
    """Return only files ending with <suffix>"""
    slist = [x for x in dlist if x.endswith( suffix )]
    return slist

def patch_cover( dlist, prefix, suffix='png'):
    """Return image (of type <suffix>) that begin with <prefix>"""
    a = prefix_files( dlist, prefix )
    return suffix_files( a, suffix ) 

def run(files, output):

    command = ('mencoder',
               files,
               '-mf',
               'type=png:h=800:w=400:fps=8',
               '-ovc',
               'lavc',
               '-lavcopts',
               'vcodec=avi',
               '-oac',
               'copy',
               '-o',
               output)

    print command

    print "Making movie from " + prefix + " barcode images..."
    #os.spawnvp(os.P_WAIT, 'mencoder', command)
    subprocess.call( command )


if __name__ == "__main__":

    # prefix for cover files
    try:
        # file name prefix
        prefix = sys.argv[1]
    except IndexError:
        prefix = "cantorMF4_2D_0.9"

    dlist = os.listdir( '.' )
    flist = patch_cover( dlist, prefix )
    flist= sort_files( flist )

    fname = prefix + ".clist"
    with open( fname, 'w' ) as fh:
        for x in flist:
            fh.write( x + "\n" )
    
    #files = 'mf://' + prefix + '*.png'
    files = 'mf://@' + fname
    output = prefix + '.avi'
    #output = prefix + '.mpeg'

    run( files, output )
