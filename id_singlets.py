import tifffile as tf
import numpy as np
import scipy as sp
import scipy.ndimage as img
import scipy.spatial as spatial
import argparse
import pandas as pd
from annotate import *

def id_singlets(image_handle, channel, threshold, min_size, max_size, min_distance):
    orig = tf.TiffFile(image_handle)
    red, green = orig
    spots = (orig[channel]).asarray() > threshold
    labels, n_spots = img.label(spots)
    print "Found {} spots.".format(n_spots)

    # keep singlets (i.e. cells separated by at least min_distance)
    centers = np.array(img.center_of_mass(labels, labels, xrange(1, n_spots+1)))
    distances = spatial.distance_matrix(centers, centers)
    doublets = np.any((distances > 0) & (distances < min_distance), axis=1)
    slices = img.find_objects(labels)
    for i in np.flatnonzero(doublets):
        this_slice = labels[slices[i]]
        this_slice[this_slice == i+1] = 0
    labels, n_singlets = img.label(labels)
    print "Retaining {} spots past distance filter.".format(n_singlets)

    # apply size filter
    sizes = img.sum(spots, labels, xrange(1, n_singlets+1))
    cells = np.flatnonzero((sizes >= min_size) & (sizes < max_size)) + 1
    slices = img.find_objects(labels)
    # iterate over the indices of the spots; blank rejected spots
    for i in xrange(1, n_singlets+1):
        if i in cells: continue
        this_slice = labels[slices[i-1]]
        this_slice[this_slice == i] = 0
    labels, n_cells = img.label(labels)

    indices = range(1, n_cells+1)
    centers = np.array(img.center_of_mass(labels, labels, indices))
    green_scores = img.sum(green.asarray(), labels, indices)
    red_scores = img.sum(red.asarray(), labels, indices)
    areas = img.sum(spots, labels, indices)
    scores = pd.DataFrame({'x': centers[:,1],
                           'y': centers[:,0],
                           'area': areas,
                           'green_intensity': green_scores,
                           'red_intensity': red_scores})
    return scores


def annotate_ld(livedead_fn, annotate_fn, channel, singlets, cell_radius):
    im = tf.TiffFile(livedead_fn)[channel].asarray()
    surface = surface_from_array(im)
    rows = len(singlets)
    foo = lambda x: np.asarray(x).reshape(-1,1)
    cell_matrix = np.hstack([foo(singlets['x']), foo(singlets['y']), np.ones((rows, 1))])
    surface = annotate_cells(surface, cell_matrix, cell_radius, (1.0, 0.0, 0.0), linewidth=2)
    img = PIL_from_surface(surface)
    img.save(annotate_fn)

def main():
    parser = argparse.ArgumentParser(description='Process live/dead images.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('livedead', metavar='stitched_livedead.tif',
                        help='Stitched live/dead image, in TIFF format.')
    parser.add_argument('--threshold', '-t', required=True, type=int,
                        help='Threshold on [0, 4096) for live/dead image.')
    parser.add_argument('--channel', '-c', type=int, default=1,
                        help='Which image channel to analyze (0-based index)')
    parser.add_argument('--cell-size', '-s', nargs=2, default=(10, 32), type=int,
                        metavar=('min', 'max'),
                        help='Recognize cells with areas on [min, max), in live/dead pixel units.')
    parser.add_argument('--distance', '-d', default=25, type=float,
                        help='Minimum acceptable separation (in pixels) between cells.')
    parser.add_argument('--annotate', '-a', type=str,
                        help='Filename for annotated live/dead image.')
    parser.add_argument('--output', '-o', default='livedead_results.txt',
                        help='Filename to save information about the live/dead image.')
    args = parser.parse_args()

    print 'Identifying singlets...'
    singlets = id_singlets(args.livedead, threshold=args.threshold, min_distance=args.distance,
                           min_size=args.cell_size[0], max_size=args.cell_size[1],
                           channel=args.channel)
    singlets.to_csv(args.output, index=False)
    print 'Found {} singlets.'.format(len(singlets))

    if args.annotate:
        r = np.sqrt(args.cell_size[1]/np.pi)
        annotate_ld(args.livedead, args.annotate, args.channel, singlets, 2*r)

if __name__ == '__main__':
    main()
