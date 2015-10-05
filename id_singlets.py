from tifffile import tifffile as tf
import numpy as np
import scipy as sp
import scipy.ndimage as img
import scipy.spatial as spatial
import argparse
import pandas as pd
from annotate import *

not_singlet_filter = lambda distances, min_distance: np.any((distances > 0) & (distances < min_distance), axis=1)
not_doublet_filter = lambda distances, min_distance: np.sum((distances > 0) & (distances < min_distance), axis=1) != 1


def id_singlets(image_handle, channel, threshold, min_size, max_size, min_distance, reject=not_singlet_filter):
    orig = tf.TiffFile(image_handle)
    if orig.asarray().ndim == 3:
        is_3d = True
        red, green = orig
        spot_image = orig[channel]
    else:
        is_3d = False
        green, red = orig, None
        spot_image = orig
    spots = spot_image.asarray() > threshold
    labels, n_spots = img.label(spots)
    print "Found {} spots.".format(n_spots)

    # keep singlets (i.e. cells separated by at least min_distance)
    centers = np.array(img.center_of_mass(labels, labels, xrange(1, n_spots+1)))
    distances = spatial.distance_matrix(centers, centers)
    dx_rejects = reject(distances, min_distance)
    slices = img.find_objects(labels)
    for i in np.flatnonzero(dx_rejects):
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
    areas = img.sum(spots, labels, indices)
    data = {'x': centers[:,1] if len(centers) else [],
            'y': centers[:,0] if len(centers) else [],
            'area': areas,
            'green_intensity': green_scores}
    if is_3d:
        red_scores = img.sum(red.asarray(), labels, indices)
        data['red_intensity'] = red_scores
    scores = pd.DataFrame(data)
    return scores


def annotate_ld(livedead_fn, annotate_fn, channel, singlets, cell_radius, doublets=None):
    im = tf.TiffFile(livedead_fn)[channel].asarray()
    surface = surface_from_array(im)
    worklist = [(singlets, (1., 0., 0.))]
    if doublets is not None:
        worklist.append((doublets, (0., 1., 0.)))
    for cells, color in worklist:
        rows = len(cells)
        foo = lambda x: np.asarray(x).reshape(-1,1)
        cell_matrix = np.hstack([foo(cells['x']), foo(cells['y']), np.ones((rows, 1))])
        surface = annotate_cells(surface, cell_matrix, cell_radius, color, linewidth=2)
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
    parser.add_argument('--doublets',
                        help='Filename to save positions of doublets in live/dead image.')
    args = parser.parse_args()

    print 'Identifying singlets...'
    singlets = id_singlets(args.livedead, threshold=args.threshold, min_distance=args.distance,
                           min_size=args.cell_size[0], max_size=args.cell_size[1],
                           channel=args.channel)
    singlets.to_csv(args.output, index=False)
    print 'Found {} singlets.'.format(len(singlets))

    if args.doublets:
        print 'Identifying doublets...'
        doublets = id_singlets(args.livedead, threshold=args.threshold, min_distance=args.distance,
                               min_size=args.cell_size[0], max_size=args.cell_size[1],
                               channel=args.channel, reject=not_doublet_filter)
        doublets.to_csv(args.doublets, index=False)
        print 'Found {} doublets.'.format(len(doublets))

    if args.annotate:
        r = np.sqrt(args.cell_size[1]/np.pi)
        annotate_ld(args.livedead, args.annotate, args.channel, singlets, 2*r,
                    doublets=doublets if args.doublets else None)

if __name__ == '__main__':
    main()
