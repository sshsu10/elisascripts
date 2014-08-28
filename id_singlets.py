import tifffile as tf
import numpy as np
import scipy as sp
import scipy.ndimage as img
import scipy.spatial as spatial
import argparse
import pandas as pd

def id_singlets(image_handle, threshold, min_size, max_size, min_distance):
    orig = tf.TiffFile(image_handle)
    red, green = orig
    spots = green.asarray() > threshold
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
    scores = pd.DataFrame({'x': centers[:,0],
                           'y': centers[:,1],
                           'area': areas,
                           'green_intensity': green_scores,
                           'red_intensity': red_scores})
    return scores

def main():
    parser = argparse.ArgumentParser(description='Process live/dead images.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('livedead', metavar='livedead.tif',
                        help='Stitched live/dead image, in TIFF format.')
    parser.add_argument('--threshold', '-t', required=True, type=int,
                        help='Threshold on [0, 4096) for live/dead image.')
    parser.add_argument('--cell-size', '-s', nargs=2, default=(10, 32), type=int,
                        metavar=('min', 'max'),
                        help='Recognize cells with areas on [min, max), in live/dead pixel units.')
    parser.add_argument('--distance', '-d', default=25, type=float,
                        help='Minimum acceptable separation (in pixels) between cells.')
    parser.add_argument('--output', '-o', default='livedead_results.txt',
                        help='Filename to save information about the live/dead image.')
    args = parser.parse_args()

    print 'Identifying singlets...'
    singlets = id_singlets(args.livedead, threshold=args.threshold, min_distance=args.distance,
                           min_size=args.cell_size[0], max_size=args.cell_size[1])
    singlets.to_csv(args.output, index=False)
    print 'Found {} singlets.'.format(len(singlets))

if __name__ == '__main__':
    main()
