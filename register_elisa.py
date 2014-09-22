import numpy as np
import numpy.linalg
import tifffile as tf
import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy as sp
import scipy.ndimage as img
import readroi
import cairo
import argparse
import pandas as pd
from PIL import Image
from shapely.geometry import Polygon, Point
from annotate import annotate_cells, surface_from_array, PIL_from_surface

def tx_matrix(from_points, to_points):
    # from_points and to_points are nx3 matrices
    # the third column should be 1s
    x_tx,_,_,_ = np.linalg.lstsq(from_points, to_points[:,0])
    y_tx,_,_,_ = np.linalg.lstsq(from_points, to_points[:,1])
    tx = np.vstack((x_tx.reshape((1,-1)), y_tx.reshape((1,-1)), [[0, 0, 1]]))
    return tx

def tx_matrix_from_rois(livedead_roi, elisa_roi, plot_filename=None):
    def prep_rois(fn):
        rois = readroi.read_roi_zip(fn)
        rois = np.vstack([roi[1] for roi in rois])
        rois = np.hstack((rois, np.ones((len(rois), 1))))
        return rois
    elisa_cx = prep_rois(elisa_roi)
    cell_cx = prep_rois(livedead_roi)
    tx = tx_matrix(cell_cx, elisa_cx)
    transformed = np.dot(tx, cell_cx.T).T
    if plot_filename:
        fig, ax = plt.subplots()
        ax.plot(elisa_cx[:,0], elisa_cx[:,1], 'o', transformed[:,0], transformed[:,1], 'r.')
        fig.savefig(plot_filename)
    return tx

def in_bounds(x, radius, shape):
    return ((x[:,0] - radius >= 0) &
            (x[:,0] + radius < shape[1]) &
            (x[:,1] - radius >= 0) &
            (x[:,1] + radius < shape[0]))

def measure_cells(data, cells, radius):
    def circle(rad):
        y, x = np.ogrid[-rad:rad, -rad:rad]
        mask = x*x + y*y <= rad*rad
        return np.ones((2*rad, 2*rad)) * mask

    stats = np.zeros((cells.shape[0], 8))
    mask = circle(radius) == 0
    for i, (x, y, z) in enumerate(list(cells)):
        region = data[y-radius:y+radius, x-radius:x+radius]
        region = np.ma.array(region, mask=mask)
        mean, sd, maxv, minv, integ = region.mean(), region.std(), region.max(), region.min(), region.sum()
        valid = region.compressed()
        lq = np.percentile(valid, 0.25)
        median = np.percentile(valid, 0.5)
        uq = np.percentile(valid, 0.75)
        stats[i,:] = mean, sd, maxv, minv, integ, lq, median, uq
    return stats

def filter_excluded(cells, roiset_fn):
    rois = readroi.read_roi_zip(roiset_fn)
    shapes = []
    for roi in rois:
        shapes.append(Polygon(list([list(i) for i in roi[1]])))
    keep = []
    for i, (x, y, z) in enumerate(cells):
        good = True
        p = Point((x,y))
        for shape in shapes:
            if shape.contains(p):
                good = False
                break
        if good:
            keep.append(i)
    return keep

def main():
    parser = argparse.ArgumentParser(description='Process single-cell ELISA images.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('livedead', metavar='livedead_results.txt',
                        help='Output of id_singlets.py.')
    parser.add_argument('livedead_roi', metavar='livedead-RoiSet.zip')
    parser.add_argument('elisa', metavar='stitched_elisa.tif', type=argparse.FileType('rb'),
                        help='Stitched ELISA image, post background subtraction.')
    parser.add_argument('elisa_roi', metavar='elisa-RoiSet.zip')
    parser.add_argument('--well-radius', default=50,
                        help='ELISA well radius, in pixels.')
    parser.add_argument('--tx-plot', '-t', default='transform.png',
                        help='Filename for transformation plot.')
    parser.add_argument('--exclude', '-e', metavar='exclude-RoiSet.zip',
                        help='ROI set of regions to exclude from measurement')
    parser.add_argument('--output', '-o', default='elisa_results.txt',
                        help='Filename to save information about the live/dead image.')
    args = parser.parse_args()

    df = pd.read_csv(args.livedead)
    cells = np.vstack([df['x'], df['y'], np.ones_like(df['x'])]).T
    tx = tx_matrix_from_rois(args.livedead_roi, args.elisa_roi, args.tx_plot)
    tx_cells = np.dot(tx, cells.T).T
    print "Loaded {} cells.".format(tx_cells.shape[0])

    data = tf.TiffFile(args.elisa)[0].asarray()
    row_idx = np.nonzero(in_bounds(tx_cells, args.well_radius, data.shape))[0]
    tx_cells = tx_cells[row_idx]
    print "Kept {} in-bounds cells.".format(tx_cells.shape[0])

    if args.exclude:
        not_excluded = filter_excluded(tx_cells, args.exclude)
        row_idx = row_idx[not_excluded]
        tx_cells = tx_cells[not_excluded]
    print "Kept {} cells not excluded.".format(tx_cells.shape[0])

    surface = surface_from_array(data)
    surface = annotate_cells(surface, tx_cells, args.well_radius, (1.0, 0.0, 0.0))
    pil_image = PIL_from_surface(surface)
    pil_image.save("annotated.jpg")

    stats = measure_cells(data, tx_cells, args.well_radius)
    out = np.hstack([row_idx.reshape((-1,1)), tx_cells[:,:2], stats])
    header = ['cells_row', 'x', 'y', 'mean', 'sd', 'max', 'min', 'integrated', 'q25', 'q50', 'q75']
    np.savetxt(args.output, out, delimiter=',', fmt=['%d'] * 3 + ['%f'] * 8, header=','.join(header), comments='')


if __name__ == '__main__':
    main()
