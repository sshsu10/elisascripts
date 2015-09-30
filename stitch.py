import argparse
import operator
import sys

from tifffile import tifffile as tf
import numpy as np

_description = "Assembles an xyp image stack into a large x,y image."

def stitch(metadata, ims):
    grid_index = [(p['GridRowIndex'], p['GridColumnIndex']) for p in
            metadata['summary']['InitialPositionList']]

    # starts at (0,0)
    N_row = max([i[0] for i in grid_index]) + 1
    N_col = max([i[1] for i in grid_index]) + 1

    h = metadata['summary']['Height']
    w = metadata['summary']['Width']

    f = lambda N, dim: int((0.9*N + 0.1) * dim)

    H = f(N_row, h)
    W = f(N_col, w)

    print("Making a {}x{} image.".format(W, H))
    canvas = np.zeros((H, W), np.uint16)

    h_margin = np.rint(0.1*h)
    w_margin = np.rint(0.1*w)
    fade = np.ones((h, w), np.double)
    h_ascending = np.linspace(0, 1, h_margin+1, endpoint=False)[1:].reshape(-1, 1)
    fade[:h_margin,:] *= h_ascending
    fade[h-h_margin:,:] *= 1-h_ascending
    w_ascending = np.linspace(0, 1, w_margin+1, endpoint=False)[1:].reshape(1, -1)
    fade[:,:w_margin] *= w_ascending
    fade[:,w-w_margin:] *= 1-w_ascending

    print('Stitching...')
    for i in xrange(len(ims)):
        if (i % 100) == 0:
            print i,
            sys.stdout.flush()
        r_n, c_n = grid_index[i]
        r_px = (N_row - 1 - r_n) * (0.9 * h)
        c_px = (N_col - 1 - c_n) * (0.9 * w)
        canvas[r_px:r_px+h, c_px:c_px+w] += np.rint(ims[i] * fade)
    print
    return canvas

def main():
    parser = argparse.ArgumentParser(description=_description)
    parser.add_argument('--output', '-o', metavar='stitched_elisa.tif',
                        default='stitched_elisa.tif')
    parser.add_argument('--saturationmap', '-s', metavar='saturated_elisa.tif',
                        default='saturated_elisa.tif')
    parser.add_argument(
            'original_image',
            help="Original MicroManager stack (used for metadata)")
    parser.add_argument(
            'normalized_image',
            help="Normalized image (used for image data)")
    args = parser.parse_args()
    orig_fn = args.original_image
    bgsub_fn = args.normalized_image
    print 'Reading metadata.'
    with open(orig_fn, 'rb') as f:
        metadata = tf.read_micromanager_metadata(f)
    print('Loading original image.')
    bgsub = tf.imread(bgsub_fn)
    canvas = stitch(metadata, bgsub)
    print('Saving image.')
    tf.imsave(args.output, canvas)
    del canvas, bgsub
    print('Computing saturation map.')
    orig = tf.imread(orig_fn)
    saturated = (orig == 4095)
    canvas = (stitch(metadata, saturated) > 0).astype(np.uint8)
    tf.imsave(args.saturationmap, canvas)

if __name__ == '__main__':
    main()
