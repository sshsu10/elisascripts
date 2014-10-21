import tifffile as tf
import sys
import numpy as np
import operator


def stitch(metadata, ims):
    grid_index = [(p['GridRowIndex'], p['GridColumnIndex']) for p in
            metadata['summary']['InitialPositionList']]

    # starts at (0,0)
    N_row = max([i[0] for i in grid_index]) + 1
    N_col = max([i[1] for i in grid_index]) + 1

    h = metadata['summary']['Height']
    w = metadata['summary']['Width']

    # deleteme
    #N_row, N_col = 10, 10
    #h, w = 256, 336

    H = int((0.9*(N_row-1) + 1) * h)
    W = int((0.9*(N_col-1) + 1) * w)

    print("Making a {}x{} image.".format(W, H))
    canvas = np.zeros((H, W), np.uint16)

    #ims = np.ones((100, h, w)) * 255 # deleteme
    #grid_index = zip(range(10)*10, reduce(operator.add, [[i]*10 for i in range(10)]))

    h_margin = int(0.1*h)
    w_margin = int(0.1*w)
    fade = np.ones((h, w), np.double)
    fade[:h_margin,:] *= np.arange(h_margin, dtype=np.double).reshape(-1,1)/h_margin
    fade[h-h_margin:,:] *= np.arange(h_margin-1, -1, -1, dtype=np.double).reshape(-1,1)/h_margin
    fade[:,:w_margin] *= np.arange(w_margin, dtype=np.double).reshape(1,-1)/w_margin
    fade[:,w-w_margin:] *= np.arange(w_margin-1, -1, -1, dtype=np.double).reshape(1,-1)/w_margin

    print('Stitching...')
    for i in xrange(len(ims)):
        if (i % 100) == 0:
            print i,
            sys.stdout.flush()
        r_n, c_n = grid_index[i]
        r_px = H - int((1 + 0.9*r_n)*h)
        c_px = W - int((1 + 0.9*c_n)*w)
        canvas[r_px:r_px+h, c_px:c_px+w] += ims[i] * fade
    print
    return canvas

def main():
    orig_fn = sys.argv[1]
    bgsub_fn = sys.argv[2]
    print 'Reading metadata.'
    with open(orig_fn, 'rb') as f:
        metadata = tf.read_micromanager_metadata(f)
    print('Loading original image.')
    bgsub = tf.imread(bgsub_fn)
    canvas = stitch(metadata, bgsub)
    print('Saving image.')
    tf.imsave('stitched_elisa.tif', canvas)
    del canvas, bgsub
    print('Computing saturation image.')
    orig = tf.imread(orig_fn)
    saturated = (orig == 4095)
    canvas = (stitch(metadata, saturated) > 0).astype(np.uint8)
    tf.imsave('saturated_elisa.tif', canvas)

if __name__ == '__main__':
    main()
