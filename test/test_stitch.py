import numpy as np
import tifffile as tf

import stitch

def test_1x1(tmpdir):
    metadata = {
        'summary': {
            'InitialPositionList':
                [{'GridRowIndex': 0, 'GridColumnIndex': 0}],
            'Height': 100,
            'Width': 100,
        }
    }
    stack = np.ones(shape=(1, 100, 100), dtype=np.uint16) * 128
    canvas = stitch.stitch(metadata, stack)
    tf.imsave(str(tmpdir.join('1x1.tif')), canvas)
    # center is uniform
    assert np.all(canvas[11:-11, 11:-11] == 128)
    # vignetting is symmetric
    assert np.all(canvas[50, :] == canvas[50, ::-1])
    assert np.all(canvas[:, 50] == canvas[::-1, 50])

def test_2x2(tmpdir):
    metadata = {
        'summary': {
            'InitialPositionList':
                [{'GridRowIndex': r, 'GridColumnIndex': c} for (r, c) in zip([0, 0, 1, 1], [0, 1, 0, 1])],
            'Height': 100,
            'Width': 100,
        }
    }
    stack = np.ones(shape=(4, 100, 100), dtype=np.uint16) * 128
    canvas = stitch.stitch(metadata, stack)
    tf.imsave(str(tmpdir.join('2x2.tif')), canvas)
    delta = np.abs(canvas[11:-11, 11:-11].astype(np.int32) - 128)
    assert np.all(delta < 2)
