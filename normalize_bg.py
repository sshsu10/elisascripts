import tifffile as tf
import numpy as np
import sys
import os

infile = sys.argv[1]

print 'Reading image stack'
t = tf.TiffFile(infile)
ar = tf.stack_pages(t.pages)
n = ar.shape[0]

if os.path.exists('background.tif'):
    print 'Reading background image'
    bg = tf.imread('background.tif')
else:
    print 'Computing background image'
    sorted_ar = ar.copy()
    sorted_ar.sort(0)
    bg = sorted_ar[int(0.05*n)]
    print 'Saving background image'
    tf.imsave('background.tif', bg)
    del sorted_ar

print 'Performing background normalization'
ar = ar.astype(np.double)
for i in range(n):
    ar[i] /= bg

print 'Converting to 16-bit TIFF'
max_normed = (4095.0 / bg.min()) - 1
ar -= 1
ar *= 65535
ar /= max_normed
ar = ar.round()
ar[ar < 0] = 0
ar = ar.astype(np.uint16)

print 'Writing normalized image'
with tf.TiffWriter('normalized.tif') as out:
    for i in range(n):
        if (i % 100) == 0:
            print i,
            sys.stdout.flush()
        out.save(ar[i])
print
