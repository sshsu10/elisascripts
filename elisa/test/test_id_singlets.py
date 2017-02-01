import subprocess

import numpy as np
import pytest
from tifffile import tifffile as tf

from elisa import id_singlets


class TestIdSinglets(object):
    @pytest.fixture
    def singlet_image(self, tmpdir_factory):
        a = np.zeros((256, 256), dtype=np.uint8)
        size = 20
        # singlet
        a[20:20+size+1, 20:20+size+1] = 255
        # doublet
        a[100:100+size+1, 100:100+size+1] = 255
        offset = 100 + size + 10
        a[100:100+size+1, offset:offset+size+1] = 255

        # duplicate into second channel
        im = np.concatenate([a[np.newaxis, :], a[np.newaxis, :]])
        filename = tmpdir_factory.mktemp('data').join('img.tif')
        tf.imsave(str(filename), im)
        return filename

    def test_counts(self, singlet_image):
        def call(filter):
            return id_singlets.id_singlets(
                image_handle=str(singlet_image),
                channel=0,
                threshold=128,
                min_size=300,
                max_size=500,
                min_distance=50,
                reject=filter)
        singlets = call(id_singlets.not_singlet_filter)
        doublets = call(id_singlets.not_doublet_filter)
        assert len(singlets) == 1
        assert len(doublets) == 2

    def test_cli(self, tmpdir, singlet_image):
        with tmpdir.as_cwd():
            subprocess.check_call([
                "python",
                "-m",
                "elisa.id_singlets",
                "--threshold",
                "128",
                "--cell-size",
                "300",
                "500",
                "--distance",
                "50",
                "--doublets",
                "doublets.txt",
                "--annotate",
                str(tmpdir.join("annotated.jpg")),
                str(singlet_image)])
