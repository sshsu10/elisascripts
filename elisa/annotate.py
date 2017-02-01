import cairo
import numpy as np
from PIL import Image

from typing import Tuple  # noqa:F401


def annotate_cells(cairo_surface,  # type: cairo.ImageSurface
                   cell_xyz,       # type: np.array  # (Nx3 floats)
                   radius,         # type: float
                   rgb,            # type: Tuple[int, int, int]
                   linewidth=5     # type: int
                   ):
    # type: (...) -> cairo.ImageSurface
    """Draws circles around cells.

    Accepts a Cairo ImageSurface and a Nx3 array of cell locations; gives the
    modified ImageSurface back.
    """
    cr = cairo.Context(cairo_surface)
    cr.set_source_rgb(*rgb)
    cr.set_line_width(linewidth)
    for x, y, z in list(cell_xyz):
        cr.arc(x, y, radius, 0, 2*np.pi)
        cr.stroke()
    return cairo_surface


def surface_from_array(a):
    # type: (np.array) -> cairo.ImageSurface
    """Produces a ImageSurface from a numpy array.

    Stretches array values to an 8-bit image and yields a Cairo ImageSurface
    for drawing.
    """
    minval = a[a.nonzero()].min()
    maxval = a.max()
    image_data = a.astype(float) - minval
    value_range = (maxval-minval) or 1
    image_data *= 255.0/(value_range)
    image_data = image_data.astype(np.uint8)
    cairo_image_data = np.dstack([image_data, image_data, image_data, np.ones_like(image_data)*255])
    surface = cairo.ImageSurface.create_for_data(
        cairo_image_data,
        cairo.FORMAT_ARGB32,
        *(reversed(image_data.shape))
    )
    return surface


def PIL_from_surface(surface):
    # type: (cairo.ImageSurface) -> Image.Image
    """Creates a PIL Image from a Cairo ImageSurface."""
    w = surface.get_width()
    h = surface.get_height()
    img = Image.frombuffer("RGBA", (w, h), surface.get_data(), "raw", "BGRA", 0, 1)
    return img
