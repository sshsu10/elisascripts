import cairo
import numpy as np
from PIL import Image


def annotate_cells(cairo_surface, cell_xyz, radius, rgb, linewidth=5):
    cr = cairo.Context(cairo_surface)
    cr.set_source_rgb(*rgb)
    cr.set_line_width(linewidth)
    for x, y, z in list(cell_xyz):
        cr.arc(x, y, radius, 0, 2*np.pi)
        cr.stroke()
    return cairo_surface


def surface_from_array(a):
    minval = a[a.nonzero()].min()
    maxval = a.max()
    image_data = a.astype(float) - minval
    image_data *= 255.0/(maxval-minval)
    image_data = image_data.astype(np.uint8)
    cairo_image_data = np.dstack([image_data, image_data, image_data, np.ones_like(image_data)*255])
    surface = cairo.ImageSurface.create_for_data(
        cairo_image_data,
        cairo.FORMAT_ARGB32,
        *(reversed(image_data.shape))
    )
    return surface


def PIL_from_surface(surface):
    w = surface.get_width()
    h = surface.get_height()
    img = Image.frombuffer("RGBA", (w, h), surface.get_data(), "raw", "BGRA", 0, 1)
    return img
