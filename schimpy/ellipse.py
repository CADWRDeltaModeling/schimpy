import math
import numpy as np

def ellipse(
    x,
    y,
    z,
    pt1_xy,
    pt2_xy,
    min_depth=None,
    max_depth=None,
    major_axis_len=1200.0,
):
    """Return depth floor based on distance from an oriented ellipse center.

    The minor axis follows pt1->pt2 and the major axis is perpendicular.
    Inside the ellipse, target depth decays quadratically from max_depth at
    the center to min_depth at the boundary. Outside the ellipse, target depth
    is clamped at min_depth. The return value is never shallower than z.
    """
    pt1_x, pt1_y = pt1_xy
    pt2_x, pt2_y = pt2_xy

    if min_depth is None:
        raise ValueError("min_depth must be provided")
    if max_depth is None:
        raise ValueError("max_depth must be provided")

    dx = float(pt2_x) - float(pt1_x)
    dy = float(pt2_y) - float(pt1_y)
    seg_len = math.hypot(dx, dy)
    if seg_len == 0.0:
        return z

    cx = 0.5 * (float(pt1_x) + float(pt2_x))
    cy = 0.5 * (float(pt1_y) + float(pt2_y))

    ux = dx / seg_len
    uy = dy / seg_len
    vx = -uy
    vy = ux

    semi_minor = max(seg_len / 2.0, 1.0)
    semi_major = max(float(major_axis_len) / 2.0, 1.0)

    x_arr = np.asarray(x)
    y_arr = np.asarray(y)
    z_arr = np.asarray(z)

    relx = x_arr - cx
    rely = y_arr - cy
    minor_coord = (relx * ux + rely * uy) / semi_minor
    major_coord = (relx * vx + rely * vy) / semi_major
    radial = minor_coord * minor_coord + major_coord * major_coord
    radial = np.minimum(1.0, radial)

    min_depth = float(min_depth)
    max_depth = float(max_depth)
    target_depth = max_depth - (max_depth - min_depth) * radial
    result = np.maximum(target_depth, z_arr)
    if np.ndim(result) == 0:
        return float(result)
    return result