import numpy as np

def bound_linear(z, z_low, z_high, val_low, val_high):
    z0 = val_low + (val_high - val_low) * (z - z_low) / (z_high - z_low)

    vmin = min(val_low, val_high)
    vmax = max(val_low, val_high)

    return np.clip(z0, vmin, vmax)
