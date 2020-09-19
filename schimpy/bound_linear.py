import numpy as np
def bound_linear(z,z_low,z_high,val_low,val_high): 
    z0 = val_low + (val_high - val_low) * (z - z_low) / (z_high - z_low)
    z0 = np.maximum(z0,val_low)
    z0 = np.minimum(z0,val_high)
    return z0