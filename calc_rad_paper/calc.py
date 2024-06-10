from astropy import units as u
import numpy as np

def get_annuli_rad(r, desire_z, cosmo):
    """ 
    Function that calculates the anulii in angular size for a given radius.
    """
    annulus_radius = (
        r * np.arange(0, 5.5, 0.5) / cosmo.angularDiameterDistance(desire_z) * u.rad
    )
    return annulus_radius