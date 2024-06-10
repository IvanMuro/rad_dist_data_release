import numpy as np

def calc_ci_rings(xplot, objct_lst, bkg_lst, l_lim=0.16, u_lim=0.84):
    """Function that calculates the ci for a group of given rings"""
    up_lim = []
    low_lim = []
    for annuli_loop, annuli in enumerate(xplot):
        
        objct = np.array(objct_lst)[:, annuli_loop]
        bkg = np.mean(bkg_lst, axis=0)[annuli_loop]

        ci_up, ci_low = calc_ci(objct, bkg, l_lim=l_lim, u_lim=u_lim)

        up_lim.append(ci_up)
        low_lim.append(ci_low)

    return up_lim, low_lim


def calc_ci(objct, bkg, l_lim=0.16, u_lim=0.84):
    """Function that calculates the confidence intervals of a measurement after substracting the bkg"""
    up_lim = np.quantile(objct - bkg, u_lim)
    low_lim = np.quantile(objct - bkg, l_lim)
    return up_lim, low_lim