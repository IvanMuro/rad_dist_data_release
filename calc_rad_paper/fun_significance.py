"""
This script has the functions that are called to calculate the significance and overdensity
at different radial distances
"""
import logging

import numpy as np

logger = logging.getLogger(__name__)


def get_peak_signif(
    ring_ovrdnst, corename, bands, oneSigma, sum_area, sigma_lim, sigma_ring_lim
):
    """ """
    significance = {}

    for name_clust in corename:
        logger.info(f"{name_clust}")
        significance.update({f"{name_clust}": {}})

        for energy_band in bands:
            logger.info(f"\t{energy_band}")
            significance[f"{name_clust}"].update({f"{energy_band}": []})

            for loop_model, up_lim in enumerate(oneSigma[name_clust][energy_band]):
                logger.info(f"\t{loop_model}")
                significance_lst = [
                    get_mask_significance(
                        area_lc, ring_ovrdnst, sigma_lim, sigma_ring_lim, up_lim
                    )
                    for area_lc in sum_area[name_clust][energy_band][loop_model]
                ]
                significance[f"{name_clust}"][f"{energy_band}"].append(significance_lst)
    return significance

def get_mask_significance(
    area_lc,
    ring_ovrdnst,
    sigma_lim,
    sigma_ring_lim,
    up_lim,
    all=False,
    max_ring=8,
):
    """ """
    if all == True:
        mask_ring_previous, mask_ring_after = get_mask_sigma_ring_all(
            ring_ovrdnst, area_lc, sigma_lim, up_lim, max_ring
        )
    else:
        mask_ring_previous = get_mask_sigma_ring(
            area_lc, ring_ovrdnst - 1, sigma_lim, up_lim, how="below"
        )
        mask_ring_after = get_mask_sigma_ring(
            area_lc, ring_ovrdnst + 1, sigma_lim, up_lim, how="below"
        )
    mask_ring_overdensity = get_mask_sigma_ring(
        area_lc, ring_ovrdnst, sigma_ring_lim, up_lim, how="above"
    )

    mask_relative = area_lc[ring_ovrdnst] > area_lc[0]
    for ring in np.arange(1, ring_ovrdnst, 1):
        mask_relative = np.logical_and(
            mask_relative, area_lc[ring_ovrdnst] > area_lc[ring]
        )
        ring = {ring}
        
    for ring in np.arange(ring_ovrdnst + 1, max_ring, 1):
        mask_relative = np.logical_and(
            mask_relative, area_lc[ring_ovrdnst] > area_lc[ring]
        )

    mask_list = [
        mask_ring_previous,
        mask_ring_after,
        mask_ring_overdensity,
        mask_relative,
    ]
    mask_significance = get_mask_joined(mask_list, mask_type="and")

    return mask_significance


def get_mask_sigma_ring_all(ring_ovrdnst, area_lc, sigma_lim, up_lim, max_ring):
    for ring in np.arange(0, ring_ovrdnst, 1):
        if ring == 0:
            mask_ring_previous = get_mask_sigma_ring(
                area_lc, ring, sigma_lim, up_lim, how="below"
            )
        else:
            mask = get_mask_sigma_ring(area_lc, ring, sigma_lim, up_lim, how="below")
            mask_ring_previous = np.logical_and(mask_ring_previous, mask)
    for ring in np.arange(ring_ovrdnst + 1, max_ring, 1):
        if ring == ring_ovrdnst + 1:
            mask_ring_after = get_mask_sigma_ring(
                area_lc, ring_ovrdnst + 1, sigma_lim, up_lim, how="below"
            )
        else:
            mask = get_mask_sigma_ring(area_lc, ring, sigma_lim, up_lim, how="below")
            mask_ring_after = np.logical_and(mask_ring_after, mask)
    return mask_ring_previous, mask_ring_after


def get_mask_joined(mask_list, mask_type="and"):
    """
    Function that join using "and" or "or"
    """
    mask_joined = mask_list[0]
    if mask_type == "and":
        for mask in mask_list[1:]:
            mask_joined = np.logical_and(mask_joined, mask)
    else:
        for mask in mask_list[1:]:
            mask_joined = np.logical_or(mask_joined, mask)
    return mask_joined


def get_mask_sigma_ring(area, ring, sigma, lim, how="below"):
    """
    Function that computes the mask for a ring above/below a given sigma limit

    Input:
    ------

    Output:
    -------

    """
    if how == "below":
        mask = area[ring] < sigma * lim[ring]
    else:
        mask = area[ring] > sigma * lim[ring]
    return mask