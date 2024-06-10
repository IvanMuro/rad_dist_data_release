"""
"""
import logging

from astropy import units as u
import numpy as np

import mockUniv.calc_rad_paper.fun_features as fun_features 
import mockUniv.calc_rad_paper.ci as ci
import mockUniv.calc_rad_paper.fun_significance as fun_significance

logger = logging.getLogger(__name__)


def get_sum_bkg_subs_and_ci(
    xplot,
    lc_type,
    lc_type_f,
    corename,
    bands,
    sum_area_model_dict,
    sum_area_field_model_suffle_dict,
    mod=False,
    l_lim=0.16,
    u_lim=0.84            
):
    """
    Input:
    ------
    + xplot,
    + lc_type,
    + lc_type_f,
    + corename,
    + bands,
    + sum_area_model_dict,
    + sum_area_field_model_suffle_dict

    Output:
    -------
    """
    sum_area_model_dict_bkgSubs = {}
    sum_area_model_dict_bkgSubs_indiv = {}
    oneSigmaUp_sum_area_model_dict_bkgSubs = {}
    oneSigmaLow_sum_area_model_dict_bkgSubs = {}
    
    for name_clust in corename:
        sum_area_model_dict_bkgSubs.update({f"{name_clust}": {}})
        sum_area_model_dict_bkgSubs_indiv.update({f"{name_clust}": {}})
        oneSigmaUp_sum_area_model_dict_bkgSubs.update({f"{name_clust}": {}})
        oneSigmaLow_sum_area_model_dict_bkgSubs.update({f"{name_clust}": {}})

        for energy_band in bands:
            sum_area_model_dict_bkgSubs[f"{name_clust}"].update({f"{energy_band}": []})
            sum_area_model_dict_bkgSubs_indiv[f"{name_clust}"].update(
                {f"{energy_band}": []}
            )
            oneSigmaUp_sum_area_model_dict_bkgSubs[f"{name_clust}"].update(
                {f"{energy_band}": []}
            )
            oneSigmaLow_sum_area_model_dict_bkgSubs[f"{name_clust}"].update(
                {f"{energy_band}": []}
            )

            for loop_model, sum_area_loop in enumerate(
                sum_area_model_dict[name_clust][energy_band][lc_type][:]
            ):
                if mod == True:
                    if loop_model < 3:
                        loop_model_f = 0
                    else:
                        loop_model_f = 1
                else:
                    loop_model_f = loop_model
                up_lim, low_lim = ci.calc_ci_rings(
                    xplot,
                    sum_area_loop,
                    sum_area_field_model_suffle_dict[name_clust][energy_band][
                        lc_type_f
                    ][loop_model_f],
                    l_lim=l_lim,
                    u_lim=u_lim,
                )

                area_bkg_subs_mean, area_bkg_subs_indiv = substract_bkg(
                    sum_area_loop,
                    sum_area_field_model_suffle_dict[name_clust][energy_band][
                        lc_type_f
                    ][loop_model_f],
                )

                sum_area_model_dict_bkgSubs[f"{name_clust}"][f"{energy_band}"].append(
                    area_bkg_subs_mean
                )
                sum_area_model_dict_bkgSubs_indiv[f"{name_clust}"][
                    f"{energy_band}"
                ].append(area_bkg_subs_indiv)

                oneSigmaUp_sum_area_model_dict_bkgSubs[f"{name_clust}"][
                    f"{energy_band}"
                ].append(up_lim)
                oneSigmaLow_sum_area_model_dict_bkgSubs[f"{name_clust}"][
                    f"{energy_band}"
                ].append(low_lim)

    return (
        sum_area_model_dict_bkgSubs,
        sum_area_model_dict_bkgSubs_indiv,
        oneSigmaUp_sum_area_model_dict_bkgSubs,
        oneSigmaLow_sum_area_model_dict_bkgSubs,
    )


def substract_bkg(objcts, bkg):
    """ """
    mean_bkg = np.mean(bkg, axis=0)
    mean = np.mean(objcts, axis=0) - mean_bkg
    indiv = [obj - mean_bkg for obj in objcts]
    return mean, indiv


def sum_area_ids(
    corename, bands, lc_dict_flux, annulus_radius, models, models_f=None
):
    """
    """
    sum_area_model_dict = {}
    for name_clust in corename:
        print(f'\t{name_clust}')
        sum_area_model_dict.update( { f"{name_clust}":{} } )
        
        for energy_band in bands:
            print(f'\t\t{energy_band}')
            sum_area_model_dict[name_clust].update( { f"{energy_band}":{} } )
            
            for lc_list in lc_dict_flux:
                print(f'\t\t\t{lc_list}')
                if np.logical_and(lc_list=='f_rand', models_f!=None):                    
                    models_aux = models_f

                else:
                    models_aux = models
                    
                sum_area_model = get_area_summed_sample_lcs(lc_dict_flux[lc_list][:], 
                                                                           annulus_radius, 
                                                                           energy_band,
                                                                           models_aux, 
                                                                           name_clust,
                                                                          )
                
                sum_area_model_dict[name_clust][energy_band].update({f"{lc_list}":sum_area_model})     
    return sum_area_model_dict

def get_sum_area_field_suffle_ids(num_obs, corename, bands, sum_area_model_dict, models):
    """
    This function returns the summed area of a random subsample with repetitions of the 
    100 light cones of the field. This value is background that is substracted to the radial 
    of the cluster.
    """
    sum_area_field_model_suffle_dict = {}     
    for name_clust in corename:
        print(f'\t{name_clust}')
        sum_area_field_model_suffle_dict.update( { f"{name_clust}":{} } )
        
        for energy_band in bands:
            #print(f'\t\t{energy_band}')
            sum_area_field_model_suffle_dict[name_clust].update( { f"{energy_band}":{} } )
            
            for lc_type in ['f_rand']: 
                #print(f'\t\t\t{lc_list}')
                
                sum_area_field_model_rand_obs_suffle = []
                for loop_model in range(len(models)):
                    sum_area_field_model_rand_obs_suffle_aux = []
                    
                    #I fix a random seed to obtain exactly the same results as in the paper
                    random_state = 12823
                    rng = np.random.RandomState(random_state)
                    rand_lst = rng.choice(np.arange(0, num_obs, 1), size=100)
                    
                    for idx in rand_lst:
                        sum_area_field_model_rand_obs_suffle_aux.append(
                            sum_area_model_dict[name_clust][energy_band][lc_type][loop_model][idx]
                        )     
                    sum_area_field_model_rand_obs_suffle.append(sum_area_field_model_rand_obs_suffle_aux)
                sum_area_field_model_suffle_dict[name_clust][energy_band].update(
                        {f"{lc_type}":sum_area_field_model_rand_obs_suffle}
                    )
    return sum_area_field_model_suffle_dict

######
def get_area_anulii_models(models, flux_list_model, interpol_flux_ring):
    """
    Function that returns the area for each ring and input seeding model

    Input:
    ------
    - models
    - flux_list_model
    - interpol_flux_ring

    Output:
    -------
    - area_anulii_list_model
    """
    area_anulii_list_model = []
    for loop_model, model in enumerate(models):
        area_anulii_list_model.append(
            get_area_anulii(flux_list_model[loop_model], interpol_flux_ring)
        )
    return area_anulii_list_model


def get_area_anulii(flux_list, interpol_flux_ring):
    """
    Function that converts fluxes into probabilities of detection
    Input:
    -----
    - flux_list (list): list of fluxes for diff. rings
    - interpol_flux_ring (list): list of interpol1D(scipy object) that converts flux
    into probabilities. See https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html

    Output:
    -------
    - area_anulii_list (list): list of porbabilities of detection associted to the input fluxes

    """
    area_anulii_list = []
    for flux_annuli in flux_list:
        area_anulii_list_aux = []
        
        for loop_ring, interpol in enumerate(interpol_flux_ring[:-1]):
            if loop_ring < len(flux_annuli):
                area_anulii_list_aux.append(interpol(flux_annuli[loop_ring]))
        area_anulii_list.append(area_anulii_list_aux)
    return area_anulii_list


def _area_model(area_anulii_list_model, models):
    """
    Function that returns the total area within a ring for the different seeding models.
    It loops calls to get_sum_area(...) for the diferent seeding models.

    Input:
    ------
    - area_anulii_list_model
    - models

    Output:
    -------

    """
    print(models)
    return [
        get_sum_area(area_anulii_list_model[loop_model])
        for loop_model, model in enumerate(models)
    ]

def get_sum_area_model(area_anulii_list_model, models):
    """
    Function that returns the total area within a ring for the different seeding models.
    It loops calls to get_sum_area(...) for the diferent seeding models.

    Input:
    ------
    - area_anulii_list_model
    - models

    Output:
    -------

    """
    print(models)
    return [
        get_sum_area(area_anulii_list_model[loop_model])
        for loop_model, model in enumerate(models)
    ]
    
def get_sum_area(area_anulii_list_model):
    """
    Function that get the total area within a ring calling do_sum_area(...) which actually do the sum.

    Input:
    ------
    - area_anulii_list_model

    Output:
    -------
    - List ...
    """
    return [ do_sum_area(area_anulii) for area_anulii in area_anulii_list_model ]


def do_sum_area(area_anulii):
    """
    Function that do the sum
    """
    return [np.sum(a, axis=0) for a in area_anulii]


def _get_area_sample_lcs(
    lc_sample,
    annulus_radius,
    band,
    models,
    corename,
    get_all_feat=False,
):
    """
    Function that returns the total area and the flux and redshift distributions. It calls
    _get_features_sample_lcs(...) to get the necesary data (number og AGN, flux list...), then it
    calls get_area_anulii_models(...) to obtain the area with the input flux list calculated in
    the previous step.

    Input:
    -----
    - lc_sample,
    - annulus_radius,
    - band,
    - models,
    - corename,
    - interpol_flux_ring,
    - get_all_feat=False,

    Output:
    -------
    - area_anulii_list_model,
    - histo_flux_lc_sample_model,
    - histo_z_lc_sample_model,
    """
    (
        num_agn_annuli_lc_sample_model,
        flux_annuli_lc_sample_model,
        z_annuli_lc_sample_model,
        histo_flux_lc_sample_model,
        histo_z_lc_sample_model,
        area_anulii_list_model,
    ) = fun_features._get_features_sample_lcs(
        lc_sample, annulus_radius, band, models, corename 
    )

    if get_all_feat == False:
        return area_anulii_list_model
    else:
        return (
            area_anulii_list_model,
            histo_flux_lc_sample_model,
            histo_z_lc_sample_model,
        )


def get_area_summed_sample_lcs(
    lc_sample,
    annulus_radius,
    band,
    models,
    corename,
):
    """
    Function that calculate the total area (probability of detecting an AGN) in each annulus/ring.
    First it calls _get_area_sample_lcs(...) to get the area of each ring. Then it calls
    get_sum_area_model(...) to sum each of the components within the rings.

    Input:
    ------
    - lc_sample,
    - annulus_radius,
    - band,
    - models,
    - corename,
    - interpol_flux_ring,

    Output:
    -------
    - sum_area_model
    """

    area_anulii_list_model = _get_area_sample_lcs(
        lc_sample, annulus_radius, band, models, corename
    )
    sum_area_model = get_sum_area_model(area_anulii_list_model, models)
    return sum_area_model