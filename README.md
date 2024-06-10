# Scrutinising evidence for the triggering of Active Galactic Nuclei in the outskirts of massive galaxy clusters at z~1

This repository provides the code that reproduces the results presented in the paper "Scrutinising evidence for the triggering of Active Galactic Nuclei in the outskirts of massive galaxy clusters at z~1" by Munoz-Rodriguez et al. (2024, MNRAS accepted, XXX DOI XXX). The notebook rad_dist_data_release.ipynb demonstrates how to reproduce the plots on the radial distribution of AGN in clusters of galaxies presented in the above paper. Running this code requires the simulated light cones of massive clusters stored in the hdf5 files that are accesible via are accessible via zenodo https://zenodo.org/records/11446317. These files should be copied in the "./input" directory of the repository. 

The data in Zenodo includes 400 light cones (under the subdirectory with name lc_rad_dst), the reconstucted radial distribution of number counts for the chandra observations of clusters PLCKG266.6-27.3 and SPT-CLJ2146-4633 (under the subdirectory with name rad_dst_obs) and the list of UniverseMachine ids that are considered as infall galaxies as described in Section 5.2 of the paper (id_infall_pop_smMask_mDMmask.csv). 

The 400 light cones are organised in a set of subdirectories attending to the target cluster in the simulation (with UniverseMachine ids 7830644447 and 7793510527) if they point to the cluster (c) or the field (f), see Section 3.2 of the paper for further details. The light cones are used to estimate the effect of cosmic variance on the observed radial distribution of AGN in clusters and test claims of preferential activation of AGN in clusters in the sample of clusters of Koulouridis & Bartalucci (2019).

The README.pdf document provides information on the columns included in the light cone hdf5 tables.

Contribution to the code: Ivan Munoz Rodriguez, Antonis Georgakakis and Angel Ruiz Camunas.
