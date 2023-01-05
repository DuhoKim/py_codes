from calsagos import lagasu
from calsagos import utils
from calsagos import clumberi

# -- select cluster members
cluster_members = clumberi.clumberi(id_galaxy, ra_galaxy, dec_galaxy, redshift_galaxy, cluster_initial_redshift,
                                    ra_cluster, dec_cluster, range_cuts)

# -- defining output parameters from clumberi
id_member = cluster_members[0]
ra_member = cluster_members[1]
dec_member = cluster_members[2]
redshift_member = cluster_members[3]

# -- estimating the galaxy separation of galaxies in the cluster sample to be used as input in lagasu
knn_distance = utils.calc_knn_galaxy_distance(ra_member, dec_member, n_galaxies)

# -- determining the distance to the k-nearest neighbor of each galaxy in the cluster
knn_galaxy_distance = knn_distance[0]

typical_separation = utils.best_eps_dbscan(id_member, knn_galaxy_distance)

# -- Assign galaxies to each substructures
label_candidates = lagasu.lagasu(id_member, ra_member, dec_member, redshift_member, range_cuts, typical_separation,
                                 n_galaxies)

# -- defining output parameters from lagasu
id_candidates = label_candidates[0]
ra_candidates = label_candidates[1]
dec_candidates = label_candidates[2]
redshift_candidates = label_candidates[3]
label_zcut = label_candidates[4]
label_final = label_candidates[5]

# -- renaming the substructures identified by using lagasu in order to identify the principal halo and separate it from the substructures
id_final = utils.rename_substructures(ra_candidates, dec_candidates, redshift_candidates, label_final, ra_cluster,
                                      dec_cluster, r200_degree, flag)

# all done
