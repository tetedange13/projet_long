#!/usr/bin/env python3


import numpy as np


def calcul_distance(dico_pdb, dico_align):
    """
    Calcul distance between all residus pair by pair.
    Args: The dictionary containing the coordinate of the template's structure
    only CA and only the fragment which matched with the query
          The dictionary of one alignment
    Returns: A matrix of distance between all residus
    """
    size_query = len(dico_align['query'])
    # Create the empty matrix distance, which will contain distances
    # between all residus
    matrix_dist = np.empty( (size_query, size_query), dtype="float" )
    # Fill the matrix distance with NaN
    matrix_dist[:] = np.nan
    for i in range(0, size_query):
        # Only the half matrix is completed because it's symetric
        for j in range(i + 1, size_query):
            # A window of 5 residus is taken because it's not interesting to
            # calculate their distance, because they are too close.
            if j not in range(i, i + 5) :
                # If the two residus aren't a gap
                if (dico_align['query'][i] != "-" and
                    dico_align['query'][j] != "-" and
                    dico_align['template'][i] != "-"):
                    # calculate the distance between the two residus
                    dist = np.sqrt((dico_pdb[i]["x"] - dico_pdb[j]["x"])**2 +
                                   (dico_pdb[i]["y"] - dico_pdb[j]["y"])**2 +
                                   (dico_pdb[i]["z"] - dico_pdb[j]["z"])**2)
                    # fill the matrix with the calculated distance
                    matrix_dist[i, j] = dist
                else:
                    # -1 is added when one residu is a gap
                    matrix_dist[i, j] = -1
    return matrix_dist
    
    
def dist_to_proba(dist):1/(1+exp((x-8)/1.5)
    if dist > 15: # Cut-off, deducted from visualization, may be an prog arg ?
        return 0
    return 1/(1+exp((dist-8)/1.5)    
    

def extract_subarrays(proba_arr, m): # PI(m)
    A = proba_arr[m:n, 0:m] # There is also np.ix_()
    B = proba_arr[0:m, m:n]
    C = proba_arr[0:m, 0:m]
    
    return (A, B, C) # Faut-il retourner plutot la somme des matrices ??
    # return (A*B - C**2)/((A+C)(B+C)) 
    

def joint_proba_arr(proba_arr):



    
def calcul_energy(matrix_dist, dope_arr, index_aa, dico_align):
    """
    Calcul total energy of the query threaded on the template. A dope
    matrix is used to retrieve the energy by distance.
    Args: The matrix distance, which contains the distance between all residus
          The dope_arr is a 3D numpy array containing the associated energy by
    distance for two residus
          The index_aa is a dictionary containing the index of the dope_arr for
    each amino acid
          The dico_align is the dictionary of one alignment
    Returns: The total energy of the query threaded on the template
    """
    size_query = len(dico_align['query'])
    energy = 0
    for i in range(0, size_query):
        # Retrieve the first amino acid
        aa1 = dico_align['query'][i]
        # Only the half matrix is traveled
        for j in range(i + 1, size_query):
            # Retrieve the second amino acid
            aa2 = dico_align['query'][j]
            # Retrieve the distance calculated between the two residus
            dist = matrix_dist[i, j]
            # If the distance is -1, it's a gap so it's penalized
            if dist == -1:
                energy += cA.QUERY_GAP
            # No distance calculated because these two residus are too close
            # or they are too distant to have an impact
            elif not(np.isnan(dist) or dist > cA.CUTOFF):
                # Retrieve index of the amino acid in the dope matrix
                idx1 = index_aa[aa1]
                idx2 = index_aa[aa2]
                # Calculate the index in the list of distance in dope matrix
                # The distance begin to 0.25A and the step is 0.5
                idx_dist = int(round((dist - 0.25) / 0.5))
                # Retrieve the energy in the dope matrix between these two
                # residus with the distance calculated
                energy += dope_arr[idx_dist][idx1][idx2]
    return energy
