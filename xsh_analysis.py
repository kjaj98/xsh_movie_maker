import numpy as np
# =====================================================
# 5. XSH probabilities (ND×NA CT + exciton)
# =====================================================
def calculate_hole_site_probabilities(NA, ND, complex_coeffs):
    pops = np.abs(complex_coeffs)**2
    num_CT = ND * NA
    CT_pops = pops[:num_CT]

    hole_prob = np.zeros(ND)
    for i in range(ND):
        start = i * NA
        end   = start + NA
        hole_prob[i] = np.sum(CT_pops[start:end])
    return hole_prob

def calculate_electron_site_probabilities(NA, ND, complex_coeffs):
    pops = np.abs(complex_coeffs)**2
    num_CT = ND * NA
    CT_pops = pops[:num_CT]

    elec_prob = np.zeros(NA)
    for j in range(NA):
        elec_prob[j] = np.sum(CT_pops[j::NA])
    return elec_prob

def calculate_exciton_site_probabilities(NA, ND, complex_coeffs):
    pops = np.abs(complex_coeffs)**2

    num_CT = ND * NA
    exc_pops = pops[num_CT:]

    nsites = NA # acceptor sites only

    exc_prob = np.zeros(nsites)
    exc_prob[:len(exc_pops)] = exc_pops[:nsites]
    return exc_prob