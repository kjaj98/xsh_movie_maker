import math

import numpy as np
from numba import njit, prange
# =====================================================
# 6a. Real-space SH_p p-orbital grid (Orbital-Visualiser style)
# =====================================================

def SH_p(grid_x, grid_y, grid_z, center, orbital_resolution):
    """
    Build p-like Cartesian components centered on `center` using Å-space
    coordinates scaled by `orbital_resolution`. Decreasing orbital_resolution
    narrows the orbital; increasing it widens the lobes.
    """
    gx = np.asarray(grid_x, dtype=float)
    gy = np.asarray(grid_y, dtype=float)
    gz = np.asarray(grid_z, dtype=float)
    scale = float(orbital_resolution)
    if scale <= 0:
        raise ValueError("orbital_resolution must be positive.")

    dx = (gx[:, None, None] - center[0]) / scale
    dy = (gy[None, :, None] - center[1]) / scale
    dz = (gz[None, None, :] - center[2]) / scale

    r = np.sqrt(dx*dx + dy*dy + dz*dz)
    exp_neg_r = np.exp(-r)

    px = dx * exp_neg_r
    py = dy * exp_neg_r
    pz = dz * exp_neg_r

    return np.array([px, py, pz])


# =====================================================
# 6b. CP2K-style p-vector calculation for one molecule
# =====================================================

def calc_pvecs_1mol(mol_crds, act_ats, cutoff=3.5):
    """
    Will calculate the pvecs for 1 molecule.

    Inputs:
        * mol_crds <array> => The coordinates of each atom in the molecule
        * act_ats  <array> => Which atom to calculate the pvecs for (local indices)
    Outputs:
        <array> The pvecs, shape (len(act_ats), 3)
    """
    nearest_neighbours = []
    for iat in act_ats:
        at_crd = mol_crds[iat]
        dists = np.linalg.norm(mol_crds - at_crd, axis=1)

        count = 0
        nn = [at_crd]
        for jat, d in enumerate(dists):
            if 0 < d < cutoff:
                nn.append(mol_crds[jat])
                count += 1
            if count == 2:
                break

        # fallback: if fewer than 3 points, pad with nearest neighbours
        if len(nn) < 3:
            order = np.argsort(dists)
            for jat in order:
                if jat == iat:
                    continue
                nn.append(mol_crds[jat])
                if len(nn) == 3:
                    break

        nearest_neighbours.append(nn)

    nearest_neighbours = np.array(nearest_neighbours)

    pvecs = []
    for a1, a2, a3 in nearest_neighbours:
        v1 = a2 - a1
        v2 = a3 - a1
        pvec = np.cross(v1, v2)
        n = np.linalg.norm(pvec)
        if n < 1e-8:
            pvec = np.array([0.0, 0.0, 1.0])
        else:
            pvec /= n
        pvecs.append(pvec)

    return np.array(pvecs)

# =====================================================
# 8. Orbital builders on 3D grid (SH_p + p-vectors)
# =====================================================

def build_global_grid(coords, spacing, margin):
    mins = coords.min(axis=0) - margin
    maxs = coords.max(axis=0) + margin
    xs = np.arange(mins[0], maxs[0] + spacing, spacing)
    ys = np.arange(mins[1], maxs[1] + spacing, spacing)
    zs = np.arange(mins[2], maxs[2] + spacing, spacing)
    return xs, ys, zs

@njit(cache=True, parallel=True)
def _build_orbital_wavefunction_site_core(gx, gy, gz,
                                          coords_site,
                                          pvecs_local,
                                          fragment_aom,
                                          orbital_resolution):
    """
    Numba-accelerated core routine that mirrors the original Python logic.
    """
    nx = gx.shape[0]
    ny = gy.shape[0]
    nz = gz.shape[0]
    n_atoms = coords_site.shape[0]

    psi = np.zeros((nx, ny, nz), dtype=np.complex128)

    # Pre-normalise p-vectors
    pnorm = np.zeros((n_atoms, 3), dtype=np.float64)
    for a in range(n_atoms):
        vx = pvecs_local[a, 0]
        vy = pvecs_local[a, 1]
        vz = pvecs_local[a, 2]
        n = math.sqrt(vx * vx + vy * vy + vz * vz)
        if n < 1e-8:
            pnorm[a, 0] = 0.0
            pnorm[a, 1] = 0.0
            pnorm[a, 2] = 1.0
        else:
            inv = 1.0 / n
            pnorm[a, 0] = vx * inv
            pnorm[a, 1] = vy * inv
            pnorm[a, 2] = vz * inv

    scale = float(orbital_resolution)

    for ix in prange(nx):
        x = gx[ix]
        for iy in range(ny):
            y = gy[iy]
            for iz in range(nz):
                z = gz[iz]

                acc = 0.0 + 0.0j

                for a in range(n_atoms):
                    a_k = fragment_aom[a]
                    if abs(a_k) < 1e-8:
                        continue

                    cx = coords_site[a, 0]
                    cy = coords_site[a, 1]
                    cz = coords_site[a, 2]

                    dx = (x - cx) / scale
                    dy = (y - cy) / scale
                    dz = (z - cz) / scale

                    r = math.sqrt(dx * dx + dy * dy + dz * dz)
                    e = math.exp(-r)

                    px = dx * e
                    py = dy * e
                    pz = dz * e

                    vx = pnorm[a, 0]
                    vy = pnorm[a, 1]
                    vz = pnorm[a, 2]

                    phi_i = vx * px + vy * py + vz * pz

                    acc += a_k * phi_i

                psi[ix, iy, iz] = acc

    return psi


def build_orbital_wavefunction_site(grid_x, grid_y, grid_z,
                                    coords,
                                    atom_indices,
                                    pvecs_local,
                                    fragment_aom,
                                    orbital_resolution):
    """
    Build φ_site(r) for ONE monomer using a Numba-accelerated core. The public
    API and resulting wavefunction match the previous NumPy implementation.
    """
    gx = np.asarray(grid_x, dtype=np.float64)
    gy = np.asarray(grid_y, dtype=np.float64)
    gz = np.asarray(grid_z, dtype=np.float64)

    coords = np.asarray(coords, dtype=np.float64)
    atom_indices = np.asarray(atom_indices, dtype=np.int64)
    coords_site = np.ascontiguousarray(coords[atom_indices], dtype=np.float64)

    pvecs = np.ascontiguousarray(pvecs_local, dtype=np.float64)
    if pvecs.shape[0] != coords_site.shape[0]:
        raise ValueError("pvecs_local and atom_indices must have the same length.")

    aom = np.ascontiguousarray(fragment_aom, dtype=np.complex128)

    psi = _build_orbital_wavefunction_site_core(
        gx,
        gy,
        gz,
        coords_site,
        pvecs,
        aom,
        float(orbital_resolution),
    )

    return psi

def normalize_orbital_on_grid(phi, gx, gy, gz, eps=1e-30):
    dx = gx[1] - gx[0]
    dy = gy[1] - gy[0]
    dz = gz[1] - gz[0]
    dV = dx * dy * dz

    norm2 = np.sum(np.abs(phi)**2) * dV
    if norm2 > eps:
        phi = phi / math.sqrt(norm2)
    return phi


def build_frame_densities(coords, sites,
                          donor_homo_fragment,
                          acc_homo_fragment,
                          acc_lumo_fragment,
                          h_prob, e_prob, x_prob,
                          spacing, margin,
                          orbital_resolution,
                          exciton_sites):
    """
    Population-only approximation in site basis:

      rho_h(r) ≈ ∑_i P_h(i) |φ_D^{(i)}(r)|²
      rho_e(r) ≈ ∑_j P_e(j) |φ_A^{(j)}(r)|²
      rho_x(r) ≈ ∑_k P_x(k) |φ_exc^{(k)}(r)|²
    """
    gx, gy, gz = build_global_grid(coords, spacing, margin)
    dx = gx[1] - gx[0]
    dy = gy[1] - gy[0]
    dz = gz[1] - gz[0]
    dV = dx * dy * dz

    # Hole density on donors (HOMO of donor fragment)
    rho_h = np.zeros((len(gx), len(gy), len(gz)))
    for i in range(len(h_prob)):
        P_h = h_prob[i]
        if P_h < 1e-10:
            continue
        label = f"D{i}"
        atom_inds = sites[label]["atoms"]
        mol_crds = coords[atom_inds]
        pvecs_local = calc_pvecs_1mol(mol_crds, np.arange(len(atom_inds)))
        phi_site = build_orbital_wavefunction_site(
            gx, gy, gz,
            coords,
            atom_indices=atom_inds,
            pvecs_local=pvecs_local,
            fragment_aom=donor_homo_fragment,
            orbital_resolution=orbital_resolution
        )
        phi_site = normalize_orbital_on_grid(phi_site, gx, gy, gz)
        rho_h += P_h * np.abs(phi_site)**2

    # Electron density on acceptors (LUMO of acceptor fragment)
    rho_e = np.zeros_like(rho_h)
    for j in range(len(e_prob)):
        P_e = e_prob[j]
        if P_e < 1e-10:
            continue
        label = f"A{j}"
        atom_inds = sites[label]["atoms"]
        mol_crds = coords[atom_inds]
        pvecs_local = calc_pvecs_1mol(mol_crds, np.arange(len(atom_inds)))
        phi_site = build_orbital_wavefunction_site(
            gx, gy, gz,
            coords,
            atom_indices=atom_inds,
            pvecs_local=pvecs_local,
            fragment_aom=acc_lumo_fragment,
            orbital_resolution=orbital_resolution
        )
        phi_site = normalize_orbital_on_grid(phi_site, gx, gy, gz)
        rho_e += P_e * np.abs(phi_site)**2

    # Exciton density (here: acceptor HOMO as exciton orbital)
    rho_x = np.zeros_like(rho_h)

    for j in range(len(x_prob)):
        P_x = x_prob[j]
        if P_x < 1e-10:
            continue
        label = f"A{j}"
        atom_inds = sites[label]["atoms"]
        mol_crds = coords[atom_inds]
        pvecs_local = calc_pvecs_1mol(mol_crds, np.arange(len(atom_inds)))

        phi_xh = build_orbital_wavefunction_site(
            gx, gy, gz,
            coords,
            atom_indices=atom_inds,
            pvecs_local=pvecs_local,
            fragment_aom=acc_homo_fragment,
            orbital_resolution=orbital_resolution
        )
        phi_xh = normalize_orbital_on_grid(phi_xh, gx, gy, gz)

        phi_xe = build_orbital_wavefunction_site(
            gx, gy, gz,
            coords,
            atom_indices=atom_inds,
            pvecs_local=pvecs_local,
            fragment_aom=acc_lumo_fragment,
            orbital_resolution=orbital_resolution
        )
        phi_xe = normalize_orbital_on_grid(phi_xe, gx, gy, gz)
        
        
        transition_density = phi_xe * np.conj(phi_xh)  # if fragments are complex
        rho_site = np.abs(transition_density)**2

        I_site = np.sum(rho_site) * dV
        if I_site > 0.0:
            rho_site /= I_site

        # weight by exciton population on that site
        rho_x += P_x * rho_site
        


    sum_e = np.sum(e_prob)
    sum_h = np.sum(h_prob)
    sum_x = np.sum(x_prob)

    I_e = np.sum(rho_e) * dV
    I_h = np.sum(rho_h) * dV
    I_x = np.sum(rho_x) * dV

    if I_e > 0 and sum_e > 0:
        rho_e *= (sum_e / I_e)
    if I_h > 0 and sum_h > 0:
        rho_h *= (sum_h / I_h)
    if I_x > 0 and sum_x > 0:
        rho_x *= (sum_x / I_x)

    return gx, gy, gz, rho_h, rho_e, rho_x
