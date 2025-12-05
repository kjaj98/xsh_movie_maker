import numpy as np

# =====================================================
# 3. Match nuclear / electronic times
# =====================================================

def match_nuclear_and_electronic_times(frames_geom, frames_coeff):
    times_elec = np.array([f['time_fs'] for f in frames_coeff])

    out_frames = []
    for fg in frames_geom:
        t = fg['time_fs']
        idx = np.where(times_elec == t)[0]
        if len(idx) == 0:
            raise ValueError(f"No matching coeffs for nuclear time {t} fs")
        if len(idx) > 1:
            raise ValueError(f"More than one coeff block for time {t} fs")
        i = idx[0]
        merged = dict(fg)
        merged['coefficients'] = frames_coeff[i]['coefficients']
        out_frames.append(merged)
    return out_frames


# =====================================================
# 4. Site -> atom mapping from ranges
# =====================================================

def build_site_atom_indices_from_ranges(D_ranges, A_ranges, natoms_total):
    sites = {}

    for i, (start, end) in enumerate(D_ranges):
        if start < 0 or end < start or end >= natoms_total:
            raise ValueError(f"D{i} range ({start},{end}) out of bounds for natoms_total={natoms_total}")
        atoms = list(range(start, end + 1))
        sites[f"D{i}"] = {"type": "donor", "atoms": atoms}

    for j, (start, end) in enumerate(A_ranges):
        if start < 0 or end < start or end >= natoms_total:
            raise ValueError(f"A{j} range ({start},{end}) out of bounds for natoms_total={natoms_total}")
        atoms = list(range(start, end + 1))
        sites[f"A{j}"] = {"type": "acceptor", "atoms": atoms}

    return sites

