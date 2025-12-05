'''
3x3 donor acceptor system example:
3 donor molecules, 3 acceptor molecules


y6 y6 y6 pdi pdi pdi

donor1 donor3 donor2 acceptor1 acceptor2 acceptor3


Psi; nadiab =      12; norbitals =       1
 i = 0, time =                0.000 (fs)
   1    1         0.0000000000        0.0000000000 # ct 11 (d1 a1)
   2    1         0.0000000000        0.0000000000 # ct 12 (d1 a2)
   3    1         0.0000000000        0.0000000000 # ct 13 (d1 a3)
   4    1         1.0000000000        0.0000000000 # ct 21 (d2 a1) # this is interfacial CT FYI
   5    1         0.0000000000        0.0000000000 # ct 22 (d2 a2)
   6    1         0.0000000000        0.0000000000 # ct 23 (d2 a3)
   7    1         0.0000000000        0.0000000000 # ct 31 (d3 a1) 
   8    1         0.0000000000        0.0000000000 # ct 32 (d3 a2)
   9    1         0.0000000000        0.0000000000 # ct 33 (d3 a3)
  10    1         0.0000000000        0.0000000000 # exciton on pdi 1
  11    1         0.0000000000        0.0000000000 # exciton on pdi 2
  12    1         0.0000000000        0.0000000000 # exciton on pdi 3

donor1_atomidx = 0..43
donor2_atomidx = 88..131
donor3_atomidx = 44..87
------------------------------- <- interface 
acceptor1_atomidx = 132..171
acceptor2_atomidx = 172..211
acceptor3_atomidx = 212..251

'''
#!/usr/bin/env python3
import argparse
import os
import shutil
import uuid
from datetime import datetime
from contextlib import contextmanager, redirect_stdout, redirect_stderr
from glob import glob
from typing import Dict, List, Optional
import time
from io_helpers import (
    read_xyz_trajectory,
    read_coeffs,
    read_aom_fragment,
    write_cube,
    write_single_frame_xyz,
)
from utils import (
    match_nuclear_and_electronic_times,
    build_site_atom_indices_from_ranges,
)
from create_real_space_psi import build_frame_densities
from vmd_api import (
    write_vmd_capture_script,
    write_vmd_server_script,
    capture_vmd_view,
    VMDRenderServer,
)
from make_animation import VideoSettings, make_video_from_tga_sequence
from xsh_analysis import (
    calculate_electron_site_probabilities,
    calculate_hole_site_probabilities,
    calculate_exciton_site_probabilities,
)


def build_default_config() -> Dict:
    base_data_dir = "./xsh_files_forKit"
    legacy_tga_dir = "./tga_images"

    config = {
        "paths": {
            "base_data_dir": base_data_dir,
            "active_xyz": os.path.join(base_data_dir, "run-pos-trimmed.xyz"),
            "whole_system_xyz": os.path.join(base_data_dir, "run-pos-1.xyz"),
            "coefficients": os.path.join(base_data_dir, "run-coeff-1.xyz"),
            "donor_homo_aom": os.path.join(base_data_dir, "HOMO_6T_AOM_COEFF.include"),
            "acceptor_homo_aom": os.path.join(base_data_dir, "HOMO_PDI_AOM_COEFF.include"),
            "acceptor_lumo_aom": os.path.join(base_data_dir, "LUMO_PDI_AOM_COEFF.include"),
            "run_dir": None,
            "wf_out_dir": None,
            "tga_out_dir": None,
            "logs_dir": None,
            "python_log": None,
            "vmd_log": None,
            "ffmpeg_log": None,
            "base_output_dir": "./runs",
            "legacy_view_state": os.path.join(legacy_tga_dir, "saved_view_matrices.txt"),
            "view_state_file": os.path.join("view_states", "saved_view_matrices.txt"),
        },
        "sites": {
            "nd": 3,
            "na": 3,
            "donor_ranges": [
                (0, 43),    # donor 1
                (44, 87),   # donor 3 (interfacial)
                (88, 131),  # donor 2
            ],
            "acceptor_ranges": [
                (132, 171),
                (172, 211),
                (212, 251),
            ],
            "natoms_expected": 252,
            "exciton_allowed_sites": "acceptor",
        },
        "grid": {
            "spacing": 0.05,   # Å - finer grid spacing for better resolution on cube files
            "margin": 2.0,     # Å - padding around the grid to capture orbital spill over 
            "orbital_resolution": 0.25,  # Å scale for SH_p basis (controls lobe size)
        },
        "workflow": {
            "region": "active",   # 'active' or 'full'
            "max_frames": None,     # None for all frames
            "keep_cubes": False,
            "keep_xyz": False,
            "align_to_first_frame": True,
            "generate_mp4": True,
        },
        "render": {
            "resize": [1280,720],
            "display": {
                "ambientocclusion": True,
                "dof": True,
                "depthcue": False,
                "projection": "Orthographic",
                "background": "white",
                "axes_location": "Off",
            },
            "atom_representation": {
                "remove_default": True,
                "representation": "Licorice 0.1 150.0 150.0",
                "selection": "all",
                "material": "Opaque",
                "color": "Name",
            },
            "atom_colors": {
                "H": "white",
                "C": "black",
                "O": "red",
                "N": "blue",
                "S": "yellow",
            },
            "isosurfaces": [
                {
                    "name": "hole",
                    "selection": "all",
                    "style": "Isosurface 2e-2 0 0 0 1 50",
                    "color": "ColorID 25",
                    "material": "AOChalky",
                },
                {
                    "name": "electron",
                    "selection": "all",
                    "style": "Isosurface 2e-2 1 0 0 1 50",
                    "color": "ColorID 30",
                    "material": "AOChalky",
                },
                {
                    "name": "exciton",
                    "selection": "all",
                    "style": "Isosurface 2e-2 2 0 0 1 50",
                    "color": "ColorID 20",
                    "material": "AOChalky",
                },
            ],
            "capture_scene": {
                "resize": [1280,720],
                "display": {
                    "ambientocclusion": None,
                    "dof": None,
                    "depthcue": None,
                    "projection": "Orthographic",
                    "background": "white",
                    "axes_location": None,
                },
                "atom_representation": {
                    "remove_default": True,
                    "representation": "Licorice 0.1 20.0 20.0",
                    "selection": "all",
                    "material": "Opaque",
                    "color": "Name",
                },
                "isosurfaces": [
                    {
                        "name": "hole",
                        "selection": "all",
                        "style": "Isosurface 5e-2 0 0 0 1 10",
                        "color": "ColorID 25",
                        "material": "AOChalky",
                    },
                    {
                        "name": "electron",
                        "selection": "all",
                        "style": "Isosurface 5e-2 1 0 0 1 10",
                        "color": "ColorID 30",
                        "material": "AOChalky",
                    },
                    {
                        "name": "exciton",
                        "selection": "all",
                        "style": "Isosurface 5e-2 2 0 0 1 10",
                        "color": "ColorID 20",
                        "material": "AOChalky",
                },
            ],
            "server_start_timeout": 100,
        },
        },
        "executables": {
            "vmd": "/scratch/kjoll/software/vmd/bin/vmd",
            "ffmpeg": "/scratch/kjoll/software/ffmpeg/ffmpeg-7.0.2-amd64-static/ffmpeg",
            "tachyon": "/scratch/kjoll/software/vmd/lib/vmd/tachyon_LINUXAMD64"
        },
        "video": {
            "output_file": "xsh_density_movie.mp4",
            "fps": 12,
            "codec": "libx264",
            "crf": 1,
            "preset": "medium",
            "pix_fmt": "yuv420p",
            "resolution": None,
            "extra_args": (),
        },
    }
    return config


CONFIG = build_default_config()


def _create_unique_run_dir(base_output_dir: str) -> str:
    base_output_dir = os.path.abspath(base_output_dir)
    while True:
        now = datetime.now()
        date_folder = now.strftime("%Y%m%d")
        timestamp = now.strftime("%H%M%S")
        unique = uuid.uuid4().hex[:8]
        run_name = f"run_{timestamp}_{unique}"
        candidate = os.path.join(base_output_dir, date_folder, run_name)
        try:
            os.makedirs(candidate, exist_ok=False)
            return candidate
        except FileExistsError:
            continue


def prepare_run_environment(config: Dict, requested_run_dir: Optional[str]) -> str:
    paths = config["paths"]
    base_output_dir = paths.get("base_output_dir", "./runs")
    if requested_run_dir:
        run_dir = os.path.abspath(requested_run_dir)
        os.makedirs(run_dir, exist_ok=True)
    else:
        run_dir = _create_unique_run_dir(base_output_dir)

    logs_dir = os.path.join(run_dir, "logs")
    wf_dir = os.path.join(run_dir, "wf_cubes")
    tga_dir = os.path.join(run_dir, "tga_images")

    paths.update({
        "run_dir": run_dir,
        "logs_dir": logs_dir,
        "wf_out_dir": wf_dir,
        "tga_out_dir": tga_dir,
        "python_log": os.path.join(logs_dir, "python.log"),
        "vmd_log": os.path.join(logs_dir, "vmd.log"),
        "ffmpeg_log": os.path.join(logs_dir, "ffmpeg.log"),
    })
    return run_dir


def ensure_view_state_destination(paths: Dict) -> None:
    view_state = paths.get("view_state_file")
    if not view_state:
        return
    view_state = os.path.abspath(view_state)
    view_dir = os.path.dirname(view_state)
    if view_dir:
        os.makedirs(view_dir, exist_ok=True)
    legacy = paths.get("legacy_view_state")
    if legacy:
        legacy = os.path.abspath(legacy)
        if os.path.exists(legacy) and not os.path.exists(view_state):
            shutil.copy2(legacy, view_state)

@contextmanager
def redirect_python_output(log_path: str):
    log_dir = os.path.dirname(os.path.abspath(log_path))
    os.makedirs(log_dir, exist_ok=True)
    with open(log_path, "w") as log_file, redirect_stdout(log_file), redirect_stderr(log_file):
        yield


def ensure_output_dirs(config: Dict) -> None:
    paths = config["paths"]
    run_dir = paths.get("run_dir")
    if not run_dir:
        raise RuntimeError("Run directory not initialised. Call prepare_run_environment first.")
    os.makedirs(run_dir, exist_ok=True)
    for key in ("wf_out_dir", "tga_out_dir", "logs_dir"):
        if not paths.get(key):
            raise RuntimeError(f"Missing path configuration for {key}.")
        os.makedirs(paths[key], exist_ok=True)
    ensure_view_state_destination(paths)


def require_view_state_file(view_state_file: str) -> None:
    if not os.path.exists(view_state_file):
        raise RuntimeError(
            f"View capture file '{view_state_file}' not found. "
            "Run `python main.py --mode capture` first to record a camera view."
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate densities, rerender VMD TGAs, or assemble videos."
    )
    parser.add_argument(
        "--mode",
        choices=("full", "rerender", "video", "capture"),
        default="full",
        help="Workflow to run: full pipeline, rerender, video-only, or capture camera view.",
    )
    parser.add_argument(
        "--max-frames",
        type=int,
        default=None,
        help="Optional override for number of frames to process.",
    )
    parser.add_argument(
        "--run-dir",
        default=None,
        help=(
            "Optional path to an existing run directory. "
            "If omitted, a new dated run directory is created automatically."
        ),
    )
    return parser.parse_args()


def create_vmd_scripts(config: Dict) -> tuple[str, str]:
    paths = config["paths"]
    render_cfg = config["render"]

    capture_script_path = os.path.join(paths["tga_out_dir"], "capture_view.tcl")
    write_vmd_capture_script(capture_script_path, render_cfg)
    server_script_path = os.path.join(paths["tga_out_dir"], "render_vmd_server.tcl")
    write_vmd_server_script(
        server_script_path,
        config["executables"]["tachyon"],
        render_cfg,
        view_state_file=paths["view_state_file"],
    )
    return capture_script_path, server_script_path


def get_effective_max_frames(config: Dict, override: Optional[int]) -> Optional[int]:
    return override if override is not None else config["workflow"]["max_frames"]


def run_full_workflow(config: Dict, max_frames_override: Optional[int]) -> None:
    paths = config["paths"]
    site_cfg = config["sites"]
    grid_cfg = config["grid"]
    workflow_cfg = config["workflow"]
    render_cfg = config["render"]
    exe_cfg = config["executables"]
    max_frames = get_effective_max_frames(config, max_frames_override)

    frames_geom = read_xyz_trajectory(paths["active_xyz"])
    frames_coeff = read_coeffs(paths["coefficients"])

    full_system_frames = None
    if workflow_cfg["region"] == "full":
        full_system_frames = read_xyz_trajectory(paths["whole_system_xyz"])

    print(f"Run directory: {paths['run_dir']}")
    print(f"Read {len(frames_geom)} geometry frames from {paths['active_xyz']}")
    print(f"Read {len(frames_coeff)} coeff frames from {paths['coefficients']}")
    if not frames_coeff:
        raise RuntimeError("No coefficient frames found.")

    natoms_traj = len(frames_geom[0]["atoms"])
    if site_cfg["natoms_expected"] is not None and natoms_traj != site_cfg["natoms_expected"]:
        print(f"WARNING: expected {site_cfg['natoms_expected']} atoms but found {natoms_traj}")

    n_states = len(frames_coeff[0]["coefficients"])
    if site_cfg["exciton_allowed_sites"] == "acceptor":
        expected_states = site_cfg["nd"] * site_cfg["na"] + site_cfg["na"]
    else:
        expected_states = site_cfg["nd"] * site_cfg["na"] + site_cfg["nd"]
    if n_states != expected_states:
        print(f"WARNING: n_states={n_states}, expected {expected_states}")

    frames = match_nuclear_and_electronic_times(frames_geom, frames_coeff)
    print("Matched nuclear and electronic times.")

    sites = build_site_atom_indices_from_ranges(
        site_cfg["donor_ranges"],
        site_cfg["acceptor_ranges"],
        natoms_traj,
    )
    print("\nSite->atom mapping:")
    for key in sorted(sites.keys()):
        val = sites[key]
        print(
            f"  {key}: type={val['type']}, natoms={len(val['atoms'])}, "
            f"indices=[{val['atoms'][0]}..{val['atoms'][-1]}]"
        )

    donor_natoms = site_cfg["donor_ranges"][0][1] - site_cfg["donor_ranges"][0][0] + 1
    acceptor_natoms = site_cfg["acceptor_ranges"][0][1] - site_cfg["acceptor_ranges"][0][0] + 1

    print("\nReading fragment AOM files...")
    donor_homo_fragment = read_aom_fragment(paths["donor_homo_aom"], donor_natoms)
    acc_homo_fragment = read_aom_fragment(paths["acceptor_homo_aom"], acceptor_natoms)
    acc_lumo_fragment = read_aom_fragment(paths["acceptor_lumo_aom"], acceptor_natoms)
    print("Fragment AOMs loaded and normalised.")

    _, server_script_path = create_vmd_scripts(config)
    require_view_state_file(paths["view_state_file"])
    vmd_log_path = paths.get("vmd_log")
    if vmd_log_path:
        print(f"Logging VMD server output to {vmd_log_path}")

    render_server: Optional[VMDRenderServer] = None
    try:
        try:
            render_server = VMDRenderServer(
                exe_cfg["vmd"],
                server_script_path,
                startup_timeout=render_cfg.get("server_start_timeout", 60),
                log_path=vmd_log_path,
            )
            print("Started persistent VMD render server.")
        except Exception as exc:
            raise RuntimeError("Failed to start persistent VMD render server.") from exc

        print("Entering loop over frames", flush=True)
        for iframe, frame in enumerate(frames):
            atoms = frame["atoms"]
            coords = frame["coords"]
            t = frame["time_fs"]
            coeffs = frame["coefficients"]

            print(f"\nFrame {iframe}: time = {t:.1f} fs", flush=True)

            e_prob = calculate_electron_site_probabilities(site_cfg["na"], site_cfg["nd"], coeffs)
            h_prob = calculate_hole_site_probabilities(site_cfg["na"], site_cfg["nd"], coeffs)
            x_prob = calculate_exciton_site_probabilities(
                site_cfg["na"],
                site_cfg["nd"],
                coeffs,
            )

            print(f"  Electron site probabilities (A0..A{site_cfg['na']-1}): {e_prob}")
            print(f"  Hole site probabilities     (D0..D{site_cfg['nd']-1}): {h_prob}")
            print(f"  Exciton site probabilities  (A0..A{site_cfg['na']-1}): {x_prob}")

            gx, gy, gz, rho_h, rho_e, rho_x = build_frame_densities(
                coords,
                sites,
                donor_homo_fragment,
                acc_homo_fragment,
                acc_lumo_fragment,
                h_prob,
                e_prob,
                x_prob,
                spacing=grid_cfg["spacing"],
                margin=grid_cfg["margin"],
                orbital_resolution=grid_cfg["orbital_resolution"],
                exciton_sites=site_cfg["exciton_allowed_sites"],
            )

            tag = f"f{iframe:08d}"
            hole_file = os.path.join(paths["wf_out_dir"], f"hole_{tag}.cube")
            elec_file = os.path.join(paths["wf_out_dir"], f"elec_{tag}.cube")
            exc_file = os.path.join(paths["wf_out_dir"], f"exc_{tag}.cube")

            write_cube(hole_file, gx, gy, gz, rho_h, atoms, coords, units="bohr")
            write_cube(elec_file, gx, gy, gz, rho_e, atoms, coords, units="bohr")
            write_cube(exc_file, gx, gy, gz, rho_x, atoms, coords, units="bohr")

            print(f"  Wrote cube files for tag {tag}", flush=True)

            xyz_frame_file = os.path.join(paths["tga_out_dir"], f"geom_{tag}.xyz")
            write_single_frame_xyz(xyz_frame_file, atoms, coords, t)
            render_xyz_file = xyz_frame_file

            if workflow_cfg["region"] == "full":
                if full_system_frames is None:
                    raise RuntimeError("Requested full-system rendering but no full trajectory was loaded.")
                full_frame = full_system_frames[iframe]
                render_xyz_file = os.path.join(paths["tga_out_dir"], f"full_system_{tag}.xyz")
                write_single_frame_xyz(render_xyz_file, full_frame["atoms"], full_frame["coords"], t)

            img_file = os.path.join(paths["tga_out_dir"], f"density_{tag}.tga")

            render_server.render_frame(
                hole_file,
                elec_file,
                exc_file,
                render_xyz_file,
                img_file,
            )

            if not workflow_cfg["keep_cubes"]:
                for cube_file in (hole_file, elec_file, exc_file):
                    if os.path.exists(cube_file):
                        os.remove(cube_file)
                print("  Deleted cube files to save space.")
            if not workflow_cfg["keep_xyz"]:
                if os.path.exists(xyz_frame_file):
                    os.remove(xyz_frame_file)
                if workflow_cfg["region"] == "full" and os.path.exists(render_xyz_file):
                    os.remove(render_xyz_file)

            if max_frames is not None and iframe + 1 >= max_frames:
                print(f"Reached max_frames limit ({max_frames}); stopping.")
                break
    finally:
        if render_server is not None:
            render_server.close()


def run_capture_mode(config: Dict) -> None:
    paths = config["paths"]
    site_cfg = config["sites"]
    grid_cfg = config["grid"]
    workflow_cfg = config["workflow"]
    exe_cfg = config["executables"]

    ensure_view_state_destination(paths)

    frames_geom = read_xyz_trajectory(paths["active_xyz"])
    frames_coeff = read_coeffs(paths["coefficients"])
    if not frames_geom or not frames_coeff:
        raise RuntimeError("Need at least one geometry and coefficient frame to capture a view.")

    full_system_frames = None
    if workflow_cfg["region"] == "full":
        full_system_frames = read_xyz_trajectory(paths["whole_system_xyz"])

    frames = match_nuclear_and_electronic_times(frames_geom, frames_coeff)
    sites = build_site_atom_indices_from_ranges(
        site_cfg["donor_ranges"],
        site_cfg["acceptor_ranges"],
        len(frames_geom[0]["atoms"]),
    )

    donor_natoms = site_cfg["donor_ranges"][0][1] - site_cfg["donor_ranges"][0][0] + 1
    acceptor_natoms = site_cfg["acceptor_ranges"][0][1] - site_cfg["acceptor_ranges"][0][0] + 1

    donor_homo_fragment = read_aom_fragment(paths["donor_homo_aom"], donor_natoms)
    acc_homo_fragment = read_aom_fragment(paths["acceptor_homo_aom"], acceptor_natoms)
    acc_lumo_fragment = read_aom_fragment(paths["acceptor_lumo_aom"], acceptor_natoms)

    capture_script_path, _ = create_vmd_scripts(config)
    frame = frames[0]
    atoms = frame["atoms"]
    coords = frame["coords"]
    t = frame["time_fs"]
    coeffs = frame["coefficients"]

    e_prob = calculate_electron_site_probabilities(site_cfg["na"], site_cfg["nd"], coeffs)
    h_prob = calculate_hole_site_probabilities(site_cfg["na"], site_cfg["nd"], coeffs)
    x_prob = calculate_exciton_site_probabilities(
        site_cfg["na"],
        site_cfg["nd"],
        coeffs,
    )

    gx, gy, gz, rho_h, rho_e, rho_x = build_frame_densities(
        coords,
        sites,
        donor_homo_fragment,
        acc_homo_fragment,
        acc_lumo_fragment,
        h_prob,
        e_prob,
        x_prob,
        spacing=grid_cfg["spacing"],
        margin=grid_cfg["margin"],
        orbital_resolution=grid_cfg["orbital_resolution"],
        exciton_sites=site_cfg["exciton_allowed_sites"],
    )

    tag = f"f{0:08d}"
    hole_file = os.path.join(paths["wf_out_dir"], f"hole_{tag}.cube")
    elec_file = os.path.join(paths["wf_out_dir"], f"elec_{tag}.cube")
    exc_file = os.path.join(paths["wf_out_dir"], f"exc_{tag}.cube")

    write_cube(hole_file, gx, gy, gz, rho_h, atoms, coords, units="bohr")
    write_cube(elec_file, gx, gy, gz, rho_e, atoms, coords, units="bohr")
    write_cube(exc_file, gx, gy, gz, rho_x, atoms, coords, units="bohr")

    xyz_frame_file = os.path.join(paths["tga_out_dir"], f"geom_{tag}.xyz")
    write_single_frame_xyz(xyz_frame_file, atoms, coords, t)
    render_xyz_file = xyz_frame_file

    if workflow_cfg["region"] == "full":
        if full_system_frames is None:
            raise RuntimeError("Requested full-system capture but no full trajectory was loaded.")
        full_frame = full_system_frames[0]
        render_xyz_file = os.path.join(paths["tga_out_dir"], f"full_system_{tag}.xyz")
        write_single_frame_xyz(render_xyz_file, full_frame["atoms"], full_frame["coords"], t)

    capture_vmd_view(
        exe_cfg["vmd"],
        capture_script_path,
        hole_file,
        elec_file,
        exc_file,
        render_xyz_file,
        paths["view_state_file"],
    )
    print(f"Saved VMD camera matrices to {paths['view_state_file']}.")


def collect_existing_frame_tags(wf_dir: str) -> List[str]:
    tags: List[str] = []
    for path in sorted(glob(os.path.join(wf_dir, "hole_f*.cube"))):
        base = os.path.basename(path)
        tag = base.replace("hole_", "").split(".")[0]
        tags.append(tag)
    return tags


def resolve_render_xyz(tga_dir: str, tag: str) -> str:
    full_path = os.path.join(tga_dir, f"full_system_{tag}.xyz")
    if os.path.exists(full_path):
        return full_path
    geom_path = os.path.join(tga_dir, f"geom_{tag}.xyz")
    if os.path.exists(geom_path):
        return geom_path
    raise FileNotFoundError(f"No XYZ file found for tag {tag}")


def rerender_existing_frames(config: Dict, max_frames_override: Optional[int]) -> None:
    paths = config["paths"]
    render_cfg = config["render"]
    exe_cfg = config["executables"]
    max_frames = get_effective_max_frames(config, max_frames_override)

    tags = collect_existing_frame_tags(paths["wf_out_dir"])
    if not tags:
        print("No existing cube files found; rerender mode requires saved cubes.")
        return

    print(f"Using run directory: {paths['run_dir']}")

    _, server_script_path = create_vmd_scripts(config)
    require_view_state_file(paths["view_state_file"])
    vmd_log_path = paths.get("vmd_log")
    if vmd_log_path:
        print(f"Logging VMD server output to {vmd_log_path}")

    render_server: Optional[VMDRenderServer] = None
    try:
        try:
            render_server = VMDRenderServer(
                exe_cfg["vmd"],
                server_script_path,
                startup_timeout=render_cfg.get("server_start_timeout", 60),
                log_path=vmd_log_path,
            )
            print("Started persistent VMD render server for rerender mode.")
        except Exception as exc:
            raise RuntimeError("Failed to start persistent VMD render server.") from exc

        print(f"Rerendering {len(tags)} frame(s) with updated VMD settings.")
        for idx, tag in enumerate(tags):
            hole_file = os.path.join(paths["wf_out_dir"], f"hole_{tag}.cube")
            elec_file = os.path.join(paths["wf_out_dir"], f"elec_{tag}.cube")
            exc_file = os.path.join(paths["wf_out_dir"], f"exc_{tag}.cube")
            render_xyz_file = resolve_render_xyz(paths["tga_out_dir"], tag)
            img_file = os.path.join(paths["tga_out_dir"], f"density_{tag}.tga")

            render_server.render_frame(
                hole_file,
                elec_file,
                exc_file,
                render_xyz_file,
                img_file,
            )

            if max_frames is not None and idx + 1 >= max_frames:
                print(f"Reached max_frames limit ({max_frames}); stopping rerender.")
                break
    finally:
        if render_server is not None:
            render_server.close()


def build_video_settings(config: Dict) -> VideoSettings:
    video_cfg = config["video"]
    paths = config["paths"]
    return VideoSettings(
        input_pattern=os.path.join(paths["tga_out_dir"], "density_f%08d.tga"),
        output_file=video_cfg["output_file"],
        fps=video_cfg["fps"],
        codec=video_cfg["codec"],
        crf=video_cfg["crf"],
        preset=video_cfg["preset"],
        pix_fmt=video_cfg["pix_fmt"],
        resolution=video_cfg["resolution"],
        extra_args=video_cfg["extra_args"],
    )


def generate_video_from_tgas(config: Dict) -> None:
    paths = config["paths"]
    video_settings = build_video_settings(config)
    tga_pattern = os.path.join(paths["tga_out_dir"], "density_f*.tga")
    if not glob(tga_pattern):
        print("No TGA files found; skipping video generation.")
        return
    print("\nGenerating video from TGA sequence...")
    make_video_from_tga_sequence(
        config["executables"]["ffmpeg"],
        video_settings,
        log_path=paths.get("ffmpeg_log"),
    )
    if paths.get("ffmpeg_log"):
        print(f"FFmpeg output logged to {paths['ffmpeg_log']}")
    print("Video generation complete.")


def main() -> None:
    args = parse_args()
    prepare_run_environment(CONFIG, args.run_dir)
    ensure_output_dirs(CONFIG)

    python_log = CONFIG["paths"]["python_log"]
    vmd_log = CONFIG["paths"].get("vmd_log")
    elapsed = 0.0
    success = False
    try:
        with redirect_python_output(python_log):
            start_time = time.time()
            if args.mode == "full":
                run_full_workflow(CONFIG, max_frames_override=args.max_frames)
                if CONFIG["workflow"]["generate_mp4"]:
                    generate_video_from_tgas(CONFIG)
            elif args.mode == "rerender":
                rerender_existing_frames(CONFIG, max_frames_override=args.max_frames)
            elif args.mode == "video":
                generate_video_from_tgas(CONFIG)
            elif args.mode == "capture":
                run_capture_mode(CONFIG)
            else:
                raise ValueError(f"Unknown mode {args.mode}")
            
            elapsed = time.time() - start_time
            print(f"Total execution time: {elapsed:.2f} seconds")
            print(
                f"Mode '{args.mode}' complete | N_frames limit: {args.max_frames} | "
                f"MP4 generation: {CONFIG['workflow']['generate_mp4']}"
            )
            success = True
    finally:
        status = "completed" if success else "failed"
        message = f"Mode '{args.mode}' {status}"
        if success:
            message += f" in {elapsed:.2f}s"
        message += f" | Run dir: {CONFIG['paths']['run_dir']}"
        message += f" | Python log: {python_log}"
        if vmd_log:
            message += f" | VMD log: {vmd_log}"
        ffmpeg_log = CONFIG["paths"].get("ffmpeg_log")
        if ffmpeg_log:
            message += f" | FFmpeg log: {ffmpeg_log}"
        print(message)


if __name__ == "__main__":
    main()
