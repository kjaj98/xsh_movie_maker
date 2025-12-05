# XSH Movie Maker

## What it is
A small Python/VMD pipeline for visualising exciton/charge-transfer dynamics in a 3x3 donor/acceptor molecular stack. It builds real-space hole, electron, and exciton densities from site-population coefficients and turns them into rendered frames and an MP4 movie.

## What it does
- Reads a trimmed geometry trajectory (`xsh_files_forKit/run-pos-trimmed.xyz`) and matching electronic coefficients (`xsh_files_forKit/run-coeff-1.xyz`).
- Uses fragment AOM coefficients for donor HOMO and acceptor HOMO/LUMO to build approximate p-like orbitals on a 3D grid.
- Computes time-resolved electron, hole, and exciton site populations, converts them to real-space densities, and saves each frame as cube files.
- Renders each frame in VMD (via a persistent render server + Tachyon) to TGA images and optionally stitches them into an MP4 with ffmpeg.

## How it works (pipeline)
1. **Run setup**: `main.py` builds a dated run directory under `runs/DATE/run_<timestamp>_<id>` with `wf_cubes`, `tga_images`, and `logs` subfolders.
2. **Load data**: geometry frames and electronic coefficients are time-matched (`utils.match_nuclear_and_electronic_times`). Site-to-atom mappings are defined in `CONFIG['sites']` (3 donors, 3 acceptors) with expected atom ranges for the active region.
3. **Site probabilities**: `xsh_analysis.py` turns complex coefficients into per-site hole, electron, and exciton populations.
4. **Real-space densities**: `create_real_space_psi.py` builds a global grid, constructs fragment p-orbitals using AOM weights and local p-vectors, normalises them, and combines them with site populations to produce hole/electron/exciton densities. Cube files are written for each frame (`io_helpers.write_cube`).
5. **Rendering**: `vmd_api.VMDRenderServer` keeps VMD running in text mode. For each frame it loads the cube volumes plus the geometry (`geom_*.xyz` or `full_system_*.xyz`), applies a saved camera view, and renders an isosurface scene to `density_f########.tga` using Tachyon.
6. **Video assembly**: if enabled, `make_animation.py` calls ffmpeg on the TGA sequence to produce `xsh_density_movie.mp4` (configurable codec/FPS/CRF/resolution).

## Repository layout
- `main.py`: orchestrates the workflow and holds the default configuration.
- `io_helpers.py`: read/write helpers for XYZ, coefficient blocks, AOM fragments, and cube files.
- `create_real_space_psi.py`: builds orbital wavefunctions and densities on a grid (Numba-accelerated core).
- `xsh_analysis.py`: converts complex amplitudes to per-site populations.
- `vmd_api.py`: emits VMD TCL scripts, captures/loads camera views, and manages the persistent render server.
- `make_animation.py`: thin ffmpeg wrapper for turning TGAs into MP4s.
- `xsh_files_forKit/`: example input trajectory, coefficients, and fragment AOM files.
- `view_states/saved_view_matrices.txt`: saved VMD camera matrices (capture once, then reused).
- `runs/`: sample outputs (cube/TGA/log directories) from previous runs.
- `xsh_density_movie.mp4`: example movie produced from the sample data.

## Requirements
- Python 3 with `numpy` and `numba` installed.
- VMD with Tachyon renderer available on PATH (or update `CONFIG['executables']['vmd']`/`['tachyon']`).
- ffmpeg for MP4 assembly.
- The provided XYZ/coefficient/AOM input files (or your own, with matching site definitions).

## Usage
1. **Adjust paths/executables**: edit `build_default_config()` in `main.py` to point to your VMD, ffmpeg, Tachyon binaries, and input data. Tweak grid spacing/margins or site atom ranges if your system differs.
2. **Capture a camera view (once per system)**: run
   ```bash
   python main.py --mode capture
   ```
   This renders the first frame, lets you orient the scene in the VMD GUI, and saves matrices to `view_states/saved_view_matrices.txt` (copied from legacy path if present).
3. **Run the full pipeline** (densities + TGA renders + MP4):
   ```bash
   python main.py --mode full --max-frames 50
   ```
   - `--max-frames` is optional; omit to process all frames.
   - Outputs go to a new run folder (unless `--run-dir` is supplied), with logs in `logs/`, intermediate cubes in `wf_cubes/`, and rendered TGAs in `tga_images/`.
   - By default cube and XYZ intermediates are deleted after rendering to save space; toggle `CONFIG['workflow']['keep_cubes']`/`['keep_xyz']` if you want to keep them.
4. **Rerender with new VMD settings** (reuse saved cubes, skip recomputation):
   ```bash
   python main.py --mode rerender --run-dir runs/20251119/run_143107_6570514d
   ```
5. **Video-only** (assemble existing TGAs):
   ```bash
   python main.py --mode video --run-dir runs/20251119/run_143714_3ee87eb9
   ```

## Output and examples
- Each run creates `runs/<date>/<run_id>/` containing cube files (`wf_cubes/`), rendered images (`tga_images/density_f########.tga`), and logs (`logs/python.log`, `logs/vmd.log`, `logs/ffmpeg.log`).
- A sample movie generated from the included data is at `xsh_density_movie.mp4` (embedded below and also downloadable).
- Example frames live under `runs/20251119/…/tga_images/` showing hole/electron/exciton isosurfaces with the saved camera view.

<p align="center">
  <video controls width="720">
    <source src="xsh_density_movie.mp4" type="video/mp4">
    <!-- GitHub sometimes needs a fully-qualified raw URL for inline playback -->
    <source src="https://raw.githubusercontent.com/kitjoll/xsh_movie_maker/main/xsh_density_movie.mp4" type="video/mp4">
    Your viewer may block inline video; use the download link below.
  </video>
</p>

<p align="center">
  <a href="xsh_density_movie.mp4">Download the MP4</a>
</p>


## Tuning
- `CONFIG['sites']`: set donor/acceptor counts and atom index ranges; update if using different geometries.
- `CONFIG['grid']`: control grid spacing/margins and orbital lobe width (`orbital_resolution`).
- `CONFIG['render']`: atom representation, colors, and isosurface styles used by VMD/Tachyon; adjust before rerendering.
- `CONFIG['video']`: MP4 filename, fps, codec, and optional resolution overrides.
