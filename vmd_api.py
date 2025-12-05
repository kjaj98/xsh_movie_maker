import subprocess
import threading
from queue import Queue, Empty
import time
from typing import Optional


def _build_scene_config(render_settings, variant=None):
    """
    Merge base render settings with optional variant-specific overrides.
    """
    base_scene = {
        "resize": tuple(render_settings.get("resize", (1920, 1080))),
        "display": dict(render_settings.get("display", {})),
        "atom_representation": dict(render_settings.get("atom_representation", {})),
        "atom_colors": dict(render_settings.get("atom_colors", {})),
        "isosurfaces": list(render_settings.get("isosurfaces", [])),
    }
    if not variant:
        return base_scene

    overrides = render_settings.get(variant, {})
    merged_scene = {
        "resize": tuple(overrides.get("resize", base_scene["resize"])),
        "display": {**base_scene["display"], **overrides.get("display", {})},
        "atom_representation": {
            **base_scene["atom_representation"],
            **overrides.get("atom_representation", {}),
        },
        "atom_colors": {**base_scene["atom_colors"], **overrides.get("atom_colors", {})},
        "isosurfaces": overrides.get("isosurfaces", base_scene["isosurfaces"]),
    }
    return merged_scene


def _format_toggle(value):
    if isinstance(value, str):
        return value
    return "on" if value else "off"


def _format_display_commands(display_cfg):
    lines = []
    if not display_cfg:
        return lines
    for key in ("ambientocclusion", "dof", "depthcue"):
        if key in display_cfg and display_cfg[key] is not None:
            lines.append(f"display {key} {_format_toggle(display_cfg[key])}")
    projection = display_cfg.get("projection")
    if projection:
        lines.append(f"display projection {projection}")
    background = display_cfg.get("background")
    if background:
        lines.append(f"color Display Background {background}")
    axes_location = display_cfg.get("axes_location")
    if axes_location:
        lines.append(f"axes location {axes_location}")
    for extra_cmd in display_cfg.get("extra_commands", ()):
        lines.append(extra_cmd)
    return lines


def _format_atom_colors(atom_colors):
    return [f"color Name {atom} {color}" for atom, color in atom_colors.items()]


def _format_atom_representation(rep_cfg):
    if not rep_cfg:
        return []
    lines = []
    if rep_cfg.get("remove_default", True):
        lines.append("mol delrep 0 $molID")
    representation = rep_cfg.get("representation")
    if representation:
        lines.append(f"mol representation {representation}")
    selection = rep_cfg.get("selection")
    if selection:
        lines.append(f"mol selection {selection}")
    material = rep_cfg.get("material")
    if material:
        lines.append(f"mol material {material}")
    color_mode = rep_cfg.get("color")
    if color_mode:
        lines.append(f"mol color {color_mode}")
    lines.append("mol addrep $molID")
    for extra_cmd in rep_cfg.get("extra_commands", ()):
        lines.append(extra_cmd)
    return lines


def _build_isosurface_style(iso_cfg, fallback_volume_index):
    style = iso_cfg.get("style")
    if style:
        return style
    iso_value = iso_cfg.get("isovalue", 5e-3)
    volume_index = iso_cfg.get("volume_index", fallback_volume_index)
    smoothing = iso_cfg.get("smoothing", 0)
    draw_normals = iso_cfg.get("draw_normals", 0)
    line_style = iso_cfg.get("line_style", 0)
    two_sided = iso_cfg.get("two_sided", 1)
    resolution = iso_cfg.get("resolution", 10)
    return (
        f"Isosurface {iso_value} {volume_index} {smoothing} "
        f"{draw_normals} {line_style} {two_sided} {resolution}"
    )


def _format_isosurface_reps(isosurface_cfg, start_idx):
    if not isosurface_cfg:
        return [], start_idx
    lines = []
    rep_idx = start_idx
    for iso_cfg in isosurface_cfg:
        rep_idx += 1
        selection = iso_cfg.get("selection", "all")
        lines.append("mol addrep $molID")
        lines.append(f"mol modselect {rep_idx} $molID {selection}")
        fallback_volume = iso_cfg.get("volume_index", rep_idx - start_idx - 1)
        style = _build_isosurface_style(iso_cfg, fallback_volume)
        lines.append(f"mol modstyle {rep_idx} $molID {style}")
        color_mode = iso_cfg.get("color")
        if color_mode:
            lines.append(f"mol modcolor {rep_idx} $molID {color_mode}")
        material = iso_cfg.get("material")
        if material:
            lines.append(f"mol modmaterial {rep_idx} $molID {material}")
        for extra_cmd in iso_cfg.get("extra_commands", ()):
            lines.append(extra_cmd)
    return lines, rep_idx


def write_vmd_capture_script(path, render_settings):
    """
    Write a TCL script that launches the GUI with the capture-scene settings and
    provides a Tk dialog for saving the current view matrices to a file.
    """
    scene = _build_scene_config(render_settings, variant="capture_scene")
    resize_w, resize_h = scene["resize"]
    display_section = "\n".join(_format_display_commands(scene["display"]))
    atom_rep_section = "\n".join(_format_atom_representation(scene["atom_representation"]))
    atom_color_section = "\n".join(_format_atom_colors(scene["atom_colors"]))
    iso_lines, _ = _format_isosurface_reps(scene["isosurfaces"], start_idx=0)
    iso_section = "\n".join(iso_lines)

    script = f"""# capture_view.tcl
set hole_cube [lindex $argv 0]
set elec_cube [lindex $argv 1]
set exc_cube  [lindex $argv 2]
set xyz_file  [lindex $argv 3]
set save_file [lindex $argv 4]

display resize {resize_w} {resize_h}
{display_section}

mol new $xyz_file type xyz waitfor all
set molID [molinfo top get id]

# Atom representation
{atom_rep_section}

# Custom atom colors
{atom_color_section}

mol addfile $hole_cube type cube waitfor all molid $molID
mol addfile $elec_cube type cube waitfor all molid $molID
mol addfile $exc_cube  type cube waitfor all molid $molID

# Isosurface representations
{iso_section}

package require Tk
proc save_view {{molID outfile}} {{
    set fh [open $outfile w]
    puts $fh "# Saved VMD view (resize + matrices)"
    puts $fh "resize [display get size]"
    puts $fh "center_matrix [molinfo $molID get center_matrix]"
    puts $fh "rotate_matrix [molinfo $molID get rotate_matrix]"
    puts $fh "scale_matrix [molinfo $molID get scale_matrix]"
    puts $fh "global_matrix [molinfo $molID get global_matrix]"
    close $fh
    vmdcon -info "Saved camera matrices to $outfile"
    quit
}}

toplevel .capture
wm title .capture "Capture VMD View"
label .capture.msg -text "1) Orient the camera as desired.\\n2) Click 'Save view' to store it for future renders."
button .capture.save -text "Save view & Quit" -command {{save_view $molID $save_file}}
button .capture.cancel -text "Cancel" -command {{quit}}
pack .capture.msg -padx 12 -pady 8
pack .capture.save -side left -padx 12 -pady 8
pack .capture.cancel -side right -padx 12 -pady 8
"""
    with open(path, "w") as f:
        f.write(script)


def write_vmd_server_script(path, tachyon_exe_path, render_settings, view_state_file=None):
    """
    Write a TCL script that keeps VMD alive and renders multiple frames from stdin.
    """
    scene = _build_scene_config(render_settings)
    resize_w, resize_h = scene["resize"]
    display_section = "\n".join(_format_display_commands(scene["display"]))
    atom_rep_section = "\n".join(_format_atom_representation(scene["atom_representation"]))
    atom_color_section = "\n".join(_format_atom_colors(scene["atom_colors"]))
    iso_lines, _ = _format_isosurface_reps(scene["isosurfaces"], start_idx=0)
    iso_section = "\n".join(iso_lines)
    view_file = view_state_file or ""

    script = f"""# render_vmd_server.tcl
proc read_view_state {{filename}} {{
    set pairs [list]
    set fh [open $filename r]
    while {{[gets $fh line] >= 0}} {{
        set line [string trim $line]
        if {{$line eq ""}} {{
            continue
        }}
        if {{[string match "#*" $line]}} {{
            continue
        }}
        set key [lindex $line 0]
        set value [lrange $line 1 end]
        lappend pairs $key $value
    }}
    close $fh
    return $pairs
}}

set view_file "{view_file}"
set default_resize {{{resize_w} {resize_h}}}

{display_section}

if {{[string length $view_file] == 0}} {{
    error "No view_state_file configured. Run capture mode first."
}}
if {{![file exists $view_file]}} {{
    error "Saved view file $view_file not found. Run capture mode first."
}}

if {{[catch {{
    set view_pairs [read_view_state $view_file]
    if {{[llength $view_pairs] % 2 != 0}} {{
        error "Malformed view file."
    }}
}} err]}} {{
    error "Could not parse saved view from $view_file: $err"
}}

array set view_arr $view_pairs

proc apply_saved_view {{molID}} {{
    global view_arr default_resize
    if {{[catch {{
        molinfo $molID set center_matrix $view_arr(center_matrix)
        molinfo $molID set rotate_matrix $view_arr(rotate_matrix)
        molinfo $molID set scale_matrix $view_arr(scale_matrix)
        molinfo $molID set global_matrix $view_arr(global_matrix)
    }} err]}} {{
        error "Failed to apply matrices from saved view: $err"
    }}
    set resize_vals $default_resize
    if {{[info exists view_arr(resize)] && [llength $view_arr(resize)] >= 2}} {{
        set resize_vals $view_arr(resize)
    }}
    display resize [lindex $resize_vals 0] [lindex $resize_vals 1]
}}

proc render_job {{hole_cube elec_cube exc_cube xyz_file out_file}} {{
    mol new $xyz_file type xyz waitfor all
    set molID [molinfo top get id]

    {atom_rep_section}

    {atom_color_section}

    mol addfile $hole_cube type cube waitfor all molid $molID
    mol addfile $elec_cube type cube waitfor all molid $molID
    mol addfile $exc_cube  type cube waitfor all molid $molID

    {iso_section}

    apply_saved_view $molID

    set scene_file "$out_file"
    set tachyon_exe "{tachyon_exe_path}"
    set nthreads 8
    set tachyon_cmd [format "\\\"%s\\\" -numthreads %d %%s -format TARGA -res 3840 2160 -o %s" $tachyon_exe $nthreads $out_file]
    render Tachyon $scene_file $tachyon_cmd

    mol delete $molID
}}

puts stdout "RENDER_SERVER_READY"
flush stdout
while {{![eof stdin]}} {{
    if {{[gets stdin line] < 0}} {{
        break
    }}
    set line [string trim $line]
    if {{$line eq ""}} {{
        continue
    }}
    if {{$line eq "quit"}} {{
        quit
    }}
    set parts [split $line "|"]
    if {{[llength $parts] != 5}} {{
        puts stdout "FRAME_ERROR malformed_job"
        flush stdout
        continue
    }}
    set hole_cube [lindex $parts 0]
    set elec_cube [lindex $parts 1]
    set exc_cube  [lindex $parts 2]
    set xyz_file  [lindex $parts 3]
    set out_file  [lindex $parts 4]
    if {{[catch {{render_job $hole_cube $elec_cube $exc_cube $xyz_file $out_file}} err]}} {{
        puts stdout "FRAME_ERROR $err"
        flush stdout
    }} else {{
        puts stdout "FRAME_DONE $out_file"
        flush stdout
    }}
}}
quit
"""
    with open(path, "w") as f:
        f.write(script)


def capture_vmd_view(vmd_exe_path, capture_script_path, hole_file, elec_file, exc_file,
                     xyz_file, view_state_file):
    """
    Launch VMD with GUI to let the user orient the molecule and save the camera matrices.
    """
    cmd = [
        vmd_exe_path,
        "-e", capture_script_path,
        "-args",
        hole_file,
        elec_file,
        exc_file,
        xyz_file,
        view_state_file,
    ]
    print("Launching VMD GUI for camera capture...")
    subprocess.run(cmd, check=True)


class VMDRenderServer:
    """
    Maintain a persistent VMD process to render multiple frames without
    re-launching the application for every image.
    """

    def __init__(self, vmd_exe_path, server_script_path,
                 startup_timeout=60, frame_timeout=600, log_path: Optional[str] = None):
        self._frame_timeout = frame_timeout
        self._log_handle = open(log_path, "w") if log_path else None
        self._echo_vmd = log_path is None
        self._process = subprocess.Popen(
            [
                vmd_exe_path,
                "-dispdev", "text",
                "-e", server_script_path,
            ],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        self._stdout_queue: Queue = Queue()
        self._stdout_thread = threading.Thread(
            target=self._drain_stdout, daemon=True
        )
        self._stdout_thread.start()
        try:
            self._wait_for_message("RENDER_SERVER_READY", startup_timeout)
        except Exception:
            self.close()
            raise

    def _log_line(self, line: str) -> None:
        if self._log_handle is not None:
            self._log_handle.write(line + "\n")
            self._log_handle.flush()

    def _drain_stdout(self):
        if not self._process.stdout:
            return
        for line in self._process.stdout:
            cleaned = line.rstrip("\n")
            self._log_line(cleaned)
            self._stdout_queue.put(cleaned)
        self._stdout_queue.put(None)

    def _wait_for_line(self, timeout):
        try:
            return self._stdout_queue.get(timeout=timeout)
        except Empty:
            raise TimeoutError("Timed out waiting for VMD output.")

    def _wait_for_message(self, prefix, timeout):
        end = time.monotonic() + timeout
        while True:
            remaining = end - time.monotonic()
            if remaining <= 0:
                raise TimeoutError(f"Timed out waiting for '{prefix}'.")
            line = self._wait_for_line(remaining)
            if line is None:
                raise RuntimeError("VMD render server terminated unexpectedly.")
            if line.startswith(prefix):
                return line
            if self._echo_vmd:
                print(f"[VMD] {line}")

    def _ensure_alive(self):
        if self._process.poll() is not None:
            raise RuntimeError("VMD render server has exited.")

    def render_frame(self, hole_file, elec_file, exc_file, xyz_file, img_file):
        """
        Submit a render job and block until the frame is finished.
        """
        self._ensure_alive()
        job_line = "|".join([hole_file, elec_file, exc_file, xyz_file, img_file])
        try:
            self._process.stdin.write(job_line + "\n")
            self._process.stdin.flush()
        except Exception as exc:
            raise RuntimeError(f"Failed to dispatch frame to VMD: {exc}") from exc

        end = time.monotonic() + self._frame_timeout
        while True:
            remaining = end - time.monotonic()
            if remaining <= 0:
                raise TimeoutError("Timed out waiting for VMD frame completion.")
            line = self._wait_for_line(remaining)
            if line is None:
                raise RuntimeError("VMD render server terminated unexpectedly.")
            if line.startswith("FRAME_DONE"):
                print(f"  Rendered TGA: {img_file}")
                return
            if line.startswith("FRAME_ERROR"):
                raise RuntimeError(f"VMD render server error: {line}")
            if self._echo_vmd:
                print(f"[VMD] {line}")

    def close(self):
        """
        Shut down the persistent VMD process.
        """
        if self._process and self._process.poll() is None:
            try:
                if self._process.stdin:
                    self._process.stdin.write("quit\n")
                    self._process.stdin.flush()
            except Exception:
                pass
            try:
                self._process.wait(timeout=10)
            except subprocess.TimeoutExpired:
                self._process.kill()
        if self._stdout_thread and self._stdout_thread.is_alive():
            self._stdout_thread.join(timeout=1)
        if self._log_handle is not None:
            try:
                self._log_handle.close()
            finally:
                self._log_handle = None
