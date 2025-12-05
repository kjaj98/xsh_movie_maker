import os
import subprocess
import shlex
from dataclasses import dataclass
from typing import Optional, Tuple, Sequence


@dataclass
class VideoSettings:
    # ffmpeg input pattern, e.g. "density_f%05d_t%dp0.tga"
    input_pattern: str
    output_file: str = "output.mp4"

    # basic video settings
    fps: int = 24
    codec: str = "libx264"
    crf: int = 1          
    preset: str = "medium" # encoding speed/size tradeoff

    # optional extras
    pix_fmt: str = "yuv420p"
    resolution: Optional[Tuple[int, int]] = None  # e.g. (1920, 1080)
    extra_args: Sequence[str] = ()                # any extra ffmpeg args

def make_video_from_tga_sequence(
    ffmpeg_path: str,
    settings: VideoSettings,
    dry_run: bool = False,
    log_path: Optional[str] = None,
):
    """
    Build and run an ffmpeg command based on the given settings.

    If dry_run=True, it only prints the command and returns it instead of running.
    """

    cmd = [
        f"{ffmpeg_path}",
        "-y",  # overwrite output without asking
        "-framerate", str(settings.fps),
        "-i", settings.input_pattern,
    ]

    # Optional scaling
    if settings.resolution is not None:
        w, h = settings.resolution
        cmd += ["-vf", f"scale={w}:{h}"]

    # Encoding settings
    cmd += [
        "-c:v", settings.codec,
        "-crf", str(settings.crf),
        "-preset", settings.preset,
        "-pix_fmt", settings.pix_fmt,
    ]

    # Any extra user-specified options
    cmd += list(settings.extra_args)

    # Output
    cmd.append(settings.output_file)

    if dry_run:
        # Show the exact shell command that would be run
        print(" ".join(shlex.quote(c) for c in cmd))
        return cmd

    if log_path is not None:
        log_dir = os.path.dirname(os.path.abspath(log_path))
        os.makedirs(log_dir, exist_ok=True)
        with open(log_path, "w") as log_file:
            subprocess.run(cmd, check=True, stdout=log_file, stderr=subprocess.STDOUT)
    else:
        subprocess.run(cmd, check=True)
    return settings.output_file
