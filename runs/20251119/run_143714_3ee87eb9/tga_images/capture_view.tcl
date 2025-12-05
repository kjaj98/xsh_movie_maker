# capture_view.tcl
set hole_cube [lindex $argv 0]
set elec_cube [lindex $argv 1]
set exc_cube  [lindex $argv 2]
set xyz_file  [lindex $argv 3]
set save_file [lindex $argv 4]

display resize 1280 720
display projection Orthographic
color Display Background white

mol new $xyz_file type xyz waitfor all
set molID [molinfo top get id]

# Atom representation
mol delrep 0 $molID
mol representation Licorice 0.1 20.0 20.0
mol selection all
mol material Opaque
mol color Name
mol addrep $molID

# Custom atom colors
color Name H white
color Name C black
color Name O red
color Name N blue
color Name S yellow

mol addfile $hole_cube type cube waitfor all molid $molID
mol addfile $elec_cube type cube waitfor all molid $molID
mol addfile $exc_cube  type cube waitfor all molid $molID

# Isosurface representations
mol addrep $molID
mol modselect 1 $molID all
mol modstyle 1 $molID Isosurface 5e-2 0 0 0 1 10
mol modcolor 1 $molID ColorID 25
mol modmaterial 1 $molID AOChalky
mol addrep $molID
mol modselect 2 $molID all
mol modstyle 2 $molID Isosurface 5e-2 1 0 0 1 10
mol modcolor 2 $molID ColorID 30
mol modmaterial 2 $molID AOChalky
mol addrep $molID
mol modselect 3 $molID all
mol modstyle 3 $molID Isosurface 5e-2 2 0 0 1 10
mol modcolor 3 $molID ColorID 20
mol modmaterial 3 $molID AOChalky

package require Tk
proc save_view {molID outfile} {
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
}

toplevel .capture
wm title .capture "Capture VMD View"
label .capture.msg -text "1) Orient the camera as desired.\n2) Click 'Save view' to store it for future renders."
button .capture.save -text "Save view & Quit" -command {save_view $molID $save_file}
button .capture.cancel -text "Cancel" -command {quit}
pack .capture.msg -padx 12 -pady 8
pack .capture.save -side left -padx 12 -pady 8
pack .capture.cancel -side right -padx 12 -pady 8
