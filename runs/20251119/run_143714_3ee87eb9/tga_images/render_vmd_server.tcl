# render_vmd_server.tcl
proc read_view_state {filename} {
    set pairs [list]
    set fh [open $filename r]
    while {[gets $fh line] >= 0} {
        set line [string trim $line]
        if {$line eq ""} {
            continue
        }
        if {[string match "#*" $line]} {
            continue
        }
        set key [lindex $line 0]
        set value [lrange $line 1 end]
        lappend pairs $key $value
    }
    close $fh
    return $pairs
}

set view_file "view_states/saved_view_matrices.txt"
set default_resize {1280 720}

display ambientocclusion on
display dof on
display depthcue off
display projection Orthographic
color Display Background white
axes location Off

if {[string length $view_file] == 0} {
    error "No view_state_file configured. Run capture mode first."
}
if {![file exists $view_file]} {
    error "Saved view file $view_file not found. Run capture mode first."
}

if {[catch {
    set view_pairs [read_view_state $view_file]
    if {[llength $view_pairs] % 2 != 0} {
        error "Malformed view file."
    }
} err]} {
    error "Could not parse saved view from $view_file: $err"
}

array set view_arr $view_pairs

proc apply_saved_view {molID} {
    global view_arr default_resize
    if {[catch {
        molinfo $molID set center_matrix $view_arr(center_matrix)
        molinfo $molID set rotate_matrix $view_arr(rotate_matrix)
        molinfo $molID set scale_matrix $view_arr(scale_matrix)
        molinfo $molID set global_matrix $view_arr(global_matrix)
    } err]} {
        error "Failed to apply matrices from saved view: $err"
    }
    set resize_vals $default_resize
    if {[info exists view_arr(resize)] && [llength $view_arr(resize)] >= 2} {
        set resize_vals $view_arr(resize)
    }
    display resize [lindex $resize_vals 0] [lindex $resize_vals 1]
}

proc render_job {hole_cube elec_cube exc_cube xyz_file out_file} {
    mol new $xyz_file type xyz waitfor all
    set molID [molinfo top get id]

    mol delrep 0 $molID
mol representation Licorice 0.1 150.0 150.0
mol selection all
mol material Opaque
mol color Name
mol addrep $molID

    color Name H white
color Name C black
color Name O red
color Name N blue
color Name S yellow

    mol addfile $hole_cube type cube waitfor all molid $molID
    mol addfile $elec_cube type cube waitfor all molid $molID
    mol addfile $exc_cube  type cube waitfor all molid $molID

    mol addrep $molID
mol modselect 1 $molID all
mol modstyle 1 $molID Isosurface 2e-2 0 0 0 1 50
mol modcolor 1 $molID ColorID 25
mol modmaterial 1 $molID AOChalky
mol addrep $molID
mol modselect 2 $molID all
mol modstyle 2 $molID Isosurface 2e-2 1 0 0 1 50
mol modcolor 2 $molID ColorID 30
mol modmaterial 2 $molID AOChalky
mol addrep $molID
mol modselect 3 $molID all
mol modstyle 3 $molID Isosurface 2e-2 2 0 0 1 50
mol modcolor 3 $molID ColorID 20
mol modmaterial 3 $molID AOChalky

    apply_saved_view $molID

    set scene_file "$out_file"
    set tachyon_exe "/scratch/kjoll/software/vmd/lib/vmd/tachyon_LINUXAMD64"
    set nthreads 8
    set tachyon_cmd [format "\"%s\" -numthreads %d %%s -format TARGA -res 3840 2160 -o %s" $tachyon_exe $nthreads $out_file]
    render Tachyon $scene_file $tachyon_cmd

    mol delete $molID
}

puts stdout "RENDER_SERVER_READY"
flush stdout
while {![eof stdin]} {
    if {[gets stdin line] < 0} {
        break
    }
    set line [string trim $line]
    if {$line eq ""} {
        continue
    }
    if {$line eq "quit"} {
        quit
    }
    set parts [split $line "|"]
    if {[llength $parts] != 5} {
        puts stdout "FRAME_ERROR malformed_job"
        flush stdout
        continue
    }
    set hole_cube [lindex $parts 0]
    set elec_cube [lindex $parts 1]
    set exc_cube  [lindex $parts 2]
    set xyz_file  [lindex $parts 3]
    set out_file  [lindex $parts 4]
    if {[catch {render_job $hole_cube $elec_cube $exc_cube $xyz_file $out_file} err]} {
        puts stdout "FRAME_ERROR $err"
        flush stdout
    } else {
        puts stdout "FRAME_DONE $out_file"
        flush stdout
    }
}
quit
