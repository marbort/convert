#REMOVE DRAWN BONDS WHEN LONGER THAN CUTOFF IN A TRAJECTORY
#Usage: source remove_long_bonds_traj_mio.tcl
#       remove_bonds_traj_mio cutoff output
#Example: remove_bonds_traj_mio 2.0 output.dcd
#Note: cutoff in Angstroms
proc remove_bonds_traj_mio {cutoff} {
    set sel [atomselect top "all"]
    set nframes [molinfo top get numframes]
    for {set i 0} {$i < $nframes} {incr i} {
        $sel frame $i
        set positions [$sel get {x y z}]
        set bonds [$sel getbonds]
        # Find atoms involved in bonds longer than cutoff
        set to_remove {}
        set n [llength $positions]
        for {set j 0} {$j < $n} {incr j} {
            set pos_j [lindex $positions $j]
            for {set k [expr {$j + 1}]} {$k < $n} {incr k} {
                set pos_k [lindex $positions $k]
                set dist [veclength [vecsub $pos_j $pos_k]]
                if {$dist < $cutoff} {
                    lappend to_remove $j
                    lappend to_remove $k
                }
            }
        }
        set to_remove [lsort -unique $to_remove]
        foreach idx $to_remove {
            $sel deleteatom $idx
        }
    }
    #animate write dcd $output sel $sel
    $sel delete
}