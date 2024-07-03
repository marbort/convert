proc remove_long_bonds { max_length } {
    set n [ molinfo top get numatoms ]
    for { set i 0 } { $i < $n } { incr i } {
        set bead [ atomselect top "index $i" ]
        set bonds [ lindex [$bead getbonds] 0 ]
        if { $i % 1000 == 0 || $i == $n-1 } then {
            vmdcon -info "$i / $n"
        }
        if { [ llength bonds ] > 0 } {
            set bonds_new {}
            set xyz [ lindex [$bead get {x y z}] 0 ]
            foreach j $bonds {
                set bead_to [ atomselect top "index $j" ]
                set xyz_to [ lindex [$bead_to get {x y z}] 0 ]
                $bead_to delete
                if { [ vecdist $xyz $xyz_to ] < $max_length } {
                    lappend bonds_new $j
                }
            }
            $bead setbonds [ list $bonds_new ]
        }
        $bead delete
    }
    vmdcon -info "$n / $n"
}
