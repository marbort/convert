# User-defined parameters
set molID top                   ;# Molecule ID (assumes only one molecule loaded)
set cutoff 10.0                ;# Distance cutoff (in Ångströms)
set outputPrefix "filtered"  ;# Output filename prefix

# Load TopoTools
package require topotools

# Get number of frames
set nframes [molinfo $molID get numframes]

# Loop over each frame
for {set f 0} {$f < $nframes} {incr f} {
    puts "Processing frame $f"
    animate goto $f

    # Get current frame coordinates
    set sel [atomselect $molID "all"]
    set coords [$sel get {x y z}]

    # Clear existing bonds to start fresh
    topo clearbonds

    # Get bonds as guessed or from the topology
    set bonds [$sel getbonds]

    # Rebuild bonds, filtering out long ones
    set newbonds {}
    foreach atom1_bonds $bonds atom1_index [range 0 [expr {[llength $bonds] - 1}]] {
        set filtered_bonds {}
        foreach atom2_index $atom1_bonds {
            set pos1 [lindex $coords $atom1_index]
            set pos2 [lindex $coords $atom2_index]
            set dist [veclength [vecsub $pos1 $pos2]]
            if {$dist <= $cutoff} {
                lappend filtered_bonds $atom2_index
            }
        }
        lappend newbonds $filtered_bonds
    }

    # Apply the filtered bonds
    $sel setbonds $newbonds

    # Optionally write to a new frame/structure file
    #set outname "${outputPrefix}_frame${f}.pdb"
    #$sel writepdb $outname

    $sel delete
}

puts "Bond filtering completed."

