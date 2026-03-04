proc draw_dipole_shifted {dipole {radius 0.1}} {
    set sel [atomselect top "all"]
    set numframes [molinfo top get numframes]
    set sum_com {0 0 0}
    for {set i 0} {$i < $numframes} {incr i} {
    set com [measure center $sel weight mass]
    set sum_com [vecadd $sum_com $com]
}
    
    set center_of_mass [vecscale [expr {1.0 / $numframes}] $sum_com]
    puts "Center of Mass: $center_of_mass"
    #set center_of_mass [measure center $sel weight mass]
    set middle_point [vecscale 0.5 [vecadd $center_of_mass $dipole] ]
    set half_lgth [veclength [vecscale 0.5 $dipole]]
    set shift_vector [vecscale 0.5 [vecsub $middle_point $center_of_mass]]
    set start [vecsub $center_of_mass $shift_vector]
    set end $middle_point
    #set end [vecadd $dipole $start]
        
    # Shift dipole vector
    set dipole_shifted [vecsub $dipole $middle_point]
    # Draw dipole vector (this is a placeholder; actual drawing code depends on the visualization tool)
    set middle [vecadd $start [vecscale 0.75 [vecsub $end $start]]]

    # Draw the cylinder part of the arrow
    graphics top cylinder $start $middle radius $radius

    # Draw the cone part of the arrow (the arrowhead)
    # The cone's radius is typically larger than the cylinder's for a visible arrowhead.
    graphics top cone $middle $end radius [expr $radius * 3]
    #$sel delete
}