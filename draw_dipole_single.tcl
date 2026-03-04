proc draw_dipole {dipole {radius 0.1}} {
    set sel [atomselect top "all"]
  
    set center_of_mass [measure center $sel weight mass]
   
    set start $center_of_mass 
    set end [vecadd $center_of_mass $dipole]
    #set end [vecadd $dipole $start]
        
    
    # Draw dipole vector (this is a placeholder; actual drawing code depends on the visualization tool)
    set middle [vecadd $start [vecscale 0.75 [vecsub $end $start]]]

    # Draw the cylinder part of the arrow
    graphics top cylinder $start $middle radius $radius

    # Draw the cone part of the arrow (the arrowhead)
    # The cone's radius is typically larger than the cylinder's for a visible arrowhead.
    graphics top cone $middle $end radius [expr $radius * 3]
    #$sel delete
}