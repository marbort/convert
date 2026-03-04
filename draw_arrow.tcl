proc vmd_draw_arrow {mol start end {radius 0.1}} {
    # An arrow is made of a cylinder and a cone.
    # The 'middle' point defines where the cylinder ends and the cone begins.
    # It's set at 75% of the distance from 'start' to 'end' by default.
    set middle [vecadd $start [vecscale 0.75 [vecsub $end $start]]]

    # Draw the cylinder part of the arrow
    graphics $mol cylinder $start $middle radius $radius

    # Draw the cone part of the arrow (the arrowhead)
    # The cone's radius is typically larger than the cylinder's for a visible arrowhead.
    graphics $mol cone $middle $end radius [expr $radius * 3]
}
