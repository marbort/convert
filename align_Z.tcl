proc align_x {one two} {
    set A [atomselect top "index $one"]
    set B [atomselect top "index $two"]
    set C [vecsub [lindex [$A get {x y z}] 0] [lindex [$B get {x y z}] 0]]

    set sel [atomselect top all]
    set M [transvecinv $C]
    $sel move $M
    set M [transaxis y -90]
    $sel move $M
}