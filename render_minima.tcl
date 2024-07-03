set fp [open "structures.dat" r]
set file_data [read $fp]
close $fp
set data [split $file_data "\n"]
foreach struct $data {
     animate goto $struct
     #display resetview
     animate write pdb $struct.pdb beg $struct end $struct top
     label hide Atoms all
     label hide Bonds all
     render POV3 $struct.pov povray +W%w +H%h -I%s -O$struct.png +D +X +A +FN +UA
}

