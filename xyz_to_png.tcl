proc rep_charged {type} {
    set files [lsort [glob *.$type]]
    foreach file $files {
        mol new $file
        mol modstyle 0 top CPK 1.000000 0.300000 32.000000 32.000000
        render snapshot $file.png 
        mol off top
    }
}
