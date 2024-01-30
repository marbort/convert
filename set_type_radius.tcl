set carbon [atomselect top "type 1"]
set oxygen [atomselect top "type 2"]
set hydrogen [atomselect top "type 3"]
set magnesium [atomselect top "type 4"]
set chlorine [atomselect top "type 5"]

$carbon set radius 1.6
$oxygen set radius 1.5
$hydrogen set radius 1.1
$magnesium set radius 1.2
$chlorine set radius 2.0