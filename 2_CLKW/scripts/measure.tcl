#############################################################
## Center and dimensions                                   ##
#############################################################

set name ../Output_Build/4-4_CLKWWW_center
mol load psf $name.psf pdb $name.pdb
set all [atomselect top all]
set center [measure center $all]
set dim [measure minmax $all]

puts "$center"
puts "$dim"

