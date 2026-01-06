# Author: Jeff Comer <jeffcomer at gmail>

set psf [lindex $argv 0]
set pdb [lindex $argv 1]
set massScale [lindex $argv 2]
set selText [lindex $argv 3]
set outPsf [lindex $argv 4]

mol new $psf
mol addfile $pdb


set sel [atomselect top "($selText) and hydrogen"]
puts "Repartitioning mass of [$sel num] hydrogen atoms."

foreach bondList [$sel getbonds] hi [$sel get index] {
    # Get the current mass.
    set hydSel [atomselect top "index $hi"]
    set hydMass [$hydSel get mass]
    set newMass [expr {$massScale*$hydMass}]
    set delMass [expr {$newMass - $hydMass}]

    if {$hydMass > 1.5} {
	    puts "WARNING: Current hydrogen mass $hydMass, new hydrogen mass $newMass. Maybe repartitioning was previously done?"
    }

    # Get the new masses of heavy bonded atoms.
    set s [atomselect top "index $bondList and noh"]
    set massList {}
    set n [$s num]

    foreach m [$s get mass] {
	    lappend massList [expr {$m - $delMass/double($n)}]
    }

    # Set the new masses.
    $hydSel set mass $newMass
    $s set mass $massList
    $s delete
    $hydSel delete
}

set all [atomselect top all]
#Cambie esta linea para que guarde en el archivo de salida
$all writepsf ../build/$outPsf

$sel delete
$all delete
exit
