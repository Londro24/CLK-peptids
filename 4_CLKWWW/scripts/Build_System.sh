#!/bin/bash

name="4_CLKWWW"
surface="step1_nanomaterial"
dir="../build"

vmd -dispdev text -e Build_System.tcl -args "${name}" "${surface}" "${dir}"
vmd -dispdev text -e repartitionHydMasses.tcl -args ${dir}/final-sist-${name}_ion.psf ${dir}/final-sist-${name}_ion.pdb 3.0 "hydrogen and not water and not resname TIP3 SPCE" repart_${name}.psf
#vmd -e Build_System.tcl -args "${name}" "${surface}" "${dir}"
