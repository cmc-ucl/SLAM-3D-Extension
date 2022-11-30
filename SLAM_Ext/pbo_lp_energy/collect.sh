#!/bin/bash

TAR=$1


grep -A 3 "Eigenvalues (eV) |"  $TAR | tail -10
echo "--"
grep "LP Total Energy   (eV) :" $TAR | tail -1
grep "Total             (eV) :" $TAR
