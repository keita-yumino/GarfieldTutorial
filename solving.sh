#!/bin/bash

GeoFile=gemcell

#gmsh ${GeoFile}.geo -3 -optimize -order 2
#gmsh ${GeoFile}.geo -3 -order 2
#gmsh ${GeoFile}.geo -3 -order 2 -format msh2
#gmsh ${GeoFile}.geo -3 -order 2 -optimize_netgen
#gmsh ${GeoFile}.geo -3 -order 2 -format msh2 -optimize_netgen
gmsh ${GeoFile}.geo -3 -order 2  -optimize_netgen
ElmerGrid 14 2 ${GeoFile}.msh -autoclean
ElmerSolver ${GeoFile}.sif

