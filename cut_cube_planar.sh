#!/bin/bash

echo "this script only works for orthorhombic cells" > "/dev/stderr"
if [ $# -lt 1 ]; then
  echo "./cut_cube.awk in.cube > out.cube"
  exit
fi

awk '{
  #### SPECIFY OPTIONS #####
  # specify three points / atom indices to define the plane
  # (these are vmd indices: starting from 0)
  n1=1
  n2=181
  n3=179
  # cut +/- eps (in Angstrom) around the plane
  cuteps=0.1
  # print cube file to stdout
  output=1
  ##########################

  #### PROCESS HEADER ######
  if(output) print; getline; if(output) print 
  getline; if(output) print; N=$1; origin[0]=$2; origin[1]=$3; origin[2]=$4
  getline; if(output) print; nx=$1; dx=$2
  getline; if(output) print; ny=$1; dy=$3
  getline; if(output) print; nz=$1; dz=$4

  if(nx<0) { # angstrom units
    nx=-nx;ny=-ny;nz=-nz;
    a0 = 1.
  } else { a0=.5291772108; }; # units are bohr
  cuteps/=a0

  for(i=0;i<N;i++) {
    getline; if(output) print
    if(NR-7==n1) { x[0]=$3; x[1]=$4; x[2]=$5 }
    if(NR-7==n2) { y[0]=$3; y[1]=$4; y[2]=$5 }
    if(NR-7==n3) { z[0]=$3; z[1]=$4; z[2]=$5 }
    cm[0]+=$3/N;cm[1]+=$4/N;cm[2]+=$5/N
  }
  getline
  print "origin of system    : ",origin[0],origin[1],origin[2] > "/dev/stderr"
  print "center of mass      : ",cm[0],cm[1],cm[2] > "/dev/stderr"
  print "position of point 1 : ",x[0],x[1],x[2] > "/dev/stderr"
  print "position of point 2 : ",y[0],y[1],y[2] > "/dev/stderr"
  print "position of point 3 : ",z[0],z[1],z[2] > "/dev/stderr"

  ##########################

  # compute tangential to the plane, v = cross(z-x,z-y) = [a,b,c]
  a = (z[1]-x[1])*(z[2]-y[2]) - (z[2]-x[2])*(z[1]-y[1])
  b = (z[2]-x[2])*(z[0]-y[0]) - (z[0]-x[0])*(z[2]-y[2])
  c = (z[0]-x[0])*(z[1]-y[1]) - (z[1]-x[1])*(z[0]-y[0])
  vnorm = sqrt(a*a+b*b+c*c)
  a/=vnorm; b/=vnorm; c/=vnorm
  d=-(a*(x[0]-origin[0])+b*(x[1]-origin[1])+c*(x[2]-origin[2])) # distance to x
  # the plane obtained is defined by
  # ax + by + cz + d = 0 
  ##########################
  go=1
  ix=0;iy=0;iz=0
  while(go) {
    for(i=1;i<=NF;i++) {
      value=$i
      ##########################
      # MAIN PROGRAM FOLLOWS : # 
      ##########################
      #density[ix, iy, iz]=value # store density if needed for some other post-processing
      dist = a*(ix*dx)+b*(iy*dz)+c*(iz*dz)+d # distance from point to plane
      if(dist*dist<cuteps) print ix*dx,iy*dy,iz*dz,dist,value > "cut.dat" # output for plotting
      else value=0; # remove for output-cube and partial sum
      ######### DONE ###########
      sum+=value
      if(output) printf " % 0.5e",value;  
      sumtotal+=$i
      iz++
      if(iz>=nz) { 
        iz=0
        iy++
        if(iy>=ny) {
          iy=0
          ix++
          if(ix>=nx) go=0
        }
      }
      if(!go) break
    }
    if(output) printf "\n"
    getline
  }
  dv=dx*dy*dz # volume element for integration
  print "sumtotal = ",sumtotal*dv  > "/dev/stderr"
  print "sum      = ",sum*dv       > "/dev/stderr"
}' $1
