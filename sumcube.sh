#!/bin/bash

mode=1
p=1
usage="usage : $0 [options] file1.cube [ file2.cube [...]] > out.cube

options are:
  -m          mode: -1=difference, 1=sum; (default=$mode) 
  -v          verbose
  -p          print output (bool)
  -h          show help
" 

while getopts "m:vph" opt; do
  case $opt in 
    m ) mode="$OPTARG";;
    v ) ((verbose++));;
    p ) p=$((1-p));;
    h ) echo -e "$usage"; exit 0;;
    * ) echo "Unimplented option."; exit 1;;
  esac
done
shift $(($OPTIND-1))

if [ $# -lt 1 ]; then echo "no input file."; exit; fi

nat=`awk '{if(NR==3) { print $1; exit}}' $1`
head -n $((nat+6)) $1

awk -v verbose="$verbose" -v p="$p" -v mode="$mode" '{
  f=FILENAME
  flist[++fn]=f
  print "reading ",f > "/dev/stderr"

  #### PROCESS HEADER ######
  if(fn==1)getline; 
  getline; N=$1; origin[0]=$2; origin[1]=$3; origin[2]=$4
  getline; nx=$1; dx=$2
  getline; ny=$1; dy=$3
  getline; nz=$1; dz=$4
  ntotal = nx*ny*nz
  if(fn==1 && verbose) printf "program will use approximately %.1f Mbytes of memory\n",ntotal/1024/1024*8*24 > "/dev/stderr"
  if(nx<0) { # angstrom units
    nx=-nx;ny=-ny;nz=-nz;
    a0 = 1.
  } else { a0=.5291772108; }; # units are bohr

  cm[0]=0;cm[1]=0;cm[2]=0
  for(i=0;i<N;i++) {
    getline;
    cm[0]+=$3/N;cm[1]+=$4/N;cm[2]+=$5/N
  }
  getline
  if(verbose) print "center of mass      : ",cm[0],cm[1],cm[2] > "/dev/stderr"
  ##########################

  iter=0; sum=0; sumtotal=0
  while(iter<ntotal) {
    for(i=1;i<=NF;i++) {
      value=$i
      ##########################
      # MAIN PROGRAM FOLLOWS : # 
      ##########################
      if(fn==1) newrho[iter] = value
      else {
        newrho[iter] += mode*value 
        #if(mode==0) newrho[iter] -= value 
        #else newrho[iter] += value 
      }
      # 1D-ARRAY TO REDUCE MEMORY FOOTPRINT #
      ######### DONE ###########
      sum+=value
      sumtotal+=newrho[iter]
      iter++
    }
    getline
  }
  dv=dx*dy*dz # volume element for integration
  if(verbose) printf "integral = %f, total = %f\n",sum*dv,sumtotal*dv > "/dev/stderr"
} END {
  if(p) {
    print "writing result to output." > "/dev/stderr"
    system("sleep 20")
    i=0;iter=0
    for(ix=0;ix<nx;ix++) {
      for(iy=0;iy<ny;iy++) {
        for(iz=0;iz<nz;iz++) {
          printf "%13.5E",newrho[iter]
          if(iz%6==5) printf "\n"
          iter++
        } 
        printf "\n"
      }
    }
  }
}' $*
