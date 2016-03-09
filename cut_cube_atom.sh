#!/bin/bash

n1=0
n2=999
p=1
c=0
eps=2.0
verbose=0
usage="usage : $0 [options] in.cube > out.cube

options are:
  -i <int>    from i...
  -j <int>    to ...j
  -e <float>  radius around atom to cut
  -r <float>  see -e
  -p          print output (bool)
  -c          use center of mass of selected atoms
  -v          verbose
  -h          show help
" 
#echo "this script only works for orthorhombic cells" > "/dev/stderr"
while getopts "i:j:e:r:pcvh" opt; do
  case $opt in 
    i ) n1=$OPTARG;;
    j ) n2=$OPTARG;;
    e ) eps=$OPTARG;;
    r ) eps=$OPTARG;;
    p ) p=$((1-p));;
    c ) c=$((1-c));;
    v ) ((verbose++));;
    h ) echo -e "$usage"; exit 0;;
    * ) echo "Unimplented option."; exit 1;;
  esac
done
shift $(($OPTIND-1))

if [ $# -lt 1 ]; then
  echo "no filename given."
  exit
fi

> "/dev/stderr"
awk -v n1=$n1 -v n2=$n2 -v output=$p -v cuteps=$eps -v verbose=$verbose -v center=$c '{
  #### SPECIFY OPTIONS #####
  # specify atoms : minimal distance from points to atoms n1..n2
  # (these are vmd indices: starting from 0)
  #n1=360
  #n2=362
  # cut +/- eps (in Angstrom) around the group of atoms
  #cuteps=2.0
  # print cube file to stdout
  #output=1
  ##########################

  #### PROCESS HEADER ######
  if(output) print; getline; if(output) print 
  getline; if(output) print; N=$1; origin[0]=$2; origin[1]=$3; origin[2]=$4
  getline; if(output) print; nx=$1; dx=$2
  getline; if(output) print; ny=$1; dy=$3
  getline; if(output) print; nz=$1; dz=$4
  if(n2>=N)n2=N-1

  if(nx<0) { # angstrom units
    nx=-nx;ny=-ny;nz=-nz;
    a0 = 1.
  } else { a0=.5291772108; }; # units are bohr
  cuteps/=a0
  cuteps2=cuteps*cuteps
  lx=nx*dx;ly=ny*dy;lz=nz*dz

  if(verbose)print "COORDINATES:" > "/dev/stderr"
  for(i=0;i<N;i++) {
    getline; if(output) print
    x[i]=$3; y[i]=$4; z[i]=$5
    cm[0]+=$3/N;cm[1]+=$4/N;cm[2]+=$5/N
    if(verbose)printf "%3d   %.2f %.2f %.2f\n",i,x[i],y[i],z[i] > "/dev/stderr"
  }
  getline
  printf "origin of system : (% 8.4f, % 8.4f, % 8.4f)\n",origin[0],origin[1],origin[2] > "/dev/stderr"
  printf "center of mass   : (% 8.4f, % 8.4f, % 8.4f)\n",cm[0],cm[1],cm[2] > "/dev/stderr"
  printf "cell size        : (% 8.4f, % 8.4f, % 8.4f)\n",lx,ly,lz > "/dev/stderr"
  for(j=n1;j<=n2;j++) {
    printf "atom %d : pos = (%.4f, %.4f, %.4f)\n",j,x[j],y[j],z[j] > "/dev/stderr"
    if(center) { xc+=x[j];yc+=y[j];zc=z[j];nc++ }
  }
  xc/=nc; yc/=nc; zc/=nc
  if(center) printf "cutting %.3f around center of atoms %d to %d : (%.4f, %.4f, %.4f)\n",cuteps,n1,n2,xc,yc,zc > "/dev/stderr"
  else printf "cutting %.3f around atoms %d to %d\n",cuteps,n1,n2 > "/dev/stderr"
  
  ## for height...
  v[1] = y[n1]-cm[1]
  v[2] = z[n1]-cm[2]
  l = sqrt(v[0]**2+v[1]**2)
  v[1] /= l
  v[2] /= l

  ##########################
  go=1
  ix=0;iy=0;iz=0
  while(go) {
    for(i=1;i<=NF;i++) {
      ##########################
      # MAIN PROGRAM FOLLOWS : # 
      ##########################
      value=0
      if(center) {
          deltax = xc-ix*dx-origin[0]; if(deltax>lx/2) deltax-=lx; else if(deltax<-lx/2)deltax += lx
          deltay = yc-iy*dy-origin[1]; if(deltay>ly/2) deltay-=ly; else if(deltay<-ly/2)deltay += ly
          deltaz = zc-iz*dz-origin[2]; if(deltaz>lz/2) deltaz-=lz; else if(deltaz<-lz/2)deltaz += lz
          delta = deltax*deltax + deltay*deltay + deltaz*deltaz
          if(delta<cuteps2) value=$i
      } else {
        for(j=n1;j<=n2;j++) {
          deltax = x[j]-ix*dx-origin[0]; if(deltax>lx/2) deltax-=lx; else if(deltax<-lx/2)deltax += lx
          deltay = y[j]-iy*dy-origin[1]; if(deltay>ly/2) deltay-=ly; else if(deltay<-ly/2)deltay += ly
          deltaz = z[j]-iz*dz-origin[2]; if(deltaz>lz/2) deltaz-=lz; else if(deltaz<-lz/2)deltaz += lz
          delta = deltax*deltax + deltay*deltay + deltaz*deltaz
          if(delta<cuteps2) {value=$i; break}
        }
      }
      ### using height defined as YZ-plane to atom n1
      #height = v[1]*(iy*dy+origin[1]-cm[1]) + v[2]*(iz*dz+origin[2]-cm[2]) # dot prod of v and pos
      #if(height>11) value=0
      ######### DONE ###########
      sum+=value
      sum2+=value*value
      if(output) printf " % 0.5e",value;  
      sumtotal+=$i
      sum2total+=$i*$i
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
  if(sum2total<1.1 && sum2total>0.9) {
    print "sum2total = ",sum2total*dv  > "/dev/stderr"
    print "sum2      = ",sum2*dv       > "/dev/stderr"
  }
}' $1
