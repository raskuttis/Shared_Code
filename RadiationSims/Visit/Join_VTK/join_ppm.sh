#!/bin/bash
minuso=$1
outbasename=$2
inbasename=$3
varid=$4
layout=$5
s=$6
e=$7

usage_err="Usage:  join_ppm.sh -o <OUTBASENAME> <INBASENAME> <IDSTRING> <LAYOUT> <STARTCYCLE> <ENDCYCLE>"
layout_err="<LAYOUT> must be in the form MxN, where M and N are integers and the number of processors = M*N"

declare -a proclist
proclist=( `find id* -name "*.ppm" | cut -d '/' -f 1 | cut -d 'd' -f 2 | sort | uniq` )
#proclist=( `find id* -name "*.ppm" | cut -b 3-3 | sort | uniq` )
#proclist=( `find id* -name "*.ppm" | cut -b 3-4 | sort | uniq` )
echo ${proclist[@]}

if [[ ${layout##*x*} != "" ]]
  then echo -e "$usage_err \n $layout_err"
  exit
else
  nx=${layout#*x*}
  ny=${layout%*x*}
  nproc=${#proclist[@]}
#  if [ $((nx*ny)) != $nproc ]
#    then echo -e "$usage_err \n $layout_err"
#    exit
#  fi
fi

if [[ $minuso -ne "-o" ]]
  then echo -e $usage_err 
  exit
elif (( "$s" < "0" ))
  then echo -e $usage_err 
  exit
elif (( "$e" < "$s" ))
  then echo -e $usage_err 
  exit
fi

framelist=""
for ((c=$s; c<=$e; c++))
do
  num=$(printf "%04d" "$c")
  arglist=""

  for ((m=$ny-1; m>=0; m--))
  do
    for ((n=0; n<=$nx-1; n++))
    do
      p=$((m*nx+n))
      if [ ${proclist[p]} -eq 0 ]; then
        arglist="$arglist id0/$inbasename.$num.$varid.ppm"
      else
        arglist="$arglist id${proclist[p]}/$inbasename-id${proclist[p]}.$num.$varid.ppm"
      fi
    done
  done

  outfilename=$outbasename.$num.$varid.ppm
  echo "Cycle $c:"
  echo "  Stitching $nproc processors together to form $outfilename"
#  echo $arglist
  montage +frame +shadow +label -geometry +0+0 -tile $layout $arglist $outfilename

  framelist="$framelist $outfilename"
done

convert -loop 0 -delay 10 -scale 400% $framelist $outbasename.gif
animate $outbasename.gif &
