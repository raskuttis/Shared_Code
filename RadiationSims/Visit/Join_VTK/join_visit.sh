#!/bin/bash

usage () {
  echo "Usage:  join_visit.sh -i <INBASENAME> -v <VAR_ID> -c <START:SKIP:END>"
  echo "  e.g., join_visit.sh -i RadParGrav -v starpar -c 0:5:100"
  exit 1
}

# Specify the location of the Athena installation directory
homedir="/scratch/gpfs/raskutti/RadParGrav/UV_M5.0e4_R15.0_N256_Tf4_All/"

# Check to see if at least one argument was specified
if [ $# -lt 1 ] ; then
   echo "You must specify at least 1 argument!"
   usage
fi

# Set defaults 
inbasename=""
varid=""
cycle="0:0"

# Parse the arguments
while getopts hi:o:v:c: opt
do
   case "$opt" in
      h) usage;;
      i) inbasename=$OPTARG;;
      v) varid=$OPTARG;;
      c) cycle=$OPTARG;;
      \?) usage;;
   esac
done

# Parse in/out basenames and variable ID
if [[ $inbasename == "" ]]; then
  echo "Please specify an input basename!"
  usage
fi

# Parse the cycle variable
if [[ ${cycle##*:*} != "" ]]; then
  usage
else
  s=${cycle%%:*}
  e=${cycle##*:}
  k=${cycle#*:}
  if [[ ${k##*:*} != "" ]]; then
    k="1"
  else
    k=${k%:*}
    if [[ ${k##*:*} == "" ]]; then
      usage
    fi
  fi
fi
echo $s $e $k

if [[ $varid != "" ]]; then
  varid=${varid}
fi
#filepattern=${inbasename}*.????${varid}.vtk
filepattern=${inbasename}*.${varid}.????.vtk
echo "filepattern = $filepattern"
filename=${inbasename}.${varid}.visit

declare -a proclist
proclist=( `find id* -name "$filepattern" | cut -d '/' -f 1 | cut -d 'd' -f 2 | sort -n | uniq` )
nproc=${#proclist[@]}
echo "nproc = $nproc"

echo "!NBLOCKS $nproc" > $filename
for c in `seq $s $k $e`
do
  num=$(printf "%04d" "$c")
  for p in ${proclist[@]}
  do
    if [[ $p == 0 ]]; then
      echo id${p}/${inbasename}.${varid}.${num}.vtk >> $filename
    else
      echo id${p}/${inbasename}-id${p}.${varid}.${num}.vtk >> $filename
    fi
  done
done
#find * -name "$filepattern" | sort -t '.' -k 2 >> $filename
cat $filename

