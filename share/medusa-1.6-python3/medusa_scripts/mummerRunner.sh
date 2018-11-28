file1=$1
file2=$2
threads=$3

fname1=`echo $file1 | awk -F "/" '{print $NF}'`
fname2=`echo $file2 | awk -F "/" '{print $NF}'`
tag1=`echo $fname1 | awk -F "." '{print $1}'`
tag2=`echo $fname2 | awk -F "." '{print $1}'`

if [[ -e $threads ]]
then
	thr="--threads=$threads"
fi

prefix=`echo $tag1"_"$tag2`

nucmer $thr --prefix=$prefix $file1 $file2 
show-coords -rcl $prefix.delta > $prefix.coords
echo $prefix.coords
