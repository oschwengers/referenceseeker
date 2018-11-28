dir=$1
reference_genome=$2
directory=$3
threads=$4

inps=`ls $dir | grep -v $reference_genome`
for i in $inps; do
	mmr=`bash $directory/mummerRunner.sh $reference_genome $dir/$i $threads`
	echo $mmr
	done
