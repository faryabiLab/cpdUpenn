for i in `awk -F "\t" '{print $1}' Primers/Manifest.txt`
do
wc -l $i
done

