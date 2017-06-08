## Run this script as ./garbPick.sh <Sample.R1.btrim.summary.txt> <Sample.R2.btrim.summary.txt> <Primers.primParts File> <PrimersRC.primParts File> <PrimerManifestFile>
#!/bin/sh


R1_btrimSumFile=$1
R2_btrimSumFile=$2
PrimFile=$3
Prim_RC_File=$4
PrimManifest=$5


NameOfSample=`echo $R1_btrimSumFile |cut -d "_" -f 1`

grep -v 'no_5_found' $R2_btrimSumFile | awk '{print $1,$8}' | sort -k 1,1 > $R2_btrimSumFile.awked.sorted
grep -v 'no_5_found' $R1_btrimSumFile | awk '{print $1,$8}' | sort -k 1,1 > $R1_btrimSumFile.awked.sorted

join -1 1 -2 1 $R1_btrimSumFile.awked.sorted $R2_btrimSumFile.awked.sorted > $NameOfSample.btrim.summary.joined

grep -v '^@' $NameOfSample.novo.sam| awk '{print $1,$2}'| sort -k 1,1 > $NameOfSample.novo.sam.nameFlag

grep -v '^@' $NameOfSample.novo.sam| awk '{print $2}'| sort | uniq > $NameOfSample.novo.sam.flagValues

join -1 1 -2 1 $NameOfSample.btrim.summary.joined $NameOfSample.novo.sam.nameFlag > $NameOfSample.flagJoin

python garbPick.py $PrimFile $Prim_RC_File $NameOfSample.flagJoin $NameOfSample.garbageStats $NameOfSample.novo.sam.flagValues $PrimManifest














##join -v1 -1 1 -2 1 CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021_R1.btrim.summary.awked.sorted.txt CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021_R2.btrim.summary.awked.sorted.txt| awk '{print $1,$2" ""-"}' > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.onlyR1.txt
##join -v2 -1 1 -2 1 CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021_R1.btrim.summary.awked.sorted.txt CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021_R2.btrim.summary.awked.sorted.txt| awk '{print $1," ""-",$2}' > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.onlyR2.txt



##(cat CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.summary.joined; cat CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.onlyR1.txt;  cat CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.onlyR2.txt) > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.summary.joined.merged
##sort -k 1,1 CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.summary.joined.merged > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.summary.joined.merged.sorted

##grep -v '^@' CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.original.novo.sam| awk '{print $1,$2}'| sort -k 1,1 > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.original.novo.sam.flags

##join -1 1 -2 1 CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.btrim.summary.joined.merged.sorted CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.original.novo.sam.flags > CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.flagJoin

##python grepExtract.py PrimersFFPE185bp.preOutput primers_RC_FFPE185bp.txt CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.flagJoin CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.garbageStats CPDC130393-TSCA-13015-A704-A502-PAL-13013-SEQ-130021.original.novo.sam.flagValues

