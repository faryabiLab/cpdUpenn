##V2.1.1
##Run this script as ./WaitTillComplete.sh <Path of RunFolder on MiseqAnalysis>
##Changed in V2.0: 1) The script now greps the line 'Ending Execution for Copier' from AnalysisLog.txt to wait until the fastqs and other important items get copied to the RunFolder on MiSeqOutput.
##                 2) Script now detects the creation of RunFolder on MiSeqAnalysis drive every 4 hours. It then detects the creation of AnalysisLog.txt file (on the same drive) every 1 hour. The echo statements are changed accordingly.
## Changes in V2.1.1: 1) Changed the echo line under the 'elif [ $grepOut -eq 1 ]' condition in order to reflect the correct statement i.e. corresponding to the new workflow adopted by Pipe3_new.sh.

#!/bin/sh


name=$1
while :
do
	if [ -d $name ];
	then
		echo "Run Folder created\nSearching for AnalysisLog.txt..."

		while :
		do
			if [ -r $name/AnalysisLog.txt ];
			then
				echo "AnalysisLog.txt created\nWaiting for analysis and copy to be completed..."
				while :
				do
					grepOut=`grep 'Ending Execution for Copier' $name/AnalysisLog.txt |wc -l`
	
					if [ $grepOut -eq 0 ];
					then
						echo "Analysis and copy not completed yet, waiting for 10 more minutes..."
						sleep 10m
					elif [ $grepOut -eq 1 ];
					then
						echo "Analysis and Copying Completed"
						exit 0 
					fi
				done
		
			else
				echo "AnalysisLog.txt is not created yet, waiting for 1 more hour.."
				sleep 1h
			fi
		done
		

	else
		echo "Run Folder not created yet\nWaiting for 4 hours..."
		sleep 4h
	fi

done




