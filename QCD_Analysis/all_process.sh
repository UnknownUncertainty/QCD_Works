#!/bin/bash

folder1=/home/fero/Desktop/QCD_Proceses/Data_Roots
folder2=/home/fero/Desktop/Programs/qcd_delphes/delphes_analysis
				
		newrootfilename=(The_QCD_100-250 The_QCD_250-500 The_QCD_500-1000 The_QCD_1000-2500 The_QCD_2500-4000 The_QCD_4000-6000 The_QCD_6000-Inf)
		rootfilename=(combined_100-250.root combined_250-500.root combined_500-1000.root combined_1000-2500.root combined_2500-4000.root combined_4000-6000.root combined_6000-Inf.root)
		
		#echo ${rootfilename[1]}

		for((i=0; i < ${#rootfilename[@]}; i++))
		do 
		cd $folder2/
		sed -i "s|"${folder1}/.*"|"${folder1}/${rootfilename[$i]}"|g" $folder2/list_myfiles
		sed -i 's#OutputFileName = "'.*'"#OutputFileName = "'${newrootfilename[$i]}'"#g' $folder2/input_ttjets_lept.par
		make
		./skimming_eventsv2 input_ttjets_lept.par
		echo ${rootfilename[$i]}
		done
		cd $folder2/
		root -l << EOF
	        .x hadd.C		
		.q
		EOF

		date

		echo "============================================================================"
                echo " All events generated. "
		echo "============================================================================"

				
exit 0
