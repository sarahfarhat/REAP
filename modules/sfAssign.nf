#!/usr/bin/env nextflow

process sfAssign {
    tag "$species"
    
    input:
    tuple path(blast_file), path(transcripts_fastafile), val(species)
    val(sep)

    publishDir '02_repeatsAssignmentResults', mode: 'copy'

    output:
    tuple path("${species}.filtered.assignation.id"), path("${species}.norepeats.fa"), val(species)

    script:
    """
    #1: Get the 10 first hits of the blast (keeping only best match per couple of hits). If less than 5 hits remove lines.
    if [ -s ${blast_file} ];
    then
    	awk '{print \$1}' ${blast_file} | sort -u > ${species}.repeats.blastout6.id
    	for i in `cat ${species}.repeats.blastout6.id`; do 
        	grep \$i ${blast_file} | sort -k12,12nr | awk '{if(id[\$2]!=1){print \$0; id[\$2]=1}}'; 
    	done | awk '{if(id==\$1){if(count<10){print \$0; count=count+1}}else{print \$0; id=\$1; count=1}}' >> ${species}.repeats.blastout6.top10
    	awk '{print \$1}' ${species}.repeats.blastout6.top10 | sort | uniq -c | awk '\$1<5 {print \$2}' > ${species}.to_remove
    	grep -v -f ${species}.to_remove ${species}.repeats.blastout6.top10 > ${species}.repeats.blastout6.top10.2; mv ${species}.repeats.blastout6.top10.2 ${species}.repeats.blastout6.top10

    	#2: Create a file with all possible assignation
    	awk '{split(\$2,t,"${sep}"); print \$1,t[2]}' ${species}.repeats.blastout6.top10 | sort -u > ${species}.filtered.assignation.id
    	awk '{print \$1}' ${species}.filtered.assignation.id | sort -u > ${species}.filtered.assignation.uniq.id
    	#3: Get fasta files
    	awk 'BEGIN{while(getline < "${species}.filtered.assignation.uniq.id" > 0){id[">"\$1]=1}} {if(\$0~/^>/){if(id[\$1]){t=1}else{print \$0; t=0}}else{if(t==0){print \$0}}}' "${transcripts_fastafile}" | sed -e '/^>/! s/\\(.*\\)/\\U\\1/; s/>_R_/>/' > "${species}.norepeats.fa"
    else
        touch ${species}.filtered.assignation.id ${species}.norepeats.fa
    fi
    """
}
