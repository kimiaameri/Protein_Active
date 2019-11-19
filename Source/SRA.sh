for x in `cat SrA.Accession.txt`; do 
fastq-dump --split-files $x ;  
sleep 5;
done
