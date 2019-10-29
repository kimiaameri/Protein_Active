import csv
import sys

if len(sys.argv) < 4:
    sys.stderr.write('No Input CSV file and genomebed\n')
    sys.exit(0)
    
inputFile = sys.argv[1]
genomeBedpath = sys.argv[2]
outputPath = sys.argv[3]

outputFile = "varpos.sh"
with open(outputFile,'w') as outFile:
    with open(inputFile) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter=',')
        for row in csv_reader:
                outFile.write(f'reference_Genome <- as.matrix(read.table(paste0({GenomebedPath},"/nctc8325.bed"),header=F,sep="\t",stringsAsFactors = F))\n')
                outFile.write('length.genome<- nrow(reference_Genome) \n')
                outFile.write(f'intersections <- read.table(paste({intersectionspath},{row[0]}.bed,sep=""),header=F,sep="\t",stringsAsFactors = F)\n')
                outFile.write('length.intersection= nrow(intersections)\n')
                outFile.write('variant.matrix<- matrix(NA, ncol=5, nrow=length.intersection)\n')
                outFile.write('colnames(variant.matrix)<- c("Gene.Id","Variant.start","Variant.end","Gene.length","Chromosome.Length")\n')
                outFile.write('for (k in 1:length.intersection) \n')
                outFile.write('for (j in 1 :length.genome)\n')
                  outFile.write('if (intersections[k,2] >= reference_Genome[j,2] & intersections[k,2] <= reference_Genome[j,3]) \n')
                  outFile.write(' { \n')
                     outFile.write('variant.matrix[k,1] = as.character(reference_Genome [j,4])\n')
                     outFile.write('z<-round(abs(as.numeric(intersections[k,2]) - as.numeric(reference_Genome [j,2]))/3)\n')
                     outFile.write('variant.matrix[k,2] = z\n')
                     outFile.write('variant.matrix[k,3] = z+1\n')
                     outFile.write('m<-as.numeric(reference_Genome [j,3])-as.numeric(reference_Genome [j,2])\n')
                     outFile.write('variant.matrix[k,4] = m\n')
                     outFile.write('variant.matrix[k,5] = round(m/3)\n')
                    outFile.write(' }\n')
               outFile.write('variant.matrix<- variant.matrix[complete.cases(variant.matrix),]\n')')
               outFile.write( f'write.csv(x=variant.matrix,file = paste(paste0({outputPath},"/VariantPosition/"),paste0(row[0],".csv"),sep=""), row.names = FALSE)\n')

                
