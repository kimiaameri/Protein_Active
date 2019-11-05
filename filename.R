argv <- commandArgs(trailingOnly = TRUE)
InputPath <- argv[1]
OutputFile <- argv[2]
filename <- list.files(InputPath )
write.csv(filename , OutputFile, row.names = F, quote=FALSE)
