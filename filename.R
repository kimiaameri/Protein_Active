argv <- commandArgs(trailingOnly = TRUE)
InputPath <- argv[1]
OutputFile <- argv[2]
filename <- list.files(InputPath )
write.table(filename , OutputFile,  col.names = F, row.names = F , quote=FALSE)
