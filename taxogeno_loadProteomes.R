source("taxogeno_library.R")
library(DBI)
library(tools)
dbConn <- dbConnect(RPostgres::Postgres(),dbname="taxogeno")

dirNameVec<-c("archaeas","bacteria","bacteria2","mammals","fungi","invert","protozoa","vert_other")
#dirNameVec<-c("plants")

for(dirName in dirNameVec){
dirNamePath<-sprintf("/opt/taxogeno-src/fuentes_internas/fuentes_crudas/%s/",dirName)

baseNameVec<-gsub("^(.*)_uniref90_go_goslim[.]tsv$","\\1",list.files(path=dirNamePath,pattern=".*_uniref90_go_goslim[.]tsv"))
tsvVec<-paste(baseNameVec,"_uniref90_go_goslim.tsv",sep="")
faaVec<-paste(baseNameVec,".faa",sep="")

filesDf<-data.frame(tsv=tsvVec,faa=faaVec)

timeTakenVec<-c()
geneCountVec<-c()

for(rowNum in seq_len(nrow(filesDf))){
    tsvFilePath<-paste(dirNamePath,filesDf[rowNum,"tsv"],sep="")
    faaFilePath<-paste(dirNamePath,filesDf[rowNum,"faa"],sep="")

    if( file.exists(tsvFilePath) &&
        file.exists(faaFilePath)    )
    {
        start.time <- Sys.time()
        saveSma3sFileSet(dbConn,
                         tsvFilePath,
                         faaFilePath,
                         isUserProteome=FALSE)
        end.time <- Sys.time()
        time.taken <- end.time - start.time
        timeTakenVec<-c(timeTakenVec,time.taken)
        print(mean(timeTakenVec))
    }
}
}
