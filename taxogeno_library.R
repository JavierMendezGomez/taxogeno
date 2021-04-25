library(uuid)
library(tools)
####################################
##  readSma3sAnnotation           ##
######################################################################################
readSma3sAnnotation<-function(sma3sAnnotationPath){
    sma3sAnnotationDf<-read.csv(sma3sAnnotationPath,
                                col.names=c("fastaheader",
                                            "genename",
                                            "genedescription",
                                            "enzyme",
                                            "goc",
                                            "gocname",
                                            "gof",
                                            "gofname",
                                            "gop",
                                            "gopname",
                                            "keyword",
                                            "pathway",
                                            "goslim"),
                                colClasses=rep("character",13),
                                quote="",
                                header=FALSE,
                                sep="\t",
                                comment.char="#",
                                stringsAsFactors = FALSE)
    fastaShortHeaderVec<-gsub("^([^ ]+) .*$","\\1",sma3sAnnotationDf[,"fastaheader"])
    sma3sAnnotationDf[,"fastashortheader"]<-fastaShortHeaderVec ; rm(fastaShortHeaderVec)
    sma3sAnnotationDf
}

####################################
##  readMultifasta                ##
######################################################################################
readMultifasta<-function(multifastaPath){   
    ## Leer fichero línea a línea para no petar cargando el fichero completo en RAM
    ## https://stackoverflow.com/questions/12626637/read-a-text-file-in-r-line-by-line#35761217

    ##inicializar variables que almacenan FASTA uno a uno
    fastaHeader=""
    fastaSequence=""

    ## booleanos control de flujo
    startSequenceBody=FALSE
    insideSequenceBody=FALSE

    fastaHeaderVec<-character()
    fastaSequenceVec<-character()
    fileConn = file(multifastaPath, "r")
    while ( TRUE ) {
        oneLine = readLines(fileConn, n = 1)
        if ( length(oneLine) == 0 ) {
            if(insideSequenceBody){
                ##cerrando cosa anterior
                fastaHeaderVec<-c(fastaHeaderVec,fastaHeader)
                fastaSequenceVec<-c(fastaSequenceVec,fastaSequence)
                
                ##limpiar variables de último elemento
                fastaHeader=""
                fastaSequence=""

                ##indicar que ya no se está leyendo el cuerpo de la secuencia
                insideSequenceBody = FALSE
            }
            break
        }

        ## Cosas de FASTA
        if(startsWith(oneLine, ">")){
            ## Encontrar una cabecera FASTA tiene que interrumpir la secuencia
            if(insideSequenceBody){
                ##cerrando cosa anterior

                fastaHeaderVec<-c(fastaHeaderVec,fastaHeader)
                fastaSequenceVec<-c(fastaSequenceVec,fastaSequence)
                
                ##limpiar variables de último elemento
                fastaHeader=""
                fastaSequence=""

                ##indicar que ya no se está leyendo el cuerpo de la secuencia
                insideSequenceBody = FALSE
            }

            ## empezando cosa nueva
            
            ## Guardar id 
            fastaHeader = gsub("^>([^\\s]*)", "\\1", oneLine)
            ##En la siguiente vuelta debe comenzar el cuerpo de la secuencia
            startSequenceBody = TRUE
            next
        }
        if(startSequenceBody){
            ##Cambio de booleano
            startSequenceBody=FALSE
            insideSequenceBody=TRUE
        }
        if(insideSequenceBody){
            ##Concatenar secuencia, por si esta está en varias líneas
            fastaSequence<-paste(fastaSequence, oneLine, sep='')
        }
    }
    close(fileConn)
    data.frame(fastaheader=fastaHeaderVec,aasequence=fastaSequenceVec)
}
###################################################
##  extractColumnDfFromInsertedGeneAnnotationDf  ##
######################################################################################
extractColumnDfFromInsertedGeneAnnotationDf<-function(insertedGeneAnnotationDf,colName,newColName=colName){
    if(length(colName)>1){
        stop("length(colName)>1")
    }
    if(length(newColName)>1){
        stop("length(newColName)>1")
    }
    colVecList<-strsplit(insertedGeneAnnotationDf[,colName],";")

    geneidVec<-rep(insertedGeneAnnotationDf[,"geneid"],lengths(colVecList))
    colVec<-unlist(colVecList) ; rm(colVecList)
    outputDf<-data.frame(geneidVec,colVec,stringsAsFactors=FALSE) ; rm(geneidVec,colVec)
    colnames(outputDf)<-c("geneid",newColName)
    
    outputDf
}
#################################################
##  ncbiAssemblyInfoFromGcaIdContainingString  ##
######################################################################################
ncbiAssemblyInfoFromGcaIdContainingString<-function(dbConn,gcaIdContainingString){
    ## extractedGcaId<-gsub('^.*(GCA_[^.]*)[.].*$','\\1',gcaIdContainingString)
    ## assemblyInfoDf<-dbGetQuery(
    ##     dbConn,
    ##     "SELECT
    ##       *
    ##      FROM ncbi.assembly_summary_genbank
    ##      WHERE assembly_accession LIKE $1||'%'",
    ##     params=list(extractedGcaId))

    extractedGcaId<-gsub('^.*(GCA_[0123456789]+[.][0123456789]+).*$','\\1',gcaIdContainingString)
    print("Extracted gcaid")
    print(extractedGcaId)
    assemblyInfoList<-dbGetQuery(
        dbConn,
        "SELECT
           assembly_accession,
           bioproject,
           biosample,
           wgs_master,
           refseq_category,
           taxid,
           species_taxid,
           organism_name,
           infraspecific_name,
           isolate,
           version_status,
           assembly_level,
           release_type,
           genome_rep,
           seq_rel_date date,
           asm_name,
           submitter,
           gbrs_paired_asm,
           paired_asm_comp,
           ftp_path,
           excluded_from_refseq,
           relation_to_type_material
     FROM ncbi.assembly_summary_genbank
     WHERE assembly_accession IN($1)",
     params=list(extractedGcaId))[1,]
    print("AssemblyInfoList")
    print(assemblyInfoList)
    assemblyInfoList
}
#################################################
##  updateNcbiAssemblyInfoForProteomeIdVec     ##
######################################################################################
updateNcbiAssemblyInfoForProteomeIdVec<-function(dbConn,proteomeIdVec){
    ## sacar la información del/de los proteomas de la base de datos
    proteomeInfoDf<-dbGetQuery(dbConn,"SELECT * FROM taxogeno.proteome WHERE proteomeid IN ($1)", params=list(proteomeIdVec))

    for(rowNum in seq_along(nrow(proteomeInfoDf))){
        gcaIdContainingString=NA
        ## si se proporcionó gcaid usarlo primero
        if(!is.na(proteomeInfoDf[rowNum,][["gcaid"]])){
            gcaIdContainingString<-proteomeInfoDf[rowNum,"gcaid"]
        }
        ## si no hay gcaid, tirar del nombre del fichero
        else if(!is.na(proteomeInfoDf[rowNum,"basename"])){
            gcaIdContainingString<-proteomeInfoDf[rowNum,"basename"]
        }
        
        ## tomar la información nueva
        ncbiAssemblyInfoList<-ncbiAssemblyInfoFromGcaIdContainingString(dbConn, gcaIdContainingString)
        str(ncbiAssemblyInfoList)
        if(length(ncbiAssemblyInfoList)>0){
            print("Entra en lista de información ncbi")
            ## ^############
            ## comprobar si había un taxón antiguo almacenado
            oldNcbiTaxId<-dbGetQuery(dbConn,
                                     "SELECT ncbitaxid FROM taxogeno.taxonomy WHERE ncbitaxid = $1 AND is_ancestor=FALSE",
                                     params=list(proteomeInfoDf[rowNum,"ncbitaxid"]) ) [1,"ncbitaxid"]

            if(!is.na(oldNcbiTaxId) && oldNcbiTaxId != ncbiAssemblyInfoList[["taxid"]]){
                dbExecute(dbConn, "CALL taxogeno.delete_taxon($1)", params=list(oldNcbiTaxId))
            }
            ## $##########

            ## ^############
            ## comprobar si había un taxón antiguo almacenado igual que el que se quiere insertar
            oldNcbiTaxId<-dbGetQuery(dbConn,
                                     "SELECT ncbitaxid FROM taxogeno.taxonomy WHERE ncbitaxid = $1 AND is_ancestor=FALSE",
                                     params=list(ncbiAssemblyInfoList[["taxid"]]) ) [1,"ncbitaxid"]
            
            if(is.na(oldNcbiTaxId)){
                dbExecute(dbConn, "CALL taxogeno.insert_taxon($1)", params=list(ncbiAssemblyInfoList[["taxid"]]))
            }
            ## $############

            dbExecute(dbConn,
                      "UPDATE taxogeno.proteome SET (gcaid,ncbitaxid,sourceurl,dbname) = ($2,$3,$4,$5) WHERE proteomeid=$1",
                      params=list(proteomeInfoDf[rowNum,"proteomeid"],
                                  ncbiAssemblyInfoList[["assembly_accession"]],
                                  ncbiAssemblyInfoList[["taxid"]],
                                  ncbiAssemblyInfoList[["ftp_path"]],
                                  "NCBI:GenBank") )
        }
    }
}
#####################
##  generateGafDf  ##
######################################################################################
generateGafDf<- function(dbConn,ontologyAnnotDf) {
    gafDf<-data.frame(
        "db"                 =rep("taxogeno", nrow(ontologyAnnotDf)),
        "dbobjectid"         =ontologyAnnotDf[,"geneid"],
        "dbobjectsymbol"     =ontologyAnnotDf[,"geneid"],
        "qualifier"          =rep(NA,nrow(ontologyAnnotDf)),
        "annotkwid"         =ontologyAnnotDf[,"annotkwid"],
        "dbreference"        =rep(NA,nrow(ontologyAnnotDf)),
        "evidencecode"       =rep("IEA",nrow(ontologyAnnotDf)),
        "with"               =rep(NA,nrow(ontologyAnnotDf)),
        "aspect"             =ontologyAnnotDf[,"aspect"],
        "dbobjectname"       =rep(NA,nrow(ontologyAnnotDf)),
        "dbobjectsynonym"    =rep(NA,nrow(ontologyAnnotDf)),
        "dbobjecttype"       =rep("protein",nrow(ontologyAnnotDf)),
        "taxon"              =rep(NA,nrow(ontologyAnnotDf)),
        "date"               =rep( strftime( as.POSIXlt(Sys.time()), "%Y-%m-%dT%H:%M:%S%z"), nrow(ontologyAnnotDf) ),
        "assignedby"         =rep("sma3s",nrow(ontologyAnnotDf)),
        "annotationextension"=rep(NA,nrow(ontologyAnnotDf)),
        "geneproductformid"  =rep(NA,nrow(ontologyAnnotDf))
    )
    gafDf
}

###########################
## owltoolsMap2SlimGafDf ##
######################################################################################
owltoolsMap2SlimGafDf<-function(gafDf,owl,subset=NA){
    gafNamesVec<-c("db","dbobjectid","dbobjectsymbol","qualifier",
                   "annotkwid","dbreference","evidencecode","with",
                   "aspect","dbobjectname","bobjectsynonym","dbobjecttype",
                   "proteomeid","date","assignedby","annotationextension",
                   "geneproductformid")
    inputGafFilePath<- tempfile()
    write.table(gafDf, file = inputGafFilePath, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "", dec = ".", row.names = FALSE,
                col.names = FALSE)
    
    outputGafFilePath<- tempfile()

    owltoolsArgs<-if(!is.na(subset))
                       c(owl,
                        "--gaf",inputGafFilePath,
                        "--map2slim",
                        "--subset",subset,
                        "--write-gaf",outputGafFilePath)
                   else
                       c(owl,
                        "--gaf",inputGafFilePath,
                        "--map2slim",
                        "--write-gaf",outputGafFilePath)
    
    map2slimGafStdOut<-system2("./bin/owltools",args=owltoolsArgs,stdout=TRUE)
    #unlink(inputGafFilePath)
    map2slimGafDf<-read.csv(text=map2slimGafStdOut,header=FALSE,quote="",stringsAsFactors=FALSE,comment.char='!',sep="\t")
    outputGafDf<-read.csv(outputGafFilePath,header=FALSE,quote="",stringsAsFactors=FALSE, col.names=gafNamesVec,comment.char='!',sep="\t")
    #unlink(outputGafFilePath)
    outputGafDf
}
############################
## convertGafToSummaryDf  ##
######################################################################################
convertGafToSummaryDf<-function(gafDf,proteomeId){
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }
    summaryDf<-aggregate(dbobjectid~annotkwid,data=gafDf,length)
    summaryDf[,"proteomeid"]<-rep(proteomeId,nrow(summaryDf))
    colnames(summaryDf)[colnames(summaryDf)=="dbobjectid"]<-"annotkwcount"   
    summaryDf
}
###################################################
## saveGeneratedGoslimSummaryForGoAnnotationData ##
######################################################################################
saveGeneratedGoslimSummaryForGoStrictAnnotationData<-function(dbConn,proteomeId,goStrictDf){
        gafDf<-generateGafDf(dbConn,goStrictDf)
        mappedGafDf<-owltoolsMap2SlimGafDf(gafDf,owl="goslim_generic.owl",subset="goslim_generic")
        summaryDf<-convertGafToSummaryDf(mappedGafDf, proteomeId)
        dbWriteTable(dbConn,SQL("taxogeno.generated_goslim_summary"), summaryDf,append=TRUE,row.names=FALSE)
}
## ###################################################
## ## saveGeneratedKeywordSummaryForGoAnnotationData ##
## ######################################################################################
## saveGeneratedKeywordSummaryForGoStrictAnnotationData<-function(dbConn,proteomeId,keywordStrictDf){
##         gafDf<-generateGafDf(dbConn,keywordStrictDf)
##         mappedGafDf<-owltoolsMap2SlimGafDf(gafDf,owl="keyword.obo")
##         summaryDf<-convertGafToSummaryDf(mappedGafDf, proteomeId)
##         dbWriteTable(dbConn,SQL("taxogeno.generated_keyword_summary"), summaryDf,append=TRUE,row.names=FALSE)
## }

####################################
##                                ##
######################################################################################
chunkVectorInEqualSizeFragments <- function(x,n) {
        split(x, ceiling(seq_along(x)/n))
}
####################################################
## getNormalizedAnnotationForProteomeIdVec ##
######################################################################################
getNormalizedAnnotationForProteomeIdVec<-function(dbConn,annotationType,proteomeIdVec){
    if(length(annotationType)>1){
        stop("length(annotationType)>1")
    }

    normalizedAnnotationDf<-data.frame(proteomeid=character(),annotkwid=character(),annotkwrelval=numeric())
    if(annotationType=="generated_goslim"){
        normalizedAnnotationDf<-dbGetQuery(
           dbConn,
           "SELECT
              taxogeno.proteome.proteomeid as proteomeid,
              annotkwid                     as annotkwid,
              annotkwcount                  as annotkwcount,
              annotkwcount::float/genecount as annotkwrelval

            FROM taxogeno.generated_goslim_summary
            INNER JOIN taxogeno.proteome
              ON taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid
              AND taxogeno.generated_goslim_summary.proteomeid IN ($1)

            ORDER BY taxogeno.proteome.proteomeid, annotkwid",
           params=list(proteomeIdVec)
        )
    } else if(annotationType=="keyword"){
        proteomeIdKeywordIdCartesianDf<-dbGetQuery(
            dbConn,
            "select proteomeid as proteomeid,
                 genecount as genecount,
                   keyword as keyword,
                 keywordid as annotkwid
           from taxogeno.proteome,
                uniprot.uniprot_keyword
          where proteomeid in ($1)",
          params=list(proteomeIdVec))

        proteomeIdKeywordCountDf<-dbGetQuery(
            dbConn,
            "select tg.proteomeid      as proteomeid,
                tgk.keyword        as keyword,
                count(tgk.keyword) as annotkwcount
         from taxogeno.gene tg
         inner join taxogeno.gene_keyword_rel tgk
	 on tg.geneid=tgk.geneid and tg.proteomeid in ($1)
         where tg.proteomeid in ($1)
         group by tg.proteomeid,tgk.keyword",
         params=list(proteomeIdVec))

        summaryDf<-merge(proteomeIdKeywordIdCartesianDf,proteomeIdKeywordCountDf,by=c("proteomeid","keyword"),all.x=TRUE)
        summaryDf[,"annotkwrelval"]<-summaryDf[,"annotkwcount"]/summaryDf[,"genecount"]
        
        normalizedAnnotationDf<-summaryDf[,c("proteomeid","annotkwid","annotkwcount","annotkwrelval")]
        rm(summaryDf)
    }
    normalizedAnnotationDf
}
####################################################
## getKeywordNormalizedAnnotationForProteomeIdVec ##
######################################################################################
getKeywordNormalizedAnnotationForProteomeIdVec<-function(dbConn,proteomeIdVec){
    proteomeIdKeywordIdCartesianDf<-dbGetQuery(
        dbConn,
        "select proteomeid as proteomeid,
                 genecount as genecount,
                   keyword as keyword,
                 keywordid as annotkwid
           from taxogeno.proteome,
                uniprot.uniprot_keyword
          where proteomeid in ($1)",
        params=list(proteomeIdVec))

    proteomeIdkeywordCountDf<-dbGetQuery(
        dbConn,
        "select tg.proteomeid      as proteomeid,
                tgk.keyword        as keyword,
                count(tgk.keyword) as annotkwcount
         from taxogeno.gene tg
         inner join taxogeno.gene_keyword_rel tgk
	 on tg.geneid=tgk.geneid and tg.proteomeid in ($1)
         where tg.proteomeid in ($1)
         group by tg.proteomeid,tgk.keyword",
        params=list(proteomeIdVec))

    summaryDf<-merge(proteomeIdKeywordIdCartesianDf,proteomeIdkeywordCountDf,by=c("proteomeid","keyword"),all.x=TRUE)
    summaryDf[,"annotkwrelval"]<-summaryDf[,"annotkwcount"]/summaryDf[,"genecount"]
    summaryDf[,c("proteomeid","annotkwid","annotkwrelval")]
}
####################################
##                                ##
######################################################################################
getGeneratedGoslimNormalizedAnnotationForProteomeIdVec<-function(dbConn,proteomeIdVec){
    dbGetQuery(dbConn,
               "SELECT
                  taxogeno.proteome.proteomeid as proteomeid,
                  annotkwid                     as annotkwid,
                  annotkwcount                  as annotkwcount,
                  annotkwcount::float/genecount as annotkwrelval

                FROM taxogeno.generated_goslim_summary
                INNER JOIN taxogeno.proteome
                  ON taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid
                  AND taxogeno.generated_goslim_summary.proteomeid IN ($1)

                ORDER BY taxogeno.proteome.proteomeid, annotkwid", params=list(proteomeIdVec))
}
####################################
##                                ##
######################################################################################
distAsDf <- function(inDist) {
    ## https://stackoverflow.com/questions/23474729/convert-object-of-class-dist-into-data-frame-in-r
    if (class(inDist) != "dist") stop("wrong input type")
    A <- attr(inDist, "Size")
    B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
    if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
    if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
    data.frame(
        row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
        col = rep(B[-length(B)], (length(B)-1):1),
        value = as.vector(inDist)
    ) ## end data.frame
}
##################################
## calculateEuclideanDistances  ##
######################################################################################
calculateEuclideanDistances<-function (normalizedAnnotationDf){   
    euclideanDistancesDf<-distAsDf(
        dist(
            ## ####################
            ##            kw
            ## proteomeid a b c d
            ##          1 1 . . .
            ##          2 . 3 . .
            ##          3 . . 5 .
            ##          4 . . . 7
            xtabs( annotkwrelval~proteomeid+annotkwid, normalizedAnnotationDf, sparse=TRUE),
            ## De esa matriz dispersa la distancia euclídea entre las filas
            method="euclidean",
            diag=FALSE,
            upper=FALSE
        ) ## end dist
    ) ## end distAsDf
    colnames(euclideanDistancesDf)<-c("proteomeid_greatest","proteomeid_least","distance")
    euclideanDistancesDf
}

################################################
## saveEuclideanDistancesForProteomeId ##
######################################################################################
saveEuclideanDistancesForProteomeId<-function(dbConn,proteomeId,annotationTypeVec){
    ## integrity tests
    allowedAnnotationTypes <- c("generated_goslim","keyword")
    if(!(all(annotationTypeVec %in% allowedAnnotationTypes))){
        stop("!(annotationTypeVec %in% allowedAnnotationTypes)")
    }
    ## only one proteome in proteomeId vector
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }

    ## only proceed if proteomeId exists in databse
    existsProteomeId<-as.logical( dbGetQuery(dbConn,"SELECT EXISTS( SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid=$1 )",
                                             params=list(proteomeId))[1,1] )
    if(!existsProteomeId){
        stop("!existsProteomeId")
    }
    
    ## ## ## ## ##    ## ## ## ## ##    ## ## ## ## ##    ## ## ## ## ##    ## ## ## ## ##    ## ## ## ## ##    ## ## ## ## ##
    for(annotationType in annotationTypeVec){
        ## examining candidates for euclidean distance comparision
        otherProteomeIdVec<-dbGetQuery(dbConn,"SELECT proteomeid FROM taxogeno.proteome")[,"proteomeid",drop=TRUE]
        ## Porque no tiene sentido comparar entre sí los que ya están guardados
        ## Hay que sacar un vector con los que ya existen
        alreadyExistsVec<-dbGetQuery(dbConn,
                                     sprintf("SELECT proteomeid_least as proteomeid
                                                FROM taxogeno.euclidean_distances_%1$s
                                               WHERE proteomeid_least=$1
                                               UNION
                                              SELECT proteomeid_greatest as proteomeid
                                                FROM taxogeno.euclidean_distances_%1$s
                                               WHERE proteomeid_greatest=$1",annotationType),
                                     params=list(proteomeId))[,"proteomeid"]

        ## Y dejar únicamente los que no estén ya comparados entre sí
        otherProteomeIdVecFiltered<-setdiff(otherProteomeIdVec, alreadyExistsVec)
        
        ## Como puede petar la memoria al generar la matriz dispersa y el objeto dist, ir de 20 en 20
        for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=1)){
            ## Se añade el proteoma que se está introduciendo al trozo de ids de proteoma contra los que se va a comparar
            comparingProteomeIdVec<-c(proteomeId,chunkVec)

            ## Obtener las anotaciones normalizadas para cada uno de ellos en un dataframe en formato mapa    
            normalizedKeywordAnnotationDf<-getNormalizedAnnotationForProteomeIdVec(dbConn,annotationType,comparingProteomeIdVec)

            ## Calcular la matriz dispersa reconvertida a data frame que tiene las distancias euclidianas
            euclideanDistancesDf<-calculateEuclideanDistances(normalizedKeywordAnnotationDf)

            ## Sólo añadir las distancias calculadas para el proteoma que se está introduciendo, porque si no se viola la clave primaria
            dbWriteTable( dbConn,
                         SQL(sprintf("taxogeno.euclidean_distances_%s",annotationType)),
                         euclideanDistancesDf [euclideanDistancesDf[,"proteomeid_greatest"]==proteomeId | euclideanDistancesDf[,"proteomeid_least"]==proteomeId , ],
                         append=TRUE, row.names=FALSE )
        } ## end for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=20))
    } ## end for(annotationType in annotationTypeVec)
}

################################################
## saveKeywordEuclideanDistancesForProteomeId ##
######################################################################################
saveKeywordEuclideanDistancesForProteomeId<-function(dbConn,proteomeId){
    ## integrity tests
    ## only one proteome in proteomeId vector
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }

    ## only proceed if proteomeId exists in databse
    existsProteomeId<-as.logical( dbGetQuery(dbConn,"SELECT EXISTS( SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid=$1 )",
                                             params=list(proteomeId))[1,1] )
    if(!existsProteomeId){
        stop("!existsProteomeId")
    }


    ## examining candidates for euclidean distance comparision
    otherProteomeIdVec<-dbGetQuery(dbConn,"SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid<>$1",
                                   params=list(proteomeId))[,"proteomeid",drop=TRUE]
    ## Porque no tiene sentido comparar entre sí los que ya están guardados
    ## Hay que sacar un vector con los que ya existen
    alreadyExistsVec<-dbGetQuery(dbConn,
                                "SELECT proteomeid_least as proteomeid
                                   FROM taxogeno.euclidean_distances_keyword
                                  WHERE proteomeid_least=$1
                                  UNION
                                 SELECT proteomeid_greatest as proteomeid
                                   FROM taxogeno.euclidean_distances_keyword
                                  WHERE proteomeid_greatest=$1",
                                params=list(proteomeId))[,"proteomeid"]

    ## Y dejar únicamente los que no estén ya comparados entre sí
    otherProteomeIdVecFiltered<-setdiff(otherProteomeIdVec, alreadyExistsVec)
    
    ## Como puede petar la memoria al generar la matriz dispersa y el objeto dist, ir de 20 en 20
    for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=20)){
        ## Se añade el proteoma que se está introduciendo al trozo de ids de proteoma contra los que se va a comparar
        comparingProteomeIdVec<-c(proteomeId,chunkVec)
        ## Obtener las anotaciones normalizadas para cada uno de ellos en un dataframe en formato mapa
        normalizedKeywordAnnotationDf<-getKeywordNormalizedAnnotationForProteomeIdVec(dbConn, comparingProteomeIdVec)
        ## Calcular la matriz dispersa reconvertida a data frame que tiene las distancias euclidianas
        euclideanDistancesDf<-calculateEuclideanDistances(normalizedKeywordAnnotationDf)
        ## Sólo añadir las distancias calculadas para el proteoma que se está introduciendo, porque si no se viola la clave primaria
        dbWriteTable( dbConn,
                      SQL("taxogeno.euclidean_distances_keyword"),
                      euclideanDistancesDf [euclideanDistancesDf[,"proteomeid_greatest"]==proteomeId | euclideanDistancesDf[,"proteomeid_least"]==proteomeId , ],
                      append=TRUE, row.names=FALSE )
    }
}


########################################################
## saveGeneratedGoslimEuclideanDistancesForProteomeId ##
######################################################################################
saveGeneratedGoslimEuclideanDistancesForProteomeId<-function(dbConn,proteomeId){
    ## integrity tests
    ## only one proteome in proteomeId vector
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }

    ## only proceed if proteomeId exists in databse
    existsProteomeId<-as.logical( dbGetQuery(dbConn,"SELECT EXISTS( SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid=$1 )",
                                             params=list(proteomeId))[1,1] )
    if(!existsProteomeId){
        stop("!existsProteomeId")
    }


    ## examining candidates for euclidean distance comparision
    otherProteomeIdVec<-dbGetQuery(dbConn,"SELECT proteomeid FROM taxogeno.proteome WHERE proteomeid<>$1",
                                   params=list(proteomeId))[,"proteomeid",drop=TRUE]
    ## Porque no tiene sentido comparar entre sí los que ya están guardados
    ## Hay que sacar un vector con los que ya existen
    alreadyExistsVec<-dbGetQuery(dbConn,
                                "SELECT proteomeid_least as proteomeid
                                   FROM taxogeno.taxogeno.euclidean_distances_generated_goslim
                                  WHERE proteomeid_least=$1
                                  UNION
                                 SELECT proteomeid_greatest as proteomeid
                                   FROM taxogeno.taxogeno.euclidean_distances_generated_goslim
                                  WHERE proteomeid_greatest=$1",
                                params=list(proteomeId))[,"proteomeid"]

    ## Y dejar únicamente los que no estén ya comparados entre sí
    otherProteomeIdVecFiltered<-setdiff(otherProteomeIdVec, alreadyExistsVec)
    
    ## Como puede petar la memoria al generar la matriz dispersa y el objeto dist, ir de 20 en 20
    for(chunkVec in chunkVectorInEqualSizeFragments(otherProteomeIdVecFiltered, n=20)){
        ## Se añade el proteoma que se está introduciendo al trozo de ids de proteoma contra los que se va a comparar
        comparingProteomeIdVec<-c(proteomeId,chunkVec)
        ## Obtener las anotaciones normalizadas para cada uno de ellos en un dataframe en formato mapa
        normalizedGeneratedGoslimAnnotationDf<-getGeneratedGoslimNormalizedAnnotationForProteomeIdVec(dbConn, comparingProteomeIdVec)
        ## Calcular la matriz dispersa reconvertida a data frame que tiene las distancias euclidianas
        euclideanDistancesDf<-calculateEuclideanDistances(normalizedGeneratedGoslimAnnotationDf)
        ## Sólo añadir las distancias calculadas para el proteoma que se está introduciendo, porque si no se viola la clave primaria
        dbWriteTable( dbConn,
                     SQL("taxogeno.euclidean_distances_generated_goslim"),
                     euclideanDistancesDf [euclideanDistancesDf[,"proteomeid_greatest"]==proteomeId | euclideanDistancesDf[,"proteomeid_least"]==proteomeId , ],
                     append=TRUE, row.names=FALSE )
    }
}


## ##############################
## ## Función tocha            ##
## ######################################################################################
## ##  1. Lee las anotaciones
## ##  2. Lee el multifasta
## ##  3. Mete la información en la tabla de proteoma
## ##  4. Si es un proteoma que tiene GCA, asignarle el ncbitaxid
## ##  5. Si no es un proteoma de usuario, introducir su taxón y el de sus ancestros
## ##     en la labla de taxonomía
## ##  6. Mete la información básica de cada gen
## ##  7. A los genes les mete su secuencia sacándola del multifasta
## ##  8. Para cada gen mete cada campo de las anotaciones en cada tabla correspondiente
## ##  9. Recalcula los goslim y saca un conteo por cada combinatoria de proteína y goslim,
## ##     cosa que guarda en una tabla
## ## 10. Calcula las distancias euclídeas respecto el resto de proteomas que haya en
## ##     la base de datos partiendo de la información de los goslim recalculados
## #########################################################-#############################
## saveSma3sFileSet<-function(dbConn,
##                            sma3sAnnotationFileName,
##                            multifastaFileName,
##                            tagVec=character(),
##                            isUserProteome=FALSE,
##                            gcaId=NA,
##                            ncbiTaxId=NA,
##                            dbName=NA,
##                            sourceUrl=NA,
##                            doUpdateNcbiAssemblyInfo=TRUE){

##     sma3sAnnotationDf<-readSma3sAnnotation(sma3sAnnotationFileName)
##     sma3sAnnotationMd5Sum<- md5sum(sma3sAnnotationFileName)
##     multifastaDf<-readMultifasta(multifastaFileName)
##     multifastaMd5Sum<-md5sum(multifastaFileName)

##     geneCount<-sum(complete.cases(sma3sAnnotationDf[,"fastaheader"]), na.rm=TRUE)
##     if(geneCount==0){
##         error("Empty Sma3s tsv file. Proteome with no genes.")
##     }

##     ## ## Campos del proteoma
##     proteomeInfoObj<-data.frame(
##         ## proteomeid=NA, ## proteomeid asignado por la base de datos. No se mete como columna en el DF porque lo va a interpretar como que es null en dbWriteTable
##         basename = basename(tools::file_path_sans_ext(sma3sAnnotationFileName)),
##         sma3sannotationmd5sum = sma3sAnnotationMd5Sum,
##         multifastamd5sum = multifastaMd5Sum,
##         genecount = geneCount,
##         is_userproteome = isUserProteome,
##         gcaid = gcaId,
##         ncbitaxid = ncbiTaxId,
##         dbname = dbName,
##         sourceurl = sourceUrl,
##         do_updatencbiassemblyinfo = doUpdateNcbiAssemblyInfo
##     )
        
##     proteomeId<-dbGetQuery(dbConn,"
##       INSERT INTO taxogeno.proteome(
##         basename,
##         sma3sannotationmd5sum,
##         multifastamd5sum,
##         genecount,
##         is_userproteome,
##         gcaid,
##         ncbitaxid,
##         dbname,
##         sourceurl,
##         do_updatencbiassemblyinfo
##       ) VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10)
##       RETURNING proteomeid 
##     ",params=list(proteomeInfoObj[["basename"]],
##                   proteomeInfoObj[["sma3sannotationmd5sum"]],
##                   proteomeInfoObj[["multifastamd5sum"]],
##                   proteomeInfoObj[["genecount"]],
##                   proteomeInfoObj[["is_userproteome"]],
##                   proteomeInfoObj[["gcaid"]],
##                   proteomeInfoObj[["ncbitaxid"]],
##                   proteomeInfoObj[["dbname"]],
##                   proteomeInfoObj[["sourceurl"]],
##                   proteomeInfoObj[["do_updatencbiassemblyinfo"]])) [1,"proteomeid"] ## tomar la primera entrada de la columna proteomeid, porque sólo hay una, del RETURNING

##     ## ## si es un proteoma no de un usuario, sacar el posible GCA del proteoma que se ha introducido a partir de su nombre
##     if(!isUserProteome && doUpdateNcbiAssemblyInfo ){
##         updateNcbiAssemblyInfoForProteomeIdVec(dbConn,proteomeId)
##     }

##     ## actualizar a lo que hay en la bbdd
##     proteomeInfoObj<-dbGetQuery(dbConn,"SELECT * FROM taxogeno.proteome WHERE proteomeid=$1", params=list(proteomeId))
    
##     ## ## Etiquetas del proteoma
##     tagDf<-data.frame(proteomeid = rep(proteomeId,length(tagVec)),
##                       tag = tagVec )
##     dbWriteTable(dbConn,SQL("taxogeno.proteome_tag_rel"),tagDf,append=TRUE,row.names=FALSE)

##     ## GENES
##     ## ## Campos generales de los genes
##     geneDf<-sma3sAnnotationDf[,c("fastashortheader","genename","genedescription","fastaheader")]
##     ## ## Engancharles sus secuencias si las hay
##     geneDf<-merge(geneDf,multifastaDf,by="fastaheader",all.x=TRUE)
##     ## ## Vincular con el proteoma que se está metiendo
##     geneDf[,"proteomeid"]<-rep(proteomeId, nrow(geneDf))
##     ## ## Escribir en la base de datos
##     dbWriteTable(dbConn,SQL("taxogeno.gene"),geneDf,append=TRUE,row.names=FALSE)

##     ## ## las claves autoincrementales hay que recuperarlas de la base de datos
##     geneIdFastaShortHeaderRelDf<-dbGetQuery(dbConn,
##                                             "SELECT geneid,fastashortheader FROM taxogeno.gene WHERE proteomeid=$1",
##                                             params=list(proteomeId))
##     insertedGeneAnnotationDf<-merge(sma3sAnnotationDf,geneIdFastaShortHeaderRelDf, by="fastashortheader")
##     rm(sma3sAnnotationDf)
    
##     ## Campos 1:n gen:campos
##     enzymeDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "enzyme", newColName="ec")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_enzyme_rel"),enzymeDf,append=TRUE,row.names=FALSE)

##     keywordDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "keyword")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_keyword_rel"),keywordDf,append=TRUE,row.names=FALSE)
##     rm(keywordDf)

##     ## #########################
##     pathwayDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "pathway")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_pathway_rel"),pathwayDf,append=TRUE,row.names=FALSE)
##     rm(pathwayDf)
    
##     ## #########################
##     goslimDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "goslim", newColName="goid")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_goslim_rel"),goslimDf,append=TRUE,row.names=FALSE)
##     rm(goslimDf)
    
##     ## #########################
##     gocDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"goc", newColName="goid")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_goc_rel"),gocDf,append=TRUE,row.names=FALSE)
##     gofDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gof", newColName="goid")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_gof_rel"),gofDf,append=TRUE,row.names=FALSE)
##     gopDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gop", newColName="goid")
##     dbWriteTable(dbConn,SQL("taxogeno.gene_gop_rel"),gopDf,append=TRUE,row.names=FALSE)
##     goCategory<-c("goc"="cellular_component",
##                   "gop"="biological_process",
##                   "gof"="molecular_function")
##     gocDf[,"gocategory"]<-rep(goCategory["goc"],nrow(gocDf))
##     gofDf[,"gocategory"]<-rep(goCategory["gof"],nrow(gofDf))
##     gopDf[,"gocategory"]<-rep(goCategory["gop"],nrow(gopDf))
##     goStrictDf<-rbind(gocDf,gofDf,gopDf)
##     colnames(goStrictDf)[colnames(goStrictDf)=="goid"]<-"annotkwid"
##     colnames(goStrictDf)[colnames(goStrictDf)=="gocategory"]<-"aspect"
##     rm(gocDf)
##     rm(gofDf)
##     rm(gopDf)   
##     saveGeneratedGoslimSummaryForGoStrictAnnotationData(dbConn,proteomeId,goStrictDf)

##     ## Cálculo y cacheo de distancias ##
##     ## saveKeywordEuclideanDistancesForProteomeId(dbConn,proteomeId)
##     ## saveGeneratedGoslimEuclideanDistancesForProteomeId(dbConn,proteomeId)
##     saveEuclideanDistancesForProteomeId(dbConn,proteomeId,c("generated_goslim","keyword"))
##     geneCount
## }

###################
## Función tocha ##
######################################################################################
##  1. Lee las anotaciones
##  2. Lee el multifasta
##  3. Mete la información en la tabla de proteoma
##  4. Si es un proteoma que tiene GCA, asignarle el ncbitaxid
##  5. Si no es un proteoma de usuario, introducir su taxón y el de sus ancestros
##     en la labla de taxonomía
##  6. Mete la información básica de cada gen
##  7. A los genes les mete su secuencia sacándola del multifasta
##  8. Para cada gen mete cada campo de las anotaciones en cada tabla correspondiente
##  9. Recalcula los goslim y saca un conteo por cada combinatoria de proteína y goslim,
##     cosa que guarda en una tabla
## 10. Calcula las distancias euclídeas respecto el resto de proteomas que haya en
##     la base de datos partiendo de la información de los goslim recalculados
#########################################################-#############################
## [x] meter parámetro para que si se sube por web borre luego el fichero
## [ ] que se borre de verdad
## [x] asignar jobid auto, importar paquete uuid, que devuelva jobid
saveSma3sFileSet<-function(dbConn,
                           sma3sAnnotationFilePath,
                           multifastaFilePath,
                           tagVec=character(),
                           isUserProteome=FALSE,
                           gcaId=NA,
                           ncbiTaxId=NA,
                           dbName=NA,
                           sourceUrl=NA,
                           doUpdateNcbiAssemblyInfo=TRUE,
                           webFile=FALSE,
                           webSma3sAnnotationFileName=NA){

    sma3sAnnotationDf<-readSma3sAnnotation(sma3sAnnotationFilePath)
    sma3sAnnotationMd5Sum<- md5sum(sma3sAnnotationFilePath)
    multifastaDf<-readMultifasta(multifastaFilePath)
    multifastaMd5Sum<-md5sum(multifastaFilePath)

    geneCount<-sum(complete.cases(sma3sAnnotationDf[,"fastaheader"]), na.rm=TRUE)
    if(geneCount==0){
        error("Empty Sma3s tsv file. Proteome with no genes.")
    }

    sma3sAnnotationBaseName<-NA
    if(webFile==FALSE){
        sma3sAnnotationBaseName<-basename(tools::file_path_sans_ext(sma3sAnnotationFilePath))
    } else {
        sma3sAnnotationBaseName<-basename(tools::file_path_sans_ext(webSma3sAnnotationFileName))
    }
    
    ## ## Campos del proteoma
    proteomeInfoObj<-data.frame(
        ## proteomeid=NA, ## proteomeid asignado por la base de datos. No se mete como columna en el DF porque lo va a interpretar como que es null en dbWriteTable
        basename = sma3sAnnotationBaseName,
        sma3sannotationmd5sum = sma3sAnnotationMd5Sum,
        multifastamd5sum = multifastaMd5Sum,
        genecount = geneCount,
        is_userproteome = isUserProteome,
        gcaid = gcaId,
        ncbitaxid = ncbiTaxId,
        dbname = dbName,
        sourceurl = sourceUrl,
        do_updatencbiassemblyinfo = doUpdateNcbiAssemblyInfo,
		jobid = UUIDgenerate()
    )
        
    proteomeId<-dbGetQuery(dbConn,"
      INSERT INTO taxogeno.proteome(
        basename,
        sma3sannotationmd5sum,
        multifastamd5sum,
        genecount,
        is_userproteome,
        gcaid,
        ncbitaxid,
        dbname,
        sourceurl,
        do_updatencbiassemblyinfo,
		jobid
      ) VALUES ($1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11)
      RETURNING proteomeid 
    ",params=list(proteomeInfoObj[["basename"]],
                  proteomeInfoObj[["sma3sannotationmd5sum"]],
                  proteomeInfoObj[["multifastamd5sum"]],
                  proteomeInfoObj[["genecount"]],
                  proteomeInfoObj[["is_userproteome"]],
                  proteomeInfoObj[["gcaid"]],
                  proteomeInfoObj[["ncbitaxid"]],
                  proteomeInfoObj[["dbname"]],
                  proteomeInfoObj[["sourceurl"]],
                  proteomeInfoObj[["do_updatencbiassemblyinfo"]],
				  proteomeInfoObj[["jobid"]])) [1,"proteomeid"] ## tomar la primera entrada de la columna proteomeid, porque sólo hay una, del RETURNING

    ## ## si es un proteoma no de un usuario, sacar el posible GCA del proteoma que se ha introducido a partir de su nombre
    if(!isUserProteome && doUpdateNcbiAssemblyInfo ){
        updateNcbiAssemblyInfoForProteomeIdVec(dbConn,proteomeId)
    }

    ## actualizar a lo que hay en la bbdd
    proteomeInfoObj<-dbGetQuery(dbConn,"SELECT * FROM taxogeno.proteome WHERE proteomeid=$1", params=list(proteomeId))
    
    ## ## Etiquetas del proteoma
    tagDf<-data.frame(proteomeid = rep(proteomeId,length(tagVec)),
                      tag = tagVec )
    dbWriteTable(dbConn,SQL("taxogeno.proteome_tag_rel"),tagDf,append=TRUE,row.names=FALSE)

    ## GENES
    ## ## Campos generales de los genes
    geneDf<-sma3sAnnotationDf[,c("fastashortheader","genename","genedescription","fastaheader")]
    ## ## Engancharles sus secuencias si las hay
    geneDf<-merge(geneDf,multifastaDf,by="fastaheader",all.x=TRUE)
    ## ## Vincular con el proteoma que se está metiendo
    geneDf[,"proteomeid"]<-rep(proteomeId, nrow(geneDf))
    ## ## Escribir en la base de datos
    dbWriteTable(dbConn,SQL("taxogeno.gene"),geneDf,append=TRUE,row.names=FALSE)

    ## ## las claves autoincrementales hay que recuperarlas de la base de datos
    geneIdFastaShortHeaderRelDf<-dbGetQuery(dbConn,
                                            "SELECT geneid,fastashortheader FROM taxogeno.gene WHERE proteomeid=$1",
                                            params=list(proteomeId))
    insertedGeneAnnotationDf<-merge(sma3sAnnotationDf,geneIdFastaShortHeaderRelDf, by="fastashortheader")
    rm(sma3sAnnotationDf)
    
    ## Campos 1:n gen:campos
    enzymeDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "enzyme", newColName="ec")
    dbWriteTable(dbConn,SQL("taxogeno.gene_enzyme_rel"),enzymeDf,append=TRUE,row.names=FALSE)

    keywordDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "keyword")
    dbWriteTable(dbConn,SQL("taxogeno.gene_keyword_rel"),keywordDf,append=TRUE,row.names=FALSE)
    rm(keywordDf)

    ## #########################
    pathwayDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "pathway")
    dbWriteTable(dbConn,SQL("taxogeno.gene_pathway_rel"),pathwayDf,append=TRUE,row.names=FALSE)
    rm(pathwayDf)
    
    ## #########################
    goslimDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf, "goslim", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_goslim_rel"),goslimDf,append=TRUE,row.names=FALSE)
    rm(goslimDf)
    
    ## #########################
    gocDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"goc", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_goc_rel"),gocDf,append=TRUE,row.names=FALSE)
    gofDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gof", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_gof_rel"),gofDf,append=TRUE,row.names=FALSE)
    gopDf<-extractColumnDfFromInsertedGeneAnnotationDf(insertedGeneAnnotationDf,"gop", newColName="goid")
    dbWriteTable(dbConn,SQL("taxogeno.gene_gop_rel"),gopDf,append=TRUE,row.names=FALSE)
    goCategory<-c("goc"="cellular_component",
                  "gop"="biological_process",
                  "gof"="molecular_function")
    gocDf[,"gocategory"]<-rep(goCategory["goc"],nrow(gocDf))
    gofDf[,"gocategory"]<-rep(goCategory["gof"],nrow(gofDf))
    gopDf[,"gocategory"]<-rep(goCategory["gop"],nrow(gopDf))
    goStrictDf<-rbind(gocDf,gofDf,gopDf)
    colnames(goStrictDf)[colnames(goStrictDf)=="goid"]<-"annotkwid"
    colnames(goStrictDf)[colnames(goStrictDf)=="gocategory"]<-"aspect"
    rm(gocDf)
    rm(gofDf)
    rm(gopDf)   
    saveGeneratedGoslimSummaryForGoStrictAnnotationData(dbConn,proteomeId,goStrictDf)

    ## Cálculo y cacheo de distancias ##
    ## saveKeywordEuclideanDistancesForProteomeId(dbConn,proteomeId)
    ## saveGeneratedGoslimEuclideanDistancesForProteomeId(dbConn,proteomeId)
    saveEuclideanDistancesForProteomeId(dbConn,proteomeId,c("generated_goslim","keyword"))
    
    proteomeInfoObj[["jobid"]]
}

                                        # ################################ #
                                        # Funciones para alterar proteomas #
                                        # ################################ #
updateProteomeFieldInDB<-function(dbConn,proteomeId,fieldName, newValue){
    if(length(proteomeId)>1){
        stop("length(proteomeId)>1")
    }
    allowedFieldNames<-c("is_userproteome", "gcaid", "ncbitaxid", "scientific_name")
    if(!(all(fieldName %in% allowedFieldNames))){
        error("!(fieldName %in% allowedFieldNames)")
    }
    dbExecute(dbConn, sprintf("UPDATE taxogeno.proteome SET %s = $1 WHERE proteomeid = $2", fieldName), params(list(newValue,proteomeId)))
    
}
getLastProteomeId<-function(dbConn){
    dbGetQuery(dbConn,"SELECT MAX(proteomeid) as max_proteomeid FROM taxogeno.proteome")[1,"max_proteomeid"]
}
deleteProteome<-function(dbConn,proteomeIdVec){
    for(proteomeId in proteomeIdVec){
        dbExecute(dbConn,"CALL taxogeno.delete_proteome($1)",params=list(proteomeId))
    }
}

                                        # ################################################# #
                                        # Funciones para actualizar bases de datos externas #
                                        # ################################################# #
updateNcbiAssemblySummaryGenbank<-function(dbConn){
    ## semiconstants
    assemblySummaryGenbankUrl<-"ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt"
    assemblySummaryGenbankFilePath<-"assembly_summary_genbank.txt"
    assemblySummaryGenbankColNames<-c("assembly_accession", "bioproject", "biosample", "wgs_master",
                                      "refseq_category", "taxid", "species_taxid", "organism_name",
                                      "infraspecific_name", "isolate", "version_status", "assembly_level",
                                      "release_type", "genome_rep", "seq_rel_date", "asm_name", "submitter",
                                      "gbrs_paired_asm", "paired_asm_comp", "ftp_path", "excluded_from_refseq",
                                      "relation_to_type_material")
    
    write("Downloading assembly summary",stdout())
    system2("wget", args=c(assemblySummaryGenbankUrl,"-O",assemblySummaryGenbankFilePath), wait=TRUE)
    write("Loading assembly summary table",stdout())
    assemblySummaryGenbankDf<-read.table(assemblySummaryGenbankFilePath,
                                         comment.char="#",
                                         header=FALSE,
                                         row.names=NULL,
                                         col.names = assemblySummaryGenbankColNames,
                                         colClasses = rep("character",22),
                                         fill=TRUE,
                                         quote=NULL,
                                         sep="\t",
                                         stringsAsFactors=FALSE)
    unlink(assemblySummaryGenbankFilePath)
    write("Inserting assembly summary table into database",stdout())
    dbExecute(dbConn,"DELETE FROM ncbi.assembly_summary_genbank")
    dbWriteTable(dbConn,SQL("ncbi.assembly_summary_genbank"), assemblySummaryGenbankDf, append=TRUE, row.names=FALSE)
}

updateBiosqlNcbiTaxonomy<-function(dbConn){
    ncbiTaxonomyDumpUrl="ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz"
    ncbiTaxonomyDumpFile="taxdump.tar.gz"

    write("Downloading taxonomy dump",stdout())
    system2("wget", args=c(ncbiTaxonomyDumpUrl,"-O",ncbiTaxonomyDumpFile), wait=TRUE)

    write("Uncompressing taxonomy dump",stdout())
    taxonomyDumpFileList<-untar(ncbiTaxonomyDumpFile,list=TRUE)
    untar(ncbiTaxonomyDumpFile, exdir=dirname(ncbiTaxonomyDumpFile))
    unlink(ncbiTaxonomyDumpFile)
    
    write("Loading taxonomy dump into biosql schema of database",stdout())
    system2("perl",
            args=c("./bin/load_ncbi_taxonomy.pl","--driver","Pg","--schema","biosql","--dbname" ,"taxogeno","--directory",dirname(ncbiTaxonomyDumpFile)),
            wait=TRUE)
}

updateUniprotKeyword<-function(dbConn){
    ## semiconstants
    uniprotKeywordOboUrl<-"'http://www.uniprot.org/keywords/?query=*&format=obo'"
    uniprotKeywordOboFilePath<-"keyword.obo"
    system2("wget", args=c(uniprotKeywordOboUrl,"-O",uniprotKeywordOboFilePath), wait=TRUE)

    uniprotKeywordTsvUrl<-"'http://www.uniprot.org/keywords/?query=*&format=tab'"
    uniprotKeywordColNames<-c("keywordid","keyword","keyworddescription","keywordcategory")
    uniprotKeywordTsvFilePath<-"uniprot_keyword.tsv"
    system2("wget", args=c(uniprotKeywordTsvUrl,"-O",uniprotKeywordTsvFilePath), wait=TRUE)
    uniprotKeywordDf<-read.table(uniprotKeywordTsvFilePath,
                                 header=TRUE,
                                 row.names=NULL,
                                 col.names=uniprotKeywordColNames,
                                 quote="",
                                 sep="\t",
                                 stringsAsFactors=FALSE)
    unlink(uniprotKeywordTsvFilePath)
    dbExecute(dbConn,"DELETE FROM uniprot.uniprot_keyword")
    dbExecute(dbConn,"DELETE FROM uniprot.uniprot_keyword_category")
    dbWriteTable(dbConn,SQL("uniprot.uniprot_keyword_category"), unique(uniprotKeywordDf[,"keywordcategory",drop=FALSE]), append=TRUE,row.names=FALSE)
    dbWriteTable(dbConn,SQL("uniprot.uniprot_keyword"), uniprotKeywordDf,append=TRUE,row.names=FALSE)
}
    
## ##https://www.uniprot.org/docs/pathlist.txt
## uniprotProteomeTsvUrl<-"http://www.uniprot.org/proteomes/?query=*&format=tab&force=true&columns=id,busco,cpd,assembly"
## uniprotProteomeTsvFilePath<-"uniprot_proteome.tsv"
## system2("wget", args=c(uniprotProteomeTsvUrl,"-O",uniprotProteomeTsvFilePath), wait=TRUE)
## uniprotProteomeColNames<-c("uniprotproteomeid","busco","cpd","gcaid")
## uniprotProteomeDf<-read.table("uniprot_proteome.tsv",
##                              header=TRUE,
##                              row.names=NULL,
##                              col.names=uniprotProteomeColNames,
##                              quote="",
##                              sep="\t",
##                              stringsAsFactors=FALSE)
## dbWriteTable(dbConn,SQL("uniprot.uniprot_proteome"),uniprotProteomeDf)
## ##

## updateGeneOntology()
## [ ] dónde meter los go ¿en el directorio de trabajo?
## [x] delete cascade gene_ontology.gene_ontology
updateGeneOntology<-function(dbConn){
    insertGoOwl<-function(dbConn,owlFileName){
        dbExecute(dbConn,"DELETE from gene_ontology.gene_ontology_xml WHERE filename=$1",params=list(owlFileName))
        owlContents <- rawToChar(readBin(owlFileName, "raw", file.info(owlFileName)$size))
        dbExecute(dbConn,"")
        dbExecute(dbConn,
                  "insert into gene_ontology.gene_ontology_xml (filename, xmldata) values($1,XMLPARSE(DOCUMENT $2))",
                  params=list(owlFileName,owlContents))
        rm(owlContents)
    }
    deleteTablifiedGoOwl<-function(dbConn,tableName){
        dbExecute(dbConn,sprintf("DELETE from gene_ontology.%s",tableName))
    }
    insertTablifiedGoOwl<-function(dbConn,tableName,owlFileName){
        dbExecute(dbConn,
                  sprintf("insert into gene_ontology.%s(goid,golabel,goaspect)
                       select
                         gene_ontology_tablified.goid     as goid,
                         gene_ontology_tablified.golabel  as golabel,
                         gene_ontology_tablified.goaspect as goaspect
                       --
                        from
                          -- table with xml source column inside
                          gene_ontology.gene_ontology_xml,
                          -- generated table from xml source
                          xmltable(
                            -- namespaces
                            xmlnamespaces(
                              'http://www.w3.org/1999/02/22-rdf-syntax-ns#'   as \"rdf\",
                              'http://www.w3.org/2002/07/owl#'                as \"owl\",
                              'http://www.geneontology.org/formats/oboInOwl#' as \"oboInOwl\",
                              'http://www.w3.org/2000/01/rdf-schema#'         as \"rdfs\"
                            ),
                            -- path
                            '/rdf:RDF/owl:Class'
                            -- xml source
                            passing
                              gene_ontology.gene_ontology_xml.xmldata
                            -- columns
                            columns
                              goid     text PATH 'oboInOwl:id/text()',
                              golabel  text PATH 'rdfs:label/text()',
                              goaspect text PATH 'oboInOwl:hasOBONamespace/text()'
                          ) gene_ontology_tablified -- generated table alias
                       --
                        where
                         gene_ontology.gene_ontology_xml.filename=$1
                        and
                         gene_ontology_tablified.goid is not null",
                       tableName),
                  params=list(owlFileName))
    }

    ## Descargar los owl
    
    goOwlUrl<-"http://current.geneontology.org/ontology/go.owl"
    goOwlFileName<-"go.owl"
    system2("wget", args=c(goOwlUrl,"-O",goOwlFileName), wait=TRUE)

    goSlimGenericOwlUrl<-"http://current.geneontology.org/ontology/subsets/goslim_generic.owl"
    goSlimGenericOwlFileName<-"goslim_generic.owl"
    system2("wget", args=c(goSlimGenericOwlUrl,"-O",goSlimGenericOwlFileName), wait=TRUE)

    ## ## ## ## ##

    insertGoOwl(dbConn,"go.owl")
    insertGoOwl(dbConn,"goslim_generic.owl")
    
    deleteTablifiedGoOwl(dbConn,"goslim_generic")
    deleteTablifiedGoOwl(dbConn,"gene_ontology")

    insertTablifiedGoOwl(dbConn,"gene_ontology","go.owl")
    insertTablifiedGoOwl(dbConn,"goslim_generic","goslim_generic.owl")
}
