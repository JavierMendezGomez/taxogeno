#!/usr/bin/Rscript
source("taxogeno_library.R")
library(shiny)
library(shinyTree)
library(ggplot2)
## https://christophergandrud.github.io/networkD3/
library(networkD3)
library(ggiraph)
library(stringr)
library(parallel)
library("DBI")
library("jsonlite")

setwd("/opt/taxogeno-clean")

##https://atrebas.github.io/post/2019-06-08-lightweight-dendrograms/

proteomeGoAnnotationSummary <- function(conn, proteomeid){
    proteomeGoAnnotationSummaryResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.generated_goslim_summary.proteomeid,
        taxogeno.generated_goslim_summary.goid,
        taxogeno.generated_goslim_summary.annotkwcount,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.generated_goslim_summary
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.generated_goslim_summary.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
    ")
    dbBind(proteomeGoAnnotationSummaryResult, list(proteomeid))
    proteomeGoAnnotationSummaryDF <- dbFetch(proteomeGoAnnotationSummaryResult)
    dbClearResult(proteomeGoAnnotationSummaryResult)
    rm(proteomeGoAnnotationSummaryResult)
    
    proteomeGoAnnotationSummaryDF
}

proteomeList<-function(conn){
    proteomeListResult <- dbSendQuery(conn, "
      SELECT
        proteomeid,
        basename,
        gcaid,
        prot.ncbitaxid as ncbitaxid,
        scientific_name
      FROM taxogeno.proteome prot
      LEFT JOIN taxogeno.taxonomy tax
      ON prot.ncbitaxid=tax.ncbitaxid
      WHERE is_userproteome = FALSE
    ")
    proteomeListDF<- dbFetch(proteomeListResult)
    dbClearResult(proteomeListResult)
    rm(proteomeListResult)

    proteomeListDF
}
userProteomeList<-function(conn,jobId){
    proteomeListDf <- dbGetQuery(conn, "
      SELECT
        proteomeid,
        basename,
        gcaid,
        prot.ncbitaxid as ncbitaxid,
        scientific_name
      FROM taxogeno.proteome prot
      LEFT JOIN taxogeno.taxonomy tax
      ON prot.ncbitaxid=tax.ncbitaxid
      WHERE is_userproteome = TRUE
      AND jobid=$1
    ",params=list(jobId))
    proteomeListDf
}

proteomeInfo <- function (conn, proteomeid){
    proteomeInfoResult <- dbSendQuery(conn, "
      SELECT
        prot.proteomeid,
        prot.basename,
        prot.creationtimestamp,
        prot.is_userproteome,
        prot.genecount,
        prot.gcaid,
        prot.ncbitaxid,
        prot.dbname,
        prot.sourceurl,

        tax.scientific_name,
        tax.node_rank

      FROM taxogeno.proteome prot
        LEFT JOIN taxogeno.taxonomy tax ON (prot.ncbitaxid=tax.ncbitaxid)
      WHERE proteomeid= $1
    ")
    dbBind(proteomeInfoResult, list(proteomeid))
    proteomeInfoDF<- dbFetch(proteomeInfoResult)
    dbClearResult(proteomeInfoResult)
    rm(proteomeInfoResult)

    proteomeInfoDF[1,]
}

proteomeIdFromJobId <- function (conn, jobid){
    proteomeIdResult <- dbSendQuery(conn, "
      SELECT
        prot.proteomeid as proteomeid,
      FROM taxogeno.proteome prot
      WHERE jobid= $1
    ")
    dbBind(proteomeIdResult, list(jobid))
    proteomeIdDF<- dbFetch(proteomeIdResult)
    dbClearResult(proteomeIdResult)
    rm(proteomeIdResult)

    proteomeInfoDF[1,"proteomeid"]
}

geneList <- function(conn, proteomeid){
    geneListResult <- dbSendQuery(conn,"
      SELECT
        geneid,
        genename,
        genedescription,
        fastaheader
      FROM taxogeno.gene WHERE proteomeid= $1
    ")
    dbBind(geneListResult, list(proteomeid))
    geneListDF <- dbFetch(geneListResult)
    dbClearResult(geneListResult)
    rm(geneListResult)
    
    geneListDF
}
geneListColumns<-list("geneid","genename","genedescription","fastaheader","aasequence")

geneInfo <- function(conn, geneid){
    geneInfoResult <- dbSendQuery(conn,"
      SELECT
        geneid,
        genename,
        genedescription,
        fastaheader,
        aasequence
      FROM taxogeno.gene WHERE geneid= $1
    ")
    dbBind(geneInfoResult, list(geneid))
    geneInfoDF <- dbFetch(geneInfoResult)
    dbClearResult(geneInfoResult)
    rm(geneInfoResult)
    
    geneInfoDF[1:1,]
}
geneInfoColumns<-list("geneid","genename","genedescription","fastaheader","aasequence")

#########################################
geneGoAnnotationList <- function(conn, geneid){
    geneGoAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_goc_rel.geneid,
        taxogeno.gene_goc_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_goc_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_goc_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
UNION
      SELECT
        taxogeno.gene_gof_rel.geneid,
        taxogeno.gene_gof_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_gof_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_gof_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1
UNION
      SELECT
        taxogeno.gene_gop_rel.geneid,
        taxogeno.gene_gop_rel.goid,
        gene_ontology.gene_ontology.golabel,
        gene_ontology.gene_ontology.goaspect

      FROM taxogeno.gene_gop_rel
        INNER JOIN gene_ontology.gene_ontology on (taxogeno.gene_gop_rel.goid = gene_ontology.gene_ontology.goid)
      WHERE geneid= $1

    ")
    dbBind(geneGoAnnotationListResult, list(geneid))
    geneGoAnnotationListDF <- dbFetch(geneGoAnnotationListResult)
    dbClearResult(geneGoAnnotationListResult)
    rm(geneGoAnnotationListResult)
    
    geneGoAnnotationListDF
}

#############################################
geneKeywordAnnotationList <- function(conn, geneid){
    geneKeywordAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_keyword_rel.geneid,
        taxogeno.gene_keyword_rel.keyword,
        uniprot.uniprot_keyword.keywordid,
        uniprot.uniprot_keyword.keywordcategory,
        uniprot.uniprot_keyword.keyworddescription

      FROM taxogeno.gene_keyword_rel
        INNER JOIN uniprot.uniprot_keyword on (taxogeno.gene_keyword_rel.keyword = uniprot.uniprot_keyword.keyword)
      WHERE geneid= $1
    ")
    dbBind(geneKeywordAnnotationListResult, list(geneid))
    geneKeywordAnnotationListDF <- dbFetch(geneKeywordAnnotationListResult)
    dbClearResult(geneKeywordAnnotationListResult)
    rm(geneKeywordAnnotationListResult)
    
    geneKeywordAnnotationListDF
}

#############################################
geneEnzymeAnnotationList <- function(conn, geneid){
    geneEnzymeAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_enzyme_rel.geneid,
        taxogeno.gene_enzyme_rel.ec
      FROM taxogeno.gene_enzyme_rel
      WHERE geneid= $1
    ")
    dbBind(geneEnzymeAnnotationListResult, list(geneid))
    geneEnzymeAnnotationListDF <- dbFetch(geneEnzymeAnnotationListResult)
    dbClearResult(geneEnzymeAnnotationListResult)
    rm(geneEnzymeAnnotationListResult)
    
    geneEnzymeAnnotationListDF
}

#############################################
genePathwayAnnotationList <- function(conn, geneid){
    genePathwayAnnotationListResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.gene_pathway_rel.geneid,
        taxogeno.gene_pathway_rel.pathway
      FROM taxogeno.gene_pathway_rel
      WHERE geneid= $1
    ")
    dbBind(genePathwayAnnotationListResult, list(geneid))
    genePathwayAnnotationListDF <- dbFetch(genePathwayAnnotationListResult)
    dbClearResult(genePathwayAnnotationListResult)
    rm(genePathwayAnnotationListResult)
    
    genePathwayAnnotationListDF
}

getGoslimAnnotationForProteomeIdVector<- function(conn, proteomeIdVector, aspectVector){

    proteomeIdArrayString<-paste('{',
                                 paste(proteomeIdVector,collapse=','),
                                 '}',
                                 sep='')
    aspectArrayString<-paste('{',
                             paste(aspectVector,collapse=','),
                             '}',
                             sep='')
    goslimAnnotationResult <- dbSendQuery(conn,"
      SELECT
        taxogeno.generated_goslim_summary.proteomeid as proteomeid,
        taxogeno.generated_goslim_summary.annotkwid  as goid,
        gene_ontology.goslim_generic.goaspect  as goaspect,
        gene_ontology.goslim_generic.golabel   as golabel,

        taxogeno.generated_goslim_summary.annotkwcount as gocount,
        taxogeno.proteome.genecount                    as genecount,
        taxogeno.generated_goslim_summary.annotkwcount::float/taxogeno.proteome.genecount::float as relval
      FROM taxogeno.generated_goslim_summary
      INNER JOIN taxogeno.proteome
        ON (taxogeno.generated_goslim_summary.proteomeid=taxogeno.proteome.proteomeid)
      INNER JOIN gene_ontology.goslim_generic
        ON taxogeno.generated_goslim_summary.annotkwid = gene_ontology.goslim_generic.goid
      WHERE taxogeno.generated_goslim_summary.proteomeid = ANY($1)
      AND taxogeno.generated_goslim_summary.annotkwcount > 0
      AND gene_ontology.goslim_generic.golabel NOT IN ('molecular_function','cellular_component','biological_process')
      AND gene_ontology.goslim_generic.goaspect = ANY($2)
    ")
    dbBind(goslimAnnotationResult, list(proteomeIdArrayString, aspectArrayString))
    goslimAnnotationDF <- dbFetch(goslimAnnotationResult)
    dbClearResult(goslimAnnotationResult)
    rm(goslimAnnotationResult)
    
    goslimAnnotationDF
}

getTaxonomyNestedStructure<-function(conn, ncbitaxid){
    taxonResult <- dbSendQuery(conn,"
      SELECT taxogeno.taxonomy_jsonb_children($1) as json_result
    ")
    dbBind(taxonResult, list(ncbitaxid))
    taxonDF <- dbFetch(taxonResult)
    dbClearResult(taxonResult)
    rm(taxonResult)
    taxonList<-fromJSON(taxonDF$json_result[1], simplifyVector=FALSE)
    rm(taxonDF)
    
    childrenTreeTraversal<-function(nodeListObj, initial=FALSE){
        
        childrenObjList<-lapply(
            nodeListObj[["children"]],
            childrenTreeTraversal
        )

        namesChildrenObjList<- function (childrenObjList) {
            lapply(
                childrenObjList,
                function(childObj){
                    objListLabel<-paste(
                        "ncbitaxid:",attributes(childObj)[["ncbitaxid"]],"; ",
                        attributes(childObj)[["node_rank"]]," ",
                        attributes(childObj)[["scientific_name"]],
                        sep=""
                    ) ## end paste
                    
                    ## return:
                    objListLabel
                } ## end function
            ) ## end lapply
        } ## end function

        names(childrenObjList)<-namesChildrenObjList(childrenObjList)
        
        objList<-structure(childrenObjList,
                           ncbitaxid=nodeListObj[["ncbitaxid"]],
                           parent_ncbitaxid=nodeListObj[["parent_ncbitaxid"]],
                           scientific_name=nodeListObj[["scientific_name"]],
                           node_rank=nodeListObj[["node_rank"]])

        if (initial==TRUE){
            wrappedObjList<-list(objList)
            names(wrappedObjList) <- namesChildrenObjList(wrappedObjList)
            wrappedObjList
        } else {
            objList
        }
    }
    childrenTreeTraversal(taxonList,TRUE)
}

getTaxonomySelectedNodes <- function(tree, ancestry=NULL, vec=list()){
    if (is.list(tree)){
        for (i in 1:length(tree)){
            anc <- c(ancestry, names(tree)[i])
            vec <- getTaxonomySelectedNodes(tree[[i]], anc, vec)
        }    
    }
    a <- attr(tree, "stselected", TRUE)
    if (!is.null(a) && a == TRUE){
                                        # Get the element name
        len_anc <- length(ancestry)
        el <- ancestry[len_anc]
        vec[length(vec)+1] <- el

        ## Save some attributes
        lapply(names(attributes(tree)),function(attribute){
            if(grepl("^st",attribute)
               || attribute == "ncbitaxid"
               || attribute == "parent_ncbitaxid"
               || attribute == "scientific_name"
               || attribute == "node_rank"){
                attr(vec[[length(vec)]], attribute) <<- attr(tree,attribute)
            }
        })
    }
    return(vec)
}

getReferenceProteomeIdsForTaxonomySelectedNodes<-function(conn, tree){
    ncbiTaxIdVector<-as.integer(
        sapply(
            getTaxonomySelectedNodes(tree),
            function(node){attr(node,"ncbitaxid")}
        )
    )
    refProteomeListResult <- dbSendQuery(conn,"
      SELECT
        proteomeid
      FROM taxogeno.proteome
      WHERE ncbitaxid in ($1)
      AND is_userproteome = FALSE
    ")
    dbBind(refProteomeListResult, list(ncbiTaxIdVector))
    refProteomeListDF <- dbFetch(refProteomeListResult)
    dbClearResult(refProteomeListResult)
    rm(refProteomeListResult)
    
    refProteomeListDF$proteomeid
}

getEuclideanInterproteomeDistances <- function(conn,
                                               annotationType,
                                               proteomeId,
                                               proteomeIdVector){
    ## distanceListResult <- dbSendQuery(conn,sprintf("
    ##    SELECT
    ##      proteomeid_least,
    ##      (SELECT tax.scientific_name
    ##       FROM taxogeno.proteome AS prot
    ##       INNER JOIN taxogeno.taxonomy AS tax
    ##       ON (prot.ncbitaxid=tax.ncbitaxid)
    ##       WHERE prot.proteomeid=eucdisttab.proteomeid_least) as ref_scientific_name,

    ##      proteomeid_greatest,
    ##      (SELECT tax.scientific_name
    ##       FROM taxogeno.proteome AS prot
    ##       INNER JOIN taxogeno.taxonomy AS tax
    ##       ON (prot.ncbitaxid=tax.ncbitaxid)
    ##       WHERE prot.proteomeid=eucdisttab.proteomeid_greatest) as comp_scientific_name,
    ##      distance
    ##    FROM taxogeno.euclidean_distances_%s AS eucdisttab
    ##    WHERE 
    ##            (proteomeid_least = %2$s and proteomeid_greatest in (%3$s))
    ##            or
    ##            (proteomeid_least in (%3$s) and proteomeid_greatest = %2$s)
    ## ",annotationType,proteomeId, paste(proteomeIdVector,collapse=",")))


    distanceListResult <- dbSendQuery(conn,sprintf("
       SELECT
         proteomeid_least,
         (SELECT basename
          FROM taxogeno.proteome AS prot
          WHERE prot.proteomeid=eucdisttab.proteomeid_least) as ref_scientific_name,

         proteomeid_greatest,
         (SELECT basename
          FROM taxogeno.proteome AS prot
          WHERE prot.proteomeid=eucdisttab.proteomeid_greatest) as comp_scientific_name,
         distance
       FROM taxogeno.euclidean_distances_%s AS eucdisttab
       WHERE 
               (proteomeid_least = %2$s and proteomeid_greatest in (%3$s))
               or
               (proteomeid_least in (%3$s) and proteomeid_greatest = %2$s)
    ",annotationType,proteomeId, paste(proteomeIdVector,collapse=",")))
    
    
    distanceListDF <- dbFetch(distanceListResult)
    dbClearResult(distanceListResult)
    rm(distanceListResult)
    
    distanceListDF
##########
}

## Página con tres pestañas
## Cada pestaña contiene un sidebarLayout
ui <- navbarPage(
    title = "Taxogeno",
    tabPanel(
        title= "Step 0: Upload proteome",
        sidebarLayout(
            sidebarPanel(
                helpText("Upload here your proteome annotated with Sma3s"),
                width=2
            ),
            mainPanel(
                column(
                    width=4,
                    h3("Upload form"),
                    textInput("tags", "Insert comma separated tags for this proteome"),
                    fileInput("tsv", "Choose Sma3s output TSV file", accept = ".tsv"),
                    fileInput("multifasta", "Choose multifasta file", accept = ".fasta"),
                    actionButton("sendButton","Insert proteome")
                ),
                column(
                    width=4,
                    h3("Warnings and errors")
                    ## Errores y advertencias podrían ser que los ficheros tuvieran nombres distintos, que los ficheros no se ajustaran al formato (por ejemplo que se subiera el summary en lugar del tsv con todos los datos) o que el multifasta no fuera un multifasta, o el resumen de si cuadra todo o no (esto no está implementado, esto depende de modificar taxogeno.R)
                ),
                column(
                    width=4,
                    h3("Upload result"),
                    h4("Job ID for uploaded proteome"),
                    verbatimTextOutput("jobId",placeholder=TRUE)
                )
            )
        )
    ),
    tabPanel(
        title = "Step 1: Select a proteome to compare",
        sidebarLayout(    
            sidebarPanel(
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Proteome selection"',
                    helpText("Proteome selection.")
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "General info"',
                    helpText("General information about proteome.")
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Annotation info"',
                    helpText("Information about functional annotation of proteome.")
                ),
                conditionalPanel(
                    condition =
                        'input.proteomeSubsection === "Genes info"',
                    helpText("Genes info conditional panel."),
                    conditionalPanel(
                        condition = 'input.genesInfoTabsetPanel === "Genes table"'
                    )
                ),
                conditionalPanel(
                    condition = 'input.proteomeSubsection === "Similarity info"',
                    helpText("Display 5 records by default.")
                )
               ,width=2),
            
            mainPanel(
                tabsetPanel(
                    id = 'proteomeSubsection',
                    tabPanel(
                        title = "Proteome selection",
                        h3("Proteome list"),
                        p("If you have uploaded a proteome, insert your jobid below. Your proteome is not visible for anyone. A valid jobid is a condition for it to appear."),
                        textInput("jobid", "Insert your jobid"),
                        DT::dataTableOutput("proteomeList")
                    ),
                    tabPanel(
                        title = "General info",
                        h3("General info"),
                        column(
                            6,
                            h4("Internal proteome Identificator in this database"),
                            verbatimTextOutput("proteomeInfo.proteomeid")
                        ),
                        column(
                            6,
                            h4("Common name of uploaded files"),
                            verbatimTextOutput("proteomeInfo.basename")
                        ),
                        column(
                            6,
                            h4("Date and time of insertion in database"),
                            verbatimTextOutput("proteomeInfo.creationtimestamp")
                        ),
                        column(
                            6,
                            h4("Has been this proteome uploaded by an user?"),
                            verbatimTextOutput("proteomeInfo.is_userproteome")
                        ),
                        column(
                            6,
                            h4("Number of genes"),
                            verbatimTextOutput("proteomeInfo.genecount")
                        ),
                        column(
                            6,
                            h4("Assembly"),
                            verbatimTextOutput("proteomeInfo.gcaid")
                        ),                        
                        column(
                            6,
                            h4("NCBI Taxonomy Identificator"),
                            verbatimTextOutput("proteomeInfo.ncbitaxid")
                        ),
                        column(
                            6,
                            h4("Scientific name"),
                            verbatimTextOutput("proteomeInfo.scientific_name")
                        ),
                        column(
                            6,
                            h4("Database of origin"),
                            verbatimTextOutput("proteomeInfo.dbname")
                        ),
                        column(
                            6,
                            h4("Proteome URL resource"),
                            verbatimTextOutput("proteomeInfo.sourceurl")
                        ) 
                    ), #end tabPanel
                    ## ###################################
                    ## GENES INFO
                    ## ###################################
                    tabPanel(
                        title = "Genes info",
                        ##Los elementos de la pestaña sep. por comas
                        h3("Genes info"),
                        column(
                            12,
                            tabsetPanel(
                                id='genesInfoTabsetPanel',
                                tabPanel(
                                    title = "Genes table",
                                    ##Los elementos de la pestaña sep. por comas
                                    h4("Gene table"),
                                    DT::dataTableOutput("geneList")
                                ), #end tabPanel
                                tabPanel(
                                    title="Gene info",
                                    h4("Gene info"),
                                    h5("geneid"),
                                    verbatimTextOutput("geneInfo.geneid"),
                                    h5("genename"),
                                    verbatimTextOutput("geneInfo.genename"),
                                    h5("genedescription"),
                                    verbatimTextOutput("geneInfo.genedescription"),
                                    h5("fastaheader"),
                                    verbatimTextOutput("geneInfo.fastaheader"),
                                    h5("aasequence"),
                                    verbatimTextOutput("geneInfo.aasequence")
                                ), #end tabPanel
                                tabPanel(
                                    title="Gene annotation info",
                                    h4("Gene GO annotations"),
                                    DT::dataTableOutput("geneGoAnnotationList"),
                                    h4("Gene Keyword annotations"),
                                    DT::dataTableOutput("geneKeywordAnnotationList"),
                                    h4("Gene EC annotations"),
                                    DT::dataTableOutput("geneEnzymeAnnotationList"),
                                    h4("Gene Pathway annotations"),
                                    DT::dataTableOutput("genePathwayAnnotationList")
                                )
                            ) #end tabsetPanel
                        ) #end column
                    ) #end tabPanel
                ) #end tabsetPanel
               ,width=10) #end mainPanel
        ) #end sidebarLayout
    ), #end tabPanel
    tabPanel(
        title = "Step 2: Select a set of reference proteomes",
        sidebarLayout(
            sidebarPanel(
                helpText("Step 1: Select proteomes to compare from taxonomic tree")
            ),
            mainPanel(
                shinyTree("tree", checkbox = TRUE, search = TRUE)
            )
        )
    ),
    tabPanel(
        title = "Step 3: Check annotation info",
        sidebarLayout(
            sidebarPanel(
                helpText("Here there is only information about generated GO Slim")
            ),
            mainPanel(
                h3("Annotation info"),
                column(
                    12,
                    tabsetPanel(
                        id='annotationInfoTabsetPanel',
                        tabPanel(
                            title="GOslim nnotation",
                            h3("GOslim annotation"),
                            column(
                                12,
                                tabsetPanel(
                                    tabPanel(
                                        title="GOslim Annotation, boxplot",
                                        ##Los elementos de la pestaña sep. por comas
                                        h4("Molecular function GOslim Annotation, boxplot"),
                                        girafeOutput("molecularFunctionGoslimAnnotationBoxplot"),
                                        h4("Biological process GOslim Annotation, boxplot"),
                                        girafeOutput("biologicalProcessGoslimAnnotationBoxplot"),
                                        h4("Cellular component GOslim Annotation, boxplot"),
                                        girafeOutput("cellularComponentGoslimAnnotationBoxplot")
                                    ),
                                    tabPanel(
                                        title="GOslim Annotation, table",
                                        h3("Molecular Function GOslim Annotation, table"),
                                        DT::dataTableOutput("molecularFunctionGoslimAnnotationList"),
                                        h3("Biological process GOslim Annotation, table"),
                                        DT::dataTableOutput("biologicalProcessGoslimAnnotationList"),
                                        h3("Cellular component GOslim Annotation, table"),
                                        DT::dataTableOutput("cellularComponentGoslimAnnotationList")
                                    )
                                ) ##end tabset panel goslim boxplot/list
                            ) ## end column
                        ) ## end tabpanel goslim
                    ) ##end tabset panel goslim/kw/...
                ) ## end column
            )
        )
    ), #end tabPanel
    
    tabPanel(
        title = "Step 4: Check similarity info",
        sidebarLayout(
            sidebarPanel(
                helpText("For this funcionality to work, you must select a proteome in step 1 and a set of proteomes from taxonomic tree in step 2"),
                helpText("You can review generated data in tabulated format and download it. Also there is a graphical representation of similarities between selected proteomes."),
                ## Para hacer la similitud debe seleccionarse un proteoma en el paso 1 y un conjunto de proteomas de referencia. Sólo cuando ambos estén seleccionados deben salir los elementos gráficos
                width=2
            ),
            mainPanel(
                h3("Similarity info between selected proteome and reference proteomes"),
                tabsetPanel(
                    tabPanel(
                        title="Similarity Info, tabulated",
                        DT::dataTableOutput("similarityList"),
                        downloadButton("downloadData", label = "Download table")
                    ),
                    tabPanel(
                        title="Similarity Info, dendrogram",
                        plotOutput("euclideanDistancesDendogramPlot",height="2000px")
                    )
                )
            )
        )
    ) #end tabPanel
) #end fluidPage

server <- function(input, output) {
    options(shiny.maxRequestSize=128*1024^2) 
    conn <- dbConnect(RPostgres::Postgres(),
                      host = '/tmp/',
                      dbname = 'taxogeno')

    ##Proteomes list
    proteomeListDfReactive<-reactive({
        jobIdInput<-input$jobid
        if(is.na(jobIdInput) || is.null(jobIdInput) || jobIdInput==""){
            jobIdInput<-NA
        }

        if(is.na(jobIdInput))
            proteomeList(conn)
        else
            userProteomeList(conn,jobIdInput)
    })
    
    
    output$proteomeList<-DT::renderDataTable({
        DT::datatable (
                # proteomeListDF[, input$proteomeListColumnsShown, drop = FALSE],
                proteomeListDfReactive(),
                selection="single",
                rownames=proteomeListDfReactive()$proteomeid
            )
    })

    ## reactive para tener el proteoma seleccionado de la lista que viene de la base de datos
    dbProteomeIdSelectedReactive<-reactive({
        proteomeIdSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        proteomeIdSelected
    })

    ##proteomeId selected reactive
    ## puede venir o bien de un proteoma de usuario o bien de un proteoma de la base de datos
    userProteomeIdSelectedReactive<-reactive({
        proteomeIdSelected<-proteomeIdFromJobId(input$jobId)
        proteomeIdSelected
    })
    
    ##Proteome detail
    proteomeInfoReactive<-eventReactive(input$proteomeSubsection,{
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        proteomeInfo(conn,proteomeidSelected)
    })

    output[["proteomeInfo.proteomeid"]]<-renderText(proteomeInfoReactive()[["proteomeid"]])
    output[["proteomeInfo.basename"]]<-renderText(proteomeInfoReactive()[["basename"]])
    output[["proteomeInfo.creationtimestamp"]]<-renderText(proteomeInfoReactive()[["creationtimestamp"]])
    output[["proteomeInfo.is_userproteome"]]<-renderText(proteomeInfoReactive()[["is_userproteome"]])
    output[["proteomeInfo.genecount"]]<-renderText(proteomeInfoReactive()[["genecount"]])
    output[["proteomeInfo.ncbitaxid"]]<-renderText(proteomeInfoReactive()[["ncbitaxid"]])
    output[["proteomeInfo.scientific_name"]]<-renderText(proteomeInfoReactive()[["scientific_name"]])
    output[["proteomeInfo.gcaid"]]<-renderText(proteomeInfoReactive()[["gcaid"]])
    output[["proteomeInfo.dbname"]]<-renderText(proteomeInfoReactive()[["dbname"]])
    output[["proteomeInfo.sourceurl"]]<-renderText(proteomeInfoReactive()[["sourceurl"]])

    ## ############
    ## Gene
    ## ############
    geneListReactive<-eventReactive(input$genesInfoTabsetPanel,{
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        geneList(conn,proteomeidSelected)
    })
    output$geneList<-DT::renderDataTable({
        DT::datatable(
                geneListReactive(),
                selection="single",
                rownames=geneListReactive()$geneid
            ) #end datatable
    })

    ## ############
    ## Annotations
    ## ############
    geneGoAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneGoAnnotationList(conn,geneidSelected) ## !!
    })
    output$geneGoAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneGoAnnotationListReactive(),
                selection="single",
                rownames=geneGoAnnotationListReactive()$geneid
            )
    })
    ##
    geneKeywordAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneKeywordAnnotationList(conn,geneidSelected) ## !!
    })
    output$geneKeywordAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneKeywordAnnotationListReactive(),
                selection="single",
                rownames=geneKeywordAnnotationListReactive()$geneid
            )
    })
    ##
    geneEnzymeAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        geneEnzymeAnnotationList(conn,geneidSelected)
    })
    output$geneEnzymeAnnotationList<-DT::renderDataTable({
        DT::datatable(
                geneEnzymeAnnotationListReactive(),
                selection="single",
                rownames=geneEnzymeAnnotationListReactive()$geneid
            )
    })
    ##gene pathway annotation list
    genePathwayAnnotationListReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,"geneid"]
        genePathwayAnnotationList(conn,geneidSelected)
    })
    output$genePathwayAnnotationList<-DT::renderDataTable({
        DT::datatable(
                genePathwayAnnotationListReactive(),
                selection="single",
                rownames=genePathwayAnnotationListReactive()$geneid
            ) #end datatable
    })
    
    ##Gene detail
    geneInfoReactive<-reactive({
        geneidSelected<-geneListReactive()[input$geneList_rows_selected,][["geneid"]]
        geneInfo(conn,geneidSelected)
    })

    output[["geneInfo.geneid"]]<-renderText(geneInfoReactive()[["geneid"]])
    output[["geneInfo.genename"]]<-renderText(geneInfoReactive()[["genename"]])
    output[["geneInfo.genedescription"]]<-renderText(geneInfoReactive()[["genedescription"]])
    output[["geneInfo.fastaheader"]]<-renderText(geneInfoReactive()[["fastaheader"]])
    output[["geneInfo.aasequence"]]<-renderText(str_replace_all(geneInfoReactive()[["aasequence"]],"(.{60})", "\\1\n"))
    
    ## Tree
    output$tree <- renderTree({
        getTaxonomyNestedStructure(conn,1)
    })
    
    referenceProteomeIdsForTaxonomySelectedNodesReactive<-reactive({
        getReferenceProteomeIdsForTaxonomySelectedNodes( conn, input[["tree"]] )
    })

    ## Similarity
    similarityListReactive<-reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        referenceProteomeIds<-referenceProteomeIdsForTaxonomySelectedNodesReactive()
        getEuclideanInterproteomeDistances(conn,'generated_goslim', proteomeidSelected, referenceProteomeIds) ## !!
    })
    output$similarityList<-DT::renderDataTable({
        DT::datatable(
                similarityListReactive() ##,
                ## selection="single" ##,
                ## rownames=genePathwayAnnotationListReactive()$geneid
            ) #end datatable
    })
    ## Clust
##    output$euclideanDistancesDendogramPlot<-renderForceNetwork({
    output$euclideanDistancesDendogramPlot<-renderPlot({
        referenceProteomeIds<-referenceProteomeIdsForTaxonomySelectedNodesReactive()
        distancesDf<-dbGetQuery(conn, "SELECT * FROM taxogeno.euclidean_distances_keyword WHERE (proteomeid_greatest  IN ($1) OR proteomeid_least IN ($1))",params=list(referenceProteomeIds))
        distancesDf<-merge(distancesDf,dbGetQuery(conn,"select genecount as genecount_greatest, proteomeid as proteomeid_greatest, basename as basename_greatest from taxogeno.proteome where proteomeid in ($1)",list(referenceProteomeIds)), by="proteomeid_greatest")
        distancesDf<-merge(distancesDf,dbGetQuery(conn,"select genecount as genecount_least, proteomeid as proteomeid_least, basename as basename_least from taxogeno.proteome where proteomeid in ($1)",list(referenceProteomeIds)), by="proteomeid_least")
        distancesDf[,"basename_greatest"]<-gsub("^(.*)_GCA.*$","\\1",distancesDf[,"basename_greatest"])
        distancesDf[,"basename_least"]<-gsub("^(.*)_GCA.*$","\\1",distancesDf[,"basename_least"])

        distancesDf<-distancesDf[,c("basename_greatest","basename_least","distance")]
        
        uniqueBasenamesVec<-unique(c(distancesDf$basename_least,distancesDf$basename_greatest))
        reflexiveDistancesVec<-rep(0,length(uniqueBasenamesVec))
        reflexiveDistancesDf<-data.frame(basename_least=uniqueBasenamesVec,
                                         basename_greatest=uniqueBasenamesVec,
                                         distance=reflexiveDistancesVec)
        distancesDf<-rbind(distancesDf,reflexiveDistancesDf)
        distancesDfSpecular<-data.frame(basename_least=distancesDf$basename_greatest,
                                basename_greatest=distancesDf$basename_least,
                                distance=distancesDf$distance)

        distancesDf<-rbind(distancesDf,distancesDfSpecular)
        rm(distancesDfSpecular)

        scalar1 <- function(x) { x / sqrt(sum(x^2))}

        distancesDf$similarity<-rep(1,length(distancesDf$distance))-scalar1(distancesDf$distance)

        distMatrix<-xtabs(distance~basename_greatest+basename_least,distancesDf)
        
        heatmap(distMatrix)
        
        ## distancesdf[,"similarity"]<-1/distancesdf[,"distance"]
        ## distancesdf<-distancesdf[order(distancesdf$distance),]
        ## distmatrix<-xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=FALSE)
        ## #distobj<-as.dist(xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=TRUE))

        ## ##euclideanDistancesClust<-hclust(distobj,method="complete")
        
        
        ## uniqueProteomeIdVec<-unique(c(distancesdf[,"proteomeid_least"],distancesdf[,"proteomeid_greatest"]))
        ## seqIdProteomeIdMap<-data.frame(id=seq_along(uniqueProteomeIdVec)-1,proteomeid=uniqueProteomeIdVec)
        ## rm(uniqueProteomeIdVec)

        ## sourceNodesDf<-merge(seqIdProteomeIdMap, distancesdf, by.x="proteomeid",by.y="proteomeid_least")
        ## colnames(sourceNodesDf)[colnames(sourceNodesDf) == 'basename_least'] <- 'name'
        ## targetNodesDf<-merge(seqIdProteomeIdMap, distancesdf, by.x="proteomeid",by.y="proteomeid_greatest")
        ## colnames(targetNodesDf)[colnames(targetNodesDf) == 'basename_greatest'] <- 'name'
        
        ## d3ForceGraphLinksDf<-data.frame(source=c(sourceNodesDf[,"id"], targetNodesDf[,"id"]),
        ##                                 target=c(targetNodesDf[,"id"], sourceNodesDf[,"id"]),
        ##                                 value=c(sourceNodesDf[,"distance"],targetNodesDf[,"distance"])*50)
        ## d3ForceGraphLinksDf<-unique(d3ForceGraphLinksDf[order(d3ForceGraphLinksDf[,"source"]),])
        
        ## d3ForceGraphNodesDf<-unique(data.frame(id=c(sourceNodesDf[,"id"],targetNodesDf[,"id"]),
        ##                                        name=c(sourceNodesDf[,"name"],targetNodesDf[,"name"]),
        ##                                        group=rep(1,nrow(sourceNodesDf))))
        ## d3ForceGraphNodesDf<-d3ForceGraphNodesDf[order(d3ForceGraphNodesDf[,"id"]),]

        ## print(d3ForceGraphNodesDf)
        ## print(d3ForceGraphLinksDf)
        
        ## rm(sourceNodesDf)
        ## rm(targetNodesDf)

        ## forceNetwork(Links = d3ForceGraphLinksDf, Nodes = d3ForceGraphNodesDf, Source = "source", Target = "target", Value = "value", NodeID = "name", Group="group", linkDistance=JS("function(d){return d.value * 30}"), opacity = 0.8, zoom=TRUE)


        ## distancesdf[,"similarity"]<-distancesdf[,"distance"]
        ## distobj<-as.dist(xtabs(similarity~basename_greatest+basename_least,distancesdf,sparse=TRUE))
        ## euclideanDistancesClust<-hclust(distobj,method="complete")

        ## dendroNetwork(euclideanDistancesClust, zoom=TRUE)
    })
    
    ## Global annotations
    molecularFunctionGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF$golabel<-str_wrap(molecularFunctionGoslimAnnotationDF$golabel,width=15)
        molecularFunctionGoslimAnnotationDF
    })
    molecularFunctionGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF$golabel<-str_wrap(molecularFunctionGoslimAnnotationDF$golabel,width=15)
        molecularFunctionGoslimAnnotationDF
    })
    molecularFunctionGoslimAnnotationBoxplotDataReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        molecularFunctionGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, "molecular_function" ) ## !!
        molecularFunctionGoslimAnnotationDF <- molecularFunctionGoslimAnnotationDF[order(molecularFunctionGoslimAnnotationDF$goid),]
        goDescriptionDf<-unique(molecularFunctionGoslimAnnotationDF[,c("goid","goaspect","golabel")])
        data.frame(
            goid=goDescriptionDf$goid,
            goaspect=goDescriptionDf$goaspect,
            golabel=goDescriptionDf$golabel,
            Q1=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.25)$relval,
            Q2=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.50)$relval,
            Q3=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,quantile,0.75)$relval,
            IQR=aggregate(relval~goid,molecularFunctionGoslimAnnotationDF,IQR)$relval
        )
    })
    output$molecularFunctionGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                molecularFunctionGoslimAnnotationBoxplotDataReactive(),
                selection="single"
            ) #end datatable
    })
    output$molecularFunctionGoslimAnnotationBoxplot <- renderGirafe({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,"proteomeid"]
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=molecularFunctionGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=molecularFunctionGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })

    ##biologicalProcessGoslim
    biologicalProcessGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        biologicalProcessGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, c("biological_process") )
        biologicalProcessGoslimAnnotationDF$golabel<-sapply(biologicalProcessGoslimAnnotationDF$golabel,str_wrap,width=20)
        biologicalProcessGoslimAnnotationDF
    })

    biologicalProcessGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        biologicalProcessGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, c("biological_process") )
        ## biologicalProcessGoslimAnnotationDF$iqrByGoid<-sapply(
        ##         biologicalProcessGoslimAnnotationDF$goid,
        ##         function(goidParam){IQR(subset(biologicalProcessGoslimAnnotationDF, goid==goidParam)$relval)}
        ## )
        biologicalProcessGoslimAnnotationDF$golabel<-sapply(biologicalProcessGoslimAnnotationDF$golabel,str_wrap,width=20)
        biologicalProcessGoslimAnnotationDF
    })
    output$biologicalProcessGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                biologicalProcessGoslimAnnotationListReactive(),
                selection="single"
            ) #end datatable
    })
    output$biologicalProcessGoslimAnnotationBoxplot <- renderGirafe({
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=biologicalProcessGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=biologicalProcessGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })
    ##cellularComponentGoslim
    cellularComponentGoslimAnnotationPointReactive <- reactive({
        proteomeidSelected<-proteomeListDfReactive()[input$proteomeList_rows_selected,][["proteomeid"]]
        cellularComponentGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, proteomeidSelected, c("cellular_component") )
        cellularComponentGoslimAnnotationDF$golabel<-sapply(cellularComponentGoslimAnnotationDF$golabel,str_wrap,width=20)
        cellularComponentGoslimAnnotationDF
    })    

    cellularComponentGoslimAnnotationListReactive <- reactive({
        referenceProteomeIds<- referenceProteomeIdsForTaxonomySelectedNodesReactive()
        cellularComponentGoslimAnnotationDF <- getGoslimAnnotationForProteomeIdVector( conn, referenceProteomeIds, c("cellular_component") )
        ## cellularComponentGoslimAnnotationDF$iqrByGoid<-sapply(
        ##         cellularComponentGoslimAnnotationDF$goid,
        ##         function(goidParam){IQR(subset(cellularComponentGoslimAnnotationDF, goid==goidParam)$relval)}
        ## )
        cellularComponentGoslimAnnotationDF$golabel<-sapply(cellularComponentGoslimAnnotationDF$golabel,str_wrap,width=20)
        cellularComponentGoslimAnnotationDF
    })
    output$cellularComponentGoslimAnnotationList <- DT::renderDataTable({
        DT::datatable(
                cellularComponentGoslimAnnotationListReactive(),
                selection="single"
            ) #end datatable
    })
    output$cellularComponentGoslimAnnotationBoxplot <- renderGirafe({
        girafe(
            ggobj = ggplot(aes(y=relval,x=golabel),
                           data=cellularComponentGoslimAnnotationListReactive()) +
                geom_boxplot() +
                geom_point(data=cellularComponentGoslimAnnotationPointReactive(),aes(y=relval,x=golabel),color="red", size=5) +
                coord_flip(),
            options = list(opts_sizing(rescale = FALSE)),
            width_svg = 10, height_svg = 30
        )
    })

    ## user proteome upload
    tsvFilePath<-reactive({
        fileTsv <- input$tsv
        path <- fileTsv$datapath
        path
    })
    tsvFileName<-reactive({
        fileTsv <- input$tsv
        path <- fileTsv$name
        path
    })
    multifastaFilePath<-reactive({
        fileMultifasta <- input$multifasta
        path <- fileMultifasta$datapath
		path
	})
    tags<-reactive({
    	unlist(strsplit(input$tags,","))
    })
	
    jobIdOutputReactive<-eventReactive(input$sendButton,{        
        if(!is.null(input$tsv) && !is.null(input$multifasta)){
            jobId<-saveSma3sFileSet(conn, tsvFilePath(), multifastaFilePath(),tags(),isUserProteome=TRUE,webFile=TRUE,webSma3sAnnotationFileName=tsvFileName())
            jobId
        }else{
            "no proteome"
        }
    })
    output$jobId <- renderText(jobIdOutputReactive())

    ## Annotation info download
    
    output$downloadData <- downloadHandler(
       filename = function() {
         paste('data-', Sys.Date(), '.csv', sep='')
       },
       content = function(con) {
         write.csv(similarityListReactive(), con)
       })
    
}


runApp(shinyApp(ui, server), port=10081)
