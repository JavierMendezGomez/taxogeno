source("taxogeno.R")
library(DBI)
library(tools)

dbConn <- dbConnect(RPostgres::Postgres(), dbname="taxogeno")

updateNcbiAssemblySummaryGenbank(dbConn)
updateBiosqlNcbiTaxonomy(dbConn)
updateUniprotKeyword(dbConn)
updateGeneOntology(dbConn)

dbDisconnect(dbConn)
