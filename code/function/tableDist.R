###################################################
## Jared Studyvin
## 26 July 2016
## 998.03
## table the distribution percentage of the top choice from the simulation
###################################################

tableDist <- function(files,path){

    names(files) <- gsub('.csv','',files)

    datList <- llply(files,function(x,p){
                         dat <- read.csv(paste0(p,x))
                         if(grepl('WL',x)){
                             dat <- ddply(dat,~rep,function(x){x[1,]})
                         }
                         ##out <- with(dat,table(distn))
                         out <- ddply(dat,~distn,summarize,proportion=length(distn)/nrow(dat))
                         return(out)
                     },p=path)


    ##head(datList[[1]])

    summaryOut <- ldply(datList,rbind,.id='seedDistn')
    summaryOut <- summaryOut %>%mutate(distn=ifelse(is.na(distn),'notConverge',as.character(distn)))

    out <- ldply(dlply(summaryOut,~seedDistn,function(x){spread(x,key=distn,value=proportion)}),rbind)


    return(out)
} #end tableDist function
