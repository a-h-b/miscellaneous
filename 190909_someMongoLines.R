#getting mongo DB from dump:mongorestore --db litterDB dump/litterDB/
# run mongod: mongod --dbpath /Users/heintzbu/Documents/iDiv/Results/metaact/prok/

library(vegan)
library(RColorBrewer)
library(ape)
library(ade4)
library(gplots)
library(gclus)
library(scatterplot3d)
library(dendextend)
library(foreach)
library(robCompositions)
library(DESeq2)
library(tsne)
library(vioplot)
library(gtools)
library(MASS)
library(phyloseq)
library(car)
library(compositions)
library(colorspace)
library(VennDiagram)
library(Hmisc)
library(scales)
library(jsonlite)

source("170801_myfunctions.R")
annaEnrichPath <- function(koTab,bgList,pw=PW,pwn=pwn,koIDname="koID",pval=0.05){
  if(class(koTab)=="data.frame") sigPW <- merge(koTab,pw,by.x=koIDname,by.y="ko")
  if(class(koTab)=="matrix") sigPW <- merge(koTab,pw,by.x=0,by.y="ko")
  if(class(koTab)=="character" ) sigPW <- merge(koTab,pw,by.x=1,by.y="ko")
  sigPW <- unique(sigPW)
  colnames(sigPW)[1] <- "koID"
  colnames(sigPW)[ncol(sigPW)] <- "pw"
  
  bgSet <- merge(bgList,pw,by.x=1,by.y="ko")
  
  colnames(bgSet)[1] <- "ko"
  bgSet <- unique(bgSet)
  pwInf <- data.frame("pathway"=unique(sigPW$pw),"numberInSig","KOs","numberTotal","pval",stringsAsFactors=F)
  colnames(pwInf) <- gsub(".","",gsub("X.","",colnames(pwInf)),fixed=T)
  
  for(i in pwInf$pathway){
    q <- length(unique(sigPW$koID[sigPW$pw==i]))
    m <- length(unique(bgSet$ko[bgSet$pw==i]))
    n <- length(unique(bgSet$ko)) - m
    k <- length(unique(sigPW$koID))
    pwInf$numberInSig[pwInf$pathway==i] <- q
    pwInf$KOs[pwInf$pathway==i] <- paste(sigPW$koID[which(sigPW$pw==i)],sep=";",collapse=";")
    pwInf$numberTotal[pwInf$pathway==i] <- m
    pwInf$pval[pwInf$pathway==i] <- phyper(q-1,m,n,k,lower.tail=F)
  }
  pwInf[,6] <- p.adjust(pwInf[,5],"fdr")
  colnames(pwInf)[6] <- "adjpval"
  
  pwEn <- pwInf[ pwInf$adjpval<=pval,]
  pwEn <- merge(pwEn,pwn,by.x="pathway",by.y="pw",all.x=T)
  return(pwEn)
}


library(mongolite)
bins <- mongo("W2Ibins","litterDB")
big <- mongo("W2Icontigs","litterDB")


bins$find('{}',limit=1)

bins$index(add='{"sample":1, "MAG":1}')
bins$index(add='{"sample":1}')
bins$index(add='{"MAG":1}')
bins$index(add='{"binScore":1}')
bins$index(add='{"uniqueEss":1}')

big$find('{}',limit=1)

big$index(add='{"sample":1, "dastoolbin":1}')
big$index(add='{"dastoolbin":1}')
big$index(add='{"sample":1}')
big$index(add='{"contig":1}')
big$index(add='{"length":1}')
big$index(add='{"sample":1, "contig":1}')
big$index(add='{"coords": "2d"}')
big$index(add='{"genes.Pfam": 1}')
big$index(add='{"sample":1, "genes.gene":1}')
big$index(add='{"genes.aveCovDNA":1}')
big$index(add='{"genes.aveCovRNA":1}')
big$index(add='{"genes.KEGG":1}')
big$index(add='{"genes.dbCAN":1}')
big$index(add='{"genes.essential":1}')
big$index(add='{"krakenFamily":1}')
big$index(add='{"krakenClass":1}')
big$index(add='{"genes.readsDNA":1}')
big$index(add='{"genes.readsRNA":1}')

test <- bins$find('{"uniqueEss" : {"$gt": 90}}')
test[,c("bin","sample")]

#utility functions to run queries

toOperator <- function(criterion){
  switch(min(which(sapply(c("<=","<",">=",">","==","!=","=~","??"),
                      function(op) grepl(op,criterion)))),
  {
    ret <- list(list("$lte"=as.numeric(gsub(".+<=","",criterion))))
    names(ret)[1] <- gsub("<=.+","",criterion)
  },
  {
    ret <- list(list("$lt"=as.numeric(gsub(".+<","",criterion))))
    names(ret)[1] <- gsub("<.+","",criterion)
  },
  {
    ret <- list(list("$gte"=as.numeric(gsub(".+>=","",criterion))))
    names(ret)[1] <- gsub(">=.+","",criterion)
  },
  {
    ret <- list(list("$gt"=as.numeric(gsub(".+>","",criterion))))
    names(ret)[1] <- gsub(">.+","",criterion)
  },
  {
    ret <- list(gsub(".+==","",criterion))
    names(ret)[1] <- gsub("==.+","",criterion)
  },
  {
    ret <- list(list("$ne"=gsub(".+!=","",criterion)))
    names(ret)[1] <- gsub("!=.+","",criterion)
  },
  {
    ret <- list(list("$regex"=gsub(".+=~","",criterion),"$options"="s"))
    names(ret)[1] <- gsub("=~.+","",criterion)
  },
  {
    ret <- list(list("$exists"=1))
    names(ret)[1] <- gsub("[?]","",criterion)
  }
  )
  ret
}

toSimpleList <- function(inputV){
  l <- lapply(1:length(inputV),function(x) 1)
  names(l) <- inputV
  gsub('\\]','',
       gsub('\\[','',
            toJSON(l)))
}

toSimpleProjection <- function(inputV,noID=T){
  l <- lapply(1:length(inputV),function(x) 1)
  names(l) <- inputV
  if(noID){
    l[["_id"]] <- 0
  }else{
    l[["_id"]] <- 1
  }
  gsub('\\]','',
       gsub('\\[','',
            toJSON(l)))
}

toQuery <- function(criteria){
  if(class(criteria)!="character") stop("criteria must be entered as string")
  crits <- unlist(strsplit(criteria,split="&"))
  gsub('\\]','',
       gsub('\\[','',
                     toJSON(
                       unlist(lapply(crits,toOperator),
                              recursive = F))))
}

toOrQuery <- function(criteria){
  if(class(criteria)!="character") stop("criteria must be entered as string")
  crits <- unlist(strsplit(criteria,split="|",fixed=T))
 paste0('{"$or":[',gsub('\\]','',
      gsub('\\[','',
            toJSON(lapply(crits,toOperator))#,
                     #recursive = F))#
           )),
      ']}')
}

toAggregation <- function(stagenames,stageobjects){
  if(length(stagenames)!=length(stageobjects)) stop("Unequal lengths of arguments")
  paste0("[",
         paste(sapply(1:length(stagenames),function(i){
           if(grepl("^[{]",stageobjects[i])) {
             paste0('{"',stagenames[i],'":',stageobjects[i],'}')
           }else{
             paste0('{"',stagenames[i],'":"',stageobjects[i],'"}')
           }
         }),
         sep=",",collapse=","),
         "]")
}

matchNunwind <- function(db,query){
  db$aggregate(toAggregation(c("$match","$unwind"),
                             c(toQuery(query),
                               "$genes")))$genes
}

countSum <- function(countitem){
  paste0('{"_id": "$',
         countitem,
         '", "count": { "$sum": 1 }}')
}

binItPerTaxon <- function(taxon){
  binIterator <- bins$iterate(toQuery(paste0("taxstring=~",taxon)),
                              toSimpleList(c("MAG","sample")))
}

mgmtTablesClassified <- function(funcCat){
  contGenes <- big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                           c(toOrQuery("??dastoolbin|??krakenClass"),
                                             "$genes",toQuery(paste0("genes.",funcCat)),
                                             toSimpleProjection(c(paste0("genes.",funcCat),
                                                                  "genes.readsRNA","genes.readsDNA","sample")))))
  
  fcmgC<- tapply(contGenes$genes$readsDNA,list(contGenes$genes[[funcCat]],contGenes$sample),sum)
  fcmgC[is.na(fcmgC)] <- 0
  fcmgCN <- decostand(fcmgC,"total",2)
  fcmtC<- tapply(contGenes$genes$readsRNA,list(contGenes$genes[[funcCat]],contGenes$sample),sum)
  fcmtC[is.na(fcmtC)] <- 0
  fcmtCN <- decostand(fcmtC,"total",2)
  
  fcmgtC<- merge(fcmgC,fcmtC,by=0,all=T,suffixes=c(".mg",".mt"))
  rownames(fcmgtC) <- fcmgtC[,1]
  fcmgtC[is.na(fcmgtC)] <- 0
  fcmgtC<- as.matrix(fcmgtC[,-1])
  fcmgtCN <- decostand(fcmgtC,"total",2)
  ol <- list(fcmgCN,fcmtCN,fcmgtCN,fcmgC,fcmtC)
  names(ol) <- c("metaG.sumnorm","metaT.sumnorm","metaGT.sumnorm",
                 "metaG.counts","metaT.counts")
  ol
}


binItGeneral <- function(query){
  binIterator <- bins$iterate(query,
                              toSimpleList(c("MAG","sample")))
}



#get numbers of reads that map to any gene, to genes on kraken-contigs, bin contigs
readsDNAGenes <- vector()
readsRNAGenes <- vector()
for(s in sInf$ID){
  readsDNAGenes <- append(readsDNAGenes,
                          sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                          c(toQuery(paste0("sample==",s)),
                                                            "$genes",
                                                            toQuery("??genes.readsDNA"),
                                                            toSimpleProjection("genes.readsDNA"))))$genes$readsDNA))
  names(readsDNAGenes)[length(readsDNAGenes)] <- s
  readsRNAGenes <- append(readsRNAGenes,
                          sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                          c(toQuery(paste0("sample==",s)),
                                                            "$genes",
                                                            toQuery("??genes.readsRNA"),
                                                            toSimpleProjection("genes.readsRNA"))))$genes$readsRNA))
  names(readsRNAGenes)[length(readsRNAGenes)] <- s
}
toSimpleProjection("genes.readsRNA")
toAggregation(c("$match","$unwind","$match","$project"),
              c(toQuery(paste0("sample==",s)),
                "$genes",toQuery("??genes.readsRNA"),
                toSimpleProjection("genes.readsRNA")))

annotatedGenes <- vector()
for(s in sInf$ID){
  annotatedGenes <- append(annotatedGenes,
                           big$aggregate(toAggregation(c("$match","$unwind","$match","$count"),
              c(toQuery(paste0("sample==",s)),
                "$genes",
                toOrQuery("??genes.KEGG|??genes.Pfam|??genes.dbCAN|??genes.Resfam"),
                "count"))))
  names(annotatedGenes)[length(annotatedGenes)] <- s
}
summary(unlist(annotatedGenes))

exprKOs <- vector()
for(s in sInf$ID){
  exprKOs <- append(exprKOs,big$aggregate(toAggregation(c("$unwind","$match","$match","$group","$count"),
                                                        c("$genes",toQuery("genes.readsRNA>0"),
                                                          toOrQuery("??genes.KEGG|??genes.Pfam|??genes.dbCAN|??genes.Resfam"),
                                                          countSum("genes.KEGG"),"count")))$count)
  names(exprKOs)[length(exprKOs)] <- s
  }
summary(exprKOs)


countSum("genes.KEGG")

exprPfam2 <- vector()
for(s in sInf$ID){
  exprPfam2 <- append(exprPfam2,big$aggregate(toAggregation(c("$match","$unwind","$match","$group","$count"),
                                                            c(toQuery(paste0("sample==",s)),
                                                              "$genes",toQuery("genes.readsRNA>0&??genes.Pfam"),
                                                              countSum("genes.Pfam"),"count")))$count)
  names(exprPfam2)[length(exprPfam2)] <- s
}
summary(exprPfam2)

readsContigsDNAGenes <- vector()
readsContigsRNAGenes <- vector()
for(s in sInf$ID){
  readsContigsDNAGenes <- append(readsContigsDNAGenes,
                                 sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                                 c(toQuery(paste0("sample==",s,"&??krakenClass")),
                                                                   "$genes",
                                                                   toQuery("??genes.readsDNA"),
                                                                   toSimpleProjection("genes.readsDNA"))))$genes$readsDNA))
                                 
  names(readsContigsDNAGenes)[length(readsContigsDNAGenes)] <- s
  readsContigsRNAGenes <- append(readsContigsRNAGenes,
                                 sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                                 c(toQuery(paste0("sample==",s,"&??krakenClass")),
                                                                   "$genes",
                                                                   toQuery("??genes.readsRNA"),
                                                                   toSimpleProjection("genes.readsRNA"))))$genes$readsRNA))
  
  names(readsContigsRNAGenes)[length(readsContigsRNAGenes)] <- s
}

readsBinsDNAGenes <- vector()
readsBinsRNAGenes <- vector()
for(s in sInf$ID){
  readsBinsDNAGenes <- append(readsBinsDNAGenes,
                              sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                              c(toQuery(paste0("sample==",s,"&??dastoolbin")),
                                                                "$genes",
                                                                toQuery("??genes.readsDNA"),
                                                                toSimpleProjection("genes.readsDNA"))))$genes$readsDNA))
  readsBinsRNAGenes <- append(readsBinsRNAGenes,
                              sum(big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                                              c(toQuery(paste0("sample==",s,"&??dastoolbin")),
                                                                "$genes",
                                                                toQuery("??genes.readsRNA"),
                                                                toSimpleProjection("genes.readsRNA"))))$genes$readsRNA))
  names(readsBinsDNAGenes)[length(readsBinsDNAGenes)] <- s
  names(readsBinsRNAGenes)[length(readsBinsRNAGenes)] <- s
}

test <- big$find(toQuery("length>1000"))

test2 <- big$aggregate(toAggregation(c("$match","$unwind"),
              c(toQuery("sample==W2I01&dastoolbin==P6.2"),"$genes")))

plot(test2$genes$aveCovDNA,test2$genes$aveCovRNA,log="xy",pch=16,cex=0.3,
     col=alpha("black",rescale(log2(test2$length),0.1,0.9)))

plot(matchNunwind(big,"sample==W2I09&dastoolbin==maxbin_res.001.fasta")$aveCovDNA,
     matchNunwind(big,"sample==W2I09&dastoolbin==maxbin_res.001.fasta")$aveCovRNA,
     log="xy",pch=16,cex=0.3)

plot(matchNunwind(big,"sample==W2I09&dastoolbin==maxbin_res.001.fasta&length>10000")$aveCovDNA,
     matchNunwind(big,"sample==W2I09&dastoolbin==maxbin_res.001.fasta&length>10000")$aveCovRNA,
     log="xy",pch=16,cex=0.3)

ttax <- "Erwinia"

test <- bins$iterate(toQuery(paste0("taxstring=~",ttax)),
                     toSimpleList(c("MAG","sample")))
while(!is.null(tmp <- test$one())){
  plot(matchNunwind(big,paste0("sample==",tmp$sample,"&dastoolbin==",tmp$MAG))$aveCovDNA,
       matchNunwind(big,paste0("sample==",tmp$sample,"&dastoolbin==",tmp$MAG))$aveCovRNA,
       log="xy",pch=16,cex=0.3,
       ylab="RNA cov",xlab="DNA cov")
  mtext(paste(tmp$sample,tmp$MAG),3,1)
}
  
test <- binItPerTaxon(ttax)
while(!is.null(tmp <- test$one())){
  plot(matchNunwind(big,paste0("sample==",tmp$sample,"&dastoolbin==",tmp$MAG))$aveCovDNA,
       matchNunwind(big,paste0("sample==",tmp$sample,"&dastoolbin==",tmp$MAG))$aveCovRNA,
       log="xy",pch=16,cex=0.3,
       ylab="RNA cov",xlab="DNA cov")
  mtext(paste(tmp$sample,tmp$MAG),3,1)
}


barplot(rbind(100*readsContigsDNAGenes/readsDNAGenes,100*readsBinsDNAGenes/readsDNAGenes),beside=T,las=1,
        col=c(brewer.pal(6,"Reds")[c(4,6)]),las=2)

barplot(rbind(100*readsContigsRNAGenes/readsRNAGenes,100*readsBinsRNAGenes/readsRNAGenes),beside=T,las=1,
        col=c(brewer.pal(6,"Blues")[c(4,6)]),las=2)

barplot(rbind(100*readsContigsDNAGenes/readsDNAGenes,100*readsBinsDNAGenes/readsDNAGenes,
              100*readsContigsRNAGenes/readsRNAGenes,100*readsBinsRNAGenes/readsRNAGenes),beside=T,las=2,
        col=c(brewer.pal(6,"Reds")[c(4,6)],brewer.pal(6,"Blues")[c(4,6)]),
        space = rep(c(1,0,0.4,0),8))

pdf("mgmt_boxplot_mapping.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.2,5,0.5,0.5),mgp=c(1.1,0.2,0),tcl=-0.3)
boxplot(100*readsContigsRNAGenes/readsRNAGenes,100*readsBinsRNAGenes/readsRNAGenes,
        100*readsContigsDNAGenes/readsDNAGenes,100*readsBinsDNAGenes/readsDNAGenes,
        col=c(brewer.pal(6,"Blues")[c(4,6)],brewer.pal(6,"Reds")[c(4,6)]),horizontal=T,ylim=c(0,100),
        xlab="% mapping reads",axes=F)
box()
axis(1,cex.axis=0.7)
axis(2,las=1,labels=rep(c("classifiable\ncontigs","MAGs"),2),at=1:4,cex.axis=13/14)
Map(axis, side=2, at=3.6, col.axis=brewer.pal(6,"Reds")[5], labels="metaG", las=2,font=2,tick=F,pos=-16)
Map(axis, side=2, at=1.6, col.axis=brewer.pal(6,"Blues")[5], labels="metaT", las=2,font=2,tick=F,pos=-16)
dev.off()

dbC_list <- mgmtTablesClassified("dbCAN")


for(i in 1:nrow(sInf)){
  plot(ifelse(dbC_list[[3]][,i]>0,
              100*dbC_list[[3]][,i],1e-7),
       ifelse(dbC_list[[3]][,i+10]>0,
              100*dbC_list[[3]][,i+10],1e-6),
       pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
       col=alpha("black",0.5),cex.lab=1.4,ylim=c(1e-6,10),xlim=c(1e-7,10))
  lines(c(1e-5,100),
        c(1e-5,100),lty=3)
  axis(1,at=10^c(-7,-6:2),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
  axis(2,at=c(10^c(-6,-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
  mtext("Cazy gene rel. abundance\n",3,0.2,cex=1.4)
  mtext(paste0(" \n",sInf$ID[i]),3,0.2,cex=1.4,col=sInf$col[i],font=2)
  mtext("% metatranscriptomics",2,3.2,cex=1.4)
  box()
}

pdf("dbC_mgVsdbC_mt_CazysInContigs_scatter_mean.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,dbC_mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(dbC_mgtCN[,1:10])>min(dbC_mgCN[dbC_mgCN>0]),
            100*rowMeans(dbC_mgtCN[,1:10]),min(dbC_mgCN[dbC_mgCN>0])),
     ifelse(rowMeans(dbC_mgtCN[,11:20])>min(dbC_mtCN[dbC_mtCN>0]),
            100*rowMeans(dbC_mgtCN[,11:20]),min(dbC_mtCN[dbC_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(dbC_mgCN[dbC_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(dbC_mtCN[dbC_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("Cazy gene rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean"),3,0.2,cex=1.4,col="black",font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

pdf("dbC_mgVsdbC_mt_CazysInContigs_scatter_meanAmbient.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,dbC_mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(dbC_mgtCN[,which(sInf$treatment=="Ambient")])>min(dbC_mgCN[dbC_mgCN>0]),
            100*rowMeans(dbC_mgtCN[,which(sInf$treatment=="Ambient")]),min(dbC_mgCN[dbC_mgCN>0])),
     ifelse(rowMeans(dbC_mgtCN[,10+which(sInf$treatment=="Ambient")])>min(dbC_mtCN[dbC_mtCN>0]),
            100*rowMeans(dbC_mgtCN[,10+which(sInf$treatment=="Ambient")]),min(dbC_mtCN[dbC_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(dbC_mgCN[dbC_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(dbC_mtCN[dbC_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("Cazy gene rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean ambient"),3,0.2,cex=1.4,col=rgb(120/255,1/255,3/255),font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

pdf("dbC_mgVsdbC_mt_CazysInContigs_scatter_meanFuture_prelim.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,dbC_mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(dbC_mgtCN[,which(sInf$treatment=="Future")])>min(dbC_mgCN[dbC_mgCN>0]),
            100*rowMeans(dbC_mgtCN[,which(sInf$treatment=="Future")]),min(dbC_mgCN[dbC_mgCN>0])),
     ifelse(rowMeans(dbC_mgtCN[,10+which(sInf$treatment=="Future")])>min(dbC_mtCN[dbC_mtCN>0]),
            100*rowMeans(dbC_mgtCN[,10+which(sInf$treatment=="Future")]),min(dbC_mtCN[dbC_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(dbC_mgCN[dbC_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(dbC_mtCN[dbC_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("Cazy gene rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean future"),3,0.2,cex=1.4,col=rgb(133/255,1/255,130/255),font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

dbC_mgCRN <- decostand(t(rrarefy(t(dbC_mgC),min(colSums(dbC_mgC)))),"total",2)
dbC_mtCRN <- decostand(t(rrarefy(t(dbC_mtC),min(colSums(dbC_mtC)))),"total",2)

pcoaSuper(dbC_mgC,"dbC_mg_Cazy_cont_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(decostand(dbC_mgCRN,"pa"),"dbC_mg_Cazy_cont_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)

pcoaSuper(dbC_mtC,"dbC_mt_Cazy_cont_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(decostand(dbC_mtCRN,"pa"),"dbC_mt_Cazy_cont_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)

treat_dbC_mgC_ds <- DESeqDataSetFromMatrix(countData=dbC_mgC,
                                       colData=data.frame("T"=sInf$treatment,
                                                          stringsAsFactors=F),
                                       design=as.formula("~ T"))
treat_dbC_mgC_ds <- DESeq(treat_dbC_mgC_ds,quiet=T)
treat_dbC_mgC_dsT <- results(treat_dbC_mgC_ds,name="T_Future_vs_Ambient")
treat_dbC_mgC_dsT <- as.data.frame(treat_dbC_mgC_dsT[order(treat_dbC_mgC_dsT$padj),])


rownames(treat_dbC_mgC_dsT)[treat_dbC_mgC_dsT$padj<0.05&!is.na(treat_dbC_mgC_dsT$padj)&treat_dbC_mgC_dsT$log2FoldChange< -1]
#"CBM21.hmm"   "GH13_40.hmm" "GT69.hmm"
# ??-amylase (EC 3.2.1.1); pullulanase (EC 3.2.1.41); cyclomaltodextrin glucanotransferase (EC 2.4.1.19); cyclomaltodextrinase (EC 3.2.1.54); trehalose-6-phosphate hydrolase (EC 3.2.1.93); oligo-??-glucosidase (EC 3.2.1.10); maltogenic amylase (EC 3.2.1.133); neopullulanase (EC 3.2.1.135); ??-glucosidase (EC 3.2.1.20); maltotetraose-forming ??-amylase (EC 3.2.1.60); isoamylase (EC 3.2.1.68); glucodextranase (EC 3.2.1.70); maltohexaose-forming ??-amylase (EC 3.2.1.98); maltotriose-forming ??-amylase (EC 3.2.1.116); branching enzyme (EC 2.4.1.18); trehalose synthase (EC 5.4.99.16); 4-??-glucanotransferase (EC 2.4.1.25); maltopentaose-forming ??-amylase (EC 3.2.1.-) ; amylosucrase (EC 2.4.1.4) ; sucrose phosphorylase (EC 2.4.1.7); malto-oligosyltrehalose trehalohydrolase (EC 3.2.1.141); isomaltulose synthase (EC 5.4.99.11); malto-oligosyltrehalose synthase (EC 5.4.99.15); amylo-??-1,6-glucosidase (EC 3.2.1.33); ??-1,4-glucan: phosphate ??-maltosyltransferase (EC 2.4.99.16); 6'-P-sucrose phosphorylase (EC 2.4.1.-); amino acid transporter
#GDP-Man: ??-1,3-mannosyltransferase (EC 2.4.1.-)
#Modules of approx. 100 residues. The granular starch-binding function has been demonstrated in one case. Sometimes designated as starch-binding domains (SBD). 
rownames(treat_dbC_mgC_dsT)[treat_dbC_mgC_dsT$padj<0.05&!is.na(treat_dbC_mgC_dsT$padj)&treat_dbC_mgC_dsT$log2FoldChange> 1]
#"PL1_3.hmm"   "CBM41.hmm"   "GH13_42.hmm"
# pectate lyase, 
# ??-amylase (EC 3.2.1.1); pullulanase (EC 3.2.1.41); cyclomaltodextrin glucanotransferase (EC 2.4.1.19); cyclomaltodextrinase (EC 3.2.1.54); trehalose-6-phosphate hydrolase (EC 3.2.1.93); oligo-??-glucosidase (EC 3.2.1.10); maltogenic amylase (EC 3.2.1.133); neopullulanase (EC 3.2.1.135); ??-glucosidase (EC 3.2.1.20); maltotetraose-forming ??-amylase (EC 3.2.1.60); isoamylase (EC 3.2.1.68); glucodextranase (EC 3.2.1.70); maltohexaose-forming ??-amylase (EC 3.2.1.98); maltotriose-forming ??-amylase (EC 3.2.1.116); branching enzyme (EC 2.4.1.18); trehalose synthase (EC 5.4.99.16); 4-??-glucanotransferase (EC 2.4.1.25); maltopentaose-forming ??-amylase (EC 3.2.1.-) ; amylosucrase (EC 2.4.1.4) ; sucrose phosphorylase (EC 2.4.1.7); malto-oligosyltrehalose trehalohydrolase (EC 3.2.1.141); isomaltulose synthase (EC 5.4.99.11); malto-oligosyltrehalose synthase (EC 5.4.99.15); amylo-??-1,6-glucosidase (EC 3.2.1.33); ??-1,4-glucan: phosphate ??-maltosyltransferase (EC 2.4.99.16); 6'-P-sucrose phosphorylase (EC 2.4.1.-); amino acid transporter
#Modules of approx. 100 residues found in primarily in bacterial pullulanases. The N-terminal module from Thermotoga maritima Pul13 has been shown to bind to the ??-glucans amylose, amylopectin, pullulan, and oligosaccharide fragments derived from these polysaccharides (Lammerts van Bueren et al. (2004) Biochemistry 43:15633-42) (PMID: 15581376). 

treat_dbC_mtC_ds <- DESeqDataSetFromMatrix(countData=dbC_mtC,
                                       colData=data.frame("T"=sInf$treatment,
                                                          stringsAsFactors=F),
                                       design=as.formula("~ T"))
treat_dbC_mtC_ds <- DESeq(treat_dbC_mtC_ds,quiet=T)
treat_dbC_mtC_dsT <- results(treat_dbC_mtC_ds,name="T_Future_vs_Ambient")
treat_dbC_mtC_dsT <- as.data.frame(treat_dbC_mtC_dsT[order(treat_dbC_mtC_dsT$padj),])

rownames(treat_dbC_mtC_dsT)[treat_dbC_mtC_dsT$padj<0.05&!is.na(treat_dbC_mtC_dsT$padj)&treat_dbC_mtC_dsT$log2FoldChange>1]
#"GH19.hmm"  "GT60.hmm"  "GH132.hmm"
rownames(treat_dbC_mtC_dsT)[treat_dbC_mtC_dsT$padj<0.05&!is.na(treat_dbC_mtC_dsT$padj)&treat_dbC_mtC_dsT$log2FoldChange< -1]
#"AA1_3.hmm" "AA5_2.hmm" "CBM21.hmm" "GH12.hmm" 


dbC2_mgC<- tapply(contCAN$genes$readsDNA,list(gsub(".hmm","",gsub("_.*","",contCAN$genes$dbCAN)),contCAN$sample),sum)
dbC2_mgC[is.na(dbC2_mgC)] <- 0

summary(apply(dbC2_mgC,
              2,function(x)length(which(x>0))),
        apply(dbC2_mgC[,gsub("I0?","m",colnames(dbC2_mgC))%in%sInfo$ID[sInfo$treatment=="Future"]],
              2,function(x)length(which(x>0))))
wilcox.test(apply(dbC2_mgC[,gsub("I0?","m",colnames(dbC2_mgC))%in%sInfo$ID[sInfo$treatment=="Ambient"]],
                  2,function(x)length(which(x>0))),
            apply(dbC2_mgC[,gsub("I0?","m",colnames(dbC2_mgC))%in%sInfo$ID[sInfo$treatment=="Future"]],
                  2,function(x)length(which(x>0))))
boxplot(apply(dbC2_mgC[,gsub("I0?","m",colnames(dbC2_mgC))%in%sInfo$ID[sInfo$treatment=="Ambient"]],
              2,function(x)length(which(x>0))),
        apply(dbC2_mgC[,gsub("I0?","m",colnames(dbC2_mgC))%in%sInfo$ID[sInfo$treatment=="Future"]],
              2,function(x)length(which(x>0))))


dbC2_mgCN <- decostand(dbC2_mgC,"total",2)
dbC2_mtC<- tapply(contCAN$genes$readsRNA,list(gsub(".hmm","",gsub("_.*","",contCAN$genes$dbCAN)),contCAN$sample),sum)
dbC2_mtC[is.na(dbC2_mtC)] <- 0
dbC2_mtCN <- decostand(dbC2_mtC,"total",2)

dbC2_mgtC<- merge(dbC2_mgC,dbC2_mtC,by=0,all=T,suffixes=c(".dbC2_mg",".dbC2_mt"))
rownames(dbC2_mgtC) <- dbC2_mgtC[,1]
dbC2_mgtC[is.na(dbC2_mgtC)] <- 0
dbC2_mgtC<- as.matrix(dbC2_mgtC[,-1])
dbC2_mgtCN <- decostand(dbC2_mgtC,"total",2)

for(i in 1:nrow(sInf)){
  plot(ifelse(dbC2_mgtCN[,i]>0,
              100*dbC2_mgtCN[,i],1e-7),
       ifelse(dbC2_mgtCN[,i+10]>0,
              100*dbC2_mgtCN[,i+10],1e-6),
       pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
       col=alpha("black",0.5),cex.lab=1.4,ylim=c(1e-6,10),xlim=c(1e-7,10))
  lines(c(1e-5,100),
        c(1e-5,100),lty=3)
  axis(1,at=10^c(-7,-6:2),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
  axis(2,at=c(10^c(-6,-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
  mtext("Cazy gene rel. abundance\n",3,0.2,cex=1.4)
  mtext(paste0(" \n",sInf$ID[i]),3,0.2,cex=1.4,col=sInf$col[i],font=2)
  mtext("% metatranscriptomics",2,3.2,cex=1.4)
  box()
}

dbC2_mgCRN <- decostand(t(rrarefy(t(dbC2_mgC),min(colSums(dbC2_mgC)))),"total",2)
dbC2_mtCRN <- decostand(t(rrarefy(t(dbC2_mtC),min(colSums(dbC2_mtC)))),"total",2)

vFac <- as.numeric(as.factor(sInf$treatment))
vcol <- c(rgb(120/255,1/255,3/255),rgb(133/255,1/255,130/255))
vchar <- rep(16,2)
vid <- levels(as.factor(sInf$treatment))

pcoaSuper(dbC2_mgC,"dbC2_mg_Cazy_cont_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(dbC2_mgCRN,"dbC2_mg_Cazy_cont_rar",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(dbC2_mgCRN,"dbC2_mg_Cazy_cont_rar_byTreat",vcol,
          vchar,vFac,vid)
pcoaSuper(mgCRN,"dbC2_mg_KEGG_cont_rar_byTreat",vcol,
          vchar,vFac,vid)
pcoaSuper(decostand(dbC2_mgCRN,"pa"),"dbC2_mg_Cazy_cont_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)

pcoaSuper(dbC2_mtC,"dbC2_mt_Cazy_cont_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(dbC2_mtCRN,"dbC2_mt_Cazy_cont_rar",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(dbC2_mtCRN,"dbC2_mt_Cazy_cont_rar_byTreat",vcol,
          vchar,vFac,vid)
pcoaSuper(decostand(dbC2_mtCRN,"pa"),"dbC2_mt_Cazy_cont_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)



treat_dbC2_mgC_ds <- DESeqDataSetFromMatrix(countData=dbC2_mgC,
                                           colData=data.frame("T"=sInf$treatment,
                                                              stringsAsFactors=F),
                                           design=as.formula("~ T"))
treat_dbC2_mgC_ds <- DESeq(treat_dbC2_mgC_ds,quiet=T)
treat_dbC2_mgC_dsT <- results(treat_dbC2_mgC_ds,name="T_Future_vs_Ambient")
treat_dbC2_mgC_dsT <- as.data.frame(treat_dbC2_mgC_dsT[order(treat_dbC2_mgC_dsT$padj),])
rownames(treat_dbC2_mgC_dsT)[treat_dbC2_mgC_dsT$padj<0.05&!is.na(treat_dbC2_mgC_dsT$padj)&treat_dbC2_mgC_dsT$log2FoldChange< -1]
#"CBM21.hmm" 
#Modules of approx. 100 residues. The granular starch-binding function has been demonstrated in one case. Sometimes designated as starch-binding domains (SBD). 
rownames(treat_dbC2_mgC_dsT)[treat_dbC2_mgC_dsT$padj<0.05&!is.na(treat_dbC2_mgC_dsT$padj)&treat_dbC2_mgC_dsT$log2FoldChange> 1]
 "CBM41.hmm" 
#Modules of approx. 100 residues found in primarily in bacterial pullulanases. The N-terminal module from Thermotoga maritima Pul13 has been shown to bind to the ??-glucans amylose, amylopectin, pullulan, and oligosaccharide fragments derived from these polysaccharides (Lammerts van Bueren et al. (2004) Biochemistry 43:15633-42) (PMID: 15581376). 

treat_dbC2_mtC_ds <- DESeqDataSetFromMatrix(countData=dbC2_mtC,
                                           colData=data.frame("T"=sInf$treatment,
                                                              stringsAsFactors=F),
                                           design=as.formula("~ T"))
treat_dbC2_mtC_ds <- DESeq(treat_dbC2_mtC_ds,quiet=T)
treat_dbC2_mtC_dsT <- results(treat_dbC2_mtC_ds,name="T_Future_vs_Ambient")
treat_dbC2_mtC_dsT <- as.data.frame(treat_dbC2_mtC_dsT[order(treat_dbC2_mtC_dsT$padj),])

rownames(treat_dbC2_mtC_dsT)[treat_dbC2_mtC_dsT$padj<0.05&!is.na(treat_dbC2_mtC_dsT$padj)&treat_dbC2_mtC_dsT$log2FoldChange>1]
#"GH19.hmm"  "GT60.hmm"  "GH132.hmm"

pdf("GH19_boxplot.pdf",width=7.2/2.54,height=8/2.54,pointsize=14)
par(mar=c(1.5,3,1.2,0.5),mgp=c(1.8,0.5,0),tcl=-0.3)
boxplot(counts(treat_dbC2_mtC_ds,normalize=T)[rownames(treat_dbC2_mtC_ds)=="GH19"]~
          sInf$treatment,
        #ylim=c(0,max(metaD$lossPerc)*120),
        las=1,ylab="normalized counts",cex.lab=1.2,outline=F)
points(jitter(as.numeric(as.factor(sInf$treatment)),amount=0.1),counts(treat_dbC2_mtC_ds,normalize=T)[rownames(treat_dbC2_mtC_ds)=="GH19"],
       col=sInf$col,pch=16)
#mtext("decomposition",3,0.1,cex=1.2)
dev.off()

pdf("PKS_boxplot.pdf",width=7.2/2.54,height=8/2.54,pointsize=14)
par(mar=c(1.5,3,1.2,0.5),mgp=c(1.8,0.5,0),tcl=-0.3)
boxplot(counts(treat_mgC_ds,normalize=T)[rownames(treat_mgC_ds)=="K05551"]~
          sInf$treatment,
        #ylim=c(0,max(metaD$lossPerc)*120),
        las=1,ylab="normalized counts",cex.lab=1.2,outline=F)
points(jitter(as.numeric(as.factor(sInf$treatment)),amount=0.1),
       counts(treat_mgC_ds,normalize=T)[rownames(treat_mgC_ds)=="K05551"],
       col=sInf$col,pch=16)
#mtext("decomposition",3,0.1,cex=1.2)
dev.off()


essC_list <- mgmtTablesClassified("essential")

ess_mgC <- essC_list[[4]]
ess_mgCN <- essC_list[[1]]
ess_mtC  <- essC_list[[5]]
ess_mtCN <- essC_list[[2]]

ess_mgtCN <- essC_list[[3]]

pdf("mgVsmt_EssInContigs_scatter_mean.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(ess_mgtCN[,1:10])>min(ess_mgCN[ess_mgCN>0]),
            100*rowMeans(ess_mgtCN[,1:10]),min(ess_mgCN[ess_mgCN>0])),
     ifelse(rowMeans(ess_mgtCN[,11:20])>min(ess_mtCN[ess_mtCN>0]),
            100*rowMeans(ess_mgtCN[,11:20]),min(ess_mtCN[ess_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(ess_mgCN[ess_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(ess_mtCN[ess_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("essential genes rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean"),3,0.2,cex=1.4,col="black",font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

pdf("mgVsmt_dbEssInContigs_scatter_meanAmbient.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(ess_mgtCN[,which(sInf$treatment=="Ambient")])>min(ess_mgCN[ess_mgCN>0]),
            100*rowMeans(ess_mgtCN[,which(sInf$treatment=="Ambient")]),min(ess_mgCN[ess_mgCN>0])),
     ifelse(rowMeans(ess_mgtCN[,10+which(sInf$treatment=="Ambient")])>min(ess_mtCN[ess_mtCN>0]),
            100*rowMeans(ess_mgtCN[,10+which(sInf$treatment=="Ambient")]),min(ess_mtCN[ess_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(ess_mgCN[ess_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(ess_mtCN[ess_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("essential genes rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean ambient"),3,0.2,cex=1.4,col=rgb(120/255,1/255,3/255),font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

pdf("mgVsmt_dbEssInContigs_scatter_meanFuture_prelim.pdf",width=11.7/2.54,height=9.9/2.54,pointsize=10)
par(mar=c(3.1,4.7,3,0.5),tcl=-0.3,mgp=c(1.7,0.4,0))
plot(ifelse(rowMeans(ess_mgtCN[,which(sInf$treatment=="Future")])>min(ess_mgCN[ess_mgCN>0]),
            100*rowMeans(ess_mgtCN[,which(sInf$treatment=="Future")]),min(ess_mgCN[ess_mgCN>0])),
     ifelse(rowMeans(ess_mgtCN[,10+which(sInf$treatment=="Future")])>min(ess_mtCN[ess_mtCN>0]),
            100*rowMeans(ess_mgtCN[,10+which(sInf$treatment=="Future")]),min(ess_mtCN[ess_mtCN>0])),
     pch=16,cex=0.4,log="xy",xlab="% metagenomics",ylab="",axes=F,
     col=alpha("black",0.2),cex.lab=1.4)
lines(c(1e-5,100),
      c(1e-5,100),lty=3)
axis(1,at=c(min(ess_mgCN[ess_mgCN>0]),10^c(-6:2)),labels=c(0,format(10^c(-6:-1),scientific=T),10^c(0:2)))
axis(2,at=c(min(ess_mtCN[ess_mtCN>0]),10^c(-5:2)),las=1,labels=c(0,format(10^c(-5:-1),scientific=T),10^c(0:2)))
mtext("essential genes rel. abundance\n",3,0.2,cex=1.4)
mtext(paste0(" \nmean future"),3,0.2,cex=1.4,col=rgb(133/255,1/255,130/255),font=2)
mtext("% metatranscriptomics",2,3.2,cex=1.4)
box()
dev.off()

ess_mgCRN <- decostand(t(rrarefy(t(ess_mgC),min(colSums(ess_mgC)))),"total",2)
ess_mtCRN <- decostand(t(rrarefy(t(ess_mtC),min(colSums(ess_mtC)))),"total",2)

pcoaSuper(ess_mgC,"mg_Ess_Contigs_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(decostand(ess_mgCRN,"pa"),"mg_Ess_Contigs_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)

pcoaSuper(ess_mtC,"mt_Ess_Contigs_nonorm",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)
pcoaSuper(decostand(ess_mtCRN,"pa"),"mt_Ess_Contigs_pa",sInf$col,sInf$symbol,1:nrow(sInf),sInf$ID)

#get numbers of reads mapping to all essential genes on good contigs in the samples
ess_gfac <- colSums(ess_mgC)
ess_tfac <- colSums(ess_mtC)
ess_tgfac <- 1/(ess_tfac/ess_gfac)



#get coverages for essential genes in interesting bins
first <- T
binIt <- binItGeneral(toQuery("binScore>0.3"))
while(!is.null(it <- binIt$one())){
  tempr <- big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                       c(toQuery(paste0("sample==",it$sample,"&dastoolbin==",it$MAG)),
                                         "$genes",
                                         toQuery("genes.essential??"),
                                         toSimpleProjection(c("genes.essential","genes.aveCovDNA","genes.aveCovRNA","sample","dastoolbin")))))
  if(first){
    binEssCovs <- data.frame(tempr$genes,"sample"=tempr$sample,
                             "dastoolbin"=tempr$dastoolbin,stringsAsFactors = F)
    first <- F
  }else{
    binEssCovs <- data.frame(rbind(binEssCovs[,1:3],tempr$genes),
                             "sample"=c(binEssCovs$sample,tempr$sample),
                             "dastoolbin"=c(binEssCovs$dastoolbin,tempr$dastoolbin),
                             stringsAsFactors = F)
  }
}


binEssCovs$covRat <- binEssCovs$aveCovRNA*sapply(binEssCovs$sample,function(x) ess_tgfac[x])/binEssCovs$aveCovDNA

boxplot(log10(binEssCovs$covRat[binEssCovs$covRat>0])~binEssCovs$essential[binEssCovs$covRat>0])

allEss <- aggregate(binEssCovs$covRat[binEssCovs$aveCovRNA>0],
                    list(binEssCovs$sample[binEssCovs$aveCovRNA>0],binEssCovs$dastoolbin[binEssCovs$aveCovRNA>0]),mean)
colnames(allEss) <- c("sample","MAG","covRat")

plotTab <- merge(binTab,allEss,by.x=c("sample","MAG"),by.y=c(1:2))

plot(plotTab$GRiD[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],
     plotTab$covRat[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],
     pch=16,col=alpha(sapply(plotTab$sample[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab$binScore[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5]),
     xlim=c(1,1.7))
cor.test(plotTab$GRiD[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
    plotTab$covRat[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5])

pdf("growthVsMT_scatter_goodBins.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.3,2.2,0.5,0.5),tcl=-0.3,mgp=c(1.2,0.4,0))
plot(plotTab$GRiD[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],
     plotTab$covRat[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],
     pch=16,col=alpha(sapply(plotTab$sample[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab$binScore[plotTab$binScore>0.5&plotTab$GRiD<5&plotTab$aveCov>5]),
     xlim=c(1,1.7), xlab="growth score", ylab="essential metaT/metaG",las=1,
      cex.axis=10/14
)
dev.off()

plotTabTax <- plotTab[gsub(";s__.*","",gsub(".+g__","",plotTab$taxstring))!="",]

pdf("mTvsmG_essential_by_genus.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.3,6.9,0.5,0.5),tcl=-0.4,mgp=c(1.2,0.4,0))
boxplot(plotTabTax$covRat[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]]~
          gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]])),
        xlab="essential metaT/metaG",las=2,
        #ylim=c(1,1.7),
        varwidth=T,
        axes=F,horizontal=T,outline=F
)
points(plotTabTax$covRat[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]],
       as.numeric(as.factor(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]])))),
       pch=16,
       col=alpha(sapply(plotTabTax$sample[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]],function(x) sInf$col[sInf$ID==x]),
                 plotTabTax$binScore[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]])
)
axis(2,at=1:length(unique(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]])))),
     font=3,las=2,
     labels=levels(as.factor(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0&gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring))%in%names(table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$aveCov>3&plotTabTax$covRat>0]))))[table(gsub(";s__.*","",gsub(".+g__","",plotTabTax$taxstring[plotTabTax$binScore>0.4&plotTabTax$covRat>0&plotTabTax$aveCov>3])))>2]])))))
axis(1,cex.axis=10/14,las=1,mgp=c(1,0.2,0))
box()
dev.off()


sepEss <- aggregate(binEssCovs$covRat[binEssCovs$aveCovRNA>0],
                    list(binEssCovs$sample[binEssCovs$aveCovRNA>0],binEssCovs$dastoolbin[binEssCovs$aveCovRNA>0],
                         binEssCovs$essential[binEssCovs$aveCovRNA>0]),mean)
colnames(sepEss) <- c("sample","MAG","essential","covRat")

corE <- vector()
for(g in unique(binEssCovs$essential)){
  plotTab <- merge(binTab,sepEss[sepEss$essential==g,],by.x=c("sample","MAG"),by.y=c(1:2))
  if(length(which(plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5))>4){
  plot(plotTab$GRiD[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
       plotTab$covRat[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
       pch=16,col=alpha(sapply(plotTab$sample[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                        plotTab$binScore[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5]),
       xlab="GRiD score",ylab="MT/MG")
  mtext(g,3,0.1)
  corE <- append(corE,cor(plotTab$GRiD[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
                          plotTab$covRat[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5]))
  names(corE)[length(corE)] <- g
  }
}

allTIGR01021 <- aggregate(binEssCovs$covRat[binEssCovs$aveCovRNA>0&binEssCovs$essential=="TIGR01021"],
                          list(binEssCovs$sample[binEssCovs$aveCovRNA>0&binEssCovs$essential=="TIGR01021"],
                               binEssCovs$dastoolbin[binEssCovs$aveCovRNA>0&binEssCovs$essential=="TIGR01021"]),mean)
colnames(allTIGR01021) <- c("sample","MAG","covRat")
plotTab2 <- merge(binTab,allTIGR01021,by.x=c("sample","MAG"),by.y=c(1:2))

pdf("growthVsMT_TIGR01021_scatter_goodBins.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.3,2.2,0.5,0.5),tcl=-0.3,mgp=c(1.2,0.4,0))
plot(plotTab2$GRiD[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     plotTab2$covRat[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     pch=16,col=alpha(sapply(plotTab2$sample[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab2$binScore[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5]),
     xlim=c(1,1.7), xlab="growth score", ylab="essential metaT/metaG",las=1,
     cex.axis=10/14
)
dev.off()

allRibosomal_S11 <- aggregate(binEssCovs$covRat[binEssCovs$aveCovRNA>0&binEssCovs$essential=="Ribosomal_S11"],
                          list(binEssCovs$sample[binEssCovs$aveCovRNA>0&binEssCovs$essential=="Ribosomal_S11"],
                               binEssCovs$dastoolbin[binEssCovs$aveCovRNA>0&binEssCovs$essential=="Ribosomal_S11"]),mean)
colnames(allRibosomal_S11) <- c("sample","MAG","covRat")
plotTab2 <- merge(binTab,allRibosomal_S11,by.x=c("sample","MAG"),by.y=c(1:2))

pdf("growthVsMT_Ribosomal_S11_scatter_goodBins.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.3,2.2,0.5,0.5),tcl=-0.3,mgp=c(1.2,0.4,0))
plot(plotTab2$GRiD[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     plotTab2$covRat[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     pch=16,col=alpha(sapply(plotTab2$sample[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab2$binScore[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5]),
     xlim=c(1,1.7), xlab="growth score", ylab="essential metaT/metaG",las=1,
     cex.axis=10/14
)
dev.off()

pdf("growthVsMT_expr_scatter_goodBins.pdf",width=11.7/2.54,height=8/2.54,pointsize=14)
par(mar=c(2.3,2.2,0.5,0.5),tcl=-0.3,mgp=c(1.2,0.4,0))
plot(binTab$GRiD[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     100*binTab$expressedGenes[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5]/binTab$CDS[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
     pch=16,col=alpha(sapply(plotTab2$sample[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab2$binScore[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5]),
     xlim=c(1,1.7), xlab="growth score", ylab="% expressed genes",las=1,
     cex.axis=10/14
)
dev.off()


cor.test(binTab$GRiD[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5],
100*binTab$expressedGenes[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5]/binTab$CDS[plotTab2$binScore>0.5&plotTab2$GRiD<5&plotTab2$aveCov>5])

plot(merge(corE,cors,by=0,suffixes=c("gridCor","ratCor"))[,2:3])

hiEss <- aggregate(binEssCovs$covRat[binEssCovs$aveCovRNA>0&binEssCovs$essential %in% names(corE)[corE>0.66]],
                    list(binEssCovs$sample[binEssCovs$aveCovRNA>0&binEssCovs$essential %in% names(corE)[corE>0.66]],
                         binEssCovs$dastoolbin[binEssCovs$aveCovRNA>0&binEssCovs$essential %in% names(corE)[corE>0.66]]),mean)
colnames(hiEss) <- c("sample","MAG","covRat")

plotTab <- merge(binTab,hiEss,by.x=c("sample","MAG"),by.y=c(1:2))

plot(plotTab$GRiD[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
     plotTab$covRat[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
     pch=16,col=alpha(sapply(plotTab$sample[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],function(x) sInf$col[sInf$ID==x]),
                      plotTab$binScore[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5]))
cor.test(plotTab$GRiD[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5],
         plotTab$covRat[plotTab$binScore>0.7&plotTab$GRiD<5&plotTab$aveCov>5])




#get coverages for KEGG in Sacchari bins
first <- T
binIt <- binItGeneral(toQuery("taxstring==d__Bacteria;p__Patescibacteria;c__Saccharimonadia;o__Saccharimonadales;f__Saccharimonadaceae;g__UBA1547;s__&binScore>0.7&aveCov>5"))
while(!is.null(it <- binIt$one())){
  tempr <- big$aggregate(toAggregation(c("$match","$unwind","$match","$project"),
                                       c(toQuery(paste0("sample==",it$sample,"&dastoolbin==",it$MAG)),
                                         "$genes",toQuery("genes.KEGG??"),
                                         toSimpleProjection(c("genes.KEGG","genes.readsDNA","genes.readsRNA","sample","dastoolbin")))))
                         
  if(first){
    saccBinKOs <- data.frame(tempr$genes,"sample"=tempr$sample,"dastoolbin"=tempr$dastoolbin,
                             stringsAsFactors = F)
    first <- F
  }else{
    saccBinKOs <- data.frame(rbind(saccBinKOs[,1:3],tempr$genes),
                             "sample"=c(saccBinKOs$sample,tempr$sample),
                             "dastoolbin"=c(saccBinKOs$dastoolbin,tempr$dastoolbin),
                             stringsAsFactors = F)
  }
}
saccBinKOs$sb <- apply(saccBinKOs[,c("sample","dastoolbin")],1,function(x) paste(x,sep=" ",collapse=" "))
ko_mg_Sac <- tapply(saccBinKOs$readsDNA,list(saccBinKOs$KEGG,saccBinKOs$sb),sum)
ko_mg_Sac[is.na(ko_mg_Sac)] <- 0
ko_mt_Sac <- tapply(saccBinKOs$readsRNA,list(saccBinKOs$KEGG,saccBinKOs$sb),sum)
ko_mt_Sac[is.na(ko_mt_Sac)] <- 0

