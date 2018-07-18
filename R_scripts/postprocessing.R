source("internal.R")

makeplot=FALSE
skeleton=FALSE

para = commandArgs()
res = pmatch("--folders", para)
if (is.na(res)) { print("no folders provided"); q("no") }
foldernames = strsplit(strsplit(para[res], "=")[[1]][2], ":")[[1]]

getBestLabels <- function(clusterscale, thresh, foldername) {
    zipfile = paste(foldername, "all_labels.txt.gz", sep="/")
    gz = gzfile(zipfile, "rt")
    labelData = scan(gz, what="")
    close(gz)
    scaleString = paste("clusterscale", clusterscale, sep="")
    threshString = paste("thresh", thresh, sep="")
    exploded = strsplit(labelData, ",")
    bestLabels = list()
    for (condition in exploded) {
        if (condition[1] == scaleString && condition[2] == threshString) {
            bestLabels = condition
            break
        }
    }
    l = length(bestLabels)
    bestLabels[3:l]
}

sapply(foldernames, function(expname){
  nexpname=expname
  if (skeleton){
    nexpname=paste("R_", expname, sep="")
    dir.create(file.path(nexpname))
    file.copy(file.path(paste(expname, "/config.txt", sep="")), paste(nexpname, "/config.txt", sep=""))
    file.copy(file.path(paste(expname, "/sim_params.txt", sep="")), paste(nexpname, "/sim_params.txt", sep=""))
  }

  r = readLines(con=file.path(paste(nexpname, "/config.txt", sep="")))
  get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
  as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
  model=get("model")
  {if (model=="Gaussian(prec)"){
    xlim = as.v(get("xlim"))
    ylim = as.v(get("ylim"))
    histbins = as.v(get("histbins"))
    histvalues = as.v(get("histvalues"))
    if (length(grep("pbackground",r))==0 | length(grep("alpha",r))==0){
      useplabel=FALSE; pb=NULL; alpha=NULL
    }
    else {
      useplabel=TRUE;
      pb=as.numeric(get("pbackground"))
      alpha=as.numeric(get("alpha"))
    }
  }
  else {stop("Haven't implemented anything else!")}}

  all=list.files(expname)
  dirnames=all[file.info(file.path(paste(expname, "/", all, sep="")))$isdir]

res=sapply(dirnames, function(dirname){
  foldername=file.path(paste(expname, "/", dirname, sep=""))
  nfoldername=file.path(paste(nexpname, "/", dirname, sep=""))
  if (skeleton){
    dir.create(nfoldername)
    file.copy(file.path(paste(foldername, "/data.txt", sep="")), file.path(paste(nfoldername, "/data.txt", sep="")))
  }
  data=read.csv(file.path(paste(nfoldername, "/data.txt", sep="")))

  pts = data[,1:2]; sds = data[,3];
  if (skeleton){
    file.copy(file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(paste(nfoldername, "/r_vs_thresh.txt", sep="")))
  }
  r = read.csv(file.path(paste(nfoldername, "/r_vs_thresh.txt",sep="")), header=FALSE, sep="\t")

  m = as.matrix(r)
  cs=(m[1,])[-1]
  thr=(m[,1])[-1]
  m = m[2:length(m[,1]),2:length(m[1,])]
  which.maxm <- function(mat){
    indcol <- rep(1:ncol(mat), each=nrow(mat))[which.max(mat)]
    indrow <- rep(1:nrow(mat), ncol(mat))[which.max(mat)]
    c(indrow, indcol)
  }
  best=which.maxm(m)
  bestcs=cs[best[2]]
  bestthr=thr[best[1]]
  bfile=file.path(paste(foldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
  nbfile=bfile
  if (skeleton){
    dir.create(paste(nfoldername, "/labels", sep=""))
    nbfile=file.path(paste(nfoldername, "/labels/clusterscale", bestcs, " thresh", bestthr, "labels.txt", sep=""))
    file.copy(bfile, nbfile)
  }
  #labelsbest = strsplit(readLines(nbfile),",")[[1]]
  labelsbest = getBestLabels(bestcs, bestthr, foldername)
  #cat(length(labelsbest), "\n")


##Some summaries
  wfile=file.path(paste(nfoldername, "/summary.txt", sep=""))
  cat("The best: clusterscale", bestcs, " thresh", bestthr, "labels.txt\nNumber of clusters:", nClusters(labelsbest),
        "\nPercentage in clusters: ", percentageInCluster(labelsbest),
        "%\nMean number of molecules per cluster: ", nMolsPerCluster(labelsbest),
        "\nMean radius: ", mean(clusterRadii(pts, labelsbest)), "\n", sep="", file=wfile)
  if (makeplot){
    if (dim(data)[2]==4){
      labelstrue=sapply(as.numeric(data[,4]), function(n){if (n==0) paste(runif(1)) else {paste(n)}})
      pdf(file.path(paste(nfoldername, "/plot.pdf", sep="")))
      par(pty="s")
      par(mfrow=c(1,2))
      par(mar=c(4,4,.5, .5))
      par(oma=c(1,1,1,1))
      plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelstrue), sub="True labels", xlab="",ylab="")
                                        #X11()
      par(mar=c(4,4,.5, .5))
      plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Estimated",xlab="",ylab="")
      dev.off()
    }
    else {
      pdf(file.path(paste(nfoldername, "/plot.pdf", sep="")))
      par(pty="s")
      par(mar=c(4,4,.5, .5))
      plot(pts, xlim=xlim, ylim=ylim, col=mkcols(labelsbest), sub="Clustering",xlab="",ylab="")
      dev.off()
    }
  }

list(radii=clusterRadii(pts, labelsbest), nmols=molsPerCluster(labelsbest), nclusters=nClusters(labelsbest), percentageclustered=percentageInCluster(labelsbest))
})


nmols=c()
for (i in 1:(length(res)/4)){
  nmols=c(nmols, res[[4*(i-1)+2]])
}
#h=hist(nmols, breaks=seq(0,max(nmols)*1.5, length=10), plot=FALSE)
#pdf(file.path(paste(nexpname, "/nmols.pdf", sep="")))
#plot(h, xlab="Number of molecules", ylab="Number of clusters", main="")
#dev.off()
f=file.path(paste(nexpname, "/nmols.txt", sep="")); cat(nmols, file=f, sep=","); cat("\n", file=f, append=TRUE)

radii=c()
for (i in 1:(length(res)/4)){
  radii=c(radii, res[[4*(i-1)+1]])
}
#h=hist(radii, breaks=seq(0,max(radii)*1.5, length=10), plot=FALSE)
#pdf(file.path(paste(nexpname, "/radii.pdf", sep="")))
#plot(h, xlab="Cluster radius", ylab="Number of clusters", main="")
#dev.off()
f=file.path(paste(nexpname, "/radii.txt", sep="")); cat(radii, file=f, sep=","); cat("\n", file=f, append=TRUE)

})
