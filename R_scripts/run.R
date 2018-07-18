source("internal.R")
para = commandArgs()
res = pmatch("--target", para)
if (is.na(res)) { print("no target provided"); q("no") }
foldername = strsplit(para[res], "=")[[1]][2]
res = pmatch("--config", para)
if (is.na(res)) { print("no target provided"); q("no") }
con = strsplit(para[res], "=")[[1]][2]
#foldername="Example"
#con=file.path(paste(foldername, "/config.txt", sep=""))
r = readLines(con)
get <- function(type){i = grep(type,r); strsplit(r[i], "=")[[1]][2]}
as.v <- function(ch){as.numeric(strsplit(ch,",")[[1]])}
model=get("model")
{if (model=="Gaussian(prec)"){
  xlim = as.v(get("xlim"))
  ylim = as.v(get("ylim"))
  histbins = as.v(get("histbins"))
  histvalues = as.v(get("histvalues"))
  rseq = as.v(get("rseq"))
  thseq = as.v(get("thseq"))
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

o = order(histbins); histbins=histbins[o]; histvalues=histvalues[o]
f = approxfun(histbins, histvalues, yleft=histvalues[1], yright=histvalues[length(histvalues)])
cst=integrate(f, lower=histbins[o],upper=histbins[length(histbins)])$value
psd <- function(sd){
  log(f(sd))-log(cst)
}
minsd = histbins[1]; maxsd = histbins[length(histbins)]

ld=list.dirs(foldername, recursive=FALSE)
ld=ld[ld!=foldername]

sapply(file.path(ld), function(foldername){

rsequence = seq(rseq[1], rseq[2], by=rseq[3])
thsequence = seq(thseq[1], thseq[2], by=thseq[3])

data= read.csv(file.path(paste(foldername, "/data.txt", sep="")))
pts = data[,1:2]; sds = data[,3];
labdir = file.path(foldername)
dir.create(labdir, showWarnings = FALSE)

res=Kclust(pts=pts, sds=sds, xlim=xlim, ylim=ylim, psd=psd, minsd=minsd, maxsd=maxsd,
    rseq=rsequence, thseq=thsequence, labdir=labdir,
    useplabel=useplabel, alpha=alpha, pb=pb, score=TRUE, rlabel=TRUE, report=TRUE)
writeRes(res, file.path(paste(foldername, "/r_vs_thresh.txt", sep="")), file.path(foldername))

})
