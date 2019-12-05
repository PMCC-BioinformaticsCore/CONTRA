#USAGE Rscript wgcnv_withBAFplot.R [contraOutF] [variantOutF] [chrLenF] [outF] 


###############################
# Variant File Column Labels
###############################
varfreq.suffix=":PMCFREQ"
canonical.lbl="CANONICAL"

x=commandArgs()
args=x[match("--args",x)+1:length(x)]

f=args[1]
vf=args[2]
chrlen.f=args[3]
fout=args[4]

initial.options <- commandArgs(trailingOnly = FALSE)
file.arg.name <- "--file="
script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
script.basename <- dirname(script.name)
genomeChrLen.f <- paste(sep="/", script.basename, chrlen.f)






#setwd("C:/Jason Li 2011/Current Projects/TRWholeGeneCNV")
#f="CNATable.10rd.10bases.20bins.txt"
#vf="1308T_Ensembl_annotated.tsv"
#fout="testout"
#genomeChrLen.f="hg19_chrlen.txt"


print (f)
print (vf)
print (fout)
print(genomeChrLen.f)



dat=read.delim(f,as.is=T)

dat$Gene.Sym=sub("\\(SPLIT\\)","",dat$Gene.Sym)

gl.varfreq.col=c()
if (vf != "NA") {
	nskip=grep("^#CHROM",scan(vf,what="character",n=150,sep="\n")) -1
	v.dat=read.delim(vf,as.is=T,skip=nskip)
	samp=sub("_.*$","",basename(vf))
#	t.varfreq.col.lbl=make.names(paste(samp,"_Variant_frequency",sep=""))
	t.varfreq.col.lbl=make.names(paste(samp,varfreq.suffix,sep=""))
	t.varfreq.col=which(colnames(v.dat)==t.varfreq.col.lbl)
#	gl.varfreq.col=setdiff(grep("_Variant_frequency",colnames(v.dat)),t.varfreq.col)
	gl.varfreq.col=setdiff(grep(sub("^X","",make.names(varfreq.suffix)),colnames(v.dat)),t.varfreq.col)
}
if (length(gl.varfreq.col)>0) {
	het.snps.wh = v.dat[,canonical.lbl] == "YES" & v.dat$ID!="." & v.dat[,gl.varfreq.col] > 0.35 & v.dat[,gl.varfreq.col] < 0.65
	baf.dat=v.dat[het.snps.wh,]
	print (paste("No. of het. snps: ",nrow(baf.dat),sep=""))
	hasGL=T
} else {
	if (vf != "NA") {
		print ("WARNING: GermLine varfreq column not found")
	}
	hasGL=F
}




genes.tmp1=unique(dat$Gene.Sym)
genes.tmp2=strsplit(genes.tmp1,";")
genes.tmp3=unlist(genes.tmp2)

genes.all=unique(genes.tmp3)


getStat=function(dat.ext) {#,note,gene) {

	lrs=dat.ext$Mean.of.LogRatio
	lrs.mean=mean(lrs)
	lrs.median=median(lrs)
	lrs.sd=sd(lrs)
	n=length(lrs)
	n2=length(unique(lrs))
	if (n>2 && lrs.mean>mu && n2>1) {
		pval=t.test(lrs,alternative="greater",mu=0.2)$p.value
		gainloss="gain"
		percentage.probes.beyond.threshold=sum(lrs>=mu) / n
	} else if (n>2 && lrs.mean < -mu && n2>1) {
		pval=t.test(lrs,alternative="less",mu=-0.2)$p.value
		gainloss="loss"
		percentage.probes.beyond.threshold=sum(lrs<= -mu) / n
	} else {
		pval=1
		gainloss="."
		percentage.probes.beyond.threshold="."
	}
	
	#u.test.pval=wilcox.test(lrs)$p.value
	ch=dat.ext$Chr[1]
	s1=min(dat.ext$OriStCoordinate)
	e1=max(dat.ext$OriEndCoordinate)
	
	#outdf=data.frame(Chr=ch,Start=s1,End=e1,mean.LR=lrs.mean,median.LR=lrs.median,sd.LR=lrs.sd,GainLoss=gainloss,pval,percentage.probes.beyond.threshold,u.test.pval,stringsAsFactors=F)
	outdf=data.frame(Chr=ch,Start=s1,End=e1,mean.LR=lrs.mean,median.LR=lrs.median,sd.LR=lrs.sd,GainLoss=gainloss,pval,n,percentage.probes.beyond.threshold,stringsAsFactors=F)
	
}


#outdf=data.frame()

app=F
mu=0.2
fo=paste(fout,"_mu",mu,".txt",sep="")

for ( i in 1:length(genes.all)) {

	#print(i)
	g=genes.all[i]
	wh=grep(paste("\\b",g,"\\b",sep=""),dat$Gene.Sym)
	dat.ext=dat[wh,]
	
	tmpdf=getStat(dat.ext)
	tmpdf1=data.frame(Gene=g,INFO="WholeGene",tmpdf,stringsAsFactors=F)
		
	
	#outdf=rbind(outdf,tmpdf1)
	write.table(tmpdf1,fo,col.names=!app,append=app,sep="\t",row.names=F,quote=F)
	app=T

}


outdf=read.delim(fo,as.is=T)
outdf$adj.pval=p.adjust(outdf$pval,method="BH")
write.table(outdf,fo,col.names=T,append=F,sep="\t",row.names=F,quote=F)


chrlen=read.delim(genomeChrLen.f,header=F,as.is=T)

chrs=chrlen$V1


dat$chr=sub("chr","",dat$Chr)
outdf$chr=sub("chr","",outdf$Chr)



plotBAF=function(ch,ch.len) {

### [CONTHERE] 130228 ----  PLOT baf.dat[,t.varfreq.col]  versus baf.dat$POS
	baf.dat.ch=baf.dat[baf.dat$X.CHROM==ch,]
	
	xs=baf.dat.ch$POS/ch.len
	ys=baf.dat.ch[,t.varfreq.col]
	plot(xs,ys,xlim=c(-0.05,1.05),ylim=c(0,1),type="p",pch=20,col=rgb(0.2,0.2,0.2,0.25),xlab="",ylab="BAF",axes=F)
	axis(2)
	
}

plotLR=function(ch,ch.len) {

	x.at= seq(0,1,0.1)
	x.lbls=paste(round(x.at*ch.len/1000000,0),"M",sep="")

	dat.ch=dat[dat$chr==ch,]
	xs=((dat.ch$OriEndCoordinate+dat.ch$OriStCoordinate)/2)/ch.len
	ys=dat.ch$Mean.of.LogRatio

	plot(xs,ys,xlim=c(-0.05,1.05),ylim=c(-3,3),type="p",pch=20,col=rgb(0.2,0.2,0.2,0.05),main=paste("chr",ch,sep=""),xlab="",ylab="LogFC",axes=F)
	axis(2)
	axis(1,at=x.at,labels=x.lbls,las=2,cex.axis=0.8)


	outdf.ch=outdf[outdf$chr==ch,]
	outdf.ch.sig=outdf.ch[outdf.ch$adj.pval < 0.05 & outdf.ch$n > 10 & outdf.ch$percentage.probes.beyond.threshold > 0.5,]

	if (nrow(outdf.ch.sig) > 0) {
	#xs=outdf.ch.sig$Start / ch.len
	#ys=outdf.ch.sig$End / ch.len
	xs=(outdf.ch.sig$Start + outdf.ch.sig$End)/2 / ch.len
	ys=outdf.ch.sig$mean.LR
	points(xs,ys,col="red",bg="red",pch=22)

	text(xs,-3,labels=outdf.ch.sig$Gene,cex=0.8,col="red",adj=0,srt=90)
	}
	abline(h=0,col=rgb(0.5,0.5,0.5,0.3))
	abline(h=-0.3,col=rgb(0.7,0.7,0.2,0.3))
	abline(h=0.3,col=rgb(0.7,0.7,0.2,0.3))


}
plotAll=function(doPlotBAF=F) {
	for (i in 1:length(chrs)) {

		ch=chrs[i]
		ch.len=chrlen$V2[i]
		
		if (doPlotBAF && hasGL) {
			nf1=layout(matrix(1:2,nrow=2,byrow=T),heights=c(3,2))
			m=par("mar")
			m[1]=4.1
			m[3]=3.1
			par("mar"=m)
		}

		plotLR(ch,ch.len)
		
		
		if (doPlotBAF && hasGL) {
			m[1]=1.1
			m[3]=0.1
			par("mar"=m)

			plotBAF(ch,ch.len)
		}
	}
}
png(paste(sub(".txt$","",fo),".png",sep=""),width=1200,height=800,bg="white",type="cairo")

nf=layout(matrix(1:25,nrow=5,byrow=T))
par("mai"=c(0.4,0.25,0.1,0.05))

plotAll()
dev.off()


png(paste(sub(".txt$","",fo),"_byChr_%02d.png",sep=""),width=1200,height=800,bg="white",type="cairo")
plotAll(doPlotBAF=T)
dev.off()



pdf(paste(sub(".txt$","",fo),"_sigGenes.pdf",sep=""))
outdf.sig=outdf[outdf$adj.pval < 0.05 & outdf$n > 10 & outdf$percentage.probes.beyond.threshold > 0.5,]

for ( i in 1:min(200,nrow(outdf.sig))) {
	g=outdf.sig$Gene[i]
	wh=grep(paste("\\b",g,"\\b",sep=""),dat$Gene.Sym)
	dat.ext=dat[wh,]

	xlim=c(min(dat.ext$OriStCoordinate),max(dat.ext$OriStCoordinate))
	ylim.lr=c(min(dat.ext$Mean.of.LogRatio),max(dat.ext$Mean.of.LogRatio))
	ylim.doc=c(min(dat.ext[,c("tumour.rd","normal.rd")]),max(dat.ext[,c("tumour.rd","normal.rd")]))

	nf=layout(matrix(c(1,2),nrow=2))

	plot(0,type="n",xlim=xlim,ylim=ylim.lr,xlab="coord",ylab="Mean of LogRatio",main=g)
	xs=(dat.ext$OriStCoordinate+dat.ext$OriEndCoordinate)/2
	ys=dat.ext$Mean.of.LogRatio
	lines(xs,ys,lwd=2,col="blue")

	plot(0,type="n",xlim=xlim,ylim=ylim.doc,xlab="coord",ylab="Read Depth (adjusted)")
	#xs=(dat.ext$OriStCoordinate+dat.ext$OriEndCoordinate)/2
	ys=dat.ext$tumour.rd
	lines(xs,ys,lwd=1,col="red")
	ys=dat.ext$normal.rd
	lines(xs,ys,lwd=1,col="green")

}
dev.off()
