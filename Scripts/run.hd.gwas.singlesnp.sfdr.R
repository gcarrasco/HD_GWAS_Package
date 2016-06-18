args <- commandArgs(trailingOnly = TRUE)


gwas.results <- args[1]
print (paste ("GWAS input file:", gwas.results))
colnames=c('chr','start','end','marker','beta','se','p', 'a1','a0','af','imputation','info')

prioritized.regions <- args[2]
print (paste ("Priority input file:", prioritized.regions ))

out_dir <-args[3]
print (paste ("Output Directory:", out_dir))

prefix <-args[4]
print (paste("Prefix:", prefix))

program.directory<- args[5]
rprograms<-paste(program.directory, "r_scripts", sep="/")

task <- args[6]
print (paste("Run HD-GWAS SNP:", task))


weight.col <- args [7]
weight.col=as.integer(weight.col)
print (paste ("Weight column:",weight.col))

gwas.dim <- args [8]

setwd(out_dir)

system(paste('intersectBed -wb -a ', gwas.results, ' -b ', prioritized.regions ,'| cut -f 4 > snpList.txt', sep="" ))

file.rename('snpList.txt',paste(prefix,'snpList.txt',sep="_"))


	snpList <- read.table(paste(prefix,'snpList.txt',sep="_"), h=F)
	colnames(snpList) <- c('snp')
	snpList$prior <- 1		

	gwas.results <- read.table(gwas.results, h=F, stringsAsFactors=F)
	colnames(gwas.results) <- c('chr','start','end','marker','beta','se','p', 'a1','a0','af','imputation','info')


	run.sfdr <- function (gwas.results, snpList.peaks,prefix=prefix) {

			snpList.res <-merge(snpList, gwas.results, by.x='snp', by.y='marker', all.y=T)
			snpList.res[is.na(snpList.res$prior), 'prior'] <- 0 
            print("Running analysis")
			snpList.res$q_fdr <- p.adjust(snpList.res$p, method="fdr")
			snpList.res[snpList.res$prior==1, 'q_sfdr'] <- p.adjust(snpList.res[snpList.res$prior==1, 'p'], method='fdr')
			snpList.res[snpList.res$prior==0, 'q_sfdr'] <- p.adjust(snpList.res[snpList.res$prior==0, 'p'], method='fdr')
			snpList.res <- snpList.res[order(snpList.res$p),]
			snpList.res$rank_fdr <- rank(snpList.res$q_fdr, ties='first')	
			snpList.res$rank_sfdr <- rank(snpList.res$q_sfdr, ties='first')
			snpList.res <- snpList.res[order(snpList.res$rank_sfdr),]
			write.table(snpList.res, file=(paste(prefix,'singlesnp.sfdr.txt',sep="_")), row.names=F, col.names=T, quote=F, sep='\t')
			write.table(snpList.res[1:1000,], file=(paste(prefix,'singlesnp.sfdr.top1000.txt',sep="_")), row.names=F, col.names=T, quote=F, sep='\t')
            print("Single SNP sFDR analysis: FINISHED")

	}

	run.sfdr(gwas.results, snpList.peaks,prefix=prefix)
		
		
	