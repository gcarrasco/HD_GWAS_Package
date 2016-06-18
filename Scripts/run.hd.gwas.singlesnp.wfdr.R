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

    gwas.dim=as.integer(gwas.dim)
    weight.col=as.integer(weight.col)

    table_col=gwas.dim + weight.col
    print (paste ("Weight column in intersected files:", table_col))
    print('WFDR of TFPI analysis starting')


	system(paste('intersectBed -wb -a ', gwas.results, ' -b ', prioritized.regions ,'| cut -f 4,',table_col, ' > snpList.peaks.txt', sep="" ))
    file.rename('snpList.peaks.txt',paste(prefix,'snpList.peaks.txt',sep="_"))


	snpList.peaks <- read.table(paste(prefix,'snpList.peaks.txt',sep="_"), h=F)
	colnames(snpList.peaks) <- c('snp','peak.p')

	gwas.results <- read.table(gwas.results, h=F, stringsAsFactors=F)
	colnames(gwas.results) <- c('chr','start','end','marker','beta','se','p', 'a1','a0','af','imputation','info')


	run.sfdr.wfdr <- function (gwas.results, snpList.peaks,prefix=prefix) {
			
			snpList.peaks$wgt <- sqrt(snpList.peaks$peak.p)		
	
			snpList.res <-merge(snpList.peaks, gwas.results, by.x='snp', by.y='marker', all.y=T)
			snpList.res[is.na(snpList.res$wgt), 'wgt'] <- 1 
            print("Running analysis")
			weights <- snpList.res$wgt
			weights=weights/(sum(weights)/length(weights))
			snpList.res$weight <- weights
			snpList.res$pwt <- snpList.res$p / snpList.res$weight

			snpList.res$q_fdr <- p.adjust(snpList.res$p, method="fdr")
			snpList.res$q_wfdr <- p.adjust(snpList.res$pwt, method='fdr')
			snpList.res[snpList.res$wgt>1, 'q_sfdr'] <- p.adjust(snpList.res[snpList.res$wgt>1, 'p'], method='fdr')
			snpList.res[snpList.res$wgt==1, 'q_sfdr'] <- p.adjust(snpList.res[snpList.res$wgt==1, 'p'], method='fdr')

			snpList.res <- snpList.res[order(snpList.res$p),]
			snpList.res$rank_fdr <- rank(snpList.res$q_fdr, ties='first')	
			snpList.res$rank_wfdr <- rank(snpList.res$q_wfdr, ties='first')	
			snpList.res$rank_sfdr <- rank(snpList.res$q_sfdr, ties='first')	
 
			snpList.res <- snpList.res[order(snpList.res$rank_sfdr),]
			write.table(snpList.res, file=paste(prefix,'singlesnp.wfdr.txt',sep="_"), row.names=F, col.names=T, quote=F, sep='\t')
			write.table(snpList.res[1:1000,], file=paste(prefix,'singlesnp.wfdr.top1000.txt',sep="_"), row.names=F, col.names=T, quote=F, sep='\t')
            print("Single SNP wFDR analysis: FINISHED")

	}

run.sfdr.wfdr(gwas.results, snpList.peaks,prefix=prefix)

	