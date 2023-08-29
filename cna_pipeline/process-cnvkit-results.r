# process cnvkit results

library(reshape2)
library(DNAcopy)

args = commandArgs(T)

cnvkit.out.dir = args[1] #Contains target, antitarget,fixed, segmented files for sample
sample = args[2] #Sample name
ref.file = args[3] #Reference cnn
genome = args[4] #hg19 or hg38
target.genes = args[5]


##### chromosome size
# prepare chromosome position for plotting all chromsomes contiguously in one plot
# Currently doesn't allow user input chr sizes
# Defaults to hg19 unless hg38 is specified specifically
if(genome == 'hg38'){
    size = read.table('/usr/bin/hg38.chr.size.txt', stringsAsFactors=F)
}else{
    size = read.table('/usr/bin/hg19.chr.size.txt', stringsAsFactors=F) 
}
colnames(size) = c('chr', 'len')
chrs = paste0('chr', c(1:22, 'X', 'Y'))
size = size[size$chr %in% chrs,]
rownames(size) = size$chr; size = size[chrs,]
n = length(chrs)
size$start = 1.0
for (i in 1:(n-1)){
    chrAfter = (i+1):n
    chlen = size$len[size$chr == chrs[i]]
    size$start[chrAfter] = size$start[chrAfter] + chlen
}
size$end = size$start + size$len - 1;
size$pos = (size$start + size$end)/2


# read targeted genes
# targeted gene infomation
g = read.table(target.genes, header=T, stringsAsFactors=F, sep='\t')
colnames(g) = c("gene.chr", "gene.start", "gene.stop", "gene", "description")

#key.cn.genes = unlist(strsplit('AR,PTEN,TP53,MYC,RB1,ERG,TMPRSS2,APC,CDKN1B,CHD1', ','))
#g$geneid = paste0(ifelse(g$gene %in% key.cn.genes, '*', ''), g$gene)


# loop thru cnvkit output and aggregate results
cnres = NULL
res = NULL

    plasma = sample
    control = ref.file
    resDir = cnvkit.out.dir
    plaCovDir = resDir
    ctlCovDir = resDir
    pooled = T

    # read plasma cov and log2r
    pla.tar.cov.file = file.path(plaCovDir, paste0(plasma, '.targetcoverage.cnn'))
    pla.offtar.cov.file = file.path(plaCovDir, paste0(plasma, '.antitargetcoverage.cnn'))
    x = read.table(pla.tar.cov.file, header=T, stringsAsFactors=F, sep='\t')
    y = read.table(pla.offtar.cov.file, header=T, stringsAsFactors=F, sep='\t')
    x = rbind(x,y)
    colnames(x)[(ncol(x)-1):ncol(x)] = paste0('plasma.', colnames(x)[(ncol(x)-1):ncol(x)])
    pla.cov = x #Plasma coverage
    
    # read reference coverage
    x = read.table(ref.file, header=T, stringsAsFactors=F, sep='\t')
    cols = colnames(x) %in% c('depth', 'log2', 'weight')
    colnames(x)[cols] = paste0('ref.', colnames(x)[cols])
    ref.cov = x #Reference coverage

    # read segmentation call
    cns.file = file.path(resDir, paste0(plasma, '.cns'))
    x = read.table(cns.file, header=T, stringsAsFactors=F, sep='\t')
    x$name = paste0(plasma, '/', x$chromosome, '_', x$start, '_', x$end, ':', x$probes)
    cns = x

    # read cna ratio and normalized/bias-corrected depth by cnvkit
    cnr.file = file.path(resDir, paste0(plasma, '.cnr'))
    x = read.table(cnr.file, header=T, stringsAsFactors=F, sep='\t')
    cols = colnames(x) %in% c('depth', 'log2', 'weight')
    colnames(x)[cols] = paste0('cnr.', colnames(x)[cols])
    pla.cnr = x
    z = merge(ref.cov, pla.cov)
    z = merge(z, pla.cnr)

    z = merge(z, g, all.x=T)
    z$desc[which(z$gene == 'Antitarget')] = 'Antitarget'
    z$desc[is.na(z$desc)] = 'Unknown'
    z$targeted = 'U'
    z$targeted[which(!(z$gene %in% c('Antitarget', '-')))] = 'Yes'
    z$targeted[which(z$gene == 'Antitarget')] = 'No'
    z$targeted[which(z$gene == '-')] = '?'

    # rescale targeted bins' log2r using CN control genes
    z$cnr.log2r.ctl = z$cnr.log2
    cn.control = z$description == 'CN-control' & !is.na(z$description)
    sex.chrs = c('chrX', 'chrY')
    # non-sex chrs
    ctrl.median.log2r = median(z$cnr.log2[cn.control & !(z$chromosome %in% sex.chrs)])
    sel = z$gene != 'Antitarget' & !(z$chromosome %in% sex.chrs)
    z$cnr.log2r.ctl[sel] = z$cnr.log2r.ctl[sel] - ctrl.median.log2r
    # sex chrs
    for (cch in sex.chrs){
        if(cch == 'chrY'){next} # no target on chrY, skip
        ctrl.median.log2r = median(z$cnr.log2[cn.control & z$chromosome==cch])
        sel = z$gene != 'Antitarget' & z$chromosome==cch
        z$cnr.log2r.ctl[sel] = z$cnr.log2r.ctl[sel] - ctrl.median.log2r
    }

    #resegmentation using cn control corrected
    z$pos = (z$start + z$end)/2
    cna.obj = CNA(cbind(z$cnr.log2r.ctl), z$chromosome, z$pos,
                  data.type='logratio', sampleid=plasma)
    smoothed.cna.obj = smooth.CNA(cna.obj)
    #cna.seg = segment(smoothed.cna.obj)
    cna.seg = segment(cna.obj)
    cna.seg = cna.seg$output
    cna.seg$chrom = as.character(cna.seg$chrom)
    cseg = cna.seg[, c('chrom', 'loc.start', 'loc.end', 'num.mark', 'seg.mean')]
    cseg$loc.start = as.integer(cseg$loc.start)
    cseg$loc.end = as.integer(cseg$loc.end)
    cseg$name = paste0(plasma, '/', cseg$chrom, '_', cseg$loc.start, '_',
        cseg$loc.end, ':', cseg$num.mark)
    cseg = cseg[, c('chrom', 'loc.start', 'loc.end', 'name', 'seg.mean')]
    write.table(cseg, file=paste0(plasma, '.bed'), sep='\t',
        quote=F, row.names=F, col.names=F)
    write.table(z, file=paste0(plasma, '.cov.tsv'), sep='\t',
        quote=F, row.names=F, col.names=T)

    # prepare seg data for plotting
    #cs = cseg
    #colnames(cs)[1:3] = c('chromosome', 'start', 'end')
    #cs = prepChr(cs)


    #key.genes = unlist(strsplit('AR,PTEN,TP53,MYC,RB1,ERG,TMPRSS2,APC,ATM,AXL,BRAF,CDKN1B,CHD1,ETV1,ETV2,ETV4,ETV5,FOXA1,KRAS,PIK3CA', ','))
    #z$color='black'; z$color[z$gene %in% key.genes] = 'red'
    #z$gene[!(z$gene %in% c('Antitarget', key.genes))] = '0ther'
    #z$rmask[is.na(z$rmask)] = 0
    #colrs = c('AR'='red', 'TP53'='darkgreen', 'PTEN'='blue', 'APC'='purple',
    #          'Antitarget'='gray', '0ther'='cyan')

    # plottttt
    #r1 = with(z[z$gene == 'Antitarget',], cor(cnr.log2, ref.depth, method='spearman'))
    #r2 = with(z[z$gene != 'Antitarget',], cor(cnr.log2, ref.depth, method='spearman'))
    #x1 = mean(log2(z[z$gene == 'Antitarget', 'ref.depth']))
    #x2 = mean(log2(z[z$gene != 'Antitarget', 'ref.depth']))
    #ttl = paste0('plasma = ', plasma, ' / control = ', control)
    #resi = data.frame(plasma=plasma, control=control, offtar.correl=r1,
    #                  ontar.correl=r2, stringsAsFactors=F)
    #if(is.null(res)){res = resi}else{res = rbind(res,resi)}

    #pdf(paste0(out.qc.dir, '/', plasma, '.pdf'), width=8, height=5)
    #par(mfrow=c(1,1))
    #plot(log2(z$ref.depth), z$cnr.log2, col=colrs[z$gene], main=ttl,
    #     xlab='log2(control depth)', ylab='log2ratio(plasma vs. control)')
    #abline(h=0, col='red', lty='dotted')
    #text(x=x1,y=1.5, paste0('Off-target\ncor = ', sprintf('%0.3g', r1)))
    #text(x=x2,y=1.5, paste0('On-target\ncor = ', sprintf('%0.3g', r2)))
    #plot(z$cnr.weight, z$cnr.log2, col=colrs[z$gene], main=ttl,
    #     xlab='bin weight', ylab='log2ratio(plasma vs. control)')
    #abline(h=0, col='red', lty='dotted')
    #plot(z$gc, log2(z$ref.depth), col=colrs[z$gene], main=ttl,
    #     xlab='GC content', ylab='log2(control depth)')
    #plot(z$rmask, log2(z$ref.depth), col=colrs[z$gene], main=ttl,
    #     xlab='Repeat masked content', ylab='log2(control depth)')
    #par(mfrow=c(1,4))
    #boxplot(log2(ref.depth) ~ gene, data=z, horizontal=T, las=2, cex.axis=0.5,
    #        main='Control', xlab='log2(control depth)')
    #boxplot(log2(cnr.depth) ~ gene, data=z, horizontal=T, las=2, cex.axis=0.5,
    #        main=plasma, xlab='log2(plasma depth)')
    #boxplot(cnr.log2 ~ gene, data=z, horizontal=T, las=2, cex.axis=0.5,
    #        main=plasma, ylim=c(-2,2), xlab='log2r.cnvkit')
    #abline(v=0, col='red', lty='dotted')
    #dev.off()

    # plot log2r QC in png format
    #png(paste0(out.qc.dir, '/', plasma, '.correl.png'))
    #plot(log2(z$ref.depth), z$cnr.log2, col=colrs[z$gene], main=ttl,
    #     xlab='log2(control depth)', ylab='log2ratio(plasma vs. control)')
    #abline(h=0, col='red', lty='dotted')
    #text(x=x1,y=1.5, paste0('Off-target\ncor = ', sprintf('%0.3g', r1)))
    #text(x=x2,y=1.5, paste0('On-target\ncor = ', sprintf('%0.3g', r2)))
    #dev.off()

    # plot all chromosomes
    #pchr = (ggplot(size)
    #    + geom_vline(xintercept=c(size$start,size$end[nrow(size)]), color='darkgray')
    #    + scale_x_continuous(breaks=size$pos, labels=gsub('chr', '', size$chr))
        #+ scale_y_continuous(breaks=seq(-10,10,1))
    #    + theme_bw(base_size=8) + ylab('log2 (plasma/normal)') + xlab(NULL)
    #    + geom_hline(yintercept=0, color='red', size=0.5)
    #    + geom_vline(xintercept=c(ar.start2, ar.stop2), color='blue', size=0.2)
    #    + geom_vline(xintercept=c(e.start2, e.stop2), color='red', size=0.2)
    #    + geom_vline(data=g2, aes(xintercept=pos), color='orange',
    #                 size=0.1, linetype='dotted')
    #    + geom_text(data=g2, aes(x=pos, y = -4, label=gene), hjust=0,
    #                angle=90, color='darkred', size=1.5)
    #    + theme(legend.position='none')
    #)

    #z = cbind(sample=plasma, z, stringsAsFactors=F)
    #if (is.null(cnres)){cnres = z}else{cnres = rbind(cnres, z)}

    #z = prepChr(z)
    #cns = prepChr(cns)
    #p = (pchr + geom_point(data=z, aes(x=pos, y = cnr.log2, color=targeted),
    #                        shape=1,size=0.1, alpha=0.75)
    #       + geom_segment(data=cns, aes(x=start, xend=end, y=log2, yend=log2),
    #                      color='black')
    #       + ggtitle(ttl)# + theme(legend.position='top')
    #       + scale_y_continuous(limits=c(-4,4))
           #+ scale_color_manual(values=c('No'='gray', 'Yes'='purple', 'U'='cyan',
           #                              '?'='blue'))
    #)
    #pX = (p + coord_cartesian(xlim=c(size$start[size$chr == 'chrX'],
    #                        size$end[size$chr == 'chrX'])) + ggtitle(NULL))

    # AR locus
    #pAR.l2r = (p + coord_cartesian(xlim=c(xmin, xmax))
    #           + theme(legend.position='none')
    #           )

    #pAR.depth = (pchr + geom_point(data=z, aes(x=pos, y=ref.depth),
    #                               color='darkgray', size=0.2)
    #           + geom_point(data=z, aes(x=pos, y=cnr.depth), color='red', size=0.2)
    #           + coord_cartesian(xlim=c(xmin, xmax))
    #           + scale_y_log10()
    #           + ylab('depth')
    #      )


    #p.l2r.ctl = (pchr + geom_point(data=z, aes(x=pos, y=cnr.log2r.ctl, color=targeted),
    #                               size=0.1, shape=1, alpha=0.75)
    #           + geom_segment(data=cs, aes(x=start, xend=end, y=seg.mean, yend=seg.mean),
    #                      color='black')
    #           + scale_y_continuous(limits=c(-4,4))
    #           + ylab('log2r.ctlnorm')
    #)
    #pX.ctrl = (p.l2r.ctl + coord_cartesian(xlim=c(size$start[size$chr == 'chrX'],
    #                            size$end[size$chr == 'chrX'])) + ggtitle(NULL))
#
    #pAR.l2r.ctl = p.l2r.ctl + coord_cartesian(xlim=c(xmin, xmax)) + ggtitle(ttl)



    #pAR.weight = (pchr + geom_point(data=z, aes(x=pos, y=cnr.weight),
    #                               color='darkgray', size=0.2)
    #           + geom_point(data=z, aes(x=pos, y=gc),
    #                               color='black', size=0.2)
    #           + coord_cartesian(xlim=c(xmin, xmax))
    #           + scale_y_continuous(limit=c(0,1))
    #           + xlab('bin weight')
    #)




    #D = 10^6
    ## TMRPSS2 locus
    #pt2 = p + coord_cartesian(xlim=c(g2$start[g2$gene=='TMPRSS2'] - D,
    #                                g2$end[g2$gene=='TMPRSS2'] + D))
    #pmyc = p + coord_cartesian(xlim=c(g2$start[g2$gene=='MYC'] - D,
    #                                g2$end[g2$gene=='MYC'] + D))
    #ppten = p + coord_cartesian(xlim=c(g2$start[g2$gene=='PTEN'] - D,
    #                                g2$end[g2$gene=='PTEN'] + D))
 
  
    #png(paste0(out.plot.dir, '/', plasma, '.log2r-allchrs.png'), unit='in', res=150,
    #    width=8, height=11)
    #grid.arrange(p,p.l2r.ctl, pX.ctrl, pX,pAR.l2r,pAR.l2r.ctl, #pAR.depth, #ncol=1)
    #             pt2,pmyc,ppten,ncol=1)
    #grid.arrange(p,pX,pAR.l2r,ncol=1)
    #dev.off()
    #print(pAR.l2r.ctl)

#}
#dev.off()

#system('for f in `ls plots/QC/png/*png`; do convert $f $f.pdf; done; pdftk plots/QC/png/correl_*pdf output plots/QC/log2r-vs-control-depth.all-patients.pdf;  pdftk plots/QC/png/allchrs_*pdf output plots/QC/log2r.allchrs.pdf')



