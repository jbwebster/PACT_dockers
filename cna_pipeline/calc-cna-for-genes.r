library(reshape2)

args = commandArgs(T)

gene.overlaps = args[1]
samples.info = args[2]
target.genes = args[3]
nsd = as.numeric(args[4]) # Number of sd above/below the mean required to make a CNA call


x = read.table(gene.overlaps,
        header=T, stringsAsFactors=F)
x$sample = sub('(^.*)/.*', '\\1', x$seg)

x$sum.log2r = x$log2r*x$overlap

x1 = aggregate(sum.log2r ~ gene + sample, data=x, FUN='sum')
x2 = aggregate(overlap ~ gene + sample, data=x, FUN='sum')
x = merge(x1,x2)
x$log2r = x$sum.log2r/x$overlap

z = x

m = read.table(samples.info,header=T)
m = m[m$sample %in% z$sample,]

g = read.table(target.genes, header=T,
    stringsAsFactors=F, sep='\t')
g = g[, c('gene', 'description')]


z1 = merge(g, z) #Merge target gene list and gene_overlap results

calc.log2r.mean.sd <- function(z1){
    m = aggregate(log2r ~ gene, data=z1, FUN='mean')
    colnames(m) = c('gene', 'log2r.mean')
    s = aggregate(log2r ~ gene, data=z1, FUN='sd')
    colnames(s) = c('gene', 'log2r.sd')
    v = merge(m,s)
    v = merge(g, v)
    v = v[order(v$log2r.sd, decreasing=T),]
    return(v)
}

v = calc.log2r.mean.sd(z1)
v = v[order(v$log2r.sd, decreasing=T),]
rownames(v) = NULL
z1 = z1[order(z1$description, -match(z1$gene, v$gene)),]


# calc log2 stat
z1 = merge(z1, m)
log2r.stat = NULL
for (type in unique(z1$type)){
    sel = which(z1$type == type)
    if (length(sel) == 0){next}
    l2r.stat = calc.log2r.mean.sd(z1[sel,])
    l2r.stat$type = type
    log2r.stat = rbind(log2r.stat, l2r.stat)
}


gx = log2r.stat[log2r.stat$type == 'plasma',]
gx = gx[order(gx$log2r.mean),]
log2r.stat$geneid = factor(log2r.stat$gene, levels=gx$gene)

log2r.stat$pt.log2r.mean = log2r.stat$log2r.mean
log2r.stat$pt.log2r.sd = log2r.stat$log2r.sd
ctrl.log2r.mean = mean(log2r.stat$pt.log2r.mean[log2r.stat$description == 'CN-control'], na.rm=T)
ctrl.log2r.sd = mean(log2r.stat$pt.log2r.sd[log2r.stat$description == 'CN-control'], na.rm=T)

#nsd based on input parameter
gain = ctrl.log2r.mean + nsd*ctrl.log2r.sd
loss = ctrl.log2r.mean - nsd*ctrl.log2r.sd

z1$cna.call = 'neutral'
z1$cna.call[z1$log2r > gain] = 'gain'
z1$cna.call[z1$log2r < loss] = 'loss'
write.table(z1, file='cna-call-targeted-genes.tsv', row.names=F,
    sep='\t', quote=F)
