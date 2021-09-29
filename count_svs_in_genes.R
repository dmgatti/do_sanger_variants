################################################################################
# Count structural variants (insertions & deletoins) that intersect with
# genes. We don't have detailed information about SVs that overlap or occur
# within each other. But this will give us a rough count.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-09-29
################################################################################

### LIBRARIES ###

library(AnnotationHub)
library(ensembldb)
library(GenomicRanges)
library(Rsamtools)

### VARIABLES ###

base_dir = '/projects/compsci/dgatti/projects/do_sanger_variants'

sanger_dir = '/projects/omics_share/temp_reorg/Mouse/Sanger/REL_1606_SV'
del_file = file.path(sanger_dir, 'mgpv5.SV_deletions.bed.gz')
ins_file = file.path(sanger_dir, 'mgpv5.SV_insertions.bed.gz')

del_exon_file = file.path(base_dir, 'deletions_in_exons')
del_cds_file  = file.path(base_dir, 'deletions_in_cds')
ins_exon_file = file.path(base_dir, 'insertions_in_exons')
ins_cds_file  = file.path(base_dir, 'insertions_in_cds')

# The README says that they aligned to GRCm38 and they didn't use an Ensembl
# build to intersect with genes. I'm going to use Ensembl 97 because that's
# what I used for the SNPs & indels.
ensembl_ver = 97
hub = AnnotationHub()
ensembl_title = paste0('Ensembl ', ensembl_ver, ' EnsDb for Mus musculus')
hub_id = names(hub)[hub$title == ensembl_title]
ensembl = hub[[hub_id]]

### FUNCTIONS ###

parse_bed = function(data, col_names) {

  data = strsplit(data, '\t')
  data = matrix(unlist(data), nrow = length(data), byrow = T,
                dimnames = list(NULL, col_names))
  data = data[,c('chr', 'start', 'end', 'type', do_founders)]
  data = data.frame(data, stringsAsFactors = FALSE)
  data$start = as.numeric(data$start)
  data$end   = as.numeric(data$end)
  
  return(data)

} # parse_bed()


filter_svs = function(data) {

  tmp = data[,5:ncol(data)]
  
  return(subset(data, rowMeans(tmp == '.;.;.;.') < 1.0))

} # filter_svs()


### MAIN ###

###########
# Deletions
tabix = TabixFile(del_file)
hdr = headerTabix(tabix)
col_names = strsplit(hdr$header, '\t')[[1]]
col_names = sub('^#', '', col_names)

# Get columns with DO founders.
do_founders = col_names[c(10, 7, 31, 33, 21, 35, 40)]

# Get the chromosmes in the file.
all_chr = hdr$seqnames

# For each chromosome, get the SVs and intersect with exons and CDSs.
for(chr in all_chr) {

  print(paste('CHR', chr))

  current_gr = GRanges(seqnames = chr, IRanges(start = 0, end = 200e6))
  svs = scanTabix(tabix, param = current_gr)[[1]]
  svs = parse_bed(data = svs, col_names = col_names)
  svs = filter_svs(svs)
  
  svs = GRanges(seqnames = svs$chr, IRanges(start = svs$start, 
                end = svs$end))

  # Get genes on current chromosome.
  genes = as.data.frame(genes(ensembl, filter = SeqNameFilter(chr)))
  genes = genes[,c('seqnames', 'start', 'end', 'strand', 'gene_id',
                   'gene_name', 'gene_biotype')]
  colnames(genes)[1] = 'chr'

  vars = data.frame(col_name    = c('num_exons',   'num_cds'),
                    file_prefix = c(del_exon_file, del_cds_file))

  # Intersect exons of CDS with SVs.
  for(j in 1:2) {
  
    exons = NULL
    if(j == 1) {
       exons = exonsBy(ensembl, by = 'gene', filter = SeqNameFilter(chr))
    } else {
       exons = cdsBy(ensembl,   by = 'gene', filter = SeqNameFilter(chr))
    }
    exons = unlist(exons)
    exons$gene_id = names(exons)
    ol = subsetByOverlaps(exons, svs)
    gene_ct = table(ol$gene_id)
    gene_ct = data.frame(gene_ct)
    colnames(gene_ct) = c('gene_id', vars$col_name[j])
    gene_ct = merge(genes, gene_ct, by = 'gene_id', all.x = TRUE)
    gene_ct[,vars$col_name[j]][is.na(gene_ct[,vars$col_name[j]])] = 0
    write.table(gene_ct, file = paste0(vars$file_prefix[j], '.csv'), sep = ',', 
                append = chr != all_chr[1], col.names = chr == all_chr[1],
                quote = FALSE, row.names = FALSE) 
  } # for(j)
} # for(chr)

# Convert files from CSV to RDS.
exons = read.csv(paste0(del_exon_file, '.csv'))
saveRDS(exons, file = paste0(del_exon_file, '.rds'))
file.remove(paste0(del_exon_file, '.csv'))

cds = read.csv(paste0(del_cds_file, '.csv'))
saveRDS(cds, file = paste0(del_cds_file, '.rds'))
file.remove(paste0(del_cds_file, '.csv'))

###########
# Insertions
tabix = TabixFile(ins_file)
hdr = headerTabix(tabix)
col_names = strsplit(hdr$header, '\t')[[1]]
col_names = sub('^#', '', col_names)

# Get columns with DO founders.
do_founders = col_names[c(10, 7, 31, 33, 21, 35, 40)]

# Get the chromosmes in the file.
all_chr = hdr$seqnames

# For each chromosome, get the SVs and intersect with exons and CDSs.
for(chr in all_chr) {

  print(paste('CHR', chr))

  current_gr = GRanges(seqnames = chr, IRanges(start = 0, end = 200e6))
  svs = scanTabix(tabix, param = current_gr)[[1]]
  svs = parse_bed(data = svs, col_names = col_names)
  svs = filter_svs(svs)
  
  svs = GRanges(seqnames = svs$chr, IRanges(start = svs$start, 
                end = svs$end))

  # Get genes on current chromosome.
  genes = as.data.frame(genes(ensembl, filter = SeqNameFilter(chr)))
  genes = genes[,c('seqnames', 'start', 'end', 'strand', 'gene_id',
                   'gene_name', 'gene_biotype')]
  colnames(genes)[1] = 'chr'

  vars = data.frame(col_name    = c('num_exons',   'num_cds'),
                    file_prefix = c(ins_exon_file, ins_cds_file))

  # Intersect exons of CDS with SVs.
  for(j in 1:2) {
  
    exons = NULL
    if(j == 1) {
       exons = exonsBy(ensembl, by = 'gene', filter = SeqNameFilter(chr))
    } else {
       exons = cdsBy(ensembl,   by = 'gene', filter = SeqNameFilter(chr))
    }
    exons = unlist(exons)
    exons$gene_id = names(exons)
    ol = subsetByOverlaps(exons, svs)
    gene_ct = table(ol$gene_id)
    gene_ct = data.frame(gene_ct)
    colnames(gene_ct) = c('gene_id', vars$col_name[j])
    gene_ct = merge(genes, gene_ct, by = 'gene_id', all.x = TRUE)
    gene_ct[,vars$col_name[j]][is.na(gene_ct[,vars$col_name[j]])] = 0
    write.table(gene_ct, file = paste0(vars$file_prefix[j], '.csv'), sep = ',', 
                append = chr != all_chr[1], col.names = chr == all_chr[1],
                quote = FALSE, row.names = FALSE) 
  } # for(j)
} # for(chr)

# Convert files from CSV to RDS.
exons = read.csv(paste0(ins_exon_file, '.csv'))
saveRDS(exons, file = paste0(ins_exon_file, '.rds'))
file.remove(paste0(ins_exon_file, '.csv'))

cds = read.csv(paste0(ins_cds_file, '.csv'))
saveRDS(cds, file = paste0(ins_cds_file, '.rds'))
file.remove(paste0(ins_cds_file, '.csv'))








