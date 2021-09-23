################################################################################
# Read in the Sanger SNPs for the DO founders and count the total number of 
# SNPs, number of SNPs in transcripts, number of SNPs in CDS.
# Note that the snps_by_type will have duplicate counts because some variants
# have multiple consequences. The total_snps should be unique SNPs & indels.
# snps_by_gene should have accurate SNP & indel counts by gene.
#
# Daniel Gatti
# dan.gatti@jax.org
# 2021-09-21
################################################################################

### LIBRARIES ###
library(AnnotationHub)
library(BiocParallel)
library(tidyverse)
library(ensembldb)
library(VariantAnnotation)

### VARIABLES ###

results_dir    = '/projects/compsci/dgatti/projects/snp_count'

sanger_dir = '/projects/omics_share/temp_reorg/Mouse/Sanger/REL_2004'
sanger_file = file.path(sanger_dir, 'mgp_REL2005_snps_indels.vcf.gz')

# 10 Mb chunk size for SNP processing.
chunk_size = 10e6

# BiocParallel/foreach setup
num_cores = 5
BPPARAM = MulticoreParam(workers = num_cores)
register(BPPARAM, default = TRUE)

### FUNCTIONS ###

# Filter the snps or indels to retain high quality, polymorphic variants.
# variants: VaraintAnnotation::CollapsedVCF containing either SNPs or Indels.
filter_variants = function(variants) {
  
    if(nrow(variants) == 0) {
      return(variants)
    }
  
    # Filter to retain high quality SNPs & indels with PASS in fixed$FILTER column.
    variants = subset(variants, fixed(variants)$FILTER == 'PASS')
    
    if(nrow(variants) == 0) {
      return(variants)
    }

    # Filter to retain SNPs & indels that don't have a missing genotype. ('./.' in GT) 
    variants = subset(variants, rowSums(geno(variants)$GT == './.') == 0)
    
    if(nrow(variants) == 0) {
      return(variants)
    }

    # Filter to retain SNPs & indels with variation in the CC/DO founders.
    variants = subset(variants, rowSums(geno(variants)$GT == '0/0') < ncol(variants))

    return(variants)

} # filter_variants()


# Intersect the variants with the exons of genes.
# variants: VaraintAnnotation::CollapsedVCF containing either SNPs or Indels.
# gene_table: data.frame containing a 'gene_id' column with the ENSMUG ID for
#             all genes in the current Ensembl build.
# ensembl: EnsemblDB containing gene information from AnnotationHub.
# gr: GRanges with current chromosome and range.
# type: string containing the type of exons to intersect. Default = 'exon'. Exon
#       uses all exons. Cds uses only coding sequences.
intersect_var_with_genes = function(variants, gene_table, ensembl, gr, 
                                    type = c('exon', 'cds')) {

  type = match.arg(type)

  # Get exons in current range.
  exons = NULL
  if(type == 'exon') {
    exons = exonsBy(ensembl, by = 'gene', filter = GRangesFilter(gr))
  } else {
    exons = cdsBy(ensembl,   by = 'gene', filter = GRangesFilter(gr))
  }
  exons = unlist(exons)
  exons$gene_id = names(exons)
  # Overlap exons with variants, retaining any kind of overlap.
  ol = countOverlaps(exons, variants, type = 'any')
  # Keep exons with at least one overlap.
  ol = ol[ol > 0]
  ol = data.frame(gene_id = names(ol), num_variants = ol)
  
  # Aggregate overlap counts by gene.
  if(nrow(ol) > 0) {
    var_ct = aggregate(ol$num_variants, by = list(gene_id = ol$gene_id), FUN = sum,
                        simplify = TRUE)
    m = match(var_ct$gene_id, gene_table$gene_id)
    gene_table$num_variants[m] = gene_table$num_variants[m] + var_ct$x
  }

  return(gene_table)

} # intersect_var_with_genes()


# Count the number of each type of variants. Note that one variant may
# have mulitple consequences, so we can't just sum this table to get
# the total number of variants.
# variants: VaraintAnnotation::CollapsedVCF containing either SNPs or Indels.
# var_table: data.frame containing a 'gene_id' column with the ENSMUG ID for
#             all genes in the current Ensembl build. This will also contain
#             columns for variant types. 
tally_vars_by_type = function(variants, var_table) {
    # Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|VARIANT_CLASS|SYMBOL_SOURCE|HGNC_ID|SIFT
    csq = info(variants)$CSQ
    csq = lapply(csq, strsplit, split = '\\|')
    csq = lapply(csq, sapply, '[', c(2, 5))
    csq = mapply(function(nm, cs) { colnames(cs) = rep(nm, ncol(cs)); cs }, 
                 rownames(variants), csq)
    csq = sapply(csq, t)
    csq = sapply(csq, unique)
    csq = do.call(rbind, csq)
    colnames(csq) = c('consequence', 'gene_id')

    # Using tidyverse...
    # I expect warnings about separate() filling in empty columns.
    csq = data.frame(variant_id = rownames(csq), csq) %>%
            separate(consequence, into = LETTERS[1:10], sep = '&') %>%
            pivot_longer(cols = LETTERS[1:10], names_to = 'LET', 
                         values_to = 'consequence') %>%
            drop_na() %>%
            dplyr::select(-LET) %>%
            dplyr::count(gene_id, consequence)
    
    stopifnot(nchar(csq$gene_id) <= 18)
    
    # Gather results by gene. Not that some SNPs may be counted more than
    # once because they have multiple consequences.    
    csq = csq %>%
            pivot_wider(names_from = consequence, values_from = n)

    var_table = full_join(var_table, csq)
    
    return(var_table)

} # tally_vars_by_type()


### MAIN ###

# Get the Ensembl version from the Sanger VCF header.
# Using Rsamtools::headerTabix() because I can't find the same info
# in VariantAnnotation::scanVcfHeader().
header = Rsamtools::headerTabix(sanger_file)$header
header = header[grep('ensembl', header)]
header = strsplit(header, ' ')[[1]]
header = header[grep('^ensembl=', header)]
ensembl_ver = as.numeric(gsub('^ensembl=|\\.[0-9]+db[0-9]+$', '', header))

###########################################################
### NOTE: Ensembl 97 GTF is missing from AnnotationHub. ###
###########################################################

# Download the correct Ensembl version GTF.
#hub = AnnotationHub()
#ensembl_title = paste0('Mus_musculus.GRCm38.', ensembl_ver, '.gtf')
#hub_id = names(hub)[hub$title == ensembl_title]

# Get the unique chromosomes in the Sanger VCF.
header = scanVcfHeader(sanger_file)

unique_chr = seqnames(seqinfo(header))
chrlength  = seqlengths(seqinfo(header))

# Get the strains.
strains = samples(header)[c(4, 2, 34, 37, 19, 40, 51)]

main_fxn = function(chr) {

  # Get the ensembldb from AnnotationHub instead.
  # I had to move this into the main_fxn because the different workers
  # couldn't use the sam DB connection.
  hub = AnnotationHub()
  ensembl_title = paste0('Ensembl ', ensembl_ver, ' EnsDb for Mus musculus')
  hub_id = names(hub)[hub$title == ensembl_title]
  ensembl = hub[[hub_id]]

  # Set up results.
  # Get genes on this chromosome.
  current_ensid = genes(ensembl, filter = SeqNameFilter(chr))$gene_id

  snps_by_gene = data.frame(gene_id = current_ensid, num_variants = 0)
  snps_by_cds  = data.frame(gene_id = current_ensid, num_variants = 0)
  snps_by_type = data.frame(gene_id = current_ensid)

  indels_by_gene = data.frame(gene_id = current_ensid, num_variants = 0)
  indels_by_cds  = data.frame(gene_id = current_ensid, num_variants = 0)
  indels_by_type = data.frame(gene_id = current_ensid)

  # Set up chunks (chromosome range to query).
  chunk = 0
  start = chunk * chunk_size
  end   = start + chunk_size
  
  print(paste('CHR', chr))

  while(end < chrlength[chr]) {

    print(paste('   ', chunk))

    # TBD: Get more information and filter more stringently?
    param = ScanVcfParam(fixed = names(fixed(header)), 
                         info  = c('INDEL', 'CSQ'), 
                         geno  = 'GT',
                         samples = strains, 
                         trimEmpty = TRUE,
                         which = GRanges(seqnames = chr, IRanges(start = start, end = end)))
    # Get the SNPs & indels in this region.
    snps = readVcf(file = sanger_file, param = param)
    
    # Split out indels.
    indels = subset(snps, info(snps)$INDEL)
    snps   = subset(snps, !info(snps)$INDEL)

    # Filter to keep high quality, polymorphic variants.
    snps   = filter_variants(snps)
    indels = filter_variants(indels)

    # Perform a raw intersection of the SNPs & indels with exons & CDSs.
    current_gr = GRanges(chr, IRanges(start, end))
    snps_by_gene   = intersect_var_with_genes(snps,   snps_by_gene,   ensembl, current_gr,
                                              type = 'exon')
    indels_by_gene = intersect_var_with_genes(indels, indels_by_gene, ensembl, current_gr,
                                              type = 'exon')
    snps_by_cds    = intersect_var_with_genes(snps,   snps_by_cds,    ensembl, current_gr,
                                            type = 'cds')
    indels_by_cds  = intersect_var_with_genes(indels, indels_by_cds,  ensembl, current_gr,
                                            type = 'cds')

    # Use the consequence column to count number of each type of variant.
    snps_by_type   = tally_vars_by_type(snps,   snps_by_type) 
    indels_by_type = tally_vars_by_type(indels, indels_by_type)

    # Increment chunking.
    chunk = chunk + 1
    start = chunk * chunk_size
    end   = start + chunk_size

  } # while()

  # Add gene locations to SNP & Indels summaries.
  genes = genes(ensembl, filter = SeqNameFilter(chr))
  genes = as.data.frame(genes)
  genes = genes[,c('gene_id', 'seqnames', 'start', 'end', 'strand', 
                   'gene_name', 'gene_biotype')]

  snps_by_gene = full_join(genes, snps_by_gene)
  snps_by_cds  = full_join(genes, snps_by_cds)

  indels_by_gene = full_join(genes, indels_by_gene)
  indels_by_cds  = full_join(genes, indels_by_cds)

  saveRDS(snps_by_gene, file = file.path(results_dir, 
                        paste0('snps_by_gene_chr', chr, '.rds')))
  saveRDS(snps_by_cds,  file = file.path(results_dir, 
                        paste0('snps_by_cds_chr', chr, '.rds')))
  saveRDS(snps_by_type, file = file.path(results_dir, 
                        paste0('snps_by_type_chr', chr, '.rds')))

  saveRDS(indels_by_gene, file = file.path(results_dir, 
                          paste0('indels_by_gene_chr', chr, '.rds')))
  saveRDS(indels_by_cds,  file = file.path(results_dir, 
                          paste0('indels_by_cds_chr', chr, '.rds')))
  saveRDS(indels_by_type, file = file.path(results_dir, 
                          paste0('indels_by_type_chr', chr, '.rds')))

} # main_fxn(chr)

# Run the main code.
x = bplapply(unique_chr, FUN = main_fxn)


# Combine the results files.
df = data.frame(pattern = c('snps_by_gene_chr', 'indels_by_gene_chr', 
                            'snps_by_type_chr', 'indels_by_type_chr',
                            'snps_by_cds_chr',  'indels_by_cds_chr'),
                outfile = c('snps_by_gene.rds', 'indels_by_gene.rds',
                            'snps_by_type.rds', 'indels_by_type.rds',
                            'snps_by_cds.rds',  'indels_by_cds.rds'))
for(i in 1:nrow(df)) {
  files = dir(results_dir, pattern = df$pattern[i], full.names = TRUE)
  sbg = lapply(files, readRDS)
  if(i <= 2 | i >= 5) {
    # We can rbind the 'by_gene' & 'by_cds' files because they have the same columns. 
    sbg = do.call(rbind, sbg)
  } else {
    # We have to join the 'by_type' files becuase they have different columns.
    tmp = sbg[[1]]
    for(j in 2:length(sbg)) {
      tmp = full_join(tmp, sbg[[j]])
    } # for(j)
    sbg = tmp
  } # else
  saveRDS(sbg, file = file.path(results_dir, df$outfile[i]))
} # for(i)

