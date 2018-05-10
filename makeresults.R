library(tidyverse)

herv_fam <- read.table('results/herv_families.txt', sep='\t', header=T, stringsAsFactors = F)
herv <- read.table('results/HERV_rmsk.hg38.v2.tsv', 
                   sep='\t', header=T, stringsAsFactors = F) %>%
    mutate(
        length=end - start + 1
        ) %>%
    separate(locus, c("family"), sep='_', extra='drop', remove=F) %>%
    left_join(herv_fam, by='family') %>%
    mutate(
        family=factor(family, levels=herv_fam$family),
        group=factor(group, levels=unique(herv_fam$group))
    ) %>%
    select(locus, chrom, start, end, length, family, group, category)


samples <- read.table('samples.txt', stringsAsFactors=F)[,1]
# samples <- samples[samples != 'SB72']

lapply(samples, function(s){
    lines <- readLines(file.path(s, 'bt2_multi.summary.txt'))
    c(s,
      as.numeric(gsub('^(\\d+) .*', '\\1',  lines[1])),
      as.numeric(gsub('^(\\d+\\.\\d+)% .*', '\\1',  lines[length(lines)]))
    )
}) %>% do.call(rbind, .) %>% data.frame(stringsAsFactors=F) -> metx
names(metx) <- c('sample', 'total_reads', 'alnrate')
metx$total_reads <- as.integer(metx$total_reads)
metx$alnrate <- as.double(metx$alnrate) * 1e-2
metx <- metx %>% mutate(naln=floor(total_reads*alnrate))
row.names(metx) <- metx$sample
metx$sample <- NULL 

write.table(metx, 'results/metrics.tsv', sep='\t', quote=F)


lapply(samples, function(s){
    tmp <- read.table(file.path(s, 'inform-telescope_report.tsv'),
                       sep='\t', header=T, stringsAsFactors=F)
    ret <- data.frame(transcript=herv$locus, stringsAsFactors=F) %>%
        left_join(tmp, by='transcript') %>%
        mutate(
            gene_id = transcript,
            count = final_count
        ) %>%
        select(gene_id, count)
    ret[is.na(ret)] <- 0
    stopifnot(all(ret$gene_id == herv$locus))
    ret$gene_id <- NULL
    names(ret) <- c(s)
    ret
}) %>%
    do.call(cbind, .) %>% 
    data.frame -> herv.counts
row.names(herv.counts) <- herv$locus

write.table(herv.counts, 'results/telescope_counts.tsv', sep='\t', quote=F)

