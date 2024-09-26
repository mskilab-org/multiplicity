#' @name snvplicity
#' @title Converts counts to copies.
#'
#' @param somatic_snv Path to somatic SNV file
#' @param germline_snv Path to germline SNV file
#' @param jabba_rds Path to jabba file
#' @param tumor_name Expected name of tumor as annotated in the VCF
#' @param normal_name Expected name of normal as annotated in the VCF
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @return Returns a GRanges with counts and converted copies
#' @export 
snvplicity = function(somatic_snv,
                      germline_snv,
                      jabba_rds,
                      snpeff_path,
                      tumor_name = NULL,
                      normal_name = NULL){

  gg = gG(jabba = jabba_rds)

  somatic.snv = parsesnpeff(vcf = somatic_snv,
                            snpeff_path = snpeff_path,
                            coding_alt_only = FALSE,
                            filterpass = FALSE,
                            tumor_id = tumor_name,
                            normal_id = normal_name,
                            keepfile = FALSE,
                            altpipe = TRUE)

  germline.snv = parsesnpeff(vcf = germline_snv,
                             snpeff_path = snpeff_path,
                             coding_alt_only = FALSE,
                             filterpass = FALSE,
                             tumor_id = tumor_name,
                             normal_id = normal_name,
                             keepfile = FALSE,
                             altpipe = TRUE)

  message("SOMATIC SNV PROCESSING")
  message("ALL VARIANTS FOR SOMATIC SNV REGARDLESS OF FILTER == PASS STATUS")

  names(somatic.snv) <- NULL
  somatic.snv.filtered = gr2dt(somatic.snv)[, .SD[1], by = c("seqnames", "start", "end")]
                                        #browser()
  somatic.snv.filtered <- somatic.snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr

  #normalization
  unique.somatic.snv = somatic.snv.filtered[,c('ref', 'alt')] %>%
    gr.nochr %Q%
    (seqnames %in% c(1:22, "X")) %>%
    unique
  unique.somatic.snv$major.count <- pmax(unique.somatic.snv$ref, unique.somatic.snv$alt)
  unique.somatic.snv$minor.count <- pmin(unique.somatic.snv$ref, unique.somatic.snv$alt)
  somatic.m = length(unique.somatic.snv)
  somatic.sf = sum(unique.somatic.snv$major.count + unique.somatic.snv$minor.count) / somatic.m
  unique.somatic.snv$major.count <- unique.somatic.snv$major.count / somatic.sf
  unique.somatic.snv$minor.count <- unique.somatic.snv$minor.count / somatic.sf

  message("GERMLINE SNV PROCESSING")
  message("ONLY FILTER == PASS VARIANTS FOR GERMLINE SNVS")

  names(germline.snv) <- NULL
  germline.snv.filtered = gr2dt(germline.snv %Q% (FILTER == "PASS"))[, .SD[1], by = c("seqnames", "start", "end")]
  germline.snv.filtered <- germline.snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr
    
  #normalization
  unique.germline.snv = germline.snv.filtered[,c('ref', 'alt')] %>%
    gr.nochr %Q%
    (seqnames %in% c(1:22, "X")) %>%
    unique
  unique.germline.snv$major.count <- pmax(unique.germline.snv$ref, unique.germline.snv$alt)
  unique.germline.snv$minor.count <- pmin(unique.germline.snv$ref, unique.germline.snv$alt)
  germline.m = length(unique.germline.snv)
  germline.sf = sum(unique.germline.snv$major.count + unique.germline.snv$minor.count) / germline.m
  unique.germline.snv$major.count <- unique.germline.snv$major.count / germline.sf
  unique.germline.snv$minor.count <- unique.germline.snv$minor.count / germline.sf

    
  #constitutional_cn assignment
  #c_subj == 1 for major allele
  #c_subj == 1 for autosomes and X chromosome in females, 0 for X and Y in males
  ncn.x = gg$nodes$dt[(seqnames == "X" | seqnames == "chrX"),
                      weighted.mean(cn,
                      w = end - start + 1,
                      na.rm = TRUE)]
    
  message("mean CN of X: ", ncn.x)
  somatic.ncn.vec = rep(2, length(unique.somatic.snv))
  germline.ncn.vec = rep(2, length(unique.germline.snv))
  if (ncn.x < 1.4) { # if male
    
    message("Adjusting ncn for XY")
    somatic.ncn.vec[which(as.character(seqnames(unique.somatic.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
    germline.ncn.vec[which(as.character(seqnames(unique.germline.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
  }
    #values(unique.somatic.snv)[, "ncn"] = ncn.vec
  values(unique.somatic.snv)[, "minor_constitutional_cn"] = somatic.ncn.vec - 1
  values(unique.germline.snv)[, "minor_constitutional_cn"] = germline.ncn.vec - 1

  jab = gg2jab(gg)
  ss.p = jab$segstats[ as.logical( strand(jab$segstats)=='+' ) ]

  unique.somatic.snv = gr.val(unique.somatic.snv, ss.p, "cn", na.rm = T)
  unique.germline.snv = gr.val(unique.germline.snv, ss.p, "cn", na.rm = T)

  ###### SOMATIC APPLICATION OF FORMULA #####
    
  tau_hat = mean(unique.somatic.snv$cn)
  alpha = jab$purity
  beta = alpha / (alpha * tau_hat + 2*(1 - alpha))
  gamma = 2*(1 - alpha) / (alpha * tau_hat + 2*(1 - alpha))

  mcols(unique.somatic.snv)$major_snv_copies =
    (2 * unique.somatic.snv$major.count - gamma) / (2 * beta)

  mcols(unique.somatic.snv)$minor_snv_copies =
    (2 * unique.somatic.snv$major.count - gamma * unique.somatic.snv$minor_constitutional_cn)/ 2 * beta

  somatic.variants = gr.val(somatic.snv.filtered %>% gr.nochr,
                            unique.somatic.snv[,c('major.count',
                                                  'minor.count',
                                                  'major_snv_copies',
                                                  'minor_snv_copies')],
                            c('major.count', 'minor.count',
                            'major_snv_copies', 'minor_snv_copies')) %>% #%Q%
                                        #(FILTER == "PASS") %>%
    gr2dt

  somatic.variants[alt >= ref, total_copies := major_snv_copies]
  somatic.variants[alt < ref, total_copies := minor_snv_copies]
  somatic.variants[, VAF := alt / (alt + ref)]
  somatic.variants = somatic.variants %>% dt2gr


  ###### GERMLINE APPLICATION OF FORMULA #####
    
  tau_hat = mean(unique.germline.snv$cn)
  alpha = jab$purity
  beta = alpha / (alpha * tau_hat + 2*(1 - alpha))
  gamma = 2*(1 - alpha) / (alpha * tau_hat + 2*(1 - alpha))

  mcols(unique.germline.snv)$major_snv_copies =
  (2 * unique.germline.snv$major.count - gamma) / (2 * beta)

  mcols(unique.germline.snv)$minor_snv_copies =
  (2 * unique.germline.snv$major.count - gamma * unique.germline.snv$minor_constitutional_cn)/ 2 * beta

  germline.variants = gr.val(germline.snv.filtered %>% gr.nochr,
                             unique.germline.snv[,c('major.count',
                                                    'minor.count',
                                                    'major_snv_copies',
                                                    'minor_snv_copies')],
                              c('major.count', 'minor.count',
                                'major_snv_copies', 'minor_snv_copies')) %>% #%Q%
      #(FILTER == "PASS") %>%
      gr2dt

  germline.variants[alt >= ref, total_copies := major_snv_copies]
  germline.variants[alt < ref, total_copies := minor_snv_copies]
  germline.variants[, VAF := alt / (alt + ref)]
  germline.variants = germline.variants %>% dt2gr

  return(list(somatic.variants, germline.variants))
}

#' @name parsesnpeff
#' @title how is this not a paragraph?
#'
#' @param vcf path to snpeff vcf
#' @param pad Exposed argument to skitools::ra.overlaps()
#' @param tumor_id Tumor name as annotated in the VCF
#' @param normal_id Normal name as annotated in the VCF
#' @param filterpass lorem impsum
#' @param coding_alt_only lorem impsum
#' @param geno lorem impsum
#' @param gr lorem impsum
#' @param kepefile lorem impsum
#' @param altpipe lorem impsum
#' @param debug lorem impsum
#' @return GRangesList of breakpoint pairs with junctions that overlap removed
#' @export
parsesnpeff = function (vcf,
                        snpeff_path,
                        tumor_id = NULL,
                        normal_id = NULL,
                        filterpass = TRUE,
                        coding_alt_only = TRUE, 
                        geno = NULL,
                        gr = NULL,
                        keepfile = FALSE,
                        altpipe = FALSE, 
                        debug = FALSE) 
{
    if (debug) 
        browser()
    out.name = paste0("tmp_", rand.string(), ".vcf.gz")
    tmp.path = paste0(tempdir(), "/", out.name)
    if (!keepfile) 
        on.exit(unlink(tmp.path))
    try2({
        catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
        onepline = paste0(snpeff_path, "/scripts/vcfEffOnePerLine.pl")
        if (coding_alt_only) {
          filt = paste0("java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar ",
                        snpeff_path, "/SnpSift.jar ",
                        "filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\"")
            if (filterpass)
                cmd = sprintf(paste(catcmd, "%s | %s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), 
                              vcf, onepline, filt, tmp.path)
            else cmd = sprintf("cat %s | %s | %s | bcftools norm -Ov -m-any | bgzip -c > %s", 
                               vcf, onepline, filt, tmp.path)
        }
        else {
            filt = ""
            if (filterpass) 
                cmd = sprintf(paste(catcmd, "%s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), 
                              vcf, onepline, tmp.path)
            else cmd = sprintf(paste(catcmd, "%s | %s | bcftools norm -Ov -m-any | bgzip -c > %s"), 
                               vcf, onepline, tmp.path)
        }
        system(cmd)
    })
    if (!altpipe) 
        out = grok_vcf(tmp.path, long = TRUE, geno = geno, gr = gr)
    else {
        vcf = readVcf(tmp.path)
        vcf = S4Vectors::expand(vcf)
        rr = within(rowRanges(vcf), {
            REF = as.character(REF)
            ALT = as.character(ALT)
        })
        ann = as.data.table(tstrsplit(unlist(info(vcf)$ANN), 
                                      "\\|"))[, 1:15, with = FALSE, drop = FALSE]
        fn = c("allele", "annotation", "impact", "gene", "gene_id", 
               "feature_type", "feature_id", "transcript_type", 
               "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos", 
               "protein_pos", "distance")
        data.table::setnames(ann, fn)
        #browser()
        if ("AD" %in% names(geno(vcf))) {
          vcf.ncol <- ncol(geno(vcf)$AD)
          vcf.names <- colnames(geno(vcf)$AD)
          ## grab the last item in the array if ids not specified... presumably tumor
          if(vcf.ncol > 1){
            vcf.order <- if(!is.na(tumor_id) && !is.na(normal_id) &&
                            !is.null(tumor_id) && !is.null(normal_id)){
                           c(base::match(c(tumor_id, normal_id), vcf.names))
                         } else {vcf.ncol} # just grab the last column!
            adep = setnames(as.data.table(geno(vcf)$AD[, vcf.order[1], 1:2]), 
                            c("ref", "alt"))
            try2(expr = {
              adep.n = setnames(as.data.table(geno(vcf)$AD[, vcf.order[2], 1:2]), 
                                c("normal.ref", "normal.alt"))
              adep = adep %>% cbind(adep.n)
            })
            gt = geno(vcf)$GT
          }
          else {
            adep = setnames(as.data.table(geno(vcf)$AD[, vcf.ncol, 1:2]), 
                            c("ref", "alt"))
            gt = geno(vcf)$GT
          }
        }
        else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in% 
                     c(names(geno(vcf))))) {
            this.col = dim(geno(vcf)[["AU"]])[2]
            d.a = geno(vcf)[["AU"]][, , 1, drop = F][, this.col, 
                                                     1]
            d.g = geno(vcf)[["GU"]][, , 1, drop = F][, this.col, 
                                                     1]
            d.t = geno(vcf)[["TU"]][, , 1, drop = F][, this.col, 
                                                     1]
            d.c = geno(vcf)[["CU"]][, , 1, drop = F][, this.col, 
                                                     1]
            mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
            rm("d.a", "d.g", "d.t", "d.c")
            refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
            refid = ifelse(!isSNV(vcf), NA_integer_, refid)
            altid = match(as.character(VariantAnnotation::fixed(vcf)$ALT), colnames(mat))
            altid = ifelse(!isSNV(vcf), NA_integer_, altid)
            refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
            altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
            this.icol = dim(geno(vcf)[["TAR"]])[2]
            refindel = d.tar = geno(vcf)[["TAR"]][, , 1, drop = F][, 
                                                                   this.icol, 1]
            altindel = d.tir = geno(vcf)[["TIR"]][, , 1, drop = F][, 
                                                                   this.icol, 1]
            adep = data.table(ref = coalesce(refsnv, refindel), 
                              alt = coalesce(altsnv, altindel))
            gt = NULL
        }
        else {
            message("ref and alt count columns not recognized")
            adep = NULL
            gt = NULL
        }
        ## mcols(rr) = BiocGenerics::cbind(mcols(rr), ann, adep, 
        ##                                 gt = gt[, 1])
        mcols(rr) = cbind(as.data.table(mcols(rr)),
                          ann,
                          adep, 
                          gt = gt[, 1])
        out = rr
    }
    this.env = environment()
    return(this.env$out)
}

#' @name rand.string
#' @title make a random string
#'
#' @return random string
#' @author Someone from Stackoverflow
#' @export rand.string
rand.string = function(n=1, length=12)
{
  randomString <- c(1:n)                  # initialize vector
  for (i in 1:n)
  {
    randomString[i] <- paste(sample(c(0:9, letters, LETTERS),
                                    length, replace=TRUE),
                             collapse="")
  }
  return(randomString)
}

#' @name try2
#' @title wrapper around tryCatch - robust to parallel:: functions; A slightly more robust version of try that works within the parallel:: set of functions that pre-deploy a cluster.
#'
#' @export
try2 = function(expr, ..., finally) {
  tryCatch(expr,
           error = function(e) {
             msg = structure(paste(conditionMessage(e), conditionCall(e), sep = "\n"), class = "err")
             cat("Error: ", msg, "\n\n")
             return(msg)
           },
           finally = finally,
           ... = ...)
}
