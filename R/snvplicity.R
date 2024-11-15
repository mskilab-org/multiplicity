#' @name snvplicity
#' @title Converts counts to copies.
#'
#' @param somatic_snv Path to somatic SNV file
#' @param germline_snv Path to germline SNV file
#' @param het_pileups_wgs Path to het_pileups_wgs file
#' @param jabba_rds Path to jabba file
#' @param tumor_name Expected name of tumor as annotated in the VCF; highly suggested this is provided for accurate VCF parsing
#' @param normal_name Expected name of normal as annotated in the VCF; highly suggested this is provided for accurate VCF parsing
#' @param purity Purity of inputted jabba_rds (optional if metadata of gGraph contains purity)
#' @param ploidy Ploidy of inputted jabba_rds (optional if metadata of gGraph contains ploidy)
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @return Returns a GRanges with counts and converted copies
#' @export 
snvplicity = function(somatic_snv = NULL,
                      germline_snv = NULL,
                      het_pileups_wgs,
                      jabba_rds,
                      snpeff_path,
                      tumor_name = NULL,
                      normal_name = NULL,
                      purity = NULL, 
                      ploidy = NULL,
                      verbose = FALSE){

  if(verbose){message("reading in jabba file")}

  if(is.null(somatic_snv) & is.null(germline_snv))
    stop("Somatic VCF and/or Germline VCF must be provided.")
  
  gg = gG(jabba = jabba_rds)
  hets = gr.val(fread(het_pileups_wgs) %>% dt2gr,
                gg$nodes$gr, 'cn') %>%
    gr.nochr() %>% gr2dt()
  hets[seqnames %in% c(1:22, "X")] %>% dt2gr -> hets
  
  #constitutional_cn assignment
  #c_subj == 1 for major allele
  #c_subj == 1 for autosomes and X chromosome in females, 0 for X and Y in males
  ncn.x = gg$nodes$dt[(seqnames == "X" |
                       seqnames == "chrX" |
                       seqnames == "23" |
                       seqnames == "chr23"), #' MSK-FACETS calls chrX as chr23
                      weighted.mean(cn,
                                    w = end - start + 1,
                                    na.rm = TRUE)]
  message("mean CN of X: ", ncn.x)
  if(ncn.x < 1.4) message("sex determination: XY")
  else message("sex determination:: XX")

  jab = gg$nodes$gr %>% gr.nochr()
  jab = dt2gr(gr2dt(jab)[seqnames == 23, seqnames := "X"])
  ss.p = jab[as.logical(strand(jab) == "+")]
  
  somatic.variants = NULL
  germline.variants = NULL

  if(verbose){message("succesfully read in JaBbA graph!")}

  if(!is.na(somatic_snv) && !is.null(somatic_snv)){

    somatic.snv = parsesnpeff(vcf = somatic_snv,
                              snpeff_path = snpeff_path,
                              coding_alt_only = FALSE,
                              filterpass = FALSE,
                              tumor_id = tumor_name,
                              normal_id = normal_name,
                              keepfile = FALSE,
                              altpipe = TRUE,
                              verbose = verbose)

    names(somatic.snv) <- NULL
    somatic.snv.filtered = gr2dt(somatic.snv)[, .SD[1], by = c("seqnames", "start", "end")]
    somatic.snv.filtered <- somatic.snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr
    
    #'normalization
     unique.somatic.snv = somatic.snv.filtered[,c('ref', 'alt')] %>%
      gr.nochr %Q%
      (!seqnames == c("Y")) %>%
      unique
    unique.somatic.snv$major.count <- pmax(unique.somatic.snv$ref, unique.somatic.snv$alt, na.rm = T)
    unique.somatic.snv$minor.count <- pmin(unique.somatic.snv$ref, unique.somatic.snv$alt, na.rm = T)
    somatic.m = ((unique.somatic.snv$ref + unique.somatic.snv$alt) %>% length) * 2
    somatic.sf = sum(unique.somatic.snv$major.count + unique.somatic.snv$minor.count, na.rm = T) / somatic.m
    unique.somatic.snv$major.count <- unique.somatic.snv$major.count / somatic.sf
    unique.somatic.snv$minor.count <- unique.somatic.snv$minor.count / somatic.sf

    somatic.ncn.vec = rep(2, length(unique.somatic.snv))

    if (ncn.x < 1.4) { # if male 
      if(verbose) message("Adjusting ncn for XY")
      somatic.ncn.vec[which(as.character(seqnames(unique.somatic.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
    }

    values(unique.somatic.snv)[, "minor_constitutional_cn"] = somatic.ncn.vec - 1
    unique.somatic.snv = gr.val(unique.somatic.snv, ss.p, "cn", na.rm = T)

    ###### SOMATIC APPLICATION OF FORMULA #####
    
    tau_hat = mean(hets$cn)
    tau = gg$meta$ploidy
    alpha = if(!is.null(gg$meta$purity)) {
              gg$meta$purity
            } else if (!is.null(purity)){
              purity
            } else {
              unique(jab$purity)
            }
    beta = alpha / (alpha * tau + 2*(1 - alpha))
    gamma = 2*(1 - alpha) / (alpha * tau_hat + 2*(1 - alpha))

    if(verbose) message(paste0("average total CN of somatic loci: " , tau_hat))
    if(verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
    if(verbose) message("applying transformation")
    
    mcols(unique.somatic.snv)$major_snv_copies =
                              (2 * unique.somatic.snv$major.count - gamma) / (2 * beta)

    mcols(unique.somatic.snv)$minor_snv_copies =
                              (2 * unique.somatic.snv$minor.count - (gamma * unique.somatic.snv$minor_constitutional_cn))/ (2 * beta)

    somatic.variants = gr.val(somatic.snv.filtered %>% gr.nochr,
                              unique.somatic.snv[,c('major.count',
                                                    'minor.count',
                                                    'major_snv_copies',
                                                    'minor_snv_copies')],
                              c('major.count', 'minor.count',
                                'major_snv_copies', 'minor_snv_copies')) %>% gr2dt()

    somatic.variants[alt >= ref, total_copies := major_snv_copies]
    somatic.variants[alt < ref, total_copies := minor_snv_copies]
    somatic.variants[total_copies <= 0, total_copies := 0]
    somatic.variants[, VAF := alt / (alt + ref)]
    somatic.variants = somatic.variants %>% dt2gr
  }

  if(!is.na(germline_snv) && !is.null(germline_snv)){

    germline.snv = parsesnpeff(vcf = germline_snv,
                              snpeff_path = snpeff_path,
                              coding_alt_only = FALSE,
                              filterpass = FALSE,
                              tumor_id = tumor_name,
                              normal_id = normal_name,
                              keepfile = FALSE,
                              altpipe = TRUE,
                              verbose = verbose)

    names(germline.snv) <- NULL
    germline.snv.filtered = gr2dt(germline.snv %Q% (FILTER == "PASS"))[, .SD[1], by = c("seqnames", "start", "end")]
    germline.snv.filtered <- germline.snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr

    #'normalization
    unique.germline.snv = germline.snv.filtered[,c('ref', 'alt')] %>%
      gr.nochr %Q%
      (!seqnames == c("Y")) %>%
      unique
    unique.germline.snv$major.count <- pmax(unique.germline.snv$ref, unique.germline.snv$alt)
    unique.germline.snv$minor.count <- pmin(unique.germline.snv$ref, unique.germline.snv$alt)
    germline.m = length(unique.germline.snv)
    germline.sf = sum(unique.germline.snv$major.count + unique.germline.snv$minor.count) / germline.m
    unique.germline.snv$major.count <- unique.germline.snv$major.count / germline.sf
    unique.germline.snv$minor.count <- unique.germline.snv$minor.count / germline.sf

    germline.ncn.vec = rep(2, length(unique.germline.snv))

    if (ncn.x < 1.4) { # if male 
      message("Adjusting ncn for XY")
      germline.ncn.vec[which(as.character(seqnames(unique.germline.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
    }

    values(unique.germline.snv)[, "minor_constitutional_cn"] = germline.ncn.vec - 1
    unique.germline.snv = gr.val(unique.germline.snv, ss.p, "cn", na.rm = T)
    
		###### GERMLINE APPLICATION OF FORMULA #####
    tau_hat = mean(hets$cn)
    alpha = gg$meta$purity
    beta = alpha / (alpha * tau_hat + 2*(1 - alpha))
    gamma = 2*(1 - alpha) / (alpha * tau_hat + 2*(1 - alpha))

    if(verbose) message(paste0("average total CN of germline loci: " , tau_hat))
    if(verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
    if(verbose) message("applying transformation")

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
  }
    
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
parsesnpeff = function (
  vcf,
  snpeff_path,
  tumor_id = NULL,
  normal_id = NULL,
  filterpass = TRUE,
  coding_alt_only = TRUE, 
  geno = NULL,
  gr = NULL,
  keepfile = FALSE,
  altpipe = FALSE, 
  debug = FALSE,
  verbose = FALSE,
  bcftools = "/gpfs/data/imielinskilab/Software/miniforge3/envs/vcftools/bin/bcftools" # FIXME: hardcoded for now.
  ) {
  if (debug)
    browser()
  tmp.path = tempfile(pattern = "tmp_", fileext = ".vcf.gz")
  if (!keepfile)
    on.exit(unlink(tmp.path))
  try2({
    catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
    onepline = paste0(snpeff_path, "/scripts/vcfEffOnePerLine.pl")
    if (verbose)(message(paste0("applying SnpSift to VCF: ", vcf)))
    if (coding_alt_only) {
      if(verbose)(message("Coding alterations only."))
      filt = paste0("java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar ",
        snpeff_path, "/SnpSift.jar ",
        "filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\"")
      if (filterpass) {
        if(verbose)(message("Coding alterations only and FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, filt, bcftools, tmp.path)  } else {if(verbose)(message("Coding alterations only."))
          cmd = sprintf("cat %s | %s | %s | bgzip -c > %s", vcf, onepline, filt, tmp.path)
        }
    } else {
      filt = ""
      if (filterpass){
        if(verbose)(message("FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, bcftools, tmp.path)
      } else {
        if(verbose)(message("All alterations included."))
        cmd = sprintf(paste(catcmd, "%s | %s | bgzip -c > %s"), vcf, onepline, tmp.path)
      }
    }
    if(verbose)(message("Performing command."))
    system(cmd)
    if(verbose)(message("SnpSift successfully applied!"))
  })
  if (!altpipe)
    out = grok_vcf(tmp.path, long = TRUE, geno = geno, gr = gr)
  else {
    if (verbose)(message(paste0("reading in SnpSift VCF.")))
    vcf = readVcf(tmp.path)
    rr = rowRanges(vcf)
    rr$REF = as.character(rr$REF)
    ann = as.data.table(tstrsplit(unlist(info(vcf)$ANN),
      "\\|"))[, 1:15, with = FALSE, drop = FALSE]
    fn = c("allele", "annotation", "impact", "gene", "gene_id",
      "feature_type", "feature_id", "transcript_type", 
      "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos", 
      "protein_pos", "distance")
    data.table::setnames(ann, fn)
    if ("AD" %in% names(geno(vcf))) {
      if (verbose)(message(paste0("parsing AD field in VCF.")))
      vcf.ncol <- ncol(geno(vcf)$AD)
      if (verbose)(message(paste0(vcf.ncol, " columns found in the VCF.")))
      vcf.names <- colnames(geno(vcf)$AD)
      if(vcf.ncol > 1) {
        if (verbose)(message(paste0("parsing tumor and normal alt/ref counts.")))
        ## grab the last item in the array if ids not specified... presumably tumor
        tumor_col_id = vcf.ncol ## assume last column bydefault
        normal_col_id = NA_integer_
        if (!is.null(tumor_id) && !is.na(tumor_id)) {
          tumor_col_id = base::match(tumor_id, vcf.names)
        }
        if (!is.null(normal_id) && !is.na(normal_id)) {
          normal_col_id = base::match(normal_id, vcf.names)
        }
        adep = data.table::transpose(geno(vcf)$AD[,tumor_col_id])
        adep = as.data.table(adep)
        adep = base::subset(
          adep,
          select = 1:2
        )
        data.table::setnames(adep, c("ref", "alt"))
        adep$normal.ref = NA_integer_
        adep$normal.alt = NA_integer_
        if (!is.na(normal_col_id)) {
          adep.n = data.table::transpose(geno(vcf)$AD[,normal_col_id])
          adep$normal.ref = adep.n[[1]]
          adep$normal.alt = adep.n[[2]]
        }
        gt = geno(vcf)$GT
      } else {
        if (verbose)(message(paste0("parsing only tumor alt/ref counts.")))
        adep = setnames(as.data.table(geno(vcf)$AD[, vcf.ncol, 1:2]),
          c("ref", "alt"))
        gt = geno(vcf)$GT
      }
    } else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in%
                     c(names(geno(vcf))))) {
      if(verbose)(message("parsing AU/GU/CU/TAR/TIR fields in VCF."))
      this.col = dim(geno(vcf)[["AU"]])[2]
      d.a = geno(vcf)[["AU"]][, , 1, drop = F][, this.col, 1]
      d.g = geno(vcf)[["GU"]][, , 1, drop = F][, this.col, 1]
      d.t = geno(vcf)[["TU"]][, , 1, drop = F][, this.col, 1]
      d.c = geno(vcf)[["CU"]][, , 1, drop = F][, this.col, 1]
      mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
      #rm("d.a", "d.g", "d.t", "d.c")
      refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
      refid = ifelse(!isSNV(vcf), NA_integer_, refid)
      altid = match(as.character(VariantAnnotation::fixed(vcf)$ALT), colnames(mat))
      altid = ifelse(!isSNV(vcf), NA_integer_, altid)
      refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
      altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
      this.icol = dim(geno(vcf)[["TAR"]])[2]
      refindel = d.tar = geno(vcf)[["TAR"]][, , 1, drop = F][, this.icol, 1]
      altindel = d.tir = geno(vcf)[["TIR"]][, , 1, drop = F][, this.icol, 1]
      try2({
        n.d.a = geno(vcf)[["AU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.g = geno(vcf)[["GU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.t = geno(vcf)[["TU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.c = geno(vcf)[["CU"]][, , 1, drop = F][, this.col - 1, 1]
        n.mat = cbind(A = n.d.a, G = n.d.g, T = n.d.t, C = n.d.c)
        #rm("n.d.a", "n.d.g", "n.d.t", "n.d.c")
        n.refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(n.mat))
        n.refid = ifelse(!isSNV(vcf), NA_integer_, refid)
        n.altid = match(as.character(VariantAnnotation::fixed(vcf)$ALT), colnames(n.mat))
        n.altid = ifelse(!isSNV(vcf), NA_integer_, altid)
        n.refsnv = n.mat[cbind(seq_len(nrow(n.mat)), refid)]
        n.altsnv = n.mat[cbind(seq_len(nrow(n.mat)), altid)]
        n.refindel = n.d.tar = geno(vcf)[["TAR"]][, , 1, drop = F][, this.icol - 1, 1]
        n.altindel = n.d.tir = geno(vcf)[["TIR"]][, , 1, drop = F][, this.icol - 1, 1]
      })
      adep = data.table(ref = coalesce(refsnv, refindel),
        alt = coalesce(altsnv, altindel))
      try2({
        adep.n = data.table(normal.ref = coalesce(n.refsnv, n.refindel),
          normal.alt = coalesce(n.altsnv, n.altindel))
        adep =  adep %>% cbind(adep.n)
      })
      gt = NULL
    } else {
      message("ref and alt count columns not recognized")

      adep = NULL
      gt = NULL
    }
    ## cbinding S4Vector DataFrame objects works much faster with dataframes
    ## Need to figure out problems with data.table conversions as some point
    mcols(rr) = cbind(mcols(rr),
      data.table::setDF(ann),
      data.table::setDF(adep),
      gt = gt[, 1])
    rr = S4Vectors::expand(rr, "ALT")
    rr$ALT = as.character(rr$ALT)
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
  randomString <- c(1:n) #initialize vector
  for (i in 1:n){
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
