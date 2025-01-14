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
                      het_pileups_wgs = NULL,
                      tumor_cbs = NULL,
                      tumor_dryclean = NULL,
                      dryclean_field, 
                      jabba_rds,
                      mask = NULL,
                      inferred_sex = NA,
                      read_size = 151,
                      snpeff_path,
                      tumor_name = NULL,
                      normal_name = NULL,
                      normfactor = 2,
                      filterpass = FALSE,
                      tau_in_gamma = FALSE,
                      purity = NULL, 
                      ploidy = NULL,
                      verbose = FALSE){

  if(verbose){message("reading in jabba file")}

  if(is.null(somatic_snv) & is.null(germline_snv) & is.null(het_pileups_wgs))
    stop("Somatic VCF, Germline VCF, and/or HetSNPs file must be provided.")
  
  gg = gG(jabba = jabba_rds)

  read_size = as.numeric(read_size)

  if(!(is.na(tumor_dryclean) || is.null(tumor_dryclean))){
    dryclean.cov <- readRDS(tumor_dryclean)
    cov.vector = mcols(dryclean.cov)[names(mcols(dryclean.cov)) %in% dryclean_field][, 1]
    mcols(dryclean.cov) <- NULL
    mcols(dryclean.cov)$bincov <- cov.vector
    mcols(dryclean.cov)$avg_basecov <- dryclean.cov$bincov * 2 * read_size / width(dryclean.cov)
  } else {
    dryclean.cov <- NULL
  }

  #' CBS WILL OVERRIDE DRYCLEAN INPUTS
  #' CBS reduces variance relative to dryclean rescaling alone
  if(!(is.na(tumor_cbs) || is.null(tumor_cbs))){
    cbs.cov <- readRDS(tumor_cbs)
    cbs.vector = mcols(cbs.cov)[names(mcols(cbs.cov)) %in% "seg.mean"][, 1]
    mcols(cbs.cov) <- NULL
    mcols(cbs.cov)$bincov <- exp(cbs.vector)
    cbs.cov <- gr.tile(seqlengths(cbs.cov), 1000) %>% gr.val(cbs.cov, "bincov")
    cbs.cov[which(is.na(cbs.cov$bincov))]$bincov <- 0
    mcols(cbs.cov)$avg_basecov <- cbs.cov$bincov * 2 * read_size / width(cbs.cov)
    dryclean.cov <- cbs.cov
  } else {
    cbs.cov <- NULL
  }

  
  if(!(is.na(mask) || is.null(mask))){
    mask = readRDS(mask)     
  }
  
  if(is.na(inferred_sex)){

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
    if(ncn.x < 1.4){
      inferred_sex <- "M"
      message("sex determination: XY")
    } else {
      inferred_sex <- "F"
      message("sex determination:: XX")
    }
  } else if(grepl("^m", inferred_sex, ignore.case = T)) {
    message("sex provided: XY")
    inferred_sex <- "M"
  } else if(grepl("^f", inferred_sex, ignore.case = T)){
    message("sex determination:: XX")
    inferred_sex <- "F"
  } else {
    stop("sex not provided and could not be inferred from provided inputs.")
  }
  

  ploidy = if(!is.null(ploidy)){
             ploidy
           } else {
             gg$meta$ploidy
           }
  purity = if(!is.null(gg$meta$purity)) {
            gg$meta$purity
          } else if (!is.null(purity)){
            purity
          } else {
            unique(jab$purity)
          }

  jab = gg$nodes$gr %>% gr.nochr()
  jab = dt2gr(gr2dt(jab)[seqnames == 23, seqnames := "X"])
  cn.gr = jab[as.logical(strand(jab) == "+")]
  
  somatic.variants = NULL
  germline.variants = NULL
  het.pileups = NULL

  if(verbose){message("succesfully read in JaBbA graph!")}

  if(!is.na(somatic_snv) && !is.null(somatic_snv)){
    somatic.variants <- transform_snv(
      vcf = somatic_snv,
      cn_gr = cn.gr,
      snpeff_path = snpeff_path,
      dryclean.cov = dryclean.cov,
      basecov_field = "avg_basecov",
      inferred_sex = inferred_sex,
      tumor_id = tumor_name,
      normal_id = normal_name,
      purity = purity,
      ploidy = ploidy,
      filterpass = filterpass,
      verbose = verbose,
      normfactor = normfactor,
      tau_in_gamma = tau_in_gamma,
      mask = mask
    )    
  }

  if(!is.na(germline_snv) && !is.null(germline_snv)){

    germline.variants <- transform_snv(
      vcf = germline_snv,
      cn_gr = cn.gr,
      snpeff_path = snpeff_path,
      dryclean.cov = dryclean.cov,
      basecov_field = "avg_basecov",
      inferred_sex = inferred_sex,
      tumor_id = tumor_name,
      normal_id = normal_name,
      purity = purity,
      ploidy = ploidy,
      filterpass = TRUE,
      verbose = verbose,
      normfactor = normfactor,
      tau_in_gamma = tau_in_gamma,
      mask = mask
    )
  }

  
  if(!is.na(het_pileups_wgs) && !is.null(het_pileups_wgs)){
    
    het.pileups <- transform_hets(
      hets = het_pileups_wgs,
      cn_gr = cn.gr,
      dryclean.cov = dryclean.cov,
      basecov_field = "avg_basecov",
      inferred_sex = inferred_sex,
      purity = purity,
      ploidy = ploidy,
      verbose = verbose,
      normfactor = normfactor,
      tau_in_gamma = tau_in_gamma,
      mask = mask
    )    
  }
    
  return(list(somatic.variants, germline.variants, het.pileups))
}


#' @name transform_snv
#' @title function to inhale SNV-laden VCF and transform counts to to estimated altered copies using linear transformation via beta, gamma conversions
#'
#' @param vcf Path to snpeff vcf
#' @param dryclean.cov Tumor drycleaned coverage for denoising/rescaling of REF/ALT counts (OPTIONAL; if provided will denoise)
#' @param basecov_field field name of dryclean.cov that indicates average Base Coverage per bin
#' @param tumor_id Tumor name as annotated in the VCF (highly suggested this is provided)
#' @param normal_id Normal name as annotated in the VCF (highly suggested this is provided)
#' @param purity Inferred purity input for conversion
#' @param ploidy Inferred ploidy input for conversion
#' @param inferred_sex enum of ["M", "F"] that indicates inferred of true sex of patient sample
#' @param filterpass process only FILTER == PASS variants?
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @return GRanges of transformed read counts
#' @export
transform_snv = function(vcf,
                         cn_gr,
                         dryclean.cov = NULL,
                         basecov_field = "avg_basecov",
                         tumor_id = NA,
                         normal_id = NA,
                         ploidy,
                         purity,
                         inferred_sex,
                         filterpass = FALSE,
                         snpeff_path,
                         verbose = TRUE,
                         normfactor = 2,
                         tau_in_gamma = FALSE,
                         mask = mask){

  snv = parsesnpeff(vcf = vcf,
                    snpeff_path = snpeff_path,
                    coding_alt_only = FALSE,
                    filterpass = filterpass,
                    tumor_id = tumor_id,
                    normal_id = normal_id,
                    keepfile = FALSE,
                    altpipe = TRUE,
                    verbose = verbose)

  #browser()

  names(snv) <- NULL; fields.to.carry <- character()
  snv.filtered = gr2dt(snv)[, .SD[1], by = c("seqnames", "start", "end")]
  snv.filtered <- snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr

  #browser()

  if(!is.null(mask)){
    snv.filtered = gr.val(snv.filtered, mask, "blacklisted", na.rm = T) %Q%
      (!is.na(blacklisted)) %Q%
      (blacklisted == F)
  }
    
  #'normalization
  unique.snv = snv.filtered[,c('ref', 'alt')] %>%
    gr.nochr %Q%
    (!seqnames == c("Y")) %>%
    unique

  if(!(is.na(dryclean.cov) || is.null(dryclean.cov))){
    unique.snv = gr.val(unique.snv, dryclean.cov, basecov_field, na.rm = T)
    unique.snv$ref_denoised <- unique.snv$ref * unique.snv$avg_basecov /
      (unique.snv$ref + unique.snv$alt)
    unique.snv$alt_denoised <- unique.snv$alt * unique.snv$avg_basecov /
      (unique.snv$ref + unique.snv$alt)

    unique.snv$major.count <- pmax(unique.snv$ref_denoised, unique.snv$alt_denoised, na.rm = T)
    unique.snv$minor.count <- pmin(unique.snv$ref_denoised, unique.snv$alt_denoised, na.rm = T)
    fields.to.carry <- c("ref_denoised", "alt_denoised")
  } else {
    unique.snv$major.count <- pmax(unique.snv$ref, unique.snv$alt, na.rm = T)
    unique.snv$minor.count <- pmin(unique.snv$ref, unique.snv$alt, na.rm = T)
  }
    
  m = length(unique.snv$ref + unique.snv$alt) * normfactor
  sf = sum(unique.snv$major.count + unique.snv$minor.count, na.rm = T) / m
  unique.snv$major.count <- unique.snv$major.count / sf
  unique.snv$minor.count <- unique.snv$minor.count / sf

  ncn.vec = rep(2, length(unique.snv))

  if (inferred_sex == "M") { # if male 
    if(verbose) message("Adjusting ncn for XY")
    ncn.vec[which(as.character(seqnames(unique.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
  }

  values(unique.snv)[, "minor_constitutional_cn"] = ncn.vec - 1
  unique.snv = gr.val(unique.snv, cn_gr, "cn", na.rm = T)


  ###### APPLICATION OF FORMULA #####
    
  tau_hat = mean(unique.snv$cn, na.rm = T)
  tau = if(tau_in_gamma){
          ploidy
        } else {
          tau_hat
        }
  alpha = purity
  beta = alpha / (alpha * ploidy + 2*(1 - alpha))
  gamma = 2*(1 - alpha) / (alpha * tau + 2*(1 - alpha))

  if(verbose) message(paste0("average total CN of somatic loci: " , tau_hat))
  if(verbose) message(paste0("ploidy of tumor sample: " , ploidy))
  if(verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
  if(verbose) message("applying transformation")
    
  mcols(unique.snv)$major_snv_copies =
                    (2 * unique.snv$major.count - (gamma * 2)) / (2 * beta)

  ## mcols(unique.snv)$minor_snv_copies =
  ## (2 * unique.snv$minor.count - (gamma * unique.snv$minor_constitutional_cn))/ (2 * beta)

  mcols(unique.snv)$minor_snv_copies =
                    (2 * unique.snv$minor.count - (gamma * 0)) / (2 * beta)

  mcols(unique.snv)$total_snv_copies =
                    unique.snv$major_snv_copies + unique.snv$minor_snv_copies
    
  variants = gr.val(snv.filtered %>% gr.nochr,
                    unique.snv[,c(fields.to.carry,
                                  'major.count',
                                  'minor.count',
                                  'total_snv_copies',
                                  'major_snv_copies',
                                  'minor_snv_copies',
                                  'cn')],
                    c(fields.to.carry,
                      'major.count', 'minor.count',
                      'total_snv_copies',
                      'major_snv_copies',
                      'minor_snv_copies',
                      'cn')) %>% gr2dt()

  variants[alt >= ref, altered_copies := major_snv_copies]
  variants[alt < ref , altered_copies := minor_snv_copies]
  #variants[total_copies <= 0, total_copies := 0]
  variants[, VAF := alt / (alt + ref)]
  variants = variants %>% dt2gr

  return(variants)
}


#' @name transform_hets
#' @title function to inhale SNV-laden VCF and transform counts to to estimated altered copies using linear transformation via beta, gamma conversions
#'
#' @param vcf Path to snpeff vcf
#' @param dryclean.cov Tumor drycleaned coverage for denoising/rescaling of REF/ALT counts (OPTIONAL; if provided will denoise)
#' @param basecov_field field name of dryclean.cov that indicates average Base Coverage per bin
#' @param tumor_id Tumor name as annotated in the VCF (highly suggested this is provided)
#' @param normal_id Normal name as annotated in the VCF (highly suggested this is provided)
#' @param purity Inferred purity input for conversion
#' @param ploidy Inferred ploidy input for conversion
#' @param inferred_sex enum of ["M", "F"] that indicates inferred of true sex of patient sample
#' @param filterpass process only FILTER == PASS variants?
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @return GRanges of transformed read counts
#' @export
transform_hets = function(hets,
                         cn_gr,
                         dryclean.cov = NULL,
                         basecov_field = "avg_basecov",
                         ploidy,
                         purity,
                         inferred_sex,
                         verbose = TRUE,
                         tau_in_gamma = FALSE,
                         mask = mask,
                         normfactor = 2){

  hets.gr <- fread(hets) %>% dt2gr()

  if(!is.null(mask)){
    hets.gr = gr.val(hets.gr, mask, "blacklisted", na.rm = T) %Q%
      (!is.na(blacklisted)) %Q%
      (blacklisted == F)
  }
      
  #'normalization
  unique.snv = hets.gr[,c('ref.count.t', 'alt.count.t')] %>%
    gr.nochr %Q%
    (!seqnames == c("Y")) %>%
    unique

  fields.to.carry  <- c("ref", "alt")

  names(mcols(unique.snv)) <- fields.to.carry
  
  if(!(is.na(dryclean.cov) || is.null(dryclean.cov))){
    unique.snv = gr.val(unique.snv, dryclean.cov, basecov_field, na.rm = T)
    unique.snv$ref_denoised <- unique.snv$ref * unique.snv$avg_basecov /
      (unique.snv$ref + unique.snv$alt)
    unique.snv$alt_denoised <- unique.snv$alt * unique.snv$avg_basecov /
      (unique.snv$ref + unique.snv$alt)

    unique.snv$major.count <- pmax(unique.snv$ref_denoised, unique.snv$alt_denoised, na.rm = T)
    unique.snv$minor.count <- pmin(unique.snv$ref_denoised, unique.snv$alt_denoised, na.rm = T)
    fields.to.carry <- c(fields.to.carry, "ref_denoised", "alt_denoised")
  } else {
    unique.snv$major.count <- pmax(unique.snv$ref, unique.snv$alt, na.rm = T)
    unique.snv$minor.count <- pmin(unique.snv$ref, unique.snv$alt, na.rm = T)
  }
    
  m = length(unique.snv$ref + unique.snv$alt) * normfactor
  sf = sum(unique.snv$major.count + unique.snv$minor.count, na.rm = T) / m
  unique.snv$major.count <- unique.snv$major.count / sf
  unique.snv$minor.count <- unique.snv$minor.count / sf

  ncn.vec = rep(2, length(unique.snv))

  if (inferred_sex == "M") { # if male 
    if(verbose) message("Adjusting ncn for XY")
    ncn.vec[which(as.character(seqnames(unique.snv)) %in% c("chrX", "chrY", "X", "Y"))] = 1
  }

  values(unique.snv)[, "minor_constitutional_cn"] = ncn.vec - 1
  unique.snv = gr.val(unique.snv, cn_gr, "cn", na.rm = T)

  ###### APPLICATION OF FORMULA #####
    
  tau_hat = mean(unique.snv$cn, na.rm = T)
  tau = if(tau_in_gamma){
          ploidy
        } else {
          tau_hat
        }
  alpha = purity
  beta = alpha / (alpha * ploidy + 2*(1 - alpha))
  gamma = 2*(1 - alpha) / (alpha * tau + 2*(1 - alpha))

  if(verbose) message(paste0("average total CN of somatic loci: " , tau_hat))
  if(verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
  if(verbose) message("applying transformation")
    
  mcols(unique.snv)$major_snv_copies =
                    (2 * unique.snv$major.count - (gamma * 1)) / (2 * beta)

  mcols(unique.snv)$minor_snv_copies =
                    (2 * unique.snv$minor.count - (gamma * unique.snv$minor_constitutional_cn))/ (2 * beta)

  ## mcols(unique.snv)$minor_snv_copies =
  ##                   (2 * unique.snv$minor.count - (gamma * 1)) / (2 * beta)

  mcols(unique.snv)$total_snv_copies =
                    unique.snv$major_snv_copies + unique.snv$minor_snv_copies
    
  variants = gr.val(hets.gr %>% gr.nochr,
                    unique.snv[,c(fields.to.carry,
                                  'major.count',
                                  'minor.count',
                                  'total_snv_copies',
                                  'major_snv_copies',
                                  'minor_snv_copies',
                                  'cn')],
                    c(fields.to.carry,
                      'major.count', 'minor.count',
                      'total_snv_copies',
                      'major_snv_copies',
                      'minor_snv_copies',
                      'cn')) %>% gr2dt()

  variants[alt >= ref, altered_copies := major_snv_copies]
  variants[alt < ref , altered_copies := minor_snv_copies]
  #variants[total_copies <= 0, total_copies := 0]
  variants[, VAF := alt / (alt + ref)]
  variants = variants %>% dt2gr

  return(variants)
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
                        debug = FALSE,
                        verbose = FALSE) {
  if (debug)
    browser()
  out.name = paste0("tmp_", rand.string(), ".vcf.gz")
  tmp.path = paste0(tempdir(), "/", out.name)
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
        cmd = sprintf(paste(catcmd, "%s | %s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, filt, tmp.path)  } else {if(verbose)(message("Coding alterations only."))
            cmd = sprintf("cat %s | %s | %s | bcftools norm -Ov -m-any | bgzip -c > %s", vcf, onepline, filt, tmp.path)
          }
    } else {
      filt = ""
      if (filterpass){
        if(verbose)(message("FILTER == PASS variants only."))
        	cmd = sprintf(paste(catcmd, "%s | %s | bcftools view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, tmp.path)
      } else {
        if(verbose)(message("All alterations included."))
        cmd = sprintf(paste(catcmd, "%s | %s | bcftools norm -Ov -m-any | bgzip -c > %s"), vcf, onepline, tmp.path)
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
    vcf = S4Vectors::expand(vcf)
    rr = S4Vectors::within(rowRanges(vcf), {
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
      if (verbose)(message(paste0("parsing AD field in VCF.")))
      vcf.ncol <- ncol(geno(vcf)$AD)
      if (verbose)(message(paste0(vcf.ncol, " columns found in the VCF.")))
      vcf.names <- colnames(geno(vcf)$AD)
      ## grab the last item in the array if ids not specified... presumably tumor
        if(vcf.ncol > 1){
          if (verbose)(message(paste0("parsing tumor and normal alt/ref counts.")))
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
        } else {
          if (verbose)(message(paste0("parsing only tumor alt/ref counts.")))
        adep = setnames(as.data.table(geno(vcf)$AD[, vcf.ncol, 1:2]),
                        c("ref", "alt"))
        gt = geno(vcf)$GT
      }
    }
    else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in%
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
    ## mcols(rr) = BiocGenerics::cbind(mcols(rr), ann, gt = gt[, 1])
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
