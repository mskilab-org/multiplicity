#' @name snvplicity
#' @title Predecated name for multiplicity
#' @description predecated name for power function to run multiplicity
#' @details Please use multiplicity() instead
#' @export
snvplicity <- function(...) {
  multiplicity(...)
}

#' @name multiplicity
#' @title Converts counts to copies.
#'
#' @param somatic_snv Path to somatic SNV file.
#' @param germline_snv Path to germline SNV file.
#' @param het_pileups_wgs Path to het_pileups_wgs file.
#' @param jabba_rds Path to jabba file
#' @param tumor_cbs Path to segmented drycleaned coverage file. will override supplied tumor_dryclean file
#' @param tumor_dryclean Path to drycleaned coverage file. if provided, will rescale REF/ALT values to expected base coverage as determined by dryclean.
#' @param dryclean_field specification of field in tumor_dryclean gRanges used as rescaled binned coverage.
#' @param read_size estimated average read size (default: 151 for Illumina sequencers.)
#' @param tumor_name Expected name of tumor as annotated in the VCF; highly suggested this is provided for accurate VCF parsing
#' @param normal_name Expected name of normal as annotated in the VCF; highly suggested this is provided for accurate VCF parsing
#' @param mask if mask gRanges is provided, loci within ranges will be excluded from analysis
#' @param inferred_sex enum of ["M", "F"] that indicates inferred of true sex of patient sample (default = NA); optional input; will be inferred from coverage/JaBbA graph if not provided
#' @param tau_in_gamma Use tau(T) or tau_hat(F) in computation of gamma? tau_hat is the average copy number only of loci; tau is simply ploidy
#' @param filterpass process only FILTER == PASS variants?
#' @param purity Purity of inputted jabba_rds (optional if metadata of gGraph contains purity)
#' @param ploidy Ploidy of inputted jabba_rds (optional if metadata of gGraph contains ploidy)
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @param verbose verbose output?
#' @return Returns a GRanges with counts and converted copies
#' @export
multiplicity <- function(somatic_snv = NULL,
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
                         filterpass = FALSE,
                         tau_in_gamma = FALSE,
                         purity = NULL,
                         ploidy = NULL,
                         verbose = FALSE) {
  ### if any filepaths are /dev/null, turn them into true NULL characters.
  vars <- c("somatic_snv", "germline_snv", "het_pileups_wgs", "tumor_dryclean", "tumor_cbs")
  for (var in vars) {
    assign(var, normalize_path(get(var)))
  }

  if (verbose) {
    message("reading in jabba file")
  }

  if (is.null(somatic_snv) & is.null(germline_snv) & is.null(het_pileups_wgs)) {
    stop("Somatic VCF, Germline VCF, and/or HetSNPs file must be provided.")
  }

  gg <- gG(jabba = jabba_rds)

  read_size <- as.numeric(read_size)

  is_cov_character = is.character(tumor_dryclean)
  is_cov_len_one = NROW(tumor_dryclean) == 1
  is_cov_null = is.null(tumor_dryclean)
  is_cov_na = is_cov_len_one && (is.na(tumor_dryclean) || identical(tumor_dryclean, "NA"))
  is_cov_existent = is_cov_character && is_cov_len_one && file.exists(tumor_dryclean)
  is_cov_nonexistent = is_cov_character && is_cov_len_one && !file.exists(tumor_dryclean)
  is_cov_invalid = is_cov_character && !is_cov_len_one
  if (!is_cov_null && !is_cov_na) {
    if (is_cov_nonexistent) {
      stop("Path to coverage provided to 'tumor_dryclean' does not exist")
    } else if (is_cov_invalid) {
      stop("Path to coverage provided to 'tumor_dryclean' is a character but not length 1")
    }
    is_cov_rds = is_cov_existent && grepl("rds$", tumor_dryclean)
    is_cov_possibly_txt = is_cov_existent && !is_cov_rds
    if (is_cov_rds) {
      dryclean.cov = readRDS(tumor_dryclean) 
    } else if (is_cov_possibly_txt) {
      dryclean.cov = fread(tumor_dryclean)
    } else {
      dryclean.cov = tumor_dryclean
    }
    is_cov_table = inherits(dryclean.cov, "data.frame")
    if (is_cov_table) {
      dryclean.cov = gUtils::dt2gr(dryclean.cov)
      GenomeInfoDb::seqlevelsStyle(dryclean.cov) = "NCBI" ## better than gr.nochr()
    }
    is_cov_granges = inherits(dryclean.cov, "GRanges")
    if (!is_cov_granges) {
      stop("Provided coverage must be a path to a file coercible to GRanges, or a GRanges/tabular ranged format")
    }
    cov.vector <- mcols(dryclean.cov)[names(mcols(dryclean.cov)) %in% dryclean_field][, 1]
    mcols(dryclean.cov) <- NULL
    mcols(dryclean.cov)$bincov <- cov.vector
    mcols(dryclean.cov)$avg_basecov <- dryclean.cov$bincov * 2 * read_size / width(dryclean.cov)
  } else {
    dryclean.cov <- NULL
  }

  # CBS WILL OVERRIDE DRYCLEAN INPUTS
  # CBS reduces variance relative to dryclean rescaling alone
  is_seg_character = is.character(tumor_cbs)
  is_seg_len_one = NROW(tumor_cbs) == 1
  is_seg_null = is.null(tumor_cbs)
  is_seg_na = is_seg_len_one && (is.na(tumor_cbs) || identical(tumor_cbs, "NA"))
  is_seg_existent = is_seg_character && is_seg_len_one && file.exists(tumor_cbs)
  is_seg_nonexistent = is_seg_character && is_seg_len_one && !file.exists(tumor_cbs)
  is_seg_invalid = is_seg_character && !is_seg_len_one
  if (!is_seg_null && !is_seg_na) {
    if (is_seg_nonexistent) {
      stop("Path to segmentation provided to 'tumor_cbs' does not exist")
    } else if (is_seg_invalid) {
      stop("Path to segmentation provided to 'tumor_cbs' is invalid (length == 0)")
    }
    is_seg_rds = is_seg_existent && grepl("rds$", tumor_cbs)
    is_seg_possibly_txt = is_seg_existent && !is_seg_rds
    if (is_seg_rds) {
      cbs.cov = readRDS(tumor_cbs) 
    } else if (is_seg_possibly_txt) {
      cbs.cov = fread(tumor_cbs)
    } else {
      cbs.cov = tumor_cbs
    }
    is_seg_table = inherits(cbs.cov, "data.frame")
    if (is_seg_table) {
      cbs.cov = gUtils::dt2gr(cbs.cov)
      GenomeInfoDb::seqlevelsStyle(cbs.cov) = "NCBI" ## better than gr.nochr()
    }
    is_seg_granges = inherits(cbs.cov, "GRanges")
    if (!is_seg_granges) {
      stop("Provided segmentation must be a path to a file coercible to GRanges, or a GRanges/tabular ranged format")
    }
    cbs.vector <- mcols(cbs.cov)[names(mcols(cbs.cov)) %in% "seg.mean"][, 1]
    mcols(cbs.cov) <- NULL
    if (verbose) message("Assuming segmentation seg.mean is log-scaled, converting to linear")
    mcols(cbs.cov)$bincov <- exp(cbs.vector)
    cbs.cov <- gr.tile(seqlengths(cbs.cov), 1000) %>% gr.val(cbs.cov, "bincov")
    cbs.cov[which(is.na(cbs.cov$bincov))]$bincov <- 0
    mcols(cbs.cov)$avg_basecov <- cbs.cov$bincov * 2 * read_size / width(cbs.cov)
    if (exists("dryclean.cov")) {
      message("Provided segmentation superseding drycleaned binned coverage")
    }
    dryclean.cov = cbs.cov
  } else {
    cbs.cov <- NULL
  }


  if (!(is.na(mask) || is.null(mask))) {
    mask <- readRDS(mask)
  }

  if (is.na(inferred_sex)) {
    ## constitutional_cn assignment
    ## c_subj == 1 for major allele
    ## c_subj == 1 for autosomes and X chromosome in females, 0 for X and Y in males
    ncn.x <- gg$nodes$dt[
      (seqnames == "X" |
        seqnames == "chrX" |
        seqnames == "23" |
        seqnames == "chr23"), # MSK-FACETS calls chrX as chr23
      weighted.mean(cn,
        w = end - start + 1,
        na.rm = TRUE
      )
    ]
    message("mean CN of X: ", ncn.x)
    if (ncn.x < 1.4) {
      inferred_sex <- "M"
      message("sex determination: XY")
    } else {
      inferred_sex <- "F"
      message("sex determination:: XX")
    }
  } else if (grepl("^m", inferred_sex, ignore.case = T)) {
    message("sex provided: XY")
    inferred_sex <- "M"
  } else if (grepl("^f", inferred_sex, ignore.case = T)) {
    message("sex determination:: XX")
    inferred_sex <- "F"
  } else {
    stop("sex not provided and could not be inferred from provided inputs.")
  }


  ploidy <- if (!is.null(ploidy)) {
    ploidy
  } else {
    gg$meta$ploidy
  }
  purity <- if (!is.null(purity)) {
    purity
  } else if (!is.null(gg$meta$purity)) {
    gg$meta$purity
  }

  jab <- dt2gr(gr2dt(gg$nodes$gr)[seqnames == 23, seqnames := "X"])
  cn.gr <- jab[as.logical(strand(jab) == "+")]

  somatic.variants <- NULL
  germline.variants <- NULL
  het.pileups <- NULL

  if (verbose) {
    message("succesfully read in JaBbA graph!")
  }

  if (!is.na(somatic_snv) && !is.null(somatic_snv)) {
    message("Processing somatic variants.")

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
      tau_in_gamma = tau_in_gamma,
      mask = mask,
      major_gamma_coeff = 2,
      minor_gamma_coeff = 0
    )
  }

  if (!is.na(germline_snv) && !is.null(germline_snv)) {
    message("Processing germline variants.")
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
      tau_in_gamma = tau_in_gamma,
      mask = mask,
      major_gamma_coeff = 1,
      minor_gamma_coeff = 1
    )
  }

  if (!is.na(het_pileups_wgs) && !is.null(het_pileups_wgs)) {
    message("Processing heterozygous SNPs.")
    het.pileups <- transform_hets(
      hets = het_pileups_wgs,
      cn_gr = cn.gr,
      dryclean.cov = dryclean.cov,
      basecov_field = "avg_basecov",
      inferred_sex = inferred_sex,
      purity = purity,
      ploidy = ploidy,
      verbose = verbose,
      tau_in_gamma = tau_in_gamma,
      mask = mask
    )
  }

  return(list(somatic.variants, germline.variants, het.pileups))
}


#' @name transform_snv
#' @title function to inhale SNV-laden VCF and transform counts to to estimated altered copies using linear transformation via beta, gamma conversions
#'
#' @param vcf Path to SnpEff-annotated VCF
#' @param dryclean.cov Tumor drycleaned coverage for denoising/rescaling of REF/ALT counts (OPTIONAL; if provided will denoise)
#' @param basecov_field field name of dryclean.cov that indicates average Base Coverage per bin
#' @param tumor_id Tumor name as annotated in the VCF (highly suggested this is provided)
#' @param normal_id Normal name as annotated in the VCF (highly suggested this is provided)
#' @param purity Inferred purity input for conversion
#' @param ploidy Inferred ploidy input for conversion
#' @param inferred_sex enum of ["M", "F"] that indicates inferred of true sex of patient sample
#' @param filterpass process only FILTER == PASS variants?
#' @param tau_in_gamma Use tau(T) or tau_hat(F) in computation of gamma? tau_hat is the average copy number only of loci; tau is simply ploidy
#' @param mask if mask gRanges is provided, loci within ranges will be excluded from analysis
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @param verbose verbose output?
#' @return GRanges of transformed read counts
#' @export
transform_snv <- function(vcf,
                          cn_gr,
                          dryclean.cov = NULL,
                          basecov_field = "avg_basecov",
                          major_gamma_coeff = 1,
                          minor_gamma_coeff = 1,
                          tumor_id = NA,
                          normal_id = NA,
                          ploidy,
                          purity,
                          inferred_sex,
                          filterpass = FALSE,
                          snpeff_path,
                          verbose = TRUE,
                          tau_in_gamma = FALSE,
                          mask = mask) {
  snv <- parsesnpeff(
    vcf = vcf,
    snpeff_path = snpeff_path,
    coding_alt_only = FALSE,
    filterpass = filterpass,
    tumor_id = tumor_id,
    normal_id = normal_id,
    keepfile = FALSE,
    altpipe = TRUE,
    verbose = verbose
  )

  names(snv) <- NULL
  fields.to.carry <- character()
  snv.filtered <- gr2dt(snv)[, .SD[1], by = c("seqnames", "start", "end")]
  snv.filtered <- snv.filtered[, variant.g := paste0(REF, ">", ALT)] %>% dt2gr()

  if (!is.null(mask)) {
    snv.filtered <- gr.val(snv.filtered, mask, "blacklisted", na.rm = T) %Q%
      (!is.na(blacklisted)) %Q%
      (blacklisted == F)
  }

  # normalization
  unique.snv <- snv.filtered[, c("ref", "alt")] %>%
    gr.nochr() %Q%
    (!seqnames == c("Y")) %>%
    unique()

  if (!(is.null(dryclean.cov))) {
    unique.snv <- gr.val(unique.snv, dryclean.cov, basecov_field, na.rm = T)
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

  m <- length(unique.snv$ref + unique.snv$alt) 
  sf <- sum(unique.snv$major.count + unique.snv$minor.count, na.rm = T) / m
  unique.snv$major.count <- unique.snv$major.count / sf
  unique.snv$minor.count <- unique.snv$minor.count / sf

  ncn.vec <- rep(2, length(unique.snv))

  if (inferred_sex == "M") { # if male
    if (verbose) message("Adjusting ncn for XY")
    ncn.vec[which(as.character(seqnames(unique.snv)) %in% c("chrX", "chrY", "X", "Y"))] <- 1
  }

  values(unique.snv)[, "minor_constitutional_cn"] <- ncn.vec - 1
  unique.snv <- gr.val(unique.snv, cn_gr, "cn", na.rm = T)


  ###### APPLICATION OF FORMULA #####

  tau_hat <- mean(unique.snv$cn, na.rm = T)
  tau <- if (tau_in_gamma) {
    ploidy
  } else {
    tau_hat
  }
  alpha <- purity
  beta <- alpha / (alpha * ploidy + 2 * (1 - alpha))
  gamma <- 2 * (1 - alpha) / (alpha * tau + 2 * (1 - alpha))

  if (verbose) message(paste0("average total CN of loci: ", tau_hat))
  if (verbose) message(paste0("ploidy of tumor sample: ", ploidy))
  if (verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
  if (verbose) message("applying transformation")

  if (major_gamma_coeff > 2) {
    ncn.add <- major_gamma_coeff - 1
  } else {
    ncn.add <- 0
  }

  mcols(unique.snv)$major_snv_copies <-
    (2 * unique.snv$major.count -
      (gamma * (unique.snv$minor_constitutional_cn + ncn.add))) / (2 * beta)

  mcols(unique.snv)$minor_snv_copies <-
    (2 * unique.snv$minor.count -
      (gamma * unique.snv$minor_constitutional_cn * minor_gamma_coeff)) / (2 * beta)

  mcols(unique.snv)$total_snv_copies <-
    unique.snv$major_snv_copies + unique.snv$minor_snv_copies

  variants <- gr.val(
    snv.filtered %>% gr.nochr(),
    unique.snv[, c(
      fields.to.carry,
      "major.count",
      "minor.count",
      "total_snv_copies",
      "major_snv_copies",
      "minor_snv_copies",
      "cn"
    )],
    c(
      fields.to.carry,
      "major.count", "minor.count",
      "total_snv_copies",
      "major_snv_copies",
      "minor_snv_copies",
      "cn"
    )
  ) %>% gr2dt()

  variants[alt >= ref, altered_copies := major_snv_copies]
  variants[alt < ref, altered_copies := minor_snv_copies]
  # variants[total_copies <= 0, total_copies := 0]
  variants[, VAF := alt / (alt + ref)]
  variants <- variants %>% dt2gr()

  return(variants)
}


#' @name transform_hets
#' @title function to inhale SNP-laden gRanges and transform counts to to estimated altered copies using linear transformation via beta, gamma conversions
#'
#' @param hets heterozygous SNPs GRanges
#' @param dryclean.cov Tumor drycleaned coverage for denoising/rescaling of REF/ALT counts (OPTIONAL; if provided will denoise)
#' @param basecov_field field name of dryclean.cov that indicates average `Base Coverage` per bin
#' @param tumor_id Tumor name as annotated in the VCF (highly suggested this is provided; otherwise inferences are made about ordering in the VCF)
#' @param normal_id Normal name as annotated in the VCF (highly suggested this is provided; otherwise inferences are made about ordering in the VCF)
#' @param purity Inferred purity input for conversion
#' @param ploidy Inferred ploidy input for conversion
#' @param inferred_sex enum of ["M", "F"] that indicates inferred of true sex of patient sample
#' @param tau_in_gamma Use tau(T) or tau_hat(F) in computation of gamma? tau_hat is the average copy number only of loci; tau is simply ploidy
#' @param mask if mask gRanges is provided, loci within ranges will be excluded from analysis
#' @return GRanges of transformed read counts
#' @export
transform_hets <- function(hets,
                           cn_gr,
                           dryclean.cov = NULL,
                           basecov_field = "avg_basecov",
                           ploidy,
                           purity,
                           inferred_sex,
                           verbose = TRUE,
                           tau_in_gamma = FALSE,
                           mask = mask) {
  hets.gr <- fread(hets) %>% dt2gr()

  if (!is.null(mask)) {
    hets.gr <- gr.val(hets.gr, mask, "blacklisted", na.rm = T) %Q%
      (!is.na(blacklisted)) %Q%
      (blacklisted == F)
  }

  # normalization
  unique.snv <- hets.gr[, c("ref.count.t", "alt.count.t")] %>%
    gr.nochr() %Q%
    (!seqnames == c("Y")) %>%
    unique()

  fields.to.carry <- c("ref", "alt")
  names(mcols(unique.snv)) <- fields.to.carry

  ## get rid of lower peak values mostly for targeted panel values
  ## (unique.snv$ref + unique.snv$alt) %>% density(bw = 100, adjust = 1.5) -> total.read.density
  ## data.table(x = total.read.density$x, y = total.read.density$y) -> density.data
  ## peaks <- density.data[x %in% density.data[y %in% sort(y, decreasing = TRUE)[1:2], x]]$x
  ## peaks <- sort(peaks)
  ## valley <- density.data[x < peaks[1] & x > peaks[2]]


  if (!(is.null(dryclean.cov))) {
    unique.snv <- gr.val(unique.snv, dryclean.cov, basecov_field, na.rm = T)
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

  m <- length(unique.snv$ref + unique.snv$alt) 
  sf <- sum(unique.snv$major.count + unique.snv$minor.count, na.rm = T) / m
  unique.snv$major.count <- unique.snv$major.count / sf
  unique.snv$minor.count <- unique.snv$minor.count / sf

  ncn.vec <- rep(2, length(unique.snv))

  if (inferred_sex == "M") { # if male
    if (verbose) message("Adjusting ncn for XY")
    ncn.vec[which(as.character(seqnames(unique.snv)) %in% c("chrX", "chrY", "X", "Y"))] <- 1
  }

  values(unique.snv)[, "minor_constitutional_cn"] <- ncn.vec - 1
  unique.snv <- gr.val(unique.snv, cn_gr, "cn", na.rm = T)
  unique.snv$cn <- as.integer(round(unique.snv$cn))

  ###### APPLICATION OF FORMULA #####

  tau_hat <- mean(unique.snv$cn, na.rm = T)
  tau <- if (tau_in_gamma) {
    ploidy
  } else {
    tau_hat
  }
  alpha <- purity
  beta <- alpha / (alpha * ploidy + 2 * (1 - alpha))
  gamma <- 2 * (1 - alpha) / (alpha * tau + 2 * (1 - alpha))

  if (verbose) message(paste0("average total CN of heterozygous SNPs: ", tau_hat))
  if (verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
  if (verbose) message("applying transformation")

  mcols(unique.snv)$major_snv_copies <-
    (2 * unique.snv$major.count - (gamma * 1)) / (2 * beta)

  mcols(unique.snv)$minor_snv_copies <-
    (2 * unique.snv$minor.count - (gamma * unique.snv$minor_constitutional_cn)) / (2 * beta)

  ## mcols(unique.snv)$minor_snv_copies =
  ##                   (2 * unique.snv$minor.count - (gamma * 1)) / (2 * beta)

  mcols(unique.snv)$total_snv_copies <-
    unique.snv$major_snv_copies + unique.snv$minor_snv_copies

  variants <- gr.val(
    hets.gr %>% gr.nochr(),
    unique.snv[, c(
      fields.to.carry,
      "major.count",
      "minor.count",
      "total_snv_copies",
      "major_snv_copies",
      "minor_snv_copies",
      "cn"
    )],
    c(
      fields.to.carry,
      "major.count", "minor.count",
      "total_snv_copies",
      "major_snv_copies",
      "minor_snv_copies",
      "cn"
    )
  ) %>% gr2dt()

  variants[alt >= ref, altered_copies := major_snv_copies]
  variants[alt < ref, altered_copies := minor_snv_copies]
  # variants[total_copies <= 0, total_copies := 0]
  variants[, VAF := alt / (alt + ref)]
  variants <- variants %>% dt2gr()

  return(variants)
}

#' @title Parse SnpEff vcfs.
#' @description Ingests and processes a snpeff annotated vcf. #' VCF annotations from a snpeff file are parsed,  with optional functionality provided to select protein coding only mutations or filter = PASS events. AD and DP fields are also collected for downstream applications.
#' @param vcf path to snpeff vcf
#' @param snpeff_path Path to unzipped SnpEff toolkit
#' @param tumor_id Tumor name as annotated in the VCF
#' @param normal_id Normal name as annotated in the VCF
#' @param filterpass process only FILTER == PASS variants?
#' @param coding_alt_only process only coding alterations?
#' @param gr todo
#' @param keepfile todo
#' @param altpipe parsesnpeff via khtools::grok_vcf?
#' @param geno todo
#' @param debug debug function?
#' @param verbose verbose output?
#' @param bcftools Path to bcftools toolkit
#' @return GRangesList of breakpoint pairs with junctions that overlap removed
#' @export
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
  snpeff_path = "~/modules/SnpEff/source/snpEff",
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
        if(verbose) (message("Coding alterations only and FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, filt, bcftools, tmp.path)  
      } else {
        if(verbose) (message("Coding alterations only."))
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
    vcf = VariantAnnotation::readVcf(tmp.path)
    rr = MatrixGenerics::rowRanges(vcf)
    rr$REF = as.character(rr$REF)
    vcf_geno_lst = VariantAnnotation::geno(vcf)
    ann = as.data.table(data.table::tstrsplit(unlist(VariantAnnotation::info(vcf)$ANN),
      "\\|"))[, 1:15, with = FALSE, drop = FALSE]
    fn = c("allele", "annotation", "impact", "gene", "gene_id",
      "feature_type", "feature_id", "transcript_type", 
      "rank", "variant.c", "variant.p", "cdna_pos", "cds_pos", 
      "protein_pos", "distance")
    data.table::setnames(ann, fn)
    if ("AD" %in% names(vcf_geno_lst)) {
      if (verbose)(message(paste0("parsing AD field in VCF.")))
      vcf.ncol <- ncol(vcf_geno_lst$AD)
      if (verbose)(message(paste0(vcf.ncol, " columns found in the VCF.")))
      vcf.names <- colnames(vcf_geno_lst$AD)
      if(vcf.ncol > 1) {
        if (verbose)(message(paste0("parsing tumor and normal alt/ref counts.")))
        ## grab the last item in the array if ids not specified... presumably tumor
        tumor_col_id = vcf.ncol ## assume last column bydefault
        normal_col_id = NA_integer_
        is_id_matchable = function(id) {
          is_null_or_na = is.null(id) || is.na(id) || identical(id, "NA")
          return(!is_null_or_na)
        }
        if (is_id_matchable(tumor_id) && any(tumor_id %in% vcf.names)) {
          tumor_col_id = base::match(tumor_id, vcf.names)
        }
        if (is_id_matchable(normal_id) && any(normal_id %in% vcf.names)) {
          normal_col_id = base::match(normal_id, vcf.names)
        }
        adep = data.table::transpose(vcf_geno_lst$AD[,tumor_col_id])
        adep = as.data.table(adep)
        adep = base::subset(
          adep,
          select = 1:2
        )
        data.table::setnames(adep, c("ref", "alt"))
        adep$normal.ref = NA_integer_
        adep$normal.alt = NA_integer_
        if (!is.na(normal_col_id)) {
          adep.n = data.table::transpose(vcf_geno_lst$AD[,normal_col_id])
          adep$normal.ref = adep.n[[1]]
          adep$normal.alt = adep.n[[2]]
        }
        gt = vcf_geno_lst$GT
      } else {
        if (verbose)(message(paste0("parsing only tumor alt/ref counts.")))
        adep = vcf_geno_lst$AD[, vcf.ncol]
        adep = data.table::transpose(adep)
        adep = setnames(as.data.table(adep), c("ref", "alt"))
        gt = vcf_geno_lst$GT
      }
    } else if (all(c("AU", "GU", "CU", "TU", "TAR", "TIR") %in%
                     c(names(vcf_geno_lst)))) {
      if(verbose)(message("parsing AU/GU/CU/TAR/TIR fields in VCF."))
      this.col = dim(vcf_geno_lst[["AU"]])[2]
      d.a = vcf_geno_lst[["AU"]][, , 1, drop = F][, this.col, 1]
      d.g = vcf_geno_lst[["GU"]][, , 1, drop = F][, this.col, 1]
      d.t = vcf_geno_lst[["TU"]][, , 1, drop = F][, this.col, 1]
      d.c = vcf_geno_lst[["CU"]][, , 1, drop = F][, this.col, 1]
      mat = cbind(A = d.a, G = d.g, T = d.t, C = d.c)
      refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(mat))
      refid = ifelse(!isSNV(vcf), NA_integer_, refid)
      ## Note this assumes that multiallelic splitting via bcftools norm -m -any is run upstream (this breaks some fields if run after vcfOneLine.pl)
        ## otherwise you can't just unlist the ALT field if multiallelics are present
      altid = match(as.character(unlist(VariantAnnotation::fixed(vcf)$ALT)), colnames(mat))
      altid = ifelse(!isSNV(vcf), NA_integer_, altid)
      refsnv = mat[cbind(seq_len(nrow(mat)), refid)]
      altsnv = mat[cbind(seq_len(nrow(mat)), altid)]
      this.icol = dim(vcf_geno_lst[["TAR"]])[2]
      refindel = d.tar = vcf_geno_lst[["TAR"]][, , 1, drop = F][, this.icol, 1]
      altindel = d.tir = vcf_geno_lst[["TIR"]][, , 1, drop = F][, this.icol, 1]
      try2({
        n.d.a = vcf_geno_lst[["AU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.g = vcf_geno_lst[["GU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.t = vcf_geno_lst[["TU"]][, , 1, drop = F][, this.col - 1, 1]
        n.d.c = vcf_geno_lst[["CU"]][, , 1, drop = F][, this.col - 1, 1]
        n.mat = cbind(A = n.d.a, G = n.d.g, T = n.d.t, C = n.d.c)
        n.refid = match(as.character(VariantAnnotation::fixed(vcf)$REF), colnames(n.mat))
        n.refid = ifelse(!isSNV(vcf), NA_integer_, refid)
        ## Note this assumes that multiallelic splitting via bcftools norm -m -any is run upstream (this breaks some fields if run after vcfOneLine.pl)
        ## otherwise you can't just unlist the ALT field if multiallelics are present
        n.altid = match(as.character(unlist(VariantAnnotation::fixed(vcf)$ALT)), colnames(n.mat))
        n.altid = ifelse(!isSNV(vcf), NA_integer_, altid)
        n.refsnv = n.mat[cbind(seq_len(nrow(n.mat)), refid)]
        n.altsnv = n.mat[cbind(seq_len(nrow(n.mat)), altid)]
        n.refindel = n.d.tar = vcf_geno_lst[["TAR"]][, , 1, drop = F][, this.icol - 1, 1]
        n.altindel = n.d.tir = vcf_geno_lst[["TIR"]][, , 1, drop = F][, this.icol - 1, 1]
      })
      adep = data.table(ref = coalesce(refsnv, refindel),
        alt = coalesce(altsnv, altindel))
      try2({
        adep.n = data.table(normal.ref = coalesce(n.refsnv, n.refindel),
          normal.alt = coalesce(n.altsnv, n.altindel))
        adep =  adep %>% cbind(adep.n)
      })
      gt = data.frame(GT = rep_len("", NROW(vcf)))
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
rand.string <- function(n = 1, length = 12) {
  randomString <- c(1:n) # initialize vector
  for (i in 1:n) {
    randomString[i] <- paste(
      sample(c(0:9, letters, LETTERS),
        length,
        replace = TRUE
      ),
      collapse = ""
    )
  }
  return(randomString)
}

#' @name try2
#' @title wrapper around tryCatch - robust to parallel:: functions; A slightly more robust version of try that works within the parallel:: set of functions that pre-deploy a cluster.
try2 <- function(expr, ..., finally) {
  tryCatch(expr,
    error = function(e) {
      msg <- structure(paste(conditionMessage(e), conditionCall(e), sep = "\n"), class = "err")
      cat("Error: ", msg, "\n\n")
      return(msg)
    },
    finally = finally,
    ... = ...
  )
}

#' @name normalize_path
#' @title easy function that returns NULL if file path set to /dev/null
normalize_path <- function(path) {
  is_character = is.character(path)
  is_length_one = NROW(path) == 1
  is_a_potential_path = is_character && is_length_one
  is_invalid_input = is_character && !is_length_one
  out = path
  if (is_a_potential_path && !file.exists(path)) stop("Provided path does not exist")
  if (is_invalid_input) stop("Provided invalid paths (not length == 1)")
  if (!(is_a_potential_path)) return(out)
  if (path == "/dev/null") return(NULL)
  ## After all this.. this will return variable "path" (could be single length character or NULL)
  return(out)
}

preprocess_snps <- function(gr) {

}
