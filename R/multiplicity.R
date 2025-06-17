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
#' @param modeltype Model type to use for multiplicity calculations. Options are "unified" (default) or "separate". "unified" uses a single model for all SNVs, while "separate" uses different models for somatic and germline SNVs.
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
             snpeff_path = system.file("extdata", "snpeff_scripts", package = "multiplicity"),
             tumor_name = NULL,
             normal_name = NULL,
             filterpass = FALSE,
             tau_in_gamma = FALSE,
             purity = NULL,
             ploidy = NULL,
             modeltype = "unified",
             verbose = FALSE) {

              preprocessed_inputs <- preprocess_multiplicity_inputs(
                somatic_snv = somatic_snv,
                germline_snv = germline_snv,
                het_pileups_wgs = het_pileups_wgs,
                tumor_cbs = tumor_cbs,
                tumor_dryclean = tumor_dryclean,
                dryclean_field = dryclean_field,
                jabba_rds = jabba_rds,
                mask = mask,
                inferred_sex = inferred_sex,
                read_size = read_size,
                ploidy = ploidy,
                purity = purity,
                verbose = verbose
              )

              somatic_variants <- NULL
              germline_variants <- NULL
              het_pileups <- NULL

              cn.gr <- preprocessed_inputs$cn.gr
              dryclean.cov <- preprocessed_inputs$dryclean.cov
              inferred_sex <- preprocessed_inputs$inferred_sex
              purity <- preprocessed_inputs$purity
              ploidy <- preprocessed_inputs$ploidy
              mask <- preprocessed_inputs$mask
              browser()

              if (modeltype == "unified") {
                ##### UNIFIED MODEL #####
                if (verbose) message("Using unified processing model.")
                
                if (!is.na(somatic_snv_path) && !is.null(somatic_snv_path)) {
                  if (verbose) message("Processing somatic variants (unified model using parsesnpeff).")
                  somatic_variants <- parsesnpeff(
                    vcf = somatic_snv,
                    snpeff_path = snpeff_path,
                    coding_alt_only = FALSE, # As per original commented-out code
                    filterpass = filterpass, # Uses the main filterpass parameter
                    tumor_id = tumor_name,
                    normal_id = normal_name,
                    keepfile = FALSE,
                    altpipe = TRUE,
                    verbose = verbose
                  )
                }
                
                if (!is.na(germline_snv_path) && !is.null(germline_snv_path)) {
                  if (verbose) message("Processing germline variants (unified model using parsesnpeff).")
                  germline_variants <- parsesnpeff(
                    vcf = germline_snv,
                    snpeff_path = snpeff_path,
                    coding_alt_only = FALSE, # As per original commented-out code
                    filterpass = filterpass, # Uses the main filterpass parameter
                    tumor_id = tumor_name,
                    normal_id = normal_name,
                    keepfile = FALSE,
                    altpipe = TRUE,
                    verbose = verbose
                  )
                }
                
                if (!is.na(het_pileups_wgs) && !is.null(het_pileups_wgs)) {
                  if (verbose) message("Processing heterozygous SNPs (unified model - basic load).")
                  if (file.exists(het_pileups_wgs_path)) {
                    # Assuming het_pileups_wgs_path is a file to be read by fread and converted
                    het_data <- data.table::fread(het_pileups_wgs)
                    if (inherits(het_data, "data.frame")) {
                      het_pileups <- gUtils::dt2gr(het_data) # Ensure gUtils is available
                    } else {
                      if (verbose) message("Loaded het_pileups data is not a data.frame, cannot convert to GRanges.")
                      het_pileups <- NULL
                    }
                    if (!is.null(het_pileups) && !inherits(het_pileups, "GRanges")) {
                      warning("Failed to convert het_pileups to GRanges in unified model.")
                      het_pileups <- NULL
                    }
                  } else {
                    if (verbose) message("Het pileups file not found for unified model: ", het_pileups_wgs_path)
                    het_pileups <- NULL
                  }
                }
                
                names(somatic_variants) <- NULL
                names(germline_variants) <- NULL
                names(het_pileups) <- NULL
                
                mcols(somatic_variants)$class <- "SOMATIC"
                mcols(germline_variants)$class <- "GERMLINE"
                mcols(het_pileups)$class <- "HET"
                
                somatic_dt <- gr2dt(somatic_variants)[, .SD[1], by = c("seqnames", "start", "end")][, c(
                  "seqnames", "start", "end", "strand", "width", "class", "ref", "alt",
                  intersect(c("normal.ref", "normal.alt"), names(gr2dt(somatic_variants)))
                ), with = FALSE]
                germline_dt <- gr2dt(germline_variants)[, .SD[1], by = c("seqnames", "start", "end")][, c(
                  "seqnames", "start", "end", "strand", "width", "class", "ref", "alt",
                  intersect(c("normal.ref", "normal.alt"), names(gr2dt(germline_variants)))
                ), with = FALSE]
                het_pileups_dt <- gr2dt(het_pileups)[, .SD[1], by = c("seqnames", "start", "end")][, .(
                  seqnames, start, end, strand, width, class,
                  ref = ref.count.t,
                  alt = alt.count.t,
                  normal.ref = if ("ref.count.n" %in% names(gr2dt(het_pileups))) ref.count.n else NA,
                  normal.alt = if ("alt.count.n" %in% names(gr2dt(het_pileups))) alt.count.n else NA
                )]
                
                variants <- rbind(somatic_dt, germline_dt, het_pileups_dt, fill = TRUE) %>% dt2gr()
                
                if (!is.null(mask)) {
                  variants <- gr.val(variants, mask, "blacklisted", na.rm = TRUE) %Q%
                    (!is.na(blacklisted)) %Q%
                    (blacklisted == FALSE) # Filter out blacklisted variants
                }
                
                unique.variants <- variants[, c("class", "ref", "alt")] %>% gr.nochr() %Q%
                  (!seqnames %in% c("Y", "MT")) %>% unique
                
                if (!is.null(dryclean.cov)) {
                  unique.variants <- gr.val(unique.variants, dryclean.cov, "avg_basecov", na.rm = TRUE)
                  unique.variants$ref_denoised <- unique.variants$ref * unique.variants$avg_basecov / (unique.variants$ref + unique.variants$alt)
                  unique.variants$alt_denoised <- unique.variants$alt * unique.variants$avg_basecov / (unique.variants$ref + unique.variants$alt)
                  
                  unique.variants$major.count <- pmax(unique.variants$ref_denoised, unique.variants$alt_denoised, na.rm = TRUE)
                  unique.variants$minor.count <- pmin(unique.variants$ref_denoised, unique.variants$alt_denoised, na.rm = TRUE)
                  fields.to.carry <- c("ref_denoised", "alt_denoised")
                } else {
                  unique.variants$major.count <- pmax(unique.variants$ref, unique.variants$alt, na.rm = TRUE)
                  unique.variants$minor.count <- pmin(unique.variants$ref, unique.variants$alt, na.rm = TRUE)
                }
                
                m <- length(unique.variants$ref + unique.variants$alt)
                sf <- sum(unique.variants$major.count + unique.variants$minor.count, na.rm = TRUE) / m
                unique.variants$major.count <- unique.variants$major.count / sf
                unique.variants$minor.count <- unique.variants$minor.count / sf
                
                ncn.vec <- rep(2, length(unique.variants))
                
                if (inferred_sex == "M") {
                  if (verbose) message("Adjusting ncn for XY")
                  ncn.vec[as.character(seqnames(unique.variants)) %in% c("X", "chrX", "23", "chr23", "Y", "chrY", "y")] <- 1
                }
                
                values(unique.variants)$minor_constitutional_cn <- ncn.vec - 1
                unique.variants <- gr.val(unique.variants, cn.gr, "cn", na.rm = TRUE)
                
                tau_hat <- mean(unique.variants$cn, na.rm = TRUE)
                tau <- if (tau_in_gamma) ploidy else tau_hat
                alpha <- purity
                beta <- alpha / (alpha * ploidy + 2 * (1 - alpha))
                gamma <- 2 * (1 - alpha) / (alpha * tau + 2 * (1 - alpha))
                
                if (verbose) message(paste0("average total CN of loci: ", tau_hat))
                if (verbose) message(paste0("ploidy of tumor sample: ", ploidy))
                if (verbose) message(paste0("purity: ", alpha, " beta: ", beta, " gamma: ", gamma))
                if (verbose) message("applying transformation")
                
                mcols(unique.variants)$major_gamma_coeff <- NA
                mcols(unique.variants)$major_gamma_coeff[unique.variants$class == "SOMATIC"] <- 2
                mcols(unique.variants)$major_gamma_coeff[unique.variants$class == "GERMLINE"] <- 1
                mcols(unique.variants)$major_gamma_coeff[unique.variants$class == "HET"] <- 1
                mcols(unique.variants)$minor_gamma_coeff <- NA
                mcols(unique.variants)$minor_gamma_coeff[unique.variants$class == "SOMATIC"] <- 0
                mcols(unique.variants)$minor_gamma_coeff[unique.variants$class == "GERMLINE"] <- 1
                mcols(unique.variants)$minor_gamma_coeff[unique.variants$class == "HET"] <- 1
                
                mcols(unique.variants)$ncn.add <- ifelse(
                  mcols(unique.variants)$major_gamma_coeff > 2,
                  mcols(unique.variants)$major_gamma_coeff - 1,
                  0
                )
                
                mcols(unique.variants)$major_snv_copies <-
                  (2 * unique.variants$major.count -
                     (gamma * (unique.variants$minor_constitutional_cn + unique.variants$ncn.add))) / (2 * beta)
                
                mcols(unique.variants)$minor_snv_copies <-
                  (2 * unique.variants$minor.count -
                     (gamma * unique.variants$minor_constitutional_cn * unique.variants$minor_gamma_coeff)) / (2 * beta)
                
                mcols(unique.variants)$total_snv_copies <-
                  unique.variants$major_snv_copies + unique.variants$minor_snv_copies
                
                unique.variants$altered_copies <- ifelse(unique.variants$alt >= unique.variants$ref,
                  unique.variants$major_snv_copies,
                  unique.variants$minor_snv_copies
                )
                unique.variants$VAF <- unique.variants$alt / (unique.variants$ref + unique.variants$alt)
                
                somatic_variants <- gr.val(somatic_variants, unique.variants %Q% (class == "SOMATIC"), c(fields.to.carry, "cn", "major_snv_copies", "minor_snv_copies", "total_snv_copies", "altered_copies", "VAF"))
                
                germline_variants <- gr.val(germline_variants, unique.variants %Q% (class == "GERMLINE"), c(fields.to.carry, "cn",
                                                                                                      "major_snv_copies", "minor_snv_copies", "total_snv_copies", "altered_copies", "VAF"))
                
                het_pileups <- gr.val(het_pileups, unique.variants %Q% (class == "HET"), c(fields.to.carry, "cn", "major_snv_copies", "minor_snv_copies", "total_snv_copies", "altered_copies", "VAF"))
                
              } else if (modeltype == "separate") {
                
                if (verbose) message("Using separate processing model.")
                
                if (!is.na(somatic_snv_path) && !is.null(somatic_snv_path)) {

                  if (verbose) message("Processing somatic variants (separate model using transform_snv).")
                  
                  somatic_variants <- transform_snv(
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
                    filterpass = filterpass, # Uses the main filterpass parameter
                    verbose = verbose,
                    tau_in_gamma = tau_in_gamma,
                    mask = mask,
                    major_gamma_coeff = 2, # As per original active code
                    minor_gamma_coeff = 0  # As per original active code
                  )
                }
                
                if (!is.na(germline_snv_path) && !is.null(germline_snv_path)) {

                  if (verbose) message("Processing germline variants (separate model using transform_snv).")

                  germline_variants <- transform_snv(
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
                    mask = mask,
                    filterpass = TRUE, # Hardcoded TRUE as per original active code for germline
                    verbose = verbose,
                    tau_in_gamma = tau_in_gamma,
                    major_gamma_coeff = 1, # As per original active code
                    minor_gamma_coeff = 1  # As per original active code
                  )
                }
                
                if (!is.na(het_pileups_wgs_path) && !is.null(het_pileups_wgs_path)) {

                  if (verbose) message("Processing heterozygous SNPs (separate model using transform_hets).")

                  het_pileups <- transform_hets(
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
              } else {
                stop(paste0("Invalid model_type specified: '", model_type, "'. Please choose 'unified' or 'separate'."))
              }

              return(list(
                somatic_variants = somatic_variants,
                germline_variants = germline_variants,
                het_pileups = het_pileups
              ))

#' @title Preprocess inputs for multiplicity
#' @description Helper function to preprocess and validate inputs for the main multiplicity function.
#' @param somatic_snv Path to somatic SNV file.
#' @param germline_snv Path to germline SNV file.
#' @param het_pileups_wgs Path to het_pileups_wgs file.
#' @param tumor_cbs Path to segmented drycleaned coverage file.
#' @param tumor_dryclean Path to drycleaned coverage file.
#' @param dryclean_field Specification of field in tumor_dryclean gRanges.
#' @param jabba_rds Path to jabba file or gGraph object.
#' @param mask Path to mask gRanges RDS file or gRanges object.
#' @param inferred_sex Enum of ["M", "F"].
#' @param read_size Estimated average read size.
#' @param ploidy Ploidy of the sample.
#' @param purity Purity of the sample.
#' @param verbose Verbose output?
#' @return A list containing preprocessed inputs:
#'         somatic_snv, germline_snv, het_pileups_wgs, cn.gr,
#'         dryclean.cov, inferred_sex, purity, ploidy, mask.
#' @keywords internal
preprocess_multiplicity_inputs <- function(somatic_snv,
                       germline_snv,
                       het_pileups_wgs,
                       tumor_cbs,
                       tumor_dryclean,
                       dryclean_field,
                       jabba_rds,
                       mask,
                       inferred_sex,
                       read_size,
                       ploidy,
                       purity,
                       verbose) {
  ### if any filepaths are /dev/null, turn them into true NULL characters.
  vars_to_normalize <- c("somatic_snv", "germline_snv", "het_pileups_wgs", "tumor_dryclean", "tumor_cbs")
  for (var_name in vars_to_normalize) {
  assign(var_name, normalize_path(get(var_name)))
  }

  if (verbose) {
    message("Starting preprocessing of inputs...")
  }

  if (is.null(somatic_snv) && is.null(germline_snv) && is.null(het_pileups_wgs)) {
    stop("Somatic VCF, Germline VCF, and/or HetSNPs file must be provided.")
  }

  gg <- NULL
  cn.gr <- NULL
  jab <- NULL

  is_jabba_character <- is.character(jabba_rds)
  is_jabba_len_one <- NROW(jabba_rds) == 1
  is_jabba_null <- is.null(jabba_rds) || identical(jabba_rds, base::nullfile())
  is_jabba_na <- is_jabba_len_one && (is.na(jabba_rds) || identical(jabba_rds, "NA"))

  if (!is_jabba_null && !is_jabba_na) {
  if (verbose) message("Processing JaBbA input...")
  is_jabba_existent <- is_jabba_character && is_jabba_len_one && file.exists(jabba_rds)
  is_jabba_nonexistent <- is_jabba_character && is_jabba_len_one && !file.exists(jabba_rds)
  is_jabba_invalid <- is_jabba_character && !is_jabba_len_one

  if (is_jabba_nonexistent) {
    stop("Path to JaBbA provided to 'jabba_rds' does not exist.")
  } else if (is_jabba_invalid) {
    stop("Path to jabba provided to 'jabba_rds' is a character but not length one.")
  }

  is_jabba_rds <- is_jabba_existent && grepl("rds$", jabba_rds, ignore.case = TRUE)
  is_jabba_list_like <- is.list(jabba_rds) && all(c("segstats", "junctions", "gtrack") %in% names(jabba_rds)) # Simplified check
  is_jabba_ggraph <- inherits(jabba_rds, "gGraph")

  if (is_jabba_rds || is_jabba_list_like) {
    gg <- gGnome::gG(jabba = jabba_rds)
  } else if (is_jabba_ggraph) {
    gg <- jabba_rds
  } else if (is_jabba_existent) { # Catch other file types if existent but not RDS
    stop("Path to jabba provided to 'jabba_rds' is an unsupported file type. Expecting an RDS, a gGraph object, or a list-like Jabba structure.")
  } else if (!is_jabba_character) { # jabba_rds is an object but not gGraph or recognized list
    stop("'jabba_rds' is an object of unrecognized type. Expecting a file path, a gGraph object, or a list-like Jabba structure.")
  }


  if (is.null(gg)) {
    stop("Failed to load or process jabba_rds into a gGraph object.")
  }

  jab <- gUtils::dt2gr(gUtils::gr2dt(gg$nodes$gr)[seqnames == 23, seqnames := "X"][seqnames == 24, seqnames := "Y"])
  GenomeInfoDb::seqlevelsStyle(jab) <- "NCBI"
  cn.gr <- jab[as.logical(strand(jab) == "+")]
  if (verbose) message("Successfully processed JaBbA graph!")
  }
  is_jabba_obj_present <- !is.null(gg) && !is.null(cn.gr)

  if (is.null(ploidy)) {
  if (!is_jabba_obj_present || is.null(gg$meta)) {
    stop("ploidy not provided, and JaBbA object (gg) or its metadata (gg$meta) is not available to infer ploidy. Please provide ploidy.")
  }
  fetched_ploidy <- try(gg$meta[["ploidy"]], silent = TRUE) # More direct for lists
  if (inherits(fetched_ploidy, "try-error") || is.null(fetched_ploidy)) {
    stop("ploidy not provided, and could not be found or is NULL in gGraph 'meta' field. You must provide a ploidy solution.")
  }
  ploidy <- fetched_ploidy
  }
  if (is.null(purity)) {
  if (!is_jabba_obj_present || is.null(gg$meta)) {
    stop("purity not provided, and JaBbA object (gg) or its metadata (gg$meta) is not available to infer purity. Please provide purity.")
  }
  fetched_purity <- try(gg$meta[["purity"]], silent = TRUE)
  if (inherits(fetched_purity, "try-error") || is.null(fetched_purity)) {
    stop("purity not provided, and could not be found or is NULL in gGraph 'meta' field. You must provide a purity solution.")
  }
  purity <- fetched_purity
  }

  read_size <- as.numeric(read_size)

  is_cov_character <- is.character(tumor_dryclean)
  is_cov_len_one <- NROW(tumor_dryclean) == 1
  is_cov_null <- is.null(tumor_dryclean) || identical(tumor_dryclean, base::nullfile())
  is_cov_na <- is_cov_len_one && (is.na(tumor_dryclean) || identical(tumor_dryclean, "NA"))
  dryclean.cov <- NULL

  if (!is_cov_null && !is_cov_na) {
  if (verbose) message("Processing tumor dryclean coverage...")
  is_cov_existent <- is_cov_character && is_cov_len_one && file.exists(tumor_dryclean)
  is_cov_nonexistent <- is_cov_character && is_cov_len_one && !file.exists(tumor_dryclean)
  is_cov_invalid <- is_cov_character && !is_cov_len_one

  if (is_cov_nonexistent) {
    stop("Path to coverage provided to 'tumor_dryclean' does not exist.")
  } else if (is_cov_invalid) {
    stop("Path to coverage provided to 'tumor_dryclean' is a character but not length 1.")
  }

  current_dryclean_obj <- NULL
  if (is_cov_existent && grepl("rds$", tumor_dryclean, ignore.case = TRUE)) {
    current_dryclean_obj <- readRDS(tumor_dryclean)
  } else if (is_cov_existent) { # Assumed text file if not RDS
    current_dryclean_obj <- data.table::fread(tumor_dryclean)
  } else if (!is_cov_character) { # tumor_dryclean is already an object
    current_dryclean_obj <- tumor_dryclean
  } else { # Should not be reached if previous checks are correct
    stop("Could not load tumor_dryclean input.")
  }

  if (inherits(current_dryclean_obj, "data.frame")) {
    current_dryclean_obj <- gUtils::dt2gr(current_dryclean_obj)
    GenomeInfoDb::seqlevelsStyle(current_dryclean_obj) <- "NCBI"
  }
  if (!inherits(current_dryclean_obj, "GRanges")) {
    stop("Provided tumor_dryclean must be a path to a file coercible to GRanges, or a GRanges/tabular ranged format.")
  }
  dryclean.cov <- current_dryclean_obj

  if (!dryclean_field %in% names(mcols(dryclean.cov))) {
    stop(paste0("dryclean_field '", dryclean_field, "' not found in tumor_dryclean metadata columns. Available columns: ", paste(names(mcols(dryclean.cov)), collapse = ", ")))
  }
  cov.vector <- mcols(dryclean.cov)[[dryclean_field]]
  mcols(dryclean.cov) <- NULL # Clear mcols before adding new ones
  mcols(dryclean.cov)$bincov <- cov.vector
  mcols(dryclean.cov)$avg_basecov <- dryclean.cov$bincov * 2 * read_size / width(dryclean.cov)
  }
  is_cov_obj_present <- !is.null(dryclean.cov)

  is_seg_character <- is.character(tumor_cbs)
  is_seg_len_one <- NROW(tumor_cbs) == 1
  is_seg_null <- is.null(tumor_cbs) || identical(tumor_cbs, base::nullfile())
  is_seg_na <- is_seg_len_one && (is.na(tumor_cbs) || identical(tumor_cbs, "NA"))
  cbs.cov.gr <- NULL # Use a distinct name for CBS GRanges
  cbs.vector <- NULL

  is_both_jabba_and_cov_present <- is_jabba_obj_present && is_cov_obj_present

  if (is_both_jabba_and_cov_present) {
  if (verbose) message("Using JaBbA segmentation and binned coverage to get segment level mean coverage.")
  if (!is.null(cn.gr) && !is.null(dryclean.cov) && "bincov" %in% names(mcols(dryclean.cov))) {
    seg_stats_result <- JaBbA:::segstats(cn.gr, dryclean.cov, field = "bincov")
    cbs.cov.gr <- cn.gr
    mcols(cbs.cov.gr) <- cbind(mcols(cn.gr), mcols(seg_stats_result)) # Combine metadata carefully
    cbs.vector <- cbs.cov.gr$mean
  } else {
    if (verbose) message("Inputs cn.gr or dryclean.cov (with bincov) are not suitable for JaBbA:::segstats. Skipping this step.")
  }
  }

  if (!is_both_jabba_and_cov_present && (!is_seg_null && !is_seg_na)) {
  if (verbose) message("Processing tumor CBS segmentation...")
  is_seg_existent <- is_seg_character && is_seg_len_one && file.exists(tumor_cbs)
  is_seg_nonexistent <- is_seg_character && is_seg_len_one && !file.exists(tumor_cbs)
  is_seg_invalid <- is_seg_character && !is_seg_len_one

  if (is_seg_nonexistent) {
    stop("Path to segmentation provided to 'tumor_cbs' does not exist.")
  } else if (is_seg_invalid) {
    stop("Path to segmentation provided to 'tumor_cbs' is invalid (not length one character).")
  }

  current_cbs_obj <- NULL
  if (is_seg_existent && grepl("rds$", tumor_cbs, ignore.case = TRUE)) {
    current_cbs_obj <- readRDS(tumor_cbs)
  } else if (is_seg_existent) {
    current_cbs_obj <- data.table::fread(tumor_cbs)
  } else if (!is_seg_character) {
    current_cbs_obj <- tumor_cbs
  } else {
    stop("Could not load tumor_cbs input.")
  }

  if (inherits(current_cbs_obj, "data.frame")) {
    current_cbs_obj <- gUtils::dt2gr(current_cbs_obj)
    GenomeInfoDb::seqlevelsStyle(current_cbs_obj) <- "NCBI"
  }
  if (!inherits(current_cbs_obj, "GRanges")) {
    stop("Provided tumor_cbs must be a path to a file coercible to GRanges, or a GRanges/tabular ranged format.")
  }
  cbs.cov.gr <- current_cbs_obj

  seg_mean_val <- NULL
  if ("seg.mean" %in% names(mcols(cbs.cov.gr))) {
    seg_mean_val <- mcols(cbs.cov.gr)[["seg.mean"]]
  }

  if (!is.null(seg_mean_val)) {
    cbs.vector <- seg_mean_val
    if (verbose) message("Assuming segmentation seg.mean is log-scaled, converting to linear.")
    cbs.vector <- exp(cbs.vector)
  } else if (is_cov_obj_present) { # seg.mean not found, but binned coverage exists
    if (verbose) message("seg.mean not found in tumor_cbs, recomputing segment means using provided binned coverage on CBS segments.")
    if (!is.null(cbs.cov.gr) && !is.null(dryclean.cov) && "bincov" %in% names(mcols(dryclean.cov))) {
    seg_stats_result <- JaBbA:::segstats(cbs.cov.gr, dryclean.cov, field = "bincov")
    mcols(cbs.cov.gr) <- cbind(mcols(cbs.cov.gr), mcols(seg_stats_result))
    cbs.vector <- cbs.cov.gr$mean
    } else {
    if (verbose) message("Inputs cbs.cov.gr or dryclean.cov (with bincov) are not suitable for JaBbA:::segstats. Cannot recompute segment means for CBS.")
    }
  } else { # seg.mean not found, and no binned coverage to recompute
    stop("No seg.mean values available in tumor_cbs, and no binned coverage (tumor_dryclean) provided to recompute.")
  }

  if (!is_jabba_obj_present && !is.null(cbs.cov.gr) && !is.null(cbs.vector)) {
    if (verbose) message("Deriving cn.gr from CBS segmentation as JaBbA was not provided/processed.")
    temp_mcols <- mcols(cbs.cov.gr) # Preserve original mcols
    mcols(cbs.cov.gr)$cov_val <- cbs.vector # Temporary column for rel2abs

    # rel2abs returns a GRanges with new 'cn' column
    abs_cn_gr <- skitools::rel2abs(cbs.cov.gr, "cov_val", purity = purity, ploidy = ploidy)

    cn.gr <- abs_cn_gr # This now has 'cn'
    jab <- cn.gr # jab is used for consistency if Jabba was primary source
    # mcols(cbs.cov.gr) <- temp_mcols # Restore original mcols if needed, or select specific ones
    is_jabba_obj_present <- TRUE # Mark as present because cn.gr is now derived
  }
  }

  is_seg_obj_processed <- !is.null(cbs.cov.gr) && !is.null(cbs.vector)
  if (is_seg_obj_processed) {
  if (verbose) message("Overriding dryclean.cov with processed CBS data.")
  # Create a new GRanges for the tiled coverage from CBS means
  tiled_cbs_gr <- NULL
  valid_seqlengths <- !is.null(seqlengths(cbs.cov.gr)) && all(!is.na(seqlengths(cbs.cov.gr))) && all(is.finite(seqlengths(cbs.cov.gr)))

  if (valid_seqlengths) {
    tiled_cbs_gr <- gUtils::gr.tile(seqlengths(cbs.cov.gr), 1000)

    # For gr.val, ensure cbs.cov.gr has only the 'bincov' (which is cbs.vector)
    temp_val_gr <- cbs.cov.gr
    mcols(temp_val_gr) <- S4Vectors::DataFrame(bincov = cbs.vector) # Use segment means as bincov

    tiled_cbs_gr <- gUtils::gr.val(tiled_cbs_gr, temp_val_gr, "bincov")
    tiled_cbs_gr$bincov[is.na(tiled_cbs_gr$bincov)] <- 0
    tiled_cbs_gr$avg_basecov <- tiled_cbs_gr$bincov * 2 * read_size / width(tiled_cbs_gr)
  } else {
    if (verbose) message("Cannot tile CBS coverage due to missing or invalid seqlengths. Using segment-level CBS data directly for dryclean.cov override.")
    # Fallback: use the segment-level cbs.cov.gr, ensuring it has bincov and avg_basecov
    tiled_cbs_gr <- cbs.cov.gr
    mcols(tiled_cbs_gr)$bincov <- cbs.vector # Ensure bincov is the segment mean
    mcols(tiled_cbs_gr)$avg_basecov <- tiled_cbs_gr$bincov * 2 * read_size / width(tiled_cbs_gr)
  }

  if (is_cov_obj_present && verbose) {
    message("Provided CBS segmentation is superseding binned coverage from tumor_dryclean.")
  }
  dryclean.cov <- tiled_cbs_gr
  is_cov_obj_present <- TRUE # dryclean.cov is now set from CBS
  }

  if (is.null(cn.gr)) {
    stop("cn.gr could not be derived from jabba_rds or tumor_cbs. This is required for multiplicity calculations.")
  }
  if (!is_cov_obj_present && verbose && (!is.null(somatic_snv) || !is.null(germline_snv) || !is.null(het_pileups_wgs))) {
    message("dryclean.cov is not available. SNV/Het count denoising/rescaling will be skipped.")
  }

  if (!(is.na(mask) || is.null(mask))) {
    if (is.character(mask) && file.exists(mask)) {
      if (verbose) message("Loading mask from RDS file: ", mask)
      mask <- readRDS(mask)
    } else if (inherits(mask, "GRanges")) {
      if (verbose) message("Using provided GRanges object as mask.")
    } else {
      warning("Mask is provided but is not a valid file path to an RDS or a GRanges object. It will be ignored.")
      mask <- NULL
    }
  }

  if (is.na(inferred_sex)) {
    if (!is.null(cn.gr) && "cn" %in% names(mcols(cn.gr))) {
      ncn.x <- gUtils::gr2dt(cn.gr)[
        (seqnames == "X" | seqnames == "chrX" | seqnames == "23" | seqnames == "chr23"),
        stats::weighted.mean(cn, w = width, na.rm = TRUE) # Use width calculation robustly
      ]
      if (verbose) message("Mean CN of X for sex inference: ", round(ncn.x, 2))
      if (is.na(ncn.x)) {
        warning("Could not calculate mean CN of X for sex inference (e.g., no X chr data or all CN values NA). Please provide 'inferred_sex'. Inferred_sex remains NA.")
      } else if (ncn.x < 1.4) {
        inferred_sex <- "M"
        if (verbose) message("Sex determination: XY")
      } else {
        inferred_sex <- "F"
        if (verbose) message("Sex determination: XX")
      }
    } else {
      warning("cn.gr (with 'cn' column) not available for sex inference. Please provide 'inferred_sex'. Inferred_sex remains NA.")
    }
  } else if (grepl("^m", inferred_sex, ignore.case = TRUE)) {
    inferred_sex <- "M"
    if (verbose) message("Sex provided: XY")
  } else if (grepl("^f", inferred_sex, ignore.case = TRUE)) {
    inferred_sex <- "F"
    if (verbose) message("Sex provided: XX")
  } else {
    if (verbose && !is.na(inferred_sex)) { # Provided but not M/F/NA
      message(paste0("Provided 'inferred_sex' (", inferred_sex, ") is not 'M', 'F', or NA. Attempting to proceed, but this may cause issues."))
    } else if (verbose && is.na(inferred_sex)) { # Was NA and could not be inferred
      message("Sex could not be inferred and was not explicitly M/F. Proceeding with 'inferred_sex' as NA.")
    }
  }
  if (verbose) message("Preprocessing of inputs complete.")

  return(
    list(
      somatic_snv = somatic_snv,
      germline_snv = germline_snv,
      het_pileups_wgs = het_pileups_wgs,
      cn.gr = cn.gr,
      dryclean.cov = dryclean.cov,
      inferred_sex = inferred_sex,
      purity = purity,
      ploidy = ploidy,
      mask = mask
    )
  )
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
  snpeff_path = system.file("extdata", "snpeff_scripts", package = "multiplicity"),
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
  bcftools = "bcftools" # FIXME: hardcoded for now.
  ) {
  if (debug)
    browser()
  tmp.path = tempfile(pattern = "tmp_", fileext = ".vcf.gz")
  if (!keepfile)
    on.exit(unlink(tmp.path))
  try2({
    # catcmd = if (grepl("(.gz)$", vcf)) "zcat" else "cat"
    catcmd = paste(bcftools, "view -Ov")
    ## Redundant logic is for backwards compatibility with workflows
    ## in which a snpeff path is already provided.
    ## Fallback option encoded here in case it doesn't.
    onepline_path1 = paste0(snpeff_path, "/vcfEffOnePerLine.pl")
    onepline_path2 = paste0(snpeff_path, "/scripts/vcfEffOnePerLine.pl")
    onepline_path3 = system.file("extdata", "snpeff_scripts", "vcfEffOnePerLine.pl", package = "multiplicity")
    onepline = onepline_path3
    is_path1_same_as_path3 = identical(onepline_path1, onepline_path3)
    if (!is_path1_same_as_path3 && file.exists(onepline_path1)) {
      onepline = onepline_path1
    } else if (!is_path1_same_as_path3 && file.exists(onepline_path2)) {
      onepline = onepline_path2
    }
    # onepline = paste0(snpeff_path, "/vcfEffOnePerLine.pl")
    snpsift_path1 = paste0(snpeff_path, "/SnpSift.jar")
    snpsift_path2 = paste0(snpeff_path, "/scripts/SnpSift.jar")
    snpsift_path3 = system.file("extdata", "snpeff_scripts", "SnpSift.jar", package = "multiplicity")
    is_path1_same_as_path3 = identical(snpsift_path1, snpsift_path3)
    snpsift = snpsift_path3
    if (!is_path1_same_as_path3 && file.exists(snpsift_path1)) {
      snpsift = snpsift_path1
    } else if (!is_path1_same_as_path3 && file.exists(snpsift_path2)) {
      snpsift = snpsift_path2
    }
    if (verbose)(message(paste0("applying SnpSift to VCF: ", vcf)))
    if (coding_alt_only) {
      if(verbose)(message("Coding alterations only."))
      filt = paste0("java -Xmx20m -Xms20m -XX:ParallelGCThreads=1 -jar ",
        # snpeff_path, "/SnpSift.jar ",
        snpsift, " ",
        "filter \"( ANN =~ 'chromosome_number_variation|exon_loss_variant|rare_amino_acid|stop_lost|transcript_ablation|coding_sequence|regulatory_region_ablation|TFBS|exon_loss|truncation|start_lost|missense|splice|stop_gained|frame' )\"")
      if (filterpass) {
        if(verbose) (message("Coding alterations only and FILTER == PASS variants only."))
        cmd = sprintf(paste(catcmd, "%s | %s | %s | %s view -i 'FILTER==\"PASS\"' | bgzip -c > %s"), vcf, onepline, filt, bcftools, tmp.path)  
      } else {
        if(verbose) (message("Coding alterations only."))
          cmd = sprintf("%s %s | %s | %s | bgzip -c > %s", catcmd, vcf, onepline, filt, tmp.path)
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
