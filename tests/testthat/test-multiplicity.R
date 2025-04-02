library(testthat)
library(multiplicity)

context("multiplicity function tests")

test_that("multiplicity fails when no SNV inputs are provided", {
  expect_error(
    multiplicity(),
    "Somatic VCF, Germline VCF, and/or HetSNPs file must be provided."
  )
})

# Get file paths from the inst/testdata folder
cov_file    <- system.file("testdata", "sample_cov.rds", package = "multiplicity")
cbs_file    <- system.file("testdata", "sample_cbs.rds", package = "multiplicity")
snpeff_vcf_file  <- system.file("testdata", "sample_chr22_revised.vcf", package = "multiplicity")
jabba_file  <- system.file("testdata", "sample_jabba.rds", package = "multiplicity")
het_pileups_wgs <- system.file("testdata", "sample_het_pileups.txt", package = "multiplicity")

# Test that the test data files exist
test_that("Test data files exist", {
  expect_true(file.exists(cov_file))
  expect_true(file.exists(cbs_file))
  expect_true(file.exists(snpeff_vcf_file))
  expect_true(file.exists(jabba_file))
  expect_true(file.exists(het_pileups_wgs))
})

test_that("multiplicity test of somatic variants transformation with CBS rescaling; recommended solution for running", {
  result1 <- suppressWarnings(multiplicity(
    somatic_snv     = snpeff_vcf_file,
    germline_snv    = "/dev/null",
    het_pileups_wgs = "/dev/null",
    tumor_cbs       = cbs_file,
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  
  result1[[1]] -> list1
  expect_true(is.list(result1))
  expect_equal(length(result1), 3)
  expect_equal(length(result1[[2]]), 0)
  expect_equal(length(result1[[3]]), 0)

  list1  %Q% (seqnames == 22 & start == 37492267) -> variant
  expect_equal(round(variant$altered_copies, 3), 0.551)
  expect_equal(round(variant$total_snv_copies, 3), 1.848)
})

test_that("multiplicity test of somatic variants transformation with DRYCLEAN rescaling", {
  result2 <- suppressWarnings(multiplicity(
    somatic_snv     = snpeff_vcf_file,
    germline_snv    = "/dev/null",
    het_pileups_wgs = "/dev/null",
    dryclean_field  = "foreground",
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  result2[[1]] -> list2
  expect_true(is.list(result2))
  expect_equal(length(result2), 3)
  expect_equal(length(result2[[2]]), 0)
  expect_equal(length(result2[[3]]), 0)
  
  list2  %Q% (seqnames == 22 & start == 37492267) -> variant
  expect_equal(round(variant$altered_copies, 3), 0.757)
  expect_equal(round(variant$total_snv_copies, 3), 2.538)
})

test_that("multiplicity test of somatic variants transformation with NO rescaling", {
  result3 <- suppressWarnings(multiplicity(
    somatic_snv     = snpeff_vcf_file,
    germline_snv    = "/dev/null",
    het_pileups_wgs = "/dev/null",
    tumor_dryclean  = cov_file,
    dryclean_field  = "foreground",
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  result3[[1]] -> list3
  expect_true(is.list(result3))
  expect_equal(length(result3), 3)
  expect_equal(length(result3[[2]]), 0)
  expect_equal(length(result3[[3]]), 0)

  list3 %Q% (seqnames == 22 & start == 37492267) -> variant
  expect_equal(round(variant$altered_copies, 3), 0.549)
  expect_equal(round(variant$total_snv_copies, 3), 1.84)
})

round(variant$altered_copies, 3)
round(variant$total_snv_copies, 3)

test_that("multiplicity test of hetsnps transformation with CBS rescaling; recommended solution for running", {
  result4 <- suppressWarnings(multiplicity(
    somatic_snv     = "/dev/null",
    germline_snv    = "/dev/null",
    het_pileups_wgs = het_pileups_wgs,
    tumor_cbs       = cbs_file,
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  result4[[3]] -> list4
  
  expect_true(is.list(result4))
  expect_equal(length(result4), 3)
  expect_equal(length(result4[[1]]), 0)
  expect_equal(length(result4[[2]]), 0)

  list4[10000]-> variant
  expect_equal(round(variant$altered_copies, 3), 1.365)
  expect_equal(round(variant$total_snv_copies, 3), 2.623)
})

test_that("multiplicity test of hetsnps transformation with DRYCLEAN rescaling", {
  result5 <- suppressWarnings(multiplicity(
    somatic_snv     = "/dev/null",
    germline_snv    = "/dev/null",
    het_pileups_wgs = het_pileups_wgs,
    dryclean_field  = "foreground",
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  result5[[3]] -> list5
  expect_true(is.list(result5))
  expect_equal(length(result5), 3)
  expect_equal(length(result5[[1]]), 0)
  expect_equal(length(result5[[2]]), 0)

  list5[10000] -> variant
  expect_equal(round(variant$altered_copies, 3), 1.612)
  expect_equal(round(variant$total_snv_copies, 3), 3.097)
})

test_that("multiplicity test of hetsnps transformation with NO rescaling", {
  result6 <- suppressWarnings(multiplicity(
    somatic_snv     = "/dev/null",
    germline_snv    = "/dev/null",
    het_pileups_wgs = het_pileups_wgs,
    tumor_dryclean  = cov_file,
    dryclean_field  = "foreground",
    jabba_rds       = jabba_file,
    snpeff_path     = "~/modules/SnpEff/source/snpEff",
    tumor_name      = "TUMOR",
    normal_name     = "NORMAL",
    filterpass      = TRUE,
    tau_in_gamma    = FALSE,
    read_size       = 151,
    verbose         = TRUE
  ))
  result6[[3]] -> list6
  expect_true(is.list(result6))
  expect_equal(length(result6), 3)
  expect_equal(length(result6[[1]]), 0)
  expect_equal(length(result6[[2]]), 0)

  list6[10000] -> variant
  expect_equal(round(variant$altered_copies, 3), 3.451)
  expect_equal(round(variant$total_snv_copies, 3), 6.629)
})

