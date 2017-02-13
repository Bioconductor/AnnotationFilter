context("expression translation")

test_that("translation of expression works for single filter/condition", {
    ## Check for some character filter.
    ## exon_id
    flt <- ExonIdFilter("EX1", condition = "==")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(exon_id == "EX1")
    expect_equal(flt, flt2)
    flt <- ExonIdFilter(c("EX1", "EX2"), condition = "!=")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(exon_id != c("EX1", "EX2"))
    expect_equal(flt, flt2)
    ## seq_name
    flt <- SeqNameFilter(c("chr3", "chrX"), condition = "==")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(seq_name == c("chr3",
                                                                           "chrX"))
    expect_equal(flt, flt2)
    flt <- SeqNameFilter(1:3, condition = "==")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(seq_name %in% 1:3)
    expect_equal(flt, flt2)
    ## Check IntegerFilter
    flt <- GeneStartFilter(123, condition = ">")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(gene_start > 123)
    expect_equal(flt, flt2)
    flt <- TxStartFilter(123, condition = "<")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(tx_start < 123)
    expect_equal(flt, flt2)
    flt <- GeneEndFilter(123, condition = ">=")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(gene_end >= 123)
    expect_equal(flt, flt2)
    flt <- ExonEndFilter(123, condition = "<=")
    flt2 <- AnnotationFilter:::.expressionToAnnotationFilter(exon_end <= 123)
    expect_equal(flt, flt2)
    ## Test exceptions/errors.
    expect_error(AnnotationFilter:::.expressionToAnnotationFilter(not_existing == 1:3))
    ## Throws an error, but is not self-explanatory.
    expect_error(AnnotationFilter:::.expressionToAnnotationFilter(gene_id * 3))
})

