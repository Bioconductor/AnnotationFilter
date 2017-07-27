context("expression translation")

test_that("translation of expression works for single filter/condition", {
    ## Check for some character filter.
    ## exon_id
    flt <- ExonIdFilter("EX1", condition = "==")
    flt2 <- AnnotationFilter(~ exon_id == "EX1")
    expect_equal(flt, flt2)
    flt <- ExonIdFilter(c("EX1", "EX2"), condition = "!=")
    flt2 <- AnnotationFilter(~ exon_id != c("EX1", "EX2"))
    expect_equal(flt, flt2)
    ## seq_name
    flt <- SeqNameFilter(c("chr3", "chrX"), condition = "==")
    flt2 <- AnnotationFilter(~ seq_name == c("chr3", "chrX"))
    expect_equal(flt, flt2)
    flt <- SeqNameFilter(1:3, condition = "==")
    flt2 <- AnnotationFilter(~ seq_name %in% 1:3)
    expect_equal(flt, flt2)
    ## Check IntegerFilter
    flt <- GeneStartFilter(123, condition = ">")
    flt2 <- AnnotationFilter(~ gene_start > 123)
    expect_equal(flt, flt2)
    flt <- TxStartFilter(123, condition = "<")
    flt2 <- AnnotationFilter(~ tx_start < 123)
    expect_equal(flt, flt2)
    flt <- GeneEndFilter(123, condition = ">=")
    flt2 <- AnnotationFilter(~ gene_end >= 123)
    expect_equal(flt, flt2)
    flt <- ExonEndFilter(123, condition = "<=")
    flt2 <- AnnotationFilter(~ exon_end <= 123)
    expect_equal(flt, flt2)
    ## Test exceptions/errors.
    expect_error(AnnotationFilter(~ not_existing == 1:3))
    ## Throws an error, but is not self-explanatory.
    expect_error(AnnotationFilter(~ gene_id * 3))
})

test_that("translation of combined expressions works", {
    res <- AnnotationFilter(~ exon_id == "EX1" & genename == "BCL2")
    cmp <- AnnotationFilterList(ExonIdFilter("EX1"), GenenameFilter("BCL2"))
    expect_equal(res, cmp)
    res <- AnnotationFilter(~ exon_id == "EX1" | genename != "BCL2")
    cmp <- AnnotationFilterList(ExonIdFilter("EX1"),
                                GenenameFilter("BCL2", "!="), logicOp = "|")
    expect_equal(res, cmp)
    ## 3 filters.
    res <- AnnotationFilter(~ exon_id == "EX1" & genename == "BCL2" |
                                seq_name != 3)
    ## Expect an AnnotationFilterList of length 3.
    expect_equal(length(res), 3)
    cmp <- AnnotationFilterList(ExonIdFilter("EX1"), GenenameFilter("BCL2"),
                                SeqNameFilter(3, "!="), logicOp = c("&", "|"))
    expect_equal(res, cmp)
    ## 4 filters.
    res <- AnnotationFilter(~ exon_id == "EX1" & genename == "BCL2" |
                                seq_name != 3 | seq_name == "Y")
    expect_equal(length(res), 4)
    cmp <- AnnotationFilterList(ExonIdFilter("EX1"), GenenameFilter("BCL2"),
                                SeqNameFilter(3, "!="), SeqNameFilter("Y"),
                                logicOp = c("&", "|",  "|"))
    expect_equal(res, cmp)
})

test_that("translation works from within other functions", {
    simpleFun <- function(x)
        AnnotationFilter(x)
    expect_equal(simpleFun(~ gene_id == 4), AnnotationFilter(~ gene_id == 4))
    filter_expr <- ~ gene_id == 4
    expect_equal(simpleFun(filter_expr),
                 AnnotationFilter(~ gene_id == 4))
})

## This might be a test if we get the nesting working.
## test_that("translation of nested expressions works" {
##     res <- convertFilterExpression((exon_id == "EX1" & gene_id == "BCL2") |
##                                    (exon_id == "EX3" & gene_id == "BCL2L11"))
##     expect_equal(logicOp(res), "|")
##     expect_true(is(res[[1]], "AnnotationFilterList"))
##     expect_equal(res[[1]][[1]], ExonIdFilter("EX1"))
##     expect_equal(res[[1]][[2]], GeneIdFilter("BCL2"))
##     expect_equal(logicOp(res[[1]]), "&")
##     expect_true(is(res[[2]], "AnnotationFilterList"))
##     expect_equal(res[[2]][[1]], ExonIdFilter("EX3"))
##     expect_equal(res[[2]][[2]], GeneIdFilter("BCL2L11"))
##     expect_equal(logicOp(res[[2]]), "&")
##     ##
##     res <- convertFilterExpression(seq_name == "Y" |
##                                    (exon_id == "EX1" & gene_id == "BCL2") &
##                                    (exon_id == "EX3" & gene_id == "BCL2L11"))
##     ## Expect: length 3, first being a SeqNameFilter, second an
##     ## AnnotationFilterList, third a AnnotationFilterList.
##     expect_equal(res[[1]], SeqNameFilter("Y"))
##     expect_equal(logicOp(res), "|")
##     expect_true(is(res[[2]], "AnnotationFilterList"))
##     expect_equal(res[[1]][[1]], ExonIdFilter("EX1"))
##     expect_equal(res[[1]][[2]], GeneIdFilter("BCL2"))
##     expect_equal(logicOp(res[[1]]), "&")
##     expect_true(is(res[[2]], "AnnotationFilterList"))
##     expect_equal(res[[2]][[1]], ExonIdFilter("EX3"))
##     expect_equal(res[[2]][[2]], GeneIdFilter("BCL2L11"))
##     expect_equal(logicOp(res[[2]]), "&")

##     expect_true(is(res[[1]], "AnnotationFilterList"))
##     expect_true(is(res[[2]], "AnnotationFilterList"))

##     convertFilterExpression((gene_id == 3) ()
## })

