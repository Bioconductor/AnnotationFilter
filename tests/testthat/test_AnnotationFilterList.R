context("AnnotationFilterList")

test_that("AnnotationFilterList() works", {
    f1 <- GeneIdFilter("somegene")
    f2 <- SeqNameFilter("chr3")
    f3 <- GeneBiotypeFilter("protein_coding", "!=")

    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2)
    expect_true(length(fL) == 2)
    expect_equal(fL[[1]], f1)
    expect_equal(fL[[2]], f2)
    expect_true(all(logicOp(fL) == "&"))
    
    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2, f3,
                                                  logicOp = c("&", "|"))
    expect_true(length(fL) == 3)
    expect_equal(fL[[1]], f1)
    expect_equal(fL[[2]], f2)
    expect_equal(fL[[3]], f3)
    expect_equal(logicOp(fL), c("&", "|"))

    ## A AnnotationFilterList with and AnnotationFilterList
    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2, logicOp = "|")
    fL2 <- AnnotationFilter:::AnnotationFilterList(f3, fL, logicOp = "&")
    expect_true(length(fL) == 2)
    expect_true(length(fL2) == 2)
    expect_true(is(value(fL2)[[1]], "GeneBiotypeFilter"))
    expect_true(is(value(fL2)[[2]], "AnnotationFilterList"))
    expect_equal(value(fL2)[[2]], fL)
    expect_equal(fL2[[2]], fL)
    expect_equal(logicOp(fL2), "&")
    expect_equal(logicOp(fL2[[2]]), "|")
})

test_that("empty elements in AnnotationFilterList", {
    ## empty elements should be removed from the AnnotationFilterList.
    empty_afl <- AnnotationFilterList()
    afl <- AnnotationFilterList(empty_afl)
    expect_true(length(afl) == 0)
    afl <- AnnotationFilterList(GeneIdFilter(4), empty_afl)
    expect_true(length(afl) == 1)
    afl <- AnnotationFilterList(GeneIdFilter(4),
        AnnotationFilter(~ gene_id == 3 | seq_name == 4),empty_afl)
    expect_true(length(afl) == 2)
    ## Check validate.
    afl@.Data <- c(afl@.Data, list(empty_afl))
    ## Fix also the logOp.
    afl@logOp <- c(afl@logOp, "|")
    expect_error(validObject(afl))
})

test_that("convertFilter works", {
    smbl <- SymbolFilter("ADA")
    txid <- TxIdFilter(1000)
    gr <- GRangesFilter(GenomicRanges::GRanges("chr15:25062333-25065121"))

    expect_identical(convertFilter(AnnotationFilter(~smbl | txid)),
        "symbol == 'ADA' | tx_id == '1000'")
    expect_identical(convertFilter(AnnotationFilter(~smbl & (smbl | txid))),
        "symbol == 'ADA' & (symbol == 'ADA' | tx_id == '1000')")
    expect_identical(convertFilter(AnnotationFilter(~smbl & !(smbl | txid))),
        "symbol == 'ADA' & !(symbol == 'ADA' | tx_id == '1000')")
    expect_error(convertFilter(AnnotationFilter(smbl | (txid & gr))))
    
})

test_that("distributeNegation works", {
    afl <- AnnotationFilter(~!(symbol == 'ADA' | symbol %startsWith% 'SNORD'))
    afl2 <- AnnotationFilter(~!symbol == 'ADA' & !symbol %startsWith% 'SNORD')
    expect_identical(distributeNegation(afl), afl2)
})
