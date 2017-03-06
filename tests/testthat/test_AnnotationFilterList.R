context("AnnotationFilterList")

test_that("AnnotationFilterList() works", {
    logOp <- AnnotationFilter:::.logOp
    f1 <- GeneIdFilter("somegene")
    f2 <- SeqNameFilter("chr3")
    f3 <- GeneBiotypeFilter("protein_coding", "!=")

    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2)
    expect_true(length(fL) == 2)
    expect_equal(fL[[1]], f1)
    expect_equal(fL[[2]], f2)
    expect_true(all(logOp(fL) == "&"))
    
    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2, f3,
                                                  logOp = c("&", "|"))
    expect_true(length(fL) == 3)
    expect_equal(fL[[1]], f1)
    expect_equal(fL[[2]], f2)
    expect_equal(fL[[3]], f3)
    expect_equal(fL@logOp, c("&", "|"))

    ## A AnnotationFilterList with and AnnotationFilterList
    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2, logOp = "|")
    fL2 <- AnnotationFilter:::AnnotationFilterList(f3, fL, logOp = "&")
    expect_true(length(fL) == 2)
    expect_true(length(fL2) == 2)
    expect_true(is(value(fL2)[[1]], "GeneBiotypeFilter"))
    expect_true(is(value(fL2)[[2]], "AnnotationFilterList"))
    expect_equal(value(fL2)[[2]], fL)
    expect_equal(fL2[[2]], fL)
    expect_equal(logOp(fL2), c("&"))
    expect_equal(logOp(fL2[[2]]), c("|"))
})
