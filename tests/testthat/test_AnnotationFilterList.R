context("AnnotationFilterList")

test_that("AnnotationFilterList() works", {
    f1 <- GeneIdFilter("somegene")
    f2 <- SeqNameFilter("chr3")
    f3 <- GeneBiotypeFilter("protein_coding", "!=")

    fL <- AnnotationFilter:::AnnotationFilterList(f1, f2)
    expect_true(length(fL) == 2)
    expect_equal(fL[[1]], f1)
    expect_equal(fL[[2]], f2)
    expect_true(all(fL@logOp == "&"))
    
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
})
