box::use(testthat,
         tra=./tra)


testthat$context("Test binarization of expression across tissues")

expression <- c(rep(1, 5), 3, 3, 4, rep(1, 2))
names(expression) <- paste0("Tissue", 1:length(expression))
binarized <- c(rep(0, 5), rep(1, 3), rep(0,2))
names(binarized) <- names(expression)

testthat$expect_equal(binarized, tra$find_tra_tissues(expression))


