symdiff <- function(x, y) {
    dif <- NULL
    for (i in 1:dim(x)[1]) {
        dif[i] <- sum(abs(x[i, ] - y))
    }
    return(dif)
}
