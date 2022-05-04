SL.mean.wrong <- function (Y, X, newX, family, obsWeights, id, ...)
{
    # temp <- weighted.mean(y - offset, weights) + 0.05
    # if (temp > 0.9) temp <- 0.9

    meanY <- weighted.mean(Y, w = obsWeights) + 0.05
    if (meanY > 0.99) meanY <- 0.99
    if (meanY < 0.01) meanY <- 0.01

    pred <- rep.int(meanY, times = nrow(newX))
    fit <- list(object = meanY)
    out <- list(pred = pred, fit = fit)
    class(out$fit) <- c("SL.mean")
    return(out)
}

predict.SL.mean.wrong <- function (object, newdata, family, X = NULL, Y = NULL, ...)
{
    pred <- rep.int(object$object, times = nrow(newdata))
    return(pred)
}
