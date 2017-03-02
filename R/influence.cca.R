### influence statistics for cca objects


## hat values need adjustment, because QR ignores Intercept

`hatvalues.cca` <-
    function(model, ...)
 {
     rowSums(qr.Q(model$CCA$QR)^2) + 1/nrow(model$CCA$QR$qr)
 }
