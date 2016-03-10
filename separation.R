library('fastICA')
library('tuneR')
library('Metrics')
library('signal')
library('ppls')

generateLinearMixture <- function(osList) {
    print("Generating linear mixture")
    A <- matrix(runif(4, 1, 5), 2, 2)
    return(osList %*% A)
}

separateUsingICA <- function(sample, use_original=TRUE, sample_original=TRUE) {
    originalSignalDir <- paste('sound_files/original_sources', sample, sep='/')
    mixedSignalDir <- paste('sound_files/mixed_sources', sample, sep='/')
    outputDir <- paste('sound_files/ica_output', sample, sep='/')

    if(sample_original) {
        os1 <- readWave(paste(originalSignalDir, '1.wav', sep='/'))
        os2 <- readWave(paste(originalSignalDir, '2.wav', sep='/'))

        osList <- cbind(os1@left, os2@left)

        sample_rate <- os1@samp.rate
        bitrate <- os1@bit
    }

    if(use_original) {
        msList <- generateLinearMixture(osList)
        ms1 <- msList[,1]
        ms2 <- msList[,2]
    } else {
        ms1 <- readWave(paste(mixedSignalDir, '1.wav', sep='/'))
        ms2 <- readWave(paste(mixedSignalDir, '2.wav', sep='/'))
        msList <- cbind(ms1@left, ms2@left)

        if (! sample_original) {
            sample_rate <- ms1@samp.rate
            bitrate <- ms1@bit
        }
    }
    startTime <- Sys.time()

    a <- fastICA(msList, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
                 method = "R", row.norm = FALSE, maxit = 200,
                 tol = 0.0001, verbose = TRUE)

    endTime <- Sys.time()
    print('fastICA running time:')
    print(endTime - startTime)

    out1 <- Wave(left = a$S[,1], samp.rate = sample_rate, bit = bitrate)
    out2 <- Wave(left = a$S[,2], samp.rate = sample_rate, bit = bitrate)

    par(mfcol = c(2, 3))

    if(sample_original) {
        plot(1:nrow(osList), osList[,1], type = "l", main = "Original Signals",
             xlab = "", ylab = "")
        plot(1:nrow(osList), osList[,2], type = "l", xlab = "", ylab = "")
    }

    plot(1:nrow(msList), msList[,1 ], type = "l", main = "Mixed Signals", xlab = "", ylab = "")
    plot(1:nrow(msList), msList[,2 ], type = "l", xlab = "", ylab = "")

    plot(1:length(out1@left), out1@left, type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
    plot(1:length(out2@left), out2@left, type = "l", xlab = "", ylab = "")

    writeWave(normalize(out1, unit = '8'), filename = paste(outputDir, '1.wav', sep='/'))
    writeWave(normalize(out2, unit = '8'), filename = paste(outputDir, '2.wav', sep='/'))

    if(sample_original) {
        os1 <- os1@left / sqrt(sum(os1@left ^ 2))
        os2 <- os2@left / sqrt(sum(os2@left ^ 2))
        out1 <- out1@left / sqrt(sum(out1@left ^ 2))
        out2 <- out2@left / sqrt(sum(out2@left ^ 2))

        e1 <- abs((abs(os1) - abs(out1)) / abs(os1))
        e1 <- e1[is.finite(e1)]
        e1 <- sqrt(mean(e1 ^ 2))
        print(e1)
        e2 <- abs((abs(os2) - abs(out2)) / abs(os2))
        e2 <- e2[is.finite(e2)]
        e2 <- sqrt(mean(e2 ^ 2))
        print(e2)
        ef <- abs(mean(c(e1, e2)))

        e1 <- abs((abs(os1) - abs(out2)) / abs(os1))
        e1 <- e1[is.finite(e1)]
        e1 <- sqrt(mean(e1 ^ 2))
        print(e1)
        e2 <- abs((abs(os2) - abs(out1)) / abs(os2))
        e2 <- e2[is.finite(e2)]
        e2 <- sqrt(mean(e2 ^ 2))
        print(e2)
        es <- abs(mean(c(e1, e2)))

        err <- min(c(ef, es))
        print('ERROR:')
        print(err)
    }
}

## Plot Signals Distribution

# par(mfcol = c(1, 1))

# In case of using original sources:
# plot(S[,1], S[,2], main = "Original Signals Distribution", xlab = "S1", ylab = "S2", cex = 0.1)

# plot(X[,1], X[,2], main = "Mixed Signals Distribution", xlab = "X1", ylab = "X2", cex = 0.1)

## Calculate error in case of using original sources:

separateUsingICA('2', use_original=FALSE)
