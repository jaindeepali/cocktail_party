library('fastICA')
library('tuneR')

## Input Sources

# In case of using mixed sources:
# S1 <- readWave('sound files/mixed sources/mix2a.wav')
# S2 <- readWave('sound files/mixed sources/mix2b.wav')
# S1 <- normalize(S1, unit = '8')
# S2 <- normalize(S2, unit = '8')
# X <- cbind(S1@left, S2@left)

# In case of using original sources:
# S1 <- readWave('sound files/original sources/source1a.wav')
# S2 <- readWave('sound files/original sources/source1b.wav')

S1 <- sine(2)
S2 <- sine(10)

S1 <- normalize(S1, unit = '8')
S2 <- normalize(S2, unit = '8')
S <- cbind(S1@left, S2@left)
A <- matrix(runif(4, 1, 5), 2, 2)
X <- S %*% A

## Fast ICA

a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)

## Plot Signals

out1 <- Wave(left = a$S[,1], samp.rate = S1@samp.rate, bit = S1@bit)
out2 <- Wave(left = a$S[,2], samp.rate = S1@samp.rate, bit = S1@bit)
out1 <- normalize(out1, unit = '8')
out2 <- normalize(out2, unit = '8')

## Plot Signals

par(mfcol = c(2, 3))

# In case of using original sources:
plot(1:nrow(S), S[,1 ], type = "l", main = "Original Signals",
     xlab = "", ylab = "")
plot(1:nrow(S), S[,2 ], type = "l", xlab = "", ylab = "")


plot(1:nrow(X), X[,1 ], type = "l", main = "Mixed Signals",
     xlab = "", ylab = "")
plot(1:nrow(X), X[,2 ], type = "l", xlab = "", ylab = "")

plot(1:length(out1@left), out1@left, type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:length(out2@left), out2@left, type = "l", xlab = "", ylab = "")

## Plot Signals Distribution

par(mfcol = c(1, 1))

# In case of using original sources:
plot(S[,1], S[,2], main = "Original Signals Distribution", xlab = "S1", ylab = "S2", cex = 0.1)

plot(X[,1], X[,2], main = "Mixed Signals Distribution", xlab = "X1", ylab = "X2", cex = 0.1)

## Save output signals

writeWave(out1, filename = 'sound files/output/outSinea.wav')
writeWave(out2, filename = 'sound files/output/outSineb.wav')

## Calculate error in case of using original sources:

S1 <- normalize(S1, unit = '1')
S2 <- normalize(S2, unit = '1')
out1 <- normalize(out1, unit = '1')
out2 <- normalize(out2, unit = '1')

err1 <- min(mean((out1**2 - S1**2)@left), mean((out1**2 - S2**2)@left))
print(err1)
err2 <- min(mean((out2**2 - S1**2)@left), mean((out2**2 - S2**2)@left))
print(err2)