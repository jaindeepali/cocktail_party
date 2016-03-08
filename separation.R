library('fastICA')
library('tuneR')

## Input Sources

# In case of using mixed sources:

# S1 <- readWave('sound files/mixed sources/X1_linear.wav')
# S2 <- readWave('sound files/mixed sources/X2_linear.wav')
# X <- cbind(S1@left, S2@left)

# In case of using original sources:

S1 <- readWave('sound files/original sources/source1.wav')
S2 <- readWave('sound files/original sources/source2.wav')
S <- cbind(S1@left, S2@left)
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A

## Fast ICA

a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)

## Plot Signals

par(mfcol = c(2, 3))

# In case of using original sources:

plot(1:nrow(S), S[,1 ], type = "l", main = "Original Signals",
     xlab = "", ylab = "")
plot(1:nrow(S), S[,2 ], type = "l", xlab = "", ylab = "")


plot(1:nrow(X), X[,1 ], type = "l", main = "Mixed Signals",
     xlab = "", ylab = "")
plot(1:nrow(X), X[,2 ], type = "l", xlab = "", ylab = "")

plot(1:nrow(a$S), a$S[,1 ], type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:nrow(a$S), a$S[,2], type = "l", xlab = "", ylab = "")

## Save output signals

out1 <- Wave(left = a$S[,1], samp.rate = S1@samp.rate, bit = S1@bit)
out2 <- Wave(left = a$S[,2], samp.rate = S1@samp.rate, bit = S1@bit)
out1 <- normalize(out1, unit = '8')
out2 <- normalize(out2, unit = '8')

writeWave(out1, filename = 'sound files/output/out1.wav')
writeWave(out2, filename = 'sound files/output/out2.wav')
