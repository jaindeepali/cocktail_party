library('fastICA')
library('tuneR')
library('Metrics')

## Input Sources

# In case of using mixed sources:
# S1 <- readWave('sound_files/mixed_sources/2/1.wav')
# S2 <- readWave('sound_files/mixed_sources/2/2.wav')
# S1 <- normalize(S1, unit = '8')
# S2 <- normalize(S2, unit = '8')
# X <- cbind(S1@left, S2@left)

# In case of using original sources:
# S1 <- readWave('sound_files/original_sources/3/1.wav')
# S2 <- readWave('sound_files/original_sources/3/2.wav')

S1 <- sine(2)
S2 <- sine(9)

# S1 <- 3*sine(2) + 4*sine(10)
# S2 <- 5*sine(3) + 7*sine(9)

S1 <- normalize(S1, unit = '8')
S2 <- normalize(S2, unit = '8')
S <- cbind(S1@left, S2@left)
A <- matrix(runif(4, 1, 5), 2, 2)
X <- S %*% A

## Fast ICA

startTime <- Sys.time()

a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)

endTime <- Sys.time()
print('fastICA running time:')
print(endTime - startTime)


## Output Signals

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


plot(1:nrow(X), X[,1 ], type = "l", main = "Mixed Signals", xlab = "", ylab = "")
plot(1:nrow(X), X[,2 ], type = "l", xlab = "", ylab = "")

plot(1:length(out1@left), out1@left, type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:length(out2@left), out2@left, type = "l", xlab = "", ylab = "")

## Plot Signals Distribution

par(mfcol = c(1, 1))

# In case of using original sources:
# plot(S[,1], S[,2], main = "Original Signals Distribution", xlab = "S1", ylab = "S2", cex = 0.1)

plot(X[,1], X[,2], main = "Mixed Signals Distribution", xlab = "X1", ylab = "X2", cex = 0.1)

## Save output signals

writeWave(out1, filename = 'sound_files/output/singleSine_1.wav')
writeWave(out2, filename = 'sound_files/output/singleSine_2.wav')

## Calculate error in case of using original sources:

S1 <- normalize(S1, unit = '1')
S2 <- normalize(S2, unit = '1')
out1 <- normalize(out1, unit = '1')
out2 <- normalize(out2, unit = '1')

err1 <- min(rmse(abs(S1@left), abs(out1@left)), rmse(abs(S2@left), abs(out1@left)))
err2 <- min(rmse(abs(S1@left), abs(out2@left)), rmse(abs(S2@left), abs(out2@left)))
err <- mean(c(err1,err2))
print(err)