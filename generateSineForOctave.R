library('fastICA')
library('tuneR')

S1 <- 3*sine(59)
S2 <- 5*sine(231)

S1 <- normalize(S1, unit = '8')
S2 <- normalize(S2, unit = '8')

writeWave(S1, filename = 'sound_files/original_sources/singleSine/1.wav')
writeWave(S2, filename = 'sound_files/original_sources/singleSine/2.wav')
