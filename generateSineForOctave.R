library('fastICA')
library('tuneR')

S1 <- sine(50)
S2 <- sine(100)

writeWave(S1, filename = 'sound_files/original_sources/singleSine/1.wav')
writeWave(S2, filename = 'sound_files/original_sources/singleSine/2.wav')
