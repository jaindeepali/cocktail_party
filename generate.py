import numpy as np
import scipy.io.wavfile as wavfile

f = "sound_files/original_sources/linear_mono_human/2.wav"

fs, s = wavfile.read(f)

s = np.concatenate([s, s[:50000-44100]])
print np.shape(s)

wavfile.write(f, 8000, s)

