from struct import unpack_from
import wave
from pylab import *
from sklearn.decomposition import dict_learning_online


def wavread(f):
    wav = wave.open(f, 'r')
    nchan, _, fs, nsamps, _, _ = wav.getparams()
    frames = wav.readframes(nsamps * nchan)
    wav.close()
    out = asarray(unpack_from('{}h'.format(nsamps * nchan), frames))
    return (out[::2], out[1::2], fs) if nchan == 2 else out, fs


close('all')
o, fs = wavread('/home/cpcloud/bird-stereo.wav')
ax = subplot(2, 1, 1)
p, freqs, bins, _ = specgram(o, Fs=fs)
axis('tight')
colorbar(orientation='horizontal', spacing='proportional')
title('Time - Frequency Plot')
ylabel('Frequency (Hz)')

subplot(2, 1, 2, sharex=ax)
plot(linspace(1, o.size, o.size) / fs, o)
axis('tight')
xlabel('Time (s)')
ylabel('Amplitude')

figure()
# subplot(2, 1, 1)
imshow(p)
axis('tight')
# subplot(2, 1, 2)
# hist(p)
show()
