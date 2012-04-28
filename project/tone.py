from __future__ import division
import numpy as np
import pyaudio

from pylab import psd, show, figure, subplot, colorbar

PA_TYPES = {
    np.int8: pyaudio.paInt8,
    np.uint8: pyaudio.paUInt8,
    np.int16: pyaudio.paInt16,
    np.int32: pyaudio.paInt32,
    np.float32: pyaudio.paFloat32
}

class PureTone(object):
    """Model a pure tone."""
    def __init__(self, hz, secs, fs=44100, nchans=1,
                 dtype=np.float32, output=1, precompute=False):
        self.hz = hz
        self.fs = fs
        self.secs = secs

        self.data = self.create() if precompute else None
        self.output = output
        
        assert nchans in (1, 2), 'only 1 or 2 channels allowed'
        self.nchans = nchans
        
        assert dtype in PA_TYPES, 'provided dtype not compatible with PyAudio types.'
        self.dtype = dtype

    def create(self):
        if self.data is None:
            secs, fs, hz, pi = self.secs, self.fs, self.hz, np.pi
            self.data = np.sin(
                np.arange(secs * fs) * hz * 2 * pi / fs
                ).astype(self.dtype)

    def play(self):
        with self as stream:
            stream.write(self.data.tostring())

    def __enter__(self):
        self.create()
        self.__p = pyaudio.PyAudio()
        self.__stream = self.__p.open(format=PA_TYPES[self.dtype],
                                      channels=self.nchans,
                                      rate=self.fs,
                                      output=self.output)
        return self.__stream

    def __exit__(self, tp, val, tb):
        self.__stream.close()
        self.__p.terminate()
        del self.__stream, self.__p

class NoisyTone(PureTone):
    """A tone with some additive noise."""
    def __init__(self, hz, secs, fs=44100, nchans=1, dtype=np.float32, output=1,
                 precompute=False, gen=np.random.randn):
        super(NoisyTone, self).__init__(hz, secs, fs=44100, nchans=1,
                                        dtype=np.float32, output=1, precompute=False)
        self.gen = gen
    
    def create(self):
        super(NoisyTone, self).create()
        self.data += self.gen(self.data.size).astype(self.dtype)

def plot_tones(*tones):
    n = len(tones)
    figure()
    for i, tone in enumerate(tones):
        subplot(n, 2, 2 * i + 1)
        specgram(tone.data, Fs=tone.fs)
        colorbar()
        xlabel('Time')
        ylabel('Frequency')
        axis('tight')
        
        subplot(n, 2, 2 * (i + 1))
        plot(tone.data[:500])
        axis('tight')
    show()

def play_tones(*tones):
    for tone in tones:
        tone.play()

def play_and_plot_tones(*tones):
    play_tones(*tones)
    plot_tones(*tones)
    
    
if __name__ == '__main__':
    a = 440
    d = 0.5
    pure = PureTone(a, d)
    noisy_unif = NoisyTone(a, d, gen=np.random.rand)
    noisy_norm = NoisyTone(a, d, gen=np.random.randn)
    noisy_beta = NoisyTone(a, d, gen=lambda x: np.random.beta(1, 2, size=(x,)))
    noisy_exp = NoisyTone(a, d, gen=lambda x: np.random.exponential(abs(randn()), size=(x,)))
    noisy_pois = NoisyTone(a, d, gen=lambda x: np.random.poisson(abs(5 * randn()), size=(x,)))
    noisy_tri = NoisyTone(a, d, gen=lambda x: np.random.triangular(size=(x,)))

    close('all')
    play_and_plot_tones(noisy_exp, noisy_pois)
    # play_and_plot_tones(pure, noisy_unif, noisy_norm, noisy_beta, noisy_exp, noisy_pois)
