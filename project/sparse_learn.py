from pylab import *


class Bunch(dict):
    """Class from sklearn for use of a dict as a MATLAB-style
    struct.
    """
    def __init__(self, **kwargs):
        super(Bunch, self).__init__(kwargs)
        self.__dict__ = self


def wavread(f):
    from struct import unpack_from
    import wave

    wav = wave.open(f, 'r')
    nchan, _, fs, nsamps, _, _ = wav.getparams()
    frames = wav.readframes(nsamps * nchan)
    wav.close()
    out = asarray(unpack_from('{}h'.format(nsamps * nchan), frames))
    return (out[::2], out[1::2], fs) if nchan == 2 else out, fs


def load_all_wav(d='sounds'):
    import glob
    import os
    import operator
    import itertools
    import re
    import tokenize
    import keyword
    
    fullnames = glob.glob('{}/*.wav'.format(d))
    
    if not fullnames:
        raise ValueError('no files found')

    n = len(fullnames)
    splitextfiles = itertools.imap(os.path.splitext, fullnames)
    sansextfiles = itertools.imap(operator.getitem, splitextfiles, [0] * n)
    filenames = map(os.path.basename, sansextfiles)
    data = Bunch(**dict(itertools.izip(filenames, [None] * n)))
    for filename, fullname in itertools.izip(filenames, fullnames):
        if not (re.match(tokenize.Name, filename) or
                keyword.iskeyword(filename)):
            raise SyntaxError(
                'fields must be valid, non-reserved Python identifier.'
                )
        d, fs = wavread(fullname)
        data[filename] = Bunch(data=d, fs=fs)
    return data

    
if __name__ == '__main__':
    close('all')
    datasets = load_all_wav()

    fig1, fig2 = figure(), figure()
    
    for i, dataset in enumerate(datasets):
        dsname = ' '.join(dataset.split('_'))
        fig1.add_subplot(4, 3, i + 1)
        ax1 = fig1.axes[i]
        p, freqs, bins, _ = ax1.specgram(datasets[dataset].data,
                                         Fs=datasets[dataset].fs)
        ax1.axis('tight')
        ax1.set_title(dsname)
        ax1.axhline(y=2e4, color='k', linestyle='--', linewidth=3)
        
        fig2.colorbar()
        fig2.add_subplot(4, 3, i + 1)
        ax2 = fig2.axes[i]
        ax2.imshow(p)
        ax2.axis('tight')
        fig2.colorbar()
        ax2.set_title(dsname)
