import ctypes
import argparse
import six
import json
from collections import Counter


lib = ctypes.CDLL("target/release/libgundam.dylib")

lib.read_kmers.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
lib.read_kmers.restype = ctypes.c_void_p

lib.get_dyad.argtypes = [ctypes.c_void_p, ctypes.c_uint]
lib.get_dyad.restype = ctypes.c_void_p

lib.get_len.argtypes = [ctypes.c_void_p]
lib.get_len.restype = ctypes.c_uint

lib.release_str.argtypes = [ctypes.c_void_p]
lib.release_str.restype = None

lib.simple_mean.argtypes = [ctypes.c_void_p, ctypes.c_uint]
lib.simple_mean.restype = ctypes.c_uint

lib.info_content.argtypes = [ctypes.c_void_p, ctypes.c_uint]
lib.info_content.restype = ctypes.c_float

class Gundam(object):
    names = None
    state = None

    def __init__(self, kmers, pos_fname, neg_fname):
        """
        kmers - list of kmer indices, format [(i, j, k, p)] where i and j are indices,
                k is the gap, and p is a float p-val
        pos_fname - name of file containing matching sequences
        neg_fname - name of file containing non-matching sequences
        """
        self.names = {}
        self.state = lib.read_kmers( json.dumps(kmers), pos_fname, neg_fname )

    def __len__(self):
        return lib.get_len(self.state)

    def __getitem__(self, idx):
        siz = len(self)
        if isinstance(idx, six.string_types):
            idx = self.names[ idx ]
        if idx >= siz:
            raise IndexError("{} > {}".format(idx,siz))
        p = lib.get_dyad(self.state, idx)
        # if we treat the char * as such, Python will try to release the memory,
        #   so instead we use this void * hack
        s = ctypes.cast(p, ctypes.c_char_p).value
        lib.release_str(p)
        return json.loads(s)

    def add_mean(self, idx, name=None):
        siz = len(self)
        if isinstance(idx, six.string_types):
            idx = self.names[ idx ]
        if idx >= siz:
            raise IndexError("{} > {}".format(idx,siz))

        new_idx = lib.simple_mean( self.state, idx )
        if name is not None:
            self.names[ name ] = new_idx
        return new_idx

    def info_content(self, idx):
        siz = len(self)
        if isinstance(idx, six.string_types):
            idx = self.names[ idx ]
        if idx >= siz:
            raise IndexError("{} > {}".format(idx,siz))

        return lib.info_content( self.state, idx )


def ctr(dyad):
    siz = dyad['motif']['scores']['dim'][0]
    subseqs = [tuple(seq[data['loc'] : data['loc'] + siz]) for  seq, data in dyad['pos_seqs']]
    return Counter(subseqs)

