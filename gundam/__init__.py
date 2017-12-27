import ctypes
import argparse
import six
import json


lib = ctypes.CDLL("target/release/libgundam.dylib")

lib.read_kmers.argtypes = [ctypes.c_char_p, ctypes.c_char_p, ctypes.c_char_p]
lib.read_kmers.restype = ctypes.c_void_p

lib.get_dyad.argtypes = [ctypes.c_void_p, ctypes.c_uint]
lib.get_dyad.restype = ctypes.c_void_p

lib.get_len.argtypes = [ctypes.c_void_p]
lib.get_len.restype = ctypes.c_uint

class Gundam(object):
    state = None

    def __init__(self, kmers, pos_fname, neg_fname):
        """
        kmers - list of kmer indices, format [(i, j, k, p)] where i and j are indices,
                k is the gap, and p is a float p-val
        pos_fname - name of file containing matching sequences
        neg_fname - name of file containing non-matching sequences
        """
        self.state = lib.read_kmers( json.dumps(kmers), pos_fname, neg_fname )

    def __len__(self):
        return lib.get_len(self.state)

    def __getitem__(self, idx):
        siz = len(self)
        if idx >= siz:
            raise IndexError("{} > {}".format(idx,siz))
        s = ctypes.cast(lib.get_dyad(self.state, idx), ctypes.c_char_p).value
        return json.loads(s)
