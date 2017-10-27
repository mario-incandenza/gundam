
use std::usize;
use std::str;
use std::ops::Index;
use std::convert::AsRef;
use ndarray::Array3;
use pssm::BasePos;

#[derive(Debug, Clone)]
pub struct GappedKmerCtr {
    pub ctr: Array3<usize>,
    pub kmer_len: usize,
    pub min_gap: usize,
    pub max_gap: usize,
}

impl GappedKmerCtr {
    /// simple bitshifting conversion of sequences to integer
    /// note that this only makes sense for fixed-width kmers, as, eg,
    /// kmer_to_int(b"A") == kmer_to_int(b"AA") == kmer_to_int(b"AAA") == 0
    pub fn kmer_to_int(kmer: &[u8]) -> usize {
        let mut val: usize = 0;
        for b in kmer.iter() {
            val <<= 2;
            val |= BasePos::get(*b);
        }
        val
    }

    ///
    pub fn int_to_kmer(kmer_len: usize, i: usize) -> Vec<u8> {
        let mut kmer: Vec<u8> = Vec::new();
        let mut val = i;
        for _ in 0..kmer_len {
            kmer.push(BasePos::put(val & 0b11));
            val >>= 2;
        }
        kmer.reverse();
        kmer
    }

    pub fn new(kmer_len: usize, min_gap: usize, max_gap: usize) -> GappedKmerCtr {
        let len: usize = 4usize.pow(kmer_len as u32);
        GappedKmerCtr {
            ctr: Array3::zeros((len, len, max_gap - min_gap + 1)),
            kmer_len: kmer_len,
            min_gap: min_gap,
            max_gap: max_gap,
        }
    }

    /// eg, foo.get(b"ATGC", b"TTCA", 0) -> 10
    pub fn get(&self, kmer1: &[u8], kmer2: &[u8], gap: usize) -> usize {
        self.ctr[[
            GappedKmerCtr::kmer_to_int(kmer1),
            GappedKmerCtr::kmer_to_int(kmer2),
            gap,
        ]]
    }

    /// increment cell and return new value
    pub fn incr(&mut self, kmer1: &[u8], kmer2: &[u8], gap: usize) -> usize {
        let idx = [
            GappedKmerCtr::kmer_to_int(kmer1),
            GappedKmerCtr::kmer_to_int(kmer2),
            gap,
        ];
        self.ctr[idx] += 1;
        self.ctr[idx]
    }

    /// given a sequence, increment all kmers
    pub fn update_with_seq(&mut self, seq: &[u8]) {
        // annoying - work around borrow-checking
        let kmer_len = self.kmer_len;
        let min_gap = self.min_gap;
        let max_gap = self.max_gap;

        for gap in min_gap..max_gap + 1 {
            for i in 0..1 + seq.len() - (2 * kmer_len + gap) {
                self.incr(
                    &seq[i..i + kmer_len],
                    &seq[i + kmer_len + gap..i + 2 * kmer_len + gap],
                    gap - min_gap,
                );
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const SEQ: &'static [u8] = b"TGCACGGT";
    #[test]
    fn kmer_to_int() {
        assert_eq!(GappedKmerCtr::kmer_to_int(b"A"), 0);
        assert_eq!(GappedKmerCtr::kmer_to_int(b"T"), 1);
        assert_eq!(GappedKmerCtr::kmer_to_int(b"G"), 2);
        assert_eq!(GappedKmerCtr::kmer_to_int(b"C"), 3);
        assert_eq!(GappedKmerCtr::kmer_to_int(b"AA"), 0);
        assert_eq!(GappedKmerCtr::kmer_to_int(b"TA"), 4);

        assert_eq!(GappedKmerCtr::int_to_kmer(1, 0), b"A");
        assert_eq!(GappedKmerCtr::int_to_kmer(1, 1), b"T");
        assert_eq!(GappedKmerCtr::int_to_kmer(1, 2), b"G");
        assert_eq!(GappedKmerCtr::int_to_kmer(1, 3), b"C");
        assert_eq!(GappedKmerCtr::int_to_kmer(2, 0), b"AA");
        assert_eq!(GappedKmerCtr::int_to_kmer(2, 4), b"TA");
    }

    #[test]
    fn index() {
        // TGC -> 27
        // CGG -> 58
        // GCA -> 44
        // GGT -> 41

        let mut ctr = GappedKmerCtr::new(3, 1, 1);
        ctr.update_with_seq(SEQ);

        assert_eq!(ctr.get(b"TGC", b"CGG", 0), 1);

        for i in 0..64 {
            for j in 0..64 {
                assert_eq!(
                    ctr.ctr[[i, j, 0]],
                    match (i, j) {
                        (27, 58) => 1,
                        (44, 41) => 1,
                        _ => 0,
                    }
                );
            }
        }
    }
}

/*

impl AsRef<[usize; 3]> for &(&[u8], &[u8], usize) {
    fn as_ref(&self) -> &T;
}

impl Index<KmerCtrIdx> for GappedKmerCtr {
    type Output = usize;

    fn index(&self, idx: KmerCtrIdx) -> &usize {
        self.ctr[[0,0,0]]
    }
}

impl From<
*/
