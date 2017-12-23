#[macro_use(s)]
extern crate ndarray;
extern crate pssm;
extern crate darwin_rs;
extern crate bio;
extern crate rand;
extern crate jobsteal;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate fishers_exact;
extern crate num_cpus;
#[macro_use]
extern crate lazy_static;

use jobsteal::{make_pool, BorrowSpliteratorMut, Spliterator, Pool};
use std::f64;
use std::str;
use darwin_rs::{Individual, SimulationBuilder, Population, PopulationBuilder};
//use darwin_rs::select::MaximizeSelector;

use pssm::{Motif, BasePos, ScoredPos};
use bio::io::fasta;
use ndarray::prelude::{Array, Array2};
use rand::Rng;

pub mod ctr;
use ctr::*;

pub mod dyad;
pub use dyad::*;

pub use dyad::find_motifs;

const KMER_LEN: usize = 5;
const MIN_GAP: usize = 0;
const MAX_GAP: usize = 20;
const MUT_INCR: f32 = 0.2;
const MIN_SCORE: f32 = 0.9;

lazy_static! {
    pub static ref CPU_COUNT: usize = num_cpus::get();
}

/*
type ScoredSeqs = Vec<(String,f32)>;

pub fn export_scores(dyad: &DyadMotif) -> (ScoredSeqs, ScoredSeqs)

pub fn get_dyad(*Vec<DyadMotif>, idx) -> DyadMotif (serde)

pub fn refine_GA(*Vec..., idx) -> idx of new

pub fn refine_mean(*Vec..., idx) -> idx of new

pub fn winnow_seqs(*Vec, idx, Vec<idx>)
*/

#[cfg(test)]
mod tests {
    use super::*;
    use fishers_exact::{fishers_exact, TestTails};
    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";
    const POS_FNAME: &'static str = "pos.fa";
    const NEG_FNAME: &'static str = "neg.fa";

    #[test]
    #[ignore]
    fn kmers_to_m() {
        let m = kmers_to_matrix(b"ATGC", 1, b"ATGC");
        let expected = Array::from_vec(vec![
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
            0.25,
            0.25,
            0.25,
            0.25,
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
            0.05,
            0.05,
            0.05,
            0.05,
            0.8,
        ]).into_shape((9, 4))
            .unwrap();
        println!("diff: {:?}", m.clone() - expected.clone());
        assert_eq!(m, expected);
    }

    #[test]
    #[ignore]
    fn test_one() {
        let motif = Motif::from(kmers_to_matrix(b"ATAGG", MAX_GAP, b"CCATG"));
        println!("score for present: {:?}", motif.score(b"GGAACGAAGTCCGTAGGGTCCATAGGAAAACCACTATGGGGCAGGATAATCATTAAAGGTCACTCGGTCGAGGCACAGATTGTGAGGAAGATGTAGGGGACCGTCGTTAAACCTAACGGACGGCTACACGGTTGTTGAAATGTCCCCCCCTTTTGCATTTTTCCTATGGGCGGCGACATAAAACTCGCAGACGAAGTTGGATATCTCCCGAATACGTGGACCGGCAGCATAACCAGACAAACGGGTAACTAACGTATGAGTGTGTCCAGCCACCATCCATAGGAAGTCCCATGAGTGAGCTTGATGATGTGAGGGCATGACATGTGCGGAAAACGAAGAACTAGGACCATAATGCAGGGCGACCTGCGCTCGAAACTCTGGATTACCATTTCCGCGGCCTAATATGGATCTCCTGTGTCTCGGATCCTTCAGGTCGACGTTCGGATCATACATGGGACTACAACGTGTCGATAGACCGCCAGACCTACACAAAGCATGCA"));
    }


    #[test]
    #[ignore]
    fn test_find() {
        let v = DyadMotif::passing_kmers(POS_FNAME, NEG_FNAME);
        let dyads = DyadMotif::motifs(v, POS_FNAME, NEG_FNAME, choose);
        let (score, new_dyad) = dyads[0].refine(100);

    }

    #[test]
    #[ignore]
    fn print_kmers() {
        for i in 0..MOTIF.len() - KMER_LEN {
            println!(
                "@@ from motif, kmer {} -> {}",
                str::from_utf8(&MOTIF[i..i + KMER_LEN]).unwrap(),
                GappedKmerCtr::kmer_to_int(&MOTIF[i..i + KMER_LEN])
            );
        }
    }
    #[test]
    fn fisher() {
        println!(
            "fisher: {:?}",
            fishers_exact(&[100, 200, 5000, 10000], TestTails::One)
        );
    }
}
