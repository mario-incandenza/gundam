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
use dyad::*;

pub use dyad::find_motifs;

const KMER_LEN: usize = 5;
const MIN_GAP: usize = 0;
const MAX_GAP: usize = 6;
const MUT_INCR: f32 = 0.2;
const MIN_SCORE: f32 = 0.9;



#[cfg(test)]
mod tests {
    use super::*;
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

    fn choose(
        pos_v: &mut Vec<(Vec<u8>, f64)>,
        neg_v: &mut Vec<(Vec<u8>, f64)>,
    ) -> Option<(Vec<Vec<u8>>, Vec<Vec<u8>>)> {
        pos_v.sort_by(|&(_, score_a), &(_, score_b)| {
            score_b.partial_cmp(&score_a).expect("float sort")
        });
        neg_v.sort_by(|&(_, score_a), &(_, score_b)| {
            score_b.partial_cmp(&score_a).expect("float sort")
        });
        Some((
            pos_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(100)
                .collect(),
            neg_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(100)
                .collect(),
        ))
    }

    #[test]
    #[ignore]
    fn test_find() {
        let indiv_ct = 100;
        let dyads = DyadMotif::motifs(POS_FNAME, NEG_FNAME, MOTIF.len(), indiv_ct, choose);
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
}
