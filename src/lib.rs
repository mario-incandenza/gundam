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
#[macro_use]
extern crate serde_derive;

extern crate serde;
extern crate serde_json;

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_void};
use std::mem;

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

type ScoredSeqs = Vec<(String, f32)>;

/*

interface ideas:

pub fn export_scores(dyad: &DyadMotif) -> (ScoredSeqs, ScoredSeqs)

pub fn get_dyad(*Vec<DyadMotif>, idx) -> DyadMotif (serde)

pub fn refine_GA(*Vec..., idx) -> idx of new

pub fn refine_mean(*Vec..., idx) -> idx of new

pub fn winnow_seqs(*Vec, idx, Vec<idx>)

pub fn shannon_entropy( Vec<seq> ) -> f32
*/


#[no_mangle]
pub extern "C" fn release_str(somestr: *mut c_char) {
    unsafe {
        CString::from_raw(somestr);
    }
}

/// wrapper for DyadMotif::motifs
#[no_mangle]
pub extern "C" fn read_kmers(
    _kmers_s: *const c_char,
    _pos_fname: *const c_char,
    _neg_fname: *const c_char,
) -> *mut c_void {
    let kmers_s = unsafe { CStr::from_ptr(_kmers_s).to_string_lossy().into_owned() };
    let pos_fname = unsafe { CStr::from_ptr(_pos_fname).to_string_lossy().into_owned() };
    let neg_fname = unsafe { CStr::from_ptr(_neg_fname).to_string_lossy().into_owned() };
    let kmers: Vec<(usize, usize, usize, f64)> =
        serde_json::from_str(&kmers_s).expect("deser kmers");
    let motifs = Box::new(DyadMotif::motifs(kmers, &pos_fname, &neg_fname, choose));
    Box::into_raw(motifs) as *mut c_void
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    DyadMotif at position idx, encoded as JSON
#[no_mangle]
pub extern "C" fn get_dyad(_dyads: *const c_void, idx: u32) -> *const c_char {
    let dyads: &Vec<DyadMotif> = unsafe { mem::transmute(_dyads) };

    CString::new(serde_json::to_string(&dyads[idx as usize]).expect(
        "get_dyad - ser",
    )).expect("get_dyad")
        .into_raw()
}

/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
/// output:
///    length of vec
#[no_mangle]
pub extern "C" fn get_len(_dyads: *const c_void) -> u32 {
    let dyads: &Vec<DyadMotif> = unsafe { mem::transmute(_dyads) };
    dyads.len() as u32
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    information content of motif
#[no_mangle]
pub extern "C" fn info_content(_dyads: *const c_void, idx: u32) -> f32 {
    let dyads: &Vec<DyadMotif> = unsafe { mem::transmute(_dyads) };

    dyads[idx as usize].motif.info_content()
}


/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    information content of motif
#[no_mangle]
pub extern "C" fn show_motif(_dyads: *const c_void, idx: u32) -> *const c_char {
    let dyads: &Vec<DyadMotif> = unsafe { mem::transmute(_dyads) };

    CString::new(dyads[idx as usize].show_motif())
        .expect("get_dyad")
        .into_raw()
}


/// creates a new DyadMotif and appends it to the vec, returning
/// indes of new motif
/// input:
///    _dyads - Vec<DyadMotif> as returned by read_kmers
///    idx - u32
/// output:
///    index of new motif
#[no_mangle]
pub extern "C" fn simple_mean(_dyads: *const c_void, idx: u32) -> u32 {
    let dyads: &mut Vec<DyadMotif> = unsafe { mem::transmute(_dyads) };

    let new = dyads[idx as usize].refine_mean();
    let new_idx = dyads.len();
    dyads.push(new);
    new_idx as u32
}



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
        let new_dyad = dyads[0].refine(100);

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
