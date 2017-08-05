extern crate ndarray;
extern crate pssm;
extern crate darwin_rs;
extern crate bio;
extern crate rand;
extern crate jobsteal;

use jobsteal::{make_pool, BorrowSpliteratorMut, Spliterator};


use std::f64;
//use darwin_rs::{Individual, SimulationBuilder, PopulationBuilder};
use darwin_rs::{Individual, SimulationBuilder, Population, PopulationBuilder};

use pssm::{ Motif, BasePos };
use bio::io::fasta;
use ndarray::prelude::{Array, Array2};
use rand::Rng;

pub mod ctr;
use ctr::*;

const KMER_LEN: usize = 5;
const GAP_LEN: usize = 5;
const MUT_INCR: f32 = 0.08;
const MIN_SCORE: f32 = 0.9;


#[derive(Debug, Clone)]
struct DyadMotif<'a> {
    /// initial state based on kmers
    init:       Motif,
    /// weights updated by GA
    motif:      Motif,
    /// kmer counts
    ctr:        GappedKmerCtr,
    /// sequences containing the motif
    pos_fname:  &'a str,
    /// sequences without the motif
    neg_fname:  &'a str,
}

fn fasta_to_ctr( fname: &str ) -> GappedKmerCtr {
    let mut ctr = GappedKmerCtr::new( KMER_LEN, GAP_LEN, GAP_LEN);
    
    for _rec in fasta::Reader::from_file(fname).unwrap().records() {
        let rec = _rec.expect("couldn't unwrap record");
        ctr.update_with_seq(rec.seq());
    };
    ctr
}

fn kmers_to_matrix( kmer1: &[u8], gap_len: usize, kmer2: &[u8] ) -> Array2<f32> {
    let mut m = Array2::from_elem((kmer1.len() + gap_len + kmer2.len(), 4), 0.05);
    for i in 0 .. kmer1.len() {
        m[[i, BasePos::get(kmer1[i])]] = 0.8;
    }
    // set gap to N, ie, equal weights
    for i in  0 .. gap_len + 1 {
        for j in 0 .. 4 {
            m[[kmer1.len() + i, j]] = 0.25;
        }
    }
    for i in 0 .. kmer2.len() {
        m[[kmer1.len() + gap_len + i, BasePos::get(kmer2[i])]] = 0.8;
    }
    m
}

impl<'a> DyadMotif<'a> {
    pub fn population( pos_fname: &'a str, neg_fname: &'a str, width: usize, count: usize)
                       -> Vec<DyadMotif<'a>> {

        println!("processing pos");
        let pos = fasta_to_ctr(pos_fname);
        println!("processing neg");
        let neg = fasta_to_ctr(neg_fname);

        let (width, height, _) = pos.ctr.dim();
        let mut dyads = Vec::new();
        for i in 0 .. width {
            for j in 0 .. height {
                if pos.ctr[[i, j, 0]] > neg.ctr[[i, j, 0]]
                    && pos.ctr[[i, j, 0]] - neg.ctr[[i, j, 0]] >= 100 /*200*/ {

                        let init = Motif::from( kmers_to_matrix(
                            GappedKmerCtr::int_to_kmer(KMER_LEN, i).as_slice(),
                            GAP_LEN,
                            GappedKmerCtr::int_to_kmer(KMER_LEN, j).as_slice()) );

                        if init.min_score == init.max_score {
                            println!("{}.....{}",
                                     String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, i)).expect("AA"),
                                     String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, j)).expect("BB"));
                            continue
                        }

                        let copy = init.clone();

                        println!("{}.....{}",
                                 String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, i)).expect("AA"),
                                 String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, j)).expect("BB"));
                        
                        dyads.push( DyadMotif {
                            init:        init,
                            motif:       copy,
                            ctr:         GappedKmerCtr::new(KMER_LEN, GAP_LEN, GAP_LEN),
                            pos_fname:   pos_fname,
                            neg_fname:   neg_fname,
                        });
                    }
            }
        }
        dyads
    }
}

/// normalize scores in-place by summing each column and dividing each value
fn normalize_scores( scores: &mut Array2<f32> ) {
    let (width, bases) = scores.dim();

    for i in 0..width {
        let mut tot = 0.0;
        for j in 0..4 {
            tot += scores[[i,j]];
        }
        for j in 0..4 {
            scores[[i,j]] = scores[[i,j]] / tot;
        }
    }
}

impl<'a> Individual for DyadMotif<'a> {

    /// shift value at each position
    fn mutate(&mut self) {
        // Mutate the scores
        for x in self.motif.scores.iter_mut() {
            // r is on (0,1)
            let r = rand::random::<f32>();
            // by subtracting 0.5, we allow for negative random numbers
            // by scaling by 0.02, we limit changes to (-0.01,0.01)
            *x += MUT_INCR * (r - 0.5);
        }
        normalize_scores( &mut self.motif.scores );
        self.motif.calc_minmax();
    }

    fn calculate_fitness(&mut self) -> f64 {
        // Calculate how good the data values are compared to the perfect solution

        let mut pool = make_pool(3).unwrap();

        let mut v = Vec::new();
        for _rec in fasta::Reader::from_file(self.pos_fname)
            .expect(format!("couldn't open {}", &self.pos_fname).as_str())
            .records() {
                let rec = _rec.expect("unwrap record");
                v.push( (rec.seq().to_vec(), 0.0) );
            };

        v.split_iter_mut()
            .for_each(&pool.spawner(), |p| {

                //pos_ct += 1;

                match self.motif.score( &p.0 ) {
                    //Some((_, score)) if score >= MIN_SCORE => pos_pass += 1,
                    Some((_, score)) => { p.1 = score as f64; },
                    _ => (),
                }
            });

        let pos_score: f64 = v.iter().map(|&(_,score)| score).sum();
        let pos_ct = v.len();
        println!("score: {}", pos_score);

        let mut neg_ct: usize = 0;
        let mut neg_pass: usize = 0;
        let mut neg_score: f64 = 0.0;
        for _rec in fasta::Reader::from_file(self.neg_fname)
            .expect(format!("couldn't open {}", &self.neg_fname).as_str())
            .records() {
                neg_ct += 1;

                let rec = _rec.expect("unwrap record");
                match self.motif.score( rec.seq() ) {
                    //Some((_, score)) if score >= MIN_SCORE => neg_pass += 1,
                    Some((_, score)) => neg_score += score as f64,
                    _ => (),
                }

            }

        if pos_score == 0.0 {
            0.0
        } else if neg_score == 0.0 {
            f64::INFINITY
        } else {
            //let x = (pos_pass as f64 / pos_ct as f64) - (neg_pass as f64 / neg_ct as f64);
            //println!("({}/{}) - ({}/{}) =  {}", pos_pass, pos_ct, neg_pass, neg_ct, x);

            let x = (pos_score / pos_ct as f64) / (neg_score / neg_ct as f64);
            println!("({}/{}) / ({}/{}) =  {}", pos_score, pos_ct, neg_score, neg_ct, x);
            x
        }
    }

    /// initialize array with random values, then normalize
    /// so each position sums to 1.0
    fn reset(&mut self) {
        // bases == 4
        self.motif = self.init.clone();
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";

    #[test]
    fn kmers_to_m() {
        let m = kmers_to_matrix(b"ATGC", 1, b"ATGC");
        assert_eq!(m,
                   Array::from_vec(vec![0.8, 0.05, 0.05, 0.05,
                                        0.05, 0.8, 0.05, 0.05,
                                        0.05, 0.05, 0.8, 0.05,
                                        0.05, 0.05, 0.05, 0.8,
                                        0.25, 0.25, 0.25, 0.25,
                                        0.8, 0.05, 0.05, 0.05,
                                        0.05, 0.8, 0.05, 0.05,
                                        0.05, 0.05, 0.8, 0.05,
                                        0.05, 0.05, 0.05, 0.8,])
                   .into_shape((9,4))
                   .unwrap());
    }

    #[test]
    fn test_one() {
        let motif = Motif::from( kmers_to_matrix( b"ATAGG", GAP_LEN, b"CCATG" ));
        println!("score for present: {:?}", motif.score(b"GGAACGAAGTCCGTAGGGTCCATAGGAAAACCACTATGGGGCAGGATAATCATTAAAGGTCACTCGGTCGAGGCACAGATTGTGAGGAAGATGTAGGGGACCGTCGTTAAACCTAACGGACGGCTACACGGTTGTTGAAATGTCCCCCCCTTTTGCATTTTTCCTATGGGCGGCGACATAAAACTCGCAGACGAAGTTGGATATCTCCCGAATACGTGGACCGGCAGCATAACCAGACAAACGGGTAACTAACGTATGAGTGTGTCCAGCCACCATCCATAGGAAGTCCCATGAGTGAGCTTGATGATGTGAGGGCATGACATGTGCGGAAAACGAAGAACTAGGACCATAATGCAGGGCGACCTGCGCTCGAAACTCTGGATTACCATTTCCGCGGCCTAATATGGATCTCCTGTGTCTCGGATCCTTCAGGTCGACGTTCGGATCATACATGGGACTACAACGTGTCGATAGACCGCCAGACCTACACAAAGCATGCA"));
    }

    #[test]
    fn find_motif() {

        let motifs = DyadMotif::population("pos.fa", "neg.fa", MOTIF.len(), 10);

        let population1 = PopulationBuilder::<DyadMotif>::new()
            .set_id(1)
            .initial_population(&motifs)
            .increasing_exp_mutation_rate(1.03)
            .reset_limit_increment(100)
            .reset_limit_start(100)
            .reset_limit_end(1000)
            .finalize().unwrap();        

        let my_builder = SimulationBuilder::<DyadMotif>::new()
            .factor(0.34)
            .threads(2)
            .add_population(population1)
            .finalize();

        match my_builder {
            Err(_) => println!("more than 10 iteratons needed"),
            Ok(mut my_simulation) => {
                my_simulation.run();

                println!("total run time: {} ms", my_simulation.total_time_in_ms);
                println!("improvement factor: {}", my_simulation.simulation_result.improvement_factor);
                println!("number of iterations: {}", my_simulation.simulation_result.iteration_counter);

                my_simulation.print_fitness();
            }
        }
    }
}

