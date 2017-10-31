#[macro_use(s)]
extern crate ndarray;
extern crate pssm;
extern crate darwin_rs;
extern crate bio;
extern crate rand;
extern crate jobsteal;

use jobsteal::{make_pool, BorrowSpliteratorMut, Spliterator, Pool};


use std::f64;
use std::str;
use darwin_rs::{Individual, SimulationBuilder, Population, PopulationBuilder};

use pssm::{Motif, BasePos, ScoredPos};
use bio::io::fasta;
use ndarray::prelude::{Array, Array2};
use rand::Rng;

pub mod ctr;
use ctr::*;

const KMER_LEN: usize = 5;
const GAP_LEN: usize = 4;
const MUT_INCR: f32 = 0.08;
const MIN_SCORE: f32 = 0.9;


#[derive(Debug, Clone)]
struct DyadMotif<'a> {
    /// initial state based on kmers
    init: Motif,
    /// weights updated by GA
    motif: Motif,
    /// kmer counts
    ctr: GappedKmerCtr,
    /// sequences containing the motif
    pos_fname: &'a str,
    /// sequences without the motif
    neg_fname: &'a str,
}

fn fasta_to_ctr(fname: &str) -> GappedKmerCtr {
    let mut ctr = GappedKmerCtr::new(KMER_LEN, 0, GAP_LEN);

    for _rec in fasta::Reader::from_file(fname)
        .expect(format!("trouble opening {}", fname).as_str())
        .records()
    {
        let rec = _rec.expect("couldn't unwrap record");
        ctr.update_with_seq(rec.seq());
    }

    println!("@@ created ctr w/ dimensions: {:?}", ctr.ctr.shape());
    ctr
}

fn kmers_to_matrix(kmer1: &[u8], gap_len: usize, kmer2: &[u8]) -> Array2<f32> {
    let mut m = Array2::from_elem((kmer1.len() + gap_len + kmer2.len(), 4), 0.05);
    for i in 0..kmer1.len() {
        m[[i, BasePos::get(kmer1[i])]] = 0.8;
    }
    // set gap to N, ie, equal weights
    for i in 0..gap_len + 1 {
        for j in 0..4 {
            m[[kmer1.len() + i, j]] = 0.25;
        }
    }
    for i in 0..kmer2.len() {
        m[[kmer1.len() + gap_len + i, BasePos::get(kmer2[i])]] = 0.8;
    }
    m
}

impl<'a> DyadMotif<'a> {
    pub fn motifs(
        pos_fname: &'a str,
        neg_fname: &'a str,
        width: usize,
        count: usize,
    ) -> Vec<DyadMotif<'a>> {

        println!("processing pos");
        let pos = fasta_to_ctr(pos_fname);
        println!("processing neg");
        let neg = fasta_to_ctr(neg_fname);

        let (width, height, gap) = pos.ctr.dim();
        let mut dyads = Vec::new();
        for i in 0..width {
            for j in 0..height {
                for k in 0..gap {
                    if pos.ctr[[i, j, k]] > neg.ctr[[i, j, k]] &&
                        pos.ctr[[i, j, k]] - neg.ctr[[i, j, k]] >= 100
                    /*200*/
                    {

                        let init = Motif::from(kmers_to_matrix(
                            GappedKmerCtr::int_to_kmer(KMER_LEN, i).as_slice(),
                            k,
                            GappedKmerCtr::int_to_kmer(KMER_LEN, j).as_slice(),
                        ));

                        if init.min_score == init.max_score {
                            println!(
                                "{}.....{}",
                                String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, i))
                                    .expect("AA"),
                                String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, j))
                                    .expect("BB")
                            );
                            continue;
                        }

                        let copy = init.clone();

                        println!(
                            "{}.....{}",
                            String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, i)).expect("AA"),
                            String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, j)).expect("BB")
                        );

                        dyads.push(DyadMotif {
                            init: init,
                            motif: copy,
                            ctr: GappedKmerCtr::new(KMER_LEN, 0, GAP_LEN),
                            pos_fname: pos_fname,
                            neg_fname: neg_fname,
                        });
                    }
                }
            }
        }
        dyads
    }
}

/// normalize scores in-place by summing each column and dividing each value
fn normalize_scores(scores: &mut Array2<f32>) {
    let (width, bases) = scores.dim();

    for i in 0..width {
        let mut tot = 0.0;
        for j in 0..4 {
            tot += scores[[i, j]];
        }
        for j in 0..4 {
            scores[[i, j]] = scores[[i, j]] / tot;
        }
    }
}


fn crossover_motifs(
    mine: &mut Motif,
    theirs: &mut Motif,
    pos_fname: &str,
    neg_fname: &str,
) -> Motif {
    assert_eq!(mine.len(), theirs.len());

    // store sequences in memory, as IO doesn't play nice w/ parallelism
    let mut pos_seqs = Vec::new();
    for _rec in fasta::Reader::from_file(pos_fname)
        .expect(format!("couldn't open {}", pos_fname).as_str())
        .records()
    {
        let rec = _rec.expect("unwrap record");
        pos_seqs.push((rec.seq().to_vec(), ScoredPos::nil(), ScoredPos::nil()));
    }
    let mut neg_seqs = Vec::new();
    for _rec in fasta::Reader::from_file(neg_fname)
        .expect(format!("couldn't open {}", neg_fname).as_str())
        .records()
    {
        let rec = _rec.expect("unwrap record");
        neg_seqs.push((rec.seq().to_vec(), ScoredPos::nil(), ScoredPos::nil()));
    }

    // step one: reduce input to match the width of motif
    // this necessarily means chosing between the best match position from @mine or @theirs
    fn max_slice(
        pwm_len: usize,
        mine: &Motif,
        theirs: &Motif,
        seq: &[u8],
    ) -> (ScoredPos, ScoredPos) {
        let s = mine.score(seq).expect(
            format!("self couldn't score: {:?}", seq)
                .as_str(),
        );
        let o = theirs.score(seq).expect(
            format!("other couldn't score: {:?}", seq)
                .as_str(),
        );

        if s.loc == o.loc {
            (s, o)
        } else if o.sum > s.sum {
            let x = mine.score(&seq[o.loc..o.loc + pwm_len]).expect(
                "couldn't score slice (1)",
            );
            (x, o)
        } else {
            let x = theirs.score(&seq[s.loc..s.loc + pwm_len]).expect(
                "couldn't score slice (2)",
            );
            (s, x)
        }
    }

    let mut pool = make_pool(3).unwrap();
    let pwm_len = mine.len();
    pos_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
        println!(
            "pos ({}): {:?} -> {}",
            str::from_utf8(mine.degenerate_consensus().as_slice()).unwrap(),
            str::from_utf8(&t.0[scores.0.loc..scores.0.loc + pwm_len]).unwrap(),
            scores.0.sum
        );
        t.1 = scores.0;
        t.2 = scores.1;
    });
    neg_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
        println!(
            "neg ({}): {:?} -> {}",
            str::from_utf8(mine.degenerate_consensus().as_slice()).unwrap(),
            str::from_utf8(&t.0[scores.0.loc..scores.0.loc + pwm_len]).unwrap(),
            scores.0.sum
        );
        t.1 = scores.0;
        t.2 = scores.1;
    });

    // step 2: create new PWM by choosing the best base at each position
    let mut new_m = Array2::zeros((pwm_len, 4));
    for i in 0..pwm_len {
        // positive set
        let mut s_ptally = 0.0;
        let mut o_ptally = 0.0;
        for &(_, ref sscore, ref oscore) in pos_seqs.iter() {
            s_ptally += sscore.scores[i];
            o_ptally += oscore.scores[i];
        }

        // negative set
        let mut s_ntally = 0.0;
        let mut o_ntally = 0.0;
        for &(_, ref sscore, ref oscore) in neg_seqs.iter() {
            s_ntally += sscore.scores[i];
            o_ntally += oscore.scores[i];
        }

        let mut __c: usize = 0;
        for b in 0..4 {
            new_m[[i, b]] = if (o_ptally / o_ntally) > (s_ptally / s_ntally) {
                __c += 1;
                theirs.scores[[i, b]]
            } else {
                mine.scores[[i, b]]
            }
        }
        if __c != 0 && __c != 4 {
            println!("@@ weird mixed");
        }
    }

    Motif::from(new_m)
}

impl<'a> Individual for DyadMotif<'a> {
    const CAN_CROSSOVER: bool = true;

    /// shift value at each position
    fn mutate(&mut self) {
        // Mutate the scores
        for x in self.motif.scores.iter_mut() {
            // r is on (0,1)
            let r = rand::random::<f32>();
            // by subtracting 0.5, we allow for negative random numbers
            // by scaling by 0.02, we limit changes to (-0.01,0.01)
            let new_x = *x + MUT_INCR * (r - 0.5);
            *x = if new_x < 0.0 { 0.0 } else { new_x };
        }
        normalize_scores(&mut self.motif.scores);
        self.motif.calc_minmax();
    }

    fn calculate_fitness(&mut self) -> f64 {
        // Calculate how good the data values are compared to the perfect solution

        let mut pool = make_pool(3).unwrap();

        fn eval_file(pool: &mut Pool, motif: &Motif, fname: &str) -> f64 {

            // FIXME: b/c we wind up re-analyzing these files again and again,
            // we should probably just read into memory once and be done w/ it
            let mut v = Vec::new();
            for _rec in fasta::Reader::from_file(fname)
                .expect(format!("couldn't open {}", fname).as_str())
                .records()
            {
                let rec = _rec.expect("unwrap record");
                v.push((rec.seq().to_vec(), 0.0));
            }

            if v.len() == 0 {
                panic!("empty file: {}", fname);
            }

            v.split_iter_mut().for_each(&pool.spawner(), |p| {

                match motif.score(&p.0) {
                    //Some((_, score)) if score >= MIN_SCORE => pos_pass += 1,
                    Some(ScoredPos { ref sum, .. }) => {
                        p.1 = *sum as f64;
                    }
                    _ => (),
                }
            });

            let pos_score: f64 = v.iter().map(|&(_, score)| score).sum();
            let pos_ct = v.len();

            (pos_score / pos_ct as f64)
        }

        let pos = eval_file(&mut pool, &self.motif, self.pos_fname);
        let neg = eval_file(&mut pool, &self.motif, self.neg_fname);

        //println!("fitness (motif[0]={:?}): {}", [self.motif.scores[[0,0]], self.motif.scores[[0,1]]], neg / pos);

        if pos == 0.0 {
            f64::INFINITY
        } else if neg == 0.0 {
            0.0
        } else {
            neg / pos
        }
    }

    /// initialize array with random values, then normalize
    /// so each position sums to 1.0
    fn reset(&mut self) {
        println!("-- reset");
        // bases == 4
        self.motif = self.init.clone();
    }



    fn crossover(&mut self, other: &mut Self) -> Self {
        let new_motif = crossover_motifs(
            &mut self.motif,
            &mut other.motif,
            self.pos_fname,
            self.neg_fname,
        );

        DyadMotif {
            init: self.motif.clone(),
            motif: new_motif,
            ctr: self.ctr.clone(),
            pos_fname: self.pos_fname,
            neg_fname: self.neg_fname,
        }
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use darwin_rs::select::MaximizeSelector;
    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";

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
    fn test_one() {
        let motif = Motif::from(kmers_to_matrix(b"ATAGG", GAP_LEN, b"CCATG"));
        println!("score for present: {:?}", motif.score(b"GGAACGAAGTCCGTAGGGTCCATAGGAAAACCACTATGGGGCAGGATAATCATTAAAGGTCACTCGGTCGAGGCACAGATTGTGAGGAAGATGTAGGGGACCGTCGTTAAACCTAACGGACGGCTACACGGTTGTTGAAATGTCCCCCCCTTTTGCATTTTTCCTATGGGCGGCGACATAAAACTCGCAGACGAAGTTGGATATCTCCCGAATACGTGGACCGGCAGCATAACCAGACAAACGGGTAACTAACGTATGAGTGTGTCCAGCCACCATCCATAGGAAGTCCCATGAGTGAGCTTGATGATGTGAGGGCATGACATGTGCGGAAAACGAAGAACTAGGACCATAATGCAGGGCGACCTGCGCTCGAAACTCTGGATTACCATTTCCGCGGCCTAATATGGATCTCCTGTGTCTCGGATCCTTCAGGTCGACGTTCGGATCATACATGGGACTACAACGTGTCGATAGACCGCCAGACCTACACAAAGCATGCA"));
    }

    #[test]
    fn find_motif() {

        for dyad in DyadMotif::motifs("pos.fa", "neg.fa", MOTIF.len(), 1000) {

            println!("-- motif: {:?}", &dyad.motif);

            // make an initial population of 100 copies of the motif
            let mut init = (0..100).map(|_| dyad.clone()).collect::<Vec<DyadMotif>>();
            for ind in init.iter_mut() {
                ind.mutate();
            }


            let population1 = PopulationBuilder::<DyadMotif>::new()
                .set_id(1)
                .initial_population(&init)
                .increasing_exp_mutation_rate(1.03)
                .reset_limit_increment(100)
                .reset_limit_start(100)
                .reset_limit_end(1000)
                .finalize()
                .unwrap();
            println!("-- population built");
            let my_builder = SimulationBuilder::<DyadMotif>::new()
                .factor(0.34)
                .threads(2)
                .add_population(population1)
                .finalize();
            println!("-- builder");
            let selector = MaximizeSelector::new(80);

            match my_builder {
                Err(_) => println!("more than 10 iteratons needed"),
                Ok(mut my_simulation) => {
                    my_simulation.run(&selector);

                    println!("total run time: {} ms", my_simulation.total_time_in_ms);
                    println!(
                        "very fittest: {}",
                        my_simulation.simulation_result.fittest[0].fitness
                    );
                    println!(
                        "improvement factor: {}",
                        my_simulation.simulation_result.improvement_factor
                    );
                    println!(
                        "number of iterations: {}",
                        my_simulation.simulation_result.iteration_counter
                    );

                    my_simulation.print_fitness();
                }
            }
        }
    }

    #[test]
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
