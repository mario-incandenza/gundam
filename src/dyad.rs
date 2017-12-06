use jobsteal::{make_pool, BorrowSpliterator, Spliterator, Pool};
use std::f64;
use std::str;
use std::cmp::max;
use darwin_rs::{Individual, SimulationBuilder, Population, PopulationBuilder};
//use darwin_rs::select::MaximizeSelector;
use fishers_exact::{fishers_exact, TestTails};

use pssm::{Motif, BasePos, ScoredPos};
use bio::io::fasta;
use ndarray::prelude::{Array, Array2};
use rand;
use rand::Rng;

use ctr::*;
use super::*;

const P_CUTOFF: f64 = 0.001;

#[derive(Debug, Clone)]
pub struct DyadMotif {
    /// initial state based on kmers
    init: Motif,
    /// weights updated by GA
    pub motif: Motif,
    /// kmer len
    kmer_len: usize,
    /// gap len
    gap_len: usize,
    /// kmer counts
    ctr: GappedKmerCtr,
    /// sequences matching our motif
    pos_seqs: Vec<Vec<u8>>,
    /// sequences representing background
    neg_seqs: Vec<Vec<u8>>,
}

fn fasta_to_ctr(fname: &str) -> (GappedKmerCtr, usize) {
    let mut ctr = GappedKmerCtr::new(KMER_LEN, MIN_GAP, MAX_GAP);
    let mut tot = 0;

    for _rec in fasta::Reader::from_file(fname)
        .expect(format!("trouble opening {}", fname).as_str())
        .records()
    {
        let rec = _rec.expect("couldn't unwrap record");
        ctr.update_with_seq(rec.seq());
        tot += 1;
    }

    (ctr, tot)
}

pub fn kmers_to_matrix(kmer1: &[u8], gap_len: usize, kmer2: &[u8]) -> Array2<f32> {
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

impl DyadMotif {
    /// P-values returned by the Fisher exact test don't change meaningfully as the values get larger.
    /// eg, both [100, 200, 10_000, 10_000] and [1_000, 2_000, 100_000, 100_000] yield a P-value well
    /// below our cutoff.  therefore, we can safely scale the values down if they're above some arbitrary
    /// threshold.
    fn scaled_fisher(_ct1: usize, _tot1: usize, _ct2: usize, _tot2: usize) -> f64 {
        let (ct1, tot1) = if _tot1 as f64 <= 1e4 {
            (_ct1 as i32, _tot1 as i32)
        } else {
            let scale = 1e4 / _tot1 as f64;
            (max(1, (scale * _ct1 as f64) as i32), 10_000)
        };

        let (ct2, tot2) = if _tot2 as f64 <= 1e4 {
            (_ct2 as i32, _tot2 as i32)
        } else {
            let scale = 1e4 / _tot2 as f64;
            (max(1, (scale * _ct2 as f64) as i32), 10_000)
        };

        fishers_exact(&[ct1, ct2, tot1, tot2], TestTails::One)
    }

    /// generate kmers, tablulate, and apply Fisher exact test
    pub fn passing_kmers(pos_fname: &str, neg_fname: &str) -> Vec<(usize, usize, usize, f64)> {
        let (pos, pos_ct) = fasta_to_ctr(pos_fname);
        let (neg, neg_ct) = fasta_to_ctr(neg_fname);

        let (width, height, gap) = pos.ctr.dim();
        let mut dyads = Vec::new();
        for i in 0..width {
            info!("i={}", i);
            for j in 0..height {
                for k in 0..gap {
                    if pos.ctr[[i, j, k]] > neg.ctr[[i, j, k]] {
                        let p = DyadMotif::scaled_fisher(
                            pos.ctr[[i, j, k]],
                            pos_ct,
                            neg.ctr[[i, j, k]],
                            neg_ct,
                        );
                        if p < P_CUTOFF {
                            dyads.push((i, j, k, p));
                        };
                    }
                }
            }
        }
        dyads
    }

    pub fn motifs<F>(
        chosen: Vec<(usize, usize, usize, f64)>,
        pos_fname: &str,
        neg_fname: &str,
        chooser: F,
    ) -> Vec<DyadMotif>
    where
        F: Fn(&mut Vec<(Vec<u8>, f64)>, &mut Vec<(Vec<u8>, f64)>) -> Option<(Vec<Vec<u8>>, Vec<Vec<u8>>)>,
    {
        info!("using {} cpus", *CPU_COUNT);
        let mut pool = make_pool(*CPU_COUNT).unwrap();
        let mut dyads = Vec::new();
        for (idx, &(i, j, k, _)) in chosen.iter().enumerate() {
            if idx % 500 == 0 {
                info!("creating dyad #{}", idx);
            }

            let init = Motif::from(kmers_to_matrix(
                GappedKmerCtr::int_to_kmer(KMER_LEN, i).as_slice(),
                k,
                GappedKmerCtr::int_to_kmer(KMER_LEN, j).as_slice(),
            ));

            if init.min_score == init.max_score {
                info!("skipping motif: {}<gap={}>{}", 
                                String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, i))
                                    .expect("DyadMotif::motifs - A"),
k,
                                String::from_utf8(GappedKmerCtr::int_to_kmer(KMER_LEN, j))
                                    .expect("DyadMotif::motifs - B"),
                            );

                continue;
            }

            let copy = init.clone();

            let mut pos_v = DyadMotif::eval_file(&mut pool, &copy, pos_fname);
            let mut neg_v = DyadMotif::eval_file(&mut pool, &copy, neg_fname);
            let (pos_seqs, neg_seqs) =
                chooser(&mut pos_v, &mut neg_v).expect("motifs found bad one (1)");

            info!("DyadMotif::motifs - {} / {} seqs used", pos_seqs.len(), pos_v.len());

            dyads.push(DyadMotif {
                init: init,
                motif: copy,
                ctr: GappedKmerCtr::new(KMER_LEN, MIN_GAP, MAX_GAP),
                kmer_len: KMER_LEN,
                gap_len: MIN_GAP + k,
                pos_seqs: pos_seqs,
                neg_seqs: neg_seqs,
            });
        }

        dyads
    }

    /// apply motif to sequences in a FASTA file, returning sequences and scores
    fn eval_file(pool: &mut Pool, motif: &Motif, fname: &str) -> Vec<(Vec<u8>, f64)> {

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
        v
    }

    ///
    pub fn refine(&self, mut_ct: usize) -> (f64, DyadMotif) {

        // make an initial population of 100 copies of the motif
        let mut init = (0..mut_ct)
            .map(|_| self.clone())
            .collect::<Vec<DyadMotif>>();
        for ind in init.iter_mut() {
            ind.mutate();
        }


        let population1 = PopulationBuilder::<DyadMotif>::new()
                .set_id(1)
                .initial_population(&init)
                .increasing_exp_mutation_rate(1.03)
                .reset_limit_end(0) // disable resetting
                .finalize()
                .expect("PopulationBuilder");
        let mut sim = SimulationBuilder::<DyadMotif>::new()
            .iterations(11)        //.factor(0.34)
            .threads(*CPU_COUNT)
            .add_population(population1)
            .finalize()
            .expect("some problem making builder");

        sim.run();
        sim.print_fitness();

        (
            sim.simulation_result.fittest[0].fitness,
            sim.simulation_result.fittest[0].individual.clone(),
        )
    }

    /// stringified self.motif.degenerate_consensus()
    pub fn show_motif(&self) -> String {
        String::from_utf8(self.motif.degenerate_consensus()).expect("show_motif")
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

        let t = if s.loc == o.loc {
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
        };

        t
    }

    let mut pool = make_pool(*CPU_COUNT).unwrap();
    let pwm_len = mine.len();
    pos_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
        t.1 = scores.0;
        t.2 = scores.1;
    });
    neg_seqs.split_iter_mut().for_each(&pool.spawner(), |t| {
        let scores = max_slice(pwm_len, mine, theirs, t.0.as_slice());
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
        //println!("@ i={}, o > s? {:?}", i, (o_ptally / o_ntally) > (s_ptally / s_ntally));
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

impl Individual for DyadMotif {
    //const CAN_CROSSOVER: bool = false;

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

        let mut pool = make_pool(*CPU_COUNT).unwrap();

        let pos_sum: f64 = self.pos_seqs
            .clone()
            .split_iter()
            .map(|seq| self.motif.score(&seq).expect("score?").sum as f64)
            .collect::<Vec<f64>>(&pool.spawner())
            .iter()
            .sum();
        let pos: f64 = pos_sum / self.pos_seqs.len() as f64;


        let neg_sum: f64 = self.neg_seqs
            .clone()
            .split_iter()
            .map(|seq| self.motif.score(&seq).expect("score?").sum as f64)
            .collect::<Vec<f64>>(&pool.spawner())
            .iter()
            .sum();
        let neg: f64 = neg_sum / self.neg_seqs.len() as f64;
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



    /* FIXME: switched from crossover fork to regular darwin-rs

    fn crossover(&mut self, other: &mut Self) -> Self {
        info!("DyadMotif::crossover");
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
        self.clone()
    }*/
}

pub fn find_motifs(
    chosen: Vec<(usize, usize, usize, f64)>,
    pos_fname: &str,
    neg_fname: &str,
) -> Vec<DyadMotif> {
    let indiv_ct = 100;

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
        let mut cutoff = 0;
        for (i, &(_, score)) in pos_v.iter().enumerate() {
            if score <= 0.9 {
                break;
            }
            cutoff = i;
        }
        if cutoff == 0 {
            return None;
        }
        Some((
            pos_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(cutoff)
                .collect(),
            neg_v
                .iter()
                .map(|&(ref s, _)| s.clone())
                .take(cutoff)
                .collect(),
        ))
    }

    let mut pool = make_pool(*CPU_COUNT).unwrap();
    let motifs = DyadMotif::motifs(chosen, pos_fname, neg_fname, choose);
    info!("got {} motifs", motifs.len());
    motifs
        .iter()
        .enumerate()
        .map(|(idx, dyad)| {

            // dyad wraps the sequences chosen by our method
            info!(
                "unrefined motif #{}: {} (gap_len {})",
                idx,
                dyad.show_motif(),
                dyad.gap_len
            );
            let (score, mut new_dyad) = dyad.refine(100);
            info!(
                "motif #{} after refine: {}, score={}",
                idx,
                dyad.show_motif(),
                score
            );

            // now that the GA has [hopefully] improved our PWM, we need to choose new seqs
            let mut pos_v = DyadMotif::eval_file(&mut pool, &new_dyad.motif, pos_fname);
            let mut neg_v = DyadMotif::eval_file(&mut pool, &new_dyad.motif, neg_fname);
            let (pos_seqs, neg_seqs) =
                choose(&mut pos_v, &mut neg_v).expect("motifs found bad one (2)");
            info!("find_motifs - motif #{} - {} / {} seqs used", idx, pos_seqs.len(), pos_v.len());
            new_dyad.pos_seqs = pos_seqs;
            new_dyad.neg_seqs = neg_seqs;

            //new_dyad.repopulate_seqs(pos_fname, neg_fname, choose);
            let (score2, new2) = new_dyad.refine(100);
            info!(
                "motif #{} after second refine: {}, score={}",
                idx,
                new2.show_motif(),
                score2
            );
            new2
        })
        .collect()
}



#[cfg(test)]
mod tests {
    use super::*;
    const MOTIF: &'static [u8] = b"GGCCTAGCCATG";
    //const POS_FNAME: &'static str = "pos.fa"; // all contain MOTIF
    const POS_FNAME: &'static str = "test.fa"; // various motifs at various frequencies
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
    fn test_find_one_motif() {
        println!("dyad::test_find");
        let v = DyadMotif::passing_kmers(POS_FNAME, NEG_FNAME);
        let dyads = DyadMotif::motifs(v, POS_FNAME, NEG_FNAME, choose);
        dyads[0].refine(100);
    }

    #[test]
    fn test_find_motifs() {
        env_logger::init();
        let v = DyadMotif::passing_kmers(POS_FNAME, NEG_FNAME);
        find_motifs(v, POS_FNAME, NEG_FNAME);
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
