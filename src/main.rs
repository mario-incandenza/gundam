extern crate gundam;
#[macro_use]
extern crate log;
extern crate env_logger;
extern crate chrono;
extern crate darwin_rs;

use std::env;
use std::process::exit;
use std::io::{BufReader, BufRead};
use std::fs::File;
use chrono::Local;
use env_logger::LogBuilder;
use gundam::*;
use darwin_rs::individual::Individual;


fn main() {
    //env_logger::init();
    let _ = LogBuilder::new()
        .format(|record| {
            format!(
                "{} [{}] - {}",
                Local::now().format("%Y-%m-%dT%H:%M:%S"),
                record.level(),
                record.args()
            )
        })
        .parse(&env::var("RUST_LOG").unwrap_or_default())
        .init()
        .expect("log init");

    let args = env::args().collect::<Vec<String>>();
    if args.len() != 4 {
        println!(
            "invalid # of args: {}.  usage: gundam <indices.txt> <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    let file = File::open(&args[1]).expect("can't open index file");
    let idx_file = BufReader::new(&file);
    let indices: Vec<(usize, usize, usize, f64)> = idx_file
        .lines()
        .map(|line| {
            let a = line.as_ref()
                .expect("no line?")
                .split(",")
                .collect::<Vec<&str>>();
            (
                a[0].parse::<usize>().expect("first"),
                a[1].parse::<usize>().expect("second"),
                a[2].parse::<usize>().expect("third"),
                a[3].parse::<f64>().expect("fourth"),
            )
        })
        .collect();

    info!("got {} indices", indices.len());

    for (idx, mut d) in find_motifs(indices, &args[2], &args[3]).into_iter().enumerate() {

        println!("{}: {}", idx, d.show_motif());
        println!("{}: {}", idx, d.calculate_fitness());
        println!("{}: {:?}", idx, d.history);
        println!("{}: {:?}", idx, d.motif.scores);
    }

}
