extern crate gundam;
#[macro_use]
extern crate log;
extern crate env_logger;

use std::env;
use std::process::exit;
use std::io::{BufReader, BufRead};
use std::fs::File;
use gundam::*;


fn main() {
    env_logger::init();
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
    let indices: Vec<(usize, usize, usize)> = idx_file
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
            )
        })
        .collect();

    info!("got {} indices", indices.len());

    for d in find_motifs(indices, &args[2], &args[3]) {
        println!("{}", d.show_motif());
    }

}
