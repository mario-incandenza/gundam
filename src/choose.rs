extern crate gundam;
#[macro_use]
extern crate log;
extern crate env_logger;

use std::env;
use std::process::exit;
use gundam::*;

fn main() {
    let _ = env_logger::init();
    let args = env::args().collect::<Vec<String>>();
    if args.len() != 3 {
        println!(
            "invalid # of args: {}.  usage: gundam <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    for (i, j, k, f) in DyadMotif::passing_kmers(args[1].as_str(), args[2].as_str()) {
        println!("{},{},{},{:e}", i, j, k, f);
    }
}
