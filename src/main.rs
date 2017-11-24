extern crate gundam;

use std::env;
use std::process::exit;
use gundam::*;


fn main() {
    let args = env::args().collect::<Vec<String>>();
    if args.len() != 3 {
        println!(
            "invalid # of args: {}.  usage: gundam <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    find_motifs(args[1].as_str(), args[2].as_str(), 16);
}
