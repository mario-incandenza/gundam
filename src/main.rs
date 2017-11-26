extern crate gundam;
#[macro_use]
extern crate log;
extern crate env_logger;

use std::env;
use std::process::exit;
use gundam::*;

fn main() {
    env_logger::init();
    error!("start");
    let args = env::args().collect::<Vec<String>>();
    if args.len() != 3 {
        println!(
            "invalid # of args: {}.  usage: gundam <pos.fa> <neg.fa>",
            args.len()
        );
        exit(1);
    }

    info!("hello world");
    for d in find_motifs(args[1].as_str(), args[2].as_str(), 16) {
        println!("{}", d.show_motif());
    }
}