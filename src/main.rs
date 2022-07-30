#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]
use std::process::exit;

use log::error;
use singular_demux_lib::run::{run, setup};

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(not(tarpaulin_include))]
#[allow(clippy::too_many_lines)]
fn main() {
    let opts = setup();

    if let Err(err) = run(opts) {
        error!("{:#}", err);
        exit(1);
    }
}
