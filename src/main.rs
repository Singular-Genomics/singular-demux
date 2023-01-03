#![allow(clippy::missing_errors_doc, clippy::missing_panics_doc)]
use std::process::exit;

use log::error;
use sgdemux_lib::opts::setup;
use sgdemux_lib::opts::DOC_LINK_AND_SUPPORT_EMAIL;
use sgdemux_lib::run::run;

#[global_allocator]
static GLOBAL: mimalloc::MiMalloc = mimalloc::MiMalloc;

#[cfg(not(tarpaulin_include))]
#[allow(clippy::too_many_lines)]
fn main() {
    let opts = setup();

    if let Err(err) = run(opts) {
        eprintln!("{}", DOC_LINK_AND_SUPPORT_EMAIL);
        error!("{:#}", err);
        exit(1);
    }
}
