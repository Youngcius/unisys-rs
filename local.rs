//! Dependencies can be specified in the script file itself as follows:
//!
//! ```cargo
//! [dependencies]
//! rust-script-ext = "0.1.0"
//! ```
//! 
use rust_script_ext::prelude::*;

fn main() -> Result<()> {
    std::fs::read("foo.txt")
        .into_diagnostic()
        .wrap_err("failed to read file `foo.txt`")?;
    Ok(())
}
