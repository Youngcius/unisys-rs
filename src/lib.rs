pub mod basic;
pub mod decompose;
pub mod mapping;
pub mod models;
pub mod phoenix;
pub mod transforms;
pub mod utils;

#[macro_export]
macro_rules! c {
    ($re:expr, $im:expr) => {
        c64::new($re, $im)
    };
}

#[macro_export]
macro_rules! r {
    ($re:expr) => {
        c64::new($re, 0.0)
    };
}

#[macro_export]
macro_rules! i {
    ($im:expr) => {
        c64::new(0.0, $im)
    };
}
