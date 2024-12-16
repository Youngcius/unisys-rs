// Extended matrix traits for ndarray

use ndarray::{Array2, ArrayView2};
use ndarray_linalg::c64;


// conjugation of a complex matrix
pub trait Conj {
    fn conj(&self) -> Array2<c64>;
}

impl Conj for Array2<c64> {
    fn conj(&self) -> Array2<c64> {
        self.mapv(|i| i.conj())
    }
}

impl Conj for ArrayView2<'_, c64> {
    fn conj(&self) -> Array2<c64> {
        self.mapv(|i| i.conj())
    }
}

// dagger (hermitian) of a complex matrix
pub trait Dagger {
    fn dagger(&self) -> Array2<c64>;
}

impl Dagger for Array2<c64> {
    fn dagger(&self) -> Array2<c64> {
        self.mapv(|i| i.conj()).reversed_axes()
    }
}

impl<'a> Dagger for ArrayView2<'a, c64> {
    fn dagger(&self) -> Array2<c64> {
        self.mapv(|i| i.conj()).reversed_axes()
    }
}

// Kronecker product of two complex matrices

