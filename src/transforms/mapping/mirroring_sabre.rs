use std::collections::HashMap;

use crate::basic::{circuits::Circuit, gates::Gate, operations::Operation};


pub fn mirroring_sabre(
    circ: Circuit,
    device: ...,
    num_pass_periods: usize,
    init_mapping: &HashMap<usize, usize>
) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
    panic!("Not implemented")
}

pub fn mirroring_sabre_one_pass(
    circ: Circuit,
    device: ...,
    init_mapping: &HashMap<usize, usize>
) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
    panic!("Not implemented")
}
