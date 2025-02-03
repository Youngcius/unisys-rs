use crate::basic::circuits::Circuit;
use crate::basic::operations::Operation;

// TODO: 将这些特性都迁移到 circuits.rs 中

pub fn circuit_to_dag(circ: &Circuit) {
    panic!("Not implemented");
}

pub fn front_layer(circ: &Circuit) -> Vec<Operation> {
    let mut front_layer: Vec<Operation> = Vec::new();
    // TODO:
    front_layer
}

pub fn last_layer(circ: &Circuit) -> Vec<Operation> {
    let mut last_layer: Vec<Operation> = Vec::new();

    last_layer
}

pub fn front_full_width_circuit(circ: &Circuit) -> Circuit {
    let mut ffwc = Circuit::new();

    ffwc
}

pub fn last_full_width_circuit(circ: &Circuit) -> Circuit {
    let mut lfwc = Circuit::new();

    lfwc
}
