use crate::basic::gates::{Gate, GateType};
use crate::basic::operations::Operation;
use ndarray::Array2;
use ndarray_linalg::c64;
use ndarray_rand::rand::distributions::DistMap;
use std::collections::HashMap;

fn is_executable(op: Operation, mapping: &HashMap<usize, usize>, dist_mat: Array2<usize>) -> bool {
    let qregs = op.qregs();
    if qregs.len() == 1 {
        return true;
    }
    if qregs.len() > 2 {
        panic!("Only 1 or 2 qubit gates are supported")
    }
    dist_mat[[mapping[&qregs[0]], mapping[&qregs[1]]]] == 1
}

fn update_mapping(mapping: HashMap<usize, usize>, swap: Operation) -> HashMap<usize, usize> {
    if swap.gate.gate_type != GateType::Swap {
        panic!("Input Operation should be a Swap gate")
    }
    let mut updated_mapping = mapping.clone();
    updated_mapping.insert(swap.tqs[0], mapping[&swap.tqs[1]]);
    updated_mapping.insert(swap.tqs[1], mapping[&swap.tqs[0]]);
    updated_mapping
}

// fn is_executable_on_coupling_map(op: Operation, mapping: &HashMap<usize, usize>, device: ...) -> bool {
//     // TODO: device -> dist_mat
// return device.has_edge(mapping[gate.qregs[0]], mapping[gate.qregs[1]])

//     // is_executable(op, mapping, dist_mat)
// }

// def is_executable(gate: Gate, mapping: Dict[int, int], dist_mat: ndarray = None, device: rx.PyGraph = None):
//     """
