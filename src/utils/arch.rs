use crate::basic::circuits::Circuit;
use crate::basic::gates::{Gate, GateType};
use crate::basic::operations::Operation;
use ndarray::Array2;
use ndarray_linalg::c64;
use ndarray_rand::rand::distributions::DistMap;
use rand::seq::SliceRandom;
// use rustworkx_core::connectivity;
// use rustworkx_core::petgraph::graph::NodeIndex;
// use rustworkx_core::petgraph::prelude::StableUnGraph;
use std::collections::HashMap;

// pub fn is_executable(op: &Operation, mapping: &HashMap<usize, usize>, dist_mat: &Array2<usize>) -> bool {
//     let qregs = op.qregs();
//     if qregs.len() == 1 {
//         return true;
//     }
//     if qregs.len() > 2 {
//         panic!("Only 1 or 2 qubit gates are supported")
//     }
//     dist_mat[[mapping[&qregs[0]], mapping[&qregs[1]]]] == 1
// }

// pub fn is_executable(
//     op: &Operation,
//     mapping: &HashMap<usize, usize>,
//     device: &StableUnGraph<usize, ()>,
// ) -> bool {
//     let qregs = op.qregs();
//     if qregs.len() == 1 {
//         return true;
//     }
//     if qregs.len() > 2 {
//         panic!("Only 1 or 2 qubit gates are supported")
//     }
//     device.contains_edge(
//         NodeIndex::new(mapping[&qregs[0]]),
//         NodeIndex::new(mapping[&qregs[1]]),
//     )
// }

// pub fn gene_trivial_init_mapping(
//     circ: &Circuit,
//     device: &StableUnGraph<usize, ()>,
// ) -> HashMap<usize, usize> {
//     let n = circ.num_qubits();
//     let logic_qubits = circ.qubits();

//     let mut subgraph = device.clone();
//     subgraph.retain_nodes(|_, idx| idx.index() < n);
//     let connected = connectivity::connected_components(&subgraph).len() == 1;
//     if !connected {
//         panic!("The subgraph of the first n physical-qubits is not connected")
//     }
//     let phys_qubits = subgraph
//         .node_indices()
//         .map(|idx| idx.index())
//         .collect::<Vec<_>>();

//     // let mut mapping = HashMap::new();
//     // for (i, q) in logic_qubits.iter().enumerate() {
//     //     mapping.insert(*q, phys_qubits[i]);
//     // }
//     let mapping = logic_qubits
//         .iter()
//         .zip(phys_qubits.iter())
//         .map(|(l, p)| (*l, *p))
//         .collect();

//     mapping
// }

// pub fn gene_random_init_mapping(
//     circ: &Circuit,
//     device: &StableUnGraph<usize, ()>,
// ) -> HashMap<usize, usize> {
//     let n = circ.num_qubits();
//     let logic_qubits = circ.qubits();
//     let mut subgraph = device.clone();
//     let mut indices = subgraph.node_indices().collect::<Vec<_>>();
//     indices.shuffle(&mut rand::rng());
//     indices.truncate(n);
//     subgraph.retain_nodes(|_, idx| indices.contains(&idx));
//     let mut connected = connectivity::connected_components(&subgraph).len() == 1;
//     while !connected {
//         subgraph = device.clone();
//         indices = subgraph.node_indices().collect::<Vec<_>>();
//         indices.shuffle(&mut rand::rng());
//         indices.truncate(n);
//         subgraph.retain_nodes(|_, idx| indices.contains(&idx));
//         connected = connectivity::connected_components(&subgraph).len() == 1;
//     }
//     let phys_qubits = subgraph
//         .node_indices()
//         .map(|idx| idx.index())
//         .collect::<Vec<_>>();
//     let mapping = logic_qubits
//         .iter()
//         .zip(phys_qubits.iter())
//         .map(|(l, p)| (*l, *p))
//         .collect();
//     mapping
// }

pub fn update_mapping(mapping: HashMap<usize, usize>, swap: Operation) -> HashMap<usize, usize> {
    if swap.gate.gate_type != GateType::Swap {
        panic!("Input Operation should be a Swap gate")
    }
    let mut updated_mapping = mapping.clone();
    updated_mapping.insert(swap.tqs[0], mapping[&swap.tqs[1]]);
    updated_mapping.insert(swap.tqs[1], mapping[&swap.tqs[0]]);
    updated_mapping
}

pub fn unify_mapped_circuit(
    circ_with_swaps: &Circuit,
    mappings: &Vec<HashMap<usize, usize>>,
) -> Circuit {
    let mut mapped_circ = Circuit::new();
    let mut mapping_idx = 0;
    let mut mapping = &mappings[mapping_idx];
    for op in &circ_with_swaps.ops {
        if op.gate.gate_type == GateType::Swap {
            mapped_circ.append(op.reapply(vec![mapping[&op.tqs[0]], mapping[&op.tqs[1]]], None));
            mapping_idx += 1;
            mapping = &mappings[mapping_idx];
        } else {
            mapped_circ.append(
                op.reapply(
                    op.tqs.iter().map(|tq| mapping[tq]).collect(),
                    op.cqs
                        .as_ref()
                        .map(|cqs| cqs.iter().map(|cq| mapping[cq]).collect()),
                ),
            );
        }
    }
    let inverse_mapping: HashMap<usize, usize> =
        mappings[0].iter().map(|(k, v)| (*v, *k)).collect();
    mapped_circ.rewire(&inverse_mapping);
    mapped_circ
}
