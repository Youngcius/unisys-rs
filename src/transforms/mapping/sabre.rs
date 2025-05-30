// use std::collections::HashMap;
// use std::ptr::NonNull;
// use rustworkx_core::petgraph::graph::UnGraph;

// use crate::basic::{gates, gates::Gate, gates::GateType, circuits::Circuit};
// use crate::basic::{operations, operations::Operation};
// use crate::mapping;
// use crate::utils::{passes, arch};

// const INIT_DECAY: f64 = 1.0;
// const DECAY_STEP: f64 = 0.001;
// const NUM_SEARCHES_TO_RESET: usize = 5;
// const EXT_WEIGHT: f64 = 0.5;
// const EXT_SIZE: usize = 20;

// pub fn sabre(
//     circ: &Circuit,
//     device: & UnGraph<usize, ()>,
//     num_pass_periods: usize,
//     init_mapping: Option<&HashMap<usize, usize>>
// ) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
//     let (mapped_circ, init_mapping, final_mapping) = sabre_one_pass(circ, device, init_mapping);

//     // acquire the reverse of the circuit
//     let ops = &circ.ops;
//     let mut ops = ops.clone();
//     ops.reverse();
//     let circ_rev = Circuit { ops };

//     let mut results = vec![(mapped_circ, init_mapping.clone(), final_mapping.clone())];

//     // let mut mut_init_mapping = init_mapping.clone();
//     let mut mut_final_mapping = final_mapping.clone();
//     for _ in 0..num_pass_periods {
//         // init_mapping = final_mapping.clone();
//         let (_, _, final_mapping) = sabre_one_pass(&circ_rev, device, Some(&mut_final_mapping));
//         // init_mapping = final_mapping.clone();
//         let (mapped_circ, init_mapping, final_mapping) = sabre_one_pass(circ, device, Some(&final_mapping));
//         // mut_init_mapping = init_mapping.clone();
//         mut_final_mapping = final_mapping.clone();
//         results.push((mapped_circ.clone(), init_mapping.clone(), final_mapping.clone()));
//     }

//     let idx_min_num_swaps = results.iter().enumerate().min_by_key(|(_, res)| res.0.gate_count(GateType::Swap)).unwrap().0;
//     results[idx_min_num_swaps].clone()
// }

// pub fn sabre_one_pass(
//     circ: &Circuit,
//     device: &UnGraph<usize, ()>,
//     init_mapping: Option<&HashMap<usize, usize>>
// ) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
//     let mut dag = circ.dag();
//     let qubits = circ.qubits();
//     let mut mappings = Vec::new();
//     let mut front_layer = circ.front_layer();
//     let mut exe_gates = Vec::new();
//     let circ_with_swaps = Circuit::new();
//     let mut num_searches = 0;
//     let mut decay_params = HashMap::new();
//     for q in qubits {
//         decay_params.insert(q, INIT_DECAY);
//     }

//     panic!("Not implemented yet");

//     if let Some(init_map) = init_mapping {
//         mappings.push(init_map.clone());
//     } else {
//         mappings.push(arch::gene_init_mapping(circ, device, "random"));
//     }

//     while !front_layer.is_empty() {
//         exe_gates.clear();

//         for op in front_layer {
//             if arch::is_executable(op, &mappings[mappings.len() - 1], device) {
//                 exe_gates.push(op);
//             }
//         }

//         if !exe_gates.is_empty() {
//             for op in exe_gates {
//                 circ.ops.re
//             }
//         }
//     }

//     (arch::unify_mapped_circuit(&circ_with_swaps, &mappings), mappings[0].clone(), mappings[mappings.len() - 1].clone())

// }

//  #[cfg(test)]
//  mod tests {
//     use super::*;

//  }
