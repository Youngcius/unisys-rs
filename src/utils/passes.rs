use std::collections::{HashMap, HashSet};
use rustworkx_core::petgraph::graph::NodeIndex;
use rustworkx_core::petgraph::prelude::StableDiGraph;

use crate::basic::circuits::Circuit;
use crate::basic::operations::Operation;

// TODO: 将这些特性都迁移到 circuits.rs 中


/*
    def to_dag(self, backend='rustworkx') -> Union[rx.PyDiGraph, nx.DiGraph]:
        """
        Convert a circuit into a Directed Acyclic Graph (DAG) according to dependency of each gate's qubits.

        Args:
            backend: 'networkx' or 'rustworkx'
        """
        all_gates = self.gates
        dag = rx.PyDiGraph(multigraph=False)
        dag.add_nodes_from(all_gates)
        gate_to_node_idx = {dag[idx]: idx for idx in dag.node_indices()}
        while all_gates:
            g = all_gates.pop(0)
            qregs = set(g.qregs)
            for g_opt in all_gates:  # traverse the subsequent optional gates
                qregs_opt = set(g_opt.qregs)
                if dependent_qubits := qregs_opt & qregs:
                    dag.add_edge(gate_to_node_idx[g], gate_to_node_idx[g_opt], {'qubits': list(dependent_qubits)})
                    qregs -= qregs_opt
                if not qregs:
                    break
        if backend == 'networkx':
            return rx_to_nx_graph(dag)
        return dag
 */

pub fn circuit_to_dag(circ: &Circuit) -> StableDiGraph<Operation, ()> {
    let mut dag = StableDiGraph::new();
    // let mut gate_to_node_idx: HashMap<&Operation, &NodeIndex> = HashMap::new();
    for (idx, op) in circ.ops.iter().enumerate() {
        let node_idx = dag.add_node(op.clone());
        // gate_to_node_idx.insert(node_idx, idx);
    }
    // let w = dag.node_weight(NodeIndex::new(0)).unwrap();
    // let gate_to_node_idx:HashMap<&Operation, &NodeIndex> = dag.node_indices().map(|idx| {
    //     let w = dag.node_weight(idx).unwrap();
    //     (w, &idx)
    // }).collect();

    let mut all_gates = circ.ops.clone();
    while let Some(g) = all_gates.pop() { // TODO: should pop the first element
        let qregs = g.qregs().iter().cloned().collect::<HashSet<_>>();
        for g_opt in &all_gates {
            let qregs_opt = g_opt.qregs().iter().cloned().collect::<HashSet<_>>();
            let dependent_qubits = qregs.intersection(&qregs_opt).cloned().collect::<HashSet<_>>();
            if !dependent_qubits.is_empty() {
                dag.add_edge(gate_to_node_idx[&g], gate_to_node_idx[&g_opt], ());
                qregs.difference(&qregs_opt);
            }
            if qregs.is_empty() {
                break;
            }
        }
    }
    dag
}

pub fn front_layer(circ: &Circuit) -> Vec<&Operation> {
    let mut front_layer: Vec<&Operation> = Vec::new();
    let mut visited_qubits = HashSet::new();
    let n = circ.num_qubits();
    for op in &circ.ops {
        if !op.qregs().iter().any(|q| visited_qubits.contains(q)) {
            front_layer.push(op);
        }
        visited_qubits.extend(op.qregs());
        if visited_qubits.len() == n {
            break;
        }
    }
    front_layer
}

pub fn last_layer(circ: &Circuit) -> Vec<&Operation> {
    let mut last_layer: Vec<&Operation> = Vec::new();
    let mut visited_qubits = HashSet::new();
    let n = circ.num_qubits();
    for op in circ.ops.iter().rev() {
        if !op.qregs().iter().any(|q| visited_qubits.contains(q)) {
            last_layer.push(op);
        }
        visited_qubits.extend(op.qregs());
        if visited_qubits.len() == n {
            break;
        }
    }
    last_layer
}


pub fn front_full_width_circuit(circ: &Circuit) -> Circuit {
    let n = circ.num_qubits();
    let mut ffwc = Circuit::new();
    let mut ffwc_qubits = HashSet::new();
    for op in &circ.ops {
        ffwc.append(op.clone());
        ffwc_qubits.extend(op.qregs());
        if ffwc_qubits.len() == n {
            break;
        }
    }
    ffwc
}

pub fn last_full_width_circuit(circ: &Circuit) -> Circuit {
    let n = circ.num_qubits();
    let mut lfwc = Circuit::new();
    let mut lfwc_qubits = HashSet::new();
    for op in circ.ops.iter().rev() {
        lfwc.append(op.clone());
        lfwc_qubits.extend(op.qregs());
        if lfwc_qubits.len() == n {
            break;
        }
    }
    lfwc.ops.reverse();
    lfwc
}


pub fn mirror_near_identity(
    circ: &Circuit,
    init_mapping: &HashMap<usize, usize>,
    kak_norm_thresh: Option<f64>,
) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
    // if kak_norm_thresh is None, set it to 1e-10
    let kak_norm_thresh = kak_norm_thresh.unwrap_or(0.3);


    panic!("Not implemented")
}
