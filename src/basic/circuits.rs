use super::gates::GateType;
use super::operations::Operation;
use crate::utils::ops;
use crate::utils::passes;
use core::panic;
use ndarray::Array2;
use ndarray_linalg::c64;
use rustworkx_core::petgraph::prelude::StableDiGraph;
use std::collections::{HashMap, HashSet};

#[derive(Clone, Debug)]
pub struct Circuit {
    pub ops: Vec<Operation>,
}

impl Circuit {
    pub fn new() -> Self {
        Circuit { ops: Vec::new() }
    }

    pub fn append(&mut self, op: Operation) {
        self.ops.push(op);
    }

    pub fn prepend(&mut self, op: Operation) {
        self.ops.insert(0, op);
    }

    pub fn compose(&mut self, other: &Circuit) {
        self.ops.extend(other.ops.clone());
    }

    pub fn num_ops(&self) -> usize {
        self.ops.len()
    }

    pub fn num_nonlocal_ops(&self) -> usize {
        self.ops.iter().filter(|op| op.qregs().len() > 1).count()
    }

    pub fn max_op_weight(&self) -> usize {
        let mut max_weight = 0;
        for op in &self.ops {
            if op.qregs().len() > max_weight {
                max_weight = op.qregs().len();
            }
        }
        max_weight
    }

    pub fn qubits(&self) -> Vec<usize> {
        let mut qubits = HashSet::new();
        for op in &self.ops {
            for tq in &op.tqs {
                qubits.insert(*tq);
            }
            if let Some(cqs) = &op.cqs {
                for cq in cqs {
                    qubits.insert(*cq);
                }
            }
        }
        let mut qubits: Vec<usize> = qubits.into_iter().collect();
        qubits.sort();
        qubits
    }

    pub fn num_qubits(&self) -> usize {
        self.qubits().len()
    }

    pub fn num_qubits_with_dummy(&self) -> usize {
        self.qubits().iter().max().unwrap() + 1
    }

    pub fn gate_count(&self, gate_type: GateType) -> usize {
        self.ops.iter().filter(|op| op.gate.gate_type == gate_type).count()
    }

    pub fn gate_stats(&self) -> HashMap<String, usize> {
        // TODO: review this implementation
        let mut stats = HashMap::new();
        for op in &self.ops {
            let gate_name = op.gate.to_string();
            let count = stats.entry(gate_name).or_insert(0);
            *count += 1;
        }
        stats
    }

    pub fn unitary(&self) -> Array2<c64> {
        let u = Array2::<c64>::eye(1 << self.num_qubits());
        self.ops.iter().fold(u, |acc, op| {
            ops::tensor_slots(&op.matrix(), self.num_qubits(), &op.qregs()).dot(&acc)
        })
    }

    pub fn rewire(&self, mapping: &HashMap<usize, usize>) -> Self {
        let mut rewired_circ = Circuit::new();
        for op in &self.ops {
            let tqs = op.tqs.iter().map(|tq| mapping[tq]).collect();
            let cqs = op
                .cqs
                .as_ref()
                .map(|cqs| cqs.iter().map(|cq| mapping[cq]).collect());
            rewired_circ.append(op.reapply(tqs, cqs))
        }
        rewired_circ
    }

    pub fn inverse(&self) -> Self {
        let mut ops = Vec::new();
        for op in self.ops.iter().rev() {
            ops.push(op.hermitian());
        }
        Circuit { ops }
    }

    pub fn dag(&self) -> StableDiGraph<Operation, ()> {
        passes::circuit_to_dag(self);
    }

    pub fn qasm(&self) -> String {
        let mut qasm_str = String::new();

        // add header
        qasm_str.push_str(format!("// Date: {}\n", chrono::Utc::now()).as_str());
        qasm_str.push_str("OPENQASM 2.0;\n");
        qasm_str.push_str("include \"qelib1.inc\";\n\n");

        // add qubits
        qasm_str.push_str(format!("qreg q[{}];\n\n", self.num_qubits()).as_str());

        // add special gate definitions
        let defined_gates: HashSet<_> = self
            .ops
            .iter()
            .filter_map(|op| op.gate.qasm_def())
            .collect();
        for gate_def in defined_gates {
            qasm_str.push_str(&format!("{}\n", gate_def));
        }

        // add operations
        for op in &self.ops {
            qasm_str.push_str(&format!("{};\n", op.qasm_cmd()));
        }

        qasm_str
    }

    pub fn depth(&self) -> usize {
        let mut wire_lengths = vec![0; self.num_qubits()];
        for op in &self.ops {
            let qubits = op.qregs();
            let current_depth = qubits.iter().map(|q| wire_lengths[*q]).max().unwrap() + 1;
            for q in qubits {
                wire_lengths[q] = current_depth;
            }
        }
        wire_lengths.into_iter().max().unwrap()
    }

    pub fn depth_nonlocal(&self) -> usize {
        let mut wire_lengths = vec![0; self.num_qubits()];
        for op in &self.ops {
            if op.qregs().len() == 1 {
                continue;
            }
            let qubits = op.qregs();
            let current_depth = qubits.iter().map(|q| wire_lengths[*q]).max().unwrap() + 1;
            for q in qubits {
                wire_lengths[q] = current_depth;
            }
        }
        wire_lengths.into_iter().max().unwrap()
    }

    pub fn layer(&self) -> Vec<Vec<&Operation>> {
        panic!("Not implemented yet")
    }

    pub fn front_layer(&self) -> Vec<&Operation> {
        passes::front_layer(self)
    }

    pub fn last_layer(&self) -> Vec<&Operation> {
        passes::last_layer(self)
    }

    pub fn front_full_width_circuit(&self) -> Circuit {
        passes::front_full_width_circuit(self)
    }

    pub fn last_full_width_circuit(&self) -> Circuit {
        passes::last_full_width_circuit(self)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::basic::gates::Gate;
    use crate::basic::matrices::{Imag, Real};

    #[test]
    fn test_demo_circuit() {
        // h q[0];
        // h q[2];
        // h q[5];
        // z q[0];
        // x q[1] q[2];
        // x q[4] q[5];
        // x q[0] q[1];
        // x q[2] q[3];
        // h q[2];
        // h q[3];
        // x q[1] q[2];
        // x q[3] q[5];
        // z q[3];
        // x q[4] q[3];
        // x q[3] q[0];

        let mut circ = Circuit::new();

        circ.append(Operation::new(Gate::h(), vec![0], None));
        circ.append(Operation::new(Gate::h(), vec![2], None));
        circ.append(Operation::new(Gate::h(), vec![5], None));
        circ.append(Operation::new(Gate::z(), vec![0], None));
        circ.append(Operation::new(Gate::x(), vec![2], Some(vec![1])));
        circ.append(Operation::new(Gate::x(), vec![5], Some(vec![4])));
        circ.append(Operation::new(Gate::x(), vec![1], Some(vec![0])));
        circ.append(Operation::new(Gate::x(), vec![3], Some(vec![2])));
        circ.append(Operation::new(Gate::h(), vec![2], None));
        circ.append(Operation::new(Gate::h(), vec![3], None));
        circ.append(Operation::new(Gate::x(), vec![2], Some(vec![1])));
        circ.append(Operation::new(Gate::x(), vec![5], Some(vec![3])));
        circ.append(Operation::new(Gate::z(), vec![3], None));
        circ.append(Operation::new(Gate::x(), vec![3], Some(vec![4])));
        circ.append(Operation::new(Gate::x(), vec![0], Some(vec![3])));

        println!("{}", circ.qasm());

        println!(
            "The first row of this unitary's real part is {}",
            circ.unitary().real().row(0)
        );
        println!(
            "The first row of this unitary's imaginary part is {}",
            circ.unitary().imag().row(0)
        );
    }
}
