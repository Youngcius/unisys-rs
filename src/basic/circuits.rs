use super::operations::Operation;
use ndarray_linalg::c64;
use std::collections::HashSet;

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

    pub fn num_ops(&self) -> usize {
        self.ops.len()
    }

    pub fn max_ops_weight(&self) -> usize {
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

    pub fn inverse(&self) -> Self {
        let mut ops = Vec::new();
        for op in self.ops.iter().rev() {
            ops.push(op.hermitian());
        }
        Circuit { ops }
    }
}
