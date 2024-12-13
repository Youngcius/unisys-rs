use super::operations::Operation;

pub struct Circuit {
    pub ops: Vec<Operation>,
}

impl Circuit {
    pub fn new() -> Self {
        Circuit { ops: Vec::new() }
    }

    // pub fn add_op(&mut self, op: Operation) {
    //     self.ops.push(op);
    // }

    // pub fn n_qubits(&self) -> usize {
    //     let mut n_qubits = 0;
    //     for op in &self.ops {
    //         for tq in &op.tqs {
    //             if *tq >= n_qubits {
    //                 n_qubits = *tq + 1;
    //             }
    //         }
    //         if let Some(cqs) = &op.cqs {
    //             for cq in cqs {
    //                 if *cq >= n_qubits {
    //                     n_qubits = *cq + 1;
    //                 }
    //             }
    //         }
    //     }
    //     n_qubits
    // }
}
