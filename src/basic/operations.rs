use super::gates::Gate;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub struct Operation {
    pub gate: Gate,
    pub tqs: Vec<usize>,
    pub cqs: Option<Vec<usize>>,
}

impl Operation {
    pub fn new(gate: Gate, tqs: Vec<usize>, cqs: Option<Vec<usize>>) -> Self {
        if tqs.is_empty() {
            panic!("Operation must designate target qubit(s)");
        }
        if tqs.len() > 1 && tqs.len() != gate.n_qubits {
            panic!(
                "{} must have act on {} target qubits",
                gate.to_string(),
                gate.n_qubits
            );
        }
        if let Some(cqs) = &cqs {
            for tq in &tqs {
                if cqs.contains(tq) {
                    panic!("Target qubit(s) and control qubit(s) should not be overlapped");
                }
            }
        }
        Operation {
            gate: gate,
            tqs,
            cqs,
        }
    }

    pub fn tq(&self) -> Result<usize, &str> {
        if self.tqs.len() != 1 {
            return Err("Operation must have exactly one target qubit to return tq");
        }
        Ok(self.tqs[0])
    }

    pub fn cq(&self) -> Result<usize, &str> {
        if let Some(cqs) = &self.cqs {
            if cqs.len() != 1 {
                return Err("Operation must have exactly one control qubit to return cq");
            }
            Ok(cqs[0])
        } else {
            Err("Operation does not have control qubit")
        }
    }

    pub fn qregs(&self) -> Vec<usize> {
        let mut qregs = self.tqs.clone();
        if let Some(cqs) = &self.cqs {
            qregs.extend(cqs.clone());
        }
        qregs
    }

    pub fn hermitian(&self) -> Self {
        Operation::new(self.gate.hermitian(), self.tqs.clone(), self.cqs.clone())
    }
}

impl std::fmt::Display for Operation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Format target qubits
        let tqs_str = if self.tqs.len() == 1 {
            self.tqs[0].to_string()
        } else {
            self.tqs
                .iter()
                .map(|tq| tq.to_string())
                .collect::<Vec<_>>()
                .join("|")
        };

        let qregs_str = match &self.cqs {
            Some(cqs) => {
                // Format control qubits
                let cqs_str = if cqs.len() == 1 {
                    cqs[0].to_string()
                } else {
                    cqs.iter()
                        .map(|cq| cq.to_string())
                        .collect::<Vec<_>>()
                        .join("|")
                };
                format!("{}←{}", tqs_str, cqs_str)
            }
            None => format!("{}", tqs_str),
        };

        // Format angle or angles
        let prefix = self.gate.to_string();
        // if let Ok(angles) = self.gate.angles() {
        //     let angles_str = angles
        //         .iter()
        //         .map(|a| format!("{:.2}π", a / PI))
        //         .collect::<Vec<_>>()
        //         .join(",");
        //     prefix += &format!("({})", angles_str);
        // }

        // if let Ok(angle) = self.gate.angle() {
        //     prefix += &format!("({:.2}π)", angle / PI);
        // }

        write!(f, "{}({})", prefix, qregs_str)
    }
}

// write code to test Display trait for Operation
#[cfg(test)]
mod tests {
    use crate::basic::gates;

    use super::*;

    #[test]
    fn test_display() {
        let op = Operation::new(Gate::h(), vec![0], None);
        assert_eq!(format!("{}", op), "H(0)");
        println!("{:?}", op);
        println!("{}", op);
        println!("tq: {}", op.tq().unwrap());
        println!("cq: {:?}", op.cq());
        println!();

        let op = Operation::new(Gate::x(), vec![0], Some(vec![1]));
        assert_eq!(format!("{}", op), "X(0←1)");
        println!("{:?}", op);
        println!("{}", op);
        println!("tq: {}", op.tq().unwrap());
        println!("cq: {:?}", op.cq().unwrap());
        println!();

        let op = Operation::new(Gate::rx(1.12313), vec![0], None);
        println!("{:?}", op);
        println!("{}", op);
        println!("tq: {}", op.tq().unwrap());
        println!("cq: {:?}", op.cq());

        let op = Operation::new(Gate::can(1.1, 2.2, 3.3), vec![0, 1], vec![2, 3, 4].into());
        println!("{:?}", op);
        println!("{}", op);
    }

    #[test]
    fn test_hermitian() {
        let op = Operation::new(Gate::s(), vec![0], None);
        let op_h = op.hermitian();
        println!("{}, {}", op, op_h);

        let op = Operation::new(Gate::x(), vec![0], Some(vec![1]));
        let op_h = op.hermitian();
        println!("{}, {}", op, op_h);

        let op = Operation::new(Gate::rx(1.12313), vec![0], None);
        let op_h = op.hermitian();
        println!("{}, {}", op, op_h);
    }
}
