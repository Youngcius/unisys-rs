
use super::gates::Gate;

#[derive(Clone, Debug)]
pub struct Operation {
    pub gate: Gate,
    pub tqs: Vec<usize>,
    pub cqs: Option<Vec<usize>>,
}


impl Operation {
    pub fn new(gate: Gate, tqs: Vec<usize>, cqs: Option<Vec<usize>>) -> Self {
        Operation { gate, tqs, cqs }
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
    
    pub fn hermitian(&self) -> Self {
        Operation::new(self.gate.hermitian(), self.tqs.clone(), self.cqs.clone())
    }
}


impl std::fmt::Display for Operation {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match &self.cqs {
            Some(cqs) => write!(f, "{}({:?} | {:?})", self.gate.name(), self.tqs, cqs),
            None => write!(f, "{}({:?})", self.gate.name(), self.tqs),
        }
    }
}



// write code to test Display trait for Operation
#[cfg(test)]
mod tests {
    use crate::basic::gates;

    use super::*;

    #[test]
    fn test_display() {
        let op = Operation::new(Gate::H(gates::HGate{}), vec![0], None);
        assert_eq!(format!("{}", op), "H([0])");
        println!("{:?}", op);
        println!("{}", op);
        println!("tq: {}", op.tq().unwrap());
        println!("cq: {:?}", op.cq());

        let op = Operation::new(Gate::X(gates::XGate{}), vec![0], Some(vec![1]));
        assert_eq!(format!("{}", op), "X([0] | [1])");
        println!("{:?}", op);
        println!("{}", op);
        println!("tq: {}", op.tq().unwrap());
        println!("cq: {:?}", op.cq().unwrap());

        // TODO: print parametrized gates
    }

    #[test]
    fn test_hermitian() {
        let op = Operation::new(Gate::S(gates::SGate{}), vec![0], None);
        let op_h = op.hermitian();
        assert_eq!(format!("{}", op_h), "SDG([0])");
        println!("{}", op_h);

        let op = Operation::new(Gate::X(gates::XGate{}), vec![0], Some(vec![1]));
        let op_h = op.hermitian();
        assert_eq!(format!("{}", op_h), "X([0] | [1])");
        println!("{}", op_h);

        // TODO: test parametrized gates
    }
}
