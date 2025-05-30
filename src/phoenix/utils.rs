use crate::basic::circuits::Circuit;
use crate::basic::operations::Operation;
use crate::models::paulis::BSF;

#[derive(Debug, Clone)]
pub enum SimpItem {
    BSFItem(BSF),
    UCGItem(Operation),
}

impl std::fmt::Display for SimpItem {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            SimpItem::BSFItem(bsf) => write!(f, "{:?}", bsf.paulis()),
            SimpItem::UCGItem(ucg) => write!(f, "{}", ucg),
        }
    }
}

pub fn config_to_circuit(config: Vec<SimpItem>) -> Circuit {
    // TODO: make this function parallel
    let mut circ = Circuit::new();
    for item in config {
        match item {
            SimpItem::BSFItem(bsf) => {
                circ.compose(&bsf.as_cnot_circuit());
            }
            SimpItem::UCGItem(ucg) => {
                circ.append(ucg);
            }
        }
    }
    circ
}
