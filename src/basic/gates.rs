use std::f64::consts::PI;
use ndarray_linalg::c64;
use ndarray::{array, Array2};

#[derive(Clone, Debug)]
pub enum Gate {
    // Constant(ConstantGate),
    H(HGate),
    T(TGate),
    TDG(TDGGate),
    S(SGate),
    SDG(SDGGate),
    X(XGate),
    Y(YGate),
    Z(ZGate),
    SWAP(SWAPGate),
    CX(CXGate),
    CZ(CZGate),
    U1(U1Gate),
    U2(U2Gate),
    U3(U3Gate),
    RX(RXGate),
    RY(RYGate),
    RZ(RZGate),
    RXX(RXXGate),
    RYY(RYYGate),
    RZZ(RZZGate),
    Can(CanonicalGate)
}

#[derive(Clone, Debug)]
pub struct HGate;

#[derive(Clone, Debug)]
pub struct TGate;

#[derive(Clone, Debug)]
pub struct TDGGate;

#[derive(Clone, Debug)]
pub struct SGate;

#[derive(Clone, Debug)]
pub struct SDGGate {}

#[derive(Clone, Debug)]
pub struct XGate {}

#[derive(Clone, Debug)]
pub struct YGate {}

#[derive(Clone, Debug)]
pub struct ZGate {}

#[derive(Clone, Debug)]
pub struct SWAPGate {}

#[derive(Clone, Debug)]
pub struct CXGate {}

#[derive(Clone, Debug)]
pub struct CZGate {}

#[derive(Clone, Debug)]
pub struct U1Gate {
    pub lambda: f64,
}

#[derive(Clone, Debug)]
pub struct U2Gate {
    pub phi: f64,
    pub lambda: f64,
}

#[derive(Clone, Debug)]
pub struct U3Gate {
    pub theta: f64,
    pub phi: f64,
    pub lambda: f64,
}

#[derive(Clone, Debug)]
pub struct RXGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct RYGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct RZGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct RXXGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct RYYGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct RZZGate {
    pub theta: f64,
}

#[derive(Clone, Debug)]
pub struct CanonicalGate {
    pub theta1: f64,
    pub theta2: f64,
    pub theta3: f64,
}


impl Gate {
    pub fn name(&self) -> &str {
        match self {
            Gate::H(_) => "H",
            Gate::T(_) => "T",
            Gate::TDG(_) => "TDG",
            Gate::S(_) => "S",
            Gate::SDG(_) => "SDG",
            Gate::X(_) => "X",
            Gate::Y(_) => "Y",
            Gate::Z(_) => "Z",
            Gate::SWAP(_) => "SWAP",
            Gate::CX(_) => "CX",
            Gate::CZ(_) => "CZ",
            Gate::U1(_) => "U1",
            Gate::U2(_) => "U2",
            Gate::U3(_) => "U3",
            Gate::RX(_) => "RX",
            Gate::RY(_) => "RY",
            Gate::RZ(_) => "RZ",
            Gate::RXX(_) => "RXX",
            Gate::RYY(_) => "RYY",
            Gate::RZZ(_) => "RZZ",
            Gate::Can(_) => "Can",
        }
    }

    pub fn n_qubits(&self) -> usize {
        match self {
            Gate::H(_) => 1,
            Gate::T(_) => 1,
            Gate::TDG(_) => 1,
            Gate::S(_) => 1,
            Gate::SDG(_) => 1,
            Gate::X(_) => 1,
            Gate::Y(_) => 1,
            Gate::Z(_) => 1,
            Gate::SWAP(_) => 2,
            Gate::CX(_) => 2,
            Gate::CZ(_) => 2,
            Gate::U1(_) => 1,
            Gate::U2(_) => 1,
            Gate::U3(_) => 1,
            Gate::RX(_) => 1,
            Gate::RY(_) => 1,
            Gate::RZ(_) => 1,
            Gate::RXX(_) => 2,
            Gate::RYY(_) => 2,
            Gate::RZZ(_) => 2,
            Gate::Can(_) => 2,
        }
    }

    pub fn angle(&self) -> Result<f64, &str> {
        match self {
            Gate::H(_) => Err("H gate has no angle"),
            Gate::T(_) =>Err("T gate has no angle"),
            Gate::TDG(_) => Err("TDG gate has no angle"),
            Gate::S(_) => Err("S gate has no angle"),
            Gate::SDG(_) => Err("SDG gate has no angle"),
            Gate::X(_) => Err("X gate has no angle"),
            Gate::Y(_) => Err("Y gate has no angle"),
            Gate::Z(_) => Err("Z gate has no angle"),
            Gate::SWAP(_) => Err("SWAP gate has no angle"),
            Gate::CX(_) => Err("CX gate has no angle"),
            Gate::CZ(_) => Err("CZ gate has no angle"),
            Gate::U1(gate) => Ok(gate.lambda),
            Gate::U2(_) => Err("U2 gate has multiple angles"),
            Gate::U3(_) => Err("U3 gate has multiple angles"),
            Gate::RX(gate) => Ok(gate.theta),
            Gate::RY(gate) => Ok(gate.theta),
            Gate::RZ(gate) => Ok(gate.theta),
            Gate::RXX(gate) => Ok(gate.theta),
            Gate::RYY(gate) => Ok(gate.theta),
            Gate::RZZ(gate) => Ok(gate.theta),
            Gate::Can(_) => Err("Can gate has multiple angles"),
        }
    }

    pub fn angles(&self) -> Result<Vec<f64>, &str> {
        match self {
            Gate::H(_) => Err("H gate has no angles"),
            Gate::T(_) => Err("T gate has no angles"),
            Gate::TDG(_) => Err("TDG gate has no angles"),
            Gate::S(_) => Err("S gate has no angles"),
            Gate::SDG(_) => Err("SDG gate has no angles"),
            Gate::X(_) => Err("X gate has no angles"),
            Gate::Y(_) => Err("Y gate has no angles"),
            Gate::Z(_) => Err("Z gate has no angles"),
            Gate::SWAP(_) => Err("SWAP gate has no angles"),
            Gate::CX(_) => Err("CX gate has no angles"),
            Gate::CZ(_) => Err("CZ gate has no angles"),
            Gate::U1(_) => Err("U1 gate has no angles"),
            Gate::U2(gate) => Ok(vec![gate.phi, gate.lambda]),
            Gate::U3(gate) => Ok(vec![gate.theta, gate.phi, gate.lambda]),
            Gate::RX(_) => Err("RX gate has no angles"),
            Gate::RY(_) => Err("RY gate has no angles"),
            Gate::RZ(_) => Err("RZ gate has no angles"),
            Gate::RXX(_) => Err("RXX gate has no angles"),
            Gate::RYY(_) => Err("RYY gate has no angles"),
            Gate::RZZ(_) => Err("RZZ gate has no angles"),
            Gate::Can(gate) => Ok(vec![gate.theta1, gate.theta2, gate.theta3]),
        }
    }

    pub fn data(&self) -> Array2<c64> {
        match self {
            Gate::H(_) => array![
                [c64::new(1.0, 0.0), c64::new(1.0, 0.0)],
                [c64::new(1.0, 0.0), c64::new(-1.0, 0.0)]
            ] / 2.0_f64.sqrt(),
            Gate::T(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(1.0 / 2.0_f64.sqrt(), 1.0 / 2.0_f64.sqrt())]
            ],
            Gate::TDG(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(1.0 / 2.0_f64.sqrt(), -1.0 / 2.0_f64.sqrt())]
            ],
            Gate::S(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 1.0)]
            ],
            Gate::SDG(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, -1.0)]
            ],
            Gate::X(_) => array![
                [c64::new(0.0, 0.0), c64::new(1.0, 0.0)],
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)]
            ],
            Gate::Y(_) => array![
                [c64::new(0.0, 0.0), c64::new(0.0, -1.0)],
                [c64::new(0.0, 1.0), c64::new(0.0, 0.0)]
            ],
            Gate::Z(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(-1.0, 0.0)]
            ],
            Gate::SWAP(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(1.0, 0.0)]
            ],
            Gate::CX(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(1.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0)]
            ],
            Gate::CZ(_) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(-1.0, 0.0)]
            ],
            Gate::U1(gate) => array![
                [c64::new(1.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(gate.lambda.cos(), gate.lambda.sin())]
            ],
            Gate::U2(gate) => array![
                [c64::new(1.0, 0.0), c64::new(- gate.lambda.cos(), - gate.lambda.sin())],
                [c64::new(gate.phi.cos(), gate.phi.sin()), c64::new((gate.phi + gate.lambda).cos(), (gate.phi + gate.lambda).sin())]
            ] / 2.0_f64.sqrt(),
            Gate::U3(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), -gate.lambda.sin() / 2.0), c64::new(-gate.lambda.cos() * (gate.theta / 2.0).sin(), gate.lambda.cos() * (gate.theta / 2.0).cos())],
                [c64::new(gate.phi.sin() * gate.theta / 2.0, gate.phi.cos() * gate.theta / 2.0), c64::new((gate.phi + gate.lambda).sin() * gate.theta / 2.0, (gate.phi + gate.lambda).cos() * gate.theta / 2.0)]
            ],
            Gate::RX(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, -(gate.theta / 2.0).sin())],
                [c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new((gate.theta / 2.0).cos(), 0.0)]
            ],
            Gate::RY(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(-(gate.theta / 2.0).sin(), 0.0)],
                [c64::new((gate.theta / 2.0).sin(), 0.0), c64::new((gate.theta / 2.0).cos(), 0.0)]
            ],
            Gate::RZ(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), (gate.theta / 2.0).sin())]
            ],
            Gate::RXX(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, -(gate.theta / 2.0).sin())],
                [c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), 0.0)]
            ],
            Gate::RYY(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, -(gate.theta / 2.0).sin())],
                [c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new((gate.theta / 2.0).cos(), 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), 0.0)]
            ],
            Gate::RZZ(gate) => array![
                [c64::new((gate.theta / 2.0).cos(), -(gate.theta / 2.0).sin()), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), (gate.theta / 2.0).sin()), c64::new(0.0, 0.0), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), (gate.theta / 2.0).sin()), c64::new(0.0, 0.0)],
                [c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new(0.0, 0.0), c64::new((gate.theta / 2.0).cos(), -(gate.theta / 2.0).sin())]
            ],
            Gate::Can(gate) => {
                let zero = c64::new(0.0, 0.0);
                let cosm = (gate.theta1 / 2.0 - gate.theta2 / 2.0).cos();
                let cosp = (gate.theta1 / 2.0 + gate.theta2 / 2.0).cos();
                let sinm = (gate.theta1 / 2.0 - gate.theta2 / 2.0).sin();
                let sinp = (gate.theta1 / 2.0 + gate.theta2 / 2.0).sin();
                let eim = c64::new(0.0, -gate.theta3 / 2.0).exp();
                let eip = c64::new(0.0, gate.theta3 / 2.0).exp();
                array![
                    [eim * cosm, zero, zero, -1.0 * eim * sinm],
                    [zero, eip * cosp, -1.0 * eip * sinp, zero],
                    [zero, -1.0 * eip * sinp, eip * cosp, zero],
                    [-1.0 * eim * sinm, zero, zero, eim * cosm]
                ]
            }
        }
    }

    pub fn hermitian(&self) -> Self {
        match self {
            Gate::H(_) => Gate::H(HGate{}),
            Gate::T(_) => Gate::T(TGate{}),
            Gate::TDG(_) => Gate::TDG(TDGGate{}),
            Gate::S(_) => Gate::SDG(SDGGate{}),
            Gate::SDG(_) => Gate::S(SGate{}),
            Gate::X(_) => Gate::X(XGate{}),
            Gate::Y(_) => Gate::Y(YGate{}),
            Gate::Z(_) => Gate::Z(ZGate{}),
            Gate::SWAP(_) => Gate::SWAP(SWAPGate{}),
            Gate::CX(_) => Gate::CX(CXGate{}),
            Gate::CZ(_) => Gate::CZ(CZGate{}),
            Gate::U1(gate) => Gate::U1(U1Gate{lambda: -gate.lambda}),
            Gate::U2(gate) => Gate::U3(U3Gate{theta: - PI / 2.0, phi: -gate.lambda, lambda: -gate.phi}),
            Gate::U3(gate) => Gate::U3(U3Gate{theta: -gate.theta, phi: -gate.lambda, lambda: -gate.phi}),
            Gate::RX(gate) => Gate::RX(RXGate{theta: -gate.theta}),
            Gate::RY(gate) => Gate::RY(RYGate{theta: -gate.theta}),
            Gate::RZ(gate) => Gate::RZ(RZGate{theta: -gate.theta}),
            Gate::RXX(gate) => Gate::RXX(RXXGate{theta: -gate.theta}),
            Gate::RYY(gate) => Gate::RYY(RYYGate{theta: -gate.theta}),
            Gate::RZZ(gate) => Gate::RZZ(RZZGate{theta: -gate.theta}),
            Gate::Can(gate) => Gate::Can(CanonicalGate{theta1: -gate.theta1, theta2: -gate.theta2, theta3: -gate.theta3}),
        }
    }


    // fn name(&self) -> &str;
    // fn on(&self, tqs: Vec<usize>, cqs: Vec<usize>) -> Self;
    // fn n_qubits(&self) -> usize {
    //     (self.data().shape()[0] as f64).log2() as usize
    // }
    // fn tqs(&self) -> Vec<usize>;
    // fn cqs(&self) -> Vec<usize>;
    // fn tq(&self) -> Result<usize, &str> {
    //     if self.tqs().len() > 1 {
    //         Err("tqs length is greater than 1")
    //     } else {
    //         Ok(self.tqs()[0])
    //     }
    // }
    // fn cq(&self) -> Result<usize, &str> {
    //     if self.cqs().len() > 1 {
    //         Err("cqs length is greater than 1")
    //     } else {
    //         Ok(self.cqs()[0])
    //     }
    // }
    // fn qregs(&self) -> Vec<usize> {
    //     let mut qregs = self.tqs();
    //     qregs.extend(self.cqs());
    //     qregs
    // }
    // fn data(&self) -> Array2<f64>;
    // fn hermitian(&self) -> Self;
    // fn num_qregs(&self) -> usize;
}






// impl Gate {
//     pub fn new(name: &str, tqs: Vec<usize>, cqs: Vec<usize>, n_qubits: usize, angle: f64, angles: Vec<f64>, params: Vec<f64>, exponent: f64) -> Self {
//         Gate {
//             name: name.to_string(),
//             tqs,
//             cqs,
//             n_qubits,
//             angle,
//             angles,
//             params,
//             exponent,
//         }
//     }
    
// }
