use crate::{c, i, r};
use ndarray::{array, Array2};
use super::matrices::Dagger;
use ndarray_linalg::c64;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub enum Gate {
    Univ(UnivGate),
    H(HGate),
    T(TGate),
    TDG(TDGGate),
    S(SGate),
    SDG(SDGGate),
    X(XGate),
    Y(YGate),
    Z(ZGate),
    SWAP(SWAPGate),
    U1(U1Gate),
    U2(U2Gate),
    U3(U3Gate),
    RX(RXGate),
    RY(RYGate),
    RZ(RZGate),
    RXX(RXXGate),
    RYY(RYYGate),
    RZZ(RZZGate),
    Can(CanonicalGate),
}

#[derive(Clone, Debug)]
pub struct UnivGate {
    pub untry: Array2<c64>,
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
            Gate::Univ(_) => "Univ",
            Gate::H(_) => "H",
            Gate::T(_) => "T",
            Gate::TDG(_) => "TDG",
            Gate::S(_) => "S",
            Gate::SDG(_) => "SDG",
            Gate::X(_) => "X",
            Gate::Y(_) => "Y",
            Gate::Z(_) => "Z",
            Gate::SWAP(_) => "SWAP",
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
            Gate::Univ(univ) => {
                let n = univ.untry.shape()[0];
                (n as f64).log2() as usize
            }
            Gate::H(_) => 1,
            Gate::T(_) => 1,
            Gate::TDG(_) => 1,
            Gate::S(_) => 1,
            Gate::SDG(_) => 1,
            Gate::X(_) => 1,
            Gate::Y(_) => 1,
            Gate::Z(_) => 1,
            Gate::SWAP(_) => 2,
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
            Gate::Univ(_) => Err("Univ gate has no angle"),
            Gate::H(_) => Err("H gate has no angle"),
            Gate::T(_) => Err("T gate has no angle"),
            Gate::TDG(_) => Err("TDG gate has no angle"),
            Gate::S(_) => Err("S gate has no angle"),
            Gate::SDG(_) => Err("SDG gate has no angle"),
            Gate::X(_) => Err("X gate has no angle"),
            Gate::Y(_) => Err("Y gate has no angle"),
            Gate::Z(_) => Err("Z gate has no angle"),
            Gate::SWAP(_) => Err("SWAP gate has no angle"),
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
            Gate::Univ(_) => Err("Univ gate has no angles"),
            Gate::H(_) => Err("H gate has no angles"),
            Gate::T(_) => Err("T gate has no angles"),
            Gate::TDG(_) => Err("TDG gate has no angles"),
            Gate::S(_) => Err("S gate has no angles"),
            Gate::SDG(_) => Err("SDG gate has no angles"),
            Gate::X(_) => Err("X gate has no angles"),
            Gate::Y(_) => Err("Y gate has no angles"),
            Gate::Z(_) => Err("Z gate has no angles"),
            Gate::SWAP(_) => Err("SWAP gate has no angles"),
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
            Gate::Univ(gate) => gate.untry.clone(),
            Gate::H(_) => array![[r!(1.0), r!(1.0)], [r!(1.0), r!(-1.0)]] / 2.0_f64.sqrt(),
            Gate::T(_) => array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), (r!(1.0) + i!(1.0)) / 2.0_f64.sqrt()]
            ],
            Gate::TDG(_) => array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), (r!(1.0) - i!(1.0)) / 2.0_f64.sqrt()]
            ],
            Gate::S(_) => array![[r!(1.0), r!(0.0)], [r!(0.0), i!(1.0)]],
            Gate::SDG(_) => array![[r!(1.0), r!(0.0)], [r!(0.0), i!(-1.0)]],
            Gate::X(_) => array![[r!(0.0), r!(1.0)], [r!(1.0), r!(0.0)]],
            Gate::Y(_) => array![[r!(0.0), i!(-1.0)], [i!(1.0), r!(0.0)]],
            Gate::Z(_) => array![[r!(1.0), r!(0.0)], [r!(0.0), r!(-1.0)]],
            Gate::SWAP(_) => array![
                [r!(1.0), r!(0.0), r!(0.0), r!(0.0)],
                [r!(0.0), r!(0.0), r!(1.0), r!(0.0)],
                [r!(0.0), r!(1.0), r!(0.0), r!(0.0)],
                [r!(0.0), r!(0.0), r!(0.0), r!(1.0)]
            ],
            Gate::U1(gate) => array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), c!(gate.lambda.cos(), gate.lambda.sin())]
            ],
            Gate::U2(gate) => {
                array![
                    [r!(1.0), c!(-gate.lambda.cos(), -gate.lambda.sin())],
                    [
                        c!(gate.phi.cos(), gate.phi.sin()),
                        c!(
                            (gate.phi + gate.lambda).cos(),
                            (gate.phi + gate.lambda).sin()
                        )
                    ]
                ] / 2.0_f64.sqrt()
            }
            Gate::U3(gate) => array![
                [
                    r!((gate.theta / 2.0).cos()),
                    c!(
                        -gate.lambda.cos() * (gate.theta / 2.0).sin(),
                        -gate.lambda.sin() * (gate.theta / 2.0).sin()
                    )
                ],
                [
                    c!(
                        gate.phi.cos() * (gate.theta / 2.0).sin(),
                        gate.phi.sin() * (gate.theta / 2.0).sin()
                    ),
                    c!(
                        (gate.phi + gate.lambda).cos() * (gate.theta / 2.0).cos(),
                        (gate.phi + gate.lambda).sin() * (gate.theta / 2.0).cos()
                    )
                ]
            ],
            Gate::RX(gate) => array![
                [r!((gate.theta / 2.0).cos()), i!(-(gate.theta / 2.0).sin())],
                [i!(-(gate.theta / 2.0).sin()), r!((gate.theta / 2.0).cos())]
            ],
            Gate::RY(gate) => array![
                [r!(gate.theta / 2.0).cos(), r!(-(gate.theta / 2.0).sin())],
                [r!((gate.theta / 2.0).sin()), r!((gate.theta / 2.0).cos())]
            ],
            Gate::RZ(gate) => array![
                [
                    c!((gate.theta / 2.0).cos(), -(gate.theta / 2.0).sin()),
                    r!(0.0)
                ],
                [
                    r!(0.0),
                    c!((gate.theta / 2.0).cos(), (gate.theta / 2.0).sin())
                ]
            ],
            Gate::RXX(gate) => {
                let cos = r!((gate.theta / 2.).cos());
                let sin = i!(-(gate.theta / 2.).sin());
                let zero = r!(0.0);
                array![
                    [cos, zero, zero, sin],
                    [zero, cos, sin, zero],
                    [zero, sin, cos, zero],
                    [sin, zero, zero, cos]
                ]
            }
            Gate::RYY(gate) => {
                let cos = r!((gate.theta / 2.).cos());
                let nsin = i!(-1.0) * (gate.theta / 2.).sin();
                let psin = i!(1.0) * (gate.theta / 2.).sin();
                let zero = r!(0.0);
                array![
                    [cos, zero, zero, psin],
                    [zero, cos, nsin, zero],
                    [zero, nsin, cos, zero],
                    [psin, zero, zero, cos]
                ]
            }
            Gate::RZZ(gate) => {
                let pos = (i!(1.) * gate.theta / 2.).exp();
                let neg = (i!(-1.) * gate.theta / 2.).exp();
                let zero = r!(0.0);
                array![
                    [neg, zero, zero, zero],
                    [zero, pos, zero, zero],
                    [zero, zero, pos, zero],
                    [zero, zero, zero, neg]
                ]
            }
            Gate::Can(gate) => {
                let zero = r!(0.0);
                let cosm = (gate.theta1 / 2.0 - gate.theta2 / 2.0).cos();
                let cosp = (gate.theta1 / 2.0 + gate.theta2 / 2.0).cos();
                let sinm = (gate.theta1 / 2.0 - gate.theta2 / 2.0).sin();
                let sinp = (gate.theta1 / 2.0 + gate.theta2 / 2.0).sin();
                let eim = c!(0.0, -gate.theta3 / 2.0).exp();
                let eip = c!(0.0, gate.theta3 / 2.0).exp();
                array![
                    [eim * cosm, zero, zero, i!(-1.0) * eim * sinm],
                    [zero, eip * cosp, i!(-1.0) * eip * sinp, zero],
                    [zero, i!(-1.0) * eip * sinp, eip * cosp, zero],
                    [i!(-1.0) * eim * sinm, zero, zero, eim * cosm]
                ]
            }
        }
    }

    pub fn hermitian(&self) -> Self {
        match self {
            Gate::Univ(univ) => Gate::Univ(UnivGate {untry: univ.untry.dagger()}),
            Gate::H(_) => Gate::H(HGate {}),
            Gate::T(_) => Gate::T(TGate {}),
            Gate::TDG(_) => Gate::TDG(TDGGate {}),
            Gate::S(_) => Gate::SDG(SDGGate {}),
            Gate::SDG(_) => Gate::S(SGate {}),
            Gate::X(_) => Gate::X(XGate {}),
            Gate::Y(_) => Gate::Y(YGate {}),
            Gate::Z(_) => Gate::Z(ZGate {}),
            Gate::SWAP(_) => Gate::SWAP(SWAPGate {}),
            Gate::U1(gate) => Gate::U1(U1Gate {
                lambda: -gate.lambda,
            }),
            Gate::U2(gate) => Gate::U3(U3Gate {
                theta: -PI / 2.0,
                phi: -gate.lambda,
                lambda: -gate.phi,
            }),
            Gate::U3(gate) => Gate::U3(U3Gate {
                theta: -gate.theta,
                phi: -gate.lambda,
                lambda: -gate.phi,
            }),
            Gate::RX(gate) => Gate::RX(RXGate { theta: -gate.theta }),
            Gate::RY(gate) => Gate::RY(RYGate { theta: -gate.theta }),
            Gate::RZ(gate) => Gate::RZ(RZGate { theta: -gate.theta }),
            Gate::RXX(gate) => Gate::RXX(RXXGate { theta: -gate.theta }),
            Gate::RYY(gate) => Gate::RYY(RYYGate { theta: -gate.theta }),
            Gate::RZZ(gate) => Gate::RZZ(RZZGate { theta: -gate.theta }),
            Gate::Can(gate) => Gate::Can(CanonicalGate {
                theta1: -gate.theta1,
                theta2: -gate.theta2,
                theta3: -gate.theta3,
            }),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::utils::ops::allclose;
    use ndarray::Array;

    #[test]
    fn test_data() {
        // use ndarray generate three random f64
        use ndarray_rand::rand_distr::Uniform;
        use ndarray_rand::RandomExt;

        // test for U2 gate
        let angles = Array::random((2,), Uniform::new(0.0, 3.0));
        let phi = angles[0];
        let lambda = angles[1];
        let u2 = Gate::U2(U2Gate {
            phi: phi,
            lambda: lambda,
        });
        let data1 = u2.data();
        let u2_from_u3 = Gate::U3(U3Gate {
            theta: PI / 2.0,
            phi: phi,
            lambda: lambda,
        });
        let data2 = u2_from_u3.data();
        assert!(allclose(&data1, &data2));

        // test for Canonical gate
        let angles = Array::random((3,), Uniform::new(0.0, 3.0));
        let theta1 = angles[0];
        let theta2 = angles[1];
        let theta3 = angles[2];
        let can = Gate::Can(CanonicalGate {
            theta1: theta1,
            theta2: theta2,
            theta3: theta3,
        });
        let rxx = Gate::RXX(RXXGate { theta: theta1 });
        let ryy = Gate::RYY(RYYGate { theta: theta2 });
        let rzz = Gate::RZZ(RZZGate { theta: theta3 });
        let data1 = can.data();
        let data2 = rxx.data().dot(&ryy.data()).dot(&rzz.data());
        assert!(allclose(&data1, &data2));
    }
}
