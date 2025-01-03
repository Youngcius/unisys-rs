use super::matrices::Dagger;
use crate::{c, i, r};
use ndarray::{array, Array2};
use ndarray_linalg::c64;
use std::f64::consts::PI;

#[derive(Clone, Debug)]
pub enum GateType {
    Univ,
    H,
    T,
    Tdg,
    S,
    Sdg,
    X,
    Y,
    Z,
    Swap,
    U1,
    U2,
    U3,
    Rx,
    Ry,
    Rz,
    Rxx,
    Ryy,
    Rzz,
    Can,
}

impl std::fmt::Display for GateType {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        match self {
            GateType::Univ => write!(f, "Univ"),
            GateType::H => write!(f, "H"),
            GateType::T => write!(f, "T"),
            GateType::Tdg => write!(f, "Tdg"),
            GateType::S => write!(f, "S"),
            GateType::Sdg => write!(f, "Sdg"),
            GateType::X => write!(f, "X"),
            GateType::Y => write!(f, "Y"),
            GateType::Z => write!(f, "Z"),
            GateType::Swap => write!(f, "Swap"),
            GateType::U1 => write!(f, "U1"),
            GateType::U2 => write!(f, "U2"),
            GateType::U3 => write!(f, "U3"),
            GateType::Rx => write!(f, "Rx"),
            GateType::Ry => write!(f, "Ry"),
            GateType::Rz => write!(f, "Rz"),
            GateType::Rxx => write!(f, "Rxx"),
            GateType::Ryy => write!(f, "Ryy"),
            GateType::Rzz => write!(f, "Rzz"),
            GateType::Can => write!(f, "Can"),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Gate {
    pub gate_type: GateType,
    pub n_qubits: usize,
    pub data: Array2<c64>,
    pub params: Option<Vec<f64>>,
}

impl std::fmt::Display for Gate {
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        let mut prefix = self.gate_type.to_string();
        if let Some(params) = &self.params {
            let params_str = params
                .iter()
                .map(|a| format!("{:.3}", a))
                .collect::<Vec<_>>()
                .join(",");
            prefix += &format!("({})", params_str);
        }
        write!(f, "{}", prefix)
    }
}

impl Gate {
    pub fn new(gate_type: GateType, n_qubits: usize, data: Array2<c64>) -> Self {
        Self {
            gate_type,
            n_qubits,
            data,
            params: None,
        }
    }

    pub fn from_gate(gate: Gate) -> Self {
        Self {
            gate_type: gate.gate_type,
            n_qubits: gate.n_qubits,
            data: gate.data,
            params: gate.params,
        }
    }

    pub fn univ(data: Array2<c64>) -> Self {
        Self {
            gate_type: GateType::Univ,
            n_qubits: (data.shape()[0] as f64).log2() as usize,
            data,
            params: None,
        }
    }

    pub fn h() -> Self {
        Self {
            gate_type: GateType::H,
            n_qubits: 1,
            data: array![[r!(1.0), r!(1.0)], [r!(1.0), r!(-1.0)]] / 2.0_f64.sqrt(),
            params: None,
        }
    }

    pub fn t() -> Self {
        Self {
            gate_type: GateType::T,
            n_qubits: 1,
            data: array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), (r!(1.0) + i!(1.0)) / 2.0_f64.sqrt()]
            ],
            params: None,
        }
    }

    pub fn tdg() -> Self {
        Self {
            gate_type: GateType::Tdg,
            n_qubits: 1,
            data: array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), (r!(1.0) - i!(1.0)) / 2.0_f64.sqrt()]
            ],
            params: None,
        }
    }

    pub fn s() -> Self {
        Self {
            gate_type: GateType::S,
            n_qubits: 1,
            data: array![[r!(1.0), r!(0.0)], [r!(0.0), i!(1.0)]],
            params: None,
        }
    }

    pub fn sdg() -> Self {
        Self {
            gate_type: GateType::Sdg,
            n_qubits: 1,
            data: array![[r!(1.0), r!(0.0)], [r!(0.0), i!(-1.0)]],
            params: None,
        }
    }

    pub fn x() -> Self {
        Self {
            gate_type: GateType::X,
            n_qubits: 1,
            data: array![[r!(0.0), r!(1.0)], [r!(1.0), r!(0.0)]],
            params: None,
        }
    }

    pub fn y() -> Self {
        Self {
            gate_type: GateType::Y,
            n_qubits: 1,
            data: array![[r!(0.0), i!(-1.0)], [i!(1.0), r!(0.0)]],
            params: None,
        }
    }

    pub fn z() -> Self {
        Self {
            gate_type: GateType::Z,
            n_qubits: 1,
            data: array![[r!(1.0), r!(0.0)], [r!(0.0), r!(-1.0)]],
            params: None,
        }
    }

    pub fn swap() -> Self {
        Self {
            gate_type: GateType::Swap,
            n_qubits: 2,
            data: array![
                [r!(1.0), r!(0.0), r!(0.0), r!(0.0)],
                [r!(0.0), r!(0.0), r!(1.0), r!(0.0)],
                [r!(0.0), r!(1.0), r!(0.0), r!(0.0)],
                [r!(0.0), r!(0.0), r!(0.0), r!(1.0)]
            ],
            params: None,
        }
    }

    pub fn u1(lambda: f64) -> Self {
        Self {
            gate_type: GateType::U1,
            n_qubits: 1,
            data: array![
                [r!(1.0), r!(0.0)],
                [r!(0.0), c!(lambda.cos(), lambda.sin())]
            ],
            params: Some(vec![lambda]),
        }
    }

    pub fn u2(phi: f64, lambda: f64) -> Self {
        Self {
            gate_type: GateType::U2,
            n_qubits: 1,
            data: array![
                [r!(1.0), c!(-lambda.cos(), -lambda.sin())],
                [
                    c!(phi.cos(), phi.sin()),
                    c!((phi + lambda).cos(), (phi + lambda).sin())
                ]
            ] / 2.0_f64.sqrt(),
            params: Some(vec![phi, lambda]),
        }
    }

    pub fn u3(theta: f64, phi: f64, lambda: f64) -> Self {
        Self {
            gate_type: GateType::U3,
            n_qubits: 1,
            data: array![
                [
                    r!((theta / 2.0).cos()),
                    c!(
                        -lambda.cos() * (theta / 2.0).sin(),
                        -lambda.sin() * (theta / 2.0).sin()
                    )
                ],
                [
                    c!(
                        phi.cos() * (theta / 2.0).sin(),
                        phi.sin() * (theta / 2.0).sin()
                    ),
                    c!(
                        (phi + lambda).cos() * (theta / 2.0).cos(),
                        (phi + lambda).sin() * (theta / 2.0).cos()
                    )
                ]
            ],
            params: Some(vec![theta, phi, lambda]),
        }
    }

    pub fn rx(theta: f64) -> Self {
        Self {
            gate_type: GateType::Rx,
            n_qubits: 1,
            data: array![
                [r!((theta / 2.0).cos()), i!(-(theta / 2.0).sin())],
                [i!(-(theta / 2.0).sin()), r!((theta / 2.0).cos())]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn ry(theta: f64) -> Self {
        Self {
            gate_type: GateType::Ry,
            n_qubits: 1,
            data: array![
                [r!(theta / 2.0).cos(), r!(-(theta / 2.0).sin())],
                [r!((theta / 2.0).sin()), r!((theta / 2.0).cos())]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn rz(theta: f64) -> Self {
        Self {
            gate_type: GateType::Rz,
            n_qubits: 1,
            data: array![
                [c!((theta / 2.0).cos(), -(theta / 2.0).sin()), r!(0.0)],
                [r!(0.0), c!((theta / 2.0).cos(), (theta / 2.0).sin())]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn rxx(theta: f64) -> Self {
        let cos = r!((theta / 2.).cos());
        let sin = i!(-(theta / 2.).sin());
        let zero = r!(0.0);
        Self {
            gate_type: GateType::Rxx,
            n_qubits: 2,
            data: array![
                [cos, zero, zero, sin],
                [zero, cos, sin, zero],
                [zero, sin, cos, zero],
                [sin, zero, zero, cos]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn ryy(theta: f64) -> Self {
        let cos = r!((theta / 2.).cos());
        let nsin = i!(-1.0) * (theta / 2.).sin();
        let psin = i!(1.0) * (theta / 2.).sin();
        let zero = r!(0.0);
        Self {
            gate_type: GateType::Ryy,
            n_qubits: 2,
            data: array![
                [cos, zero, zero, psin],
                [zero, cos, nsin, zero],
                [zero, nsin, cos, zero],
                [psin, zero, zero, cos]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn rzz(theta: f64) -> Self {
        let pos = (i!(1.) * theta / 2.).exp();
        let neg = (i!(-1.) * theta / 2.).exp();
        let zero = r!(0.0);
        Self {
            gate_type: GateType::Rzz,
            n_qubits: 2,
            data: array![
                [neg, zero, zero, zero],
                [zero, pos, zero, zero],
                [zero, zero, pos, zero],
                [zero, zero, zero, neg]
            ],
            params: Some(vec![theta]),
        }
    }

    pub fn can(x: f64, y: f64, z: f64) -> Self {
        let zero = r!(0.0);
        let cosm = ((x - y) * PI / 2.0).cos();
        let cosp = ((x + y) * PI / 2.0).cos();
        let sinm = ((x - y) * PI / 2.0).sin();
        let sinp = ((x + y) * PI / 2.0).sin();
        let eim = c!(0.0, -z * PI / 2.0).exp();
        let eip = c!(0.0, z * PI / 2.0).exp();
        Self {
            gate_type: GateType::Can,
            n_qubits: 2,
            data: array![
                [eim * cosm, zero, zero, i!(-1.0) * eim * sinm],
                [zero, eip * cosp, i!(-1.0) * eip * sinp, zero],
                [zero, i!(-1.0) * eip * sinp, eip * cosp, zero],
                [i!(-1.0) * eim * sinm, zero, zero, eim * cosm]
            ],
            params: Some(vec![x, y, z]),
        }
    }

    pub fn hermitian(&self) -> Self {
        match self.gate_type {
            GateType::Univ => Gate::univ(self.data.dagger()),
            GateType::H => self.clone(),
            GateType::T => Gate::tdg(),
            GateType::Tdg => Gate::t(),
            GateType::S => Gate::sdg(),
            GateType::Sdg => Gate::s(),
            GateType::X => self.clone(),
            GateType::Y => self.clone(),
            GateType::Z => self.clone(),
            GateType::Swap => self.clone(),
            GateType::U1 => Gate::u1(-self.params.as_ref().unwrap()[0]),
            GateType::U2 => Gate::u2(
                -self.params.as_ref().unwrap()[1] - PI,
                -self.params.as_ref().unwrap()[0] + PI,
            ),
            GateType::U3 => Gate::u3(
                -self.params.as_ref().unwrap()[0],
                -self.params.as_ref().unwrap()[1],
                -self.params.as_ref().unwrap()[2],
            ),
            GateType::Rx => Gate::rx(-self.params.as_ref().unwrap()[0]),
            GateType::Ry => Gate::ry(-self.params.as_ref().unwrap()[0]),
            GateType::Rz => Gate::rz(-self.params.as_ref().unwrap()[0]),
            GateType::Rxx => Gate::rxx(-self.params.as_ref().unwrap()[0]),
            GateType::Ryy => Gate::ryy(-self.params.as_ref().unwrap()[0]),
            GateType::Rzz => Gate::rzz(-self.params.as_ref().unwrap()[0]),
            GateType::Can => Gate::can(
                -self.params.as_ref().unwrap()[0],
                -self.params.as_ref().unwrap()[1],
                -self.params.as_ref().unwrap()[2],
            ),
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
        let u2 = Gate::u2(phi, lambda);
        let data1 = u2.data;
        let u2_from_u3 = Gate::u3(PI / 2.0, phi, lambda);
        let data2 = u2_from_u3.data;
        assert!(allclose(&data1, &data2));

        // test for Canonical gate
        let angles = Array::random((3,), Uniform::new(0.0, 3.0));
        let theta1 = angles[0];
        let theta2 = angles[1];
        let theta3 = angles[2];
        let can = Gate::can(theta1 / PI, theta2 / PI, theta3 / PI);

        let (rxx, ryy, rzz) = (Gate::rxx(theta1), Gate::ryy(theta2), Gate::rzz(theta3));
        let data1 = can.data;
        let data2 = rxx.data.dot(&ryy.data).dot(&rzz.data);
        assert!(allclose(&data1, &data2));
    }

    #[test]
    fn test_aaa() {
        use ndarray::Array;
        use ndarray::Array2;
        let a: Array2<c64> = Array::zeros((3, 2));
        let b: Array2<c64> = Array2::eye(3);
        println!("{:?}", a);
        println!("{}", b);
    }
}
