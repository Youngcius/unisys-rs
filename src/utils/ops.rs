use crate::{c, i, r};
use ndarray::Array2;
use ndarray_linalg::c64;

pub fn allclose(data1: &Array2<c64>, data2: &Array2<c64>) -> bool {
    let diff = data1 - data2;
    let diff = diff.mapv(|x| x.norm());
    let tol = 1e-8;
    diff.iter().all(|&x| x < tol)
}

// use ndarray Array2<c64> to compute its matrix trace:
// pub fn trace(data: &Array2<c64>) -> c64 {
//     data.diag().sum()
// }

pub fn match_global_phase(a: Array2<c64>, b: Array2<c64>) -> (Array2<c64>, Array2<c64>) {
    // If the shapes are different or either is empty, return copies
    if a.shape() != b.shape() || a.len() == 0 {
        return (a.clone(), b.clone());
    }

    // Find the index k of the element with the largest magnitude in b
    let mut max_mag = 0.0;
    let mut k = (0, 0);
    let shape = b.shape();
    for i in 0..shape[0] {
        for j in 0..shape[1] {
            let mag = b[(i, j)].norm();
            if mag > max_mag {
                max_mag = mag;
                k = (i, j);
            }
        }
    }

    // Calculate the multiplier required to zero the phase of a complex number
    fn dephase(v: c64) -> c64 {
        let re = v.re;
        let im = v.im;
        if im == 0.0 {
            // 避免浮点误差: 实数轴上，如果 r<0 则乘 -1，否则乘 +1
            if re < 0.0 {
                r!(-1.0)
            } else {
                r!(1.0)
            }
        } else if re == 0.0 {
            // English: If i < 0, multiply by 1j; otherwise, multiply by -1j
            if im < 0.0 {
                i!(1.0)
            } else {
                i!(-1.0)
            }
        } else {
            // Multiply by e^{-i * arg(v)} to zero the phase
            let angle = im.atan2(re); // i.e., arg(v)
            i!(-angle).exp()
        }
    }

    // Calculate the multiplier required to zero the phase at index k
    let phase_factor_a = dephase(a[k]);
    let phase_factor_b = dephase(b[k]);

    // Apply the phase correction to the entire matrices a and b
    let a_prime = a.mapv(|val| val * phase_factor_a);
    let b_prime = b.mapv(|val| val * phase_factor_b);

    (a_prime, b_prime)
}

pub fn equiv_unitary(data1: &Array2<c64>, data2: &Array2<c64>) -> bool {
    let (data1, data2) = match_global_phase(data1.clone(), data2.clone());
    allclose(&data1, &data2)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic::gates;
    use ndarray::Array;
    use ndarray_rand::rand_distr::Uniform;
    use ndarray_rand::RandomExt;

    #[test]
    fn test_equiv_unitary() {
        let angles = Array::random((3,), Uniform::new(0.0, 3.0));
        let theta1 = angles[0];
        let theta2 = angles[1];
        let theta3 = angles[2];
        let u3 = gates::Gate::u3(theta1, theta2, theta3);
        let mat1 = u3.data;
        let phase = i!(1.234).exp();
        let mat2 = &mat1 * phase;

        println!("mat1: {}", mat1);
        println!("mat2: {}", mat2);
        assert!(equiv_unitary(&mat1, &mat2));
    }
}
