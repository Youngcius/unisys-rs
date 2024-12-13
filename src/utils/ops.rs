use crate::{c, i};
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
    // 如果形状不同或任意一个为空，则直接返回拷贝
    if a.shape() != b.shape() || a.len() == 0 {
        return (a.clone(), b.clone());
    }

    // 在 b 中找到模最大的元素下标 k
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

    // 定义 dephase 函数，用于计算将某个复数 "相位归零" 所需的乘子
    fn dephase(v: c64) -> c64 {
        let r = v.re;
        let i = v.im;
        // 与 Python 版本等价的逻辑
        if i == 0.0 {
            // 避免浮点误差: 实数轴上，如果 r<0 则乘 -1，否则乘 +1
            if r < 0.0 {
                c!(-1.0, 0.0)
            } else {
                c!(1.0, 0.0)
            }
        } else if r == 0.0 {
            // 纯虚数轴上，如果 i<0 则乘 1j，否则乘 -1j
            // 注意 Rust 中的虚数单位为 c!(0.0, 1.0)
            if i < 0.0 {
                c!(0.0, 1.0)
            } else {
                c!(0.0, -1.0)
            }
        } else {
            // 一般情况，乘上 e^{-i * arg(v)} 使其相位归零
            let angle = i.atan2(r); // 即 arg(v)
            i!(-angle).exp()
        }
    }

    // 计算在索引 k 处相位归零所需的乘子
    let phase_factor_a = dephase(a[k]);
    let phase_factor_b = dephase(b[k]);

    // 对 a 和 b 全矩阵进行相位修正
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
        let u3 = gates::Gate::U3(gates::U3Gate {
            theta: theta1,
            phi: theta2,
            lambda: theta3,
        });
        let mat1 = u3.data();
        let phase = i!(1.234).exp();
        let mat2 = &mat1 * phase;

        println!("mat1: {}", mat1);
        println!("mat2: {}", mat2);
        assert!(equiv_unitary(&mat1, &mat2));
    }
}
