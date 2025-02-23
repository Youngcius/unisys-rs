use ndarray::{Array1, Array2, Zip};
use ndarray_linalg::c64;

use crate::basic::circuits::Circuit;

/*
if num_qubits is None:
    num_qubits = circ.num_qubits_with_dummy
left_end = np.full(num_qubits, -1)
circ_part = front_full_width_circuit(circ, lambda g: g.num_qregs > 1)
for num_layer, layer in enumerate(circ_part.layer()):
    for q in reduce(add, [g.qregs for g in layer]):
        if left_end[q] < 0:
            left_end[q] = num_layer
    if np.all(left_end >= 0):
        break
left_end[left_end == -1] = left_end.max() + 1
return left_end
 */
pub fn left_end_empty_layer(circ: &Circuit, num_qubits: Option<usize>) -> Array1<usize> {
    let mut n = circ.num_qubits_with_dummy();
    if let Some(num_qubits) = num_qubits {
        n = num_qubits;
    }

    let left_end = Array1::from_elem(n, -1);
    let ffwc = circ.front_full_width_circuit();
    // for (num_layer, layer) in ffwc.layer().iter().enumerate() {
    //     for q in layer.iter().flat_map(|g| g.qregs()) {
    //         if left_end[q] < 0 {
    //             left_end[q] = num_layer;
    //         }
    //     }
    //     if left_end.iter().all(|&x| x >= 0) {
    //         break;
    //     }
    // }

    left_end.mapv(|x| x as usize)
}

pub fn depth_overhead(lhs: &Array1<usize>, rhs: &Array1<usize>) -> usize {
    let lhs_bool: Array1<bool> = lhs.mapv(|x| x == 0); // herein we use x == 0 instead of x != 0
    let rhs_bool: Array1<bool> = rhs.mapv(|x| x == 0);
    let xor: Array1<bool> = lhs_bool.clone() ^ rhs_bool.clone();
    let condition = lhs_bool | rhs_bool;
    let filtered_xor = Zip::from(&xor)
        .and(&condition)
        .map_collect(|&x, &c| if c { x } else { false });
    if filtered_xor.iter().all(|&x| x) {
        (lhs + rhs).sum() - xor.len()
    } else {
        (lhs + rhs).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_depth_overhead() {
        let lhs = Array1::from_vec(vec![4, 1, 2, 0, 5, 0, 3]);
        let rhs = Array1::from_vec(vec![3, 0, 0, 0, 0, 1, 2]);
        let cost = depth_overhead(&lhs, &rhs);
        assert_eq!(cost, 21);
        println!("{}", cost);
    }
}
