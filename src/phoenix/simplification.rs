use crate::basic::gates::{Clifford2Q, Gate};
use crate::basic::operations::Operation;
use crate::models::paulis::BSF;
use itertools::Itertools;
use ndarray::Array2;
use ndarray::Axis;
use rayon::prelude::*;

pub fn simplify_bsf(bsf: &BSF) -> (BSF, Vec<(Operation, BSF)>) {
    let mut bsf = bsf.clone();
    let mut cliffs_with_locals: Vec<(Operation, BSF)> = Vec::new();
    let mut avoid: (usize, usize) = (bsf.num_qubits(), bsf.num_qubits());
    while bsf.total_weight() > 2 {
        let local_bsf = bsf.pop_local_paulis();
        let (t, c, cliff) = search_clifford2q(&bsf, avoid);
        avoid = (cliff.tqs[0], cliff.tqs[1]);
        cliffs_with_locals.push((cliff, local_bsf));
        bsf = t;
    }

    (bsf, cliffs_with_locals)
}

pub fn search_clifford2q(bsf: &BSF, avoid: (usize, usize)) -> (BSF, f64, Operation) {
    let cliff_variants: Vec<Clifford2Q> = Clifford2Q::variants().collect();
    let qubits_with_ops = bsf.qubits_with_ops();
    let pair_indices: Vec<_> = qubits_with_ops
        .iter()
        .combinations(2)
        .filter(|pair| (*pair[0], *pair[1]) != avoid) // && (*pair[1], *pair[0]) != avoid
        .collect();

    // Following is sequential version
    // let mut best_bsf = bsf.clone();
    // let mut min_cost = heuristic_bsf_cost(&best_bsf);
    // let mut best_cliff_idx: usize = 0;
    // let mut best_pair: Vec<usize> = vec![*pair_indices[0][0], *pair_indices[0][1]];
    // for (cliff_idx, cliff) in cliff_variants.iter().enumerate() {
    //     for pair in pair_indices.iter() {
    //         let (i, j) = (*pair[0], *pair[1]);
    //         let bsf_ = bsf.apply_clifford2q(cliff, i, j);
    //         let cost = heuristic_bsf_cost(&bsf_);
    //         if cost < min_cost {
    //             best_bsf = bsf_;
    //             min_cost = cost;
    //             best_cliff_idx = cliff_idx;
    //             best_pair[0] = i;
    //             best_pair[1] = j;
    //         }
    //     }
    // }

    // Following is parallel version
    let results: Vec<_> = cliff_variants
        .par_iter()
        .enumerate()
        .flat_map(|(cliff_idx, cliff)| {
            let bsf_clone = bsf.clone();
            pair_indices
                .par_iter()
                .map(move |pair| {
                    let (i, j) = (*pair[0], *pair[1]);
                    let bsf_ = bsf_clone.apply_clifford2q(cliff, i, j);
                    let cost = heuristic_bsf_cost(&bsf_);
                    (bsf_, cost, cliff_idx, vec![i, j])
                })
                .collect::<Vec<_>>() // Note: Collect into a Vec to avoid borrow issues
        })
        .collect();

    let (best_bsf, min_cost, best_cliff_idx, best_pair) = results
        .into_iter()
        .min_by(|(_, cost1, _, _), (_, cost2, _, _)| cost1.partial_cmp(cost2).unwrap())
        .unwrap();

    let cliff = &cliff_variants[best_cliff_idx];
    (
        best_bsf,
        min_cost,
        Operation::new(Gate::ucg(cliff.p0(), cliff.p1()), best_pair, None),
    )
}

pub fn heuristic_bsf_cost(bsf: &BSF) -> f64 {
    let mut cost = 0.0;
    let which_nonlocal_paulis = bsf.which_nonlocal_paulis();
    if which_nonlocal_paulis.len() > 1 {
        let row_combs: Vec<_> = which_nonlocal_paulis.iter().combinations(2).collect();
        let row_combs: Vec<usize> = row_combs.into_iter().flatten().copied().collect();
        let row_combs = Array2::from_shape_vec((row_combs.len() / 2, 2), row_combs)
            .unwrap()
            .reversed_axes();

        // let row_indices_1 = row_combs.slice(s![0, ..]).to_owned().into_raw_vec_and_offset().0;
        // let row_indices_2 = row_combs.slice(s![1, ..]).to_owned().into_raw_vec_and_offset().0;
        let row_indices_1 = row_combs.row(0).to_vec();
        let row_indices_2 = row_combs.row(1).to_vec();

        let with_ops = bsf.with_ops();
        let x = bsf.x();
        let z = bsf.z();

        cost += (with_ops.select(Axis(0), &row_indices_1)
            | with_ops.select(Axis(0), &row_indices_2))
        .mapv(|v| v as usize)
        .sum() as f64;
        cost += (x.select(Axis(0), &row_indices_1) | x.select(Axis(0), &row_indices_2))
            .mapv(|v| v as usize)
            .sum() as f64
            * 0.5;
        cost += (z.select(Axis(0), &row_indices_1) | z.select(Axis(0), &row_indices_2))
            .mapv(|v| v as usize)
            .sum() as f64
            * 0.5;
    }
    cost += bsf.total_weight() as f64 * bsf.num_nonlocal_paulis().pow(2) as f64;
    cost
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array1;
    use std::time::Instant;

    fn large_paulis() -> Vec<String> {
        [
            "ZYZZXZZYIZ",
            "IYZYXYZIYX",
            "ZXIYIZYXXX",
            "XYIIZXIZXI",
            "XYXIYIXXXZ",
            "YYZZIZYXXI",
            "ZYXIZZZZYZ",
            "XIIIYXZXZY",
            "IZXZZIXIIZ",
            "IYIZYIYZXX",
            "YXZYYIXYYX",
            "XIYZIIIYYI",
            "YZIIIIXZZY",
            "IZZIYZXZIZ",
            "YYZZIXIYIY",
            "XYXYZZIIXY",
            "YYXIYYIXZI",
            "YZXYXZXIXI",
            "IZXZYXXIXX",
            "IXXYXXXZIZ",
        ]
        .iter()
        .map(|&x| x.to_string())
        .collect()
    }

    #[test]
    fn test_bsf_cost() {
        let paulis = large_paulis();
        println!("{:?}", paulis);

        let start = Instant::now();
        let bsf = BSF::new(paulis, None, None);
        let cost = heuristic_bsf_cost(&bsf);
        let duration = start.elapsed();
        assert_eq!(cost, 7192.0);

        println!(
            "qubits_with_nonlocal_ops: {:?}",
            bsf.qubits_with_nonlocal_ops()
        );

        println!("Cost: {}", cost);
        println!("Duration: {:?}", duration);
    }

    #[test]
    fn test_simplify_bsf() {
        let paulis = [
            "ZYZZ", "XZZY", "ZYXY", "ZIYX", "ZXIY", "IZYX", "XXXY", "XYXI", "YIXX", "XZYY", "ZZIZ",
            "YXXI",
        ]
        .iter()
        .map(|&x| x.to_string())
        .collect();
        let coeffs = Array1::from_vec(vec![
            0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2,
        ]);
        let start = Instant::now();
        let bsf = BSF::new(paulis, Some(coeffs), None);
        let (bsf1, cliffs_with_locals) = simplify_bsf(&bsf);
        let duration = start.elapsed();
        println!("now cost is {}", heuristic_bsf_cost(&bsf1));
        println!(
            "paulis: {:?}, coeffs: {}",
            bsf1.paulis(),
            bsf1.coeffs.unwrap()
        );
        for (cliff, local_bsf) in cliffs_with_locals.iter() {
            // println!("cliff: {}, local_bsf: {:?}", cliff, local_bsf);
            if local_bsf.num_paulis() > 0 {
                println!("{}, {:?}", cliff, local_bsf.paulis());
            } else {
                println!("{}", cliff);
            }
        }
        println!("Duration: {:?}", duration);
    }

    #[test]
    #[ignore]
    fn test_simplify_bsf_large() {
        let start = Instant::now();
        let paulis = large_paulis();

        println!("{:?}", paulis);
        let bsf = BSF::new(paulis, None, None);
        let (bsf1, cliffs_with_locals) = simplify_bsf(&bsf);
        let duration = start.elapsed();
        println!("Duration: {:?}", duration);
        println!("now cost is {}", heuristic_bsf_cost(&bsf1));
        println!("paulis: {:?}", bsf1.paulis());
        for (cliff, local_bsf) in cliffs_with_locals.iter() {
            // println!("cliff: {}, local_bsf: {:?}", cliff, local_bsf);
            if local_bsf.num_paulis() > 0 {
                println!("{}, {:?}", cliff, local_bsf.paulis());
            } else {
                println!("{}", cliff);
            }
        }
    }
}
