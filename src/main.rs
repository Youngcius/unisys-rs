// use ndarray::Array2;
// use ndarray_linalg::svd::SVD;
// use ndarray_linalg::Norm;

// pub fn norm_spec(data: &Array2<f64>) -> f64 {
//     data.svd(false, false).unwrap().1[0]
// }

// fn main() {
//     // 创建一个 2x2 矩阵
//     let matrix = Array2::from_shape_vec((2, 2), vec![1.0, 2.0, 3.0, 4.0]).unwrap();

//     println!("矩阵:\n{}", matrix);
//     println!("norm: {}", matrix.norm());
//     println!("norm_l1: {}", matrix.norm_l1());
//     println!("norm_l2: {}", matrix.norm_l2());
//     println!("norm_max: {}", matrix.norm_max());

//     // 计算 SVD
//     let (u, sigma, vt) = matrix.svd(true, true).unwrap();

//     println!("矩阵 A:\n{}", matrix);
//     println!("左奇异向量 U:\n{}", u.unwrap());
//     println!("奇异值 Σ:\n{}", sigma);
//     println!("右奇异向量 V^T:\n{}", vt.unwrap());

//     let (u, sigma, vt) = matrix.svd(true, true).unwrap();

//     let reconstructed_matrix = u.unwrap().dot(&Array2::from_diag(&sigma)).dot(&vt.unwrap());
//     println!("重构矩阵:\n{}", reconstructed_matrix);

//     println!("{}", matrix.svd(false, false).unwrap().1);
//     // println!("{}", matrix.svd(false, false).unwrap().1);

//     println!("norm_spec is {}", norm_spec(&matrix));
// }

use ndarray::{s, Array1, Array2, Zip};
use ndarray_rand::rand;
use rand::seq::SliceRandom;
use rand::Rng;
use std::time::Instant;
use unisys::{basic::gates::Clifford2Q, models::paulis::BSF, phoenix::simplification};

fn test_for_loop() {
    let mut matrix = Array2::from_shape_vec((10000, 9), (0..90000).collect()).unwrap();
    let start = Instant::now();
    let (mut col_a, mut col_b) = matrix.multi_slice_mut((s![.., 1], s![.., 3]));
    for (a, b) in col_a.iter_mut().zip(col_b.iter_mut()) {
        std::mem::swap(a, b);
    }
    println!("After swapping col-1 and col-3 (for-loop):");
    println!("Duration: {:?}", Instant::now() - start);
}

fn test_for_loop_2() {
    let mut matrix = Array2::from_shape_vec((10000, 9), (0..90000).collect()).unwrap();
    let start = Instant::now();
    for mut row in matrix.outer_iter_mut() {
        // std::mem::swap(&mut row[1], &mut row[3]);
        row.swap(1, 3);
    }
    println!("After swapping col-1 and col-3 (for-loop):");
    println!("Duration: {:?}", Instant::now() - start);
}

fn test_zip() {
    let mut matrix = Array2::from_shape_vec((10000, 9), (0..90000).collect()).unwrap();
    let start = Instant::now();
    let (mut col_a, mut col_b) = matrix.multi_slice_mut((s![.., 1], s![.., 3]));
    Zip::from(col_a)
        .and(col_b)
        .for_each(|a, b| std::mem::swap(a, b));
    println!("After swapping col-1 and col-3 (Zip):");
    println!("Duration: {:?}", Instant::now() - start);
}

fn test_bsf_cost() {
    // def random_pauli(n: int) -> str:
    // return ''.join(np.random.choice(['X', 'Y', 'Z', 'I'], n, replace=True))

    fn random_pauli(n: usize) -> String {
        let mut rng = rand::thread_rng();
        let paulis = ['X', 'Y', 'Z', 'I'];
        (0..n).map(|_| paulis[rng.gen_range(0..4)]).collect()
    }

    // let num_qubits = 10;
    // let num_paulis = 20;
    // let paulis: Vec<String> = (0..num_paulis).map(|_| random_pauli(num_qubits)).collect();

    /*['ZYZZXZZYIZ', 'IYZYXYZIYX', 'ZXIYIZYXXX', 'XYIIZXIZXI', 'XYXIYIXXXZ', 'YYZZIZYXXI', 'ZYXIZZZZYZ', 'XIIIYXZXZY', 'IZXZZIXIIZ', 'IYIZYIYZXX', 'YXZYYIXYYX', 'XIYZIIIYYI', 'YZIIIIXZZY', 'IZZIYZXZIZ', 'YYZZIXIYIY', 'XYXYZZIIXY', 'YYXIYYIXZI', 'YZXYXZXIXI', 'IZXZYXXIXX', 'IXXYXXXZIZ'] */

    let paulis = [
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
    .collect();

    println!("{:?}", paulis);

    let start = Instant::now();
    let bsf = BSF::new(paulis, None, None);
    let cost = simplification::heuristic_bsf_cost(&bsf);
    let duration = start.elapsed();
    assert_eq!(cost, 7192.0);

    println!(
        "qubits_with_nonlocal_ops: {:?}",
        bsf.qubits_with_nonlocal_ops()
    );

    println!("Cost: {}", cost);
    println!("Duration: {:?}", duration);
}

fn test_simplify_bsf_large() {
    let start = Instant::now();
    let paulis = [
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
    .collect();

    println!("{:?}", paulis);
    let bsf = BSF::new(paulis, None, None);
    let (bsf1, cliffs_with_locals) = simplification::simplify_bsf(&bsf);
    let duration = start.elapsed();
    println!("Duration: {:?}", duration);
    println!("now cost is {}", simplification::heuristic_bsf_cost(&bsf1));
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
    let (bsf1, cliffs_with_locals) = simplification::simplify_bsf(&bsf);
    let duration = start.elapsed();
    println!("Duration: {:?}", duration);
    println!("now cost is {}", simplification::heuristic_bsf_cost(&bsf1));
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
}

fn test_vec_to_ndarray() {
    fn vec_to_ndarray(data: Vec<Vec<usize>>) -> Array2<usize> {
        // 将 Vec<Vec<usize>> 转换为 Vec<usize>
        let flattened: Vec<usize> = data.into_iter().flatten().collect();
        // 获取行数
        let rows = flattened.len() / 2;
        // 使用 from_shape_vec 创建 ndarray::Array2
        Array2::from_shape_vec((rows, 2), flattened).unwrap()
    }

    let data = vec![
        vec![0, 1],
        vec![0, 2],
        vec![0, 3],
        vec![0, 4],
        vec![0, 5],
        vec![0, 6],
        vec![0, 7],
        vec![0, 8],
        vec![0, 9],
        vec![0, 10],
        vec![0, 11],
        vec![0, 12],
        vec![0, 13],
        vec![0, 14],
        vec![0, 15],
        vec![0, 16],
        vec![0, 17],
        vec![0, 18],
        vec![0, 19],
        vec![1, 2],
        vec![1, 3],
        vec![1, 4],
        vec![1, 5],
        vec![1, 6],
        vec![1, 7],
        vec![1, 8],
        vec![1, 9],
        vec![1, 10],
        vec![1, 11],
        vec![1, 12],
        vec![1, 13],
        vec![1, 14],
        vec![1, 15],
        vec![1, 16],
        vec![1, 17],
        vec![1, 18],
        vec![1, 19],
        vec![2, 3],
        vec![2, 4],
        vec![2, 5],
        vec![2, 6],
        vec![2, 7],
        vec![2, 8],
        vec![2, 9],
        vec![2, 10],
        vec![2, 11],
        vec![2, 12],
        vec![2, 13],
        vec![2, 14],
        vec![2, 15],
        vec![2, 16],
        vec![2, 17],
        vec![2, 18],
        vec![2, 19],
        vec![3, 4],
        vec![3, 5],
        vec![3, 6],
        vec![3, 7],
        vec![3, 8],
        vec![3, 9],
        vec![3, 10],
        vec![3, 11],
        vec![3, 12],
        vec![3, 13],
        vec![3, 14],
        vec![3, 15],
        vec![3, 16],
        vec![3, 17],
        vec![3, 18],
        vec![3, 19],
        vec![4, 5],
        vec![4, 6],
        vec![4, 7],
        vec![4, 8],
        vec![4, 9],
        vec![4, 10],
        vec![4, 11],
        vec![4, 12],
        vec![4, 13],
        vec![4, 14],
        vec![4, 15],
        vec![4, 16],
        vec![4, 17],
        vec![4, 18],
        vec![4, 19],
        vec![5, 6],
        vec![5, 7],
        vec![5, 8],
        vec![5, 9],
        vec![5, 10],
        vec![5, 11],
        vec![5, 12],
        vec![5, 13],
        vec![5, 14],
        vec![5, 15],
        vec![5, 16],
        vec![5, 17],
        vec![5, 18],
        vec![5, 19],
        vec![6, 7],
        vec![6, 8],
        vec![6, 9],
        vec![6, 10],
        vec![6, 11],
        vec![6, 12],
        vec![6, 13],
        vec![6, 14],
        vec![6, 15],
        vec![6, 16],
        vec![6, 17],
        vec![6, 18],
        vec![6, 19],
        vec![7, 8],
        vec![7, 9],
        vec![7, 10],
        vec![7, 11],
        vec![7, 12],
        vec![7, 13],
        vec![7, 14],
        vec![7, 15],
        vec![7, 16],
        vec![7, 17],
        vec![7, 18],
        vec![7, 19],
        vec![8, 9],
        vec![8, 10],
        vec![8, 11],
        vec![8, 12],
        vec![8, 13],
        vec![8, 14],
        vec![8, 15],
        vec![8, 16],
        vec![8, 17],
        vec![8, 18],
        vec![8, 19],
        vec![9, 10],
        vec![9, 11],
        vec![9, 12],
        vec![9, 13],
        vec![9, 14],
        vec![9, 15],
        vec![9, 16],
        vec![9, 17],
        vec![9, 18],
        vec![9, 19],
        vec![10, 11],
        vec![10, 12],
        vec![10, 13],
        vec![10, 14],
        vec![10, 15],
        vec![10, 16],
        vec![10, 17],
        vec![10, 18],
        vec![10, 19],
        vec![11, 12],
        vec![11, 13],
        vec![11, 14],
        vec![11, 15],
        vec![11, 16],
        vec![11, 17],
        vec![11, 18],
        vec![11, 19],
        vec![12, 13],
        vec![12, 14],
        vec![12, 15],
        vec![12, 16],
        vec![12, 17],
        vec![12, 18],
        vec![12, 19],
        vec![13, 14],
        vec![13, 15],
        vec![13, 16],
        vec![13, 17],
        vec![13, 18],
        vec![13, 19],
        vec![14, 15],
        vec![14, 16],
        vec![14, 17],
        vec![14, 18],
        vec![14, 19],
        vec![15, 16],
        vec![15, 17],
        vec![15, 18],
        vec![15, 19],
        vec![16, 17],
        vec![16, 18],
        vec![16, 19],
        vec![17, 18],
        vec![17, 19],
        vec![18, 19],
    ];
    let array = vec_to_ndarray(data);

    println!("{:?}", array);
}

fn test_aaa() {
    let data = vec![vec![1, 2], vec![3, 4], vec![5, 6]];
    let rows = 3;
    let cols = 2;

    // 将 Vec 转换为 Array2

    let array = Array2::from_shape_vec((3, 2), data.into_iter().flatten().collect()).unwrap();
    println!("{:?}", array);
}

fn main() {
    // 行列交换 by for-loop
    test_for_loop();

    // 行列交换 by for-loop 2
    test_for_loop_2();

    // 行列交换 by Zip
    test_zip();

    // test_vec_to_ndarray();

    // test_aaa();
}
