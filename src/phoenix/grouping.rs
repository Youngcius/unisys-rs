use crate::models::paulis;
use indexmap::IndexMap;
use std::collections::HashMap;

fn least_overlap(indices: &Vec<Vec<usize>>, existing_indices: &Vec<Vec<usize>>) -> Vec<usize> {
    let mut overlaps: Vec<usize> = Vec::new();
    for idx in indices.iter() {
        let mut overlap = 0;
        for eidx in existing_indices.iter() {
            overlap += idx.iter().filter(|i| eidx.contains(i)).count();
        }
        overlaps.push(overlap);
    }
    let min_overlap = *overlaps.iter().min().unwrap();
    indices[overlaps.iter().position(|&x| x == min_overlap).unwrap()].clone()
}

pub fn group_paulis(paulis: &Vec<String>) -> Vec<(Vec<usize>, Vec<String>)> {
    assert!(
        paulis.len()
            == paulis
                .iter()
                .collect::<std::collections::HashSet<&String>>()
                .len(),
        "Pauli strings must be unique"
    );

    let mut groups: HashMap<Vec<usize>, Vec<String>> = HashMap::new();

    let nontrivial = paulis
        .iter()
        .map(|pauli| paulis::nontrivial_indices(pauli))
        .collect::<Vec<Vec<usize>>>();

    for (idx, pauli) in nontrivial.iter().zip(paulis.iter()) {
        if !groups.contains_key(idx) {
            groups.insert(idx.clone(), vec![pauli.clone()]);
        } else {
            groups.get_mut(idx).unwrap().push(pauli.clone());
        }
    }

    // sort "groups" according to the length of the keys and keys themselves
    let mut sorted_groups: Vec<(Vec<usize>, Vec<String>)> = groups.clone().into_iter().collect();
    sorted_groups.sort_by(|a, b| {
        if a.0.len() == b.0.len() {
            a.0.cmp(&b.0)
        } else {
            b.0.len().cmp(&a.0.len()) // a.len=5, b.len=4
        }
    });
    let mut groups = sorted_groups;

    // reorder items to reduce overall length when organizing as circuit
    let mut groups_on_length: IndexMap<usize, Vec<(Vec<usize>, Vec<String>)>> = IndexMap::new(); // {length: [(idx: [pauli, ...]), ...]}}
    for (idx, paulis) in groups.iter() {
        let length = idx.len();
        if !groups_on_length.contains_key(&length) {
            let mut inner_vec: Vec<(Vec<usize>, Vec<String>)> = Vec::new();
            inner_vec.push((idx.clone(), paulis.clone()));
            groups_on_length.insert(length, inner_vec);
        } else {
            groups_on_length
                .get_mut(&length)
                .unwrap()
                .push((idx.clone(), paulis.clone()));
        }
    }

    groups.clear();

    /*
    Python code:
        for equal_len_groups in groups_on_length.values():
        selected_indices = []
        while equal_len_groups:
            idx = least_overlap(list(equal_len_groups.keys()), selected_indices)
            selected_indices.append(idx)
            groups[idx] = equal_len_groups.pop(idx)
     */
    for equal_len_groups in groups_on_length.values_mut() {
        let mut selected_indices: Vec<Vec<usize>> = Vec::new();
        while !equal_len_groups.is_empty() {
            let idx = least_overlap(
                &equal_len_groups
                    .iter()
                    .map(|(idx, _)| idx.clone())
                    .collect::<Vec<_>>(),
                &selected_indices,
            );
            selected_indices.push(idx.clone());
            let pos = equal_len_groups
                .iter()
                .position(|(i, _)| i == &idx)
                .unwrap();
            groups.push(equal_len_groups.remove(pos));
        }
    }
    groups
}

pub fn group_paulis_and_coeffs(
    paulis: &Vec<String>,
    coeffs: &Vec<f64>,
) -> Vec<(Vec<usize>, (Vec<String>, Vec<f64>))> {
    assert!(
        paulis.len() == coeffs.len(),
        "The number of Pauli strings must be equal to the number of coefficients"
    );
    let pauli_groups = group_paulis(paulis);
    let mut groups: Vec<(Vec<usize>, (Vec<String>, Vec<f64>))> = Vec::new();
    for (idx, paulis) in pauli_groups.iter() {
        let mut group_coeffs: Vec<f64> = Vec::new();
        for pauli in paulis {
            let pos = paulis.iter().position(|p| p == pauli).unwrap();
            group_coeffs.push(coeffs[pos]);
        }
        groups.push((idx.clone(), (paulis.clone(), group_coeffs)));
    }
    groups
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_group_paulis_1() {
        /*
        E.g.,

            ['XXIII', 'YYIII', 'ZZIII', 'IXXII', 'IYYII', 'IZZII', 'IIXXI', 'IIYYI', 'IIZZI', 'IIIXX', 'IIIYY', 'IIIZZ', 'ZIIII', 'IZIII', 'IIZII', 'IIIZI', 'IIIIZ']

        will be grouped as

           {(0, 1): ['XXIII', 'YYIII', 'ZZIII'],
            (2, 3): ['IIXXI', 'IIYYI', 'IIZZI'],
            (3, 4): ['IIIXX', 'IIIYY', 'IIIZZ'],
            (1, 2): ['IXXII', 'IYYII', 'IZZII'],
            (0,): ['ZIIII'],
            (1,): ['IZIII'],
            (2,): ['IIZII'],
            (3,): ['IIIZI'],
            (4,): ['IIIIZ']}
        */

        let paulis = vec![
            "XXIII".to_string(),
            "YYIII".to_string(),
            "ZZIII".to_string(),
            "IXXII".to_string(),
            "IYYII".to_string(),
            "IZZII".to_string(),
            "IIXXI".to_string(),
            "IIYYI".to_string(),
            "IIZZI".to_string(),
            "IIIXX".to_string(),
            "IIIYY".to_string(),
            "IIIZZ".to_string(),
            "ZIIII".to_string(),
            "IZIII".to_string(),
            "IIZII".to_string(),
            "IIIZI".to_string(),
            "IIIIZ".to_string(),
        ];
        let groups = group_paulis(&paulis);
        // print groups
        for (idx, paulis) in groups.iter() {
            println!("{:?}: {:?}", idx, paulis);
        }
    }

    #[test]
    fn test_group_paulis_and_coeffs_1() {
        let paulis = vec![
            "XXIII".to_string(),
            "YYIII".to_string(),
            "ZZIII".to_string(),
            "IXXII".to_string(),
            "IYYII".to_string(),
            "IZZII".to_string(),
            "IIXXI".to_string(),
            "IIYYI".to_string(),
            "IIZZI".to_string(),
            "IIIXX".to_string(),
            "IIIYY".to_string(),
            "IIIZZ".to_string(),
            "ZIIII".to_string(),
            "IZIII".to_string(),
            "IIZII".to_string(),
            "IIIZI".to_string(),
            "IIIIZ".to_string(),
        ];
        let coeffs = vec![
            0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 0.1, 0.2, 0.3, 1.0, 1.0, 1.0, 1.0, 1.0,
        ];
        let groups = group_paulis_and_coeffs(&paulis, &coeffs);
        // print groups and coeffs
        for (idx, (paulis, coeffs)) in groups.iter() {
            println!("{:?}: {:?}, {:?}", idx, paulis, coeffs);
        }
    }
    #[test]
    fn test_group_paulis_2() {
        /*
        E.g.,
            ['ZYZZ', 'XZZY', 'ZYXY', 'ZIYX', 'ZXIY', 'IZYX', 'XXXY', 'XYXI', 'YIXX', 'XZYY', 'ZZIZ', 'YXXI']

        will be grouped as

           {(0, 1, 2, 3): ['ZYZZ', 'XZZY', 'ZYXY', 'XXXY', 'XZYY'],
            (0, 1, 2): ['XYXI', 'YXXI'],
            (0, 1, 3): ['ZXIY', 'ZZIZ'],
            (0, 2, 3): ['ZIYX', 'YIXX'],
            (1, 2, 3): ['IZYX']}
        */
        let paulis = vec![
            "ZYZZ".to_string(),
            "XZZY".to_string(),
            "ZYXY".to_string(),
            "ZIYX".to_string(),
            "ZXIY".to_string(),
            "IZYX".to_string(),
            "XXXY".to_string(),
            "XYXI".to_string(),
            "YIXX".to_string(),
            "XZYY".to_string(),
            "ZZIZ".to_string(),
            "YXXI".to_string(),
        ];
        let groups = group_paulis(&paulis);
        // print groups
        for (idx, paulis) in groups.iter() {
            println!("{:?}: {:?}", idx, paulis);
        }
    }
}
