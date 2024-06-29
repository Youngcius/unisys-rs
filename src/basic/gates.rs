pub struct Gate {
    pub name: String,
    pub tqs: Vec<usize>,
    pub cqs: Vec<usize>,
    pub n_qubits: usize,
    pub angle: f64,
    pub angles: Vec<f64>,
    pub params: Vec<f64>,
    pub exponent: f64,
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
