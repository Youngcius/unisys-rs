use std::collections::HashMap;
use rustworkx_core::petgraph::graph::UnGraph;

use crate::basic::{circuits::Circuit, gates::Gate, operations::Operation};


pub fn mirroring_sabre(
    circ: &Circuit,
    device: & UnGraph<usize, ()>,
    num_pass_periods: usize,
    init_mapping: &HashMap<usize, usize>
) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
    panic!("Not implemented")
}

pub fn mirroring_sabre_one_pass(
    circ: &Circuit,
    device: &UnGraph<usize, ()>,
    init_mapping: &HashMap<usize, usize>
) -> (Circuit, HashMap<usize, usize>, HashMap<usize, usize>) {
    panic!("Not implemented")
}


#[cfg(test)]
mod tests {
    use rustworkx_core::{connectivity::connected_components, petgraph::{graph::{NodeIndex, UnGraph}, prelude::StableUnGraph}};

    use super::*;

    #[test]
    fn test_basic_rustworkx_core() {
        // let mut graph = UnGraph::<usize, ()>::default();
        let mut graph = StableUnGraph::<usize, ()>::default();

        let a = graph.add_node(10);
        let b = graph.add_node(20);
        let c = graph.add_node(30);
        graph.add_edge(a, b, ());
        graph.add_edge(b, c, ());
        graph.add_edge(c, a, ());
        assert_eq!(graph.node_count(), 3);
        assert_eq!(graph.edge_count(), 3);

        // drop b node
        graph.remove_node(b);
        assert_eq!(graph.node_count(), 2);

        for idx in graph.node_indices() {
            println!("{:?}", NodeIndex::from(idx.index()));
            println!("{:?}, {}", idx.index(), graph.node_weight(idx).unwrap());
        }
        println!();
        // println!("{:?}", graph.node_indices().collect::<Vec<_>>());
        // println!("{:?}", graph.node_weights().collect::<Vec<_>>());
        // println!("{:?}", graph.edge_indices().collect::<Vec<_>>());
        let components = connected_components(&graph);
        println!("{:?}", components);

        let mut subgraph = graph.clone();
        subgraph.retain_nodes(|_, idx| idx.index() < 2);
        println!("{:?}", subgraph.node_indices().collect::<Vec<_>>());

    }

}
