// import numpy as np
// import rustworkx as rx
// from typing import List, Dict, Tuple
// from functools import reduce
// from operator import add
// from numpy import ndarray
// from regulus.basic.gates import Gate, SWAPGate, SWAP
// from regulus.basic.circuits import Circuit
// from regulus.utils.passes import obtain_front_layer
// from regulus.utils.arch import gene_init_mapping, is_executable, update_mapping, unify_mapped_circuit
// from regulus.utils.arch import obtain_logical_neighbors
// from rich.console import Console

use crate::basic::{gates, gates::Gate, gates::GateType};
use crate::basic::{operations, operations::Operation};
use crate::utils::passes;

const INIT_DECAY: f64 = 1.0;
const DECAY_STEP: f64 = 0.001;
const NUM_SEARCHES_TO_RESET: usize = 5;
const EXT_WEIGHT: f64 = 0.5;
const EXT_SIZE: usize = 20;

/*
def sabre_search_one_pass(circ: Circuit, device: rx.PyGraph, init_mapping: Dict[int, int] = None,
                          return_circ_with_swaps: bool = False,
                          gene_init_mapping_type: str = 'random',
                          seed: int = None) -> Tuple[Circuit, Dict[int, int], Dict[int, int]]:
    assert _has_decomposed_completely(circ), "The input circuit should be decomposed into 1Q + 2Q gates completely"
    qubits = circ.qubits
    rng = np.random.default_rng(seed)
    circ = circ.clone()
    mappings = []

    INIT_DECAY_PARAMS = {q: INIT_DECAY for q in qubits}
    decay_params = INIT_DECAY_PARAMS.copy()

    # find the front layer first
    front_layer = obtain_front_layer_from_circuit(circ)

    # begin SWAP searching until front_layer is empty
    dist_mat = rx.floyd_warshall_numpy(device)
    mappings.append(gene_init_mapping(circ, device, gene_init_mapping_type) if init_mapping is None else init_mapping)
    circ_with_swaps = Circuit()
    exe_gates = []
    num_searches = 0
    while front_layer:
        exe_gates.clear()

        for g in front_layer:
            if is_executable(g, mappings[-1], dist_mat):
                exe_gates.append(g)

        if exe_gates:
            for g in exe_gates:
                circ.remove(g)
            circ_with_swaps.append(*exe_gates)
            front_layer = obtain_front_layer_from_circuit(circ)

        else:  # find suitable SWAP gates
            # reset decay_params every 5 rounds
            num_searches += 1
            if num_searches % NUM_SEARCHES_TO_RESET == 0:
                decay_params.update(INIT_DECAY_PARAMS)

            swap_candidates = obtain_swap_candidates(front_layer, mappings[-1], device)

            # get extending_layer including only 2Q gates with maximal number EXT_SIZE
            extending_layer = []
            num = 0
            for g in circ:
                if g in front_layer:
                    continue
                if g.num_qregs == 2:
                    extending_layer.append(g)
                    num += 1
                if num == EXT_SIZE:
                    break


            scores = np.array([heuristic_score(front_layer, extending_layer,
                                               mappings[-1], swap, dist_mat, decay_params) for swap in swap_candidates])

            # find the SWAP with minimal score
            swap = swap_candidates[rng.choice(np.where(np.isclose(scores, scores.min()))[0])]
            circ_with_swaps.append(swap)
            decay_params[swap.tqs[0]] += DECAY_STEP
            decay_params[swap.tqs[1]] += DECAY_STEP

            mappings.append(update_mapping(mappings[-1], swap))
    if return_circ_with_swaps:
        return circ_with_swaps, mappings
    return unify_mapped_circuit(circ_with_swaps, mappings), mappings[0], mappings[-1]
 */
