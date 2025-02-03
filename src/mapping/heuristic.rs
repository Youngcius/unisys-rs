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
