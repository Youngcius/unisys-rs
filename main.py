import time
import numpy as np
from regulus.models import BSF
from regulus.phoenix import simplification


def random_pauli(n: int) -> str:
    return ''.join(np.random.choice(['X', 'Y', 'Z', 'I'], n, replace=True))

def test_bsf_cost():
    np.random.seed(123)
    num_qubits = 10
    num_paulis = 20

    paulis = [random_pauli(num_qubits) for _ in range(num_paulis)]
    print(paulis)

    start = time.perf_counter()
    tab = BSF(paulis)
    cost = simplification.heuristic_bsf_cost(tab)
    end = time.perf_counter()

    print('qubits_with_nonlocal_ops = ', tab.qubits_with_nonlocal_ops)
    print('Cost: {}'.format(cost))
    print('Time elapsed: {:.2f} us'.format((end - start) * 1e6))

def test_simplify_bsf():
    start = time.perf_counter()
    tab = BSF(
        ['ZYZZ', 'XZZY', 'ZYXY', 'ZIYX', 'ZXIY', 'IZYX', 'XXXY', 'XYXI', 'YIXX', 'XZYY', 'ZZIZ', 'YXXI'],
        [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2]
    )
    tab1, cliffs_with_locals = simplification.simplify_bsf(tab)
    end = time.perf_counter()
    print('Time elapsed: {:.2f} ms'.format((end - start) * 1e3))
    print('now cost is ', simplification.heuristic_bsf_cost(tab1))
    print(tab1.paulis, tab1.coeffs)
    for cliff, local_bsf in cliffs_with_locals:
        if local_bsf.num_paulis > 0:
            print(cliff, local_bsf.paulis)
        else:
            print(cliff)

def test_simplify_bsf_large():
    start = time.time()
    np.random.seed(123)
    num_qubits = 10
    num_paulis = 20

    paulis = [random_pauli(num_qubits) for _ in range(num_paulis)]
    print(paulis)

    tab = BSF(paulis)
    tab1, cliffs_with_locals = simplification.simplify_bsf(tab)
    end = time.time()
    print('Time elapsed: {:.2f} s'.format((end - start)))
    print('now cost is ', simplification.heuristic_bsf_cost(tab1))
    print(tab1.paulis, tab1.coeffs)
    for cliff, local_bsf in cliffs_with_locals:
        if local_bsf.num_paulis > 0:
            print(cliff, local_bsf.paulis)
        else:
            print(cliff)

def test_large_paulis_reconfigure():
    from regulus.models import HamiltonianModel
    start = time.time()
    np.random.seed(123)
    num_qubits = 10
    num_paulis = 20

    paulis = [random_pauli(num_qubits) for _ in range(num_paulis)]
    coeffs = np.repeat(0.1, num_paulis)
    ham = HamiltonianModel(paulis, coeffs)
    config = ham.reconfigure()


    for item in config:
        print(item)

    print('Time elapsed: {:.2f} s'.format(time.time() - start))
    print('number of items in config: ', len(config))

    circ = ham.reconfigure_and_generate_circuit()
    print()
    print('circ.num_gates = ', circ.num_gates)
    print('circ.num_nonlocal_gates = ', circ.num_nonlocal_gates)
    print('circ.depth = ', circ.depth)
    print('circ.depth_nonlocal = ', circ.depth_nonlocal)
    print(circ.gate_stats())
    print('Time elapsed: {:.2f} s'.format(time.time() - start))




if __name__ == '__main__':
    # test_simplify_bsf_large()
    # test_simplify_bsf()

    test_large_paulis_reconfigure()
