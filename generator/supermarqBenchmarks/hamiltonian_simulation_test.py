from .benchmarks.hamiltonian_simulation import HamiltonianSimulation


def test_hamiltonian_simulation_circuit() -> None:
    hs = HamiltonianSimulation(4, 1, 1)
    assert len(hs.circuit().all_qubits()) == 4


def test_hamiltonian_simulation_score() -> None:
    hs = HamiltonianSimulation(4, 1, 1)
    assert hs._average_magnetization({"1111": 1}, 1) == -1.0
    assert hs._average_magnetization({"0000": 1}, 1) == 1.0
    assert hs.score(hs._get_ideal_counts(hs.circuit())) > 0.99
