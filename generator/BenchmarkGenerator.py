from qiskit import QuantumCircuit
from .supermarqBenchmarks.qaoa_vanilla_proxy import QAOAVanillaProxy
import random
import os
def generateQAOA(qubit_num):
    vqe = QAOAVanillaProxy(qubit_num)
    return QuantumCircuit.from_qasm_str(vqe.circuit().to_qasm())

def generateHamiltonSimulation(qubit_num):
    vqe = hamiltonian_simulation.HamiltonianSimulation(qubit_num)
    return QuantumCircuit.from_qasm_str(vqe.circuit().to_qasm())

def generateQFT(qubit_num):
    circuit = QuantumCircuit(qubit_num)
    for i in range(qubit_num):
        circuit.h(i)
        for j in range(i + 1, qubit_num):
            p = 1
            for k in range(j - i + 1):
                p /= 2
            circuit.cu(0, 0, p, 0, j, i)
    return circuit

def generateRandom(qubit_num):
    circuit = QuantumCircuit(qubit_num)
    gate_num = qubit_num * qubit_num
    for _ in range(gate_num):
        uv = random.sample(range(qubit_num), 2)
        circuit.cx(uv[0], uv[1])
    return circuit

def generateQASM(func, qubit_num_list, name, regenerate = False):
    print("=====================")
    print(f"Check benchmark: {name}")
    path = f"benchmark/{name}"
    if not os.path.exists(path):
        os.makedirs(path)
    for qubit_num in qubit_num_list:
        if regenerate or not os.path.exists(f"{path}/{qubit_num}.qasm"):
            print(f"{path}/{qubit_num}.qasm")
            circuit = func(qubit_num)
            circuit.qasm(filename = f"{path}/{qubit_num}.qasm")
    print(f"Benchmark {name} ready")
    print("=====================")

