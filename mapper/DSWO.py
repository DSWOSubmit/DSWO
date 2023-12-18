import os.path

from qiskit import QuantumCircuit
from qiskit.circuit.library.standard_gates import SwapGate

import ctypes
from ctypes import *


def DSWO(circuit : QuantumCircuit, couplingMap, initialLayout, limited = 15, name = "DSWO0.so") -> QuantumCircuit:
    qubit_num = circuit.num_qubits
    gate_num = circuit.size()

    LAYOUT = [False for _ in range(qubit_num * qubit_num)]
    for edge in couplingMap:
        LAYOUT[edge[0] * qubit_num + edge[1]] = True
        LAYOUT[edge[1] * qubit_num + edge[0]] = True

    gate_passin = []
    gateid = 0
    gates_mapping = {}

    for gate in circuit.data:
        if len(gate.qubits) == 1:
            gate_passin.append(gate.qubits[0].index)
            gate_passin.append(-1)
        elif len(gate.qubits) == 2:
            gate_passin.append(gate.qubits[0].index)
            gate_passin.append(gate.qubits[1].index)
        else:
            assert False
        gateid += 1
        gate_passin.append(gateid)
        gates_mapping[gateid] = gate

    if not os.path.exists(f"./mapper/{name}"):
        command = f"g++ -std=c++17 -O2 -o ./mapper/{name} -shared -fPIC ./mapper/DSWO.cpp"
        os.system(command)

    ll = ctypes.cdll.LoadLibrary
    lib = ll(f"./mapper/{name}")

    # int *gateContent, bool *layout, int *l2p, int gateLen, int qubitNumber
    lib.DSWO.argtypes = [POINTER(c_int * (3 * gateid)), POINTER(c_bool * (qubit_num * qubit_num)),
                              POINTER(c_int * qubit_num),
                                c_int, c_int, c_int]
    # lib.DSWO.restype = ctypes.POINTER(ctypes.c_int * (gateid * 1000))
    lib.DSWO.restype = c_void_p
    # print(LAYOUT, qubit_num, gate_passin)
    gate_passin = (c_int * (3 * gateid))(*(i for i in gate_passin))
    LAYOUT = (c_bool * (qubit_num * qubit_num))(*(i for i in LAYOUT))
    initialLayout = (c_int * qubit_num)(*(i for i in initialLayout))

    data_out_str = lib.DSWO(gate_passin, LAYOUT, initialLayout,
                              gate_num, qubit_num, limited)

    data_out = ctypes.cast(data_out_str, ctypes.c_char_p).value.decode()

    lib.freeDSWO.argtypes = [c_void_p]
    lib.freeDSWO.restype = None
    lib.freeDSWO(data_out_str)

    new_circuit = QuantumCircuit(qubit_num)

    anastr = data_out.split(".")

    gate_num = int(anastr[0])
    for i in range(gate_num):
        gate = anastr[i + 1]
        uvid = gate.split(",")
        if len(uvid) == 3:
            g = gates_mapping[int(uvid[0])]
            u = int(uvid[1])
            v = int(uvid[2])
            _g = g.copy()
            if v == -1:
                _g.qubits = [new_circuit.qubits[u]]
            else:
                _g.qubits = [new_circuit.qubits[u], new_circuit.qubits[v]]
            # print(u, v, _g.qubits)
            new_circuit.append(_g)
        else:
            new_circuit.swap(int(uvid[0]), int(uvid[1]))

    return new_circuit
