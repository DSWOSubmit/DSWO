import pandas as pd
import sys
import json

from generator.BenchmarkGenerator import *
from generator.LayoutGenerator import *

from qiskit.transpiler.passes import SabreSwap, SabreLayout
from qiskit.transpiler import CouplingMap, PassManager


from mapper.DSWO import DSWO

import numpy as np

class EXPResult:
    def __init__(self):
        self.SABRESize = []
        self.SABREDepth = []
        self.DSWOSize = []
        self.DSWODepth = []

def run_round(index, couple_graph, qiskit_qc, limitL, qasm_out_path):
    coupling_map = CouplingMap(couple_graph)
    coupling_map.make_symmetric()

    passmanager = PassManager()
    sabre_layout = SabreLayout(coupling_map, max_iterations=1, skip_routing=True)
    passmanager.append(sabre_layout)

    passmanager.run(qiskit_qc)

    SLR = sabre_layout.property_set["layout"]

    # rebuild circuit with coupling_map size
    new_circuit = QuantumCircuit(coupling_map.size())
    for gate in qiskit_qc.data:
        if len(gate.clbits) > 0:
            continue
        _g = gate.copy()
        # _g.qubits = [new_circuit.qubits[qiskit_qc.find_bit(qubit).index] for qubit in _g.qubits]
        # print([SLR[q] for q in _g.qubits])
        _g.qubits = [SLR[q] for q in _g.qubits]
        new_circuit.append(_g)

    # run SABRE
    passmanager = PassManager()
    sabre_swap = SabreSwap(coupling_map, heuristic="decay", trials=1)
    passmanager.append(sabre_swap)
    sabre_result = passmanager.run(new_circuit)

    # run DSWO
    initial_mapping = [qubit for qubit in range(new_circuit.num_qubits)]
    dswo_result = DSWO(new_circuit, couple_graph, initial_mapping, limitL, f"DSWP{index}.so")

    dswo_result.qasm(filename = f"{qasm_out_path}/{qiskit_qc.num_qubits}_{index}.qasm")
    return sabre_result.size(), sabre_result.depth(), dswo_result.size(), dswo_result.depth()

def run_serveral_round(couple_graph, qiskit_qc, limitL, _round, qasm_out_path):
    result = EXPResult()

    for index in range(_round):
        run_data = run_round(index, couple_graph, qiskit_qc, limitL, qasm_out_path)
        result.SABRESize.append(run_data[0])
        result.SABREDepth.append(run_data[1])
        result.DSWOSize.append(run_data[2])
        result.DSWODepth.append(run_data[3])

    return np.min(result.SABRESize), np.min(result.SABREDepth), np.min(result.DSWOSize), np.min(result.DSWODepth)

def runEXP(couple_graph, graph_name, benchmark_set_path, limitL_set, exp_out_path, fit_scale, _round):
    # run exp

    device_qubit_num = couple_graph[0]
    device_graph = couple_graph[1]

    qasm_out_path = f"{exp_out_path}{graph_name}{device_qubit_num}"
    exp_out_path = f"{qasm_out_path}.csv"
    print(f"run exp in {graph_name}, qasm will output in {qasm_out_path}")
    if not os.path.exists(f"{qasm_out_path}"):
        os.makedirs(f"{qasm_out_path}")

    # couple_graph, benchmark_qasm, limitL, _round
    df = pd.DataFrame(columns=[
        "qubit_number", "limitL", "origin_size", "origin_depth", "SABRE_size", "SABRE_depth", "DSWO_size", "DSWP_depth"
    ])
    for file in os.listdir(benchmark_set_path):
        # print(file, not file.startswith(".") and file.endswith(".qasm"))
        if not file.startswith(".") and file.endswith(".qasm"):
            qiskit_qc = QuantumCircuit.from_qasm_file(f"{benchmark_set_path}{file}")
            # print(fit_scale, qiskit_qc.num_qubits, device_qubit_num, benchmark_set_path)
            if fit_scale and qiskit_qc.num_qubits != device_qubit_num:
                continue
            if qiskit_qc.num_qubits > device_qubit_num:
                continue
            for limitL in limitL_set:
                line = [qiskit_qc.num_qubits, limitL, qiskit_qc.size(), qiskit_qc.depth()]
                ans = run_serveral_round(device_graph, qiskit_qc, limitL, _round, qasm_out_path)
                line.extend(list(ans))
                df.loc[len(df)] = line
                df.to_csv(exp_out_path, index=False)

    print("===============")

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Expected: python main.py <exp_config.json>")
        sys.exit(1)

    with open(sys.argv[1]) as f:
        config = json.load(f)

    name = config["benchmark"]
    benchmark_path = config["benchmark_path"]
    out_path = config["out_path"]
    fit_scale = config["fit_scale"]
    limitLength = config["limitLength"]
    regenerate = config["regenerate"]
    _round = config["round"]
    exp_name = config["exp_name"]

    print("start run exp", exp_name)

    if not os.path.exists(out_path):
        os.makedirs(out_path)

    # generate benchmark
    if exp_name == "test":
        qubit_list = [5, 6, 7, 8, 9, 10]
    elif fit_scale:
        qubit_list = [209, 311, 433, 575, 737, 919]
    else:
        qubit_list = [100, 127, 200, 250, 300, 350, 400, 433]

    generate_func = generateQFT
    custom_benchmark = False
    if name == "qaoa":
        generate_func = generateQAOA
    elif name == "qft" or name == "qfttest" or name == "qftfit":
        generate_func = generateQFT
    elif name == "random" or name == "randomfit":
        generate_func = generateRandom
    elif name == "hamilton":
        generate_func = generateHamiltonSimulation
    else:
        # custom benchmark
        custom_benchmark = True

    if not custom_benchmark:
        generateQASM(generate_func, qubit_list, name, regenerate)

    if exp_name == "test":
        layouts = [[generateGridLayout(5, 5), "Grid"]]
    elif fit_scale:
        layouts = []
        for i in range(18, 39, 4):
            layouts.append([generateIBMQLayout(i), "Heavy Hex"])
    else:
        layouts = [[generateGridLayout(20, 20), "Grid"], [generateIBMQLayout(26), "Heavy Hex"]]

    for layout in layouts:
        runEXP(layout[0], layout[1], benchmark_path, limitLength, out_path, fit_scale, _round)
