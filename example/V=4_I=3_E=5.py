from GCPToHUBO import *

if __name__ == '__main__':
    vertex_comb = x(1, 2) + x(1, 4) + x(2, 3) + x(2, 4) + x(3, 4)
    objfunc, costfunc, penaltyfunc = GCPToHUBO(vertex_comb, 4, 3, binary_desc=True)
    qubo_qubit, hubo_qubit, hubo_or_qubit = GCP_Qubit_Count(5, 4, 3)
    #qubo_qubit, hubo_qubit, hubo_or_qubit = GCP_Qubit_Count(len(vertex_comb), 4, 3)
    qubo_CNOT, hubo_asc_CNOT, hubo_pf_CNOT, hubo_or_CNOT = GCP_CNOT_Count(5, 4, 3)
    qubo_TGate, hubo_asc_TGate, hubo_pf_TGate, hubo_or_TGate = GCP_TGate_Count(5, 4, 3)
    qubo_totalTGate, hubo_asc_totalTGate, hubo_pf_totalTGate, hubo_or_totalTGate = GCP_TotalTGate_Count(5, 4, 3)