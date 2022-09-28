from __future__ import annotations
import numpy as np
import pandas as pd

from fractions import Fraction
from math import gcd
from qiskit.circuit import ClassicalRegister, QuantumCircuit, QuantumRegister
from qiskit.circuit.controlledgate import ControlledGate
from qiskit.circuit.library import QFT


class ShorAlgorithm:
    """Shor's algorithm class"""
    
    def __init__(self) -> ShorAlgorithm:
        pass
    
    def _initialize_circuit(self, num_qubits: int, num_qubits_qft: int) -> QuantumCircuit:
        """Initialize Shor's algorithm circuit.

        Args:
            number (int): The number that we want to factorize.
            num_qubits_qft (int): The number of qubits for QFT inverse circuit.

        Returns:
            QuantumCircuit: A initialized Shor's algorithm circuit.
        """
        
        qubits_qft = QuantumRegister(size=num_qubits_qft)
        qubits = QuantumRegister(size=num_qubits)
        bits = ClassicalRegister(size=num_qubits_qft)
        
        qc = QuantumCircuit(qubits_qft, qubits, bits)
        
        qc.h(qubit=qubits_qft)
        qc.x(qubit=qubits[-1])
        qc.barrier()
        
        return qc
    
    def _c_amodN(self, number: int, a: int, power: int) -> ControlledGate:
        """Return a ControlledGate for QPE.

        Args:
            number (int): The number that we want to factorize.
            a (int): An integer number that must belong to set of numbers defined by gcd with the input number.
            power (int): The number that will be applied as power of a.

        Raises:
            ValueError: If a doest not belong to finite field defined by input number.

        Returns:
            ControlledGate: The controlled gate for QPE in Shor's algorithm.
        """
        
        gcd_results = [(i,gcd(i,number)) for i in range(2,number)]
        list_of_numbers = [tp[0] for tp in gcd_results if tp[1] == 1]
        
        if a not in list_of_numbers:
            raise ValueError("a is not in list of numbers!")
        
        for _ in range(power):
            num_qubits = len(bin(number)[2:])
            U = QuantumCircuit(num_qubits)
            a_bin = bin(a)[2:]
            
            if len(a_bin) < num_qubits:
                new_a_bin = ''
                for _ in range(num_qubits-len(a_bin)):
                    new_a_bin += '0'
                new_a_bin += a_bin
            else:
                new_a_bin = a_bin
            
            for i in range(num_qubits-1,num_qubits-len(a_bin), -1):
                U.swap(qubit1=i, qubit2=i-1)
            
            aux = 0
            for j, char in enumerate(new_a_bin):
                if char == '1':
                    aux = j
                    break
                
            for k, char in enumerate(new_a_bin):
                if char == '1' and k != aux:
                    U.x(qubit=k)
        
        U = U.to_gate()
        U.name = f"${a}^{power} mod {number}$"
        controlled_U = U.control()
                    
        return controlled_U
        
    def quantum_circuit(self, number: int, a: int, num_qubits_qft: int) -> QuantumCircuit:
        """Return a full Shor's algorithm circuit.

        Args:
            number (int): The number that we want to factorize.
            a (int): An integer number that must belong to set of numbers defined by gcd with the input number.
            num_qubits_qft (int): The number of qubits for QFT inverse circuit.

        Returns:
            QuantumCircuit: A full quantum circuit for Shor's algorithm.
        """
        
        num_qubits = len(bin(number)[2:])
        qc = self._initialize_circuit(num_qubits=num_qubits, num_qubits_qft=num_qubits_qft)
        for i in range(num_qubits_qft):
            controlled_U = self._c_amodN(number=number, a=a, power=2**i)
            qc.append(controlled_U, [qc.qregs[0][i]] + [qc.qregs[1][j] for j in range(qc.qregs[1].size)])
        qft_dagger = QFT(num_qubits=qc.qregs[0].size).inverse()
        qc.append(qft_dagger, qc.qregs[0])
        qc.measure(qubit=qc.qregs[0], cbit=qc.cregs[0])
        
        return qc
    
    def get_number_prime_factors(self, number: int, counts: dict) -> tuple:
        """Return a DataFrame and a Set with information about prime factors.

        Args:
            number (int): The number that we want to factorize.
            counts (dict): The result of the experiments done with
                           Shor's algorithm quantum circuit.

        Returns:
            tuple: A DataFrame with information about prime factors
                   and a Set with prime factors found.
        """
        
        guesses_list = []
        rows = []
        columns = ["register_output", "phase", "fraction", "guess_for_r"]
        keys_list = list(counts.keys())
        num_qubits = len(keys_list[0])
        
        for output in counts:
            decimal = int(output, 2)
            phase = decimal/(2**num_qubits)
            frac = Fraction(phase).limit_denominator(15)
            rows.append([f"{output}(bin) = {decimal:>3}(dec)", 
                        phase, 
                        f"{frac.numerator}/{frac.denominator}",
                        frac.denominator])
            
        df = pd.DataFrame(rows, columns=columns)

        for phase, r in zip(df.phase.values, df.guess_for_r.values):
            if not np.isclose(0.0, phase):
                if r%2 == 0:
                    guesses = [gcd(2**(r//2)-1, number), gcd(2**(r//2)+1, number)]
                    for guess in guesses:
                        if guess not in [1,number] and (number % guess) == 0:
                            guesses_list.append(guess)
                        
        guesses_set = set(guesses_list)
        
        return df, guesses_set