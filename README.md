# Shor's algorithm

## Description:

A implementation of Shor's algorithm that works for an arbitrary number. 

## How to use:

### Installing requirements:

`git clone https://github.com/jvscursulim/shor_algorithm`

`pip install -r requirements.txt`

### Code example:

Factoring `number=15` and using `a=2` as an element from the finite field defined by 15.

```python
from qiskit import execute
from qiskit_aer import Aer

from shor import ShorAlgorithm

shor = ShorAlgorithm()

qc = shor.quantum_circuit(number=15, a=2, num_qubits_qft=2)

backend = Aer.get_backend("qasm_simulator")
counts = execute(qc, backend=backend, shots=8192).result().get_counts()

df, factors = shor.get_number_prime_factors(number=15, counts=counts)

print(factors)
```

Output: `{3,5}`

### Quantum Circuit of code example:

![image](figures/qc_example.png)

## References:

[Qiskit Textbook - Shor's algorithm](https://qiskit.org/textbook/ch-algorithms/shor.html)
