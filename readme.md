In this work, we use LP bounds to rule out the existence of quantum LDPC codes.

The code is meant to be executed on matlab. To run, see the file scan_ldpc.

It searches for simple types of LDPC codes where a code C is associated with both X- and Z-type operators.
Furthermore, the code is assumed to be regular.
Cperp maps to the stabilizers of both.
In the file scan_ldpc, you can set dv, dc, the degrees of the qubits and checks respectively.
The program will try to rule out codes with block sizes between n_min and n_max that cannot encode even a single qubit (tweaked using the k parameter).
Note that if the program does NOT rule out a code with some parameters, it does not mean the code exists.

This is joint work with Jean-Pierre Tillich.
