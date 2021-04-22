*.itp files *must* conform to the following standards:

1. All heavy-atom *names* in the [ atoms ] directive (column 5) must match the PDB standard. 

2. All heavy atoms in the [ atoms ] directive must appear in the same order as the PDB standard.

3. All hydrogen atoms in the [ atoms ] directive must appear *immediately after* the heavy atom 
to which they are bound.

4. All hydrogen atoms in the [ atoms ] directive must be named (in column 5) "H1", "H2", etc.
*in sequence* beginning from the first hydrogen in the structure. 