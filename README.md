# Generator of rankmap for 4-dim simulation such as LQCD

Copyright(c) 2020-2023, Issaku Kanamori <kanamori-i@riken.jp>

License: GPL v3

## Description

This program generates a 4-dim rankmap for supercomputer Fugaku (or any other FX-1000 machine that uses Tofu network), as a
 (3 dim node torus) x  (4 intra-node ranks).
The number of intra-node ranks must be 4.

Edit the following files if needed
  calc_rankid.c  : defines the map from 4 dim process coordinate to 1 dim rank id (default: lexical)
  config.h       : defines the filename of the rankmap

process coordinate (p1,p2,p3,p4)  [process size: (P1,P2,P3,P4)]
with the following rankid gives the optimal rankmap 
  rankmap_4d_lex            rankid = p1 + p2 P1 + p3 P1 P2 + p4 P1 P2 P3 ( = Fortran style)
  rankmap_4d_lex_reversed   rankid = p4 + p3 P4 + p2 P4 P3 + p1 P4 P3 P2 (= C style)

It requires Fujitsu MPI.

## Usage (for rankmap_4d_general)

Suppose that your application (=./a.out) uses 4-dim lexical process rankmap, 
of which process size is (P1,P2,P3,P4), and intra-node 4 processes are divided 
as (in1, in2, in3, in4)  [in1 * in2 * in3 * in4 = 4].

```
#PJM --rsc-list "node=PP1xPP2xPP3"
mpirun ./rankmap_4d_general_lex P1 P2 P3 P4 in1 in2 in3 in4  # or mpirun ./rankmap_4d_general_lex_reversed ....
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

Here, (in1,in2,in3,in4) must divide (P1,P2,P3,P4).
Furthermore, one of (P1/in1, P2/in2, P3/in3, P4/in4) must be 1.
Let us suppose P3/in3 is 1, and the remaining 3-dim (P1/in1 ,P2/in2, P4/in4)
should be a permutation of (PP1,PP2,PP3) for node="...".

### example

```
// the process size: 4x3x4x2, the intra-node division is 1x1x2x2
#PJM --rsc-list "node=4x3x2"
mpirun ./rankmap_4d_lex 4 3 4 2 1 1 2 2
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

```
// the process size: 4x3x4x2, the intra-node division is 1x1x4x1
#PJM --rsc-list "node=4x3x2"
mpirun ./rankmap_4d_lex 4 3 4 2 1 1 4 1
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

### Todo

to allow 1ppn and 2ppn.


## Usage (for rankmap_4d)

Suppose that your application (=./a.out) uses 4-dim lexical process rankmap, 
of which process size is (P1,P2,P3,P4). Here, one of P* must be 4.

```
#PJM --rsc-list "node=PP1xPP2xPP3"
mpirun ./rankmap_4d_lex P1 P2 P3 P4 [dir]  # or mpirun ./rankmap_4d_lex_reversed ....
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

the parameter dir=[1234] specifies the intra-node direction.
If dir is not given, the first "4" in the process size becomes the inner node direction.  After removing the "4", the remaining 3-dim shape must a permutation
of (PP1,PP2,PP3).

### examples:

```
// the process size: 8x6x4x10, the intra-node direction is the 3rd
#PJM --rsc-list "node=8x6x10"
mpirun ./rankmap_4d_lex 8 6 4 10 3
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

```
// the process size: 4x4x4x6, the first direction is for the intra-node
#PJM --rsc-list "node=4x4x6"
mpirun ./rankmap_4d_lex 4 4 4 6
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```

```
// the process size: 4x4x4x6, the second direction is for the intra-node
#PJM --rsc-list "node=4x4x6"
mpirun ./rankmap_4d_lex 4 4 4 6 2
mpirun --vcoordfile ./rankmap_4d_list.txt ./a.out
```




## ACKNOWLEDGMENTS

I.K. acknowledges co-design working group for the lattice QCD
supported by MEXT's programs for the Development and Improvement for the Next
Generation Ultra High-Speed Computer System, under its Subsidies for Operating the
Specific Advanced Large Research Facilities, and Priority Issue 9 
(Elucidation of the Fundamental Laws and Evolution of the Universe) to be tackled by
using the Supercomputer Fugaku.
He also thanks the MEXT as ``Program for Promoting Researches on the Supercomputer Fugaku''
 (Simulation for basic science: from fundamental laws of particles to creation of nuclei),
JSPS KAKENHI (grant id: 20K03961) and JICFuS.

