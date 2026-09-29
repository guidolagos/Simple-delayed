# Test systems: Yeh et al. (2024)

Source: W.-C. Yeh, G. Yang, C.-L. Huang, "A New Node-Based Algorithm for
Identifying the Complete Minimal Cut Set", *Electronics* 13 (2024) 603
(`\cite{yeh2024new}`), Figure 7 (test networks 1–22) and Table 9.

All systems are two-terminal networks: the components are the edges, and the
system works iff the source (`Start`) and the target (`End`) are connected by
working edges.

## Files

- `yeh2024_raw.txt`: `SistemasYeh_corrected.txt`, the edge lists used for the
  subsystems paper (22 networks). Kept unchanged, for provenance.
- `yeh2024_systems.csv`, `yeh2024_edges.csv`: the systems used here. Component
  `e` of a system is its `e`-th edge in the raw file.

## Differences with the raw file

1. **Networks 21 and 22 are dropped.** Section 6.3 of Yeh et al. states that
   Figure 7(21) and 7(22) are Figure 7(13) and 7(18) redrawn with renumbered
   nodes, and their Table 9 lists 20 networks (13 and 18 appear only through
   21 and 22). In the raw file, network 21 is isomorphic to network 13
   (terminals included), and network 22 is network 18 plus one edge (the drawing
   of Figure 7(22) has 23 edges, but Table 9 reports 22 arcs for it). The 20
   systems kept are networks 1–20.
2. **Networks 6 and 7 are undirected.** The raw file labels them `6D` and `7D`
   and lists their edges as arcs. In Yeh et al.'s Figure 7 they are mixed
   graphs (some edges have arrows), but their algorithm and Table 9 treat all
   networks as undirected: network 7 has 133 minimal cuts as undirected (Table 9:
   133) and 110 with the arcs of the raw file. With the arcs of the raw file,
   network 6 also had four irrelevant components (the edges at node 4 of the
   raw file). Read as undirected, every component of every system is relevant.

## Check against Table 9 of Yeh et al. (2024)

`n` = nodes, `m` = edges, `c` = number of minimal cuts (edge sets), computed by
exhaustive enumeration in C++ and by the Julia tests (systems with m <= 16).

| system | ours n, m, c | Table 9 n, m, c | |
|---|---|---|---|
| 1 | 4, 5, 4 | 4, 5, 4 | ok |
| 2 | 6, 8, 9 | 6, 8, 9 | ok |
| 3 | 5, 8, 8 | 5, 8, 8 | ok |
| 4 | 6, 9, 9 | 6, 9, 9 | ok |
| 5 | 9, 12, 28 | 9, 12, 28 | ok |
| 6 | 7, 14, 26 | 7, 14, 25 | c differs by one (see below) |
| 7 | 11, 21, 133 | 11, 21, 133 | ok (undirected) |
| 8 | 9, 13, 24 | 9, 13, 24 | ok |
| 9 | 8, 12, 19 | 8, 12, 19 | ok |
| 10 | 8, 12, 20 | 9, 14, 25 | Table 9 row repeats row 15; Figure 7(10) has 8 nodes, 12 edges |
| 11 | 7, 12, 20 | 7, 12, 20 | ok |
| 12 | 8, 13, 29 | 8, 13, 29 | ok |
| 13 | 16, 30, 644 | 16, 30, 644 (row 21) | ok |
| 14 | 21, 26, 528 | 21, 26, 528 | ok |
| 15 | 9, 14, 25 | 9, 14, 25 | ok |
| 16 | 10, 21, 104 | 10, 21, 104 | ok |
| 17 | 18, 27, 1249 | 18, 27, 1249 | ok |
| 18 | 13, 22, 219 | 13, 22, 214 (row 22) | c differs (see below) |
| 19 | 20, 30, 7376 | 20, 30, 7376 | ok |
| 20 | 16, 24, 105 | 16, 24, 105 | ok |

Unresolved: for networks 6 and 18 our edge lists coincide with Figure 7(6) and
7(18) of the paper (checked by graph isomorphism with a transcription of the
figure), yet the minimal-cut counts differ from Table 9 by 1 and 5. For network
6, no reading of the arrows gives 25 (all directed: 13; arrows as drawn,
other edges undirected: 18 or 8; undirected: 26). We attribute these to the
table, which has at least one clear error (row 10).
