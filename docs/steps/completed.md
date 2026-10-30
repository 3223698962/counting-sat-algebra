# Completed milestones

Historical checklist moved out of `ROADMAP.md` in package 0.22.1 so the roadmap lists open work only. These items are finite witnesses. They do not prove Theorem 1.1, refute ETH, or provide a general \(2^{o(n)}\) algorithm. `theorem_1_1_complete`, `lean_proof_complete`, and `recursive_contraction` stay false.

## Done

### Step 0 — counting polynomial and interpolation

- [x] Clause polynomial Φ_e and brute-force #F
- [x] Product identity on disjoint variable groups
- [x] Random replacements of degree ≤ R with value X_v at T = 0; Lagrange recovery of H(0)
- [x] Self-reduction demo that searches for one satisfying assignment
- Cost: each H(t), t ≠ 0, still enumerates 2^n assignments

### Step 1 — bounded-occurrence testbed

- [x] Generator and checker for 3-CNF with at most 3 occurrences of each variable
- [x] Degree-cap experiment at D = m·k·R and D = 3·n·R, compared with brute-force #F
- [x] Fixed-seed formulas under `benches/` and pytest coverage
- See `docs/PROGRESS.md` for the numbers from the fixed seed

### Step 2 — simplified (r, τ, D)-design and blocked sum

- [x] On a hand-built replacement family, group the variables and sum inside each group after the shared measurements are fixed
- [x] Check blocked sum = full-space sum H(t)
- [x] Enumeration of measurements may stay exponential; record how time grows with the number of groups
- Success criterion met: agreement on small n. The popcount tuple count is (r+1)^b, exponential in the number of groups. A subexponential runtime was not required and is not claimed
- The design here is a block partition plus τ explicit linear measurements, with D carried and unused. It is not the matrix M0

### Step 3 — small-parameter replacement matrix

- [x] One machine-built slice of M0 = V0 f_I(W0) and the rows P_u(T), at side lengths 4, 9, 16, and 64. The 4×4 exchange and the row ζ_(00,11) − (T/2) ζ_(10,01) − T ζ_(01,10) are also written out by hand
- [x] T = 0 restores the auxiliary coordinate and, on four explicit placements, restores X_v. H(0) matches brute force. The measured degree sits inside M·(R−1) and inside m·k·R
- [x] A legal design (partition, absolute row sum at most the side length, kernel condition) holds at every t ∈ {1, …, Q} on those placements, including one 3-clause that meets three groups
- [x] Full M0 is absent. The parameter list with R ≥ 32 and 6R ≤ Q forces side length 2^384 at r = 2, above the cap of 64. Schedules are hand-written. The point R = Q = 2 is an exact identity and sits outside the hypothesis R < Q
- Not claimed: Theorem 1.1, a subexponential runtime, or 2^{o(n)}

### Step 4 — Hermite / multi-coefficient interpolation

- [x] Exact rational Hermite interpolation from Hasse coefficients of H. The coefficient of index k is the k-th Taylor coefficient divided by k!. Both were checked on three closed forms
- [x] The jet at a sample t is the formal derivative of the clause product, summed over {0,1}^n in one pass. Derivative orders are not separate passes. The condition count stays D+1, packed into fewer distinct points than the Lagrange path: at most two coefficients per point, at most three, or every coefficient at t = 1
- [x] H(0) matches brute-force #F and `interpolate_H0` on the closed forms, on Step 0 random formulas, and on the Step 1 bench row n6-m4-r1 (7 paired points against 13, #F = 44). On these draws the measured degree equalled the cap
- Not claimed: fewer than D+1 conditions, a subexponential evaluation, Theorem 1.1, or 2^{o(n)}. A full jet of order D does quadratic arithmetic per assignment. Each sample point still costs 2^n

### Step 5 — sparsification and a reduction from general small k-CNF

- [x] Implication-cycle copying at B = 3. A variable of occurrence d > 3 becomes d copies, one per occurrence, plus (¬y_1 ∨ y_2) ∧ … ∧ (¬y_d ∨ y_1). Each copy occurs three times. A variable that already occurs at most three times stays one variable. #F' = #F, with no extra factor of 2
- [x] Brute force agrees on the hand formulas inside n' ≤ 20, including a free bit, an empty clause, one heavy width-3 variable, and a width-2 formula whose clauses sit in two blocks. The two-block row goes through the Step 2 blocked sum. The rows with n' ≤ 7 go through `interpolate_H0`, and on those draws the measured degree equalled M·R. Sixteen units are counted; their interpolation exceeds the assignment budget and is refused
- [x] Recorded, and not counted, when n' > 20. A seeded 3-CNF with n = 5 and m = 12 has n' = 36; the cap-3 bound n' ≥ 3m − 2n is already 26. Sixteen units fit in the brute window and miss the interpolation budget
- The output is mixed width. The exact-width-3 checker rejects the cycle clauses. Step 3 is not applied. B = 3 is the smallest cap this gadget meets
- Not claimed: Theorem 1.1, a subexponential runtime, 2^{o(n)}, or a reduction that keeps n' = n + o(n)

### Step 6 — parameter scan

- [x] For n in a machine-limited range, time brute force against the popcount blocked sum. The grid runs blocks of size 2 for b = 1..7 (n = 2..14), plus rows at r = 1, r = 3, and R = 2. On every row that runs, blocked H(0) matches brute force and the full-cube sum
- [x] Lagrange, paired Hermite, and the Step 1 degree check run on that family where the assignment budget allows, and are refused where it does not. At n = 14 the blocked sum still runs; Lagrange and paired Hermite do not
- [x] Occurrence copying is recorded on the popcount row (the map is the identity), on the seeded 3-CNF with n = 5 and m = 12 (n' = 36, cube not enumerated), on n' = 21, and on a clause that crosses blocks. Step 3 is not applied to that clause. Full M0 at R = 32 and Q = 192 is refused before allocation
- [x] The reproducible record is the operation counts: 2^n, (r+1)^b, the local-evaluation count, and the assignment budget. Wall times are one machine's clock. A refused pipeline has no timing
- Not claimed: Theorem 1.1, a subexponential runtime, 2^{o(n)}, or that the scan closes the distance to a full counting algorithm. That distance is the list in `docs/PROGRESS.md`

### Step 7 — larger slice and a balance-check scaffold

- [x] The slice cap is 256. `auxiliary_problems` passes on two points the cap of 64 refused: r = 3, Q = 2, R = 1, D = 81, and r = 2, Q = 4, R = 3, D = 256. Both have R < Q. T = 0 restoration, degree ≤ R−1, support, the measurement rows, and the kernel condition were checked, and on one explicit placement each, H(0) matches brute force inside M·(R−1) and m·k·R
- [x] `schedule_placement` checks 2(b+w)D² ≤ M and 2Qr² ≤ n, with n ≤ M ≤ bn. Those inequalities need n on the order of D². They fail for every formula this repository compares with brute force, and the call then refuses. It does not draw a placement
- [x] Caller-supplied schedules, including the four Step 3 placements and the two new points, are placed and `design_problems_at` is run. Choosing schedules when the inequalities hold is Step 8
- [x] `build_auxiliary_slice(2, 192, 32)` and `schedule_placement` at those parameters both refuse before a matrix is allocated. The next pure powers, 729 and 1024, stay above the cap. `verify --quick` runs the D = 81 point and the refusals. `verify --step7` also runs D = 256
- Not claimed: Theorem 1.1, a refutation of ETH, a subexponential runtime, or 2^{o(n)}. A passing balance check is not a schedule

### Step 8 — original balancer and legacy Boolean budget

Historical scaffold from 0.9.0. The explicit Boolean mode reproduces its budget stops; Step 9 implements the exact reference scheduler and default rational check.

- [x] Successive choice runs only when the Lemma 4.5 inequalities hold and the slice is under the cap. At each variable and each clause it tries every symbol in {0, …, r−1}^Q and keeps the symbol that minimizes a quadratic load potential. Ties break toward the lexicographically first symbol. Nothing is sampled. The schedule cube is not enumerated
- [x] In the explicit `design_method="boolean"` mode, on 64 disjoint unit clauses, r = 2, Q = 1, D = 4, the search meets label multiplicity 16 ≤ 2M/D = 32 and largest group 32, with 32·r ≤ 2n. `design_problems_at` is not called: an outside group has size 32, and that walk enumerates 2^32 Boolean points, above the budget 4096. `place_formula` is not called. The same kind of refusal, with 1 ≤ R < Q holding, is what 1024 disjoint units give at r = 2, Q = 2, R = 1, D = 16 (outside group 512, multiplicity 64 ≤ 128)
- [x] Formulas that fail the inequalities, including the brute-comparable ones, still refuse with no placement and without running the search. Full M0 and side length 1024 still refuse before allocation. A direct call on eight units at r = 2, Q = 2, which is not admitted, passes `design_problems_at` and the T = 0 restoration. `schedule_placement` does not make that call
- [x] `verify --step8` and `verify --quick` record these rows. Refusals are not timed. The direct check on eight units took about 0.01 seconds on one run. That clock is not a complexity result
- Not claimed: Theorem 1.1, a refutation of ETH, a subexponential runtime, or 2^{o(n)}. The search is not a proof of Lemma 4.5. The binomial tail is not checked. A budget refusal is not a placement

### Step 9 — reference conditional expectation and rational design

- [x] Minimize the exact reference potential Phi/(M(b+w)) + Xi/(Qn), including shared-incidence dependencies
- [x] Check each prefix expectation against all completions on tiny inputs; independently check the final potential
- [x] Eliminate each group's measurement matrix once and check the entire rational kernel through row-span membership
- [x] Complete a real admitted n=1024, Q=2, R=1 placement without a Boolean cube check
- [x] Check Lemma 4.7's exact binomial tail. The n=1024 point fails it; n=16384, Q=3, R=2, tau=1/2 meets every numerical hypothesis

### Step 10 — Section 3 common counting routine

- [x] Stream integer measurement tuples and scan each local Boolean cube twice, with no attained-measurement table
- [x] Recover truncated Taylor jets using U subset E and fixed-bit restrictions, then restore H(0) by Hermite interpolation
- [x] Compare all requested coefficients and counts on the four Step 3 examples with independent cube oracles
- [x] Check second-order recovery, negative literals, empty clauses, free variables, invalid designs and pre-enumeration budgets
- This counts supplied finite designs. Their measurement space can still be exponential; no general-formula complexity bound follows from these runs

### Step 11 — structured auxiliary operator

- [x] Query M0 and replacement-polynomial entries by Fourier/Krawtchouk sums over swap orbits
- [x] Compare every entry with dense slices of side 4, 9, 16, 64 and 81
- [x] Query coordinates at r=2, Q=192, R=32, D=2^384 with no dense allocation
- [x] Describe individual integer measurement modes and refuse oversized enumeration
- Step 12 connects this descriptor to finite sparse placements and the common counting entry point

### Step 12 — sparse admitted placement

- [x] Store occupied polynomial and measurement coordinates without M×n dense arrays
- [x] Bound replacement support, mode-coordinate queries, pullback candidates and sparse elimination before resource overuse
- [x] Check restoration, common denominator, coefficient norms and exact rational kernels
- [x] Complete n=16384, r=2, Q=3, R=2 placement with all Lemma 4.7 numerical hypotheses satisfied
- [x] Record 59392 polynomial entries and 14336 measurement entries versus 536870912 dense polynomial slots
- [x] Plan the common count and refuse before enumeration: q=5461, J=5462, zero Boolean vectors examined

### Step 13 — DIMACS and independent components

- [x] Parse strict DIMACS with declared free variables, multiline clauses and empty clauses
- [x] Normalize clauses and identify independent components; place each component afresh
- [x] Count each nontrivial component through the sparse shared-measurement backend; cap optional memoization
- [x] Add actual count CLI, JSON reporting, explicit budgets and exact hex output for huge counts
- [x] Count the eight-variable example as 196 and the 32769-variable disjoint 3-CNF as 7^10923
- [x] Verify no cube-oracle backend, free-bit factors, disabled memoization and refusals without partial counts
- This completes a finite supplied-CNF counter within component budgets. General large connected instances and the reference's asymptotic bounds remain open

### Step 14 — count the same global placement by replacement supports

- [x] Form components from every nonzero coefficient of the actual replaced factors, including dependencies introduced by replacement
- [x] Sum each small support cube at T=1 in a Taylor algebra covering the component's degree; recover each local H(0) and multiply with the free-variable factor
- [x] Count the admitted n=16384, Q=3, R=2 global placement as 1, with 10240 support components of at most four variables and degree at most four
- [x] Verify the same placement with and without memoization: 96 and 49152 actual local Boolean vectors, both giving count 1; the uncached run attains the planned bound
- [x] Check local nonconstant jets against an independent cube oracle, replacements coupling original CNF components, free bits and refusals before enumeration
- [x] Pass the Step 14 milestone of 181 regression tests and 17 reference checks, including export of the executed checks and failure status
- Large connected replacement supports still refuse within the configured budgets. This finite contraction does not implement the general Appendix D algorithm or establish its asymptotic bound

### Step 15 — finite Taylor tensor elimination on connected supports

- [x] Eliminate Boolean variables of the actual replaced factors over \(\mathbb Q[z]\), keeping a Taylor jet at T=1 through each component's degree bound
- [x] Plan every component before building a numerical table; a budget refusal returns no count and creates no table. This backend does not memoize
- [x] Count the 13-variable connected path as 610 after Step 14 refuses that same placement; bucket width 2 and 98 local Boolean vectors
- [x] Count a 1024-variable connected path plus two free bits as 4·Fibonacci(1026), with bucket width 2, 8186 local Boolean vectors and peak table size 4094; `count --backend tensor` matches that count
- [x] Check the coupled-replacement jet 7/4, 1/2, −1/4 against the independent full-cube Hasse jet, and refuse a width-13 factor before enumeration
- [x] Expose DIMACS `--backend common|tensor`, with tensor width and degree caps. The current suite is 205 pytest tests and 21 reference checks at version 0.16.0
- Step 15 is finite variable elimination on connected supports whose intermediate tables fit explicit budgets. It is not Appendix D's general modular contraction, and it does not establish that asymptotic bound. Independent-component counts and disjoint-support counts are not a general connected-instance bound. Theorem 1.1 is not implemented, ETH is not refuted, and there is no general \(2^{o(n)}\) algorithm. `theorem_1_1_complete` and `lean_proof_complete` stay false

### Step 16 — choose and run a finite member of the asymptotic parameter family

- [x] Compute Lemma 4.8's integer schedule from n and the factor scopes: λ = ⌈log₂ n⌉, ν = ⌈log₂ λ⌉, r the greatest integer with r²ν⁴ ≤ λ, Q = ⌊λ/(128ν)⌋, R = (β+3)⌈√λ⌉. β is the least integer ≥ 1 with b, w ≤ λ^β
- [x] Form τ, ε, D and the centered interval only after the proof's preliminary guard r ≥ 2 and Q > R. At 16384 alternating unit clauses the schedule is r=0, Q=0, R=16, so D is not formed. At λ=2¹⁹ the integers are r=2, Q=215, R=2900, and Q ≤ R, so D is not formed there either
- [x] Inside explicit choice budgets, keep a tuple only when `check_reference_parameters` accepts every Lemma 4.7 numerical hypothesis. The interval is the centered one from the proof. τ is the least value in (0, 1/2] that meets the exact binomial tail
- [x] On 8192 alternating unit clauses, no tuple is accepted. r=2, Q=2, R=1 has tail 7/16, above every τ ≤ 1/2, and r=2, Q=3 needs 2(b+w)D² = 16384 > M. Label `CHECKED`
- [x] On 16384 alternating unit clauses the chooser selects r=2, Q=3, R=2, τ=5/16. The earlier hand-picked point used τ=1/2. The tail is still 5/32, and the measurement allowance 2τn falls from 16384 to 10240. This tuple does not meet R ≥ 32 or 6R ≤ Q. Label `CHECKED`
- [x] Sparse placement at that tuple reproduces the Step 12 storage record: 59392 polynomial entries, 14336 measurement entries, 536870912 dense polynomial slots, multiplicity 256, largest group 8192, and group ranks (2048, 2048) at each of three nodes. The independent-measurement count is 4096, which is within 10240
- [x] Counting goes through the existing backends with explicit budgets, and a refusal returns no partial count. The common routine refuses before enumeration: q=5461, J=5462, a Boolean-vector bound of 36697 bits, and zero vectors scanned. Memoized support recovery counts 1, with 96 local Boolean vectors, 13 evaluations and 10227 memo hits. Finite Taylor elimination counts 1, with bucket width 4, 176128 local Boolean vectors, 262144 jet products and peak table size 72
- [x] A candidate-trial budget of 1 refuses before the schedule search. `build_auxiliary_slice(2, 192, 32)` still refuses before allocation; the dense side-length cap is unchanged
- [x] Running these finite tuples is not a proof of Lemma 4.8's asymptotic time or space bound. Step 16 does not implement Sections 5–6 or Appendix D. Step 15 remains finite Taylor elimination. Independent-component counts and disjoint-support counts are not a general connected-instance bound. Theorem 1.1 is not implemented, ETH is not refuted, and there is no general 2^{o(n)} algorithm. `theorem_1_1_complete` and `lean_proof_complete` stay false

### Step 17 — finite sparse systems and Lemma 5.1

- [x] Represent each equation as square-free monomials. Duplicate monomials cancel modulo 2. The empty monomial is 1, and the empty list is 0. Degree and the exact bound \((log_2 n)^c\) are checked on the supplied list
- [x] The integer lift is the sum of the monomial products. On every Boolean point of the checked systems, the lift is even exactly when the equation is zero over \(\mathbb F_2\)
- [x] Lemma 5.1 builds \(p_\ell(Y)=(1-Y)^\ell\sum_{j<\ell}\binom{\ell+j-1}{j}Y^j\) in \(\mathbb Z[Y]\) and divides by the monic \(h(Y)=Y^2-Y\). For \(\ell=2\) the coefficients are \(1-3Y^2+2Y^3=(1-Y)+(2Y-1)h(Y)\). Even integers are congruent to 1 and odd integers to 0 modulo \(2^\ell\)
- [x] The truncated product multiplies one equation at a time and keeps only terms whose total number of \(h\) factors is less than \(\ell\). Each retained term has at most \(m_{test}\) linear factors and at most \(\ell-1\) quadratic factors. Term count, coefficient bits, accumulator bits, and the Boolean cube are refused before the corresponding allocation. A refusal returns no count
- [x] On a 4-variable system the exhaustive residue modulo \(2^5\) is 4, matching the independent \(\mathbb F_2\) scan of all 16 points. An empty system is \(2^n\) and a constant equation \(1=0\) is 0, both without a scan
- [x] This is Sections 5.1–5.2 on a budgeted cube. Sections 5.3–5.5 replacement forms, Section 6 equation selection, and Appendix D are not implemented. The scan is not Theorem 2.2's running time. Theorem 1.1 is not implemented, ETH is not refuted, and there is no general \(2^{o(n)}\) algorithm. `theorem_1_1_complete` and `lean_proof_complete` stay false

### Step 18 — finite replacement sum in \(A_\ell\)

- [x] Build one factor \(B_i=\sum_{j<\ell}e_j(L_i)h(L_i)^j V^j\) per non-constant equation, and one neutral factor equal to 1 with scope \(\{v\}\) per declared variable. The input-degree bound is \(d=k(m_{test}+2(\ell-1))\) and \(d_H=kR(m_{test}+2(\ell-1))\)
- [x] Sum the product with the Section 3 measurement scan, now over \(A_\ell[z]/(z^J)\). Hermite interpolation is applied to each coefficient of \(V\). The least nonnegative residue of \(S_\ell(H(0))\) is the count when \(\ell\ge n+1\)
- [x] The rows are an explicit finite design: one full block, one empty block, and \(P_{f,v}(T)=(1+T/Q)e_v\). They meet Definition 3.1 on these scopes. They are not the Section 4 auxiliary array
- [x] Section 5.5 fixes variables that occur in more than \(b_0=\lceil(M/n)\lceil\sqrt{\lambda}\rceil\rceil\) scopes, substitutes the bits, and adds the branch counts. On the checked branch Lemma 4.8 has \(r=1\), \(Q=0\), \(R=8\), so \(\tau\) and \(D\) are not formed
- [x] On \(n=4\), \(x_0x_1+x_2=0\), \(x_3+1=0\), the residue modulo \(2^5\) is 4. The unreduced coefficient sum is 2308, with \(V\)-coefficients \((-2,6,0,-192,2496)\). The value at \(T=1\) is not 4. A split partition and a jet of length 2 recover the same \(n=2\) count 2. High-incidence branches return \((4,0)\)
- [x] A Boolean-vector budget, a branch budget, and a coefficient-bit budget each return no count. The scan or the branch loop has not started, and the jet table is not allocated. An invalid kernel returns no count before enumeration
- [x] This is a finite sum. It does not implement Section 6 equation selection, Appendix D, or Theorem 2.2's bit complexity. Theorem 1.1 is not implemented, ETH is not refuted, and there is no general \(2^{o(n)}\) algorithm. `theorem_1_1_complete` and `lean_proof_complete` stay false

### Step 19 — finite Section 6 equation selection and disjoint sparsification

- [x] Build clause violation polynomials over \(\mathbb F_2\): positive literal \(1+x\), negative literal \(x\), tautology 0, empty clause 1. Expand a binary row of \((e,c)\) pairs into square-free monomials. Duplicate clause indices cancel modulo 2
- [x] Compute Lemma 6.1's integers \(\mu=\lfloor\log_2 m\rfloor+1\), \(P=20000 B^2\mu^2\), \(D_0=320000 B^2\mu^3\). On \(n=m=3\), \(B=1\) the values are \(\mu=2\), \(P=80000\), \(D_0=2560000\). \(q_G=1\) and \(h\le 0\), so \(A(0)\) is the two even-parity checks. Every nonzero accepted vector has weight 2, which is at least \(m/P\); there are 2 rows, each of weight 3, and \(2\le 3n\). Label `CHECKED`
- [x] Refuse the local-code search for \(n=m=2048\), \(B=1\) before any matrix is allocated: the search is 524288 bits. Label `REFUSED`. No kernel vectors are scanned
- [x] Greedy-append rows of a budgeted list length, breaking ties by enumeration order. On two unit clauses in two variables, list length 1 uses 6 rounds and 24 counts. Final \(N(A)=1=\#F\), and `count_common_zeros` agrees. Label `COUNTED`. The paper-scale \(P=80000\), a candidate budget of 1, and a cube budget of 1 each return no count and enumerate nothing. Label `REFUSED`
- [x] Lemma 6.3 / Appendix A: normalize without unit propagation and without copying variables. Branch on the deterministic shared-literal rule. Children keep the declared \(n\), and satisfying assignments partition across the leaves. On the 3-variable 2-CNF the leaf counts are \((3,1)\), summing to \(\#F=4\), with occurrence at most \(B=6\). Label `COUNTED`. A node budget of 1 returns no leaves and no count. Label `REFUSED`. Step 5 implication-cycle copying increases \(n\); this tree does not
- [x] Proposition 6.4 integers: at \(n=20\), \(\lambda_n=5\), \(\nu_n=3\), \(s=0\), so the schedule is the exhaustive-search case. The integer witness \(s=2\) is only \(\lambda=521284\), \(\nu=19\). For \(k=2\) and \(\sigma=1\) the Appendix A power of two is \(M=64\). Label `CHECKED`. The \(O_k(\rho(n))\) bound is not claimed
- [x] This is a finite, budgeted Section 6 milestone. Appendix D is not implemented. Theorem 1.1 is not implemented, ETH is not refuted, and there is no general \(2^{o(n)}\) algorithm. `theorem_1_1_complete` and `lean_proof_complete` stay false

### Step 20 — Appendix D finite modular networks and local recovery

- [x] Implement \(B_{p,\ell}\) and its unramified extension \(A_p\), with exact integer polynomial reduction, deterministic irreducibility checks, unit inversion and residue-field interpolation. Choose \(e=4\lceil\log_2s\rceil\) within explicit search budgets
- [x] Preserve each implicit tensor's creation precision. Canonically lift entries separately, keeping their original descriptions; lifting is not a ring homomorphism
- [x] Encode CNFs as D.2 networks with trivalent equality trees and one subdivision per wire. Preserve ordered ports, parallel edges and declared free variables. The 8-variable example has 12 edges, 14 vertices and count **196**, with 4096 edge assignments and 57344 entry evaluations. Four copies of a unit clause on three declared variables give count **4** through a 10-edge equality-tree network
- [x] Compute integral Cauchy local coordinates by canceling known powers of two before inverting odd denominators. On \(a=1,t=1,x=0,y=1\), the coefficients are \((-3,-22,-48,-32)\). Check \(W(0)=I\)
- [x] Evaluate D.8 on bounded port-pair cubes, interpolate over the residue field, and apply D.9 with D.5's tapered precisions. A two-edge network gives \(9+4V\) modulo 32 with \(V^2=0\). The one-level check uses 224 port pairs; the two-level check uses 1440 pairs and precisions 101/149. A degree-6 nonlinear deformation gives \(8+5V\), with an independent rational polynomial oracle in pytest
- [x] Refuse before network allocation or before kernel/jet allocation and entry evaluation when the relevant budget fails. Export `verify --step20`, include its 14 checks in `--reference` and `--quick`, and keep the theorem/proof flags false
- These are finite D.1–D.2 and D.4.1, D.4.3–D.4.4 objects. The direct leaf routine enumerates \(2^M\) edge assignments; D.8 enumerates \(2^{2M}\) port pairs per value. The supplied packet layout does not establish Lemma D.3. No general recursive contraction or asymptotic bound is claimed. See [APPENDIX_D.md](../APPENDIX_D.md)

### Step 21 — collision labels, packet layouts and child topology

- [x] Compute D.4's integer schedule without floating point, including H, Delta, u and cumulative a. Use a popcount identity for logarithmic witnesses without constructing s or iterating through huge depths. On s=16, M=20, u=a=0; the active guard stops before forming L or D_arr
- [x] Check D.6 exactly: M >= 2*D_arr^2*Delta, M >= L*D_arr and 2*L <= Delta. The first condition contains D_arr squared. Keep numerical admission separate from actual finite conclusions and from membership in the fixed D.4 schedule
- [x] Select labels by exact conditional expectations using incidence histograms. Include shared-endpoint and parallel-edge correlations; preserve vertex order and break ties by the lowest label. Independently enumerate remaining labels to check every candidate expectation on small multigraphs
- [x] Pad every pair-label class uniformly and construct packets with a fixed orientation convention. Equal-coordinate slot bits select distinct copies. Check actual loads, potentials, group ports, packet counts and degree bounds; reject inconsistent certificates
- [x] Construct D.4.6 port topology with original port order and isolated vertices preserved. This produces no packet tensor and evaluates no entries. Child edge and degree bounds are checked on the finite topology; recursive invariants remain open
- [x] On a 128-edge matching, r=2, a=1, Delta=4 passes D.6 and yields 64 packets; its 128-edge child does not shrink to M/8. On an 8192-edge matching, r=64, a=1, Delta=4 fails D.6 and is refused by default. Explicit finite mode yields 4096 packets and a largest child topology with 256 edges and degree 4. These are not D.4 schedule members
- [x] Feed an automatic layout into Step 20 recovery: two repeated positive units on two declared variables yield count 2, including the free variable. This scans 6656 port pairs and evaluates 33280 entries; it does not run recursive contraction
- [x] Preflight label trials, score updates, optional trace and conservative padded storage before label search. Check child budgets before allocating its topology. Export `verify --step21`; add 33 pytest tests and 14 reference rows. Full v0.22.0 verification passes 322 pytest tests and 84 reference checks, with `--quick` and both CI validators
- See [PACKET_PACKING.md](../PACKET_PACKING.md) for APIs, exact formulas, finite witnesses and refusal semantics

### Step 22 — finite D.4.5 weighted modes

- [x] Build an explicit mode table for one packet and one active node. Weight is \(b+1\). Source and target are the two symbols of that coordinate. Real columns are recorded and not evaluated
- [x] Refuse before allocation when the mode count, the \(6M/D\) weight class, the exponent-tuple product, or the partial-choice bound does not fit. The refusal has no table, no tuple list, and no count
- [x] Keep only exponent tuples whose weighted degree is strictly below the precision. An inadmissible degree record, \(\deg X > \deg U + \deg Z\), is refused before it is stored. That inequality is a checklist shape, not the SqrtLog incidence polynomial
- [x] Recover one crossing Cauchy entry by the truncated geometric series and compare it with Step 20's local polynomial. On \(a=2,t=1,x=0,y=2,\zeta=1\), precision 4, the residue is **8**. Coefficient valuations are \((3,4,5)\)
- [x] After the mode table is admitted, recover the two-edge network through the existing Step 20 routine: **9+4V**, 224 port pairs, 448 entry evaluations. Mode-budget refusal does not call that routine. Packet-tensor evaluations stay 0
- [x] Export `verify --step22`. The frozen description does not carry cache, rows, or origins. `recursive_contraction` stays false
- See [WEIGHTED_MODES.md](../WEIGHTED_MODES.md) and [LEAN_SPEC_CHECKLIST.md](../LEAN_SPEC_CHECKLIST.md)

### Step 23 — Lean-compatible syntax and streamed disjoint counting

- [x] Encode and decode Lean's framed binary CNF syntax, preserve declared n and implement canonical little-endian natural count words of length at most n+1
- [x] Choose the first dyadic M that passes the exact integer parameter test. For k=3, sigma=1, select M=512 at exponent 9 and check that all smaller choices fail
- [x] Normalize clause sets and implement the priority-core disjoint split. Independently check each leaf's models, their disjointness and their union on finite inputs
- [x] Preflight and replay a depth-first frontier without storing completed leaves or histories. Count the four-variable star as 9 and three independent stars as 125 while retaining the original variable set
- [x] Connect common/tensor finite leaves to `--backend lean`, report the completed work on refusal and return no partial total. Export 14 `verify --step23` rows

### Step 24 — native ternary indicator and exact query algebra

- [x] Construct Newton indicator coefficients modulo `3^t`, normalize by `H^t` and divide into t three-coefficient digits; encode each coefficient in 2t bits
- [x] Implement independent U/Z truncated dense arrays and complete-product marker readout; verify small-precision coefficients and convolution with independent oracles
- [x] Implement immutable squarefree support tables with support-union multiplication, exact Boolean monomial weights and live-table/scratch logical storage charges
- [x] Count accepted and marked clause-violation queries without Boolean assignment scans; preserve original clause IDs, repetitions, declared n and fixed-coordinate semantics. Verify accepted/marked counts 8/4 and nonzero higher marker sums
- [x] Stream unit history rows for CNF leaves and connect `--leaf-backend native` to the disjoint frontier. Verify example count 196, star count 9 and three-star count 125
- [x] Share monomial/scalar operation budgets across every leaf; refuse without partial total after an earlier leaf succeeds. Export 17 `verify --step24` rows
- [x] Audit the supplied Lean endpoints' static import closure and record source hashes. No missing local imports, code-level placeholders or local axiom declarations are found among 3969 modules; kernel compilation remains pending
- Full sparse-code search, greedy selection, counting extensions, placement, streamed Taylor leaves, whole-program resource refinement and Appendix D recursion remain open. See [LEAN_NATIVE.md](../LEAN_NATIVE.md) and [ROADMAP.md](../../ROADMAP.md)

### Step 25 — scheduled sparse codes and violation-mass greedy selection

- [x] Construct the scheduled ternary Tanner code from the exact SparseCodeParameters integers, including floor divisions and the least cube dimension; stream both promised searches and refuse without a partial code when a matrix, candidate, word or storage budget fails
- [x] Regenerate every row from slots (column = source vertex times degree plus port, coefficient from the local matrix) and independently check the relative-distance identity, the kernel lower bound distanceDenom\*min_weight >= m, the depth bound and hand-built accept/reject certificates
- [x] Enumerate the canonical ternary candidate family once per increasing odometer rank without storing a library; check coverage, ranks and per-tuple weights against an independent enumeration
- [x] Score every candidate with real marked queries accumulating every clause's violation count; a strict minimum keeps the earliest candidate, and the default round cap is loopCap = 2\*(n + size m + 1)
- [x] Reproduce the star's full greedy trajectory by brute force over assignments; zero mass accepts exactly the models, and the final accepted query runs separately from every trial
- [x] Connect `--backend lean --leaf-backend hash`; share marked-query budgets across the frontier and refuse without a partial total after an earlier leaf completed
- [x] Export 15 `verify --step25` rows and 18 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The Lean schedule's candidate depth `distanceDenom` is above the Python candidate budgets, so the depth is an explicit finite parameter. The high-variable loop, counting extensions, placement, streamed Taylor leaves and whole-machine resource refinement remain open

### Step 26 — high-variable frontend and the streamed rank loop

- [x] Compute EarlyThreshold's integer threshold from Asymptotic.L = floor(log2 (n+16)) + 1 and the exact ceiling cutoff (totalScope\*threshold)/ambient on the original declared n
- [x] Build the ambient-order high mask and independently recheck the membership rule, both integer bounds (cardinality\*(cutoff+1) <= totalScope, cardinality\*threshold <= ambient) and the equal-degree boundary exclusion
- [x] Stream all 2^|H| high assignments through the binary odometer; every rank runs one real conditioned query with the high bits fixed and every total coordinate retained
- [x] Verify the rank-sum identities against direct queries for unit rows, empty and explicit histories, marked clauses and base-fixed coordinates, including a two-high-coordinate case with four ranks
- [x] Route the hash leaf's marked and accepted queries through the loop with `--high-loop`; direct and routed greedy counts agree on the hub, star, forest and streamed examples
- [x] Refuse without a partial count on rank, query and algebra budgets
- [x] Export 12 `verify --step26` rows and 21 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The full parameter frontend, degree and occurrence guards, and the counting extension remain open obligations

### Step 27 — counting extension rings and Hermite placement

- [x] Reproduce LeafParametersProgram.produce exactly: L, w = floor(log2 L) + 1, s = isqrt(L)/w^2, r, layout denominator, Q, energy, radius R, t, d, N and J, all by integer arithmetic
- [x] Check the guard boundary without enumerating inputs: Q = 0 below the threshold, the K = 2^24 radius case refuses (4R > Q), K = 2^25 passes with card\*threshold <= ambient and the field cap 3^d <= 16(Q+1)^2
- [x] Search the irreducible monic lift by complete odometer enumeration with a division test, and verify the first witness against an independent multiplication oracle on every rejected prefix
- [x] Build the extension ring O = (Z/3^t)[X]/(F): commutativity, distributivity, reduced representatives, residue units, Newton-lifted inverses and non-unit refusals
- [x] Verify unique natural readback below 3^t: natCast injectivity, decodeCount round trips and the 3^t-periodic zero test
- [x] Fix the radix nodes 1..Q with unit pairwise differences and physical 2t-bit coefficient words; reject 0 and 3^d
- [x] Recover degree-below-Q\*J targets exactly from exact Taylor jets by Newton-Hermite with divisions only by node differences; detect mismatches at the condition count
- [x] Check the Hermite margin Q\*J > width\*R\*(2\*factors + 3(t-1) + 1) with the exact division-plus-one identity, and refuse on degree, candidate, division and ring budgets
- [x] Export 12 `verify --step27` rows and 26 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The full native placement into the ABI, grouping replacement and the streamed Taylor leaf remain open obligations

### Step 28 — streamed correction pairs and Taylor readback

- [x] Stream all correction pairs (|E| < J, U subset of E) in the exact fold order with strict-budget pruning; check the pair counts, the codec round trips and the enumeration against an independent mask family
- [x] Extract every squarefree coefficient by the alternating difference sum and match the support table on all masks
- [x] Verify the Boolean finite-difference identity: the streamed totals equal the direct weight sums and the brute-force point sums on signed probes, and the model indicators of the star, double star and hub give 9, 25 and 129
- [x] Emit per-node blocks of the J exact Taylor coefficients of the shifted target as Q\*J\*d little-endian 2t-bit words, and rebuild them from the stream alone
- [x] Recover the count by Newton-Hermite with Q\*J above the degree; detect short and trailing streams, reordered nodes and out-of-range readback
- [x] Refuse on order, pair, assignment, term, jet-order and stream-word budgets without a partial sum, including the first-order constant case
- [x] Export 9 `verify --step28` rows and 29 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The production native leaf still evaluates the full squarefree support table. Switching it to the streamed path with immediate U-specialization and complete pair regeneration, the native placement and the whole-machine resource refinement remain open

### Step 29 — streamed pointwise correction-pair leaf

- [x] Evaluate the native queries through the pointwise pair identity: one evaluation per correction pair, no squarefree support table, `polynomial_terms` always zero
- [x] Bound the order by all row degrees plus the marked clause degree; a fully forced-violated clause contributes a constant zero; the precision rule `3^t > 2^n` mirrors the native path
- [x] Match the native algebra exactly: accepted 8 and marked 4 with marker sums (8, 36, 45, 0) on the fixture; unit-row leaf counts 9/25/129/196 on the star, double star, hub and supplied example; randomised seeded cross-checks agree on counts, marker sums and precision
- [x] Count a star instance whose native run refuses on the polynomial-term budget, evidencing the support-table replacement
- [x] Connect `--leaf-backend stream` with shared pair and point-evaluation budgets; refusals publish no partial total
- [x] Export 9 `verify --step29` rows and 26 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The full machine's degree reduction, the native placement and the whole-machine resource refinement remain open; the default leaf stays `common`

### Step 30 — leaf assembly and placement

- [x] Place the ring constants at the recorded QueryReady ports: modulus bits, the `d^3` structure words, capacity width, scalar codes of 0/1/-1/1/2, unary degree/precision/jet/`3^d-2` exponent words, the `6t^2`-bit idempotent code, flat and framed node tables and the unary row/ring widths
- [x] Check every placed length against the recorded formula and every value against an independent oracle: power words, half-inverts-two, node digits, inverse and negative codes, and pairwise unit differences with a zero diagonal
- [x] Frame the regenerated rows exactly like `Encoding.encodeList`, decode back to the raw payloads and to the slot lists, and satisfy the `initial_encoding` identity; the append encoding matches `true :: frame (payload.reverse) ++ history`
- [x] Thread one run from the scheduled code through placement and the framed stream to the greedy rounds and the little-endian count word, with the decoded-history marked score equal to the direct rows
- [x] Count the star, hub, forest and supplied example as 9/129/25/196, equal to brute force and the direct greedy; refusals on node, degree, exponent, table and stream budgets publish no count
- [x] Export 9 `verify --step30` rows and 27 pytest checks; keep the theorem, Lean and recursive-contraction flags false
- The physical 896-stack machine, the query-derived ring target, grouping replacement and the prefix-space refinement remain open
