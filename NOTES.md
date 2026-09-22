# Erdős problem #1162 — working notes

> **Give an asymptotic formula for the number of subgroups of $S_n$. Is there a
> statistical theorem on their order?**
> — a problem of Erdős and Turán, source [Va99, 5.73]

**Status: open. Nothing here solves it.** The tracker records the problem as
`OPEN` with the note *"This is open, and cannot be resolved with a finite
computation."* That is correct, and no computation below is offered as a proof
of an asymptotic statement.

What this directory does contain:

1. a correction to an ambiguity in the problem statement as displayed
   (§1: the logarithm is base 2);
2. exact counts, far past the published ranges, of the family of subgroups
   that carries the whole first-order answer (§2–§3);
3. two limiting constants identified in closed form and confirmed to 9–12
   decimal places (§4, §5);
4. a second-order term derived from the orbit structure and confirmed by
   out-of-sample extrapolation (§6);
5. an exactly located lossy step in the one published theorem that addresses
   the *second* half of Erdős's question, and a conjecture for what replaces
   it (§7);
6. the evidence that does **not** support the conjecture, stated plainly (§8).

Every number quoted below is re-derived from the artifacts by
`src/final_check.py`, which exits non-zero if any of them is wrong.

---

## 1. The statement: which logarithm?

The problem page states Pyber's result as `log f(n) ≍ n²` and
Roney-Dougal–Tracey's as `log f(n) = (1/16 + o(1))n²` without fixing a base,
and the constant 1/16 is meaningless without one.

It is base 2. From page 1 of [RoTr25], read both from `pdftotext -raw` and
from a 150 dpi render of the page:

> "|Sub(S_n)| of subgroups of S_n is at most 2^{ξn²+o(n²)} where ξ = ⅙ log 24
> (our logarithms are to the base 2). Pyber conjectured, however, that
> |Sub(S_n)| = 2^{n²/16+o(n²)}."

Independently, the OEIS A005432 comment records Pyber's brackets as
`c^{n²(1+o(1))} ≤ a(n) ≤ d^{n²(1+o(1))}` with `c = 2^{1/16}`, `d = 24^{1/6}`,
which is the same statement with the base made explicit.

*(The default and `-layout` modes of `pdftotext` both sort text by position
and splice superscripts from the displayed math into this sentence — it comes
out as "are to the 2 2 base 2)" — so a verbatim search against them fails even
though the sentence is on the page. `src/final_check.py` uses `-raw`.)*

**All logarithms in these notes are base 2.**

### What is and is not answered

| Half of #1162 | Status |
| --- | --- |
| Asymptotic formula for $f(n) = \|\mathrm{Sub}(S_n)\|$ | Largely answered. [RoTr25] Theorem 1: $2^{n^2/16+\alpha n\log n} \le f(n) \le 2^{n^2/16+\beta n^{3/2}}$ for absolute constants $\alpha>0,\beta$. The second-order term is not pinned down: $n\log n$ vs $n^{3/2}$. |
| Statistical theorem on their order | Barely touched. [RoTr25] Theorem 6 is the only result, it is one-sided, and §7 shows its constants are not optimal. |

## 2. What was computed exactly

### 2.1 The full subgroup lattice, $n \le 8$

`src/subgroups.c` enumerates every subgroup of $S_n$ by expanding one
representative per conjugacy class and closing under conjugation.

| $n$ | subgroups $f(n)$ | classes | agrees with |
| --- | --- | --- | --- |
| 1–7 | 1, 2, 6, 30, 156, 1455, 11300 | 1, 2, 4, 11, 19, 56, 96 | A005432 / A000638 |
| 8 | 151221 | 296 | A005432(8) / A000638(8) |

$n=8$ took 883 s. This is ground truth for the *actual* object of #1162 — the
full multiset of subgroup orders — and it is the anchor that the much larger
elementary abelian computation is checked against.

**A wrong version of this program is preserved in the history of this
directory and is worth recording.** The first attempt pruned with a "covered"
set, marking `⟨H,g⟩` covered for each new generator `g`. It produced
1, 2, 6, 30, **144, 1066, 6415** — wrong from $n=5$ on, 43 % low at $n=7$ —
because if `g ∈ ⟨H,g'⟩` for an earlier `g'`, then `⟨H,g⟩` may still be a
different and *smaller* subgroup, so genuine minimal overgroups were dropped.
The rewrite prunes by right cosets instead (`⟨H,g⟩ = ⟨H,hg⟩` for `h ∈ H`, so
the whole coset `Hg` may be skipped), which is sound. The error was invisible
by inspection and was caught only by comparing against OEIS.

### 2.2 Elementary abelian 2-subgroups, $n \le 512$

Write $a(n,k) = \#\{H \le S_n : H \cong (\mathbb{Z}/2)^k\}$ and
$e(n) = \sum_k a(n,k)$. `src/elemab.py` computes these exactly:

* $\#\mathrm{Hom}(G,S_n) = n!\,[x^n]\exp\!\big(\sum_{H\le G} x^{[G:H]}/[G:H]\big)$;
  for $G=\mathbb{F}_2^k$ the index-$2^j$ subgroups are the codimension-$j$
  subspaces, so $a_{2^j} = \binom{k}{j}_2 (2^j-1)!$ and only $j \le \log_2 n$
  contribute — each coefficient costs $O(\log n)$ big-integer multiplications;
* $\#\mathrm{Inj}$ from $\#\mathrm{Hom}$ by $q$-binomial inversion at $q=2$,
  $i_k = \sum_m (-1)^{k-m} 2^{\binom{k-m}{2}} \binom{k}{m}_2 h_m$;
* $a(n,k) = i_k / |\mathrm{GL}_k(2)|$.

For example `rank_counts(8) = [1, 763, 6685, 4440, 350]`, and $e(512)$ is an
exact integer of 18 997 bits.

**Why this family.** It is not a proxy chosen for convenience: it already
carries the entire first-order answer. $\log e(n) = n^2/16 + \Theta(n\log n)$
(§6), and $\log f(n) = n^2/16 + o(n^2)$, so elementary abelian 2-subgroups
account for all of the leading term. §8 records what that does *not* imply.

### 2.3 Verification chain for $a(n,k)$

Four independent paths, all agreeing:

| check | range | status |
| --- | --- | --- |
| `src/verify_elemab.py`, brute-force permutation enumeration sharing no code with `elemab.py` | $n \le 10$ | all ranks match ($n=10$: `[1, 9495, 330015, 620325, 158220, 5670]`, 397 s) |
| `src/subgroups.c` lattice, a third implementation | $n \le 8$ | all ranks match |
| $a(n,1) = A000085(n) - 1$ (involutions, independent recurrence) | $n \le 30$ | matches |
| structural assertion $a(n,\lfloor n/2\rfloor+1) = 0$ | $n \le 512$ | holds — a sharp test of a massively cancelling alternating sum |

## 3. The rank distribution: an exact statistical theorem for this family

For each $n \le 512$ the distribution $k \mapsto a(n,k)/e(n)$ is known exactly,
so the following are *computations*, not estimates:

| $n$ | mode | mean$/n$ | sd | window holding $\ge 1-10^{-6}$ |
| --- | --- | --- | --- | --- |
| 64 | 16 | 0.249997 | 0.849353 | $\pm 4$ |
| 128 | 32 | 0.250000 | 0.849306 | $\pm 4$ |
| 256 | 64 | 0.250000 | 0.849306 | $\pm 4$ |
| 512 | 128 | 0.250000 | 0.849306 | $\pm 4$ |

The mode is exactly $n/4$ for every $n \equiv 0 \pmod 4$ tested from 32 up.
**The width does not grow with $n$**: the standard deviation converges to an
absolute constant and the $1-10^{-6}$ window stays at $\pm 4$ while $n$ grows
by a factor of 8.

The limit law is identified in closed form. For $n \equiv 0 \pmod 4$,

$$\frac{a(n,\ n/4+j)}{a(n,\ n/4)} \longrightarrow 2^{-j^2},$$

and this holds to better than $10^{-6}$ already at $n=128$:

| $j$ | $-3$ | $-2$ | $-1$ | $0$ | $+1$ | $+2$ | $+3$ |
| --- | --- | --- | --- | --- | --- | --- | --- |
| $n=512$ | 0.001953 | 0.062500 | 0.500000 | 1 | 0.500000 | 0.062500 | 0.001953 |
| $2^{-j^2}$ | 0.001953 | 0.062500 | 0.500000 | 1 | 0.500000 | 0.062500 | 0.001953 |

The standard deviation of that law is

$$\sqrt{\frac{\sum_{j\in\mathbb{Z}} j^2 2^{-j^2}}{\sum_{j\in\mathbb{Z}} 2^{-j^2}}} = 0.8493059610\ldots$$

against a measured 0.849305961022908 — and the agreement is far better than
the stored artifact suggests. `data/elemab_stats.json` keeps `sd` as a double,
which truncates at about $10^{-16}$ relative; recomputing the variance from
the exact integer rank counts at 60 decimal places gives

| $n$ | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- |
| $\|\mathrm{sd}(n) - \text{closed form}\|$ | 4.7e−5 | 8.0e−9 | 1.1e−17 | 2.7e−36 |
| decimal places | 4 | 7 | 16 | **35** |

At $n=512$ the exact standard deviation agrees with
$0.849305961022908280743697585694$ in all 30 digits shown. The error is
roughly squaring at each doubling of $n$ — the $\log_{10}$ errors are
$-4.3, -8.1, -17.0, -35.6$, ratios $1.87, 2.10, 2.10$ — i.e. it decays like
$2^{-cn}$ with $c\approx0.23$, not polynomially.

Two cautions. The deviations read off the stored doubles instead
(2.8e−10 at $n=256$, 1.2e−9 at $n=512$, and non-monotone) are **float noise
in the artifact, not convergence error**; I reported them as convergence
before checking, and they are wrong. And this is convergence of an exact
finite quantity to a closed form, checked at four values of $n$ — it is
extremely strong numerical evidence, but it is not a proof.

**Where it comes from.** $a(n, n/4+j)$ behaves like $2^{nk/2-k^2}$ near
$k=n/4$ (§6), and putting $k = n/4+j$ gives $n^2/16 - j^2$: the quadratic in
$k$ has its vertex at $n/4$, and the $2^{-j^2}$ profile is exactly the same
mechanism that makes $\binom{r}{r/2+j}_2 \approx 2^{r^2/4-j^2}$.

So, **for this family**, the answer to the second half of #1162 is sharp:
$\log_2 |H| = n/4 + O(1)$, with an explicit limiting distribution.

## 4. A single maximal elementary abelian is not enough

The classical Pyber lower bound takes one maximal elementary abelian
$E \cong \mathbb{F}_2^{\lfloor n/2\rfloor}$ and counts its subspaces, giving
the Galois number $G_{\lfloor n/2\rfloor}$. Exactly:

$$\log_2 G_r - \frac{r^2}{4} \longrightarrow \log_2\!\frac{\vartheta}{\eta} = 2.8820499654\ldots,
\qquad \vartheta=\sum_{j\in\mathbb{Z}}2^{-j^2},\quad \eta=\prod_{i\ge1}(1-2^{-i}).$$

Measured: $2.882049965$ at $r = 128, 256, 512$ (even) and $2.882046186$
(odd), stable to 9 decimals; the deviation from the closed form falls like
$2^{-r/2}$ ($5.2\times10^{-5}$ at $r=32$, $8.0\times10^{-10}$ at $r=64$).

$\eta^{-1} = 3.46275\ldots$ is exactly the constant $\zeta_2$ of [RoTr25]
Lemma 2.4 ("$P$ has at most $\zeta_p p^{k(\ell-k)}$ subgroups of order $p^k$
… $\zeta_p < 4$").

**Consequence.** One maximal $E$ contributes $n^2/16 + O(1)$ — a *bounded*
excess. So the entire $\alpha n\log n$ of [RoTr25] Theorem 1's lower bound
comes from the multiplicity of elementary abelian subgroups, not from any
single one.

## 5. An explicit lower bound for the constant $\alpha$

[RoTr25] Theorem 1 gives $f(n) \ge 2^{n^2/16+\alpha n\log n}$ with $\alpha>0$
unspecified. Since $e(n) \le f(n)$ and $e(n)$ is exact, every computed $n$
certifies a value:

| $n$ | 16 | 32 | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- | --- | --- |
| $(\log_2 e(n)-n^2/16)/(n\log_2 n)$ | 0.4487 | 0.4790 | 0.5039 | 0.5276 | 0.5487 | **0.5671** |

This is a statement about individual $n$, not a proof that $\alpha \ge 0.567$
for all $n$; the sequence is still rising.

## 6. The second-order term is $\tfrac34 n\log_2 n$

**Derivation (orbit structure).** Let $\mathbb{F}_2^k$ act on $[n]$ with every
orbit of size $m=2^j$. Exactly
$N_j(n,k) = \frac{n!}{t!}\big(\binom{k}{j}_2/m\big)^t$ with $t=n/m$, and
dividing by $|\mathrm{GL}_k(2)|$,

$$\log_2 \frac{N_j(n,k)}{|\mathrm{GL}_k(2)|} = (1-2^{-j})\,n\log_2 n \;+\; \frac{jkn}{2^j} - k^2 \;+\; O_j(n).$$

Maximising the $n^2$ part gives $k = jn/2^{j+1}$ with value $j^2n^2/2^{2j+2}$:

| $j$ | orbit size | optimal $k/n$ | $n^2$ coefficient | $n\log_2 n$ coefficient |
| --- | --- | --- | --- | --- |
| 1 | 2 | 1/4 | **1/16** | 1/2 |
| 2 | 4 | 1/4 | **1/16** | **3/4** |
| 3 | 8 | 3/16 | 9/256 | 7/8 |
| 4 | 16 | 1/8 | 1/64 | 15/16 |

Sizes 2 and 4 **tie** at the leading constant $1/16$ — the constant of Pyber
and of Theorem 1 — and size 4 wins the tie-break at order $n\log n$. Mixtures
do not help: with $c_j$ the fraction of points in orbits of size $2^j$, the
$n^2$ coefficient is $\kappa\sum_j c_j j2^{-j} - \kappa^2$ and
$\sum_j c_j j 2^{-j} \le 1/2$ with equality iff $c_j$ is supported on
$\{1,2\}$; on that face the $n\log n$ coefficient is $1/2 + c_2/4$, maximised
at $c_2=1$. Hence

$$\log_2 e(n) = \frac{n^2}{16} + \frac34 n\log_2 n + O(n).$$

**Check 1 — the closed form.** `pure_type_log2` was validated against 15
brute-forced hom counts (all $(n,k,j)$ with $n\le 8$, $k\le3$), exact match.
Its optimal rank equals the predicted $j/2^{j+1}$ *exactly* in every case
tested ($j=1,2,3,4$; $n=64,128,256$), and the measured
$\mathrm{excess}/(n\log_2 n)$ approaches $1-2^{-j}$ with an $O(n)$ correction
whose size is confirmed constant after multiplying by $\log_2 n$
($j=1$: 0.687, 0.703, 0.713; $j=2$: 1.685, 1.706, 1.717).

**Check 2 — the orbit profile.** Exactly,
$\mathbb{E}[\#\text{orbits of size }m] = \binom{n}{m}a_m h_k(n-m)/h_k(n)$.
At $k=n/4$ the expected fraction of points in orbits of size 4 is

| $n$ | 64 | 128 | 256 | 512 |
| --- | --- | --- | --- | --- |
| size 4 | 0.735008 | 0.804914 | 0.857902 | 0.897348 |
| size 2 | 0.264883 | 0.195086 | 0.142098 | 0.102652 |
| size $\ge 8$ | $1.1\times10^{-4}$ | $4.0\times10^{-9}$ | $2.2\times10^{-18}$ | $2.6\times10^{-37}$ |

rising toward 1 with the residual entirely on size 2, exactly as the mixture
argument predicts. (`src/final_check.py` verifies $\sum_m m\,\mathbb{E}[\#\text{orbits}_m]=n$
to $10^{-9}$ as an arithmetic self-check.)

**Check 3 — out-of-sample extrapolation.** A direct fit is not enough:
$A$ converges slowly, and a basis chosen by in-sample fit proves nothing when
seven candidate shapes are available. So $A_{\text{local}}(n)$ was measured by
second differences (which annihilate the $Bn+C$ terms exactly), fitted on
$n\in[96,296]$ only, and scored on $n\in[400,496]$, which no fit ever saw:

| basis | $A_\infty$ | train max-res | **held-out max-res** |
| --- | --- | --- | --- |
| $n^{-1/2}$ | **0.74962** | 3.3e−5 | **6.8e−5** |
| $1/\log_2^2 n$ | 0.75598 | 1.7e−4 | 2.9e−4 |
| $1/\log_2 n$ | 0.78532 | 1.1e−4 | 5.0e−4 |
| $n^{-1/4}$ | 0.77239 | 3.3e−4 | 1.1e−3 |
| $\log_2 n/n$ | 0.74102 | 3.9e−4 | 1.1e−3 |
| $n^{-1/2}\log_2 n$ | 0.76416 | 4.5e−4 | 1.3e−3 |
| $n^{-1}$ | 0.73823 | 7.3e−4 | 1.9e−3 |

Scanning a free exponent $p$ in $A_\infty + cn^{-p}$ the same way selects
$p = 0.48$ and $A_\infty = \mathbf{0.75057}$, within $5.7\times10^{-4}$ of the
predicted $3/4$, held-out max residual $1.6\times10^{-5}$.

**A rejected hypothesis, kept.** The first version of the model predicted
$\varepsilon \sim c/\log_2 n$ for the fraction of points off 4-orbits, hence
$A_{\text{local}} = 3/4 - c/\log_2 n$. That is **false**: measured
$\varepsilon\log_2 n$ falls 1.589 → 0.924 over $n=64..512$ while
$\varepsilon\sqrt{n}$ only drifts 2.119 → 2.323, and the $1/\log_2 n$ basis
extrapolates to 0.785, missing the target by 0.035 despite a good in-sample
fit. The exponent $1/2$ is now selected by held-out error. **I have no
derivation of why it is $1/2$** — that is an open question here.

## 7. The gap in [RoTr25] Theorem 6, located exactly

Theorem 6 is the only published result addressing the second half of #1162:

> "Let µ be a real number in [0, 1/16) and let ν be a real number in
> [0, ½ − √3⁄4). Then a random subgroup of S_n has an elementary abelian
> 2-section of order at least 2^{µn}, and a Sylow 2-subgroup of order at
> least 2^{νn}."

Both constants come from the same place. From the proof (p. 35):

> "use Lemma 2.4(i) to see that the number of R-subgroups of a Sylow
> 2-subgroup of S_n is at most 4·2^{νn(n−νn)} … Since ν(1 − ν) < 1/16, the
> result follows"

and for the section, "Since µ < 1/16 the result follows from the lower bound
in Theorem 1". Solving $\nu(1-\nu) = 1/16$ gives
$\nu = \tfrac12 - \tfrac{\sqrt3}{4} = 0.06698729810778\ldots$, reproducing
their printed constant to 25 digits — so this is the binding step, not a
by-product.

**The lossy step.** Lemma 2.4(i) is applied to a Sylow 2-subgroup of $S_n$,
of order $2^{\ell}$ with $\ell = n - s_2(n) \approx n$. The bound
$\zeta_2 2^{k(\ell-k)}$ is attained only when the ambient group is elementary
abelian of rank $\ell$. But the maximal elementary abelian 2-subgroups of
$S_n$ have rank $\lfloor n/2\rfloor$, not $n$. Replacing $\ell \approx n$ by
$n/2$ turns the threshold equation into

$$\nu\left(\tfrac12 - \nu\right) = \tfrac{1}{16},$$

whose discriminant is **exactly zero**, with a double root at
$\nu = \tfrac14$. The same substitution sends µ's threshold
$\mu < 1/16$ to $\mu(\tfrac12-\mu) < \tfrac1{16}$, again a double root at 1/4.

That the barrier is a *double* root rather than an interval endpoint is the
signature of a tight extremal problem, and $1/4$ is precisely the rank at
which §3 finds the elementary abelian family concentrated.

**Conjecture 7.1.** Theorem 6 holds for all $\mu < 1/4$ and all $\nu < 1/4$,
and 1/4 is optimal for both.

Improvement factors: $4.000$ for µ, $3.732$ for ν.

**Conjecture 7.2 (the statistical theorem Erdős asked for).** For a uniformly
random $H \le S_n$ with $n \equiv 0 \pmod 4$,
$$\log_2|H| = \left(\tfrac14 + o(1)\right)n \quad\text{in probability.}$$

Posed along $n \equiv 0 \pmod 4$ deliberately: [RoTr25] Theorem 4 disproves
Kantor's conjecture only for $n \equiv 3 \pmod 4$ ("Let n be congruent to 3
modulo 4 … the probability that a random subgroup of $S_n$ is nilpotent is
bounded away from 1"), so $n \equiv 0 \pmod 4$ is the residue class where a
random subgroup may still be nilpotent, hence a 2-group by Theorem 5 ("A
random nilpotent subgroup of $S_n$ is a 2-group").

The $\ge$ direction of 7.2 is Conjecture 7.1. The $\le$ direction is not
addressed by anything in [RoTr25].

## 8. Evidence against, and what is not established

**The elementary abelian family may not be typical.** This is the load-bearing
gap. $\log e(n)$ and $\log f(n)$ agree to $o(n^2)$, but Theorem 1's upper
bound permits $f(n)/e(n)$ as large as $2^{\beta n^{3/2}}$, so elementary
abelian subgroups could be a vanishing fraction of all subgroups forever.
Everything in §3 is then a theorem about a measure-zero family, and
Conjecture 7.2 would not follow. **Nothing here rules that out.**

**The small-$n$ data points the other way.** From the exact lattices:

| $n$ | 4 | 5 | 6 | 7 | 8 |
| --- | --- | --- | --- | --- | --- |
| $e(n)/f(n)$ | 0.4667 | 0.2949 | 0.1863 | 0.1165 | **0.0809** |
| median $\|H\|$ | 4 | 6 | 8 | 12 | 16 |
| $\mathbb{E}[\log_2\|H\|]/n$ | 0.474 | 0.491 | 0.507 | 0.500 | **0.4999** |

The elementary abelian fraction is *falling*, and the mean of $\log_2|H|$ sits
at $n/2$, not $n/4$. This is not evidence for Conjecture 7.2; if anything it
runs against it. The honest reading is that $n \le 8$ is nowhere near the
asymptotic regime — the $n^2/16$ term only overtakes competing structures for
$n$ far beyond anything enumerable — and that small-$n$ lattice data simply
cannot discriminate here. It is recorded because omitting it would be
misleading.

**Not proved.** The limits in §3, §4 and §6 are established as exact
computations at each tested $n$ and as heuristic derivations, not as theorems.
No saddle-point analysis of $a(n,k)$ was carried out. Conjectures 7.1 and 7.2
are conjectures.

**Novelty not established.** A search of arXiv (metadata, not full text),
the Semantic Scholar citation list of 2503.05416, and a web search found no
statement of the $n/4$ constant. That citation list grew from zero to three
between the first check and the final one, which is a useful reminder that
this kind of negative result has a short shelf life. The three, re-fetched at
the end, are: *Proportion of Simple Subgroups in Finite Groups and Their
Applications* (arXiv:2606.17488), which studies
$\mathcal{V}(G)=\mathrm{Simp}(G)/|L(G)|$ — a different statistic on the same
lattice, and not an asymptotic for the order of a random subgroup;
*Homological Nielsen realization for the manifolds $\#_n\mathbb{CP}^2$*
(arXiv:2605.27537), unrelated; and *Groups having 12 cyclic subgroups*
(arXiv:2210.11788), which predates 2503.05416 and is presumably a
Semantic Scholar linkage artifact. None overlaps §7. Absence of search
results is weak evidence regardless; the arXiv `all:` field does not index
full text. Re-run `refs/fetch.sh` before relying on any of this.

**Open questions raised here.** (i) Why is the approach exponent in §6 equal
to $1/2$? (ii) Is $e(n)/f(n) \to 0$, and if so how fast? (iii) Can the
$\vartheta$-law of §3 be proved by a saddle-point argument?

## 9. Reproduction

```sh
bash scripts/reproduce.sh              # everything, in order (~35 min)
FAST=1 bash scripts/reproduce.sh       # ~60 s; skips n=8 lattice, n=10 brute
                                       # force, and the n>=256 high-precision rows
python src/final_check.py              # re-derive every number above; exit 0 = pass
```

`final_check.py` runs 118 checks and was validated by planting defects rather
than by reading it. Five plants, each producing exit status 1:

| planted defect | caught by |
| --- | --- |
| one stored standard deviation perturbed by $10^{-4}$ | the theta-law checks |
| `total_subgroups` for $n=8$ incremented by 1 | A005432 comparison |
| exponent $\binom{k-m}{2}\to\binom{k-m}{2}+1$ in the $q$-binomial inversion | 10 separate checks |
| $n=512$ `decimal_places` 35 → 12 | the high-precision sd floor |
| $n=512$ sd digits `…694` → `…999` | the 30-digit string comparison |

A clean pass of an unchallenged checker would prove nothing, and `final_check`
distinguishes *skip* from *pass*: under `FAST=1` it prints the held-out
extrapolation as skipped and says so in the summary line.

Timings on an idle 11-core arm64 machine, CPython 3.14.7. Ranges are over two
runs; single figures are one run and should be read as indicative.

| step | time |
| --- | --- |
| `subgroups 8` | 883 s |
| `subgroups 1..7` | 8.8–11.6 s |
| `asymptotics.py --n-max 512` | 218 s |
| `verify_elemab.py --n-max 9` | 43.9–50.8 s |
| `verify_elemab.py --n-max 10` | 397 s and 913 s on two runs — it is very sensitive to competing load |
| `sd_precision.py` (to $n=512$, 60 dp) | 115 s |
| `final_check.py` | 0.86–0.94 s |

$n=9$ was attempted (`subgroups 9`) and deliberately abandoned. After 11
minutes it had not finished the first class, whose cost is a full
$9! = 362880$-element coset scan; scaling the $n=8$ run by the $9\times$
larger permutation group and the $1.87\times$ larger class count puts the
whole run at four hours at best. It would also have re-verified a value
already listed in A005432 rather than extending anything, and the enumerator
is already validated against OEIS at all eight of $n=1,\ldots,8$.
`final_check.py` carries the expected $n=9$ values (1694723 subgroups, 554
classes), so the check will run automatically if anyone completes it.

## References

* **[RoTr25]** C. M. Roney-Dougal and G. Tracey, *Subgroups of symmetric
  groups: enumeration and asymptotic properties*, arXiv:2503.05416v1
  [math.GR], submitted 7 March 2025. <https://arxiv.org/abs/2503.05416>
  (v1 only as of the last fetch, no journal reference; 3 Semantic Scholar
  citations, none bearing on §7 — see §8.)
* **[Py93]** L. Pyber, *Enumerating finite groups of given order*,
  Ann. of Math. 137 (1993), 203–220. [RoTr25] reference 16.
* **[Va99]** Various, *Some of Paul's favorite problems*, booklet for the
  conference "Paul Erdős and his mathematics", Budapest, July 1999 — the
  source of #1162 (5.73). Not consulted directly; not available to me.
* OEIS [A005432](https://oeis.org/A005432) (subgroups of $S_n$),
  [A000638](https://oeis.org/A000638) (conjugacy classes of subgroups),
  [A000085](https://oeis.org/A000085) (involutions),
  [A006116](https://oeis.org/A006116) (Galois numbers).
