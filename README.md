# Erdős #1162 — the number of subgroups of $S_n$, and their order

**Question.** Give an asymptotic formula for the number of subgroups of $S_n$.
Is there a statistical theorem on their order? (Erdős–Turán, [Va99, 5.73].)
**Status: open.** Not solved here.

## ELI5

Take $n$ numbered cards. A *shuffle* is any way of rearranging them, and $S_n$
is the collection of all of them — for 8 cards that is $8! = 40\,320$ shuffles.

Some sub-collections of shuffles are self-contained: do any two of them in a
row and you land on something still in the collection, and anything you do can
be undone without leaving it. Those are the **subgroups**. "Leave everything
alone" on its own is one. All 40 320 together is another. The interesting ones
sit in between — for instance, every shuffle that disturbs only the first three
cards and leaves the other five alone.

Erdős and Turán asked two questions about them.

**1. How many subgroups are there?** For 8 cards, exactly 151 221. For 18
cards it is about $7.6\times10^{18}$, so actually listing them is hopeless
well before you get there and the goal is a formula for roughly how big the
number is. This half is now largely answered by other people: it is about
$2^{n^2/16}$. For 100 cards that is a number with 189 digits.

**2. How big is a typical subgroup?** Pick one of those 151 221 at random. How
many shuffles does it contain? Erdős wanted a "statistical theorem" — a rule
for what to expect. This half is barely touched. One published theorem
addresses it, from 2025.

**What is in this repository.** Mostly work on question 2. Counting all
subgroups is out of reach, so we counted one well-behaved family of them
exactly, out to 512 cards, and found that family is rigid: its members are
overwhelmingly concentrated at a single size, and the spread around that size
stays bounded as $n$ grows rather than widening. That suggests an answer to
question 2 — a typical subgroup should contain about $2^{n/4}$ shuffles — and
it locates a specific step in the 2025 proof where a factor of 4 is given away.

**The catch, and it is a real one.** Everything above concerns one family of
subgroups, and that family might be a vanishing sliver of all of them. At
$n\le8$, where we can check directly, its share is *shrinking* — from 47 % down
to 8 % — and the average subgroup there is nearer $2^{n/2}$ than $2^{n/4}$.
Nothing we have rules out that the pattern we measured describes a corner of
the problem rather than the whole of it. So this is a conjecture with evidence
behind it, not an answer.

**Objective.** Fix the ambiguity in the displayed statement, compute exactly
— to $n=512$, far past any tabulated range — the rank-stratified family of
subgroups that carries the whole first-order answer, identify the limiting
constants in closed form, and turn the vague second half of the question into
a precise conjecture with the gap to the published theorem located exactly.
The full-lattice runs ($n\le8$) are a *validation* of the enumerator against
OEIS, not an extension: A005432 is already tabulated to $n=18$.

**Inputs.** The problem statement from
[erdosproblems.com/1162](https://www.erdosproblems.com/1162); Roney-Dougal–Tracey,
arXiv:2503.05416; OEIS A005432, A000638, A000085, A006116. Run `refs/fetch.sh`
— nothing is redistributed here.

**Findings.**

* The logarithm in the displayed statement is **base 2** (§1) — the constant
  1/16 is meaningless otherwise, and the problem page never says so.
* The full subgroup lattice of $S_n$ recomputed for $n \le 8$ (151 221
  subgroups, 296 classes at $n=8$), matching A005432 and A000638 exactly.
* Elementary abelian 2-subgroups counted exactly to $n = 512$ — this family
  already carries the entire leading term $2^{n^2/16}$. Its rank is
  $n/4 + O(1)$ with an **absolutely bounded** spread: the limiting law is
  $\mathbb{P}(n/4+j) \propto 2^{-j^2}$, and the standard deviation of the
  exact rank distribution agrees with the closed form
  $0.849305961022908\ldots$ to **35 decimal places** at $n=512$.
* A single maximal elementary abelian gives only $n^2/16 + 2.8820499654\ldots$
  — a *bounded* excess, with the constant $\log_2(\vartheta/\eta)$ in closed
  form. So the $n\log n$ gain in [RoTr25] Theorem 1 comes entirely from
  multiplicity.
* The second-order term is $\tfrac34 n\log_2 n$, derived from orbit structure
  (orbits of size 4 beat size 2 in a tie at the leading constant) and
  confirmed by **out-of-sample** extrapolation: $A_\infty = 0.75057$.
* [RoTr25] Theorem 6 — the only published result on the *order* — is not
  optimal. Its constants $1/16$ and $\tfrac12-\tfrac{\sqrt3}{4}$ both come
  from applying Lemma 2.4 with $\ell \approx n$; the elementary abelian
  picture replaces that with $n/2$, giving a threshold equation with
  **discriminant exactly zero** and a double root at $\mathbf{1/4}$.
  Improvement factors 4.000 and 3.732.

**Stated plainly:** the evidence for the last point is about elementary
abelian subgroups only, and they could be a vanishing fraction of all
subgroups. At $n\le8$ that fraction is *falling* (0.47 → 0.08) and
$\mathbb{E}[\log_2|H|]/n \approx 0.50$, not 0.25. See NOTES.md §8.

Read [`NOTES.md`](NOTES.md). Reproduce with `bash scripts/reproduce.sh`.

| path | contents |
|---|---|
| `src/elemab.py` | exact $a(n,k)$ by exponential formula + $q$-binomial inversion |
| `src/verify_elemab.py` | brute-force permutation check, shares no code with it |
| `src/subgroups.c` | full subgroup lattice of $S_n$, $n\le8$ (C, conjugacy-class expansion) |
| `src/asymptotics.py` | rank distribution, Galois numbers, the excess fit |
| `src/sd_precision.py` | the rank sd against its closed form, exact at 60 dp |
| `src/orbit_model.py` | orbit-structure model for the $n\log n$ term; held-out extrapolation |
| `src/final_check.py` | re-derives every claim; exit 0 = pass |
| `data/`, `results/` | exact lattices, rank counts, measured artifacts |

**Conventions.** All logarithms base 2. Needs a Python 3 with numpy and mpmath
(developed on 3.14.7); `scripts/reproduce.sh` takes `$PY` if set and otherwise
falls back to `python3` on `$PATH`. The C enumerator needs a C99 compiler; it
uses NEON intrinsics on arm64 and falls back to portable code elsewhere.

## Attribution

The question as displayed at the top is quoted from
[erdosproblems.com/1162](https://www.erdosproblems.com/1162), maintained by
Thomas Bloom, whose companion repository
[teorth/erdosproblems](https://github.com/teorth/erdosproblems) is licensed
Apache 2.0. The original attribution there is to Erdős and Turán via
[Va99] — *Some of Paul's favorite problems*, the booklet produced for the
conference "Paul Erdős and his mathematics", Budapest, July 1999 — problem
5.73.

No third-party paper, OEIS page, or API response is redistributed here.
`refs/fetch.sh` retrieves them; `refs/README.md` records which claim each one
supports.

The code, data, and prose in this repository are MIT licensed — see
[`LICENSE`](LICENSE). That covers this work only: the quoted question and any
material `refs/fetch.sh` downloads remain under their own terms.
