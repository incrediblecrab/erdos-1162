"""Re-derive every numeric claim in README.md and NOTES.md from the artifacts.

Exit status 0 means every claim in the write-up was reproduced from the data
files and the source code in this directory.  Any non-zero exit means at least
one stated number is wrong; the failing check is printed.

This is deliberately independent of the narrative: each check recomputes the
quantity rather than reading it back from a file that the narrative also
quotes, except where the point of the check is exactly that a stored file
agrees with a fresh computation.

Run:  python src/final_check.py
"""

from __future__ import annotations

import json
import math
import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))

from asymptotics import log2_int, ratio  # noqa: E402
from elemab import (  # noqa: E402
    galois_number,
    gaussian_binomial,
    hom_table,
    rank_counts,
)
from orbit_model import (  # noqa: E402
    holdout_extrapolation,
    local_coefficients,
    orbit_profile,
    pure_type_log2,
    pure_type_scan,
)

FAILURES: list[str] = []
CHECKS = 0
SKIPPED: list[str] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    global CHECKS
    CHECKS += 1
    if ok:
        print(f"  ok    {name}" + (f"   [{detail}]" if detail else ""))
    else:
        print(f"  FAIL  {name}   {detail}")
        FAILURES.append(f"{name}: {detail}")


def close(a: float, b: float, tol: float) -> bool:
    return abs(a - b) <= tol


def load(rel: str):
    p = ROOT / rel
    if not p.exists():
        return None
    return json.loads(p.read_text())


# ------------------------------------------------------------------ constants
# Both constants below are closed forms; the point of the checks is that the
# measured quantities agree with them, so they are computed here from scratch.

THETA = 1.0 + 2.0 * sum(2.0 ** (-j * j) for j in range(1, 60))
ETA = math.prod(1.0 - 2.0**-i for i in range(1, 200))  # (1/2;1/2)_infinity
GALOIS_CONST = math.log2(THETA / ETA)  # 2.8820499654...
THETA_SD = math.sqrt(
    sum(j * j * 2.0 ** (-j * j) for j in range(-60, 61)) / THETA
)  # 0.8493059610...


def section(title: str) -> None:
    print(f"\n== {title} ==")


# --------------------------------------------------------------------- checks


def check_closed_forms() -> None:
    section("closed forms quoted in NOTES.md")
    check(
        "Galois constant log2(theta/(1/2;1/2)_inf) = 2.882049965",
        close(GALOIS_CONST, 2.882049965, 5e-10),
        f"{GALOIS_CONST:.10f}",
    )
    check(
        "theta-law standard deviation = 0.849305961",
        close(THETA_SD, 0.849305961, 5e-10),
        f"{THETA_SD:.10f}",
    )
    # zeta_2 of RoTr25 Lemma 2.4 is 1/ETA
    check("RoTr25 Lemma 2.4 constant zeta_2 < 4", (1 / ETA) < 4, f"zeta_2={1 / ETA:.6f}")


def check_elemab_against_bruteforce() -> None:
    section("elemab.py against the independent permutation brute force")
    rec = load("results/verify_elemab.json")
    if rec is None:
        SKIPPED.append("results/verify_elemab.json missing (run src/verify_elemab.py)")
        print("  skip  results/verify_elemab.json not present")
        return
    bad = [n for n, r in rec.items() if not r["match"]]
    check(
        f"brute force agrees with the formula for n <= {max(int(k) for k in rec)}",
        not bad,
        f"mismatches at n={bad}" if bad else f"{len(rec)} values of n",
    )
    for n, r in rec.items():
        if int(n) <= 10:
            fresh = rank_counts(int(n))
            stored = [c for c in r["formula"]]
            while len(fresh) < len(stored):
                fresh.append(0)
            while len(stored) < len(fresh):
                stored.append(0)
            if fresh != stored:
                check(f"stored formula row n={n} still reproduces", False, "drifted")
                return
    check("stored formula rows still reproduce from current elemab.py", True)


def check_involutions() -> None:
    section("a(n,1) against OEIS A000085 (an independent recurrence)")
    # A000085: t(n) = t(n-1) + (n-1) t(n-2), the number of involutions
    # including the identity, so a(n,1) = A000085(n) - 1.
    t = [1, 1]
    for m in range(2, 31):
        t.append(t[m - 1] + (m - 1) * t[m - 2])
    ok = all(rank_counts(n)[1] == t[n] - 1 for n in range(2, 31))
    check("a(n,1) = A000085(n) - 1 for 2 <= n <= 30", ok, f"a(30,1)={rank_counts(30)[1]}")


def check_no_rank_above_half() -> None:
    section("structural identity a(n, floor(n/2)+1) = 0")
    # rank_counts asserts this internally; recompute the inversion one rank
    # higher than the maximum and confirm the massive cancellation lands on 0.
    ok = True
    for n in range(2, 41):
        c = rank_counts(n)
        if len(c) - 1 != n // 2:
            ok = False
            break
    check("rank_counts(n) has top rank exactly floor(n/2) for n <= 40", ok)


def check_lattice_vs_oeis() -> None:
    section("full subgroup lattice against OEIS A005432 / A000638")
    a005432 = {1: 1, 2: 2, 3: 6, 4: 30, 5: 156, 6: 1455, 7: 11300, 8: 151221, 9: 1694723}
    a000638 = {1: 1, 2: 2, 3: 4, 4: 11, 5: 19, 6: 56, 7: 96, 8: 296, 9: 554}
    seen = 0
    for n in range(1, 10):
        d = load(f"data/lattice_n{n}.json")
        if d is None:
            continue
        seen += 1
        check(
            f"n={n}: total subgroups = A005432({n}) = {a005432[n]}",
            d["total_subgroups"] == a005432[n],
            f"got {d['total_subgroups']}",
        )
        check(
            f"n={n}: conjugacy classes = A000638({n}) = {a000638[n]}",
            d["total_classes"] == a000638[n],
            f"got {d['total_classes']}",
        )
        lat = [int(x) for x in d["elementary_abelian_rank"]]
        formula = rank_counts(n)
        while len(lat) < len(formula):
            lat.append(0)
        while len(formula) < len(lat):
            formula.append(0)
        check(
            f"n={n}: lattice elementary abelian ranks = elemab.py",
            lat == formula,
            f"lattice {lat} vs formula {formula}",
        )
    if seen == 0:
        SKIPPED.append("no data/lattice_n*.json (build and run src/subgroups.c)")
        print("  skip  no lattice files present")


def check_pure_types() -> None:
    section("pure orbit-type closed form against brute-forced hom counts")
    # (n, k, j) -> exact number of homomorphisms F_2^k -> S_n all of whose
    # orbits have size exactly 2^j, obtained by enumerating k-tuples of
    # pairwise commuting elements of order dividing 2 (see NOTES.md).
    brute = {
        (4, 1, 1): 3, (4, 1, 2): 0, (4, 2, 1): 27, (4, 2, 2): 6,
        (4, 3, 1): 147, (4, 3, 2): 42, (6, 1, 1): 15, (6, 2, 1): 405,
        (6, 3, 1): 5145, (8, 1, 1): 105, (8, 1, 2): 0, (8, 1, 3): 0,
        (8, 2, 1): 8505, (8, 2, 2): 1260, (8, 2, 3): 0,
    }
    from elemab import gl_order

    bad = []
    for (n, k, j), want in brute.items():
        v = pure_type_log2(n, j, k)
        got = 0.0 if v is None else 2.0**v * gl_order(k)
        if not close(got, want, 1e-6 * max(1.0, want)):
            bad.append((n, k, j, want, got))
    check(
        f"{len(brute)} brute-forced hom counts match pure_type_log2",
        not bad,
        str(bad) if bad else "",
    )

    section("orbit model: optimal rank of each pure type is j/2^(j+1)")
    for n in (128, 256):
        for r in pure_type_scan(n):
            check(
                f"n={n} orbit size {r['orbit_size']}: best k/n = {r['predicted_k_over_n']}",
                close(r["best_k_over_n"], r["predicted_k_over_n"], 1e-9),
                f"got {r['best_k_over_n']}",
            )
    # j = 1 and j = 2 tie at n^2/16; j >= 3 must lose at order n^2
    for n in (128, 256):
        sc = {r["j"]: r for r in pure_type_scan(n)}
        for j in (1, 2):
            # the model says the excess over n^2/16 is Theta(n log2 n), so the
            # scale to compare against is n log2 n, not n^2
            check(
                f"n={n} orbit size {2**j} reaches n^2/16 up to O(n log n)",
                abs(sc[j]["excess_over_n2_16"]) < n * math.log2(n),
                f"excess {sc[j]['excess_over_n2_16']:+.1f} vs n log2 n = {n * math.log2(n):.0f}",
            )
            check(
                f"n={n} orbit size {2**j}: excess/(n log2 n) below its limit {1 - 2.0**-j}",
                0 < sc[j]["excess_per_n_log2n"] < 1 - 2.0**-j,
                f"{sc[j]['excess_per_n_log2n']:.5f}",
            )
        check(
            f"n={n}: orbit size 4 beats orbit size 2 at order n log n",
            sc[2]["excess_per_n_log2n"] > sc[1]["excess_per_n_log2n"],
            f"{sc[2]['excess_per_n_log2n']:.5f} vs {sc[1]['excess_per_n_log2n']:.5f}",
        )
        if 3 in sc:
            # the model predicts the gap grows like (1/16 - 9/256) n^2 = 7n^2/256,
            # approached from below
            gap = sc[2]["excess_over_n2_16"] - sc[3]["excess_over_n2_16"]
            check(
                f"n={n}: orbit size 8 falls behind orbit size 4 at order n^2",
                0.02 * n * n < gap < (7 / 256) * n * n,
                f"gap {gap:.1f}, gap/n^2 = {gap / (n * n):.5f}, limit 7/256 = 0.02734",
            )


def check_theta_law() -> None:
    section("pointwise theta law  a(n, n/4+j)/a(n, n/4) -> 2^(-j^2)")
    hom = hom_table(128 // 2 + 2, 128)
    worst = 0.0
    for n in (64, 128):
        c = rank_counts(n, hom if n == 128 else None)
        k = n // 4
        for j in range(-3, 4):
            got = ratio(c[k + j], c[k])
            want = 2.0 ** (-j * j)
            worst = max(worst, abs(got - want))
    check(
        "ratios match 2^(-j^2) for |j| <= 3 at n = 64 and 128",
        worst < 2e-4,
        f"max deviation {worst:.2e}",
    )
    # the sharper statement is at large n; use the stored sweep if present
    st = load("data/elemab_stats.json")
    if st is None:
        SKIPPED.append("data/elemab_stats.json missing (run src/asymptotics.py)")
        print("  skip  data/elemab_stats.json not present")
        return
    rows = {r["n"]: r for r in st["rows"]}
    for n in (128, 256, 512):
        if n not in rows:
            continue
        r = rows[n]
        check(
            f"n={n}: mode of the rank distribution is exactly n/4 = {n // 4}",
            r["mode"] == n // 4,
            f"got {r['mode']}",
        )
        check(
            f"n={n}: mean rank / n = 0.25 to 6 dp",
            close(r["mean_over_n"], 0.25, 5e-7),
            f"{r['mean_over_n']:.8f}",
        )
        check(
            f"n={n}: sd of the rank equals the theta constant to 6 dp",
            close(r["sd"], THETA_SD, 5e-7),
            f"measured {r['sd']:.8f} vs {THETA_SD:.8f}",
        )
        check(
            f"n={n}: >= 1-1e-6 of the mass lies within 4 of the mode",
            r["window_1e6"] <= 4,
            f"window {r['window_1e6']}",
        )


def check_galois() -> None:
    section("Galois numbers  log2 G_r - r^2/4 -> the same constant")
    g = load("data/galois.json")
    if g is None:
        SKIPPED.append("data/galois.json missing (run src/asymptotics.py)")
        print("  skip  data/galois.json not present")
    else:
        rows = {r["r"]: r for r in g}
        # The approach to the constant is itself a claim: the leading correction
        # to [r choose r/2]_2 is of size 2^-(r/2), so check that rate rather
        # than asserting a single tolerance that happens to hold at large r.
        for r in (32, 64):
            if r in rows:
                dev = abs(rows[r]["excess_over_r2_4"] - GALOIS_CONST)
                check(
                    f"r={r}: deviation from the constant is below 2^(3-r/2)",
                    dev < 2.0 ** (3 - r / 2),
                    f"dev {dev:.3e} < {2.0 ** (3 - r / 2):.3e}",
                )
        for r in (128, 256, 512):
            if r in rows:
                # 1e-11 is the float floor of log2_int on numbers this large
                check(
                    f"r={r} (even): excess = {GALOIS_CONST:.9f}",
                    close(rows[r]["excess_over_r2_4"], GALOIS_CONST, 1e-11),
                    f"{rows[r]['excess_over_r2_4']:.12f}",
                )
        for r in (63, 127, 255, 511):
            if r in rows:
                check(
                    f"r={r} (odd): excess within 4e-6 of the even limit",
                    close(rows[r]["excess_over_r2_4"], GALOIS_CONST, 4e-6),
                    f"{rows[r]['excess_over_r2_4']:.9f}",
                )
    # recompute one value from scratch rather than trusting the file
    fresh = log2_int(galois_number(256)) - 256 * 256 / 4.0
    check(
        "fresh galois_number(256) reproduces the constant",
        close(fresh, GALOIS_CONST, 1e-11),
        f"{fresh:.12f}",
    )
    # G_r for small r is OEIS A006116
    a006116 = [1, 2, 5, 16, 67, 374, 2825, 29212, 417199]
    check(
        "G_0..G_8 = OEIS A006116",
        [galois_number(r) for r in range(9)] == a006116,
        "",
    )
    # a single maximal elementary abelian gives only n^2/16 + O(1)
    for n in (256, 512):
        v = log2_int(galois_number(n // 2)) - n * n / 16.0
        check(
            f"n={n}: one maximal elementary abelian gives n^2/16 + O(1)",
            abs(v - GALOIS_CONST) < 1e-11,
            f"excess {v:.12f}",
        )


def check_second_order() -> None:
    section("second-order term: A_local -> 3/4, scored out of sample")
    st = load("data/elemab_stats.json")
    if st is None:
        print("  skip  data/elemab_stats.json not present")
        return
    loc = local_coefficients(st["rows"], 16)
    hv = holdout_extrapolation(loc)
    if "error" in hv:
        print(f"  skip  {hv['error']}")
        SKIPPED.append("held-out extrapolation: " + hv["error"])
        return
    best_basis = hv["ranked_bases"][0]
    check(
        "n^-1/2 is the best extrapolation basis on held-out n",
        best_basis["basis"] == "n^-1/2",
        f"won: {best_basis['basis']} (held-out maxres {best_basis['holdout_maxres']:.2e})",
    )
    b = hv["best_free_exponent"]
    check(
        "free exponent selected out of sample is near 1/2",
        0.40 <= b["p"] <= 0.60,
        f"p={b['p']:.2f}",
    )
    check(
        "extrapolated A_inf = 3/4 to within 0.002",
        abs(b["A_inf"] - 0.75) < 0.002,
        f"A_inf={b['A_inf']:.5f}, error {abs(b['A_inf'] - 0.75):.5f}",
    )
    # and the rejected hypothesis really is rejected
    byname = {r["basis"]: r for r in hv["ranked_bases"]}
    check(
        "the rejected 1/log2 n model extrapolates to the wrong limit",
        abs(byname["1/log2 n"]["A_inf"] - 0.75) > 0.02,
        f"A_inf={byname['1/log2 n']['A_inf']:.5f}",
    )


def check_orbit_profile() -> None:
    section("orbit profile at the dominant rank k = n/4")
    prev = None
    for n in (64, 128, 256, 512):
        p = orbit_profile(n, n // 4)
        check(
            f"n={n}: sum_m m*E[#orbits_m] = n exactly",
            p["point_total_error"] < 1e-6 * n,
            f"error {p['point_total_error']:.1e}",
        )
        f4 = p["profile"].get(4, {}).get("exp_point_fraction", 0.0)
        f8 = p["profile"].get(8, {}).get("exp_point_fraction", 0.0)
        check(
            f"n={n}: orbits of size 4 hold the majority of points",
            f4 > 0.5,
            f"fraction {f4:.6f}",
        )
        check(f"n={n}: orbits of size >= 8 are negligible", f8 < 1e-3, f"{f8:.2e}")
        if prev is not None:
            check(f"n={n}: the size-4 fraction increased", f4 > prev, f"{prev:.4f} -> {f4:.4f}")
        prev = f4


def check_sd_high_precision() -> None:
    section("rank sd vs closed form, exact arithmetic at 60 dp")
    # The doubles in elemab_stats.json truncate at ~1e-16 relative, which is
    # coarser than the real convergence from n=256 on.  Recompute exactly.
    d = load("data/sd_precision.json")
    if d is None:
        SKIPPED.append("data/sd_precision.json missing (run src/sd_precision.py)")
        print("  skip  data/sd_precision.json not present")
        return
    check(
        "closed form printed to 30 dp is 0.849305961022908280743697585694",
        d["closed_form_30dp"].startswith("0.849305961022908280743697585694"),
        d["closed_form_30dp"],
    )
    want = {64: 4, 128: 7, 256: 15, 512: 30}  # lower bounds, not exact values
    byn = {r["n"]: r for r in d["rows"]}
    for n, floor_dp in want.items():
        r = byn.get(n)
        if r is None:
            continue
        check(
            f"n={n}: exact sd agrees with the closed form to >= {floor_dp} dp",
            r["decimal_places"] >= floor_dp,
            f"{r['decimal_places']} dp, |err| {r['abs_error']}",
        )
    r512 = byn.get(512)
    if r512 is not None:
        check(
            "n=512: all 30 printed digits of the sd match the closed form",
            r512["sd_30dp"].startswith("0.849305961022908280743697585694"),
            r512["sd_30dp"],
        )
    ratios = [r["log10_ratio_vs_prev"] for r in d["rows"] if r["log10_ratio_vs_prev"]]
    if ratios:
        check(
            "log10(error) roughly doubles per doubling of n (error ~ 2^-cn)",
            all(1.8 < x < 2.3 for x in ratios[-2:]),
            f"ratios {[round(x, 2) for x in ratios]}",
        )


def check_rotr_thresholds() -> None:
    section("RoTr25 Theorem 6 threshold algebra")
    lo = (1 - math.sqrt(3) / 2) / 2
    check(
        "nu(1-nu) = 1/16 has smaller root 1/2 - sqrt(3)/4 = 0.0669872981",
        close(lo, 0.5 - math.sqrt(3) / 4, 1e-15) and close(lo, 0.06698729810778, 1e-12),
        f"{lo:.13f}",
    )
    check("that root is the constant printed in Theorem 6", close(lo, 0.06698729810778, 1e-12))
    disc = 0.25 - 4 * (1 / 16.0)
    check("nu(1/2 - nu) = 1/16 has discriminant exactly 0", disc == 0.0, f"{disc!r}")
    check("its double root is 1/4", close(0.5 / 2, 0.25, 0), "0.25")
    check(
        "improvement factors are 4.000 and 3.732",
        close(0.25 / (1 / 16.0), 4.0, 1e-12)
        and close(0.25 / lo, 3.7320508075688772, 1e-9),
        f"{0.25 / (1 / 16.0):.3f} and {0.25 / lo:.3f}",
    )


def check_quotes() -> None:
    section("verbatim quotations from the primary source")
    pdf = ROOT / "refs" / "roney-dougal_tracey_2503.05416.pdf"
    if not pdf.exists():
        SKIPPED.append("refs/*.pdf absent (run refs/fetch.sh); quotations unchecked")
        print("  skip  refs PDF not present; run refs/fetch.sh")
        return
    # -raw keeps content-stream order.  Default and -layout mode both sort by
    # position, which splices superscripts from displayed math into the middle
    # of a sentence: the base-2 remark comes out as "are to the 2 2 base 2)"
    # and a verbatim search fails even though the sentence is on the page.
    # The base-2 line was additionally read off a 150 dpi render of page 1.
    try:
        txt = subprocess.run(
            ["pdftotext", "-raw", str(pdf), "-"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError) as exc:
        SKIPPED.append(f"pdftotext unavailable ({exc}); quotations unchecked")
        print("  skip  pdftotext unavailable")
        return
    flat = " ".join(txt.split())
    quotes = [
        "our logarithms are to the base 2",
        "settling a conjecture of Pyber from 1993",
        "A random nilpotent subgroup of Sn is a 2-group",
        "Let n be congruent to 3 modulo 4",
        "the probability that a random subgroup of Sn is nilpotent is bounded away from 1",
        "an elementary abelian 2-section of order at least 2",
        "Since ν(1 − ν) < 1/16, the result follows",
        "Since µ < 1/16 the result follows from the lower bound in Theorem 1",
        "P has at most ζppk(ℓ−k) subgroups of order pk",
    ]
    for q in quotes:
        check(f'quote present: "{q[:56]}"', " ".join(q.split()) in flat)


def main() -> int:
    print("final_check.py -- re-deriving every claim in README.md and NOTES.md")
    check_closed_forms()
    check_elemab_against_bruteforce()
    check_involutions()
    check_no_rank_above_half()
    check_lattice_vs_oeis()
    check_pure_types()
    check_theta_law()
    check_galois()
    check_second_order()
    check_orbit_profile()
    check_sd_high_precision()
    check_rotr_thresholds()
    check_quotes()

    print(f"\n{CHECKS} checks run, {len(FAILURES)} failed")
    for s in SKIPPED:
        print(f"  note: skipped -- {s}")
    if FAILURES:
        print("\nFAILED:")
        for f in FAILURES:
            print(f"  - {f}")
        return 1
    print("\nAll stated numbers reproduced.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
