"""Automated checks: split disjointness, equivalence to v1 folds, conformal ranks, finite-sample coverage.
Run: python -m pytest -q code_v2/tests   (from the repository root)"""
import os, sys
import numpy as np
import pytest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "code_v2"))
import swap_engine_v2 as V2

V1_SRC = os.path.join(ROOT, "original_v1", "Source_Code")


@pytest.mark.parametrize("n", [692, 781, 859, 967, 1101])
@pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
def test_v2_split_disjoint_and_complete(n, seed):
    seen_test = np.zeros(n, int)
    for f, tr, ca, rc, te in V2.folds_v2(n, seed):
        parts = [set(tr), set(ca), set(rc), set(te)]
        for i in range(4):
            for j in range(i + 1, 4):
                assert not (parts[i] & parts[j]), f"overlap between partitions {i},{j}"
        assert set().union(*parts) == set(range(n))
        assert len(rc) >= 160, "recalibration pool must allow k=160"
        seen_test[te] += 1
    assert (seen_test == 1).all()


@pytest.mark.parametrize("seed", [0, 3])
def test_v1_folds_identical_to_original(seed):
    sys.path.insert(0, V1_SRC)
    import swap_engine as E1
    for (f1, tr1, ca1, te1), (f2, tr2, ca2, rc2, te2) in zip(E1.folds(859, seed), V2.folds_v1(859, seed)):
        assert f1 == f2 and np.array_equal(tr1, tr2) and np.array_equal(ca1, ca2) and np.array_equal(te1, te2)
        assert len(rc2) == 0


def test_v2_same_test_subjects_as_v1():
    for (_, _, _, _, te1), (_, _, _, _, te2) in zip(V2.folds_v1(900, 2), V2.folds_v2(900, 2)):
        assert np.array_equal(np.sort(te1), np.sort(te2))


@pytest.mark.parametrize("n,level,rank", [(19, 0.95, 19), (20, 0.95, 20), (40, 0.95, 39), (80, 0.95, 77),
                                          (160, 0.95, 153), (150, 0.90, 136), (10, 0.95, 11)])
def test_conformal_rank(n, level, rank):
    assert V2.conformal_rank(n, level) == rank


def test_conformal_quantile_is_rth_order_statistic_or_inf():
    rng = np.random.default_rng(0)
    for n in [10, 19, 20, 40, 80, 150, 160]:
        s = rng.normal(size=n)
        r = V2.conformal_rank(n, 0.95)
        q = V2.conformal_quantile(s, 0.95)
        if r > n:
            assert np.isinf(q)
        else:
            assert q == np.sort(s)[r - 1]


def test_v1_quantile_is_one_rank_higher_unless_capped():
    rng = np.random.default_rng(1)
    for n in [20, 40, 80, 139, 160, 194]:
        for level in (0.90, 0.95):
            s = rng.normal(size=n)
            r = V2.conformal_rank(n, level)
            used = V2.v1_rank_used(n, level)
            assert V2.v1_quantile(s, level) == np.sort(s)[used - 1]
            assert used == min(n, r + 1)


@pytest.mark.parametrize("k", [20, 40, 80])
@pytest.mark.parametrize("rule", ["textbook", "v1"])
def test_two_sided_recalibration_coverage_matches_theory(k, rule):
    """Exchangeable simulation: separate lower/upper corrections at alpha/2 each."""
    rng = np.random.default_rng(k)
    Q = V2.QUANT[rule]
    hits = []
    for _ in range(4000):
        e = rng.normal(size=k + 1)              # signed residuals; last one is the test point
        lo = Q(-e[:k], 0.95); hi = Q(e[:k], 0.95)
        hits.append((-e[k] <= lo) and (e[k] <= hi))
    assert abs(np.mean(hits) - V2.expected_recal_coverage(k, rule=rule)) < 0.012


def test_split_conformal_marginal_coverage():
    rng = np.random.default_rng(5)
    n = 150
    hits = [abs(x[-1]) <= V2.conformal_quantile(np.abs(x[:n]), 0.90) for x in rng.normal(size=(6000, n + 1))]
    assert abs(np.mean(hits) - V2.conformal_rank(n, 0.90) / (n + 1)) < 0.01
