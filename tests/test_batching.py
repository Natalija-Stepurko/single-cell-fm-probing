import numpy as np

from scfm.batching import length_batches, run_sharded


def test_length_batches_budget_coverage_order():
    rng = np.random.default_rng(0)
    lengths = rng.integers(1, 500, 300)
    budget = 2_000
    batches = length_batches(lengths, budget)

    flat = [i for b in batches for i in b]
    assert sorted(flat) == list(range(len(lengths)))
    for b in batches:
        assert len(b) == 1 or max(lengths[b]) * len(b) <= budget
    firsts = [lengths[b[0]] for b in batches]
    assert firsts == sorted(firsts, reverse=True)
    assert lengths[flat[0]] == lengths.max()
    assert list(lengths[flat]) == sorted(lengths, reverse=True)


def test_length_batches_single_long_item():
    batches = length_batches(np.array([5000, 10, 10]), 1000)
    assert batches[0] == [0]
    assert sorted(i for b in batches for i in b) == [0, 1, 2]


def test_run_sharded_writes_and_resumes(tmp_path):
    lengths = np.random.default_rng(1).integers(1, 50, 120)
    calls = []

    def forward(b):
        calls.append(list(b))
        return np.stack([np.full(3, i, dtype=np.float32) for i in b])

    out = run_sharded(len(lengths), lengths, 3, tmp_path / "shards", forward, 200, 25)
    assert np.array_equal(out[:, 0], np.arange(120))
    shards = sorted((tmp_path / "shards").glob("shard_*.npz"))
    assert len(shards) > 1
    assert sum(len(np.load(f)["idx"]) for f in shards) == 120

    n_calls = len(calls)
    again = run_sharded(len(lengths), lengths, 3, tmp_path / "shards", forward, 200, 25)
    assert len(calls) == n_calls
    assert np.array_equal(again, out)
