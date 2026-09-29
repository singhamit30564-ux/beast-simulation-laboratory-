import torch
from beastlab.ml.efficiency import EfficiencyMLP, features, load, predict, train


def test_shapes():
    assert features('ACGT'*5).shape == (24,)
    assert EfficiencyMLP()(torch.zeros(3, 24)).shape == (3, 1)


def test_determinism():
    a, ma = train(epochs=3)
    b, mb = train(epochs=3)
    assert ma == mb
    assert all(torch.equal(a.state_dict()[k], b.state_dict()[k]) for k in a.state_dict())


def test_range():
    model, meta = load()
    for guide in ['A'*20, 'G'*20, 'ACGT'*5, 'T'*20]:
        r = predict(model, meta, guide)
        assert 0 <= r['illustrative_band'][0] <= r['efficiency_percent'] <= r['illustrative_band'][1] <= 100
