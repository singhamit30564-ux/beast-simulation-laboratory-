import numpy as np
import pytest
from beastlab.rl.guide_env import GuideEnv, mutate


def test_env_reset():
    env = GuideEnv(lambda g: g.count('G'))
    a, ia = env.reset(seed=42)
    env.step(0)
    b, ib = env.reset(seed=42)
    assert np.array_equal(a, b) and ia == ib
    assert a.shape == (20, 4) and np.all(a.sum(axis=1) == 1)


def test_action_validity():
    for a in range(80):
        assert len(mutate('A'*20, a)) == 20
    for a in [-1, 80, 1.5, True]:
        with pytest.raises(ValueError):
            mutate('A'*20, a)


def test_reward_sign():
    env = GuideEnv(lambda g: g.count('G')*5, 'A'*20)
    env.reset()
    assert env.step(2)[1] == 5
    assert env.step(0)[1] == -5
    assert env.step(0)[1] == 0
    for _ in range(7):
        result = env.step(0)
    assert result[3]
    with pytest.raises(RuntimeError):
        env.step(0)
