"""Gymnasium-style reset/step API without an extra Gym dependency."""
import random
import numpy as np
from ..io import parse_dna

BASES = "ACGT"


def mutate(guide, action):
    if not isinstance(action, (int, np.integer)) or isinstance(action, bool) or not 0 <= action < 80:
        raise ValueError("Action must be an integer 0..79 (position*4 + base)")
    pos, base = divmod(int(action), 4)
    return guide[:pos] + BASES[base] + guide[pos+1:]


class GuideEnv:
    action_count = 80
    observation_shape = (20, 4)

    def __init__(self, score, guide=None):
        self.score = score
        self.initial = parse_dna(guide, guide=True) if guide is not None else None
        self.guide = None

    def observation(self):
        return np.eye(4, dtype=np.float32)[[BASES.index(b) for b in self.guide]]

    def reset(self, *, seed=42, options=None):
        rng = random.Random(seed)
        self.guide = parse_dna((options or {}).get("guide", self.initial or "".join(rng.choices(BASES, k=20))), guide=True)
        self.steps = 0
        self.done = False
        return self.observation(), {"guide": self.guide, "efficiency": self.score(self.guide)}

    def step(self, action):
        if self.guide is None or self.done:
            raise RuntimeError("Reset before stepping a new episode")
        new = mutate(self.guide, action)
        old_score = self.score(self.guide)
        value = self.score(new)
        self.guide = new
        self.steps += 1
        terminated, truncated = value > 90, self.steps >= 10
        self.done = terminated or truncated
        return self.observation(), value-old_score, terminated, truncated, {"guide": new, "efficiency": value}
