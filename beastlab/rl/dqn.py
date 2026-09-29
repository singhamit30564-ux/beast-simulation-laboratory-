"""Small DQN, imitation warm-start, bounded replay and honest comparisons."""
import hashlib
import json
import random
from collections import deque
from pathlib import Path
import numpy as np
import torch
from torch import nn
from .guide_env import GuideEnv, mutate
from ..ml.efficiency import features, load as load_efficiency

ROOT = Path(__file__).resolve().parents[2] / "models/guide_dqn"


class DQN(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(80, 64), nn.ReLU(), nn.Linear(64, 80))

    def forward(self, x):
        return self.net(x.flatten(start_dim=-2))


def scorer(model):
    def score(g):
        with torch.no_grad():
            return model(features(g)).item()*100
    return score


def best_action(guide, score):
    return max(range(80), key=lambda a: score(mutate(guide, a)))


def train(seed=42, episodes=120):
    torch.set_num_threads(1)
    rng = random.Random(seed)
    mlp, _ = load_efficiency()
    score = scorer(mlp)
    with torch.random.fork_rng():
        torch.manual_seed(seed)
        agent = DQN()
        optimizer = torch.optim.Adam(agent.parameters(), lr=.002)
        # Imitate greedy trajectories, ending at a heuristic high-scoring guide.
        states, labels = [], []
        good = None
        for _ in range(20):
            env = GuideEnv(score)
            state, _ = env.reset(seed=rng.randrange(100000))
            for _ in range(10):
                action = best_action(env.guide, score)
                states.append(state); labels.append(action)
                state, _, done, trunc, _ = env.step(action)
                if done or trunc:
                    break
            good = env.guide
        x = torch.tensor(np.stack(states))
        y = torch.tensor(labels)
        for _ in range(60):
            loss = nn.functional.cross_entropy(agent(x), y)
            optimizer.zero_grad(); loss.backward(); optimizer.step()
        target = DQN(); target.load_state_dict(agent.state_dict())
        memory = deque(maxlen=1000)
        returns = []
        for episode in range(episodes):
            env = GuideEnv(score, good if episode % 5 == 0 else None)
            state, _ = env.reset(seed=rng.randrange(100000))
            epsilon = max(.05, .95 * .97**episode)
            total = 0.
            for _ in range(10):
                with torch.no_grad():
                    action = rng.randrange(80) if rng.random() < epsilon else int(agent(torch.tensor(state)).argmax())
                nxt, reward, done, trunc, _ = env.step(action)
                memory.append((state, action, reward/100, nxt, done or trunc))
                total += reward
                if len(memory) >= 32:
                    batch = rng.sample(list(memory), 32)
                    s, a, r, ns, d = zip(*batch)
                    q = agent(torch.tensor(np.stack(s))).gather(1, torch.tensor(a)[:, None]).squeeze(1)
                    with torch.no_grad():
                        expected = torch.tensor(r) + .95 * target(torch.tensor(np.stack(ns))).max(1).values * (1-torch.tensor(d, dtype=torch.float32))
                    loss = nn.functional.smooth_l1_loss(q, expected)
                    optimizer.zero_grad(); loss.backward()
                    nn.utils.clip_grad_norm_(agent.parameters(), 1.)
                    optimizer.step()
                state = nxt
                if done or trunc:
                    break
            if episode % 10 == 0:
                target.load_state_dict(agent.state_dict())
            returns.append(total)
    from ..ml.efficiency import ROOT as MLP_ROOT
    meta = {"seed": seed, "episodes": episodes, "replay_capacity": 1000, "epsilon": "max(.05,.95*.97**episode)",
            "warm_start": "60 imitation epochs on 20 greedy trajectories; every fifth episode starts at heuristic guide",
            "reward": "MLP(new)-MLP(old), percentage points; divided by 100 for DQN training",
            "mlp_weights_sha256": hashlib.sha256((MLP_ROOT / "weights.pt").read_bytes()).hexdigest(),
            "returns": returns, "label": "Educational simulation; mutations do not preserve genomic target binding"}
    return agent.eval(), meta


def save(agent, meta):
    ROOT.mkdir(parents=True, exist_ok=True)
    torch.save(agent.state_dict(), ROOT / "weights.pt")
    (ROOT / "metadata.json").write_text(json.dumps(meta, indent=2))


def load():
    agent = DQN()
    agent.load_state_dict(torch.load(ROOT / "weights.pt", map_location="cpu", weights_only=True))
    return agent.eval()


def walkthrough(agent, score, guide, policy="DQN", seed=42):
    env = GuideEnv(score, guide)
    state, info = env.reset(seed=seed)
    before = info["efficiency"]
    rng = random.Random(seed)
    rows = []
    for step in range(10):
        if policy == "random":
            action = rng.randrange(80)
        elif policy == "greedy":
            action = best_action(env.guide, score)
        else:
            with torch.no_grad():
                action = int(agent(torch.tensor(state)).argmax())
        pos, base = divmod(action, 4)
        mutation = f"{pos+1}: {env.guide[pos]} → {'ACGT'[base]}"
        state, reward, done, trunc, info = env.step(action)
        rows.append({"step": step+1, "mutation": mutation, "reward (pp)": reward, "efficiency %": info["efficiency"], "guide": info["guide"]})
        if done or trunc:
            break
    return {"policy": policy, "before": before, "after": info["efficiency"], "final_guide": info["guide"], "steps": rows}


def benchmark(agent, score, guide):
    return [walkthrough(agent, score, guide, p) for p in ("DQN", "random", "greedy")]


if __name__ == "__main__":
    save(*train())
