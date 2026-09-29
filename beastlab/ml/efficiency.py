"""24→32→16→1 MLP, trained solely on synthetic labels. CPU only."""
import json
import hashlib
import random
from pathlib import Path
import torch
from torch import nn
from ..io import parse_dna

ROOT = Path(__file__).resolve().parents[2] / "models/efficiency_mlp"
LABEL = "Simulated training data — illustrative model, not clinical-grade"
FEATURE_VERSION = "v1:GC,16-overlapping-dinucleotide-fractions,4-quarter-position-weights,3-PAM-one-hot"
PAMS = ("NGG", "NAG", "OTHER")


def features(guide, pam="NGG"):
    g = parse_dna(guide, guide=True)
    if pam not in PAMS:
        raise ValueError("Unknown PAM category")
    gc = sum(b in "GC" for b in g) / 20
    di = [sum(g[i:i+2] == a+b for i in range(19))/19 for a in "ACGT" for b in "ACGT"]
    # Four position-sensitive quarter summaries: fixed base preferences and positional weights.
    weights = {"A": .2, "C": .8, "G": 1., "T": 0.}
    pos = [sum(weights[g[i]] * (i+1) for i in range(q*5, q*5+5)) /
           sum(i+1 for i in range(q*5, q*5+5)) for q in range(4)]
    return torch.tensor([gc] + di + pos + [float(pam == p) for p in PAMS], dtype=torch.float32)


class EfficiencyMLP(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(24, 32), nn.ReLU(), nn.Linear(32, 16), nn.ReLU(), nn.Linear(16, 1), nn.Sigmoid())

    def forward(self, x):
        return self.net(x)


def train(seed=42, epochs=160):
    """5,000 guides: 4,000 training / 1,000 held-out. No user inputs used."""
    torch.set_num_threads(1)
    rng = random.Random(seed)
    guides = ["".join(rng.choices("ACGT", k=20)) for _ in range(5000)]
    pams = [rng.choices(PAMS, weights=[.8, .15, .05])[0] for _ in guides]
    x = torch.stack([features(g, p) for g, p in zip(guides, pams)])
    targets = []
    for g, p in zip(guides, pams):
        gc = sum(b in "GC" for b in g)/20
        # Deliberately synthetic: broad, bounded efficiencies with irreducible noise.
        mean = .32 + .38*(1-2*abs(gc-.55)) + .09*(g[-1] == "G") - .22*("TTTT" in g)
        mean += .08*(g[3] == "C") - .18*(p == "NAG") - .35*(p == "OTHER")
        targets.append(min(.99, max(.01, mean + rng.gauss(0, .09))))
    y = torch.tensor(targets, dtype=torch.float32).reshape(-1, 1)
    with torch.random.fork_rng():
        torch.manual_seed(seed)
        model = EfficiencyMLP()
        opt = torch.optim.Adam(model.parameters(), lr=.008)
        history = []
        for epoch in range(epochs):
            model.train()
            opt.zero_grad()
            loss = nn.functional.mse_loss(model(x[:4000]), y[:4000])
            loss.backward(); opt.step()
            model.eval()
            with torch.no_grad():
                val = nn.functional.mse_loss(model(x[4000:]), y[4000:]).item()
            history.append({"epoch": epoch+1, "train_mse": loss.item(), "validation_mse": val})
        with torch.no_grad():
            radius = torch.quantile((model(x[4000:])-y[4000:]).abs(), .90).item() * 100
    meta = {"seed": seed, "n_guides": 5000, "train": 4000, "validation": 1000,
            "epochs": epochs, "features": FEATURE_VERSION, "label": LABEL,
            "band_radius_percent": radius, "band_method": "90th percentile absolute residual on simulated held-out guides; not biological confidence",
            "loss": history, "target_formula": "clipped GC/position/PAM/polyT heuristic plus Gaussian noise sd=0.09"}
    return model.eval(), meta


def save(model, meta):
    ROOT.mkdir(parents=True, exist_ok=True)
    torch.save(model.state_dict(), ROOT / "weights.pt")
    (ROOT / "metadata.json").write_text(json.dumps(meta, indent=2))


def load():
    torch.set_num_threads(1)
    model = EfficiencyMLP()
    model.load_state_dict(torch.load(ROOT / "weights.pt", map_location="cpu", weights_only=True))
    return model.eval(), json.loads((ROOT / "metadata.json").read_text())


def predict(model, meta, guide, pam="NGG"):
    with torch.no_grad():
        value = float(model(features(guide, pam)).item()*100)
    radius = meta["band_radius_percent"]
    weight_hash = hashlib.sha256(b"".join(t.detach().cpu().numpy().tobytes() for t in model.state_dict().values())).hexdigest()
    return {"model_tensor_sha256": weight_hash, "seed": meta["seed"], "features": FEATURE_VERSION,
            "efficiency_percent": value, "illustrative_band": [max(0., value-radius), min(100., value+radius)],
            "label": LABEL, "band_method": meta["band_method"]}


if __name__ == "__main__":
    save(*train())
