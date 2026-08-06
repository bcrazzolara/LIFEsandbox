"""Train the ViT to regress the 6 SBW parameters from the signal image."""
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader

from sbw import build_model, norm_targets, denorm_targets, prep_image, NAMES

NGPU = torch.cuda.device_count()
EPOCHS = 50
BATCH = 64 * max(1, NGPU)   # scale batch across GPUs
LR = 3e-4
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"

# Up-weight the faint planet params (idx1, idx2, Amplitude_planet); normalized
# to mean 1 so the loss scale (and thus LR) stays comparable.
WEIGHTS = torch.tensor([2., 2., 3., 1., 1., 1.])
WEIGHTS = (WEIGHTS / WEIGHTS.mean()).to(DEVICE)


class SBWDataset(Dataset):
    def __init__(self, npz, mean, std):
        d = np.load(npz)
        self.X, self.Y = d["X"], d["Y"]
        self.mean, self.std = mean, std

    def __len__(self):
        return len(self.X)

    def __getitem__(self, i):
        x = prep_image(self.X[i], self.mean, self.std)
        y = torch.as_tensor(norm_targets(self.Y[i]), dtype=torch.float32)
        return x, y


def main():
    stats = np.load("data/stats.npz")
    mean, std = float(stats["mean"]), float(stats["std"])
    tr = DataLoader(SBWDataset("data/train.npz", mean, std), BATCH, shuffle=True,
                    num_workers=4, pin_memory=True)
    te = DataLoader(SBWDataset("data/test.npz", mean, std), BATCH,
                    num_workers=4, pin_memory=True)

    model = build_model().to(DEVICE)
    if NGPU > 1:
        model = torch.nn.DataParallel(model)
        print(f"DataParallel over {NGPU} GPUs")
    opt = torch.optim.AdamW(model.parameters(), lr=LR)
    sched = torch.optim.lr_scheduler.CosineAnnealingLR(opt, EPOCHS)

    for ep in range(EPOCHS):
        model.train()
        tot = 0.0
        for x, y in tr:
            x, y = x.to(DEVICE), y.to(DEVICE)
            opt.zero_grad()
            loss = (WEIGHTS * (model(x) - y) ** 2).mean()   # planet-weighted MSE
            loss.backward()
            opt.step()
            tot += loss.item() * len(x)
        sched.step()
        core = model.module if NGPU > 1 else model
        torch.save({"model": core.state_dict(), "mean": mean, "std": std},
                   "model.pt")  # checkpoint every epoch
        mae = evaluate(model, te)
        line = "  ".join(f"{n}={m:.3g}" for n, m in zip(NAMES, mae))
        print(f"epoch {ep + 1}/{EPOCHS}  train={tot / len(tr.dataset):.4f} | {line}")

    print("done; model.pt holds the last epoch")


@torch.no_grad()
def evaluate(model, loader):
    """Return per-param test MAE in real units."""
    model.eval()
    preds, trues = [], []
    for x, y in loader:
        preds.append(model(x.to(DEVICE)).cpu().numpy())
        trues.append(y.numpy())
    pred = denorm_targets(np.concatenate(preds))
    true = denorm_targets(np.concatenate(trues))
    return np.abs(pred - true).mean(0)


if __name__ == "__main__":
    main()
