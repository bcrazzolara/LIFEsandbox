"""Flask debug UI: draw a random sample, show its signal image and the
model's predicted vs. true parameters."""
import base64
import io

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import torch
from flask import Flask, render_template_string

from sbw import (transmission_scenery, sample_params, build_scene,
                 signal_from_scene, build_model, prep_image, denorm_targets,
                 NAMES, lambdas, t)

app = Flask(__name__)

T_FLAT = transmission_scenery()
RNG = np.random.default_rng()
CKPT = torch.load("model.pt", map_location="cpu")
MEAN, STD = CKPT["mean"], CKPT["std"]
MODEL = build_model()
MODEL.load_state_dict(CKPT["model"])
MODEL.eval()

PAGE = """
<!doctype html><title>SBW debug</title>
<style>body{font-family:sans-serif;margin:2em}table{border-collapse:collapse}
td,th{border:1px solid #ccc;padding:4px 10px;text-align:right}
th:first-child,td:first-child{text-align:left}</style>
<h2>SBW ConvNext debug</h2>
<a href="/">↻ new random sample</a>
<p><img src="data:image/png;base64,{{img}}"></p>
<table><tr><th>param</th><th>true</th><th>predicted</th><th>|err|</th></tr>
{% for n,tv,pv in rows %}<tr><td>{{n}}</td><td>{{tv}}</td><td>{{pv}}</td>
<td>{{'%.3f'|format((tv|float-pv|float)|abs)}}</td></tr>{% endfor %}</table>
"""


@app.route("/")
def index():
    p = sample_params(RNG)
    sig = signal_from_scene(build_scene(p), T_FLAT)
    with torch.no_grad():
        pred = denorm_targets(MODEL(prep_image(sig, MEAN, STD)[None]).numpy()[0])

    fig, ax = plt.subplots(figsize=(6, 3))
    ax.imshow(sig, aspect="auto", origin="lower", cmap="inferno",
              extent=(t[0] / 3600, t[-1] / 3600, lambdas[0] * 1e6, lambdas[-1] * 1e6))
    ax.set_xlabel("time [h]"); ax.set_ylabel(r"$\lambda$ [$\mu$m]")
    ax.set_title("signal_both_sbw")
    buf = io.BytesIO(); fig.tight_layout(); fig.savefig(buf, format="png"); plt.close(fig)
    img = base64.b64encode(buf.getvalue()).decode()

    rows = [(n, f"{tv:.3f}", f"{pv:.3f}") for n, tv, pv in zip(NAMES, p, pred)]
    return render_template_string(PAGE, img=img, rows=rows)


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=True)
