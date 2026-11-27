"""Run the released ProDomino checkpoint on the four DI_screen proteins.

The published model is ESM-2 3B (layer 36, per-residue) -> Linear(2560,1280) -> ReLU ->
Linear(1280,1) -> sigmoid. The checkpoint weights are loaded directly, so no Lightning version
has to match. Nothing is trained or tuned here.

Usage: python 03_prodomino.py <DI_screen/analysis/input_data> <ProDomino/checkpoints/main_checkpoint.ckpt> <out.csv>
"""
import sys

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

PROTEINS = ["AraC", "Flp", "SigF", "TVMV"]


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]
            seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def load_head(ckpt_path):
    ck = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    sd = ck.get("state_dict", ck)
    want = {}
    for k, v in sd.items():
        kk = k.split("model.", 1)[-1] if k.startswith("model.") else k
        kk = kk.split("model.", 1)[-1] if kk.startswith("model.") else kk
        want[kk] = v
    head = nn.Sequential(nn.Linear(2560, 1280), nn.ReLU(), nn.Linear(1280, 1))
    try:
        head.load_state_dict(want, strict=True)
    except Exception as e:
        print("direct load failed:", e)
        print("checkpoint keys:", [(k, tuple(v.shape)) for k, v in sd.items() if hasattr(v, "shape")])
        raise
    head.eval()
    return head


def main(src, ckpt, out_path):
    import esm
    seqs = read_fasta(src + "/proteins.fasta")
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    print("loading ESM-2 3B (downloads ~11 GB on first run)", flush=True)
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    model = model.half().to(dev).eval()
    bc = alphabet.get_batch_converter()
    head = load_head(ckpt).to(dev).float()

    rows = []
    for p in PROTEINS:
        s = seqs[p]
        _, _, toks = bc([(p, s)])
        with torch.no_grad():
            out = model(toks.to(dev), repr_layers=[36], return_contacts=False)
            emb = out["representations"][36][0, 1:-1].float()
            score = torch.sigmoid(head(emb)).squeeze(-1).cpu().numpy()
        assert len(score) == len(s), (p, len(score), len(s))
        for i, (aa, sc) in enumerate(zip(s, score)):
            rows.append({"protein": p, "res0": i, "aa": aa, "prodomino": float(sc)})
        print(f"{p}: L={len(s)} mean={score.mean():.3f} max={score.max():.3f} "
              f"top5={np.argsort(-score)[:5].tolist()}", flush=True)
        del emb, out
        torch.cuda.empty_cache()
    pd.DataFrame(rows).to_csv(out_path, index=False)
    print("wrote", out_path)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
