"""Run the released ProDomino checkpoint on an arbitrary FASTA. See 03_prodomino.py.

Usage: python 05_prodomino_fasta.py <in.fasta> <checkpoint.ckpt> <out.csv>
"""
import sys
import numpy as np
import pandas as pd
import torch
import torch.nn as nn


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main(fa, ckpt, out):
    import esm
    seqs = read_fasta(fa)
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    model = model.half().to(dev).eval()
    bc = alphabet.get_batch_converter()
    sd = torch.load(ckpt, map_location="cpu", weights_only=False)
    sd = sd.get("state_dict", sd)
    w = {k.replace("model.model.", "").replace("model.", "", 1): v for k, v in sd.items()}
    head = nn.Sequential(nn.Linear(2560, 1280), nn.ReLU(), nn.Linear(1280, 1))
    head.load_state_dict({k: v for k, v in w.items() if k in head.state_dict()}, strict=True)
    head = head.to(dev).float().eval()
    rows = []
    for name, s in seqs.items():
        _, _, toks = bc([(name, s)])
        with torch.no_grad():
            rep = model(toks.to(dev), repr_layers=[36], return_contacts=False)["representations"][36][0, 1:-1].float()
            sc = torch.sigmoid(head(rep)).squeeze(-1).cpu().numpy()
        assert len(sc) == len(s)
        for i, (aa, v) in enumerate(zip(s, sc), start=1):
            rows.append({"name": name, "pos": i, "aa": aa, "prodomino": float(v)})
        print(name, "L=", len(s), "mean=%.3f" % sc.mean(), flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    main(*sys.argv[1:4])
