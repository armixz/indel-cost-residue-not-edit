"""Gate 8: per-position ESM-2 3B features (entropy of the residue distribution, log-likelihood of
the wild-type residue) for every protein in the benchmark. One forward pass per protein.

Usage: python 09_esm_features.py <out.csv> <fasta> [<fasta> ...]
"""
import sys

import numpy as np
import pandas as pd
import torch


def read_fasta(path):
    seqs, name = {}, None
    for line in open(path):
        line = line.strip()
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = []
        elif name:
            seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main(out, *fastas):
    import esm
    dev = "cuda" if torch.cuda.is_available() else "cpu"
    model, alphabet = esm.pretrained.esm2_t36_3B_UR50D()
    model = model.half().to(dev).eval()
    bc = alphabet.get_batch_converter()
    std = [alphabet.get_idx(a) for a in "ACDEFGHIKLMNPQRSTVWY"]
    rows = []
    for fa in fastas:
        for name, s in read_fasta(fa).items():
            _, _, toks = bc([(name, s)])
            with torch.no_grad():
                lg = model(toks.to(dev), repr_layers=[], return_contacts=False)["logits"][0, 1:-1]
            lg = lg[:, std].float()
            lp = torch.log_softmax(lg, dim=-1)
            p = lp.exp()
            ent = (-(p * lp).sum(-1)).cpu().numpy()
            wt_idx = torch.tensor([("ACDEFGHIKLMNPQRSTVWY".index(a) if a in "ACDEFGHIKLMNPQRSTVWY" else 0)
                                   for a in s], device=lp.device)
            wt_ll = lp.gather(1, wt_idx[:, None]).squeeze(1).cpu().numpy()
            for i, a in enumerate(s, start=1):
                rows.append({"name": name, "pos": i, "aa": a,
                             "esm_entropy": float(ent[i - 1]), "esm_wt_loglik": float(wt_ll[i - 1])})
            print(f"{name}: L={len(s)} entropy mean={ent.mean():.3f} wt_ll mean={wt_ll.mean():.3f}", flush=True)
    pd.DataFrame(rows).to_csv(out, index=False)
    print("wrote", out)


if __name__ == "__main__":
    main(sys.argv[1], *sys.argv[2:])
