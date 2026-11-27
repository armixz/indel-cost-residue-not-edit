#!/usr/bin/env bash
set -e
VENV=/root/envs/prodomino
if command -v uv >/dev/null 2>&1; then UV=uv; else UV=/root/.local/bin/uv; fi
"$UV" venv --python 3.11 "$VENV" 2>&1 | tail -2
VPY="$VENV/bin/python"
"$UV" pip install --python "$VPY" -q \
  torch --index-url https://download.pytorch.org/whl/cu121 2>&1 | tail -3
"$UV" pip install --python "$VPY" -q fair-esm==2.0.0 "lightning>=2.3,<3" torchmetrics biopython numpy pandas scikit-learn 2>&1 | tail -3
"$VPY" - <<'PY'
import torch, esm, lightning, torchmetrics, Bio
print("torch", torch.__version__, "cuda", torch.cuda.is_available())
print("lightning", lightning.__version__, "torchmetrics", torchmetrics.__version__, "bio", Bio.__version__)
if torch.cuda.is_available():
    p = torch.cuda.get_device_properties(0)
    print("gpu", p.name, round(p.total_memory/1e9,1), "GB")
PY
