#!/usr/bin/env bash
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
for m in Bio scipy sklearn pandas numpy torch; do
  $PY -c "import $m; print('$m', getattr($m,'__version__','?'))" 2>/dev/null || echo "$m MISSING"
done
nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo "no nvidia-smi"
df -h /root | tail -1
