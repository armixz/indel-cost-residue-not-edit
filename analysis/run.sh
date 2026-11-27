#!/usr/bin/env bash
# usage: run.sh <script.py> [args...]
PY=$(dirname "$(readlink -f /root/.local/bin/state)")/python
D=/mnt/c/Users/armin/Desktop/AXZ/insertion-lab
cd "$D/analysis" && exec "$PY" -I "$D/analysis/$1" "${@:2}"
