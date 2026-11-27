#!/usr/bin/env bash
VPY=/root/envs/prodomino/bin/python
D=/mnt/c/Users/armin/Desktop/AXZ/insertion-lab
exec "$VPY" "$D/analysis/$1" "${@:2}"
