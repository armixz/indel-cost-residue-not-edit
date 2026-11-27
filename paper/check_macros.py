"""Cross-check: every generated macro used in the manuscript is defined, and vice versa.

The build already fails loudly on an undefined control sequence, so this is a second,
cheaper guard that also catches the opposite problem: a macro that update_numbers.py still
emits but the text no longer uses, which usually means a claim was deleted and its number
left behind.
"""
import io
import re
import sys

tex = io.open("main.tex", encoding="utf-8").read()
num = io.open("numbers.tex", encoding="utf-8").read()

used = set(re.findall(r"\\([A-Z][A-Za-z]*)\{\}", tex))
defined = set(re.findall(r"\\newcommand\{\\([A-Za-z]+)\}", num))

missing = sorted(used - defined)
unused = sorted(defined - used)

print(f"used in main.tex : {len(used)}")
print(f"defined          : {len(defined)}")
print(f"USED BUT UNDEFINED: {missing if missing else 'none'}")
print(f"defined but unused: {len(unused)}")
if unused:
    print("  " + ", ".join(unused))
sys.exit(1 if missing else 0)
