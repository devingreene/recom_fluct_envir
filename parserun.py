#!/usr/bin/env python3
from pyparsing import *
import sys
import numpy as np
import itertools

def stdin_no_sharp():
    for line in sys.stdin:
        if line != "\n" and not line.strip().startswith("#"):
            yield line

stdin = stdin_no_sharp()

maybecolon = Suppress(Optional(":"))
def key(s):
    return Suppress(Literal(s)) + maybecolon

to_int = pyparsing_common.integer
to_float = pyparsing_common.fnumber
mr_pattern = key("mutation_rate") + to_float
sr_pattern = key("shift_rate") + to_float
ss_pattern = key("shift_size") + to_float
fd_pattern = key("fitness_discount") +to_float
nl_pattern = key("number_of_loci") + to_int
nai_pattern = key("number_of_asex_individuals") + to_int
nit_pattern = key("number_of_iterations") + to_int
nsi_pattern = key("number_of_sex_individuals") + to_int

try:
    mr = mr_pattern.parseString(next(stdin),parseAll=True)[0]
    sr = sr_pattern.parseString(next(stdin),parseAll=True)[0]
    ss = ss_pattern.parseString(next(stdin),parseAll=True)[0]
    fd = fd_pattern.parseString(next(stdin),parseAll=True)[0]
    nl = nl_pattern.parseString(next(stdin),parseAll=True)[0]
    nai = nai_pattern.parseString(next(stdin),parseAll=True)[0]
    nsi = nsi_pattern.parseString(next(stdin),parseAll=True)[0]
    nit = nit_pattern.parseString(next(stdin))[0]

except StopIteration:
    print("Too many lines in config file",file=sys.stderr)
    sys.exit(1)

import subprocess

exec_args = ["./exec",str(nl),str(sr),str(ss),str(fd),str(nai),str(nsi),str(mr),str(nit)]
print(" ".join(exec_args),file=sys.stderr)
res = subprocess.run(exec_args,capture_output=True)


rc = res.returncode
if rc != 0:
    sys.exit(rc)
else:
    sys.stdout.write(res.stdout.decode())
sys.stderr.write(res.stderr.decode())
