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

maybecolon = Optional(":")

to_int = pyparsing_common.integer

to_float = pyparsing_common.fnumber

mr_pattern = Literal("mutation_rate") + maybecolon + to_float

try:
    mr = mr_pattern.parseString(next(stdin),parseAll=True)[-1]
    
    assert 0 <= mr <=1
    
    tp_pattern = Literal("shift_rate") + maybecolon + to_float
    
    tp = tp_pattern.parseString(next(stdin),parseAll=True)[-1]
    
    assert tp <= 1.
    
    fd_pattern = Literal("fitness_discount") + maybecolon +to_float
    
    fd = fd_pattern.parseString(next(stdin),parseAll=True)[-1]
    
    nl_pattern = Literal("number_of_loci") + maybecolon + to_int
    
    nl = nl_pattern.parseString(next(stdin),parseAll=True)[-1]
            
    nai_pattern = Literal("number_of_asex_individuals") + maybecolon + to_int

    nai = nai_pattern.parseString(next(stdin),parseAll=True)[-1]

    nsi_pattern = Literal("number_of_sex_individuals") + maybecolon + to_int

    nsi = nsi_pattern.parseString(next(stdin),parseAll=True)[-1]

    nit_pattern = Literal("number_of_iterations") + to_int

    nit = nit_pattern.parseString(next(stdin))[-1]

except StopIteration:
    print("Too many lines in config file",file=sys.stderr)
    sys.exit(1)

import subprocess

res = subprocess.run(
        ["./exec",str(nl),str(tp),str(fd),str(nai),str(nsi),str(nit)],
        capture_output=True)

rc = res.returncode
if rc != 0:
    sys.stderr.write(res.stderr.decode())
    sys.exit(rc)
else:
    sys.stdout.write(res.stdout.decode())
