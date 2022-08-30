#!/usr/bin/env python3

#
# Parses lines in input file
#

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
sr_pattern = key("shift_rate") + to_float
ss_pattern = key("shift_size") + to_float
fd_pattern = key("fitness_discount") +to_float
nl_pattern = key("number_of_loci") + to_int
na_pattern = key("number_of_alleles") + to_int
mr_pattern = (Literal("mutation_rate") + Optional(":")).setParseAction(lambda s,loc,toks: s[len(toks[0]):])
mc_pattern = (Literal("mutation_contrib") + Optional(":")).setParseAction(lambda s,loc,toks: s[len(toks[0]):])
trts_pattern =(Literal("traits") + Optional(":")).setParseAction(lambda s,loc,toks: s[len(toks[0]):])
smr_pattern = key("sex_mutation_rate") + to_float
nai_pattern = key("number_of_asex_individuals") + to_int
nit_pattern = key("number_of_iterations") + to_int
nsi_pattern = key("number_of_sex_individuals") + to_int

try:
    sr = sr_pattern.parseString(next(stdin),parseAll=True)[0]
    ss = ss_pattern.parseString(next(stdin),parseAll=True)[0]
    fd = fd_pattern.parseString(next(stdin),parseAll=True)[0]
    nl = nl_pattern.parseString(next(stdin),parseAll=True)[0]
    na = na_pattern.parseString(next(stdin),parseAll=True)[0]
    mr = mr_pattern.parseString(next(stdin))[0].strip()
    mc = mc_pattern.parseString(next(stdin))[0].strip()
    trts = trts_pattern.parseString(next(stdin))[0].strip()
    smr = smr_pattern.parseString(next(stdin),parseAll=True)[0]
    nai = nai_pattern.parseString(next(stdin),parseAll=True)[0]
    nsi = nsi_pattern.parseString(next(stdin),parseAll=True)[0]
    nit = nit_pattern.parseString(next(stdin))[0]

except StopIteration:
    print("Couldn't parse input",file=sys.stderr)
    sys.exit(1)

m = '0'
if "-m" in sys.argv[1:]:
    m = '1'

import subprocess

exec_args = ["./exec",str(nl),
                      str(na),
                      str(sr),
                      str(ss),
                      str(fd),
                      str(nai),
                      str(nsi),
                      str(mr),
                      str(mc),
                      str(trts),
                      str(smr),
                      m,
                      str(nit)]

print(" ".join(exec_args),file=sys.stderr)
res = subprocess.run(exec_args,stdout = subprocess.PIPE,
        stderr = subprocess.PIPE)

sys.stdout.write(res.stdout.decode())

err_msg = res.stderr.decode()
if err_msg:
    print("Errors:",file=sys.stderr)
    sys.stderr.write("  "+err_msg)
