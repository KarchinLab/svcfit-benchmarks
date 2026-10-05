#!/usr/bin/env python
"""
Wrapper to run VISOR SHORtS with a controlled random seed.

VISOR SHORtS hardcodes seed=0 in all wgsim.core() calls, which maps
to time(0) — no reproducibility control. This wrapper monkey-patches
pywgsim.wgsim.core to inject a deterministic seed before VISOR runs.

Each wgsim.core() call (one per region per clone) gets a unique seed
derived from the base seed: base_seed + call_count.

Usage:
    python visor_seeded.py <SEED> <VISOR SHORtS arguments...>

Example:
    python visor_seeded.py 12345 -g ref.fa -s hack_dir -b short.bed \
        -o output --coverage 50 --error 0 --indels 0
"""

import sys

if len(sys.argv) < 2:
    print("Usage: python visor_seeded.py <SEED> <VISOR SHORtS args...>")
    sys.exit(1)

BASE_SEED = int(sys.argv[1])

# Monkey-patch pywgsim before VISOR imports it
import pywgsim.wgsim as _wgsim
_original_core = _wgsim.core
_call_count = 0

def _seeded_core(*args, **kwargs):
    global _call_count
    kwargs['seed'] = BASE_SEED + _call_count
    _call_count += 1
    return _original_core(*args, **kwargs)

_wgsim.core = _seeded_core

# Now run VISOR SHORtS with remaining arguments
sys.argv = ['VISOR', 'SHORtS'] + sys.argv[2:]

from VISOR.VISOR import main as visor_main
visor_main()
