#!/usr/bin/env python3
# -*- coding: utf8 -*-
"""
Brouwer Bounds Loader module.
Unifies bounds parsing, theoretical limits, and blockcode parameters
into a clean Singleton manager: BoundsManager.
"""
import os
import re
import json
import itertools
from cwc.utils.logging import log1, log2, log3, log4, log5, log6, is_debug, set_debug
from cwc.utils.paths import DATA_DIR
from cwc.utils.math import nextPrimePower
from core import binom

set_debug(3)

def load_hardcoded_strings():
    bounds_dir = os.path.join(DATA_DIR, 'bounds')
    main_path = os.path.join(bounds_dir, 'brouwer_main.txt')
    more_path = os.path.join(bounds_dir, 'brouwer_more.txt')
    block_path = os.path.join(bounds_dir, 'brouwer_block.txt')
        
    with open(main_path, 'r', encoding='utf-8') as f:
        boundstr = f.read()
        
    with open(more_path, 'r', encoding='utf-8') as f:
        more_bound = json.load(f)
        
    with open(block_path, 'r', encoding='utf-8') as f:
        block_bound_str = json.load(f)
        
    return boundstr, more_bound, block_bound_str

lost = [(23,6,10,2969),(23,6,11,3535),(24,6,8,1848),(24,6,10,4174),(24,6,12,5558), (25,6,8,2541),(26,6,8,3460),(27,6,8,4715),(28,6,8,6248)]

def parse_bounds_cell(cell_str):
    if not cell_str.strip():
        return None
    parts = cell_str.split('-')
    
    def get_leading_int(s):
        m = re.match(r'^\d+', s.strip())
        return int(m.group(0)) if m else None

    if len(parts) == 1:
        val = get_leading_int(parts[0])
        if val is not None:
            return (val, val)
    elif len(parts) == 2:
        val1 = get_leading_int(parts[0])
        val2 = get_leading_int(parts[1])
        if val1 is not None or val2 is not None:
            return (val1, val2)
            
    return None


class BoundsManager:
    """Singleton Manager for Constant Weight and Block Code Bounds."""
    _instance = None

    def __new__(cls, *args, **kwargs):
        if not cls._instance:
            cls._instance = super(BoundsManager, cls).__new__(cls)
            cls._instance._initialized = False
        return cls._instance

    def __init__(self):
        if getattr(self, "_initialized", False):
            return
        self._initialized = True
        self.bounds = {}
        self.block_bounds = {}
        self.original_block_bounds = {}
        self._load_all_bounds()

    def setBound(self, m, D, w, lb, ub, desc):
        if D == 2:
            return
        self.bounds[(m, D, w)] = (lb, ub, desc)

    def setBlockBound(self, w, q, d, lb, ub):
        if (w, q, d) in self.block_bounds:
            old = self.block_bounds[(w, q, d)]
            if lb is None and ub is None:
                self.block_bounds[(w, q, d)] = (old[1], old[1])
            elif ub is None:
                self.block_bounds[(w, q, d)] = (lb, old[1])
        else:
            self.block_bounds[(w, q, d)] = (lb, ub)

    def _load_all_bounds(self):
        boundstr, more_bound, block_bound_str = load_hardcoded_strings()
        
        # 1. Parse main bounds table
        for codedesc in boundstr.split("Bounds on "):
            D = 0
            for liness in codedesc.split('\n'):
                if '\u03c4' in liness or 'τ' in liness:
                    continue
                lines = liness.encode('ascii', 'ignore').decode('ascii')
                if lines.startswith('A(n,'):
                    D = [int(s) for s in lines.split(',') if s.isdigit()][0]
                elif lines.startswith('n\\w'):
                    numb = [int(s) for s in lines.split() if s.isdigit()]
                    w = numb[0]
                else:
                    first = True
                    m = D + 2
                    i = 0
                    for mmm in lines.split('\t'):
                        mm = re.findall(r'\d+', mmm)
                        if len(mm) > 0:
                            if first:
                                m = int(mm[0])
                            else:
                                if len(mm) <= 1:
                                    self.setBound(m, D, w + i, int(mm[0]), int(mm[0]), mmm)
                                if len(mm) >= 2:
                                    self.setBound(m, D, w + i, int(mm[0]), int(mm[1]), mmm)
                                i += 1
                        first = False

        # 2. Parse individual/extension tables
        for D, w, codedesc in more_bound:
            ndesc = True
            mm = []
            for liness in codedesc.split('\n'):
                if '\u03c4' in liness or 'τ' in liness:
                    ndesc = True
                    continue
                lines = liness.encode('ascii', 'ignore').decode('ascii')
                if len(lines.strip()) == 0:
                    continue
                if ndesc:
                    mm = [int(s) for s in lines.split('\t') if s.isdigit()]
                    ndesc = False
                else:
                    ndesc = True
                    i = 0
                    for nnn in lines.split('\t'):
                        if i < len(mm):
                            res = parse_bounds_cell(nnn)
                            if res:
                                lb, ub = res
                                self.setBound(mm[i], D, w, lb, ub, nnn)
                        i += 1

        # 3. Patchwork adjustments
        # This is just a lazy patch for a parsing exception where specific ranges of A(m, 14, 8) are missing/malformed in raw tables
        print("Applying quickfix patch for A(m, 14, 8) bounds...")
        for m_quickfix in [65, 66, 67, 68, 69, 70]:
            x = 72 + (m_quickfix - 64)
            self.setBound(m_quickfix, 14, 8, 72, x, str(x))
        for m_quickfix in [73, 74, 75, 76]:
            x = 89 + (m_quickfix - 72)
            self.setBound(m_quickfix, 14, 8, 89, x, str(x))

        # 4. Ingest lost bounds
        for (m, D, w, n) in lost:
            if (m, D, w) in self.bounds:
                lb, ub, desc = self.bounds[(m, D, w)]
                self.bounds[(m, D, w)] = (n, ub, 'lost:' + desc)

        # 5. Parse block code tables
        for q, codedesc in block_bound_str:
            ndesc = True
            mm = []
            for liness in codedesc.split('\n'):
                lines = liness.encode('ascii', 'ignore').decode('ascii')
                if len(lines.strip()) == 0:
                    continue
                if ndesc:
                    mm = [int(s) for s in lines.split('\t') if s.isdigit()]
                    ndesc = False
                else:
                    i = 0
                    nnn = lines.split('\t')
                    w = nnn[0]
                    for nn in nnn[1:]:
                        if nnn[0].endswith('lb'):
                            if nn == "=":
                                self.setBlockBound(int(w[:-2]), q, mm[i], None, None)
                            else:
                                self.setBlockBound(int(w[:-2]), q, mm[i], int(nn), None)
                        elif nnn[0].endswith('ub'):
                            self.setBlockBound(int(w[:-2]), q, mm[i], None, int(nn))
                        else:
                            nnnn = re.findall(r'\d+', nn)
                            if len(nn) > 0:
                                if len(nnnn) == 1:
                                    self.setBlockBound(int(w), q, mm[i], int(nn), int(nn))
                                elif len(nnnn) == 2:
                                    self.setBlockBound(int(w), q, mm[i], int(nnnn[0]), int(nnnn[1]))
                        i += 1

        # 5.5. Save original block bounds
        self.original_block_bounds = dict(self.block_bounds)

        # 6. Ingest linear block codes from linear_quary_codes.json to update lower bounds
        try:
            from cwc.utils.paths import LINEAR_QUARY_CODES_JSON
            if LINEAR_QUARY_CODES_JSON.exists():
                with open(LINEAR_QUARY_CODES_JSON, "r") as f:
                    linear_blockdata = json.load(f)
                for q_str, by_w in linear_blockdata.items():
                    q = int(q_str)
                    for w_str, by_l in by_w.items():
                        w = int(w_str)
                        for l_str, val in by_l.items():
                            l = int(l_str)
                            if isinstance(val, (list, tuple)) and len(val) >= 1:
                                d = val[0]
                                size = pow(q, l)
                                if d > 0:
                                    old_lb, old_ub = self.block_bounds.get((w, q, d), (None, None))
                                    new_lb = size
                                    if old_lb is not None:
                                        new_lb = max(old_lb, size)
                                    self.block_bounds[(w, q, d)] = (new_lb, old_ub)
        except Exception as e:
            pass


# --- Compatibility Layer ---

def loadBounds():
    mgr = BoundsManager()
    return mgr.bounds, mgr.block_bounds

def boundA(m, D, w):
    mgr = BoundsManager()
    if D == 2:
        n = binom(m, w)
        return (n, n)
    if (m, D, w) in mgr.bounds:
        return mgr.bounds[(m, D, w)][0], mgr.bounds[(m, D, w)][1]
    return None, None

def boundADesc(m, D, w):
    mgr = BoundsManager()
    if (m, D, w) in mgr.bounds:
        return mgr.bounds[(m, D, w)][2]
    return None

def boundblockA(w, q, d):
    mgr = BoundsManager()
    if q == 2 and d % 2 == 1:
        return boundblockA(w + 1, q, d + 1)
    if (w, q, d) in mgr.block_bounds:
        return mgr.block_bounds[(w, q, d)]
    elif w <= q:
        if nextPrimePower(q) == q:
            rs = pow(q, w - d + 1)
            return (rs, rs)
    return None, None

def boundblockAOriginal(w, q, d):
    mgr = BoundsManager()
    if q == 2 and d % 2 == 1:
        return boundblockAOriginal(w + 1, q, d + 1)
    if (w, q, d) in mgr.original_block_bounds:
        return mgr.original_block_bounds[(w, q, d)]
    elif w <= q:
        if nextPrimePower(q) == q:
            rs = pow(q, w - d + 1)
            return (rs, rs)
    return None, None

def __getattr__(name):
    if name == 'bounds':
        return BoundsManager().bounds
    if name == 'block_bounds':
        return BoundsManager().block_bounds
    if name == 'original_block_bounds':
        return BoundsManager().original_block_bounds
    raise AttributeError(f"module {__name__} has no attribute {name}")
