import re
import difflib
from pathlib import Path
import pandas as pd
import pytest
import datetime as dt

import schimpy.param as param_mod  # <-- this imports YOUR /mnt/data/param.py module


from schimpy import param  # uses your Params and read_params APIs
START_KEYS = {"start_year": 2025, "start_month": 12, "start_day": 31, "start_hour": 12}
KEY_PATTERN = re.compile(r"^\s*(start_year|start_month|start_day|start_hour)\s*=")

@pytest.fixture
def datadir():
    # adjust if your layout differs
    return Path(__file__).parent.resolve() / "testdata"

def _changed_body_lines(a: str, b: str):
    """Return only the '-' and '+' body lines from a unified diff (skip headers)."""
    diff = list(difflib.unified_diff(a.splitlines(), b.splitlines(), lineterm=""))
    body = []
    for line in diff:
        if line.startswith(('---', '+++', '@@')):
            continue
        if line.startswith(('-', '+')):
            body.append(line)
    return body

def test_param_writer_changes_only_start_fields(datadir):
    src = datadir / "param.nml"
    out = datadir / "param.nml.test"

    original_text = src.read_text(encoding="utf-8")

    # Prefer constructing directly; your current read_params passes file contents
    # to Params(...), which depends on nml.parse accepting raw text.
    p = param.Params(str(src))   # exercises the same API surface
    p.run_start = pd.Timestamp(2025, 12, 31, 12)
    p.write(out)

    new_text = out.read_text(encoding="utf-8")
    changes = _changed_body_lines(original_text, new_text)

    # Keep only the start_* lines on both sides
    minus = [ln[1:] for ln in changes if ln.startswith('-') and KEY_PATTERN.match(ln[1:])]
    plus  = [ln[1:] for ln in changes if ln.startswith('+') and KEY_PATTERN.match(ln[1:])]

    # There should be one '-' and one '+' for each of the four keys
    assert len(minus) == 4 and len(plus) == 4, f"Unexpected changes: {changes}"

    # Parse the '+' lines into a dict and check values
    new_vals = {}
    for ln in plus:
        m = re.match(r"^\s*(\w+)\s*=\s*([0-9]+)\s*(?:!.*)?$", ln)
        assert m, f"Unexpected assignment format: {ln}"
        k, v = m.group(1), int(m.group(2))
        new_vals[k] = v

    assert set(new_vals.keys()) == set(START_KEYS.keys())
    for k, v in START_KEYS.items():
        assert new_vals[k] == v, f"{k} expected {v}, got {new_vals[k]}"

    # Sanity: property round-trip
    assert p.run_start == pd.Timestamp(2025, 12, 31, 12)

DATA_DIR = Path(__file__).parent / "testdata"

def _load_params():
    src = DATA_DIR / "param.nml"
    assert src.exists(), f"Missing fixture: {src}"
    return param.Params(str(src))

def test_run_start_setter_updates_fields():
    p = _load_params()

    # Set via string; minutes/seconds ignored, hour precision only
    p.run_start = "2025-12-31T12:00"

    assert p["start_year"]  == 2025
    assert p["start_month"] == 12   # <-- fix: December, not January
    assert p["start_day"]   == 31
    assert p["start_hour"]  == 12

    # Getter round-trips at hour precision
    assert p.run_start == pd.Timestamp(2025, 12, 31, 12)

def test_run_start_accepts_various_input_types():
    p = _load_params()

    cases = [
        "2025-12-31 12:00",
        dt.datetime(2025, 12, 31, 12, 30),        # minutes ignored
        pd.Timestamp(2025, 12, 31, 12, 59, 59),   # seconds ignored
    ]
    for val in cases:
        p.run_start = val
        assert (p["start_year"], p["start_month"], p["start_day"], p["start_hour"]) == (2025, 12, 31, 12)
        assert p.run_start == pd.Timestamp(2025, 12, 31, 12)   