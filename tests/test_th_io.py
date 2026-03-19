import pytest
import pandas as pd
import os
import tempfile
from schimpy.th_io import read_th, write_th
from schimpy.model_time import multi_file_to_elapsed


def test_read_write_dated_th():
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "flux.th"
    )
    df_original = read_th(input_file, time_basis=None, to_timestamp=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        write_th(df_original, tmp_path, elapsed=False, source_file=input_file)
        df_reread = read_th(tmp_path, time_basis=None, to_timestamp=True)
        pd.testing.assert_frame_equal(df_original, df_reread, check_dtype=False)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_read_write_elapsed_th():
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "grantline_culvert_elapsed_20130115.th"
    )
    ref_time = "2013-01-15T00:00:00"
    df_original = read_th(input_file, time_basis=ref_time, to_timestamp=True, head=os.path.join(os.path.dirname(__file__), "testdata", "th_files", "grantline_culvert.th"))
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        write_th(df_original, tmp_path, elapsed=True, ref_time=ref_time, source_file=input_file)
        df_reread = read_th(tmp_path, time_basis=ref_time, to_timestamp=True)
        pd.testing.assert_frame_equal(df_original, df_reread, check_dtype=False, atol=1.0)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_comments_preserved():
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "montezuma_boat_lock.th"
    )
    original_lines = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            if (
                not line.strip()
                or line.startswith("datetime")
                or line.startswith("date")
                or line.startswith("time")
            ):
                continue
            original_lines.append(line)
    if not any("#" in line for line in original_lines):
        pytest.skip("No comments in input file to test")

    df = read_th(input_file, time_basis=None, to_timestamp=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        write_th(df, tmp_path, elapsed=False, source_file=input_file)
        output_lines = []
        with open(tmp_path, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if (
                    not line.strip()
                    or line.startswith("datetime")
                    or line.startswith("date")
                    or line.startswith("time")
                ):
                    continue
                output_lines.append(line)
        assert original_lines == output_lines
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def make_single_timestamp_with_inline_comments(path):
    text = """# top comment
 datetime a b
2000-01-01T00:00 1.00 2.00 # note one
2000-01-02T00:00 3.00 4.00
   # note for previous row
2000-01-03T00:00 5.00 6.00 # note three
"""
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def test_read_th_timestamped_preserves_inline_and_following_comments(tmp_path):
    infile = make_single_timestamp_with_inline_comments(tmp_path / "comments.th")
    df = read_th(str(infile))
    assert isinstance(df, pd.DataFrame)
    assert df.index.name == "datetime"
    assert "__comment__" in df.columns
    assert df["__comment__"].iloc[0] == "# note one"
    assert df["__comment__"].iloc[1] == "# note for previous row"
    assert df["__comment__"].iloc[2] == "# note three"


def test_roundtrip_with_format_options():
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "flux.th"
    )
    df_original = read_th(input_file, time_basis=None, to_timestamp=True)
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name
    try:
        write_th(df_original, tmp_path, elapsed=False, time_format="%Y-%m-%dT%H:%M", float_format="%.3f")
        df_reread = read_th(tmp_path, time_basis=None, to_timestamp=True)
        pd.testing.assert_frame_equal(df_original, df_reread, check_dtype=False, atol=0.001)
    finally:
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_multi_file_to_elapsed():
    input_files = [
        os.path.join(os.path.dirname(__file__), "testdata", "th_files", "flux.th"),
        os.path.join(os.path.dirname(__file__), "testdata", "th_files", "grantline_culvert.th"),
    ]
    ref_time = "2013-06-19T00:00"

    def name_transform(name):
        base = os.path.basename(name)
        return base.replace(".th", "_elapsed_test.th")

    with tempfile.TemporaryDirectory() as tmpdir:
        multi_file_to_elapsed(input_files, tmpdir, ref_time, name_transform=name_transform)
        for input_file in input_files:
            out_file = os.path.join(tmpdir, name_transform(input_file))
            assert os.path.exists(out_file)
            df_original = read_th(input_file, time_basis=None, to_timestamp=True)
            df_reread = read_th(out_file, time_basis=ref_time, to_timestamp=True, head=input_file)
            pd.testing.assert_frame_equal(df_original, df_reread, check_dtype=False, atol=1.0)
