import pytest
import pandas as pd
import os
import tempfile
from schimpy.model_time import read_th, write_th, multi_file_to_elapsed


def test_read_write_dated_th():
    """Test that reading and writing a dated .th file preserves the data"""
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "flux.th"
    )

    # Read the file - assuming it's a dated file
    df_original = read_th(input_file, time_basis=None, to_timestamp=True)

    # Write to a temporary file
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Write the dataframe
        write_th(df_original, tmp_path, elapsed=False, source_file=input_file)

        # Read it back
        df_reread = read_th(tmp_path, time_basis=None, to_timestamp=True)

        # Compare dataframes
        pd.testing.assert_frame_equal(df_original, df_reread, check_dtype=False)

    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_read_write_elapsed_th():
    """Test that reading and writing an elapsed .th file preserves the data"""
    input_file = os.path.join(
        os.path.dirname(__file__),
        "testdata",
        "th_files",
        "grantline_culvert_elapsed_20130115.th",
    )
    ref_time = "2013-01-15T00:00:00"
    header_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "grantline_culvert.th"
    )

    # Read the file
    df_original = read_th(
        input_file, time_basis=ref_time, to_timestamp=True, head=header_file
    )

    # Write to elapsed format
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        # Write as elapsed
        write_th(
            df_original,
            tmp_path,
            elapsed=True,
            ref_time=ref_time,
            source_file=input_file,
        )

        # Read it back
        df_reread = read_th(tmp_path, time_basis=ref_time, to_timestamp=True)

        # Compare dataframes (allow for small floating point differences)
        pd.testing.assert_frame_equal(
            df_original,
            df_reread,
            check_dtype=False,
            atol=1.0,  # Allow 1 second difference due to rounding
        )

    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_comments_preserved():
    """Test that comments are preserved when source_file is specified"""
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "montezuma_boat_lock.th"
    )

    # Read original file to extract inline comments
    original_lines = []
    with open(input_file, "r") as f:
        for line in f:
            line = line.rstrip("\n")
            # Skip header and empty lines
            if (
                not line.strip()
                or line.startswith("datetime")
                or line.startswith("date")
                or line.startswith("time")
            ):
                continue
            original_lines.append(line)

    # If there are no lines with comments, skip this test
    if not any("#" in line for line in original_lines):
        pytest.skip("No comments in input file to test")

    # Read the file
    df = read_th(input_file, time_basis=None, to_timestamp=True)

    # Write with source_file to preserve comments
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        write_th(df, tmp_path, elapsed=False, source_file=input_file)

        # Read lines from output file
        output_lines = []
        with open(tmp_path, "r") as f:
            for line in f:
                line = line.rstrip("\n")
                # Skip header and empty lines
                if (
                    not line.strip()
                    or line.startswith("datetime")
                    or line.startswith("date")
                    or line.startswith("time")
                ):
                    continue
                output_lines.append(line)

        # Compare lines (including comments)
        assert original_lines == output_lines, "Comments were not preserved"

    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_roundtrip_with_format_options():
    """Test roundtrip with custom format options"""
    input_file = os.path.join(
        os.path.dirname(__file__), "testdata", "th_files", "flux.th"
    )

    # Read the file
    df_original = read_th(input_file, time_basis=None, to_timestamp=True)

    # Write with custom formats
    with tempfile.NamedTemporaryFile(mode="w", suffix=".th", delete=False) as tmp:
        tmp_path = tmp.name

    try:
        write_th(
            df_original,
            tmp_path,
            elapsed=False,
            time_format="%Y-%m-%dT%H:%M",
            float_format="%.3f",
        )

        # Read it back
        df_reread = read_th(tmp_path, time_basis=None, to_timestamp=True)

        # Compare dataframes (with tolerance for format changes)
        pd.testing.assert_frame_equal(
            df_original, df_reread, check_dtype=False, atol=0.001
        )

    finally:
        # Clean up
        if os.path.exists(tmp_path):
            os.remove(tmp_path)


def test_multi_file_to_elapsed():
    """Test multi_file_to_elapsed produces an elapsed file and preserves data."""
    input_files = [
        os.path.join(os.path.dirname(__file__), "testdata", "th_files", "flux.th"),
        os.path.join(os.path.dirname(__file__), "testdata", "th_files", "grantline_culvert.th"),
    ]

    ref_time = "2013-06-19T00:00"

    def name_transform(name):
        base = os.path.basename(name)
        return base.replace(".th", "_elapsed_test.th")

    with tempfile.TemporaryDirectory() as tmpdir:
        # Convert to elapsed using multi_file_to_elapsed
        multi_file_to_elapsed(
            input_files, tmpdir, ref_time, name_transform=name_transform
        )

        for input_file in input_files:
            out_file = os.path.join(tmpdir, name_transform(input_file))
            assert os.path.exists(out_file), "Elapsed output file was not created"

            # Read original dated file
            df_original = read_th(input_file, time_basis=None, to_timestamp=True)

            # Read back and compare
            df_reread = read_th(out_file, time_basis=ref_time, to_timestamp=True, head=input_file)
            pd.testing.assert_frame_equal(
                df_original, df_reread, check_dtype=False, atol=1.0
            )
