from pathlib import Path

import pandas as pd
from click.testing import CliRunner

import schimpy.model_time as model_time


# ---------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------

def write_text(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8", newline="\n")
    return path


def read_nonempty_lines(path: Path):
    return [
        line.rstrip("\n")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def make_single_timestamp_th(path: Path) -> Path:
    text = """# Originally entitled: sample.th
# Top comment line 2
datetime a b
2000-01-01T00:00 1.00 2.00
2000-01-02T00:00 3.00 4.00
2000-01-03T00:00 5.00 6.00
"""
    return write_text(path, text)


def make_single_timestamp_with_inline_comments(path: Path) -> Path:
    text = """# top comment
datetime a b
2000-01-01T00:00 1.00 2.00 # note one
2000-01-02T00:00 3.00 4.00
   # note for previous row
2000-01-03T00:00 5.00 6.00 # note three
"""
    return write_text(path, text)


def make_elapsed_th_with_header(path: Path) -> Path:
    text = """# elapsed file
time a b
0.0 10.0 20.0
86400.0 11.0 21.0
172800.0 12.0 22.0
"""
    return write_text(path, text)


def make_plain_elapsed_th(path: Path) -> Path:
    text = """0.0 10.0 20.0
86400.0 11.0 21.0
172800.0 12.0 22.0
"""
    return write_text(path, text)


def make_multiindex_timestamp_th(path: Path) -> Path:
    idx = pd.Index(
        pd.to_datetime(["2000-01-01 00:00", "2000-01-02 00:00"]),
        name="datetime",
    )
    cols = pd.MultiIndex.from_tuples(
        [
            ("temperature", "delta_src_1"),
            ("temperature", "delta_src_3"),
            ("salinity", "delta_src_1"),
            ("salinity", "delta_src_3"),
        ],
        names=["variable", "location"],
    )
    df = pd.DataFrame(
        [
            [9.50, -9999.0, 0.67, 0.80],
            [9.60, -9999.0, 0.68, 0.81],
        ],
        index=idx,
        columns=cols,
    )
    df.to_csv(
        path,
        sep=" ",
        date_format="%Y-%m-%dT%H:%M",
        float_format="%.2f",
    )
    return path


# ---------------------------------------------------------------------
# CLI naming / dispatch
# ---------------------------------------------------------------------

def test_cli_help_uses_underscore_subcommands():
    runner = CliRunner()
    result = runner.invoke(model_time.model_time_cli, ["--help"])
    assert result.exit_code == 0, result.output
    assert "to_elapsed" in result.output
    assert "to_date" in result.output
    assert "to-elapsed" not in result.output
    assert "to-date" not in result.output


# ---------------------------------------------------------------------
# file_to_elapsed: single-header timestamped files
# ---------------------------------------------------------------------

def test_file_to_elapsed_single_header_basic(tmp_path):
    infile = make_single_timestamp_th(tmp_path / "single.th")
    outfile = tmp_path / "elapsed.th"

    model_time.file_to_elapsed(
        infile=str(infile),
        start="2000-01-01T00:00",
        outpath=str(outfile),
    )

    lines = read_nonempty_lines(outfile)
    assert len(lines) == 3
    assert lines[0].split() == ["0.0", "1.00", "2.00"]
    assert lines[1].split() == ["86400.0", "3.00", "4.00"]
    assert lines[2].split() == ["172800.0", "5.00", "6.00"]


def test_cli_to_elapsed_single_header(tmp_path):
    infile = make_single_timestamp_th(tmp_path / "single.th")
    outfile = tmp_path / "elapsed.th"

    runner = CliRunner()
    result = runner.invoke(
        model_time.model_time_cli,
        [
            "to_elapsed",
            str(infile),
            "--start",
            "2000-01-01T00:00",
            "--out",
            str(outfile),
        ],
    )
    assert result.exit_code == 0, result.output

    lines = read_nonempty_lines(outfile)
    assert len(lines) == 3
    assert lines[0].split() == ["0.0", "1.00", "2.00"]
    assert lines[1].split() == ["86400.0", "3.00", "4.00"]
    assert lines[2].split() == ["172800.0", "5.00", "6.00"]


def test_file_to_elapsed_start_after_file_emits_last_prior_row(tmp_path):
    infile = make_single_timestamp_th(tmp_path / "single.th")
    outfile = tmp_path / "elapsed_late.th"

    model_time.file_to_elapsed(
        infile=str(infile),
        start="2000-01-10T00:00",
        outpath=str(outfile),
    )

    lines = read_nonempty_lines(outfile)
    assert len(lines) == 1
    assert lines[0].split() == ["-604800.0", "5.00", "6.00"]


# ---------------------------------------------------------------------
# pandas MultiIndex timestamped files
# ---------------------------------------------------------------------

def test_pandas_multiindex_fixture_is_loadable(tmp_path):
    infile = make_multiindex_timestamp_th(tmp_path / "multi.th")

    df = pd.read_csv(infile, sep=r"\s+", header=[0, 1], index_col=0)
    assert df.shape == (2, 4)
    assert tuple(df.columns[0]) == ("temperature", "delta_src_1")
    assert tuple(df.columns[-1]) == ("salinity", "delta_src_3")


def test_file_to_elapsed_multiindex_timestamp_file(tmp_path):
    infile = make_multiindex_timestamp_th(tmp_path / "multi.th")
    outfile = tmp_path / "multi_elapsed.th"

    model_time.file_to_elapsed(
        infile=str(infile),
        start="2000-01-01T00:00",
        outpath=str(outfile),
    )

    lines = read_nonempty_lines(outfile)
    assert len(lines) == 2
    assert lines[0].split()[0] == "0.0"
    assert lines[1].split()[0] == "86400.0"
    assert "9.50" in lines[0]
    assert "-9999.00" in lines[0]
    assert "0.67" in lines[0]


def test_cli_to_elapsed_multiindex_timestamp_file(tmp_path):
    infile = make_multiindex_timestamp_th(tmp_path / "multi.th")
    outfile = tmp_path / "multi_elapsed.th"

    runner = CliRunner()
    result = runner.invoke(
        model_time.model_time_cli,
        [
            "to_elapsed",
            str(infile),
            "--start",
            "2000-01-01T00:00",
            "--out",
            str(outfile),
        ],
    )

    assert result.exit_code == 0, result.output
    lines = read_nonempty_lines(outfile)
    assert len(lines) == 2
    assert lines[0].split()[0] == "0.0"
    assert lines[1].split()[0] == "86400.0"


# ---------------------------------------------------------------------
# clip / fast-forward on timestamped input
# ---------------------------------------------------------------------

def test_clip_timestamped_single_header_basic(tmp_path):
    infile = make_single_timestamp_th(tmp_path / "single.th")
    outfile = tmp_path / "clipped.th"

    runner = CliRunner()
    result = runner.invoke(
        model_time.model_time_cli,
        [
            "clip",
            str(infile),
            "--start",
            "2000-01-01T00:00",
            "--clip_start",
            "2000-01-02T00:00",
            "--out",
            str(outfile),
        ],
    )

    assert result.exit_code == 0, result.output
    lines = read_nonempty_lines(outfile)
    assert len(lines) == 2
    assert lines[0].split() == ["0.0", "3.00", "4.00"]
    assert lines[1].split() == ["86400.0", "5.00", "6.00"]


def test_clip_timestamped_multiindex_file_after_patch(tmp_path):
    infile = make_multiindex_timestamp_th(tmp_path / "multi.th")
    outfile = tmp_path / "multi_clipped.th"

    runner = CliRunner()
    result = runner.invoke(
        model_time.model_time_cli,
        [
            "clip",
            str(infile),
            "--start",
            "2000-01-01T00:00",
            "--clip_start",
            "2000-01-02T00:00",
            "--out",
            str(outfile),
        ],
    )

    assert result.exit_code == 0, result.output
    lines = read_nonempty_lines(outfile)
    assert len(lines) == 1
    assert lines[0].split()[0] == "0.0"
    assert "9.60" in lines[0]
    assert "0.68" in lines[0]


# ---------------------------------------------------------------------
# to_date round-trip smoke test
# ---------------------------------------------------------------------

def test_elapsed_to_date_cli_basic(tmp_path):
    infile = make_plain_elapsed_th(tmp_path / "elapsed_plain.th")
    outfile = tmp_path / "dated.th"

    runner = CliRunner()
    result = runner.invoke(
        model_time.model_time_cli,
        [
            "to_date",
            str(infile),
            "--start",
            "2000-01-01T00:00",
            "--out",
            str(outfile),
        ],
    )

    assert result.exit_code == 0, result.output
    text = outfile.read_text(encoding="utf-8")
    assert "2000-01-01T00:00" in text
    assert "2000-01-02T00:00" in text
    assert "2000-01-03T00:00" in text