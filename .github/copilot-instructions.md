
# dms_datastore — Workspace Instructions

Follow organization standards from BayDeltaSCHISM https://raw.githubusercontent.com/CADWRDeltaModeling/BayDeltaSCHISM/refs/heads/master/AGENTS.md

Follow local project rules in AGENTS.md.

Local project rules override organization defaults.


# Build and Test

The `dms_datastore` conda environment is assumed to exist. Always activate it before running any tests or install commands.

```bash
# Install (development mode)
conda activate dms_datastore
pip install --no-deps -e .

```

pytest is configured in `pyproject.toml` (`[tool.pytest.ini_options]`): strict markers, JUnit XML output, ignores `setup.py` and `build/`.

