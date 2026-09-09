
# schimpy agent guidance

`schimpy` provides reusable SCHISM preprocessing functionality intended to be useful beyond BayDeltaSCHISM.

Keep Bay-Delta-specific policy, stations, and run assumptions in downstream packages such as `bdschism`.

## Package role

`schimpy` owns generally reusable SCHISM preprocessing, file/configuration support, and related utilities.

It may depend on `vtools`.

Its dependency on `dms_datastore` is transitional. Do not introduce new `dms_datastore` dependencies.

Do not push BayDeltaSCHISM-specific workflow policy into `schimpy` unless the behavior is genuinely reusable for SCHISM users generally.

## Dependency policy

- Assume the established scientific Python/tool stack is available.
- Do not program indirect dependencies on HEC-DSS.
- Do not add to the spatial dependency stack without a specific need.
- Do not add to the plotting dependency stack without a specific need.

## Coding practice

- Plan before coding.
- Keep functions single-purpose.
- Keep functions testable.
- Do not refactor outside the scope of the requested work. Alert the user if broader refactoring appears warranted.
- Do not contract existing documentation.
- Preserve NumPy-style documentation; repair it when interfaces change.
- Prefer explicit errors such as `ValueError` over elaborate recovery from invalid arguments.
- Avoid making inference from surrounding files the only way to use an API. Inference may be a convenience, but important inputs should also be supplyable explicitly.
- Assume the established tool stack exists; do not add elaborate defensive discovery for expected dependencies or executables.
- Prefer reusable behavior over application-specific assumptions.

## Configuration

Use `schimpy.schism_yaml` for SCHISM preprocessor configuration.

Its include behavior is an established part of the preprocessing workflow and should be reused rather than replaced by ad hoc YAML handling.

Do not introduce `schimpy.schism_yaml` into packages that do not need SCHISM-specific configuration.

OmegaConf may be used elsewhere for workflow description, but it does not replace the role of `schimpy.schism_yaml` in the SCHISM preprocessor.

## SCHISM wrappers and utilities

SCHISM wrapper functionality should use established configuration to avoid hard-wiring executable or version identities.

Prefer stable utility identities such as `combine_hotstart` rather than version-specific executable names.

Wrappers may infer simulation context from standard SCHISM files such as `param.nml`, `hgrid.gr3`, or model outputs, but important context must also be explicitly supplyable.

Assume established SCHISM utilities are available on the execution path rather than adding elaborate fallback discovery.

Link creation should follow established configuration.

## Command-line interfaces

For commands that operate on data:

- provide a workhorse function where practical
- let the CLI or wrapper validate arguments and open external inputs
- provide informative `-h` and `--help`
- use `@click.help_option("-h", "--help")`
- configure logging at the CLI entry point when logging is used

If a command participates in a higher-level BayDeltaSCHISM or `bdschism` hierarchy, ensure the owning downstream package handles its own registration.

## Logging

Use the established logging configuration patterns when present.

Configure logging at the entry point rather than inside reusable workhorse functions.

Do not invent command-specific logging frameworks when shared machinery already exists.

## Time-series conventions

Where `schimpy` uses `vtools`, use established `vtools` functionality rather than recreating common time-series operations.

Use lower-case frequency strings such as:

`min`, `h`, `d`, `s`

Do not embed Bay-Delta repository policy in reusable SCHISM preprocessing code.

## Documentation

schimpy documentation should focus on CLI, functions and a manual for the preprocessor.

## Testing

- Use `pytest`.
- Use the SCHISM environment for formal or informal testing.
- Mark network-dependent tests as `integration`.
- Normal CI should exclude integration tests where appropriate.


