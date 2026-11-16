---
name: documentation
description: "Use when setting up, building, updating, or troubleshooting schimpy Sphinx documentation: Conda runtime plus docs overlay, make html, automatic API generation, failed imports, and GitHub Pages actions."
---

# schimpy documentation

Read [contributing_docs.rst](../../../docsrc/contributing_docs.rst) for contributor
commands and publishing instructions. Run setup commands from the repository
root in a Conda-initialized terminal.
Let the user create or change environments unless execution is requested.

## Environment

- `environment.yml` supplies the runtime/test dependencies.
- `environment-doc.yml` adds documentation tools and Sphinx `>=8.2,<9`.
- Create `schimpy_doc` from the runtime YAML, then update it with the overlay
  without `--prune`. The overlay is not a standalone environment.
- `build_doc_env.bat` and `build_doc_env.sh` are convenience wrappers over those
  steps and the editable installation, not alternative dependency definitions.
- Activate `schimpy_doc` before building. Do not modify `schism_sphinx` or another
  existing environment merely to make this documentation build work.
- Install the checkout with `python -m pip install --no-deps -e .`.
  Adding `[doc]` with `--no-deps` does not install extra dependencies.
- An existing environment is not disposable. Update a compatible environment,
  or ask before replacing it.

## Build

From `docsrc`, run `make html` (`.\make.bat html` in Windows PowerShell).
For a fresh build, run `make clean` and then `make html`.
`doall.bat` is a compatibility alias, not a separate API-generation procedure.

`sphinx.ext.apidoc` generates API RST under ignored `docsrc/_api/` during every
normal build. The handwritten `schimpy.rst` includes the generated package page.
Do not edit or commit generated RST/HTML, or run a separate apidoc command into
`docsrc` that would overwrite handwritten wrappers.

The default API layout is intentional. `* - Copy.py` artifacts are excluded.
Keep intended module coverage when changing exclusions.
HTML goes to ignored `docs/html`; inspect `index.html`, affected API sections
and CLI help rather than treating successful HTML generation as completeness.

## Failures and warnings

- Missing API/CLI content: examine import warnings in the docs environment.
  Correct dependency inputs or report a source-level failure; do not add mock
  imports or suppress warnings merely to produce a green build.
- Existing mocks in `conf.py` are retained. Do not silently expand their scope.
- Existing docstring/legacy CLI warnings remain. Warning-free builds are not
  currently required; do not hide those warnings or introduce unrelated fixes.
- Ordinary HTML Mermaid rendering does not require a global npm installation.
  Do not change global npm registry or TLS settings as part of setup.

## CI and deployment

The docs workflow uses the same runtime YAML, overlay and checkout installation,
then `make html`. Do not install a released schimpy over the checkout or restore
the removed inline apidoc step.

GitHub Pages receives `docs/html` through upload/deploy actions, not a branch.
Settings -> Pages must use GitHub Actions, and the `github-pages` environment
must allow the publishing branch/tag. Pull-request checks and branch protection
are outside this change.
