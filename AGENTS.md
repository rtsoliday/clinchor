# Agent Notes for clinchor

This file contains a short tactical summary based on repository evidence. `../llm-wiki/scripts/refresh_wiki.py` rewrites only the machine-managed block.

<!-- BEGIN MACHINE:summary -->
## Quick start
- Repository-local guidance is sufficient: start with `AGENTS.md`, `README.md`, `docs/`, build/test/config files, and the source tree.
- **clinchor** is a specialized software tool designed for the calculation of growth rates of multi-bunch beam modes in particle accelerators. It can be used to analyze and predict beam instabilities—a critical aspect of designing and operating high-performance accelerator systems.
- Primary work areas: `doc`, `examples`, `rpm`, `src`.

## Read first
- `README.md`: Primary project overview and workflow notes
- `examples/MonteCarlo/README`: Supporting documentation with operational details
- `examples/sweep/README`: Supporting documentation with operational details
- `doc/Makefile`: Build system entry point or dependency manifest
- `src/Makefile`: Build system entry point or dependency manifest

## Build and test
- Documented setup/build commands: `make`, `make a file name with a different memory if one wants to expand`.
- Detected build systems: GNU Make.
- Unknown: no test workflow evidence was found in the inspected files.
- Likely run commands or operator entry points: `./clinchor`.

## Operational warnings
- Local checkout layout appears significant; avoid casual changes to sibling-repo assumptions or relative paths.
- Platform-specific dependency setup matters; do not assume one platform's build recipe carries over unchanged.

## Compatibility constraints
- Cross-platform support exists, but platform-specific dependency setup matters.
- Build and runtime behavior likely depends on neighboring core toolkit checkouts.

## Related knowledge
- Repository-local documentation should be treated as authoritative.
- If a shared `llm-wiki/` directory is present in this workspace or parent folder, consult [the matching repo page](../llm-wiki/repos/clinchor.md) for additional architectural context.
- If no shared wiki is present, continue using repository-local evidence only.
- If available, [the SDDS concept page](../llm-wiki/concepts/sdds.md) adds broader cross-repo context.
- If present in this workspace, [the cross-repo map](../llm-wiki/insights/cross-repo-map.md) helps explain related repositories.
<!-- END MACHINE:summary -->

## Human notes
Add durable repo-specific instructions here.
