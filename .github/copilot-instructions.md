## Quick context

- Language: R (scripts live under `src/`). Core analysis is script-driven, not a packaged R library.
- Data: raw TCGA downloads are under `GDCdata/TCGA-*/Transcriptome_Profiling/Gene_Expression_Quantification/`.
- Inputs/Outputs: `data/raw/`, `data/processed/`, `output/`, and top-level RDS files (`df.rds`, `results.rds`).

## Primary entry points

- `run.r` — intended pipeline entry point (currently empty in the repo). If adding runnable pipeline code, prefer adding an Rscript-compatible wrapper here that sources `src/*.r` and writes outputs to `output/`.
- `src/rnaseq.r` — RNA-seq preprocessing / helpers (currently empty but expected place for read/normalize functions).
- `src/hmm.r` — contains an R6 `ProteinDomainHMM` class. This is the most complete example of project style and idioms.

## What an AI agent should know to be productive

- Project is script-first. Modify or add `src/*.r` functions and expose a small driver in `run.r` or light CLI wrappers (use `Rscript run.r` or `R -e \"source('run.r')\"`).
- Data is large and stored under `GDCdata/` — avoid attempting to load entire directories during quick edits. Use small synthetic examples when writing or testing functions (e.g., short amino-acid strings for `ProteinDomainHMM`).
- Persist outputs to `output/` or write small RDS files at repo root only when needed (`results.rds`, `df.rds` are already used).

## Coding patterns & conventions (examples from code)

- R6 classes are used for stateful algorithms: see `ProteinDomainHMM` in `src/hmm.r`. Follow this pattern when you need object-like state and methods.
  - Constructor parameters live in `initialize(...)` and are assigned to `self$*` fields.
  - Private helper functions go into the `private` list and operate on `self`/`private` state.
- Numeric stability: `src/hmm.r` uses log-space computations (`log_sum_exp`) and small pseudocounts — preserve these practices when working on probabilistic code.
- Vectorized/base-R constructs preferred over heavy external dependencies. Add dependencies only if necessary and document them in the README.

## Workflow / dev tasks the agent may be asked to automate

- Implement a runnable pipeline: update `run.r` to source `src/rnaseq.r` and `src/hmm.r`, parse minimal CLI args (use `commandArgs(TRUE)`), and write results to `output/`.
- Add unit-like scripts: create small example scripts under `tests/` or `scripts/` that exercise functions with tiny inputs (e.g., test `ProteinDomainHMM` with a short AA sequence) — this project currently has no formal test framework.
- When changing data paths, update `data/processed/` reads/writes and keep raw files in `data/raw/` or `GDCdata/` to avoid accidental large downloads.

## Integration points & external dependencies

- R6 is required by `src/hmm.r` (guarded by `requireNamespace('R6')` in the file). If you add package dependencies, mention them in `README.md`.
- No package manager files found (no DESCRIPTION, renv, or packrat). Assume system R and `install.packages()` for ad-hoc dependencies.

## Debugging and quick run examples (concrete snippets)

- Run a quick unit of the HMM locally (small, safe example):

  source('src/hmm.r')
  hmm <- ProteinDomainHMM$new(n_states = 3, min_domain = 10, verbose = TRUE)
  res <- hmm$identify_domains('MEEPQSDPSVEPPLSQETF')
  print(res$domains)

- To add a runnable pipeline entry (suggested pattern for `run.r`): source the required `src/*.r` files, accept `commandArgs`, and write outputs under `output/` as RDS/CSV.

## When to ask humans / what not to change without confirmation

- Don't rewrite `GDCdata/` contents or delete raw TCGA directories. Those are canonical inputs.
- Large refactors that change the data layout, I/O format, or output filenames: ask for confirmation and update README + any downstream references.

## Files and locations to reference when making changes

- `src/hmm.r` — example R6 pattern, numeric stability idioms, and an algorithmic implementation.
- `src/rnaseq.r` — place for RNA-seq helpers (empty presently — put preprocessing here).
- `run.r` — pipeline entrypoint (add driver code here).
- `data/processed/`, `data/raw/`, `GDCdata/` — input data locations.
- `output/` — where analysis outputs go.

If any of the above is incomplete or you want more examples (tests, a sample `run.r`, or a requirements list), tell me which area to expand and I will update this document.
