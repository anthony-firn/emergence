# Notes for coding agents

Rust workspace in `fluid/` (crates under `fluid/crates/`); write-ups in
`analysis/`. The organising idea is in `analysis/emergence_ladder.md`: a
cheap model, a cheap validity indicator, and an unfold to a finer model
(or a closure) where the indicator fails.

- Validate every model against an exact or published result before
  claiming it works, and state the error. Report misses and negative
  results honestly; never tune a model to the number it is meant to predict.
- Cite the source of every material constant in a comment (handbook,
  paper, or "estimate").
- Tests must be fast (each well under 30 s in `cargo test`, which builds
  with opt-level 3). Never loosen a tolerance to make a test pass.
- Commands: `cd fluid && cargo test -p <crate>`; CUDA crates need
  `source fluid/scripts/cuda-env.sh`.
- Never kill processes by name or pattern (`pkill`, `killall`, `ps | grep | kill`):
  other agents and the agent's own shell carry the same words in their
  arguments. Kill only a PID you started yourself, or bound commands with
  `timeout`.
- Do not rewrite git history.
