# Which model to run pi with (Together AI, 2026-10-06)

pi 1.0.4, provider `together` (the only one configured). Each model ran
headless (`pi -p --mode json --thinking high`) in its own git worktree with
the same prompt, a 30-minute cap, and its own cargo target directory. Cost
is pi's own per-message accounting from the catalog prices. Grading was
done by hidden tests appended after the run, plus a read of every diff for
loosened tolerances and test-gaming. One run per model, so treat
differences of one check as noise.

## Task 1: debug the half-finished epidemic model (all models)

The model left by the earlier pi session had three failing tests and four
real bugs: correlated per-virion RNG streams, hash vs naive RNG
consumption, settling never applied, and a test reading the wrong `stats()`
field. Hidden checks: seed-averaged agreement with the ODE, exact
hash = naive on two new configurations (one with h not dividing 2π),
and exact settling.

| model | $ in / out per M tok | cost | time | hidden checks | notes |
|---|---|---|---|---|---|
| **DeepSeek-V4.1-Flash** | 0.30 / 1.20 | **$0.11** | 7 min | **4/4** | also found Poisson variance and periodic distance |
| Inkling | 1.00 / 4.05 | $2.80 | 16 min | 4/4 | |
| Kimi-K3 (old pi default) | 3.00 / 15.0 | $1.69 | 30 min (cap) | 4/4 | most physical fix of the front test |
| GLM-5.3 | 1.40 / 4.40 | $1.27 | 30 min (cap) | 4/4 | |
| GLM-5.2 | 1.40 / 4.40 | $1.29 | 15 min | 3/4 | missed settling |
| GLM-5.3-Flash | 0.15 / 0.50 | $0.17 | 30 min (cap) | 3/4 | hash ≠ naive when h does not divide 2π |
| DeepSeek-V4-Pro | 1.32 / 3.96 | $0.42 + $0.77 | 7 + 16 min | 4/4 on retry | 1st run killed by stream error, 2nd by credit limit |
| MiniMax-M3 | 0.30 / 1.20 | $0.33 | 30 min (cap) | 0/4 | no fix |
| DeepSeek-V4-Flash-0731 | 0.14 / 0.28 | $0.07 | 30 min (cap) | 0/4 | same loop the original session got stuck in |
| Qwen3.5-9B | 0.17 / 0.25 | $3.51 | 30 min (cap) | 0/4 | 20M uncached input tokens |

Could not be tested on this account:
- Qwen3.7-Max and Qwen3.6-Plus need third-party data sharing enabled for
  the Together organization (403).
- nemotron-3-ultra, MiniMax-M2.7 and Qwen2.5-7B need a dedicated
  endpoint (not serverless).
- Llama-3.3-70B: pi sends max_tokens = context window, which Together
  rejects once the input passes ~19k tokens (catalog entry bug).
- gpt-oss-120b: with any thinking level the second request is rejected
  ("prompt cannot be empty"; pi's replay of a reasoning-only turn). With
  thinking off it leaks its analysis channel and answers before reading
  tool output.

## Task 2: periodic geometry and a mean-field-exact kernel (finalists)

Starting from fixed task-1 code. Hidden checks: translation invariance
(exact, three hash settings) and S(t)/S₀ = exp(−(β/vol)·V·t) for uniformly
mixed virions at r_bind = 0.5, 2 and 20.

| model | cost | time | hidden | visible tests | notes |
|---|---|---|---|---|---|
| **DeepSeek-V4.1-Flash** | **$0.11** | 10 min | **6/6** | 10/10 | exact periodic-ball volume, 6 new tests; committed (e9b763a) |
| GLM-5.2 | $2.74 | 13 min | 6/6 | 11/11 | hit the credit limit on its final turn |
| GLM-5.3-Flash | $0.16 | 26 min | 6/6 | 7/9 | cut off by credit limit |
| Kimi-K3 | $1.03 | 16 min | 6/6 | 3/4 | cut off by credit limit |
| Inkling | $0.70 | 12 min | 3/6 | 3/4 | cut off by credit limit |
| GLM-5.3, DeepSeek-V4-Pro | — | — | — | — | cut off early by credit limit |

The account's credit ran out during task 2 (≈ $17.75 spent in total), so
task 3 (magnetic heat capacity of bcc iron, prepared with a hidden
heating-curve test) was not run.

## Recommendation

Use **DeepSeek-V4.1-Flash** as pi's default. It was the only model that
got everything right on both tasks, and it was the cheapest and fastest of
the models that worked: 10–25× cheaper than Kimi-K3, GLM-5.3 or Inkling
for the same or better results. GLM-5.3-Flash is the cheap backup (close
on task 2, one miss on task 1). The expensive models bought nothing
here. Kimi-K3 and GLM-5.3 were still working at the 30-minute cap.

Tool-loop reliability counts as much as model quality: stream errors,
503s and request timeouts killed or stalled runs for DeepSeek-V4-Pro,
Inkling and Kimi-K3, and pi stops the run on a stream error.
