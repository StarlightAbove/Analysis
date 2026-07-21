# docker-images — CNV reproducibility pipeline

Two-image pipeline that reproduces the CNV-calling comparison (MethylMaster,
Conumee, Sesame on FFPE tissue) from the baked-in `LabData/`.

## Why two images

MethylMaster requires a **patched** sesame (`StarlightAbove/sesame@patch-2`) at
R 4.1.2; Conumee, unpatched Sesame, and `analyticFile.R` all target
**R 4.5.2 / Bioc 3.21** (pinned in `renv.lock`). Patched and unpatched sesame
cannot share one R library, so the stages are quarantined:

| Image | Base | Runs | Emits |
|---|---|---|---|
| `cnv-methylmaster` | `mmariani123/methylmaster` (R 4.1.2) | MethylMaster arm | `Outputs/MethylMaster/**` |
| `cnv-analysis` | `rocker/tidyverse:4.5` + `renv::restore()` | Conumee + Sesame arms, then `analyticFile.R` | `Outputs/{Conumee,SeSAMe}/**`, `quarterCutoff/**` |

## Run

From this directory:

```bash
docker compose up --build
```

Stage 1 runs to completion; the compose completion gate
(`condition: service_completed_successfully`) then releases Stage 2. `Outputs/`
is a shared bind-mount handoff; final metrics and figures land in
`quarterCutoff/` on the host.

## Data flow

```
cnv-methylmaster ─▶ Outputs/MethylMaster ─┐
                                          ├─▶ [shared Outputs/] ─▶ cnv-analysis
cnv-analysis     ─▶ Outputs/{Conumee,SeSAMe} ─┘   ├─ conumee.R
                                                  ├─ sesame.R
                                                  └─ analyticFile.R ─▶ quarterCutoff/
```

## Build context

The build context is the **Analysis project root** (`../..`). Add a
`.dockerignore` there (one is provided) to keep `FinalSet.zip`, `.git`, caches,
and other bulk out of the build.

## Reproducibility hardening (for archival / publication)

- **Pin base images by digest** (`@sha256:...`) — `:latest` and `:4.5` drift.
- **Pin `StarlightAbove/sesame` to a commit SHA**, not the `patch-2` branch.
- **Pin `cnAnalysis450k`** to a commit SHA instead of `master`.
- **Freeze the Posit Package Manager snapshot date** in the analysis image so
  `renv::restore()` resolves identically over time.
- **Archive the *built* images** (`docker save` → Zenodo, DOI) — the definitive
  guarantee, since `mmariani123/methylmaster:latest` is not under your control.
