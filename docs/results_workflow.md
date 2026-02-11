# Results Workflow Runbook

This runbook is for producing presentation-ready rescue results without repeatedly rerunning the full pipeline.

## 1) Environment check

```bash
p53cad doctor
```

## 2) Run campaign

High-budget full matrix (Big-8 singles + pairs across delivery modes):

```bash
p53cad campaign-run --budget high --seed 42 --shortlist-n 30
```

Reduced dry-run:

```bash
p53cad campaign-run --budget fast --max-scenarios 6 --shortlist-n 10
```

## 3) Inspect and regenerate report artifacts

```bash
p53cad campaign-list
p53cad campaign-report --run-id <run_id> --shortlist-n 30
```

## 4) Present from saved artifacts in Streamlit

```bash
streamlit run p53cad/app/main.py
```

In **Design Studio**:
1. Set **Result Source** to `Local campaign artifacts (auto latest)`
2. Select a run id
3. Click **Load selected campaign results**

Optional: upload a zipped run folder with `manifest.json` + parquet tables.

## 5) Export top candidates for slides/tables

Use:
- `data/campaigns/<run_id>/top30.csv`
- `data/campaigns/<run_id>/summary.md`

## 6) Organize generated files (non-destructive)

Dry run:

```bash
python scripts/organize_workspace.py
```

Apply moves:

```bash
python scripts/organize_workspace.py --apply
```

This moves loose generated files into:
- `data/processed/exports/`
- `data/processed/archive/<timestamp>/`

It does not touch `_legacy/`.
