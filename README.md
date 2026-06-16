# Displacements Analysis

Computes per-atom displacement vectors (and optional Zimmerman slip vectors) relative to a reference frame.

## Install

```bash
vpm install @voltlabs/displacements-analysis
```

## CLI

```bash
displacement-analysis <input_dump> [output_base] [options]
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `<input_dump>` | yes | — | Input LAMMPS dump. |
| `[output_base]` | no | derived from input | Base path for output files. |
| `--reference <file>` | no | current frame | Reference LAMMPS dump. If omitted, the current frame is used. |
| `--mic` | no | `true` | Use the minimum image convention. |
| `--affine_mapping <mode>` | no | `noMapping` | Affine mapping: `noMapping`, `toReferenceCell`, or `toCurrentCell`. |
| `--compute_slip_vector` | no | `true` | Compute the Zimmerman slip vector per atom. |
| `--slip_cutoff <Å>` | no | `3.5` | Reference-frame neighbor cutoff distance for the slip vector. |
| `--slip_threshold <Å>` | no | `0.5` | Relative displacement above which a neighbor counts as slipped. |

## Exports

| Output file | Exposure | Exporter → artifact |
|---|---|---|
| `{output_base}_displacements.parquet` | Displacements | — |
| `{output_base}_atoms.parquet` | Displacements Model | AtomisticExporter → glb |

---

Full input contract and examples: https://docs.voltcloud.dev/docs/plugins/displacements-analysis
