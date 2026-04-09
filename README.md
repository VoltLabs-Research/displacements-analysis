# DisplacementsAnalysis

`DisplacementsAnalysis` computes atomic displacements relative to an optional reference frame.

## CLI

Usage:

```bash
displacement-analysis <lammps_file> [output_base] [options]
```

### Arguments

| Argument | Required | Description | Default |
| --- | --- | --- | --- |
| `<lammps_file>` | Yes | Input LAMMPS dump file. | |
| `[output_base]` | No | Base path for output files. | derived from input |
| `--reference <file>` | No | Reference LAMMPS dump file. If omitted, the current frame is used. | current frame |
| `--mic` | No | Use minimum image convention. | `true` |
| `--affineMapping <mode>` | No | Affine mapping mode: `noMapping`, `toReferenceCell`, `toCurrentCell`. | `noMapping` |
| `--threads <int>` | No | Maximum worker threads. | auto |
| `--help` | No | Print CLI help. | |

## Build With CoreToolkit

```bash
cd /path/to/voltlabs-ecosystem/tools/CoreToolkit
conan create . -nr

cd /path/to/voltlabs-ecosystem/plugins/DisplacementsAnalysis
conan create . -nr
```
