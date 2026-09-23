# Curved Nonogram Generator

A Python-based tool for generating **curved nonograms** from SVG input files.

The application processes vector paths from an SVG image, generates bridge and connection curves, applies optional optimization steps, and exports the resulting curved nonogram as an SVG file.

## Features

- Import SVG vector graphics
- Generate curved nonogram structures
- Automatic bridge curve generation
- Connection curve generation
- Multiple optimization strategies
- Configurable geometric constraints
- Optional filled solution output
- SVG export

---


## Project Structure

```text
.
├── main.py
├── inputs/
├── outputs/
├── Nonogram_settings.py
└── Nonogram_generation_functions.py
```

---

## Usage

### Basic Execution

```bash
python main.py
```

Default input/output paths:

- Input: `inputs/muffin.svg`
- Output: `outputs/output.svg`

### Specify Input and Output Files

```bash
python main.py \
    --input-file inputs/example.svg \
    --out-file outputs/result.svg
```

### Generate a Filled Solution

```bash
python main.py --fill
```

---

## Command Line Parameters

### Input / Output

| Parameter | Description | Default |
|------------|-------------|----------|
| `--input-file` | Input SVG file | `inputs/muffin.svg` |
| `--out-file` | Output SVG file | `outputs/output.svg` |
| `--fill` | Generate an additional filled solution image | `False` |

---

### Bridge Curve Generation

| Parameter | Description | Default |
|------------|-------------|----------|
| `--all-bridge-curves` | Consider all ECEs for bridge generation | `False` |
| `--inner-bridge-curves` | Consider only inward-facing ECEs | `True` |
| `--bridge-curves-only` | Generate bridge curves only | `False` |
| `--find-best-bridge` | Select the bridge with the lowest penalty | `True` |

---

### Optimization

| Parameter | Description | Default |
|------------|-------------|----------|
| `--do-checks` | Enable optimization checks | `True` |
| `--early-checks` | Optimize immediately after curve generation | `True` |
| `--check-angles` | Evaluate intersection angles | `True` |
| `--check-intersections` | Evaluate intersection spacing | `True` |
| `--check-distance` | Evaluate curve distances | `True` |
| `--additional-curves` | Generate additional curves | `True` |

---

### Geometric Constraints

| Parameter | Description | Default |
|------------|-------------|----------|
| `--max-loops` | Maximum optimization iterations per curve | `30` |
| `--min-angle` | Minimum allowed intersection angle | `20` |
| `--min-dist` | Minimum allowed distance between intersections | Automatically calculated as `min(viewBox width, viewBox height) / 20` |
| `--min-gap` | Minimum allowed gap between curves | Automatically calculated as `0.2 × min-dist` |
| `--max-dist` | Maximum allowed gap between border endpoints for additional curve generation | Automatically calculated as `(viewBox width + viewBox height) / 5` |
| `--control-point-dist` | Control point spacing on generated curves | `0.3` |

---

### Penalty Weights

| Parameter | Description | Default |
|------------|-------------|----------|
| `--pen-threshold` | Minimum acceptable penalty | `1` |
| `--w-angle` | Weight of angle penalty | `1` |
| `--w-intersection` | Weight of intersection penalty | `3` |
| `--w-distance` | Weight of distance penalty | `1` |

---

### Angular Parameters

| Parameter | Description | Default |
|------------|-------------|----------|
| `--alpha` | Outer angle constraint (radians) | `π/3` |
| `--beta` | Inner angle constraint (radians) | `π/20` |

---

## Example

```bash
python main.py \
    --input-file inputs/logo.svg \
    --out-file outputs/logo_nonogram.svg \
    --fill \
    --min-angle 25 \
    --max-loops 50
```

---

## Generation Pipeline

The program follows these steps:

1. Load and preprocess the SVG input.
2. Determine generation and optimization parameters.
3. Calculate ECEs (edge connection elements).
4. Generate bridge curves.
5. Generate connection curves.
6. Export the final SVG output.

---

## Output

The generator creates an SVG file containing the curved nonogram structure.

Example output:

```text
outputs/output.svg
```

Optionally, an additional filled solution image can be generated using the `--fill` option.

---

## Author

Created by Johanna Klinger.
