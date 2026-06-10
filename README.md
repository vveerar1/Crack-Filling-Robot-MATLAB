# Crack-Filling Robot — MATLAB

Path-planning and simulation for an autonomous **field robot that scans a surface for complete
coverage and fills every crack it finds**. The robot carries a 360° sensor for crack detection and
a nozzle mounted on an XY gantry within its footprint for dispensing filler.

Two planners are provided:

- **SCC (Sensor-based Complete Coverage)** — *offline*. The environment and crack locations are
  known in advance; the planner computes a low-cost coverage-and-filling path and the robot
  executes it.
- **OnlineSCC** — *online*. Crack locations are **unknown**. The robot scans the area in a
  boustrophedon (zig-zag) pattern, and whenever the sensor detects a crack it plans and fills it
  before resuming the scan — repeating until the whole area is covered.

Work areas are supplied as binary **crack maps** (PNG images), where **each pixel corresponds to a
robot pose**.

---

## How it works

Both planners share a coverage pipeline and differ only in how cracks are discovered and filled:

```
crack map (binary image)
  │
  ├─ Morse Cell Decomposition (MCD)   split the free space into monotone cells at critical points
  ├─ Reeb graph                       one node per critical point, one edge per cell
  ├─ Reeb path                        order the cells for efficient, low-overlap traversal
  └─ Boustrophedon path               generate the zig-zag coverage path at sensor spacing
                                      (with cell-connection for disjoint regions)
```

**Crack filling.**
- In **SCC**, the full crack network is known, so a crack graph is extracted up front and a
  **modified Chinese Postman** routine (Eulerian routing, with integer-programming matching of
  odd-degree nodes when needed) produces the shortest closed filling route.
- In **OnlineSCC**, the robot senses cracks within its sensor radius as it scans, extracts the
  crack skeleton and its endpoints/branch points on the fly, and plans crack-filling waypoints
  using a visibility graph and Chinese Postman routing. Covered area is removed from the free
  space, and the remaining region is re-decomposed until coverage is complete.

The planners report coverage path length, area covered, overlap, and run time.

---

## Requirements

- **MATLAB R2025b** (or a recent release).
- Toolboxes:
  - **Image Processing Toolbox** — skeletonization and morphology (`bwskel`, `bwmorph`,
    `imbinarize`).
  - **Optimization Toolbox** — integer-programming matching of odd-degree nodes (`intlinprog`).
  - **Curve Fitting Toolbox** — cubic-spline Reeb-edge curves (`csapi`).
- Polygon geometry (`polyshape`, `polybuffer`) and graph routines (`graph`, `shortestpath`) are
  part of base MATLAB.

---

## Getting started

1. **Clone** the repository and open the folder in MATLAB:
   ```matlab
   % from the MATLAB command window, in the repository folder
   addpath(genpath(pwd))
   ```
   The pipeline helpers live in `private/` and resolve automatically.
2. **Run a planner:**
   ```matlab
   OnlineSCC    % online — unknown cracks
   SCC          % offline — known cracks
   ```

Each run runs on the configured map and prints a results row. Which map, which route, and whether to
animate are all set by editing a few variables at the top of `OnlineSCC.m` / `SCC.m`, as described
next.

---

## Configuring a run

Open `OnlineSCC.m` or `SCC.m` and edit the variables in the configuration block near the top.

### Choosing a crack map

Maps live in `CrackMaps/`: uniform maps are named `myCrack<index>_<density>_<sample>.png`, Gaussian
maps `myCrackGauss_s<sigma>_<density>`. The map name is assembled from these variables:

| Variable | Meaning |
|---|---|
| `Gau` | `0` = uniform map, `1` = Gaussian map |
| `den` | the density list `[35 45 50 65 80 90 95 100]` (leave as-is) |
| `k` (in `SCC.m`) / `dd` (in `OnlineSCC.m`) | density index into `den` — e.g. `k = 4` → 65 % density |
| `mapN` | uniform map number (sample variant) |
| `Gaussb` | Gaussian source folder (`1`–`8`), used when `Gau = 1` |
| `b` | sigma index into `sig = [5 10 20]`, used when `Gau = 1` |

For example, `OnlineSCC.m` ships with `Gau = 0; dd = 8; mapN = 1;` → `myCrack8_100_1`; set
`Gau = 1; Gaussb = 6; b = 3;` for the Gaussian map `Gaussian6/myCrackGauss_s20_100`.

### Routing (SCC only)

`route` selects the coverage-route formulation: `'rpp'` (default — geometry-routed cell
Rural-Postman) or `'native'` (Reeb-graph Chinese Postman).

### Visualization

Every run always saves a styled **final-path figure** when it finishes (see Outputs). The live
animation and any saved frames/GIF are controlled by the flags in the `%% Visualization` block:

| Variable | Effect |
|---|---|
| `animate` | `true` = draw the live animation figure as the robot scans/fills (nothing written to disk); `false` = headless (default) |
| `makeGif` | `true` = stitch the animation into `Results/GIF/<map>.gif` |
| `writeFrames` | `true` = also save each step as `Results/GIF/<map>_<n>.png` |
| `fps` | GIF frame rate |

Setting `makeGif` or `writeFrames` automatically turns the animation figure on. With all three
`false` the run is fully headless (fastest) and still produces the final-path figure.

Robot dimensions (base / footprint / sensor diameters) are configured in `robot_config.json`.

---

## Outputs

| What | Where |
|---|---|
| Results row (printed) | the MATLAB command window |
| Final-path figure (always) | `Results/SCC/<map>_SCC.png` (SCC) / `Results/OnlineSCC/<map>_OnlineSCC.png` (OnlineSCC) |
| Animation GIF (when `makeGif = true`) | `Results/GIF/<map>.gif` |
| Per-step frames (when `writeFrames = true`) | `Results/GIF/<map>_<n>.png` |

The `Results/` subfolders are created automatically. The results row is:

```
res = [ num_iterations , density , runtime_s , path_length_ft , area_covered_ft2 ]
```

along with coverage and overlap percentages.

---

## Notes

- Distances are computed in image pixels; `1 px ≈ 2 mm`. The robot footprint, sensor range, and
  base size are configured in `robot_config.json` (in inches, converted to pixels) and read by both
  planners via `config_loader`.

---

## Author

Vishnu Veeraraghavan — Automated Control Systems and Robotics Lab.
For questions, contact `vveerar1@binghamton.edu`.
