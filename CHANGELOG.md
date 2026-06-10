# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project
adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-06-11

First public release — the MATLAB implementation of both planners and the shared coverage pipeline.

### Added
- **OnlineSCC** — online planner for *unknown* crack locations: boustrophedon (zig-zag) scanning
  with on-the-fly crack detection (skeleton endpoints / branch points), visibility-graph +
  Chinese Postman crack-fill planning, and iterative re-decomposition of the remaining free space
  until coverage is complete.
- **SCC** — offline planner for *known* crack maps: up-front crack-graph extraction combined with
  the coverage graph and routed with a modified Chinese Postman (Eulerian routing, with
  `intlinprog` matching of odd-degree nodes when needed).
- Shared coverage pipeline: Morse cell decomposition (MCD), Reeb graph, Reeb-path cell ordering,
  and boustrophedon path generation with cell-connection for disjoint regions.
- `private/` helpers: crack-skeleton tracers, morphology, visibility / line-of-sight tests,
  polygon decimation, and the polygon-intersection suite.
- Configurable robot dimensions via `robot_config.json` (read by `config_loader`), SCC route
  selection (`route`, default the geometry-routed `'rpp'`), and uniform / Gaussian map import.
- Visualization options (`animate` / `writeFrames` / `makeGif`): live on-screen animation of the
  robot scanning and filling, per-step frame export, in-memory animated-GIF assembly, and a
  final-path PNG written to `Results/SCC` / `Results/OnlineSCC`.
- Crack-map dataset under `CrackMaps/` (uniform and Gaussian-distributed maps across densities).

[1.0.0]: https://github.com/vveerar1/Crack-Filling-Robot-MATLAB/releases/tag/v1.0.0
