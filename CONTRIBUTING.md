# Contributing

Thank you for your interest in contributing to **scCRIPTURE**.

## Reporting Issues

Please open a GitHub Issue with:

- A clear title describing the problem
- The pipeline step/module where the issue occurs
- Relevant log output (from `logs/` directory)
- Your environment details (`conda list`, R `sessionInfo()`)

## Submitting Changes

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-improvement`)
3. Make your changes following the existing code style
4. Test your changes on a small dataset
5. Commit with clear messages (`git commit -m "Add support for X in module Y"`)
6. Push to your fork and open a Pull Request

## Code Style

- **R scripts**: Follow Tidyverse style guide; use `snake_case` for variables
- **Python scripts**: Follow PEP 8; use `snake_case`
- **SLURM scripts**: Include header comments describing inputs/outputs; use the standard helper function block
- **New modules**: Follow the existing naming convention (`XX_module_name.R`) and register in `run_pipeline.R`

## Adding a New Pipeline Module

1. Create the module script in `modules/` following the naming pattern
2. Add the module to the `run_pipeline.R` execution logic
3. Add corresponding parameters to `params.R`
4. Update the configurator HTML if applicable
5. Document the module in this README
