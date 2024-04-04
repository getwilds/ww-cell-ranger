
# ww-cell-ranger
[![Project Status: Experimental – Useable, some support, not open to feedback, unstable API.](https://getwilds.org/badges/badges/experimental.svg)](https://getwilds.org/badges/#experimental)

This WILDS WDL workflow serves as a basic starting point for a [Cell Ranger](https://www.10xgenomics.com/support/software/cell-ranger/latest) workflow and is intended to be a relatively simplistic demonstration within the context of the WILDS ecosystem.

## Usage

For Fred Hutch users that are new to WDL, we recommend using [PROOF](https://sciwiki.fredhutch.org/dasldemos/proof-how-to/) to submit this workflow directly to the on-premise HPC cluster, as it simplifies interaction with Cromwell and provides a user-friendly front-end for job submission and tracking. For users outside of Fred Hutch or more advanced users who who would like to run the workflow locally, command line execution is relatively straightforward: 
```
java -jar cromwell-86.jar run ww-cell-ranger.wdl --inputs ww-cell-ranger-inputs.json
```
Although Cromwell is demonstrated here, this pipeline is not specific to Cromwell and can be run using whichever WDL execution method you prefer ([miniwdl](https://github.com/chanzuckerberg/miniwdl), [Terra](https://terra.bio/), [HealthOmics](https://docs.aws.amazon.com/omics/latest/dev/workflows.html), etc.).

## Support

For questions, bugs, and/or feature requests, reach out to the Fred Hutch Data Science Lab (DaSL) at wilds@fredhutch.org, or open an issue on our [issue tracker](https://github.com/getwilds/ww-cell-ranger/issues).

## Contributing

If you would like to contribute to this WILDS WDL workflow, see our [contribution guidelines](.github/CONTRIBUTING.md) as well out our [WILDS Contributor Guide](https://getwilds.org/guide/) for more details.

## License

Distributed under the MIT License. See `LICENSE` for details.
