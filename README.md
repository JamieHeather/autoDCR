# autoDCR
#### Jamie Heather, MGH, 2025

[![PyPI - Version](https://img.shields.io/pypi/v/autoDCR?color=%239467bd)](https://pypi.org/project/autoDCR/) [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) ![Static Badge](https://img.shields.io/badge/experimental%20release-8A2BE2)


`autoDCR` (short for **auto**matic **D**e**c**ombinato**r**) is a python script to perform T cell receptor (TCR) gene annotation. This is inspired by and in part built upon the core functionality of [Decombinator](https://github.com/innate2adaptive/Decombinator), the TCR analysis software developed by the Chain lab at UCL. It uses a similar conceptual framework of using fast Aho-Corasick tries to search for the presence of 'tag' sequences in DNA reads, and use these to identify V and J TCR genes. However it applies that core concept in different ways, to perform several niche functions that are not well catered to in other TCR annotation pipelines.

**Note** that `autoDCR` is under development and should be considered experimental, specifically aiming to cater to specific case uses. [The documentation can be found here: https://jamieheather.github.io/autoDCR/](https://jamieheather.github.io/autoDCR/).

The 0.2.7 version used in prior publications [can be accessed via the releases page](https://github.com/JamieHeather/autoDCR/releases/tag/v0.2.7).