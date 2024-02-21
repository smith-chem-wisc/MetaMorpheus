# <img src="https://user-images.githubusercontent.com/16883585/75211541-da01c680-5749-11ea-9f6c-096dc2ec4dbc.png" width="30"/> MetaMorpheus: Free, Open-Source PTM Discovery [![Follow us on Twitter](https://img.shields.io/twitter/follow/smith_chem_wisc?label=Twitter&style=social)](https://twitter.com/smith_chem_wisc)

[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/MetaMorpheus)](https://github.com/smith-chem-wisc/MetaMorpheus/releases/latest)
[![Build status](https://ci.appveyor.com/api/projects/status/0jt31252xny5aoxt/branch/master?svg=true)](https://ci.appveyor.com/project/smith-chem-wisc/metamorpheus/branch/master)
[![codecov](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus/branch/master/graph/badge.svg)](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/MetaMorpheus/total.svg)](https://github.com/smith-chem-wisc/MetaMorpheus/releases)
[![Github All Releases](https://img.shields.io/docker/pulls/smithchemwisc/metamorpheus)](https://hub.docker.com/r/smithchemwisc/metamorpheus/tags?page=1&ordering=last_updated)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/metamorpheus/badges/version.svg)](https://anaconda.org/conda-forge/metamorpheus)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/metamorpheus/badges/platforms.svg)](https://anaconda.org/conda-forge/metamorpheus)
[![Anaconda-Server Badge](https://anaconda.org/conda-forge/metamorpheus/badges/downloads.svg)](https://anaconda.org/conda-forge/metamorpheus)

Download the current version [here](https://github.com/smith-chem-wisc/MetaMorpheus/releases/latest). For first-time Windows users, choose "MetaMorpheusInstaller.msi" and install MetaMorpheus. Check out our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[getting started playlist](https://www.youtube.com/playlist?list=PLVk5tTSZ1aWlYiTvJbRj6hjVDq4qH3w__) on YouTube 

MetaMorpheus is a bottom-up proteomics database search software with integrated post-translational modification (PTM) discovery capability.
This program combines features of [Morpheus](https://github.com/cwenger/Morpheus) and [G-PTM-D](https://github.com/smith-chem-wisc/gptmd) in a single tool.

Check out the [wiki page](https://github.com/smith-chem-wisc/MetaMorpheus/wiki) for software details!

## Major Features

* Database Search: A robust search algorithm that identifies peptides by their fragmentation spectra. Watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[Search task YouTube video](https://www.youtube.com/watch?v=sUM12UBJNuA)
* Calibration: A calibration tool that uses peptides identified by a database search to calibrate the m/z values of all peaks in the spectra. This improves the quality of any subsequent search or analysis of the data. Watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[calibration task YouTube video](https://www.youtube.com/watch?v=_LfiOqqqj8Q).
* G-PTM-D: Post-translational modification (PTM) discovery, which expands the scope of peptide identifications to include both known and unknown PTMs. Watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[GPTMD task YouTube video](https://www.youtube.com/watch?v=fXGT4XExLBo).
* [Quantification](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Quantification): Ultrafast label-free peptide quantification with FlashLFQ. MS2-identified peptides are used as "seeds" for peakfinding, including PTM-containing peptides. Watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[Label-free quantification with MetaMorpheus](https://www.youtube.com/watch?v=jgXRuExtuRI) video on YouTube.
* [O-glycopeptide Characterization](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/O-Glyco-Search-Task): O-Pair Search identifies O-glycopeptides using an ion-indexed open modification search and localizes O-glycosites using graph theory and probability-based localization.

## System Requirements

* Environment:
  * 64-bit operating system
  * .NET Core 6.0:
     * Windows: https://dotnet.microsoft.com/en-us/download/dotnet/thank-you/sdk-6.0.419-windows-x64-installer
     * macOS, x64 Intel processor: https://dotnet.microsoft.com/en-us/download/dotnet/thank-you/sdk-6.0.419-macos-x64-installer
     * macOS, ARM Apple Silicon processor: https://dotnet.microsoft.com/en-us/download/dotnet/thank-you/sdk-6.0.419-macos-arm64-installer
     * Linux: https://docs.microsoft.com/dotnet/core/install/linux-package-managers
* Note that the installer (MetaMorpheusInstaller.msi) only works on Windows. The command-line version of MetaMorpheus supports any operating system that supports .NET Core (Windows, MacOS, Linux)
* 8 GB RAM recommended

## Spectra Requirements

* One of the following formats:
   * Thermo .raw (Windows and Linux only)
   * .mzML file in centroid mode. Please watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[How to convert files to .mzML](https://www.youtube.com/watch?v=hOJ6ibCA5Pk) video on YouTube.
   * .mgf
* MS1 and MS2 scans
* If you would like to know more about the types of files that can be searched with MetaMorpheus, please watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[Mass Spectra Files Video](https://www.youtube.com/watch?v=SN6_T2JyxhA&list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb&index=3) on YouTube.

## Database Requirements

UniProt .XML or .fasta format; may be used in compressed (.gz) format. If you would like to know how to obtain a UniProt .XML databases, please watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">[Protein Databases Video](https://www.youtube.com/watch?v=LFvCj04r5kU&index=2&list=PLVk5tTSZ1aWlhNPh7jxPQ8pc0ElyzSUQb) on YouTube.

## Getting Started

Please check out our [wiki](https://github.com/smith-chem-wisc/MetaMorpheus/wiki) for useful information and guides.

Installation and typical usage is described for the on the [Getting Started](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started) page:
* [Test Installation (Windows GUI)](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-installation-windows-gui)
* [Typical Usage (Windows GUI)](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#typical-usage-windows-gui)
* [Test Installation (Windows Command Line Executable)](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-installation-windows-command-line-executable)
* [Test Installation (via .NET Core .dll - Linux, macOS, Windows)](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-installation-via-net-core-dll---linux-macos-windows)
* [Test Conda Installation (Linux, macOS, Windows)](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/Getting-Started#test-conda-installation-linux-macos-windows)


## References & Citation Guide for MetaMorpheus

MetaMorpheus:
* MetaMorpheus: [Enhanced Global Post-translational Modification Discovery with MetaMorpheus, J Proteome Res **2018**, _17_, 1844-1851.](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00873)
* Morpheus: [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra, J Proteome Res **2013**, _12_, 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)

GPTMD searches: 
  * [Global Post-Translational Modification Discovery, J Proteome Res **2017**, _16_, 1383–1390](https://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)
  * [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search, J Proteome Res, **2015**, _14_, 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

Quantification: 
  * [Ultrafast Peptide Label-Free Quantification with FlashLFQ, J Proteome Res **2018**, _17_, 386–391.](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00608)
  * If you use SILAC quantification: [An atlas of protein turnover rates in mouse tissues, Nat Communications **2021**, _12_, 6778.](https://www.nature.com/articles/s41467-021-26842-3)

Crosslinking MS (XL-MS) search: [Identification of MS-Cleavable and Noncleavable Chemically Cross-Linked Peptides with MetaMorpheus
, J. Proteome Res. **2018**, 17, 7, 2370–2376.](https://pubs.acs.org/doi/10.1021/acs.jproteome.8b00141)

Multiple protease parsimony: [Improved Protein Inference from Multiple Protease Bottom-Up Mass Spectrometry Data, J Proteome Res **2019**, _18_, 9, 3429–3438.](https://pubs.acs.org/doi/10.1021/acs.jproteome.9b00330)

Glycoproteomic searches: [O-Pair Search with MetaMorpheus for O-glycopeptide characterization, Nat Methods **2020**, _17_, 1133–1138.](https://www.nature.com/articles/s41592-020-00985-5)

Proteogenomic database searches with Spritz: [Spritz: A Proteogenomic Database Engine, J Proteome Res **2021**, _20_, 1826–1834.](https://pubs.acs.org/doi/10.1021/acs.jproteome.0c00407)

Long-read proteogenomic characterization: [Enhanced protein isoform characterization through long-read proteogenomics, Genome Biology **2022**, _23_, 69.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02624-y)

Spectral library GPTMD search: [A Hybrid Spectral Library and Protein Sequence Database Search Strategy for Bottom-Up and Top-Down Proteomic Data Analysis, J of Proteome Res **2022**, _21_, 2609-2618](https://pubs.acs.org/doi/10.1021/acs.jproteome.2c00305)

## mzLib, an all-purpose mass spectrometry toolchest implemented by MetaMorpheus

[mzLib](https://github.com/smith-chem-wisc/mzLib) is a [nuget](https://www.nuget.org/packages/mzLib/) package that we created as an all-purpose toolchest for mass-spec data analysis and many of its functions provide the tools for MetaMorpheus. mzLib is freely available for use in mass-spec applications. You do not need to download mzLib separately to run MetaMorpheus; it is already included.
