# <img src="https://user-images.githubusercontent.com/16883585/75211541-da01c680-5749-11ea-9f6c-096dc2ec4dbc.png" width="30"/> MetaMorpheus: Free, Open-Source PTM Discovery [![Follow us on Twitter](https://img.shields.io/twitter/follow/smith_chem_wisc?label=Twitter&style=social)](https://twitter.com/smith_chem_wisc)

[![Release](https://img.shields.io/github/v/release/smith-chem-wisc/MetaMorpheus)](https://github.com/smith-chem-wisc/MetaMorpheus/releases/latest)
[![Build status](https://ci.appveyor.com/api/projects/status/0jt31252xny5aoxt/branch/master?svg=true)](https://ci.appveyor.com/project/smith-chem-wisc/metamorpheus/branch/master)
[![codecov](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus/branch/master/graph/badge.svg)](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus)
[![Github All Releases](https://img.shields.io/github/downloads/smith-chem-wisc/MetaMorpheus/total.svg)](https://github.com/smith-chem-wisc/MetaMorpheus/releases)
[![Github All Releases](https://img.shields.io/docker/pulls/smithchemwisc/metamorpheus)](https://hub.docker.com/r/smithchemwisc/metamorpheus/tags?page=1&ordering=last_updated)

[![Anaconda-Server Badge](https://anaconda.org/conda-forge/metamorpheus/badges/installer/conda.svg)](https://conda.anaconda.org/conda-forge)
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
* [O-glycopeptide Characterization](https://github.com/smith-chem-wisc/MetaMorpheus/wiki/1_New-Task:-O-Glyco-Search): O-Pair Search identifies O-glycopeptides using an ion-indexed open modification search and localizes O-glycosites using graph theory and probability-based localization.

## System Requirements

* Environment:
  * 64-bit operating system
  * .NET Core 5.0:
     * Windows: https://dotnet.microsoft.com/download/dotnet/thank-you/runtime-desktop-5.0.4-windows-x64-installer
     * Mac: https://dotnet.microsoft.com/download/dotnet/thank-you/runtime-5.0.4-macos-x64-installer
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

## Test Installation (Windows GUI)

1. Download the latest MetaMorpheusInstaller.msi [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases), and install MetaMorpheus.
2. Download the example spectra and database files from [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/s/woh246eq30vtckxsfbae5pw0ltoznlca).
3. Open MetaMorpheus from the start menu, and drag and drop the .raw spectra files and the UniProt .xml database into MetaMorpheus.
4. Add a series of Tasks to make a workflow for MetaMorpheus to follow. Drag the .toml files (these files store MetaMorpheus's search parameters) (Task1 - Task5) into the application.
  * Task1-SearchTaskconfig.toml - the standard search functionality.
  * Task2-CalibrateTaskconfig.toml - will mass-calibrate the spectra file based on high scoring search results and write a new calibrated .mzML file.
  * Task3-SearchTaskconfig.toml - searches the newly calibrated data file, which demonstrates improved performance (more PSMs, lower mass errors) and allows for tighter search tolerances.
  * Task4-GPTMDTaskconfig.toml - searches the calibrated data file to find high-probability PTMs. This search task generates a new .xml protein database with annotated PTM possibilities discovered by G-PTM-D.
  * Task5-SearchTaskconfig.toml - searches the calibrated input file against the G-PTM-D .xml database. This search result is the highest confidence in terms of total PSMs and modified peptides.
5. Click "Run All Tasks!"
6. As the third task completes, open the results.txt files for the first and third tasks (before and after calibration, respectively). Observe the increase in the number of confident PSMs and identified peptides due to calibration.
7. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to discovered PTM-containing peptides.

## Typical Usage (Windows GUI)
1. Open MetaMorpheus from the start menu, and drag and drop your .raw spectra files and protein database into the GUI.
2. Select "New Calibrate Task" tab and enter appropriate search parameters, using slightly liberal mass tolerances. Then "Add the Calibration Task". Subsequent tasks (searches) will use suggested ppm tolerances automatically generated by the calibration task.
3. Select "New GPTMD Task" tab. Specify the G-PTM-D modifications that you think may be present in your sample. Many typical modifications are pre-selected. Then "Add the GPTMD Task".
4. Select "New Search Task" tab. Specify the Post-Search Parameters (e.g. protein parsimony, quantification). Then "Add the Search Task".
5. Select "Run all tasks!". This search automatically looks for PTMs uncovered in the G-PTM-D step.

## Test Installation (Windows Command Line Executable)

Please watch our <img src ="https://user-images.githubusercontent.com/16841846/40379523-eb130166-5dbb-11e8-8a03-559599cdd560.png">["How to run MetaMorpheus command line](https://www.youtube.com/watch?v=hYLe4NwZNWU) video on YouTube 
1. Download the latest [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases). Extract "MetaMorpheus_CommandLine.zip" using, for example, [7-Zip](http://www.7-zip.org/).
2. Download the example spectra and database files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/s/2u42qp0b8jllywqzeungmjj04gplw5in) to the folder with the CMD.exe executable.
3. Run the command:

```
CMD.exe -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```
4. As the third task completes, open the results.txt files for the first and third tasks (before and after calibration). Observe the increase in the number of confident PSMs identified due to calibration.
5. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## Test Installation (via .NET Core .dll - Linux, macOS, Windows)

1. Download the latest [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases). Extract files from "MetaMorpheus_CommandLine.zip".
2. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/s/2u42qp0b8jllywqzeungmjj04gplw5in) to the folder with the CMD.dll file.
3. Run the command:

* Thermo RAW files - Linux and Windows only (Thermo does not support macOS):

```
dotnet CMD.dll -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

* mzML files - Linux, macOS:

```
dotnet CMD.dll -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s mzML/04-30-13_CAST_Frac4_6uL.mzML mzML/04-30-13_CAST_Frac5_4uL.mzML -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

* mzML files - Windows

```
dotnet CMD.dll -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s mzML\04-30-13_CAST_Frac4_6uL.mzML mzML\04-30-13_CAST_Frac5_4uL.mzML -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

## Test Conda Installation (Linux, macOS, Windows)

1. Install [miniconda](https://docs.conda.io/en/latest/miniconda.html).
2. Open the terminal and enter `conda install -c conda-forge metamorpheus`.
3. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/s/2u42qp0b8jllywqzeungmjj04gplw5in).
4. Within that folder, run the command:

* Thermo RAW files - Linux and Windows only (Thermo does not support macOS):

```
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

* mzML files - Linux, macOS:

```
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s mzML/04-30-13_CAST_Frac4_6uL.mzML mzML/04-30-13_CAST_Frac5_4uL.mzML -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

* mzML files - Windows

```
metamorpheus -t Task1-SearchTaskconfig.toml Task2-CalibrateTaskconfig.toml Task3-SearchTaskconfig.toml Task4-GPTMDTaskconfig.toml Task5-SearchTaskconfig.toml -s mzML\04-30-13_CAST_Frac4_6uL.mzML mzML\04-30-13_CAST_Frac5_4uL.mzML -d uniprot-mouse-reviewed-1-24-2018.xml.gz uniprot-cRAP-1-24-2018.xml.gz
```

## mzLib


[mzLib](https://github.com/smith-chem-wisc/mzLib) is a [nuget](https://www.nuget.org/packages/mzLib/) package that we created as an all-purpose toolchest for mass-spec data analysis and many of its functions provide the tools for MetaMorpheus. mzLib is freely available for use in mass-spec applications. You do not need to download mzLib separately to run MetaMorpheus; it is already included.


## References
* [Enhanced Global Post-translational Modification Discovery with MetaMorpheus--J Proteome Res. 2018 May 4;17(5):1844-1851](https://pubs.acs.org/doi/10.1021/acs.jproteome.7b00873)
* [Global Post-translational Modification Discovery--J. Proteome Res., 2017, 16, 1383-1390](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)
