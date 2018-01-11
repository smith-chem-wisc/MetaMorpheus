# MetaMorpheus: Why Not See All Your Peptides? <a href="https://twitter.com/intent/tweet?text=Discover PTMs in bottom-up data with MetaMorpheus https://tinyurl.com/y9an55ah"> <img src="https://img.shields.io/twitter/url/http/shields.io.svg?style=social" alt="Tweet"></a>

[![Build status](https://ci.appveyor.com/api/projects/status/0kpjdrn9tn6y387k/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/metamorpheus/branch/master)
[![Build Status](https://travis-ci.org/smith-chem-wisc/MetaMorpheus.svg?branch=master)](https://travis-ci.org/smith-chem-wisc/MetaMorpheus/branches)
[![codecov](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus/branch/master/graph/badge.svg)](https://codecov.io/gh/smith-chem-wisc/MetaMorpheus/branch/master)
[![Codacy Badge](https://api.codacy.com/project/badge/Grade/a610ab7d56904e85846fae61f79ea925?branch=master)](https://www.codacy.com/app/solntsev_2/MetaMorpheus/dashboard?bid=5044926)
<a href="https://scan.coverity.com/projects/metamorpheus">
  <img alt="Coverity Scan Build Status"
       src="https://scan.coverity.com/projects/11282/badge.svg"/>
</a>

Download the current version at [https://github.com/smith-chem-wisc/MetaMorpheus/releases](https://github.com/smith-chem-wisc/MetaMorpheus/releases).
 
MetaMorpheus is a bottom-up proteomics database search software with integrated posttranslational modification (PTM) discovery capability.
This program combines features of [Morpheus](https://github.com/cwenger/Morpheus) and [G-PTM-D](https://github.com/smith-chem-wisc/gptmd) in a single tool.

Check out the [wiki page](https://github.com/smith-chem-wisc/MetaMorpheus/wiki) for software details!

## Major Features

* Database Search: A robust dictionary search algorithm that identifies peptides by their fragmentation spectra.
* Calibration: A calibration tool that uses peptides identified by a database search to calibrate the mz values of all peaks in the spectra. This improves the quality of any subsequent search or analysis of the data.
* G-PTM-D: Post-translational modification (PTM) discovery framework, which expands the scope of peptide identifications to include both known and unknown PTMs.
* Quantification: Ultrafast quantification of peptide abundances is enabled. Values reported are from the intensity of the most intense and most abundant isotopic form.

## System Requirements

* Environment:
  * For .NET Core versions: .NET Core 2.0.0 runtime
  * For .NET Framework versions: .NET Framework 4.7.1
    * The .NET Framework versions have the ability to read Thermo .RAW files. Those require [Thermo MSFileReader](https://thermo.flexnetoperations.com/control/thmo/search?query=MSFileReader) installed, and at least an x64 Windows 7. For some .RAW files, minimum required MSFileReader version is 3.1.9 (3.1 SP4).
* At least 16 GB RAM recommended


## Spectra Requirements

* One of the following formats:
   * Thermo .raw
   * .mzML file in centroid mode
* MS and MS/MS scans

## Database Requirements

UniProt .XML or .fasta format, may be used in compressed (.gz) format.

## Test Installation (Windows GUI)

1. Download the latest [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases). Extract "MetaMorpheusGUI.zip" using, for example, [7-Zip](http://www.7-zip.org).
2. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic).
3. Open MetaMorpheusGUI.exe, and drag and drop the raw spectra files and the compressed Uniprot XML database on the GUI.
4. Add search tasks that test all of the functionality of MetaMorpheus. Drag the .toml files **IN ORDER** (Task1 - Task5) onto the application. 
  * Task1SearchExample.toml - tests the standard search functionality.
  * Task2CalibrationExample.toml - will calibrate a .raw of .mzML file based on high scoring search results and write a new calibrated .mzML file.
  * Task3SearchExample.toml - searches the newly calibrated data file, which demonstrates improved performance and allows for tighter search tolerances.
  * Task4GptmdExample.toml - Searches the calibrated data file to find high probability PTMs. This search task generates a new .xml protein database.
  * Task5SearchExample.toml - Searches the calibrated input file against the G-PTM-D .xml database. This search should have the best output interms of total PSMs and total modified-peptides.
5. Click Run All Tasks!
6. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration.
7. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.


## Typical Usage (Windows GUI)
1. Open MetaMorpheusGUI.exe, and drag and drop the raw spectra files and the compressed Uniprot XML database on the GUI.
2. Select "New Calibrate Task" tab and enter the appropriate search parameters, using slightly liberal mass tolerances (10-20 ppm). Then "Add the Calibration Task".
3. Select "New GPTMD Task" tab. Use tighter parent mass tolerance than you would in your typical search. Specify the G-PTM-D modifications that you think may be present in your sample. Many typical modifications are pre-selected. Then "Add the GPTMD Task".
4. Select "New Search Task" tab. Specify the search paramters. High-resolution data that has been calibrated can frequently use a parent mass tolerance of 5ppm or less. Specify the Post-Search-Parameters (e.g. protein parsimony, quantification). Then "Add the Search Task".
5. Select "Run all tasks!". This search automatically looks for PTMs uncovered in the G-PTM-D step.

## Test Installation (Windows Command Line)

1. Download the latest [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases). Extract "MetaMorpheusCommandLine.zip" using, for example, [7-Zip](http://www.7-zip.org/).
2. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic) to the folder with MetaMorpheusCommandLine.exe executable.
3. Run the command

```
MetaMorpheusCommandLine.exe -t Task1SearchExample.toml Task2CalibrationExample.toml Task3SearchExample.toml Task4GptmdExample.toml Task5SearchExample.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-3-9-2017.xml.gz
```
4. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration.
5. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## Test Installation (Command Line with Mono)

1. Download the latest [release](https://github.com/smith-chem-wisc/MetaMorpheus/releases). Extract files from "MonoCommandLine.tar.gz".
2. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic) to the folder with MetaMorpheusCommandLine.exe executable
3. Run the command

```
mono MetaMorpheusCommandLine.exe -t Task1SearchExample.toml Task2CalibrationExample.toml Task3SearchExample.toml Task4GptmdExample.toml Task5SearchExample.toml -s example.mzML -d mouse-10.xml
```

## mzLib


[mzLib](https://github.com/smith-chem-wisc/mzLib) is a [nuget](https://www.nuget.org/packages) package that we created as an all-purpose tool chest for mass-spec data analysis and many of its functions provide the tools for MetaMorpheus. mzLib is freely available for use in mass-spec applications.


## References

* [Global Post-translational Modification Discovery--J. Proteome Res., 2017, 16, 1383-1390](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)
