# MetaMorpheus

[![Build status](https://ci.appveyor.com/api/projects/status/0kpjdrn9tn6y387k/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/metamorpheus/branch/master)
[![Coverage Status](https://coveralls.io/repos/github/smith-chem-wisc/MetaMorpheus/badge.svg?branch=master)](https://coveralls.io/github/smith-chem-wisc/MetaMorpheus?branch=master)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/11282/badge.svg)](https://scan.coverity.com/projects/metamorpheus)
[![Build Status](https://travis-ci.org/smith-chem-wisc/MetaMorpheus.svg?branch=master)](https://travis-ci.org/smith-chem-wisc/MetaMorpheus)
 
Download the current version at [https://github.com/smith-chem-wisc/MetaMorpheus/releases](https://github.com/smith-chem-wisc/MetaMorpheus/releases).

MetaMorpheus is a bottom-up proteomics database search software with integrated posttranslational modification (PTM) discovery capability.
This program combines features of [Morpheus](https://github.com/cwenger/Morpheus) and [G-PTM-D](https://github.com/smith-chem-wisc/gptmd) in a single tool.

Check out the [wiki page](https://github.com/smith-chem-wisc/MetaMorpheus/wiki) for software details!

## Major Features

* Database Search: A robust dictionary search algorithm that identifies peptides by their fragmentation spectra.

* Calibration: A calibration tool that uses peptides identified by a database search to calibrate the mz values of all peaks in the spectra. This improves the quality of any subsequent search or analysis of the data.

* G-PTM-D: Post-translational modification (PTM) discovery framework, which expands the scope of peptide identifications to include both known and unknown PTMs.

* Quantification: Ultrafast quantification of peptide abundances is enabled. Values reported are from the intensity of the most intense and most abundant isotopic form.

* TMT: TMT-tagged peptides may be searched and identified by selecting the tmt.txt modification and setting it to fixed. Reporter ion abundances are not yet reported. GPTM-D discovery may be run with TMT-tagged peptides.


## Requirements

* .raw or .mzML file in centroid mode
* MS2 resolution of 15,000
* 16 GB of RAM is recommended
* For thermo .RAW files: Need to have [Thermo MSFileReader 3.1 SP2](https://thermo.flexnetoperations.com/control/thmo/search?query=MSFileReader) installed.


## Example Usage (GUI)

1. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic).

2. Open MetaMorpheusGUI.exe, and drag and drop the raw spectra files and the compressed Uniprot XML database on the GUI.

3. Add search tasks that test all of the functionality of MetaMorpheus. Drag the .toml files **IN ORDER** (Task1 - Task5) onto the application. 
  * Task1SearchExample.toml - tests the standard search functionality.
  * Task2CalibrationExample.toml - will calibrate a .raw of .mzML file based on high scoring search results and write a new calibrated .mzML file.
  * Task3SearchExample.toml - searches the newly calibrated data file, which demonstrates improved performance and allows for tighter search tolerances.
  * Task4GptmdExample.toml - Searches the calibrated data file to find high probability PTMs. This search task generates a new .xml protein database.
  * Task5SearchExample.toml - Searches the calibrated input file against the G-PTM-D .xml database. This search should have the best output interms of total PSMs and total modified-peptides.

4. Click Run All Tasks!

5. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration.

6. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## Example Usage (Windows Command Line)

1. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic) to the folder with MetaMorpheusCommandLine.exe executable

2. Run the command

```
MetaMorpheusCommandLine.exe -t Task1SearchExample.toml Task2CalibrationExample.toml Task3SearchExample.toml Task4GptmdExample.toml Task5SearchExample.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-3-9-2017.xml.gz
```

3. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration.

4. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## Example Usage (Command Line with Mono)

1. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic) to the folder with MetaMorpheusCommandLine.exe executable

2. Run the command

```
mono MetaMorpheusCommandLine.exe -t Task1SearchExample.toml Task2CalibrationExample.toml Task3SearchExample.toml Task4GptmdExample.toml Task5SearchExample.toml -s example.mzML -d mouse-10.xml
```

## References

* [Global Post-translational Modification Discovery--J. Proteome Res., 2017, 16, 1383-1390](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)
