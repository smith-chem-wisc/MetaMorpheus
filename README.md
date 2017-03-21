# MetaMorpheus

[![Build status](https://ci.appveyor.com/api/projects/status/0kpjdrn9tn6y387k/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/metamorpheus/branch/master) 
[![Coverage Status](https://coveralls.io/repos/github/smith-chem-wisc/MetaMorpheus/badge.svg?branch=master)](https://coveralls.io/github/smith-chem-wisc/MetaMorpheus?branch=master)
[![Coverity Scan Build Status](https://scan.coverity.com/projects/11282/badge.svg)](https://scan.coverity.com/projects/metamorpheus)
[![Build Status](https://travis-ci.org/smith-chem-wisc/MetaMorpheus.svg?branch=master)](https://travis-ci.org/smith-chem-wisc/MetaMorpheus)

Download the current version at [https://github.com/smith-chem-wisc/MetaMorpheus/releases](https://github.com/smith-chem-wisc/MetaMorpheus/releases).
The Windows GUI release is fully operational, command line version is still in development.

MetaMorpheus is a bottom-up proteomics database search software with integrated posttranslational modification (PTM) discovery capability.
This program combines features of [Morpheus](https://github.com/cwenger/Morpheus) and [G-PTM-D](https://github.com/smith-chem-wisc/gptmd) in a single tool.

Check out the [wiki page](https://github.com/smith-chem-wisc/MetaMorpheus/wiki) for software details!

## Major Features

* Database Search: A robust dictionary search algorithm that identifies peptides by their fragmentation spectra.

* Calibration: A calibration tool that uses peptides identified by a database search to calibrate the mz values of all peaks in the spectra. This improves the quality of any subsequent search or analysis of the data.

* G-PTM-D: Post-translational modification (PTM) discovery framework, which expands the scope of peptide identifications to include both known and unknown PTMs. 

* Quantification: Ultrafast quantification of peptide abundances is enabled. Values reported are from the intensity of the most intense and most abundant isotopic form. Correlation with MaxQuant results are ~90%. Caveat Emptor!

* TMT: TMT-tagged peptides may be searched and identified by selecting the tmt.txt modification and setting it to fixed. Reporter ion abundances are not yet reported. GPTM-D discovery may be run with TMT-tagged peptides.

## Example Usage (GUI)

1. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic).

2. Open MetaMorpheusGUI.exe, and drag and drop the raw spectra files and the compressed Uniprot XML database on the GUI.

3. For all following tasks, have ptmlist.txt selected for localization, f.txt for fixed, and v.txt for variable. Add, in order: 
   * Combined 5 and 10 ppm precursor tolerance Search Task 
   * Calibration Task
   * Combined 5 and 10 ppm precursor tolerance Search Task
   * G-PTM-D Task with m.txt, metals.txt, pt.txt selected in the G-PTM-D column to add these plausible modifications to the augmented database
   * 5 ppm precursor tolerance Search Task with m.txt, metals.txt, pt.txt and ptmlist.txt selected for localization

4. Click Run All Tasks!

5. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration. 

6. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## Example Usage (Command Line)

1. Download the files at [https://uwmadison.box.com/v/MetaMorpheusPublic](https://uwmadison.box.com/v/MetaMorpheusPublic) to the folder with MetaMorpheusCommandLine.exe executable

2. Run the command 

```
MetaMorpheusCommandLine.exe -t Task1SearchExample.toml Task2CalibrationExample.toml Task3SearchExample.toml Task4GptmdExample.toml Task5SearchExample.toml -s 04-30-13_CAST_Frac4_6uL.raw 04-30-13_CAST_Frac5_4uL.raw -d uniprot-mouse-reviewed-3-9-2017.xml.gz
```

3. As the third task completes, open the results.txt files for the first and third tasks. Observe the increase in the number of confident PSMs identified due to calibration. 

4. As the fifth task completes, open the results.txt files for the third and fifth tasks. Observe the increase in the number of confident PSMs identified due to an addition of new plausible PTMs.

## System Requirements

* X GB of RAM is recommended

* For thermo .RAW files: Need to have [Thermo MSFileReader 3.1 SP2](https://thermo.flexnetoperations.com/control/thmo/search?query=MSFileReader) installed.

## References

* [Global Post-translational Modification Discovery--J. Proteome Res., 2016 Just Accepted](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)
