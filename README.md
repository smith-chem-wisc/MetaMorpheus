[![Build status](https://ci.appveyor.com/api/projects/status/0kpjdrn9tn6y387k/branch/master?svg=true)](https://ci.appveyor.com/project/stefanks/metamorpheus/branch/master) 
[![Coverage Status](https://coveralls.io/repos/github/smith-chem-wisc/MetaMorpheus/badge.svg?branch=master)](https://coveralls.io/github/smith-chem-wisc/MetaMorpheus?branch=master)
 [![Coverity Scan Build Status](https://scan.coverity.com/projects/11282/badge.svg)](https://scan.coverity.com/projects/metamorpheus) 


# MetaMorpheus
A bottom-up proteomics database search software with integrated posttranslational modification (PTM) discovery. This program integrates
features of [Morpheus](https://github.com/cwenger/Morpheus) and [G-PTMD](https://github.com/smith-chem-wisc/gptmd) in a single tool.


## Expanded Functionality
MetaMorpheus adds two major components to the original Morpheus functionality

* Calibration: A calibration tool that uses search results of identified peptides to adjust the mz values of all peaks in the spectra. The calibrated peaks have a more accurate value, and improves the quality of any subsequent search or analysis of the data.

* GPTM-D: We incorporated the posttranslational modification (PTM) discovery framework in MetaMorpheus, which expands the scope of peptide identifications to include both known and unknown PTMs. 

## Minor differences

### Custom Modifications
In addition to the small set of UniProt modifications, MetaMorpheus allows using an expanded set of arbitrary user-defined modifications throughout.

### Chain and Signal peptides
Some proteins are present in biological samples as subsequences of the complete sequence specified in the database. Since they are common, and Uniprot lists these protein fragments, we expanded the search functionality to look for those as well.

### Histogram Peak Analysis
The modification discovery component is enhanced by the automated peak analysis heuristic. Every database search result is analyzed, and for every freqeuntly occuring mass shift (determined by a peak-finding algorithm), an analysis is conducted. The results of the analysis are written in a separate file, and they include the total number of unique peptides associated with the mass shift, the fraction of decoys, mass match with any known entry in the unimod or uniport database, mass match to an amino acid addition/removal combination, mass match to a combination of higher frequency peaks, fraction of localizable targets, localization residues and/or termini, and presence of any modifications in the matched peptides. All of this data can then be used to determine the nature of the peak, and the characteristics of the corresponding modificaiton. 

###General Requirements
The following files must be present in the folder with the executable. If not, they are automatically downloaded (to update a file to a newer version, delete it, and the application will download a new version). This is the only network usage by the application. 

* uniprot.xml: A UniProt reference database in .xml format

  [http://www.uniprot.org](http://www.uniprot.org)

* ptmlist.txt: A PTM library
 
  [http://www.uniprot.org/docs/ptmlist.txt](http://www.uniprot.org/docs/ptmlist.txt) 

* for thermo .RAW files: [Thermo MSFileReader](https://thermo.flexnetoperations.com/control/thmo/search?query=MSFileReader)

###System Requirements and Usage
- 8 GB of RAM is recommended

###References

* [Global Post-translational Modification Discovery--J. Proteome Res., 2016 Just Accepted](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.6b00034)

* [Global Identification of Protein Post-translational Modifications in a Single-Pass Database Search--J. Proteome Res., 2015, 14 (11), pp 4714–4720](http://pubs.acs.org/doi/abs/10.1021/acs.jproteome.5b00599)

* [A Proteomics Search Algorithm Specifically Designed for High-Resolution Tandem Mass Spectra--J. Proteome Res., 2013, 12 (3), pp 1377–1386](http://pubs.acs.org/doi/abs/10.1021/pr301024c)
