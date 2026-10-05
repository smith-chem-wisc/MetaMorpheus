# FlashDeconv MS1 feature fixtures

This folder contains FlashDeconv output for three small MetaMorpheus search fixtures. The raw mzML files and the FASTA database remain in the parent `TestData` folder; the feature files here are paired by basename:

| Raw file | Feature file used for precursor deconvolution |
| --- | --- |
| `TaGe_SA_A549_3_snip.mzML` | `TaGe_SA_A549_3_snip_ms1.feature` |
| `TaGe_SA_A549_3_snip_2.mzML` | `TaGe_SA_A549_3_snip_2_ms1.feature` |
| `TaGe_SA_HeLa_04_subset_longestSeq.mzML` | `TaGe_SA_HeLa_04_subset_longestSeq_ms1.feature` |

The matching `*_ms2.feature` outputs and other FlashDeconv exports are retained as reference data, but the search test maps each raw file to its MS1 feature file. Paths are assembled from NUnit's test directory so the test does not depend on a developer-specific checkout location.

`FromFileDeconvolutionSearchTests.SearchTask_UsesFlashDeconvFeaturesForThreeRawFiles` constructs the mapped deconvolution parameters directly, verifies that each feature file loads, and runs a search across all three raw files. It deliberately bypasses task-TOML parsing while validating the backend execution path.

FlashDeconv retention-time values in these files are expressed in seconds. mzLib's from-file reader normalizes these values to minutes for comparison with scan retention times.
