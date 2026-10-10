# FlashDeconv MS1 feature fixtures

This folder contains FlashDeconv output for three small MetaMorpheus search fixtures. The raw mzML files and the FASTA database remain in the parent `TestData` folder; the feature files here are paired by basename:

| Raw file | Feature file used for precursor deconvolution |
| --- | --- |
| `TaGe_SA_A549_3_snip.mzML` | `TaGe_SA_A549_3_snip_ms1.feature` |
| `TaGe_SA_A549_3_snip_2.mzML` | `TaGe_SA_A549_3_snip_2_ms1.feature` |
| `TaGe_SA_HeLa_04_subset_longestSeq.mzML` | `TaGe_SA_HeLa_04_subset_longestSeq_ms1.feature` |

The matching `*_ms2.feature` outputs and other FlashDeconv exports are retained as reference data, but the search test maps each raw file to its MS1 feature file. Paths are assembled from NUnit's test directory so the test does not depend on a developer-specific checkout location.

`FromFileDeconvolutionSearchTests.SearchTask_UsesFlashDeconvFeaturesForThreeRawFiles` constructs the mapped deconvolution parameters directly and verifies that each feature file loads. It runs the search once, loads the task TOML written by that run into a second `SearchTask`, reruns the same inputs, and requires the two `AllPSMs.psmtsv` files to be byte-identical. This covers the backend TOML round trip without requiring GUI-generated settings.

The generated TOML stores each raw-file basename and the corresponding feature-file path relative to that raw file's directory, rather than writing machine-specific absolute paths. For example:

```toml
[CommonParameters.PrecursorDeconvolutionParameters.FeatureFileMap]
"TaGe_SA_A549_3_snip.mzML" = "file-specific-decon/TaGe_SA_A549_3_snip_ms1.feature"
"TaGe_SA_A549_3_snip_2.mzML" = "file-specific-decon/TaGe_SA_A549_3_snip_2_ms1.feature"
"TaGe_SA_HeLa_04_subset_longestSeq.mzML" = "file-specific-decon/TaGe_SA_HeLa_04_subset_longestSeq_ms1.feature"
```

FlashDeconv retention-time values in these files are expressed in seconds. mzLib's from-file reader normalizes these values to minutes for comparison with scan retention times.
