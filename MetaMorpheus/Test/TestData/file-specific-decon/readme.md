# FlashDeconv MS1 feature fixtures

This folder contains FlashDeconv output for three small MetaMorpheus search fixtures. The raw mzML files and the FASTA database remain in the parent `TestData` folder; the feature files here are paired by basename:

| Raw file | Feature file used for precursor deconvolution |
| --- | --- |
| `TaGe_SA_A549_3_snip.mzML` | `TaGe_SA_A549_3_snip_ms1.feature` |
| `TaGe_SA_A549_3_snip_2.mzML` | `TaGe_SA_A549_3_snip_2_ms1.feature` |
| `TaGe_SA_HeLa_04_subset_longestSeq.mzML` | `TaGe_SA_HeLa_04_subset_longestSeq_ms1.feature` |

The matching `*_ms2.feature` outputs and other FlashDeconv exports are retained as reference data, but the search test maps each raw file to its MS1 feature file. Paths are assembled from NUnit's test directory so the test does not depend on a developer-specific checkout location.

`FromFileDeconvolutionSearchTests.SearchTask_UsesFlashDeconvFeaturesForThreeRawFiles` constructs the task-map deconvolution parameters directly and verifies that each feature file loads. It runs the search once, loads the task TOML written by that run into a second `SearchTask`, reruns the same inputs, and requires the two `AllPSMs.psmtsv` files to be byte-identical. This covers the backend TOML round trip without requiring GUI-generated settings.

`FromFileDeconvolutionSearchTests.SearchTask_FileSpecificFromFile_UsesFlashDeconvFeaturesForThreeRawFiles` copies the three raw fixtures into an isolated temporary input directory and writes one companion TOML per raw (a native `FromFile` precursor deconvolution with the absolute MS1 feature path) next to each copied raw. It asserts every companion feature collection loads non-empty and that the search writes at least one PSM row.

`FromFileDeconvolutionSearchTests.SearchTask_FileSpecificFromFile_MatchesTaskMapOutputBytes` runs the task-map route and the per-file companion-TOML route against the same copied raw files, FASTA, raw order, common/search settings, charge bounds, polarity, `UseGenericScore` and feature pairings, then requires the two `AllPSMs.psmtsv` outputs to be byte-identical. Only the PSM files are compared; metadata and log files are ignored.

Replaying the file-specific route requires the task TOML **plus** the adjacent per-raw companion `<raw basename>.toml` files: the task TOML carries the shared search settings, while each companion TOML carries that raw file's native FromFile `FilePath`, charge range, polarity and `UseGenericScore`. The tests never create or modify TOMLs next to the repository fixtures; all companion TOMLs live under temporary test directories and are deleted afterward.

The generated TOML stores each raw-file basename and the corresponding feature-file path relative to that raw file's directory, rather than writing machine-specific absolute paths. For example:

```toml
[CommonParameters.PrecursorDeconvolutionParameters.FeatureFileMap]
"TaGe_SA_A549_3_snip.mzML" = "file-specific-decon/TaGe_SA_A549_3_snip_ms1.feature"
"TaGe_SA_A549_3_snip_2.mzML" = "file-specific-decon/TaGe_SA_A549_3_snip_2_ms1.feature"
"TaGe_SA_HeLa_04_subset_longestSeq.mzML" = "file-specific-decon/TaGe_SA_HeLa_04_subset_longestSeq_ms1.feature"
```

FlashDeconv retention-time values in these files are expressed in seconds. mzLib's from-file reader normalizes these values to minutes for comparison with scan retention times.
