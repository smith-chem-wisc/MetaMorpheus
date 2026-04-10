# Isotopic Envelope Visualization in MetaDraw

## Overview
Enable MetaDraw to display the entire isotopic envelope for fragment ion annotations, rather than just annotating the single closest peak. All peaks in the envelope will be colored, with the text annotation (e.g., "b7") appearing only on the tallest peak.

## Background
- mzLib PR 1030 already adds `MatchedFragmentIonWithEnvelope` class with optional `Envelope` property
- The search engine performs deconvolution to find fragment peaks - we will re-use this machinery in MetaDraw
- Currently, MetaDraw finds the closest peak to the theoretical m/z and annotates only that peak

## Implementation Phases

### Phase 1: TSV Reading with Envelope Placeholder
**File: `MetaMorpheus/GuiFunctions/MetaDraw/MetaDrawDataLoader.cs`**

When loading PSMs from TSV, use `SpectrumMatchParsingParameters` to ensure `MatchedFragmentIonWithEnvelope` objects are created:

```csharp
var parsingParams = new SpectrumMatchParsingParameters
{
    FragmentIonsHavePlaceholderForEnvelope = true,
    ParseMatchedFragmentIons = true
};

var psms = SpectrumMatchTsvReader.ReadPsmTsv(filePath, out warnings, parsingParams);
```

This ensures all read-in matched ions have the type that can hold the envelope (cached internally).

---

### Phase 2: MetaDraw Settings
**File: `MetaMorpheus/GuiFunctions/MetaDraw/MetaDrawSettings.cs`**

Add setting to control envelope annotation:

```csharp
public static bool AnnotateIsotopicEnvelopes { get; set; } = true;
```

Add to `MakeSnapShot()` and `LoadSettings()` for persistence.

---

### Phase 3: SpectrumMatchPlot Envelope Logic
**File: `MetaMorpheus/GuiFunctions/MetaDraw/SpectrumMatch/SpectrumMatchPlot.cs`**

#### 3.1: Helper method for lazy envelope calculation

Add method to calculate envelopes on-demand for a scan:

```csharp
private void PopulateEnvelopesForScanIfNeeded(MsDataScan scan, List<MatchedFragmentIon> ions)
{
    // Get unique products from matched ions that need envelopes
    var ionsNeedingEnvelopes = ions
        .OfType<MatchedFragmentIonWithEnvelope>()
        .Where(p => p.Envelope == null)
        .ToList();
    
    if (ionsNeedingEnvelopes.Count == 0)
        return;
    
    var products = ionsNeedingEnvelopes
        .Select(p => p.NeutralTheoreticalProduct)
        .Distinct()
        .ToList();
    
    // Get deconvolution parameters from MetaDraw settings
    var commonParams = new CommonParameters(
        productDeconParams: MetaDrawSettingsViewModel.Instance.DeconHostViewModel
            .ProductDeconvolutionParameters.Parameters,
        productMassTolerance: new PpmTolerance(20)); // TODO: get from settings
    
    var specificMass = new Ms2ScanWithSpecificMass(scan, scan.MassSpectrum, 
        scan.Ms2ScanNumber, scan.FileName, commonParams);
    
    // Run fragment ion matching - returns envelopes in MatchedFragmentIonWithEnvelope
    var allMatches = MetaMorpheusEngine.MatchFragmentIons(specificMass, products, commonParams);
    
    // Build lookup for products to envelopes
    var productToMatch = allMatches
        .OfType<MatchedFragmentIonWithEnvelope>()
        .Where(p => p.Envelope != null)
        .ToDictionary(p => p.NeutralTheoreticalProduct);
    
    // CACHE the envelope in each ion
    foreach (var ion in ionsNeedingEnvelopes)
    {
        if (productToMatch.TryGetValue(ion.NeutralTheoreticalProduct, out var match))
        {
            ion.Envelope = match.Envelope;
        }
    }
}
```

#### 3.2: Update AnnotateMatchedIons to trigger lazy population

In the `AnnotateMatchedIons` method (or before calling it), add:

```csharp
if (MetaDrawSettings.AnnotateIsotopicEnvelopes)
{
    PopulateEnvelopesForScanIfNeeded(Scan, matchedFragmentIons);
}
```

#### 3.3: Update AnnotatePeak to draw envelopes

Modify `AnnotatePeak()` method:

```csharp
protected void AnnotatePeak(MatchedFragmentIon matchedIon, bool isBetaPeptide, ...)
{
    // Get color (existing logic)...
    
    // Check for cached envelope
    if (MetaDrawSettings.AnnotateIsotopicEnvelopes && 
        matchedIon is MatchedFragmentIonWithEnvelope envIon && 
        envIon.Envelope != null)
    {
        var envelope = envIon.Envelope;
        
        // Find tallest peak in envelope
        var tallestPeak = envelope.Peaks.MaxBy(p => p.intensity);
        
        // Annotate ALL peaks in envelope with color
        foreach (var (mz, intensity) in envelope.Peaks)
        {
            bool isTallest = (mz == tallestPeak.mz && intensity == tallestPeak.intensity);
            
            // Draw peak marker (colored)
            DrawPeak(mz, intensity, MetaDrawSettings.StrokeThicknessAnnotated, 
                ionColor, peakAnnotation: null);
            
            // Draw TEXT annotation only on tallest peak
            if (isTallest)
            {
                // Create annotation text (existing logic)
                var peakAnnotationText = BuildAnnotationText(matchedIon, ...);
                
                var peakAnnotation = new TextAnnotation
                {
                    Text = peakAnnotationText,
                    TextPosition = new DataPoint(mz, intensity),
                    // ... other properties
                };
                
                DrawPeak(mz, intensity, ..., peakAnnotation);
            }
        }
        
        return;
    }
    
    // Fallback: current behavior (closest peak only)
    // ... existing code
}
```

---

### Phase 4: Optional UI Toggle
**File: `MetaMorpheus/GUI/MetaDraw/MetaDrawSettingsWindow.xaml[.cs]`**

Add checkbox to toggle `AnnotateIsotopicEnvelopes`

---

## Key Design Decisions

1. **Caching in the ion**: Once `MatchedFragmentIonWithEnvelope.Envelope` is populated, it's cached there - no separate cache needed. Subsequent annotations use the cached value.

2. **Lazy calculation**: Envelopes are calculated on first display of each PSM, not during initial load. This keeps initial load fast.

3. **Batch processing**: When envelopes are needed for a scan, all envelopes for that scan are calculated in a single `MatchFragmentIons` call, not per-ion.

4. **Text on tallest**: Only the tallest peak in the envelope gets the text annotation (e.g., "b7"); all peaks get colored markers.

5. **Backward compatible**: If `AnnotateIsotopicEnvelopes` is false (or if envelopes couldn't be calculated), falls back to current behavior.

---

## Files to Modify

| Phase | File | Change |
|-------|------|--------|
| 1 | `GuiFunctions/MetaDraw/MetaDrawDataLoader.cs` | Use `SpectrumMatchParsingParameters` with `FragmentIonsHavePlaceholderForEnvelope=true` |
| 2 | `GuiFunctions/MetaDraw/MetaDrawSettings.cs` | Add `AnnotateIsotopicEnvelopes` setting |
| 3 | `GuiFunctions/MetaDraw/SpectrumMatch/SpectrumMatchPlot.cs` | Add envelope calculation + caching, update annotation |
| 4 | `GUI/MetaDraw/MetaDrawSettingsWindow.xaml[.cs]` | Optional: add UI toggle |
