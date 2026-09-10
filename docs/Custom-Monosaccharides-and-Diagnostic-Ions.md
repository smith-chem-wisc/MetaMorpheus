# Custom Monosaccharides & Diagnostic (Oxonium) Ions

A guide to extending MetaMorpheus glyco search with **your own monosaccharides** and **your own diagnostic oxonium ions**, without recompiling.

This covers two related features:

1. **Custom monosaccharides** — register sugars MetaMorpheus does not ship with, so they can appear in your glycan databases.
2. **Diagnostic / oxonium ions** — attach observed oxonium-ion *m/z* values to a monosaccharide so they boost scoring and (optionally) act as a **strict spectrum filter**.

Both are configured in a single tab-separated file: **`MonosaccharidesCustom.tsv`**.

---

## Contents

- [Quick start](#quick-start)
- [Two ways in: the GUI dialog or the file](#two-ways-in-the-gui-dialog-or-the-file)
- [Where the file lives](#where-the-file-lives)
- [The file format](#the-file-format)
- [Built-in monosaccharides (do not redefine)](#built-in-monosaccharides-do-not-redefine)
- [Part 1 — Define a custom monosaccharide](#part-1--define-a-custom-monosaccharide)
- [Part 2 — Use it in a glycan search](#part-2--use-it-in-a-glycan-search)
- [Part 3 — Diagnostic & oxonium ions (column 4)](#part-3--diagnostic--oxonium-ions-column-4)
- [End-to-end worked example](#end-to-end-worked-example)
- [Validation & error messages](#validation--error-messages)
- [Tips & caveats](#tips--caveats)
- [FAQ](#faq)

---

## Quick start

1. In the MetaMorpheus GUI, open the **Settings** tab and click **Create new monosaccharide**. Fill in the **Name**, **Single-Char Code**, a **Chemical Formula** *and/or* **Monoisotopic Mass**, and — for this feature — the **Diagnostic Ions (m/z)** box, then **Save Monosaccharide**.
   (Prefer a text editor? Add the row to `MonosaccharidesCustom.tsv` by hand instead — see [Two ways in](#two-ways-in-the-gui-dialog-or-the-file).)
2. Reference the new monosaccharide (by name or code) in a glycan database file (`.gdb`).
3. Restart MetaMorpheus, point your Glyco Search task at that database, and run.
4. To make the diagnostic ions *filter* spectra, leave **`OxoniumIonFilt`** checked in the Glyco Search task.

```
Name	SingleCharCode	MonoisotopicMass	DiagnosticIonMasses	Description
HexA	U	176.03209	175.02482,157.01425	Hexuronic acid (GlcA, GalA)
```

---

## Two ways in: the GUI dialog or the file

Both routes write the same row to the same `MonosaccharidesCustom.tsv`, so pick whichever you prefer. The GUI is the shorter path for one sugar; the file is better for bulk edits, for version-controlling a lab-wide set, or for adding the explanatory comments the parser ignores.

### The GUI dialog

Main window → the **Settings** tab → **Create new monosaccharide**. The dialog has one field per TSV column:

| Dialog field | TSV column | Required | Notes |
|---|---|:---:|---|
| **Name** | 1 · `Name` | ✅ | e.g. `HexA`. Used in composition-format glycans. |
| **Single-Char Code** | 2 · `SingleCharCode` | ✅ | One ASCII letter, e.g. `U`. Used in structure-format glycans. |
| **Chemical Formula** | — | ⬜ | e.g. `C6H8O6`. Converted to a mass on save; **takes precedence** over the mass box if both are filled. |
| **Monoisotopic Mass** | 3 · `MonoisotopicMass` | ✅\* | Residue mass in Da, e.g. `176.03209`. \*Required unless you gave a formula. |
| **Diagnostic Ions (m/z)** | 4 · `DiagnosticIonMasses` | ⬜ | Comma-separated observed m/z, e.g. `175.02482,157.01425`. **These are the custom oxonium ions** this guide is about — see [Part 3](#part-3--diagnostic--oxonium-ions-column-4). |
| **Description** | 5 · `Description` | ⬜ | Free text; cite a source or note a formula. |

**Save Monosaccharide** validates the entry, registers it for the current session, and appends it to the file. The same validation rules and error messages apply as when the file is loaded at startup (see [Validation & error messages](#validation--error-messages)) — including the rule that a diagnostic ion may not duplicate a built-in oxonium ion.

> If the save reports that the monosaccharide is usable this session but could not be written to file, the file was not writable; fix the permissions and re-enter it, or the sugar will be gone next launch.

### The file

Same folder as everything else in the data directory — see below. **Open mods/data folder** on the same Settings tab takes you straight there.

---

## Where the file lives

`MonosaccharidesCustom.tsv` sits at the **top level of the MetaMorpheus data directory**, alongside `proteases_custom.tsv` and `rnase_custom.tsv`. (It used to live in the `Glycan_Mods` subfolder; it was moved so that an install, repair or upgrade cannot overwrite it. If you have an old copy under `Glycan_Mods`, MetaMorpheus copies it to the new location once, on first launch, and leaves the original alone.)

The quickest way there is the **Settings** tab → **Open mods/data folder**.

| Install type | Typical path |
|---|---|
| Windows GUI install | `%LOCALAPPDATA%\MetaMorpheus\MonosaccharidesCustom.tsv` |
| Built from source / portable | `<build output folder>\MonosaccharidesCustom.tsv` |
| Custom data directory | `<your data dir>\MonosaccharidesCustom.tsv` |

If the file is missing, MetaMorpheus writes a fresh documented template there on startup. It ships with documentation and commented-out examples and **no active rows**, so out of the box nothing changes. It is loaded **once at startup**, **before** any glycan database is parsed, and is **shared by both O-glycan and N-glycan searches**.

> **You must restart MetaMorpheus after editing the file by hand** for changes to take effect. A monosaccharide added through the **Create new monosaccharide** dialog is registered immediately and needs no restart.

---

## The file format

Tab-separated, five columns. Only the first three are required.

| # | Column | Required | Notes |
|---|---|---|---|
| 1 | **Name** | ✅ | Unique. Used in *composition-format* glycans, e.g. `HexNAc(2)Hex(5)HexA(1)`. Must not collide with a built-in name or another row. |
| 2 | **SingleCharCode** | ✅ | Exactly one ASCII letter (`A–Z`/`a–z`). Used in *structure-format* glycans, e.g. `(N(H)(U))`. Must not collide with a built-in code or another row. |
| 3 | **MonoisotopicMass** | ✅ | Monoisotopic mass in **Daltons** of the **residue** (the form built into the chain, i.e. after water loss). Decimal, e.g. `176.03209`. |
| 4 | **DiagnosticIonMasses** | ⬜ | Comma-separated observed **m/z** values, no spaces required, e.g. `175.02482,157.01425`. See [Part 3](#part-3--diagnostic--oxonium-ions-column-4). Leave blank for none. |
| 5 | **Description** | ⬜ | Free text, ignored by the parser. Use it to cite a source or note a formula. |

Parsing rules:

- Lines beginning with `#` are **comments** and ignored. A line that Excel re-saved as `"#…` (leading double-quote) is also treated as a comment.
- **Blank lines** are ignored. Indentation before `#` is allowed.
- A **header row** beginning with `Name<TAB>` is recognized and skipped, so you can keep the template header in place.
- A **missing file is not an error** — the loader simply registers no customs.

---

## Built-in monosaccharides (do not redefine)

These names and codes are reserved. A custom row whose **Name** *or* **SingleCharCode** collides with any of these (or with an earlier custom row) is an error and the file fails to load.

| Name | Code | Residue mass (Da) |
|---|---|---|
| Hex | H | 162.05282 |
| HexNAc | N | 203.07937 |
| NeuAc | A | 291.09542 |
| NeuGc | G | 307.09033 |
| Fuc | F | 146.05791 |
| Phospho | P | 79.96633 |
| Sulfo | S | 79.95681 |
| Na | Y | 22.98977 |
| Ac | C | 42.01056 |
| Xylose | X | 150.05282 |
| Kdn | K | 250.06897 |

> `dHex` is a permanent alias of `Fuc` (same code `F`, same slot).

Available single-char codes are therefore any ASCII letter **except** `H N A G F P S Y C X K`.

---

## Part 1 — Define a custom monosaccharide

Add one row per sugar. Example — hexuronic acid (GlcA/GalA), residue mass 176.03209 Da:

```
Name	SingleCharCode	MonoisotopicMass	DiagnosticIonMasses	Description
HexA	U	176.03209		Hexuronic acid (GlcA, GalA)
```

Rules recap:

- **Name** unique, used in composition format.
- **Code** one ASCII letter, unique, used in structure format.
- **Mass** is the **residue** monoisotopic mass (subtract one water from the free-sugar mass).

That is all that is needed to make the sugar *known* to MetaMorpheus. It does **not** yet appear in any search — for that, see Part 2.

---

## Part 2 — Use it in a glycan search

Defining a monosaccharide only adds it to the vocabulary. To search for glycans that contain it, the monosaccharide must appear in a **glycan database file** (`.gdb`) that your Glyco Search task uses.

### Glycan database formats

MetaMorpheus ships two database styles, and your custom sugar is written differently in each:

| Search | Format | Uses | Example with custom `HexA` (name) / `U` (code) |
|---|---|---|---|
| **O-glycan** | **Structure** format | single-char **codes** | `(N(H)(U))` |
| **N-glycan** | **Composition** format | **names** with counts | `HexNAc(2)Hex(5)HexA(1)` |

(For reference, the built-in `OGlycan.gdb` contains lines like `(N(H))`, `(N(A))`, `(N(H)(N))`; the built-in `NGlycan.gdb` contains lines like `HexNAc(2)Hex(1)Fuc(1)`.)

### Steps

1. Make a copy of the relevant database (e.g. `OGlycan.gdb` or `NGlycan.gdb`) and add lines that use your custom sugar in the format above. Keep one glycan per line.
2. In the MetaMorpheus GUI, open the **Glyco Search** task.
3. Point the **O-glycan database** or **N-glycan database** selector at your custom `.gdb`.
4. Make sure `MonosaccharidesCustom.tsv` defines every custom name/code your `.gdb` references (the monosaccharide file is loaded first, so unknown tokens in the `.gdb` would otherwise fail to parse).
5. Run the task.

> In a saved task `.toml`, the relevant keys are `OGlycanDatabasefile` and `NGlycanDatabasefile`.

---

## Part 3 — Diagnostic & oxonium ions (column 4)

Column 4, `DiagnosticIonMasses`, is a comma-separated list of **observed oxonium-ion m/z values** for that monosaccharide. Provide the m/z exactly as you would **see it in the spectrum** — no hydrogen-mass offset is applied.

These ions play **two roles**:

### Role 1 — Scoring (always on)

Whenever a candidate glycan contains the monosaccharide, its diagnostic ions are added to the theoretical diagnostic-ion set, matched against the spectrum, and folded into the **diagnostic-ion score**. This happens regardless of any filter setting and never rejects a candidate.

### Role 2 — Strict filter (when `OxoniumIonFilt` is enabled)

The Glyco Search task has a checkbox labelled **`OxoniumIonFilt`** (on by default; `OxoniumIonFilt = true` in the task `.toml`). When it is enabled, each custom diagnostic ion becomes a **strict gate**:

| Ion observed in spectrum? | Candidate contains the monosaccharide? | Result |
|:---:|:---:|---|
| Yes | Yes | ✅ keep |
| Yes | No | ❌ reject |
| No | Yes | ❌ reject |
| No | No | ✅ keep |

"Observed" means a peak matching the ion's m/z was found within your **product-mass tolerance** **and** rising above **2% of the HexNAc 138.055 diagnostic ion** — the same relative-intensity test the built-in sialic-acid and HexHexNAc rules use. A bare presence test would let a single noise peak at a custom ion's m/z reject every candidate lacking that monosaccharide, for the whole scan. If the 138.055 reference is itself absent from the spectrum there is nothing to measure the ion against, so the strict rule stands down for that spectrum and rejects nothing — a low-energy HCD scan with a strong 204.087 and no 138.055 is ordinary, and neither verdict would be honest there. If **any** listed ion disagrees with the candidate's composition, the candidate is rejected.

This runs **after** the built-in oxonium rules (the 138/144 ratio, the 204 requirement, the 274/292 sialic-acid and 366 HexHexNAc rules), which are unchanged.

**A custom diagnostic ion may not duplicate a built-in oxonium ion** (109.029, 115.040, 126.055, 127.040, 138.055, 144.066, 163.061, 168.066, 186.077, 204.087, 274.093, 290.088, 292.103, 308.098, 366.140, 657.235, 673.230), **or an ion already claimed by another custom monosaccharide**. Values within 0.01 Da of one of those are rejected when the file is loaded, with the offending line named. Listing something ubiquitous such as 204.087 would make the ion "observed" on essentially every glycopeptide spectrum, and the strict rule would then silently reject every candidate lacking your sugar.

> **Turn the filter off** (uncheck `OxoniumIonFilt`) and column-4 ions keep only their scoring role.

### Choosing good ions

Oxonium-ion prevalence depends on **collision energy, stepped-HCD settings, instrument type, precursor charge state, glycan branching, sialylation, and derivatization chemistry**. An ion that is reliably present under one method may be weak under another. **Validate any ion against your own acquisition conditions before relying on it to filter** — otherwise the strict rule may discard real identifications. When in doubt, list ions for scoring but leave `OxoniumIonFilt` off, or start with a single, robust ion.

---

## End-to-end worked example

**Goal:** search for N-glycans containing hexuronic acid (HexA), and require its diagnostic ions when present.

1. **Define the sugar with diagnostic ions** — **Settings** tab → **Create new monosaccharide**, with `175.02482,157.01425` in **Diagnostic Ions (m/z)**; or add the row to `MonosaccharidesCustom.tsv` by hand:

   ```
   Name	SingleCharCode	MonoisotopicMass	DiagnosticIonMasses	Description
   HexA	U	176.03209	175.02482,157.01425	Hexuronic acid (GlcA, GalA)
   ```

2. **Add glycans that use it** to a copy of the N-glycan database, e.g. `NGlycan_HexA.gdb` (composition format, by name):

   ```
   HexNAc(2)Hex(3)HexA(1)
   HexNAc(2)Hex(5)HexA(1)
   ```

3. **Restart MetaMorpheus** (only needed if you edited the TSV by hand; the dialog registers the sugar straight away).

4. In the **Glyco Search** task: set the search type to N-glycan, point the **N-glycan database** at `NGlycan_HexA.gdb`, and leave **`OxoniumIonFilt`** checked.

5. **Run.** Candidate glycopeptides whose spectra show 175.02482 / 157.01425 but lack HexA — or that contain HexA but show neither ion — are rejected; matching candidates keep these ions in their diagnostic-ion score.

---

## Validation & error messages

If a row is malformed, the file fails to load with a `MetaMorpheusException` that names the **file** and the **1-based line number**:

```
Could not parse custom monosaccharide in 'MonosaccharidesCustom.tsv' at line 7: "...".
Expected at least 3 tab-separated columns (Name, SingleCharCode, MonoisotopicMass).
```

Conditions that raise an error:

| Problem | Message contains |
|---|---|
| Fewer than 3 columns | `Expected at least 3 tab-separated columns` |
| Code not exactly one character | `SingleCharCode must be exactly one character` |
| Code is not an ASCII letter | `must be an ASCII letter` |
| Mass not a number | `MonoisotopicMass "…" is not a valid decimal number` |
| A diagnostic-ion value not a number | `DiagnosticIonMasses entry "…" is not a valid decimal number` |
| Name already exists (built-in or earlier row) | `name '…' already exists` |
| Code already exists | `code '…' already exists` |

A **missing file is not an error**.

---

## Tips & caveats

- **Restart required after hand-editing.** The file is read once at startup. The **Create new monosaccharide** dialog registers the sugar for the running session as well as writing it to the file, so it needs no restart.
- **Mass is the residue mass** (free sugar minus one water), not the free monosaccharide mass.
- **Defining ≠ using.** A custom monosaccharide does nothing until a glycan database references it and that database is selected in the task.
- **Shared across O- and N-glycan searches** — names/codes must be globally unique.
- **Excel-safe.** You can edit the file in Excel; comment lines re-saved as `"#…` are still treated as comments. Keep the file **tab-separated**.
- **Default behavior is unchanged.** With no custom rows, search output is identical to stock MetaMorpheus.

---

## FAQ

**GUI dialog or TSV file — is there any difference?**
No. **Create new monosaccharide** writes exactly the row `MonosaccharidesCustom.tsv` expects, and both go through the same validation. The dialog also registers the sugar for the current session; a hand-edit needs a restart.

**Do I have to provide diagnostic ions?**
No. Column 4 is optional. Without it, the monosaccharide still works in databases; it just contributes no extra diagnostic ions.

**Will adding diagnostic ions change results even if I don't filter?**
They can slightly affect the *diagnostic-ion score* (they are matched and scored), but they never reject a candidate unless `OxoniumIonFilt` is enabled.

**Can custom ions tie to a built-in monosaccharide?**
Yes — put a built-in name (e.g. `HexNAc`) in column 1 with ions in column 4. The strict filter then ties those ions to that built-in sugar.

**What m/z do I enter — neutral mass or observed?**
The **observed** singly-charged oxonium m/z, exactly as it appears in the spectrum. No hydrogen offset is applied.

**The strict filter is discarding good hits — what now?**
Your listed ions may be weak under your acquisition method. Remove the unreliable ions, reduce to a single robust ion, or uncheck `OxoniumIonFilt` to keep scoring only.
