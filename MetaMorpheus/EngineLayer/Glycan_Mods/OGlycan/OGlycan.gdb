# ============================================================================
# GLYCAN DATABASE FILE (structure format)
# ============================================================================
# Lines starting with '#' are comments and are ignored by the parser.
# Blank lines are also ignored. Indentation before '#' is allowed.
#
# DISCOVERY
# ----------------------------------------------------------------------------
# MetaMorpheus scans the Glycan_Mods/OGlycan/ folder (inside the install
# directory) at startup. Every file in that folder is offered in the
# Glyco Search task's "O-Glycan Database" dropdown by its filename. Pick
# this file there to search against it.
#
# ============================================================================
# LINE FORMAT (structure)
# ============================================================================
# One glycan per data line. The glycan is encoded as a parenthesized
# linkage tree: each single-character code is a monosaccharide, and
# parenthesized groups attached to it are its branches.
#
# Syntax:
#     (<root>(<branch>)(<branch>)...)
#
# Examples:
#     (N)                              single GlcNAc
#     (N(H))                           GlcNAc-Hex linear
#     (N(H)(N))                        GlcNAc branching to Hex and GlcNAc
#     (N(H(A)))                        GlcNAc-Hex-NeuAc linear
#     (N(H(A))(N(H(A))(F)))            biantennary, sialylated, fucosylated
#
# An alternative COMPOSITION format (example: HexNAc(2)Hex(5)NeuAc(1) ) is
# also accepted in glycan files, but EVERY line of a single file must use
# the SAME format. The parser auto-detects format from the first
# non-comment, non-blank line. Do not mix formats within one file.
#
# ============================================================================
# RECOGNIZED MONOSACCHARIDE CODES (structure format uses single chars)
# ============================================================================
#   Code    Monosaccharide      Composition-format name
#   ----    -----------------   -----------------------
#   H       Hexose              Hex
#   N       N-Acetylhexosamine  HexNAc
#   A       N-Acetylneuraminate NeuAc
#   G       N-Glycolylneuraminate NeuGc
#   F       Fucose              Fuc (or dHex)
#   P       Phosphate           Phospho
#   S       Sulfate             Sulfo
#   Y       Sodium              Na
#   C       Acetyl              Ac
#   X       Xylose              Xylose
#   K       Kdn                 Kdn
#
# Any character outside this set will cause a load error or a wrong mass.
#
# ============================================================================
# AUTOMATIC EXPANSION (you do not write these out yourself)
# ============================================================================
# For each structure listed below, MetaMorpheus automatically:
#   - Generates ONE entry targeting Ser (S) and a SECOND entry targeting
#     Thr (T). Do not duplicate rows for the two motifs.
#   - Enumerates all sub-tree fragments to produce Y-ion neutral losses,
#     plus diagnostic oxonium ions from the composition. You do NOT need
#     to enumerate fragments.
#
# ============================================================================
# ADDING YOUR OWN CUSTOM GLYCANS
# ============================================================================
# Two options:
#
#  (a) Append new structure lines below this header. They will be loaded
#      the next time MetaMorpheus starts.
#
#  (b) Create a new file with any name and a .gdb (or .txt) extension in
#      this same Glycan_Mods/OGlycan/ folder, then restart MetaMorpheus.
#      It will appear in the Glyco Search task dropdown alongside this
#      one. Copy this header into your new file so readers know the format.
#
# You can use '#' on any line to comment your additions.
#
# ============================================================================
# GLYCAN DEFINITIONS
# ============================================================================
(N)
(N(H))
(N(A))
(N(H)(N))
(N(H(A)))
(N(H)(N(H)))
(N(H(A))(N))
(N(H(A))(A))
(N(H(A))(N(H)))
(N(H)(N(H(A))(F)))
(N(H(A))(N(H(A))))
(N(H(A))(N(H(A))(F)))