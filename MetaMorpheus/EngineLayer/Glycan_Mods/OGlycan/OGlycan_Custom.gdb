# ============================================================================
# CUSTOM O-GLYCAN DATABASE (composition format)
# ============================================================================
# This file is a TEMPLATE for users to add their own O-glycan compositions
# without modifying the built-in OGlycan.gdb. It ships with no glycan rows
# of its own; only the example rows below (commented out) demonstrate the
# format.
#
# Lines starting with '#' are comments and are ignored by the parser.
# Blank lines are also ignored. Indentation before '#' is allowed.
#
# DISCOVERY
# ----------------------------------------------------------------------------
# MetaMorpheus scans the Glycan_Mods/OGlycan/ folder (inside the install
# directory) at startup. Every file in that folder appears in the Glyco
# Search task's "O-Glycan Database" dropdown by its filename. Select
# "OGlycan_Custom.gdb" there to search against the entries you add below.
#
# As long as this file has no uncommented data rows it is harmless to
# leave it selected: the search will simply find no glycans in it.
#
# ============================================================================
# HOW TO ADD A GLYCAN
# ============================================================================
# Add one composition per line, below the GLYCAN DEFINITIONS banner.
#
# Syntax:
#     <Monosaccharide>(<count>)<Monosaccharide>(<count>)...
#
# Order does not matter. Monosaccharides with zero count are omitted.
#
# Recognized monosaccharide names:
#     Hex, HexNAc, NeuAc, NeuGc, Fuc (or dHex), Phospho, Sulfo,
#     Na, Ac, Xylose, Kdn
#
# Any name outside this set will cause a load error.
#
# An alternative STRUCTURE format (example: (N(H(A))(N(H(A))(F))) ) is
# also accepted, but EVERY line of a single file must use the SAME format.
# Do not mix formats within one file. The built-in OGlycan.gdb is in
# structure format; this template uses composition because it is easier
# for non-experts to write.
#
# ============================================================================
# AUTOMATIC EXPANSION (you do not write these out yourself)
# ============================================================================
# For each composition listed below, MetaMorpheus automatically:
#   - Generates ONE entry targeting Ser (S) and a SECOND entry targeting
#     Thr (T). Do not duplicate rows for the two motifs.
#   - Generates neutral-loss masses, Y-ions, and diagnostic oxonium ions
#     from the composition. You do NOT need to enumerate fragments.
#
# ============================================================================
# EXAMPLES (remove the leading '#' to enable any of these)
# ============================================================================
# HexNAc(1)                             single GalNAc (Tn antigen)
# HexNAc(1)Hex(1)                       core-1 (T antigen)
# HexNAc(1)Hex(1)NeuAc(1)               sialyl-core-1
# HexNAc(2)Hex(1)                       core-2
#
# ============================================================================
# GLYCAN DEFINITIONS (add your compositions below this line)
# ============================================================================
