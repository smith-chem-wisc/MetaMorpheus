# ============================================================================
# CUSTOM N-GLYCAN DATABASE (composition format)
# ============================================================================
# This file is a TEMPLATE for users to add their own N-glycan compositions
# without modifying the built-in NGlycan.gdb. It ships with no glycan rows
# of its own; only the example rows below (commented out) demonstrate the
# format.
#
# Lines starting with '#' are comments and are ignored by the parser.
# Blank lines are also ignored. Indentation before '#' is allowed.
#
# DISCOVERY
# ----------------------------------------------------------------------------
# MetaMorpheus scans the Glycan_Mods/NGlycan/ folder (inside the install
# directory) at startup. Every file in that folder appears in the Glyco
# Search task's "N-Glycan Database" dropdown by its filename. Select
# "NGlycan_Custom.gdb" there to search against the entries you add below.
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
# Do not mix formats within one file.
#
# ============================================================================
# AUTOMATIC EXPANSION (you do not write these out yourself)
# ============================================================================
# For each composition listed below, MetaMorpheus automatically:
#   - Generates ONE entry targeting motif Asn-X-Ser (Nxs) and a SECOND
#     entry targeting Asn-X-Thr (Nxt). Do not duplicate rows for the two
#     motifs.
#   - Generates neutral-loss masses, Y-ions, and diagnostic oxonium ions
#     from the composition. You do NOT need to enumerate fragments.
#
# ============================================================================
# EXAMPLES (remove the leading '#' to enable any of these)
# ============================================================================
# HexNAc(2)Hex(5)                       core mannose-5
# HexNAc(2)Hex(5)NeuAc(1)               sialylated
# HexNAc(4)Hex(5)NeuAc(2)Fuc(1)         biantennary, fucosylated
#
# ============================================================================
# GLYCAN DEFINITIONS (add your compositions below this line)
# ============================================================================
