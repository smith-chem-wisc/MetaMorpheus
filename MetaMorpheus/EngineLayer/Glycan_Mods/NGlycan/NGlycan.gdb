# ============================================================================
# GLYCAN DATABASE FILE (composition format)
# ============================================================================
# Lines starting with '#' are comments and are ignored by the parser.
# Blank lines are also ignored. Indentation before '#' is allowed.
#
# DISCOVERY
# ----------------------------------------------------------------------------
# MetaMorpheus scans the Glycan_Mods/NGlycan/ folder (inside the install
# directory) at startup. Every file in that folder is offered in the
# Glyco Search task's "N-Glycan Database" dropdown by its filename. Pick
# this file there to search against it.
#
# ============================================================================
# LINE FORMAT (composition)
# ============================================================================
# One glycan per data line. Only the text before the first TAB is parsed,
# so additional tab-separated columns (e.g. a short name or a mass) are
# tolerated and ignored.
#
# Syntax:
#     <Monosaccharide>(<count>)<Monosaccharide>(<count>)...
#
# Order does not matter. Monosaccharides with zero count are omitted.
#
# Examples:
#     HexNAc(2)Hex(5)                       core mannose-5
#     HexNAc(2)Hex(5)NeuAc(1)               sialylated
#     HexNAc(4)Hex(5)NeuAc(2)Fuc(1)         biantennary, fucosylated
#
# An alternative STRUCTURE format (example: (N(H(A))(N(H(A))(F))) ) is
# also accepted in glycan files, but EVERY line of a single file must use
# the SAME format. The parser auto-detects format from the first
# non-comment, non-blank line. Do not mix formats within one file.
#
# ============================================================================
# RECOGNIZED MONOSACCHARIDE TOKENS
# ============================================================================
#   Composition name    Single-char code (used in structure format)
#   ------------------  -------------------------------------------
#   Hex                 H
#   HexNAc              N
#   NeuAc               A
#   NeuGc               G
#   Fuc (or dHex)       F
#   Phospho             P
#   Sulfo               S
#   Na                  Y
#   Ac                  C
#   Xylose              X
#   Kdn                 K
#
# Any token outside this set will cause a load error.
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
# ADDING YOUR OWN CUSTOM GLYCANS
# ============================================================================
# Two options:
#
#  (a) Append new composition lines below this header. They will be loaded
#      the next time MetaMorpheus starts.
#
#  (b) Create a new file with any name and a .gdb (or .txt) extension in
#      this same Glycan_Mods/NGlycan/ folder, then restart MetaMorpheus.
#      It will appear in the Glyco Search task dropdown alongside this
#      one. Copy this header into your new file so readers know the format.
#
# You can use '#' on any line to comment your additions.
#
# ============================================================================
# GLYCAN DEFINITIONS
# ============================================================================
HexNAc(1)
HexNAc(1)Fuc(1)
HexNAc(2)
HexNAc(2)Fuc(1)
HexNAc(2)Hex(1)
HexNAc(2)Hex(1)Fuc(1)
HexNAc(2)Hex(10)
HexNAc(2)Hex(11)
HexNAc(2)Hex(12)
HexNAc(2)Hex(2)
HexNAc(2)Hex(2)Fuc(1)
HexNAc(2)Hex(3)
HexNAc(2)Hex(3)Fuc(1)
HexNAc(2)Hex(4)
HexNAc(2)Hex(4)Fuc(1)
HexNAc(2)Hex(5)
HexNAc(2)Hex(5)Fuc(1)
HexNAc(2)Hex(6)
HexNAc(2)Hex(6)Fuc(1)
HexNAc(2)Hex(6)Phospho(1)
HexNAc(2)Hex(7)
HexNAc(2)Hex(7)Fuc(1)
HexNAc(2)Hex(8)
HexNAc(2)Hex(9)
HexNAc(3)Hex(3)
HexNAc(3)Hex(3)Fuc(1)
HexNAc(3)Hex(4)
HexNAc(3)Hex(4)Fuc(1)
HexNAc(3)Hex(4)Fuc(1)NeuAc(1)
HexNAc(3)Hex(4)Fuc(2)
HexNAc(3)Hex(4)Fuc(2)NeuAc(1)
HexNAc(3)Hex(4)NeuAc(1)
HexNAc(3)Hex(5)
HexNAc(3)Hex(5)Fuc(1)
HexNAc(3)Hex(5)Fuc(1)NeuAc(1)
HexNAc(3)Hex(5)Fuc(1)NeuAc(2)
HexNAc(3)Hex(5)NeuAc(1)
HexNAc(3)Hex(6)
HexNAc(3)Hex(6)Fuc(1)
HexNAc(3)Hex(6)Fuc(1)NeuAc(1)
HexNAc(3)Hex(6)NeuAc(1)
HexNAc(4)Hex(3)
HexNAc(4)Hex(3)Fuc(1)
HexNAc(4)Hex(3)NeuAc(1)
HexNAc(4)Hex(4)
HexNAc(4)Hex(4)Fuc(1)
HexNAc(4)Hex(4)Fuc(1)NeuAc(1)
HexNAc(4)Hex(4)Fuc(2)
HexNAc(4)Hex(4)NeuAc(1)
HexNAc(4)Hex(5)
HexNAc(4)Hex(5)Fuc(1)
HexNAc(4)Hex(5)Fuc(1)NeuAc(1)
HexNAc(4)Hex(5)Fuc(1)NeuAc(2)
HexNAc(4)Hex(5)Fuc(2)
HexNAc(4)Hex(5)Fuc(2)NeuAc(1)
HexNAc(4)Hex(5)Fuc(2)NeuAc(2)
HexNAc(4)Hex(5)Fuc(3)NeuAc(1)
HexNAc(4)Hex(5)Fuc(3)NeuAc(2)
HexNAc(4)Hex(5)NeuAc(1)
HexNAc(4)Hex(5)NeuAc(2)
HexNAc(4)Hex(6)
HexNAc(4)Hex(6)Fuc(1)
HexNAc(4)Hex(6)Fuc(1)NeuAc(1)
HexNAc(4)Hex(6)Fuc(2)
HexNAc(4)Hex(6)NeuAc(1)
HexNAc(4)Hex(7)
HexNAc(4)Hex(7)Fuc(1)
HexNAc(4)Hex(7)NeuAc(1)
HexNAc(5)Hex(3)
HexNAc(5)Hex(3)Fuc(1)
HexNAc(5)Hex(3)Fuc(1)NeuAc(1)
HexNAc(5)Hex(3)Fuc(2)
HexNAc(5)Hex(4)
HexNAc(5)Hex(4)Fuc(1)
HexNAc(5)Hex(4)Fuc(1)NeuAc(1)
HexNAc(5)Hex(4)Fuc(1)NeuAc(2)
HexNAc(5)Hex(4)Fuc(2)
HexNAc(5)Hex(4)NeuAc(1)
HexNAc(5)Hex(4)NeuAc(2)
HexNAc(5)Hex(5)
HexNAc(5)Hex(5)Fuc(1)
HexNAc(5)Hex(5)Fuc(1)NeuAc(1)
HexNAc(5)Hex(5)Fuc(1)NeuAc(2)
HexNAc(5)Hex(5)Fuc(2)
HexNAc(5)Hex(5)Fuc(2)NeuAc(1)
HexNAc(5)Hex(5)Fuc(3)
HexNAc(5)Hex(5)NeuAc(1)
HexNAc(5)Hex(5)NeuAc(2)
HexNAc(5)Hex(6)
HexNAc(5)Hex(6)Fuc(1)
HexNAc(5)Hex(6)Fuc(1)NeuAc(1)
HexNAc(5)Hex(6)Fuc(1)NeuAc(2)
HexNAc(5)Hex(6)Fuc(1)NeuAc(3)
HexNAc(5)Hex(6)Fuc(2)
HexNAc(5)Hex(6)Fuc(2)NeuAc(1)
HexNAc(5)Hex(6)Fuc(3)
HexNAc(5)Hex(6)Fuc(3)NeuAc(1)
HexNAc(5)Hex(6)Fuc(3)NeuAc(3)
HexNAc(5)Hex(6)Fuc(4)
HexNAc(5)Hex(6)Fuc(4)NeuAc(1)
HexNAc(5)Hex(6)NeuAc(1)
HexNAc(5)Hex(6)NeuAc(2)
HexNAc(5)Hex(6)NeuAc(3)
HexNAc(5)Hex(7)Fuc(1)NeuAc(1)
HexNAc(5)Hex(7)Fuc(1)NeuAc(2)
HexNAc(5)Hex(8)
HexNAc(5)Hex(8)Fuc(1)
HexNAc(5)Hex(9)Fuc(1)
HexNAc(6)Hex(3)
HexNAc(6)Hex(3)Fuc(1)
HexNAc(6)Hex(3)Fuc(1)NeuAc(1)
HexNAc(6)Hex(3)Fuc(1)NeuAc(2)
HexNAc(6)Hex(3)Fuc(2)
HexNAc(6)Hex(4)
HexNAc(6)Hex(4)Fuc(1)
HexNAc(6)Hex(4)Fuc(2)
HexNAc(6)Hex(4)NeuAc(1)
HexNAc(6)Hex(5)
HexNAc(6)Hex(5)Fuc(1)
HexNAc(6)Hex(5)Fuc(1)NeuAc(1)
HexNAc(6)Hex(5)Fuc(1)NeuAc(2)
HexNAc(6)Hex(5)Fuc(1)NeuAc(3)
HexNAc(6)Hex(5)Fuc(2)
HexNAc(6)Hex(5)Fuc(2)NeuAc(1)
HexNAc(6)Hex(5)Fuc(3)
HexNAc(6)Hex(6)
HexNAc(6)Hex(6)Fuc(1)
HexNAc(6)Hex(6)Fuc(1)NeuAc(1)
HexNAc(6)Hex(6)Fuc(1)NeuAc(2)
HexNAc(6)Hex(6)Fuc(2)
HexNAc(6)Hex(6)Fuc(2)NeuAc(1)
HexNAc(6)Hex(6)Fuc(2)NeuAc(2)
HexNAc(6)Hex(6)NeuAc(1)
HexNAc(6)Hex(6)NeuAc(2)
HexNAc(6)Hex(6)NeuAc(3)
HexNAc(6)Hex(7)
HexNAc(6)Hex(7)Fuc(1)NeuAc(1)
HexNAc(6)Hex(7)Fuc(1)NeuAc(2)
HexNAc(6)Hex(7)Fuc(2)
HexNAc(6)Hex(7)Fuc(2)NeuAc(1)
HexNAc(6)Hex(7)Fuc(3)
HexNAc(6)Hex(7)Fuc(3)NeuAc(1)
HexNAc(6)Hex(7)NeuAc(1)
HexNAc(6)Hex(7)NeuAc(2)
HexNAc(6)Hex(7)NeuAc(3)
HexNAc(6)Hex(7)NeuAc(4)
HexNAc(6)Hex(8)Fuc(1)NeuAc(1)
HexNAc(6)Hex(8)NeuAc(1)
HexNAc(6)Hex(9)
HexNAc(6)Hex(9)Fuc(1)NeuAc(2)
HexNAc(7)Hex(3)
HexNAc(7)Hex(3)Fuc(1)
HexNAc(7)Hex(4)
HexNAc(7)Hex(4)Fuc(1)
HexNAc(7)Hex(6)
HexNAc(7)Hex(6)Fuc(1)
HexNAc(7)Hex(7)
HexNAc(7)Hex(7)Fuc(1)
HexNAc(7)Hex(7)Fuc(1)NeuAc(3)
HexNAc(7)Hex(8)
HexNAc(7)Hex(8)Fuc(1)
HexNAc(7)Hex(8)Fuc(1)NeuAc(1)
HexNAc(7)Hex(8)Fuc(1)NeuAc(4)
HexNAc(7)Hex(8)NeuAc(1)
HexNAc(8)Hex(3)
HexNAc(8)Hex(3)Fuc(1)
HexNAc(8)Hex(4)
HexNAc(8)Hex(5)
HexNAc(8)Hex(5)Fuc(1)
HexNAc(8)Hex(7)
HexNAc(8)Hex(8)
HexNAc(8)Hex(9)
HexNAc(8)Hex(9)Fuc(1)
HexNAc(9)Hex(10)
HexNAc(9)Hex(10)Fuc(1)
HexNAc(9)Hex(3)
HexNAc(9)Hex(3)Fuc(1)
HexNAc(9)Hex(4)
HexNAc(9)Hex(4)Fuc(1)
HexNAc(9)Hex(6)
HexNAc(9)Hex(6)Fuc(1)
HexNAc(9)Hex(9)Fuc(1)
