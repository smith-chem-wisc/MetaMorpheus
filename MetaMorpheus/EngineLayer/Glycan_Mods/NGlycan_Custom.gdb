# ============================================================================
# CUSTOM N-GLYCAN DATABASE
# ============================================================================
# This file is yours. MetaMorpheus creates it once, if it is missing, and then
# never touches it again -- not on install, not on repair, not on upgrade.
#
# It is offered in the N-glycan database list of the GlycoSearch task alongside
# the databases MetaMorpheus ships, so anything you add below can be searched
# without replacing or editing a shipped file.
#
# Lines starting with '#' are comments and are ignored. Blank lines are ignored
# too, so you may space the file out however you like. Indentation before '#'
# is allowed.
#
# ============================================================================
# FILE FORMAT -- one glycan per line, and ONE format for the whole file
# ============================================================================
# MetaMorpheus infers the format from the FIRST non-comment line and then reads
# the entire file that way. Do not mix the two in one file.
#
# COMPOSITION format (what NGlycan.gdb uses, and what most published N-glycan
# lists are distributed in). Counts only, no branching:
#
#     HexNAc(2)Hex(5)
#     HexNAc(4)Hex(5)Fuc(1)NeuAc(1)
#
#   Names must be ones MetaMorpheus knows:
#
#     Hex        HexNAc     NeuAc      NeuGc      Fuc
#     Phospho    Sulfo      Na         Ac         Xylose      Kdn
#
# STRUCTURE format. Nested parentheses are the tree, and each letter is one
# monosaccharide written with its SINGLE-CHARACTER CODE:
#
#     (N(N(H(H(H)))))
#
#     H  Hex        N  HexNAc     A  NeuAc      G  NeuGc     F  Fuc
#     P  Phospho    S  Sulfo      Y  Na         C  Ac        X  Xylose
#     K  Kdn
#
#   Structure format carries the branching, so the Y-ion series is generated
#   from the tree rather than inferred from the counts. Prefer it when you know
#   the topology.
#
# Every glycan you list is searched on both N-X-S and N-X-T sequons; you do not
# write the motif.
#
# ============================================================================
# MONOSACCHARIDES MetaMorpheus DOES NOT SHIP WITH
# ============================================================================
# Declare them in MonosaccharidesCustom.tsv first -- it sits next to this file
# in the same folder and documents its own format. Once declared, that
# monosaccharide's name and single-character code are legal here, in both
# formats. A name or code this file does not recognize is an error naming the
# line, not a glycan silently given the wrong mass.
#
# ============================================================================
# EXAMPLES (remove the leading '#' to enable, or write your own)
# ============================================================================
# HexNAc(2)Hex(5)
# HexNAc(2)Hex(3)Fuc(1)
# HexNAc(4)Hex(5)Fuc(1)NeuAc(1)
#
# ============================================================================
# YOUR N-GLYCANS (add them below this line)
# ============================================================================
