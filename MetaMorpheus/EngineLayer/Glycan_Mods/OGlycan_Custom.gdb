# ============================================================================
# CUSTOM O-GLYCAN DATABASE
# ============================================================================
# This file is yours. MetaMorpheus creates it once, if it is missing, and then
# never touches it again -- not on install, not on repair, not on upgrade.
#
# It is offered in the O-glycan database list of the GlycoSearch task alongside
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
# STRUCTURE format (what OGlycan.gdb uses; preferred for O-glycans, because the
# branching is what localization and child-ion generation work from):
#
#     (N(H(A)))
#
#   Nested parentheses are the tree. Each letter is one monosaccharide, written
#   with its SINGLE-CHARACTER CODE:
#
#     H  Hex        N  HexNAc     A  NeuAc      G  NeuGc     F  Fuc
#     P  Phospho    S  Sulfo      Y  Na         C  Ac        X  Xylose
#     K  Kdn
#
# COMPOSITION format (what the shipped .txt databases use). Counts only, no
# branching:
#
#     HexNAc(1)Hex(1)NeuAc(1)
#
# Every glycan you list is searched on both serine and threonine; you do not
# write the residue.
#
# ============================================================================
# MONOSACCHARIDES MetaMorpheus DOES NOT SHIP WITH
# ============================================================================
# Declare them in MonosaccharidesCustom.tsv first -- it sits next to this file
# in the same folder and documents its own format. Once declared, that
# monosaccharide's name and single-character code are legal here, in both
# formats. A code this file does not recognize is an error naming the line, not
# a glycan silently given the wrong mass.
#
# ============================================================================
# EXAMPLES (remove the leading '#' to enable, or write your own)
# ============================================================================
# (N)
# (N(H))
# (N(H(A)))
# (N(H(A))(A))
#
# ============================================================================
# YOUR O-GLYCANS (add them below this line)
# ============================================================================
