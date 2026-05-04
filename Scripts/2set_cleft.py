import subprocess
import textwrap

# =============================
# USER INPUTS
# =============================

gro_file = "em.gro"
xtc_file = "md_noPBC.xtc"

# ---- Set 1 residues ----
res1_set1 = 194
res2_set1 = 249

# ---- Set 2 residues ----
res1_set2 = 198
res2_set2 = 256

output_file = "helix_opening_distances.txt"

vmd_path = "vmd"

# =============================
# TCL SCRIPT
# =============================

tcl_script = f"""
# =============================
# Load trajectory
# =============================

mol new {gro_file} type gro waitfor all
mol addfile {xtc_file} type xtc waitfor all

set molid top
set nframes [molinfo $molid get numframes]

# =============================
# Alignment (backbone)
# =============================

set ref_sel [atomselect $molid "backbone" frame 0]
set mob_sel [atomselect $molid "backbone"]

# =============================
# Residue selections
# =============================

# --- Set 1 ---
set sel1_set1 [atomselect $molid "resid {res1_set1}"]
set sel2_set1 [atomselect $molid "resid {res2_set1}"]

# --- Set 2 ---
set sel1_set2 [atomselect $molid "resid {res1_set2}"]
set sel2_set2 [atomselect $molid "resid {res2_set2}"]

# =============================
# Output file
# =============================

set outfile [open "{output_file}" w]
puts $outfile "#Frame  Dist_Set1(Angstrom)  Dist_Set2(Angstrom)"

# =============================
# Loop over frames
# =============================

for {{set i 0}} {{$i < $nframes}} {{incr i}} {{

    animate goto $i

    # ---- Align ----
    $mob_sel frame $i
    $mob_sel update
    set trans_mat [measure fit $mob_sel $ref_sel]
    $mob_sel move $trans_mat

    # ---- Update selections ----
    foreach sel [list $sel1_set1 $sel2_set1 $sel1_set2 $sel2_set2] {{
        $sel frame $i
        $sel update
    }}

    # =============================
    # ---- SET 1 ----
    # =============================

    set com1_set1 [measure center $sel1_set1 weight mass]
    set com2_set1 [measure center $sel2_set1 weight mass]
    set dist_set1 [veclength [vecsub $com1_set1 $com2_set1]]

    # =============================
    # ---- SET 2 ----
    # =============================

    set com1_set2 [measure center $sel1_set2 weight mass]
    set com2_set2 [measure center $sel2_set2 weight mass]
    set dist_set2 [veclength [vecsub $com1_set2 $com2_set2]]

    # ---- Write ----
    puts $outfile "$i  $dist_set1  $dist_set2"
}}

# =============================
# Cleanup
# =============================

close $outfile

foreach sel [list $ref_sel $mob_sel $sel1_set1 $sel2_set1 $sel1_set2 $sel2_set2] {{
    $sel delete
}}

quit
"""

# =============================
# Write TCL
# =============================

with open("helix_opening_analysis.tcl", "w") as f:
    f.write(textwrap.dedent(tcl_script))

# =============================
# Run VMD
# =============================

subprocess.run([
    vmd_path,
    "-dispdev", "text",
    "-e", "helix_opening_analysis.tcl"
])

print("Done. Output:", output_file)
