import subprocess
import tempfile

print("=== Ligand COM RMSD Calculator + Gnuplot Plotter (LIG301 vs all mutants) ===")

holo_pdb = "holo_fixed.pdb"
gro_file = "em.gro"
xtc_file = "md_noPBC.xtc"

# Mutant ligand residue IDs
mutant_resids = [264, 265, 266, 267, 268, 269]

# -----------------------------------------
# (1) Build TCL script for VMD computation
# -----------------------------------------

tcl_script = """\
# Load holo (mol 0)
mol new "{holo}" type pdb waitfor all

# Load mutant GRO + XTC (mol 1)
mol new "{gro}" type gro waitfor all
mol addfile "{xtc}" type xtc waitfor all

# Align trajectory (protein CA)
set sel_ref [atomselect 0 "protein and name CA"]
set sel_mut [atomselect 1 "protein and name CA"]
set all_mut [atomselect 1 "all"]
set nf [molinfo 1 get numframes]

for {{set i 0}} {{$i < $nf}} {{incr i}} {{
    $sel_mut frame $i
    set T [measure fit $sel_mut $sel_ref]
    $all_mut frame $i
    $all_mut move $T
}}

# Reference ligand (holo LIG301)
set lig_ref [atomselect 0 "resname LIG and resid 301 and noh"]
if {{[$lig_ref num] == 0}} {{
    puts "ERROR: Reference LIG301 not found!"
    quit
}}

set com_ref [measure center $lig_ref weight mass]

# Mutant ligand residue list
set mutant_resids {{{resids}}}

# Loop over each mutant ligand
foreach rid $mutant_resids {{

    set lig_mut [atomselect 1 "resname LIG and resid $rid and noh"]

    if {{[$lig_mut num] == 0}} {{
        puts "WARNING: LIG $rid not found, skipping."
        continue
    }}

    set outfile [format "COM_RMSD_vs_LIG%d.dat" $rid]
    set fp [open $outfile w]

    for {{set i 0}} {{$i < $nf}} {{incr i}} {{
        $lig_mut frame $i
        set com_mut [measure center $lig_mut weight mass]
        set d [vecdist $com_mut $com_ref]
        puts $fp "$i $d"
    }}

    close $fp
    puts "Finished LIG $rid"
}}

puts "All COM RMSD calculations completed."
quit
""".format(
    holo=holo_pdb,
    gro=gro_file,
    xtc=xtc_file,
    resids=" ".join(map(str, mutant_resids))
)

# Write TCL file
with tempfile.NamedTemporaryFile(mode="w", delete=False, suffix=".tcl") as tfile:
    tfile.write(tcl_script)
    tcl_path = tfile.name

print("\nRunning VMD...")
subprocess.run(["vmd", "-dispdev", "text", "-e", tcl_path], check=True)

print("\nCOM RMSD files created:")
for rid in mutant_resids:
    print(f"  COM_RMSD_vs_LIG{rid}.dat")

# -----------------------------------------
# (2) Generate GNUplot scripts (Å → nm)
# -----------------------------------------

for rid in mutant_resids:

    gnu_script = f"""\
set terminal pngcairo size 1600,1000 font "Helvetica,24"
set output "COM_RMSD_LIG{rid}_nm.png"

# Axes ranges
set xrange [0:250]
set yrange [0:8]

# Labels (BOLD)
set xlabel "Time (ns)" font ",28"
set ylabel "COM RMSD (nm)" font ",28"
set title "Ligand COM RMSD vs Time (LIG301 vs LIG{rid})" font ",30"

# Ticks (bold & thick)
set tics font ",22"
set border linewidth 2

# Grid (subtle)
set grid lw 1 lc rgb "#cccccc"

# Legend
set key top right font ",22" spacing 1.2

# Line style (bold)
set style line 1 lw 4 lc rgb "#1f77b4"

plot \\
    "COM_RMSD_vs_LIG{rid}.dat" using ($1*0.1):($2/10.0) \\
    with lines ls 1 title "LIG {rid}"
"""

    gnu_file = f"plot_com_rmsd_LIG{rid}.gnu"
    with open(gnu_file, "w") as gp:
        gp.write(gnu_script)

    print(f"\nRunning gnuplot for LIG {rid}...")
    subprocess.run(["gnuplot", gnu_file], check=True)
