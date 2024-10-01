Before you begin… 
Read the .md called "STAR, JUM, and DCC using 180 SRA files.md": That document should be read prior to this document. It will show you how to structure the following files: 

```
JUM_A_${prefix1}vs${prefix2}.txt 

Rscript${prefix1}vs${prefix2}.txt

JUM_B_${prefix1}vs${prefix2}.txt 

JUM_C_${prefix1}vs${prefix2}.txt 

JUM_A_${prefix1}vs${prefix2}.txt 

DCC.txt 
```
 

These loops always contain a vector called “prefixes.” Modify your prefixes to match what you want to call your replicates, cell types, and/or treatments. This will be unique to your experiment and overall scientific question. I am using 46 different cell types with varying numbers of replicates, totaling 180 unique biological replicates. 


### JUM_A_cp_and_dos2unix_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem  

module load gcc dos2unix 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Copy JUM_A.txt for files with the given prefixes 

for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir"  

      # Copy the file to the working directory 
      cp "JUM_A_${prefix1}vs${prefix2}.txt" "$working_dir" 

      # Convert the file to Unix format (if necessary):
      dos2unix "$working_dir/JUM_A_${prefix1}vs${prefix2}.txt" 
    fi 
  done 
done 
```

### JUM_A_sbatch_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
 
# Run JUM_A for files with the given prefixes 

for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 

      # Submit the job to Slurm 
      sbatch "$working_dir/JUM_A_${prefix1}vs${prefix2}.txt" 
    fi 
  done 
done 
```

### JUM_B_cp_and_dos2unix_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem  

module load gcc dos2unix 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Copy JUM_B for files with the given prefixes 

for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 
      jum_diff_dir="$source_dir/$sub_dir/JUM_diff" 

      # Copy the files to the working directory 
      cp "experiment_design${prefix1}vs${prefix2}.txt" "$jum_diff_dir" 
      cp "Rscript${prefix1}vs${prefix2}.txt" "$jum_diff_dir" 
      cp "JUM_B_${prefix1}vs${prefix2}.txt" "$jum_diff_dir" 

      # Convert the file to Unix format (if necessary) 
      dos2unix "$jum_diff_dir/experiment_design${prefix1}vs${prefix2}.txt" 
      dos2unix "$jum_diff_dir/Rscript${prefix1}vs${prefix2}.txt" 
      dos2unix "$jum_diff_dir/JUM_B_${prefix1}vs${prefix2}.txt"
    fi 
  done 
done 
```

### Rscript_sbatch_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Go to JUM_diff for files with the given prefixes 
for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 
      jum_diff_dir="$source_dir/$sub_dir/JUM_diff"  

      # Change to the working directory 
      cd "$jum_diff_dir/" 

      # sbatch Rscript step 
      sbatch "Rscript${prefix1}vs${prefix2}.txt" 
    fi 
  done 
done 
```

### JUM_B_sbatch_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Go to JUM_diff for files with the given prefixes 
for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 
      jum_diff_dir="$source_dir/$sub_dir/JUM_diff" 

      # Change to the working directory 
      cd "$jum_diff_dir/" 
 
      # sbatch Rscript step 
      sbatch "JUM_B_${prefix1}vs${prefix2}.txt" 
    fi 
  done 
done
```

### JUM_C_cp_and_dos2unix_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

module load gcc dos2unix 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Copy JUM_C for files with the given prefixes 
for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 
      jum_diff_dir="$source_dir/$sub_dir/JUM_diff" 
      final_jum_dir="$jum_diff_dir/FINAL_JUM_OUTPUT_pvalue_1" 
 
      # Copy the files to the working directory 
      cp "JUM_C_${prefix1}vs${prefix2}.txt" "$final_jum_dir" 
      cp "refFlat.txt" "$final_jum_dir" 

      # Convert the file to Unix format (if necessary) 
      dos2unix "$final_jum_dir/JUM_C_${prefix1}vs${prefix2}.txt" 
      dos2unix "$final_jum_dir/refFlat.txt"
    fi 
  done 
done 
```

### JUM_C_sbatch_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 

prefixes=("ADL" "AFD" "AIM" "AIN" "AIY" "ASEL" "ASER" "ASG" "ASI" "ASK" "AVA" "AVE" "AVG" "AVH" "AVK" "AVL" "AVM" "AWA" "AWB" "AWC" "BAG" "CAN" "DA" "DD" "DVC" "I5" "IL1" "IL2" "LUA" "NSM" "OLL" "OLQ" "PHA" "PVC" "PVD" "PVM" "RIA" "RIC" "RIM" "RIS" "RMD" "SMB" "SMD" "VB" "VC" "VD") 

# Go to FINAL_JUM_OUTPUT_pvalue_1 for files with the given prefixes 
for prefix1 in "${prefixes[@]}"; do 
  for prefix2 in "${prefixes[@]}"; do 
    if [[ "$prefix1" != "$prefix2" ]]; then 
      source_dir="/lustre/work/client/group/norrislab/Zach/STAR_results" 
      sub_dir="${prefix1}vs${prefix2}" 
      working_dir="$source_dir/$sub_dir" 
      final_jum_dir="$source_dir/$sub_dir/JUM_diff/FINAL_JUM_OUTPUT_pvalue_1" 

      # Change to the working directory 
      cd "$final_jum_dir/" 
      
      # sbatch JUM_C step 
      sbatch "JUM_C_${prefix1}vs${prefix2}.txt" 
    fi 
  done 
done
```

### DCC_sbatch_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p highmem 
#SBATCH --mem=200G  

prefixes=("ADL1" "ADL2" "ADL3" "ADL4" "AFD1" "AFD2" "AFD3" "AFD4" "AFD5" "AIM1" "AIM2" "AIM3" "AIM4" "AIN1" "AIN2" "AIN3" "AIN4" "AIN5" "AIN6" "AIY1" "AIY2" "AIY3" "ASEL1" "ASEL2" "ASEL3" "ASER1" "ASER2" "ASER3" "ASER4" "ASG1" "ASG2" "ASG3" "ASG4" "ASI1" "ASI2" "ASI3" "ASI4" "ASK1" "ASK2" "ASK3" "ASK4" "AVA1" "AVA2" "AVA3" "AVA4" "AVA5" "AVA6" "AVE1" "AVE2" "AVE3" "AVG1" "AVG2" "AVG3" "AVH1" "AVH2" "AVH3" "AVH4" "AVK1" "AVK2" "AVK3" "AVK4" "AVL1" "AVL2" "AVL3" "AVM1" "AVM2" "AVM3" "AWA1" "AWA2" "AWA3" "AWA4" "AWB1" "AWB2" "AWB3" "AWB4" "AWB5" "AWC1" "AWC2" "AWC3" "AWC4" "BAG1" "BAG2" "BAG3" "BAG4" "CAN1" "CAN2" "CAN3" "DA1" "DA2" "DA3" "DA4" "DD1" "DD2" "DD3" "DD4" "DVC1" "DVC2" "DVC3" "DVC4" "I51" "I52" "I53" "I54" "IL11" "IL12" "IL13" "IL21" "IL22" "IL23" "IL24" "LUA1" "LUA2" "LUA3" "LUA4" "NSM1" "NSM2" "NSM3" "OLL1" "OLL2" "OLQ1" "OLQ2" "OLQ3" "PHA1" "PHA2" "PHA3" "PHA4" "PVC1" "PVC2" "PVC3" "PVC4" "PVC5" "PVD1" "PVD2" "PVM1" "PVM2" "RIA1" "RIA2" "RIA3" "RIA4" "RIA5" "RIA6" "RIC1" "RIC2" "RIC3" "RIC4" "RIM1" "RIM2" "RIM3" "RIM4" "RIS1" "RIS2" "RIS3" "RMD1" "RMD2" "RMD3" "RMD4" "RMD5" "SMB1" "SMB2" "SMB3" "SMB4" "SMB5" "SMD1" "SMD2" "SMD3" "SMD4" "VB1" "VB2" "VB3" "VB4" "VC1" "VC2" "VC3" "VC4" "VC5" "VC6" "VD1" "VD2" "VD3" "VD4") 

for prefix in "${prefixes[@]}"; do 
  samplesheet_file="/lustre/work/client/group/norrislab/Zach/${prefix}Chimeric.out.junction" 
  mate1_file="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/${prefix}_1Chimeric.out.junction" 
  mate2_file="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/${prefix}_2Chimeric.out.junction" 
  output_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/${prefix}_DCC_output" 

  mkdir -p "$output_dir"  # Create the output directory if it doesn't exist 

  DCC "$samplesheet_file" -mt1 "$mate1_file" -mt2 "$mate2_file" -D -an /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.107.gtf -B /lustre/work/client/group/norrislab/Zach/${prefix}Aligned.sortedByCoord.out.bam -Pi -N -F -M -Nr 5 1 -fg -G -A /work/users/zwolfe/genome/Caenorhabditis_elegans.WBcel235.dna.toplevel.fa 

  # Move the output files to the output directory 

  mv output.err "$output_dir" 
  mv output.out "$output_dir" 
  mv LinearCount "$output_dir" 
  mv CircSkipJunctions "$output_dir" 
  mv CircCoordinates "$output_dir" 
  mv CircRNACount "$output_dir" 
  mv DCC-*.log "$output_dir" 
done 
```

### DCC_output_file_rename_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=20G 

prefixes=("ADL1" "ADL2" "ADL3" "ADL4" "AFD1" "AFD2" "AFD3" "AFD4" "AFD5" "AIM1" "AIM2" "AIM3" "AIM4" "AIN1" "AIN2" "AIN3" "AIN4" "AIN5" "AIN6" "AIY1" "AIY2" "AIY3" "ASEL1" "ASEL2" "ASEL3" "ASER1" "ASER2" "ASER3" "ASER4" "ASG1" "ASG2" "ASG3" "ASG4" "ASI1" "ASI2" "ASI3" "ASI4" "ASK1" "ASK2" "ASK3" "ASK4" "AVA1" "AVA2" "AVA3" "AVA4" "AVA5" "AVA6" "AVE1" "AVE2" "AVE3" "AVG1" "AVG2" "AVG3" "AVH1" "AVH2" "AVH3" "AVH4" "AVK1" "AVK2" "AVK3" "AVK4" "AVL1" "AVL2" "AVL3" "AVM1" "AVM2" "AVM3" "AWA1" "AWA2" "AWA3" "AWA4" "AWB1" "AWB2" "AWB3" "AWB4" "AWB5" "AWC1" "AWC2" "AWC3" "AWC4" "BAG1" "BAG2" "BAG3" "BAG4" "CAN1" "CAN2" "CAN3" "DA1" "DA2" "DA3" "DA4" "DD1" "DD2" "DD3" "DD4" "DVC1" "DVC2" "DVC3" "DVC4" "I51" "I52" "I53" "I54" "IL11" "IL12" "IL13" "IL21" "IL22" "IL23" "IL24" "LUA1" "LUA2" "LUA3" "LUA4" "NSM1" "NSM2" "NSM3" "OLL1" "OLL2" "OLQ1" "OLQ2" "OLQ3" "PHA1" "PHA2" "PHA3" "PHA4" "PVC1" "PVC2" "PVC3" "PVC4" "PVC5" "PVD1" "PVD2" "PVM1" "PVM2" "RIA1" "RIA2" "RIA3" "RIA4" "RIA5" "RIA6" "RIC1" "RIC2" "RIC3" "RIC4" "RIM1" "RIM2" "RIM3" "RIM4" "RIS1" "RIS2" "RIS3" "RMD1" "RMD2" "RMD3" "RMD4" "RMD5" "SMB1" "SMB2" "SMB3" "SMB4" "SMB5" "SMD1" "SMD2" "SMD3" "SMD4" "VB1" "VB2" "VB3" "VB4" "VC1" "VC2" "VC3" "VC4" "VC5" "VC6" "VD1" "VD2" "VD3" "VD4") 

for prefix in "${prefixes[@]}"; do 
  output_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/${prefix}_DCC_output" 

  # Rename specific output files 
  mv "$output_dir/CircSkipJunctions" "$output_dir/${prefix}CircSkipJunctions" 
  mv "$output_dir/LinearCount" "$output_dir/${prefix}LinearCount" 
  mv "$output_dir/CircRNACount" "$output_dir/${prefix}CircRNACount" 
  mv "$output_dir/CircCoordinates" "$output_dir/${prefix}CircCoordinates" 
done 
```

### DCC_output_file_cp_loop: 

```
#!/bin/bash 
#SBATCH -N 1 
#SBATCH -o output.out 
#SBATCH -e output.err 
#SBATCH -p standard-s 
#SBATCH --mem=2G 

output_files_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/all_DCC_output_files" 

prefixes=("ADL1" "ADL2" "ADL3" "ADL4" "AFD1" "AFD2" "AFD3" "AFD4" "AFD5" "AIM1" "AIM2" "AIM3" "AIM4" "AIN1" "AIN2" "AIN3" "AIN4" "AIN5" "AIN6" "AIY1" "AIY2" "AIY3" "ASEL1" "ASEL2" "ASEL3" "ASER1" "ASER2" "ASER3" "ASER4" "ASG1" "ASG2" "ASG3" "ASG4" "ASI1" "ASI2" "ASI3" "ASI4" "ASK1" "ASK2" "ASK3" "ASK4" "AVA1" "AVA2" "AVA3" "AVA4" "AVA5" "AVA6" "AVE1" "AVE2" "AVE3" "AVG1" "AVG2" "AVG3" "AVH1" "AVH2" "AVH3" "AVH4" "AVK1" "AVK2" "AVK3" "AVK4" "AVL1" "AVL2" "AVL3" "AVM1" "AVM2" "AVM3" "AWA1" "AWA2" "AWA3" "AWA4" "AWB1" "AWB2" "AWB3" "AWB4" "AWB5" "AWC1" "AWC2" "AWC3" "AWC4" "BAG1" "BAG2" "BAG3" "BAG4" "CAN1" "CAN2" "CAN3" "DA1" "DA2" "DA3" "DA4" "DD1" "DD2" "DD3" "DD4" "DVC1" "DVC2" "DVC3" "DVC4" "I51" "I52" "I53" "I54" "IL11" "IL12" "IL13" "IL21" "IL22" "IL23" "IL24" "LUA1" "LUA2" "LUA3" "LUA4" "NSM1" "NSM2" "NSM3" "OLL1" "OLL2" "OLQ1" "OLQ2" "OLQ3" "PHA1" "PHA2" "PHA3" "PHA4" "PVC1" "PVC2" "PVC3" "PVC4" "PVC5" "PVD1" "PVD2" "PVM1" "PVM2" "RIA1" "RIA2" "RIA3" "RIA4" "RIA5" "RIA6" "RIC1" "RIC2" "RIC3" "RIC4" "RIM1" "RIM2" "RIM3" "RIM4" "RIS1" "RIS2" "RIS3" "RMD1" "RMD2" "RMD3" "RMD4" "RMD5" "SMB1" "SMB2" "SMB3" "SMB4" "SMB5" "SMD1" "SMD2" "SMD3" "SMD4" "VB1" "VB2" "VB3" "VB4" "VC1" "VC2" "VC3" "VC4" "VC5" "VC6" "VD1" "VD2" "VD3" "VD4") 
 
for prefix in "${prefixes[@]}"; do 
  output_dir="/lustre/work/client/group/norrislab/Zach/STAR_results_single_mates/DCC/${prefix}_DCC_output" 

  # Copy files to the desired filepath 
 cp "$output_dir/${prefix}CircSkipJunctions" "$output_files_dir" 
 cp "$output_dir/${prefix}LinearCount" "$output_files_dir" 
 cp "$output_dir/${prefix}CircRNACount" "$output_files_dir" 
 cp "$output_dir/${prefix}CircCoordinates" "$output_files_dir" 
done 
```

Copyright 2024 The Regents of the University of California

All Rights Reserved

Created by Zachery Wolfe

Department of Biochemistry

This file is part of Differential Expression in C. elegans. \
Differential Expression in C. elegans is a free pipeline: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. \
Differential Expression in C. elegans is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. \
You should have received a copy of the GNU General Public License along with Differential Expression in C. elegans. If not, see <https://www.gnu.org/licenses/>.
