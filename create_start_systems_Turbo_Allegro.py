import numpy as np
import math
import itertools


compounds={"Me":{"Cl":["MeMgCl","MgMe2","LiCl","Li","MgCl2","Cl"],
"Br":["MeMgBr","MgMe2","LiCl","Li","MgBr2","Br"]},
"Et":{"Cl":["EtMgCl","MgEt2","LiCl","Li","MgCl2","Cl"],
"Br":["EtMgBr","MgEt2","LiCl","Li","MgBr2","Br"]},
"iPr":{"Cl":["iPrMgCl","MgiPr2","LiCl","Li","MgCl2","Cl"],
"Br":["iPrMgBr","MgiPr2","LiCl","Li","MgBr2","Br"]},
"tBu":{"Cl":["tBuMgCl","MgtBu2","LiCl","Li","MgCl2","Cl"],
"Br":["tBuMgBr","MgtBu2","LiCl","Li","MgBr2","Br"]},
"Ph":{"Cl":["PhMgCl","MgPh2","LiCl","Li","MgCl2","Cl"],
"Br":["PhMgBr","MgPh2","LiCl","Li","MgBr2","Br"]}}


combinations={"Me":{"Cl":{2:list(itertools.combinations(compounds["Me"]["Cl"], 2)),3:list(itertools.combinations(compounds["Me"]["Cl"], 3))},
                    "Br":{2:list(itertools.combinations(compounds["Me"]["Cl"], 2)),3:list(itertools.combinations(compounds["Me"]["Br"], 3))}
                    },
              "Et":{"Cl":{2:list(itertools.combinations(compounds["Et"]["Cl"], 2)),3:list(itertools.combinations(compounds["Et"]["Cl"], 3))},
                    "Br":{2:list(itertools.combinations(compounds["Et"]["Br"], 2)),3:list(itertools.combinations(compounds["Et"]["Br"], 3))}
                    },
              "iPr":{"Cl":{2:list(itertools.combinations(compounds["iPr"]["Cl"], 2)),3:list(itertools.combinations(compounds["iPr"]["Cl"], 3))},
                    "Br":{2:list(itertools.combinations(compounds["iPr"]["Br"], 2)),3:list(itertools.combinations(compounds["iPr"]["Br"], 3))}
                    },
              "tBu":{"Cl":{2:list(itertools.combinations(compounds["tBu"]["Cl"], 2)),3:list(itertools.combinations(compounds["tBu"]["Cl"], 3))},
                    "Br":{2:list(itertools.combinations(compounds["tBu"]["Br"], 2)),3:list(itertools.combinations(compounds["tBu"]["Br"], 3))}
                    },
              "Ph":{"Cl":{2:list(itertools.combinations(compounds["Ph"]["Cl"], 2)),3:list(itertools.combinations(compounds["Ph"]["Cl"], 3))},
                    "Br":{2:list(itertools.combinations(compounds["Ph"]["Br"], 2)),3:list(itertools.combinations(compounds["Ph"]["Br"], 3))}
                    }
              }


print("Combinations:")
for compound in combinations:
    for halide in combinations[compound]:
        for n in combinations[compound][halide]:
            print(f"{compound} {halide} {n}: {len(combinations[compound][halide][n])} combinations")
            for comb in combinations[compound][halide][n]:
                print("  ", " + ".join(comb))
                
                
packmol_input_template_2 = """# Packmol input file for {compound} {halide} {n}-mer
tolerance 2.0
filetype xyz

structure {structure1_file}
  number 1
  inside box 4.0 4.0 4.0 12.0 12.0 12.0
end structure

structure {structure2_file}
  number 1
  inside box 4.0 4.0 4.0 12.0 12.0 12.0
end structure


structure {solvent_file}
  number 39
  inside box 0.0 0.0 0.0 18.0 18.0 18.0
end structure

output {output_file}
"""

packmol_input_template_3 = """# Packmol input file for {compound} {halide} {n}-mer
tolerance 2.0
filetype xyz

structure {structure1_file}
  number 1
  inside box 4.0 4.0 4.0 12.0 12.0 12.0
end structure

structure {structure2_file}
  number 1
  inside box 4.0 4.0 4.0 12.0 12.0 12.0
end structure

structure {structure3_file}
  number 1
  inside box 4.0 4.0 4.0 12.0 12.0 12.0
end structure


structure {solvent_file}
  number 38
  inside box 0.0 0.0 0.0 18.0 18.0 18.0
end structure

output {output_file}
"""
idx=0
for compound in combinations:
    for halide in combinations[compound]:
        for n in combinations[compound][halide]:
            for comb in combinations[compound][halide][n]:
                if n == 2:
                    structure1_file = f"{comb[0]}.xyz"
                    structure2_file = f"{comb[1]}.xyz"
                    solvent_file = "THF.xyz"
                    output_file = f"{idx:05d}_packmol_{compound}_{halide}_{n}_mer_{comb[0]}_{comb[1]}.xyz"
                    with open(f"{idx:05d}_packmol_{compound}_{halide}_{n}_mer_{comb[0]}_{comb[1]}.inp", 'w') as f:
                        f.write(packmol_input_template_2.format(compound=compound, halide=halide, n=n,
                                                                structure1_file=structure1_file,
                                                                structure2_file=structure2_file,
                                                                solvent_file=solvent_file,
                                                                output_file=output_file))
                if n == 3:
                    structure1_file = f"{comb[0]}.xyz"
                    structure2_file = f"{comb[1]}.xyz"
                    structure3_file = f"{comb[2]}.xyz"
                    solvent_file = "THF.xyz"
                    output_file = f"{idx:05d}_packmol_{compound}_{halide}_{n}_mer_{comb[0]}_{comb[1]}_{comb[2]}.xyz"
                    with open(f"{idx:05d}_packmol_{compound}_{halide}_{n}_mer_{comb[0]}_{comb[1]}_{comb[2]}.inp", 'w') as f:
                        f.write(packmol_input_template_3.format(compound=compound, halide=halide, n=n,
                                                                structure1_file=structure1_file,
                                                                structure2_file=structure2_file,
                                                                structure3_file=structure3_file,
                                                                solvent_file=solvent_file,
                                                                output_file=output_file))
                idx += 1
print("Packmol input files created.")
                    
               