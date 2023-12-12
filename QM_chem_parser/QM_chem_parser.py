import json
import argparse
import cclib
import glob
import re
import sys
import logging
import os.path
from cclib.parser import ccData
from cclib.io import ccopen
from cclib.io import ccwrite
from cclib.scripts.ccget import ccget
from cclib.parser import logfileparser
import io
from tech_parse import parse_tech
import pandas as pd

class NumpyAwareJSONEncoder(json.JSONEncoder):
    """A encoder for numpy.ndarray's obtained from the cclib attributes.
       For all other types the json default encoder is called.
       Do Not rename the 'default' method as it is required to be implemented
       by any subclass of the json.JSONEncoder
    """
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            if obj.ndim == 1:
                nan_list = obj.tolist()
                return [None if np.isnan(x) else x for x in nan_list]
            else:
                return [self.default(obj[i]) for i in range(obj.shape[0])]
        return json.JSONEncoder.default(self, obj)


def parsed_data():
    pdata={
  "aonames": {
    "desc": "atomic orbital names",
    "units": "",
    "type": "list of strings"
  },
  "aooverlaps": {
    "desc": "atomic orbital overlap matrix",
    "units": "",
    "type": "array of rank 2"
  },
  "atombasis": {
    "desc": "indices of atomic orbitals on each atom",
    "units": "",
    "type": "list of lists"
  },
  "atomcharges": {
    "desc": "atomic partial charges",
    "units": "",
    "type": "dict of arrays of rank 1"
  },
  "atomcoords": {
    "desc": "atom coordinates",
    "units": "angstroms",
    "type": "array of rank 3"
  },
  "atommasses": {
    "desc": "atom masses",
    "units": "daltons",
    "type": "array of rank 1"
  },
  "atomnos": {
    "desc": "atomic numbers",
    "units": "",
    "type": "array of rank 1"
  },
  "atomspins": {
    "desc": "atomic spin densities",
    "units": "",
    "type": "dict of arrays of rank 1"
  },
  "ccenergies": {
    "desc": "molecular energies with Coupled-Cluster corrections",
    "units": "eV",
    "type": "array of rank 2"
  },
  "charge": {
    "desc": "net charge of the system",
    "units": "",
    "type": "integer"
  },
  "coreelectrons": {
    "desc": "number of core electrons in atom pseudopotentials",
    "units": "",
    "type": "array of rank 1"
  },
  "dispersionenergies": {
    "desc": "dispersion energy corrections",
    "units": "eV",
    "type": "array of rank 1"
  },
  "enthalpy": {
    "desc": "sum of electronic and thermal enthalpies",
    "units": "hartree/particle",
    "type": "float"
  },
  "entropy": {
    "desc": "entropy",
    "units": "hartree/particle*kelvin",
    "type": "float"
  },
  "etenergies": {
    "desc": "energies of electronic transitions",
    "units": "1/cm",
    "type": "array of rank 1"
  },
  "etoscs": {
    "desc": "oscillator strengths of electronic transitions",
    "units": "",
    "type": "array of rank 1"
  },
  "etdips": {
    "desc": "electric transition dipoles of electronic transitions",
    "units": "ebohr",
    "type": "array of rank 2"
  },
  "etveldips": {
    "desc": "velocity-gauge electric transition dipoles of electronic transitions",
    "units": "ebohr",
    "type": "array of rank 2"
  },
  "etmagdips": {
    "desc": "magnetic transition dipoles of electronic transitions",
    "units": "ebohr",
    "type": "array of rank 2"
  },
  "etrotats": {
    "desc": "rotatory strengths of electronic transitions",
    "units": "??",
    "type": "array of rank 1"
  },
  "etsecs": {
    "desc": "singly-excited configurations for electronic transitions",
    "units": "",
    "type": "list of lists"
  },
  "etsyms": {
    "desc": "symmetries of electronic transitions",
    "units": "",
    "type": "list of string"
  },
  "freeenergy": {
    "desc": "sum of electronic and thermal free energies",
    "units": "hartree/particle",
    "type": "float"
  },
  "fonames": {
    "desc": "fragment orbital names",
    "units": "",
    "type": "list of strings"
  },
  "fooverlaps": {
    "desc": "fragment orbital overlap matrix",
    "units": "",
    "type": "array of rank 2"
  },
  "fragnames": {
    "desc": "names of fragments",
    "units": "",
    "type": "list of strings"
  },
  "frags": {
    "desc": "indices of atoms in a fragment",
    "units": "",
    "type": "list of lists"
  },
  "gbasis": {
    "desc": "coefficients and exponents of Gaussian basis functions",
    "units": "",
    "type": "PyQuante format"
  },
  "geotargets": {
    "desc": "targets for convergence of geometry optimization",
    "units": "",
    "type": "array of rank 1"
  },
  "geovalues": {
    "desc": "current values for convergence of geometry optmization",
    "units": "",
    "type": "array of rank 1"
  },
  "grads": {
    "desc": "current values of forces (gradients) in geometry optimization",
    "units": "",
    "type": "array of rank 3"
  },
  "hessian": {
    "desc": "elements of the force constant matrix",
    "units": "",
    "type": "array of rank 1"
  },
  "homos": {
    "desc": "molecular orbital indices of HOMO(s)",
    "units": "",
    "type": "array of rank 1"
  },
  "metadata": {
    "desc": "various metadata about the package and computation",
    "units": "",
    "type": "dict"
  },
  "mocoeffs": {
    "desc": "molecular orbital coefficients",
    "units": "",
    "type": "list of arrays of rank 2"
  },
  "moenergies": {
    "desc": "molecular orbital energies",
    "units": "eV",
    "type": "list of arrays of rank 1"
  },
  "moments": {
    "desc": "molecular multipole moments",
    "units": "a.u.",
    "type": "list of arrays[]"
  },
  "mosyms": {
    "desc": "orbital symmetries",
    "units": "",
    "type": "list of lists"
  },
  "mpenergies": {
    "desc": "molecular electronic energies with M\u00f8ller-Plesset corrections",
    "units": "eV",
    "type": "array of rank 2"
  },
  "mult": {
    "desc": "multiplicity of the system",
    "units": "",
    "type": "integer"
  },
  "natom": {
    "desc": "number of atoms",
    "units": "",
    "type": "integer"
  },
  "nbasis": {
    "desc": "number of basis functions",
    "units": "",
    "type": "integer"
  },
  "nmo": {
    "desc": "number of molecular orbitals",
    "units": "",
    "type": "integer"
  },
  "nmrtensors": {
    "desc": "Nuclear magnetic resonance chemical shielding tensors",
    "units": "",
    "type": "dict of dicts of array of rank 2"
  },
  "nmrcouplingtensors": {
    "desc": "Nuclear magnetic resonance spin-spin coupling tensors",
    "units": "",
    "type": "dict of dicts of array of rank 2"
  },
  "nocoeffs": {
    "desc": "natural orbital coefficients",
    "units": "",
    "type": "array of rank 2"
  },
  "nooccnos": {
    "desc": "natural orbital occupation numbers",
    "units": "",
    "type": "array of rank 1"
  },
  "nsocoeffs": {
    "desc": "natural spin orbital coefficients",
    "units": "",
    "type": "list of array of rank 2"
  },
  "nsooccnos": {
    "desc": "natural spin orbital occupation numbers",
    "units": "",
    "type": "list of array of rank 1"
  },
  "optdone": {
    "desc": "flags whether an optimization has converged",
    "units": "",
    "type": "Boolean"
  },
  "optstatus": {
    "desc": "optimization status for each set of atomic coordinates",
    "units": "",
    "type": "array of rank 1"
  },
  "polarizabilities": {
    "desc": "(dipole) polarizabilities, static or dynamic",
    "units": "",
    "type": "list of arrays of rank 2"
  },
  "pressure": {
    "desc": "pressure used for Thermochemistry",
    "units": "atm",
    "type": "float"
  },
  "rotconsts": {
    "desc": "rotational constants",
    "units": "GHz",
    "type": "array of rank 2"
  },
  "scancoords": {
    "desc": "geometries of each scan step",
    "units": "angstroms",
    "type": "array of rank 3"
  },
  "scanenergies": {
    "desc": "energies of potential energy surface",
    "units": "",
    "type": "list"
  },
  "scannames": {
    "desc": "names of variables scanned",
    "units": "",
    "type": "list of strings"
  },
  "scanparm": {
    "desc": "values of parameters in potential energy surface",
    "units": "",
    "type": "list of tuples"
  },
  "scfenergies": {
    "desc": "molecular electronic energies after SCF (Hartree-Fock, DFT)",
    "units": "eV",
    "type": "array of rank 1"
  },
  "scftargets": {
    "desc": "targets for convergence of the SCF",
    "units": "",
    "type": "array of rank 2"
  },
  "scfvalues": {
    "desc": "current values for convergence of the SCF",
    "units": "",
    "type": "list of arrays of rank 2"
  },
  "temperature": {
    "desc": "temperature used for Thermochemistry",
    "units": "kelvin",
    "type": "float"
  },
  "time": {
    "desc": "time in molecular dynamics and other trajectories",
    "units": "fs",
    "type": "array of rank 1"
  },
  "transprop": {
    "desc": "all absorption and emission spectra (dictionary {name:",
    "units": "etoscs",
    "type": "etenergies"
  },
  "vibanharms": {
    "desc": "vibrational anharmonicity constants",
    "units": "1/cm",
    "type": "array of rank 2"
  },
  "vibdisps": {
    "desc": "cartesian displacement vectors",
    "units": "delta angstrom",
    "type": "array of rank 3"
  },
  "vibfreqs": {
    "desc": "vibrational frequencies",
    "units": "1/cm",
    "type": "array of rank 1"
  },
  "vibfconsts": {
    "desc": "force constants of vibrations",
    "units": "mDyne/angstrom",
    "type": "array of rank 1"
  },
  "vibirs": {
    "desc": "IR intensities",
    "units": "km/mol",
    "type": "array of rank 1"
  },
  "vibramans": {
    "desc": "Raman activities",
    "units": "A^4/Da",
    "type": "array of rank 1"
  },
  "vibrmasses": {
    "desc": "reduced masses of vibrations",
    "units": "daltons",
    "type": "array of rank 1"
  },
  "vibsyms": {
    "desc": "symmetries of vibrations",
    "units": "",
    "type": "list of strings"
  },
  "zpve": {
    "desc": "zero-point vibrational energy correction",
    "units": "hartree/particle",
    "type": "float"
  }
}
    parsed_str=[f"{x:20s} {pdata[x]['type']:40s} {pdata[x]['units']:20s} {pdata[x]['desc']}" for x in pdata]
    technical="{:20s} {:40s} {:20s} {}\n".format("Technical","None","None","Technical parameters of the calculation")
    help_str="{:20s} {:40s} {:20s} {}\n".format("Kewyword","Type","Units","Desciption")+technical+"\n".join(parsed_str)
    return(pdata,help_str)


        
        
    

def main():
    pdata,hlp=parsed_data()
    sw_list=["ADF","DALTON","FChk","GAMESS","GAMESSDAT","GAMESSUK","Gaussian","Jaguar","Molpro","Molcas","MOPAC","NWChem","ORCA","Psi4","QChem","Turbomole"]
    
    parser = argparse.ArgumentParser(description="Parse File and give Properties")
    
    prop = parser.add_argument_group('Properties')
    prop.add_argument("--prop", 
                      nargs='+',
                      default="scfenergies",
                      help='Select relevant property. Selection help print available properties')
    
    sw = parser.add_argument_group('Supported Softwares')
    sw.add_argument('--sw',help='Program(s) output to extract data from.\n \
                    Currently supported: '+"\n".join(sw_list))
    
    output_pars=parser.add_argument_group('Output Parameters')
    
    
    
    output_pars.add_argument('--outputtype',
                        choices=('json', 'cjson', 'cml', 'xyz', 'molden', 'wfx'),
                        help='the output format to write (json/cjson are identical)')
    
    output_pars.add_argument('--input',
                            nargs='+',
                            help='Input files')

    output_pars.add_argument('--full',
                            action='store_true',
                            help='print full parsed data in single json files')
    output_pars.add_argument('--recursive',
                            action='store_true',
                            help='Look recursively for input files')
    
    output_pars.add_argument('-v', '--verbose',
                        action='store_true',
                        help='more verbose parsing output (only errors by default)')

    output_pars.add_argument('-g', '--ghost',
                        type=str,
                        default=None,
                        help='Symbol to use for ghost atoms')

    output_pars.add_argument('-t', '--terse',
                        action='store_true',
                        help='CJSON by default is not indented for readability, saves space (indented for readability\'s sake)')

    output_pars.add_argument('-u', '--future',
                        action='store_true',
                        help='use experimental features (currently optdone_as_list)')

    output_pars.add_argument('-i', '--index',
                        type=int,
                        default=None,
                        nargs='+',
                        help='optional zero-based index to extract a specific element of the property array')
    output_pars.add_argument('-c', '--cols',
                        default=False,
                        action='store_true',
                        help='columns layout')
    output_pars.add_argument('-m', '--meta',
                        default=["concise"],
                        type=str,
                        nargs='+',
                        help='Metadata to print. Accepts one of the values in the dict. Full prints all metadata\
                        and concise a selection of the most common parameters')
    
        
    args = parser.parse_args()
    
    if  'help' in args.prop:
        print(hlp)
        sys.exit()
    
    #filenames="LlOQ_tolyl_ACN_S0opt_exc_TD_PBE1PBE.log"
    filenames=[]
    for file in args.input:
        if args.recursive:
            filenames+=glob.glob(file,recursive=True)
        else:
            filenames+=glob.glob(file)
    prp_print=[]
    for i,prp in enumerate(args.prop):
      try:
        prp_print.append((prp,args.index[i]))
      except:
        prp_print.append((prp,"all"))
    print(f"Properties (requested, index): {prp_print} ")
    print(f"Files Found: {len(filenames)}")
    files=dict()
    for file in filenames:
        try:
            data=ccopen(file)
            files[file]=f"{data}".split()[0]
        except:
            print(f"Cannot read file {file}")
        #for prop in args.
    #print(f"{files}".split()[0])
    #print(files)
    #types=guess_filetype(files)
    #print(types)

    outputtype = args.outputtype
    verbose = args.verbose
    terse = args.terse
    future = args.future
    index = args.index
    ghost = args.ghost
    props_all=[]
    concise_metadata=["methods","functional","basis_set","success"]

    for filename in files:
        technical=parse_tech(filename,files[filename])
        # We might want to use this option in the near future.
        ccopen_kwargs = dict()
        if future:
            ccopen_kwargs['future'] = True

        print(f"Attempting to parse {filename}")
        log = ccopen(filename, **ccopen_kwargs)

        if not log:
            print(
                f"Cannot figure out what type of computational chemistry output file '{filename}' is."
            )
            print("Report this to the cclib development team if you think this is an error.")
            sys.exit()

        if verbose:
            log.logger.setLevel(logging.INFO)
        else:
            log.logger.setLevel(logging.ERROR)
        try:
            data = log.parse()
            parsed=True
        except:
            print(f"Error parsing logfile {filename}")
            parsed=False
        #print(f"cclib can parse the following attributes from {filename}:")
        #hasattrs = [f"  {attr}" for attr in ccData._attrlist if hasattr(data, attr)]
        #print("\n".join(hasattrs))
        
        
        if parsed:
            if args.full:
                # Write out to disk.
                outputdest = '.'.join([os.path.splitext(os.path.basename(filename))[0], outputtype])
                ccwrite_kwargs = dict()
                if future:
                    ccwrite_kwargs['future'] = True
                if ghost:
                    ccwrite_kwargs['ghost'] = ghost
                # For XYZ files, write the last geometry unless otherwise
                # specified.
                if not index:
                    index = -1
                ccwrite_kwargs['jobfilename'] = filename

                # The argument terse presently is only applicable to
                # CJSON/JSON formats
                ccwrite(data, outputtype, outputdest, indices=index, terse=terse,
                        **ccwrite_kwargs)
                
                technical = {'pars':technical[0],'num_calc':technical[1]}
                with open(outputdest+"_technical",'w') as jsofile:
                    json.dump(technical,jsofile,indent=4)
                
            
            else:
                props=[filename]
                #print(args.meta)
                if args.meta[0] == 'full':
                    props.append(getattr(data,'metadata'))
                elif args.meta[0] == 'concise':
                  args.meta=concise_metadata
                  for tech in concise_metadata:
                    try:
                      props.append(getattr(data,'metadata')[tech])
                    except:
                      props.append("not found")
                else:
                  for tech in args.meta:
                    try:
                      props.append(getattr(data,'metadata')[tech])
                    except:
                      props.append("not found")
                for i,prop in enumerate(args.prop):
                    try:
                        props.append(getattr(data,prop)[args.index[i]])
                    except:
                        try:
                            props.append(getattr(data,prop))
                        except:
                            print(f"{prop} not found for {filename}")
                                
                props_all.append(props)
    #print(props_all)
    prop_string=["Name"]+ [x for x in args.meta] + [f"{x}" for x in args.prop]
    prop_string_unit=[[" "] + [x for x in args.meta] + [pdata[x]['units'] for x in args.prop]]
    #print(prop_string)
    #print(prop_string_unit)
    datf_units=pd.DataFrame(prop_string_unit,columns=prop_string)
    datf_props=pd.DataFrame(props_all,columns=prop_string)
    datf=pd.concat([datf_units,datf_props])
    #print(datf_props)
    #print(datf[datf['Name']=="2rings_DHI.log"]['etenergies'])
    pd.DataFrame.to_excel(datf,"Properties.xlsx")
    with open("Properties.csv",'w') as ofile:
      if args.cols:
        name="Name"
        ofile.write(f"{name:30s}")
        for prop in range(len(props_all[0])):
          if prop_string[prop] in list(pdata.keys()):
                    string=f"{prop_string[prop]} ({pdata[prop_string[prop]]['units']})"
                    ofile.write(f"{string:30s}")
          for comp in props_all:
            if isinstance(comp[prop],list):
                    for subprop in comp[prop]:
                        ofile.write(f"{subprop};")
            else:
                    ofile.write(f"@{comp[prop]}")
          ofile.write("\n")
        
      else:
        for line in props_all:
            for i,prop in enumerate(line):
                if prop_string[i] in list(pdata.keys()):
                    string=f"{prop_string[i]} ({pdata[prop_string[i]]['units']})"
                    ofile.write(f"{string:30s}@")
                else:
                    ofile.write(f"{prop_string[i]:30s}@")
                if isinstance(prop,list):
                    for subprop in prop:
                        ofile.write(f"{subprop};")
                else:
                    ofile.write(f"{prop}")
                ofile.write("\n")
            ofile.write("\n")
    
    #print(f"Files found: {len(files)}")        
    print(f"Files correctly parsed: {len(props_all)}")
                        
        

        
        


    
if __name__=="__main__":
    main()
    



