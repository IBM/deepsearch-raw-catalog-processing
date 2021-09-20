#
# Copyright 2020- IBM Inc. All rights reserved
# SPDX-License-Identifier: Apache2.0
#

import sys
import os

BRENDA_PARSER_DIR = os.path.abspath(
    os.path.join(
        os.path.dirname(__file__),
        'BRENDA-Parser'
    )
)
sys.path.append(BRENDA_PARSER_DIR)
from brenda.parser import BRENDAParser

import regex as re
import argparse
import json


class BrendaTXTParser():

    def __init__(self):

        self._section_pairs = {
            "AC":   "activating_compound",
            "AP":   "application",
            "CF":   "cofactor",
            "CL":   "cloned",
            "CR":   "crystallization",
            "CR":   "cas_registry_number",
            "EN":   "engineering",
            "GS":   "general_stability",
            "IC50": "ic50_value",
            "IN":   "inhibitors",
            "KI":   "ki_value",
            "KM":   "km_value",
            "LO":   "localization",
            "ME":   "metals_ions",
            "MW":   "molecular_weight",
            "NSP":  "natural_substrate_product",
            "OS":   "oxidation_stability",
            "OSS":  "organic_solvent_stability",
            "PHO":  "ph_optimum",
            "PHR":  "ph_range",
            "PHS":  "ph_stability",
            "PI":   "pi_value",
            "PM":   "posttranslational_modification",
            "PR":   "protein",
            "PU":   "purification",
            "RE":   "reaction",
            "RF":   "reference",
            "RN":   "renatured",
            "RN":   "recommended_name",
            "RT":   "reaction_type",
            "SA":   "specific_activity",
            "SN":   "systematic_name",
            "SP":   "substrate_product",
            "SS":   "storage_stability",
            "ST":   "source_tissue",
            "SU":   "subunits",
            "SY":   "synonyms",
            "TN":   "turnover_number",
            "TO":   "temperature_optimum",
            "TR":   "temperature_range",
            "TS":   "temperature_stability"
        }

    def get_mapping_name(self):
        return "cps_extdbs/resources/mapping_Brenda_txt.json"

    def _append_reactants_products(self, doc):

        if 'substrate_product' in doc:
            reactant_tmp = set()
            product_tmp = set()
            for reaction in doc['substrate_product']:
                sides = reaction.split('=')
                if len(sides) != 2:
                    continue
                else:
                    for g in sides[0].split(' + '):
                        g = g.strip().lower()
                        reactant_tmp.add(g)
                    #
                    for g in sides[1].split(' + '):
                        g = g.strip().lower()
                        product_tmp.add(g)
                    #
                #
            #
            doc['products'] = list(product_tmp)
            doc['reactants'] = list(reactant_tmp)
        #
    #

    def _append_proteins(self, doc, proteins):

        tmp = []
        identifiers = []
        for prot_id, protein in proteins.items():
            #print(prot_id) 
            #print('\tInformation:', protein.information) 
            #print('\tComment:', protein.comment) 
            #print('\tOrganism:', protein.organism) 
            #print('\tIdentifiers:', protein.identifiers) 
            #print('\tReferences:', protein.references)
            tmp.append(protein.organism)
            identifiers += protein.identifiers
        #
        if tmp != []: doc['organism'] = tmp
        if identifiers != []: doc['identifiers'] = identifiers
    #

    
    def _append_entries(self, doc, name, entries, proteins):

        tables = {
            'KM_VALUE'            : {'key_name':'substrate'},
            'TURNOVER_NUMBER'     : {'key_name':'substrate'},
            'KCAT_KM_VALUE'       : {'key_name':'substrate'},
            'PH_RANGE'            : {'key_name':None,      },
            'TEMPERATURE_OPTIMUM' : {'key_name':None,      }
        }

        tmp = []

        for i_entry,entry in enumerate(entries):
            #vwe there is a bug in entries for eg 3.2.1.22. the last entry in duplicated for kcat_km_value (6555).
            #vwe we bypass the issue by avoinding repetition in the table.

            if name in tables:

                key_name = tables[name]['key_name']

                if entry.msg == '-999': continue
                count = entry.msg.count('-')
                if   count == 0:
                    try:
                        float( entry.msg )
                    except:
                        continue
                    #
                    value = {'value':{ 'min':float( entry.msg ), 'max':float( entry.msg ) }}
                elif count == 1:
                    v_split = entry.msg.split('-')
                    vmin = v_split[0]
                    vmax = v_split[1]
                    try:
                        float( vmin )
                        float( vmax )
                    except:
                        continue
                    #
                    value = {'value':{ 'min':float( vmin ), 'max':float( vmax ) }}
                else:
                    continue
                #

                # entry.comment -> complex object
                if entry.proteins: 
                    for prot_id in entry.proteins:
                        if prot_id not in proteins: continue
                        comments = []
                        if entry.comment:
                            for com in entry.comment.msg.split(';'):
                                m=re.search(r'(?:#\d+#\s*)(.*)(?:\s*<\d+>)',com)
                                if m: comments.append(m.group(1).strip())
                            #
                        #
                        if key_name:
                            data = {key_name:str(entry.information).lower().strip(),
                                    'organism':proteins[prot_id].organism,
                                    'uniprot_id':proteins[prot_id].identifiers,
                                    'comment':comments}
                        else:
                            data = {'organism':proteins[prot_id].organism,
                                    'uniprot_id':proteins[prot_id].identifiers,
                                    'comment':comments}
                        #
                        data.update(value)
                        if data not in tmp: tmp.append(data)
                    #
                else:
                    if key_name:
                        data = {key_name:str(entry.information).lower().strip(),
                                'organism':'',
                                'uniprot_id':'',
                                'comment':comments}
                    else:
                        data = {'organism':'',
                                'uniprot_id':'',
                                'comment':comments}
                    #
                    data.update(value)
                    if data not in tmp: tmp.append(data)
                #
            else:
                if entry.msg not in tmp: tmp.append(entry.msg)
            #
        #
        #if doc['ec_number'] == '3.2.1.22' and name.lower() == 'kcat_km_value': print('<HERE2>',tmp)
        
        doc[name.lower()] = tmp
    #
    

    def parse(self, filename):

        bp = BRENDAParser
        bp._sections['KCAT_KM_VALUE']= 'KKM'
        
        with bp(filename) as parser:
            brenda = parser.parse()
        #

        for major in ['1','2','3','4','5','6','7']:

            if major not in brenda: continue

            for enzyme in brenda[major]:
                #print('e',enzyme.ec_number,enzyme)

                doc = {}
                doc['ec_number'] = enzyme.ec_number
                #
                self._append_proteins(doc, enzyme.proteins)
                #
                for name,entries in enzyme.entries.items():
                    #print('>>>>',name)
                    #print('name',name,'entries:',entries)
                    #if enzyme.ec_number == '3.2.1.22': print('<HERE0>',name,entries)
                    self._append_entries(doc, name, entries, enzyme.proteins)
                #
                self._append_reactants_products(doc)

                #print(sorted(list(doc.keys())))
                #if 'ph_range' in doc: print('>>>',doc['ph_range'])
                #if enzyme.ec_number == '3.2.1.22': print('<HERE>',doc['kcat_km_value'])
                yield doc
            #
        #
    #
#


def parse_arguments():
    argparser = argparse.ArgumentParser(description="parser package")

    argparser.add_argument(
        '--input', action='store',
        default='brenda_download.txt',
        dest='input_file',
        help='Brenda data file (brenda_download.txt)')

    argparser.add_argument(
        '--output', action='store',
        default='brenda_download.jsonl',
        dest='output_file',
        help='output collection file (jsonl)')

    pargs = argparser.parse_args()

    return pargs.input_file, pargs.output_file
#


def main():

    input_file, output_file = parse_arguments()

    bp = BrendaTXTParser()

    with open(output_file, 'w') as fod:
        for doc in bp.parse(input_file):
            fod.write(json.dumps(doc)+'\n')
        #
    #
#


if __name__ == "__main__":
    main()
#
