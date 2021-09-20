#
# Copyright 2020- IBM Inc. All rights reserved
# SPDX-License-Identifier: Apache2.0
#

import gzip
from pathlib import Path
from Bio import SeqIO
import xml.etree.ElementTree as ET
import regex as re
import argparse
import json

class UniProtParser():
    
    def __init__(self):
        pass

    def parse(self, filename):
        fileName = Path(filename)

        open_ = dict(
            gzip=gzip.open,
            gz=gzip.open,
        )

        with open_.get(fileName.suffix[1:], open)(str(fileName), 'r') as handle:
            #SeqIO.UniprotIO.UniprotIterator(handle, alphabet=ProteinAlphabet(), return_raw_comments=False)
            for record in SeqIO.UniprotIO.UniprotIterator(handle,return_raw_comments=True):
                #for record in SeqIO.parse(handle, "uniprot-xml",return_raw_comments=True):
                doc = {}
                d = record.annotations

                #doc['id'] = record.id

                doc['accessions'] = d['accessions']
                if 'recommendedName_fullName' in d:
                    doc['recommended_name'] = d['recommendedName_fullName']
                if 'recommendedName_ecNumber'  in d: 
                    doc['ec_number'] = d['recommendedName_ecNumber']
                
                if 'gene_name_ORF'             in d: 
                    doc['gene_name_ORF'] = d['gene_name_ORF']
                
                if 'gene_name_synonym'         in d: 
                    doc['gene_name_synonym'] = d['gene_name_synonym']
                
                if 'gene_name_ordered locus'   in d: 
                    doc['gene_name_ordered locus'] = d['gene_name_ordered locus']
                
                if 'organism_name'             in d: 
                    doc['organism_name'] = d['organism_name']

                doc['taxonomy'] = d['taxonomy']
                #lineage     # maybe in record.features?

                refs = []
                for ref_ in d['references']:
                    #print(ref_)

                    tmp = {
                        "title": ref_.title,
                        "authors": ref_.authors,
                        "journal": ref_.journal,
                        "pubmed_id": ref_.pubmed_id
                        }
                    refs.append(tmp)

                doc['references'] = refs

                if 'comment_function' in d: 
                    doc['comment_function'] = d['comment_function']

                if 'comment_catalyticactivity' in d:
                    doc['comment_catalyticactivity'] = d['comment_catalyticactivity']
                #

                #extract kinetics
                if 'comment_biophysicochemicalproperties_xml' in d:
                    kms = []
                    for xmlstring in d['comment_biophysicochemicalproperties_xml']:
                        root = ET.fromstring(xmlstring)
                        for child in root[0]:
                            m = re.search('\}KM$',child.tag)
                            if m:
                                #print('tag:',child.tag, ' text:', child.text)
                                km = {}
                                #4.1 uM for AMP
                                stext = child.text.split(' for ')
                                st0 = stext[0].split()

                                #> convert to mM
                                conversion = 1.0
                                if st0[1] == 'uM': 
                                    conversion = 1.0e-3
                                    st0[1] = 'mM'
                                #<
                                
                                try:
                                    km['value'] = float( st0[0] ) * conversion
                                except:
                                    continue
                                #
                                km['units'] = st0[1]
                                st1 = stext[1].split('(')
                                km['substance'] = st1[0].strip().lower()
                                if len(st1)>1:km['conditions'] = st1[1][0:-1]
                                km['text'] = child.text
                                kms.append(km)
                            #
                        #
                    #
                    if kms != []:
                        doc['km_value'] = kms
                #


                #catalytic activity references
                dbrefs = {}
                for item in record.dbxrefs:
                    kv = item.split(':', 1)
                    key = kv[0]
                    value = kv[1]

                    dbrefs.setdefault(key, []).append(value)

                doc['dbxrefs'] = dbrefs
                doc['sequence']     = str(record.seq)
                doc['url'] = f'https://www.uniprot.org/uniprot/{record.id}.txt'
                #print(doc)
                yield doc
            #
        #
    #

#

def parse_arguments():
    argparser = argparse.ArgumentParser(description="parser package")

    argparser.add_argument(
        '--input', action='store',
        default=None,
        dest='input_file',
        help='Uniprot data file (uniprot_sprot.xml.gz)')

    argparser.add_argument(
        '--output', action='store',
        default=None,
        dest='output_file',
        help='output collection file (jsonl)')

    pargs = argparser.parse_args()

    return pargs.input_file, pargs.output_file
#


def main():

    input_file, output_file = parse_arguments()

    up = UniProtParser()

    with open(output_file, 'w') as fod:
        for doc in up.parse(input_file):
            fod.write(json.dumps(doc)+'\n')
        #
    #
#


if __name__ == "__main__":
    main()
#
