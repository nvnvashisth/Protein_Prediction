# -*- coding: utf-8 -*-
"""
IMPORTANT!:
Before writing an email asking questions such as
'What does this input has to be like?' or 
'What return value do you expect?' PLEASE read our
exercise sheet and the information in this template
carefully.
If something is still unclear, PLEASE talk to your
colleagues before writing an email!

If you experience technical issues or if you find a
bug we are happy to answer your questions. However,
in order to provide quick help in such cases we need 
to avoid unnecessary emails such as the examples
shown above.
"""

from Bio import SeqIO # Tip: This module might be useful for parsing... 
import re
############ Exercise 3: SwissProt ##########
class SwissProt_Parser:
    
    #path = "P09616.xml"
    PARSER = SeqIO
    
    def __init__( self, path, frmt='uniprot-xml' ):
        '''
            Initialize every SwissProt_Parser with a path to a XML-formatted UniProt file.
            An example file is included in the repository (P09616.xml).
            Tip: Store the parsed XML entry in an object variable instead of parsing it
            again & again ...
        '''
        self.sp_anno = list(self.PARSER.parse(path, frmt))
        #self.sp_anno = None # Parse the XML file once and re-use it in the functions below
        
    # 3.2 SwissProt Identifiers
    def get_sp_identifier( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Unique SwissProt identifier for the given xml file
        '''
        for i in self.sp_anno:
            identifier = i.id
        #print(identifier)
        return identifier
        
    # 3.3 SwissProt Sequence length
    def get_sp_sequence_length( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return sequence length of the UniProt entry as an integer.
        '''
        
        for i in self.sp_anno:
            seq_len = i.annotations['sequence_length']
                # print(i)
        #print(seq_len)
        return seq_len
    
    # 3.4 Organism 
    def get_organism( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return the name of the organsim as stated in the corresponding field
                of the XML data. Return value has to be a string.
        '''
        for i in self.sp_anno:
            organism = i.annotations['organism']
        #print(organism)
        return organism
    
    # 3.5 Localizations
    def get_localization( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Return the name of the subcellular localization as stated in the 
                corresponding field.
                Return value has to be a list of strings.
        '''
        localization = []
        for i in self.sp_anno:
            localization = i.annotations['comment_subcellularlocation_location']
        #localization = ['localization_1', 'localization_2']
        #print(type(localization))
        return localization
    
    # 3.6 Cross-references to PDB
    def get_pdb_support( self ):
        '''
            Input: 
                self: Use XML entry which has been parsed & saved during object initialization 
            Return:
                Returns a list of all PDB IDs which support the annotation of the
                given SwissProt XML file. Return the PDB IDs as list.
        '''
        pdb_ids = list()
        for records in self.sp_anno:
            dref = records.dbxrefs
        #print(dref)
        for item in dref:
            if re.search(r"PDB:", item):
                pdb_ids.append(item[4:])
        return pdb_ids


#__init__("P09616.xml")

def main():
    print('SwissProt XML Parser class')
    #__init__("P09616.xml")
    return None
    
if __name__ == '__main__':
    main()
    p1 = SwissProt_Parser("P09616.xml")
    #p1.get_sp_identifier()
    #p1.get_sp_sequence_length()
    #p1.get_organism()
    #p1.get_localization()
    #p1.get_pdb_support()
    #__init__()
