#!/usr/bin/env python
#functions made throughout classes
#Python version 3.9.5

DNAbases = set('ATGCatcg')
RNAbases = set('AUGCaucg')

def validate_base_seq(seq: str,RNAflag: bool=False) -> bool:
    '''This function takes a string. Returns True if string is composed
    of only As, Ts (or Us if RNAflag), Gs, Cs. False otherwise. Case insensitive.'''
    return set(seq)<=(RNAbases if RNAflag else DNAbases)

if __name__ == "__main__":
    assert validate_base_seq("AATAGAT") == True, "Validate base seq does not work on DNA"
    assert validate_base_seq("AAUAGAU", True) == True, "Validate base seq does not work on RNA"
    assert validate_base_seq("Hi there!") == False, "Validate base seq fails to recognize nonDNA"
    assert validate_base_seq("Hi there!", True) == False, "Validate base seq fails to recognize nonDNA"
    print("Passed DNA and RNA tests")

def gc_content(DNA: str):
    '''Returns GC content of a DNA sequence as a decimal between 0 and 1.'''
    assert validate_base_seq(DNA) == True, "String contains invalid characters"
    
    DNA = DNA.upper()
    Gs = DNA.count("G")
    Cs = DNA.count("C")
    return (Gs+Cs)/len(DNA)

if __name__ == "__main__":
    assert gc_content("AAAAAATTTTAT") == 0, "gc_content does not work with no G's or C's"
    assert gc_content("GGGGGGGCCCCCGCGC") == 1, "gc_content does not work with only G's or C's in DNA"
    assert gc_content("ATATATGCGCGC") == 0.5, "gc_content does not work with mix of all nucleotides"
    print("Passed g_c_content tests!")

#convert phred score into numerical score from 0 - 42
def convert_phred(letter):
    """Converts a single character into a phred score"""
    return ord(letter) - 33

#returns average quality score
def ave_qual_score(phred_score):
    qual_sum = 0
    for letter in phred_score:
        qual_sum += convert_phred(letter)
    return qual_sum/len(phred_score)

#combines wrapped fasta lines
def fasta_line_combiner(filein, fileout):
    """
    takes in fasta file, loops through lines, and combines sequence lines to make them a single line. 
    Output is fasta file with combined sequence segments.
    """
    with open(fileout, "w") as out:
        with open(filein, "r") as fl:
            new_pep = ""
            for line in fl:
                line = line.strip()
                if line.startswith(">"):
                    if new_pep != "":
                        out.writelines(new_pep + "\n")
                        new_pep = ""
                    out.writelines(line + "\n")

                else:
                    new_pep = new_pep + line.strip()
                
            out.writelines(new_pep)

