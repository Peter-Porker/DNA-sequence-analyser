#This program allows you to input a DNA sequence in a fasta file format and analyze the length of the sequence, nucleotide ratio, and the GC content
#It also returns the complement sequence and the orf sites present across the 3 reading frames.

def read_fasta(filepath):           #Function reads the filepath entered and retrieves the sequence 
    sequence = []                   #This list will store the DNA sequence locally in the function
    reading = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            
            if not line: continue       #Omits blank lines
            if line.startswith(">"):    #Checks if line starts with '>'
                reading = True          #and starts reading the DNA sequence
                continue
            if reading:
                sequence.append(line)   #stores the sequence into the list
            
    return "".join(sequence).upper()    #returns the sequence as a string when the function is called

def validate(seq):                              #Validates if the sequence contains only the nucleotides "ATGC"
    if not seq.isalpha():
        return False
    for x in seq:
        if x not in {"A", "T", "G", "C"}:
            return False                        #Returns False if sequence contains symbols, numbers, or does not belong to "ATGC"
    return True

def analysis(seq):
    print("Length: ",len(seq))                                      #Prints the length of the DNA sequence
    num_A = seq.count("A")
    num_T = seq.count("T")
    num_G = seq.count("G")
    num_C = seq.count("C")
    GC = (num_G+num_C)/len(seq)*100                                 #Calculates the GC content using formula: [(Number of G)+(Number of C)]/length of sequence
    print('A:',num_A,' | T:',num_T,' | G:',num_G,' | C:',num_C)
    print('GC content: ',GC,'%')                                    #Prints the results from analysis

def complement(seq):                                        #Function maps the nucleotide to their complement and returns the complementary sequence
    mapping = {"A":"T","T":"A","G":"C","C":"G"}             #Dictionary stores the mapping
    return "".join(mapping[x] for x in seq)

def orf(seq,frame):                                     #Function identifies the Open Reading Frames in the sequence
    orfs =[]
    stops = ["TAG","TGA","TAA"]                         # 3 stop codons
    for i in range(frame,len(seq)-2,3):                 #For loop goes through the entire sequence
        if seq[i:i+3] == "ATG":                         #Conditional checks for the start codon 'ATG'
            for j in range(i,len(seq)-2,3):             #For loop searches from the start codon to the end of the sequence    
                if seq[j:j+3] in stops:                 #Conditional checks if codon belongs to stop codon list
                    orfs.append({                       #Entering details of the ORF into 'orfs' list
                        "start":i,
                        "end":j+3,
                        "length":j+3 - i,
                        "sequence":seq[i:j+3]
                    })
                    break
    return orfs                                         #Returns the ORFs found as a list when function is called

def main():
    filepath = input("Enter the file location: ")       #Taking the file location as input
    sequence = read_fasta(filepath)

    if not validate(sequence):                          #Displays error if validate() returns False
        print("Error in sequence provided.")
        return
    
    analysis(sequence)                                  #Function prints the results from analysis

    inp= input("Print complement sequence? [y/n] ")     #Checks if user needs complement sequence printed
    if inp == "y":
        print(complement(sequence))                     #Prints the complement sequence
    elif inp == "n":
        pass
    else: 
        print("Invalid input. Skipping...")             #Skips function if rejected or invalid input
        pass
    
    for frame in range(3):                              #Searches for ORFs in sequence in 3 frames and 
        orfs = orf(sequence,frame)                      #stores the results in 'orfs'
        if not orfs:
            print("No ORFs found")
            return
        else: print(f"\nFrame {frame+1}")

        #Printing the Start and End positions, the Length and ORF sequences found in each frame
        for x in orfs:
            print(f"ORF: Start={x['start']}, End={x['end']}, Length={x['length']}, \nSequence={x['sequence']}\n")              


main()




