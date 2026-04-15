#This program allows you to input a DNA sequence in a fasta file format and analyze the length of the sequence, nucleotide ratio, and the GC content
#It also returns the complement sequence and the orf sites present across the 3 reading frames.

def read_fasta(filepath):
    sequence = []
    reading = False

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            
            if not line: continue
            if line.startswith(">"):
                reading = True
                continue
            if reading:
                sequence.append(line)
            
    return "".join(sequence).upper()

def validate(seq):
    if not seq.isalpha():
        return False
    for x in seq:
        if x not in {"A", "T", "G", "C"}:
            return False
    return True

def analysis(seq):
    print("Length: ",len(seq))
    num_A = seq.count("A")
    num_T = seq.count("T")
    num_G = seq.count("G")
    num_C = seq.count("C")
    GC = (num_G+num_C)/len(seq)*100
    print('A:',num_A,' | T:',num_T,' | G:',num_G,' | C:',num_C)
    print('GC content: ',GC,'%')

def complement(seq):
    mapping = {"A":"T","T":"A","G":"C","C":"G"}
    return "".join(mapping[x] for x in seq)

def orf(seq,frame):
    orfs =[]
    stops = ["TAG","TGA","TAA"]
    for i in range(frame,len(seq)-2,3):
        if seq[i:i+3] == "ATG":
            for j in range(i,len(seq)-2,3):
                if seq[j:j+3] in stops:
                    orfs.append({
                        "start":i,
                        "end":j+3,
                        "length":j+3 - i,
                        "sequence":seq[i:j+3]
                    })
                    break
    return orfs

def main():
    filepath = input("Enter the file location: ")
    sequence = read_fasta(filepath)

    if not validate(sequence):
        print("Error in sequence provided.")
        return
    
    analysis(sequence)

    inp= input("Print complement sequence? [y/n] ")
    if inp == "y":
        print(complement(sequence))
    elif inp == "n":
        pass
    else: 
        print("Invalid input. Skipping...")
        pass
    
    for frame in range(3):
        orfs = orf(sequence,frame)
        if not orfs:
            print("No ORFs found")
            return
        else: print(f"\nFrame {frame}")

        for x in orfs:
            print(f"ORF: Start={x['start']}, End={x['end']}, Length={x['length']}, Sequence={x['sequence']}")


main()




