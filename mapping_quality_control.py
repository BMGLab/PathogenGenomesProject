import pysam
import os
import sys



def get_chromososmes_count(samfile):
    counter = 1
    f = open(samfile , "r")
    f.readline()
    line = f.readline().split()
    example = line[1][:5]
    line = f.readline()
    while(example in line):
        counter += 1
        line = f.readline()
    f.close()
    return counter




def get_mapped_nucleotides(reference_sequences , bamfile , result_dictionary):
    for ref in reference_sequences:
        counter = 0
        for i in bamfile.pileup(ref):
            counter+=1
        result_dictionary[ref] = counter






def get_percentages_of_mapped_nucleotides(samfile , chromosome_count , result_dictionary):
    f = open(samfile , "r")
    f.readline()
    for i in range(chromosome_count):
        line = f.readline().split()
        val = int(line[2][3:])
        chro = line[1][3:]
        result_dictionary[chro] = (result_dictionary[chro] / val)*100

def print_results(result_dictionary):
    for name , percentage in result_dictionary.items():
        print(f"{name} --> {percentage:.5f}%")
    



def main():

    
    file = sys.argv[1]
    samfile = "sorted_sam_file.sam"

    bamfile = pysam.AlignmentFile(file , "rb")
    os.system(f"samtools view -h {file} > {samfile}")

    chromosome_count = get_chromososmes_count(samfile)

    reference_sequences = bamfile.references[:chromosome_count]

    result_dictionary = {}  

    get_mapped_nucleotides(reference_sequences , bamfile , result_dictionary)

    get_percentages_of_mapped_nucleotides(samfile , chromosome_count , result_dictionary)

    print_results(result_dictionary)


main()
