"""
@ Andrew Li and Arul Prasad, PRIMES 2014
@ Script purpose: Given several objects, lengths, the location of drug resistant SNPs within the proteins, and arrays of ranges of disorder within the protein, the program will output two textfiles of distances from the SNPs to every other amino acid in each protein.
This was used to show the correlation between SNPs and disordered regions in cancer related proteins.
"""

from pymol import cmd
import re  # for splitting


def residueDist():
    names = []
    for obj in names:
        cmd.load(obj + ".pdb")
        fin_d = open("/home/andrew/pymol/" + obj + "_disorder.txt", "r")
        fin_s = open("/home/andrew/pymol/" + obj + "_snps.txt", "r")
        # split the Fasta-like input file and the SNPs file
        d_input = fin_d.read().split()
        s_input = fin_s.read()
        snps = re.split(r"[\n,\s]", s_input)
        for e in snps:
            if e == "":
                snps.remove(e)
        d_list = []
        for i in range(1, len(d_input) - 2):
            # gets the two numbers that denote the range of the disordered region
            d_range = re.split(r"[\n,\s]", d_input[i])
            for e in d_range:
                if e == "":
                    d_range.remove(e)
            d_temp = str(d_range[0])
            d_temp_array = re.findall(r"\d+", d_temp)
            d_list.append(d_temp_array)

        # length of protein and all the disordered regions
        length = int(d_input[len(d_input) - 2])
        start = int(d_input[len(d_input) - 1])
        fout = open("/home/andrew/pymol/" + obj + "_n_results.txt", "w")  # normal results
        fout_d = open("/home/andrew/pymol/" + obj + "_d_results.txt", "w")
        for i in range(len(snps)):
            for j in range(start, length + start - 1):
                if snps[i] != j:
                    temp_dist = calcDist(obj, i, j)
                    # if temp_dist is -1, it means that the residue could not be found
                    if temp_dist == -1:
                        pass
                    elif (inDisorder(j, d_list)):
                        fout_d.write(str(int(temp_dist)) + "\n")
                    else:
                        fout.write(str(int(temp_dist)) + "\n")
                else:
                    pass
        fin_d.close()
        fin_s.close()
        fout.close()
        fout_d.close()
        print obj + " has been analyzed."
        cmd.delete(obj)


def calcDist(obj, l1, l2):
    # obj is the protein in question, l1, l2 are the positions of residues.
    # delete previous selections and displayed distances
    print ""
    cmd.delete("dist")
    cmd.delete("sele")
    cmd.delete("Residue_1")
    cmd.delete("Residue_2")

    # select two residues based on input, at positions l1 and l2 in the sequence
    cmd.select("Residue_1", "resi {}".format(l1))
    cmd.select("Residue_2", "resi {}".format(l2))

    #the underscores ARE IMPORTANT
    distance = cmd.distance("dist", "Residue_1", "Residue_2")
    return distance


def inDisorder(i, disord):
    # given an array of ranges for the disordered region, determine whether given SNPs/amino acid is located within that region.
    for pair in disord:
        if i >= int(pair[0]) and i <= int(pair[1]):
            return True
    return False


cmd.extend("residueDist", residueDist)
