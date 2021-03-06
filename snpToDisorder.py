"""
@ Andrew Li and Arul Prasad, PRIMES 2014
@ Script purpose: Given several objects, lengths, the location of drug resistant SNPs within the proteins, and arrays of ranges of disorder within the protein, the program will output two textfiles of distances from the SNPs to every other amino acid in each protein.
This was used to show the correlation between SNPs and disordered regions in cancer related proteins.
"""

from pymol import cmd
import re  # for splitting
import xlrd  # manipulating xlsx files
import os
import numpy as np
import matplotlib.pyplot as plt
import pylab
from chisquare import chisquare


def snpToDisorder():

    n_avglist = []
    d_avglist = []

    if not os.path.exists("/home/andrew/Documents/CS/projects/andrewarul/histograms"):
        os.makedirs("/home/andrew/Documents/CS/projects/andrewarul/histograms")

    workbook = xlrd.open_workbook("/home/andrew/Documents/CS/projects/andrewarul/bacteria.xlsx")
    worksheet = workbook.sheet_by_name("Sheet1")
    num_rows = worksheet.nrows
    for cur_row in range(1, num_rows):
        obj = worksheet.cell_value(cur_row, 1)
        if obj != "":
            cmd.load("/home/andrew/Documents/CS/projects/andrewarul/" + obj + ".pdb")
        else:
            continue
        length = int(worksheet.cell_value(cur_row, 2))
        start = int(worksheet.cell_value(cur_row, 3))

        s_input = worksheet.cell_value(cur_row, 4)
        snps = [int(s) for s in s_input.split(", ") if s.isdigit()]
        print snps
        # remove duplicate values
        d_input = worksheet.cell_value(cur_row, 5).split()
        d_list = []
        for i in range(1, len(d_input)):
            # gets the two numbers that denote the range of the disordered region
            d_range = re.split(r"[\n,\s]", d_input[i])
            for e in d_range:
                if e == "":
                    d_range.remove(e)
            d_temp = str(d_range[0])
            d_temp_array = re.findall(r"\d+", d_temp)
            d_list.append(d_temp_array)
        print "-----Analyzing " + obj + ".-----"
        if obj == "" or length == "" or start == "" or s_input == "" or d_input == "":
            print "Protein skipped, incomplete information."
            continue
        else:
            print "Data successfully obtained."

        # check if chain is not A
        ch = "a"
        if worksheet.cell_value(cur_row, 6) != "":
            ch = worksheet.cell_value(cur_row, 6)
        # create folder for obj
        if not os.path.exists("/home/andrew/Documents/CS/projects/andrewarul/" + obj):
            os.makedirs("/home/andrew/Documents/CS/projects/andrewarul/" + obj)
        fout = open("/home/andrew/Documents/CS/projects/andrewarul/" + obj + "/n_results.txt", "w")  # normal results
        fout_d = open("/home/andrew/Documents/CS/projects/andrewarul/" + obj + "/d_results.txt", "w")  # disordered results

        # used for stat - numpy/matplotlib
        fstat = open("/home/andrew/Documents/CS/projects/andrewarul/" + obj + "/avgs.txt", "w")
        n_sum = 0  # sums and counts of all distances, for avg.
        d_sum = 0
        n_count = 0
        d_count = 0
        n_hist = []  # values to be placed in histograms
        d_hist = []

        for i in range(len(snps)):
            for j in range(start, length + start - 1):
                if snps[i] != j:
                    temp_dist = calcDist(obj, snps[i], j, ch)
                    # if temp_dist is -1, it means that the residue could not be found
                    if temp_dist == -1:
                        pass
                    elif (inDisorder(j - start + 1, d_list)):
                        fout_d.write(str(temp_dist) + "\n")
                        d_sum += temp_dist
                        d_count += 1
                        d_hist.append(int(temp_dist))
                        # print str(snps[i]) + ", " + str(j) + ", " + str(temp_dist)
                    else:
                        fout.write(str(temp_dist) + "\n")
                        n_sum += temp_dist
                        n_count += 1
                        n_hist.append(int(temp_dist))
                else:
                    pass
        fstat.write("Avg Disordered Distance: ")
        if d_count != 0:
            fstat.write(str(d_sum / d_count) + "\n")
        else:
            fstat.write("N/A")
        fstat.write("Avg Ordered Distance: ")
        if n_count != 0:
            fstat.write(str(n_sum / n_count))
        else:
            fstat.write("N/A")
        if d_count != 0 and n_count != 0:
            n_avglist.append(n_sum / n_count)
            d_avglist.append(d_sum / d_count)

        fout.close()
        fout_d.close()
        print obj + " has been analyzed."
        cmd.delete(obj)

        # generate histograms
        if len(n_hist) > 0 and len(d_hist) > 0:
            bins = np.linspace(0, 100, 100)
            plt.xlabel("Distances (angstroms)")
            plt.ylabel("Frequency of Distance")
            plt.hist(n_hist, bins, color="green", alpha=0.5, label="Ordered")
            plt.hist(d_hist, bins, color="blue", alpha=0.6, label="Disordered")
            plt.legend(loc="upper right", numpoints=1)
            pylab.savefig("/home/andrew/Documents/CS/projects/andrewarul/" + obj + "/" + obj + "_hist.png")
            pylab.savefig("/home/andrew/Documents/CS/projects/andrewarul/histograms/" + obj + "_hist.png")
            print "Histogram generated.\n"
            plt.cla()
            plt.clf()
        else:
            print "Failed to generate histogram.\n"
    print "reached"
    bins = np.linspace(0, 100, 100)
    plt.xlabel("Distances (angstroms)")
    plt.ylabel("Frequency of Distance")
    plt.hist(n_avglist, bins, color="green", alpha=0.5, label="Ordered Averages")
    plt.hist(d_avglist, bins, color="blue", alpha=0.5, label="Disordered Averages")
    plt.legend(loc="upper right", numpoints=1)
    pylab.savefig("/home/andrew/Documents/CS/projects/andrewarul/histograms/avg_hist.png")
    print "\n\nGENERAL HISTOGRAM GENERATED.\n"

    # run chisquare program (in same folder)
    # chisquare()

    balance = 0
    for i in range(len(n_avglist)):
        if n_avglist[i] > d_avglist[i]:
            pass
        else:
            balance += 1
    print balance


def calcDist(obj, l1, l2, ch):
    # obj is the protein in question, l1, l2 are the positions of residues.
    # delete previous selections and displayed distances
    cmd.delete("dist")
    cmd.delete("sele")
    cmd.delete("Residue_1")
    cmd.delete("Residue_2")

    # select two residues based on input, at positions l1 and l2 in the sequence
    cmd.select("Residue_1", "resi {} and chain {}".format(l1, ch))
    cmd.select("Residue_2", "resi {} and chain {}".format(l2, ch))

    #the underscores ARE IMPORTANT
    distance = cmd.distance("dist", "Residue_1", "Residue_2")
    # print distance
    return distance


def inDisorder(i, disord):
    # given an array of ranges for the disordered region, determine whether given SNPs/amino acid is located within that region.
    for pair in disord:
        if i >= int(pair[0]) and i <= int(pair[1]):
            return True
    return False


def stat(lst):
    mean = 0
    for i in lst:
        mean += i
    mean /= len(lst)
    stddev = 0
    for i in lst:
        stddev += (mean - i) ** 2
    stddev /= (len(lst) - 1)
    return [mean, stddev ** 2]


cmd.extend("snpToDisorder", snpToDisorder)
