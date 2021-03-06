import os


def chisquare():
    sub_dir = "/home/andrew/Documents/CS/projects/andrewarul/"
    count = 1
    ddist = []
    odist = []
    degrees = 0

    for folder in os.listdir(sub_dir):
        for fl in os.listdir(os.path.join(sub_dir, folder)):
            if count == 1:
                f = open(os.path.join(sub_dir, folder, fl), "r")
                line1 = f.readline()
                check = line1[25:28]
                if check != 'N/A':
                    ddist.append(line1[25:len(line1) - 2])
                    line2 = f.readline()
                    odist.append(line2[22:len(line2)])
                    degrees += 1
                f.close()
            if count == 3:
                count = 1
            else:
                count += 1
    ddist.remove('')
    odist.remove('')
    nddist = []
    nodist = []
    for i in ddist:
        nddist.append(float(i))
    for j in odist:
        nodist.append(float(j))
    chisquared = 0
    for i in xrange(0, len(ddist)):
        chisquared += ((nddist[i] - nodist[i]) ** 2) / (nodist[i])
    return chisquared, degrees
    if not os.path.exists("/home/andrew/Documents/CS/projects/andrewarul/chisquared_results.txt"):
        os.makedirs("/home/andrew/Documents/CS/projects/andrewarul/chisquared_results.txt")
    f = open("/home/andrew/Documents/CS/projects/andrewarul/chisquared_results.txt", 'w')
    f.write("Chi-squared value: " + str(chisquared) + "\n" + "Degrees of freedom: " + str(degrees - 1) + "\n")
    f.close()
    print("Statistical analysis completed.")
