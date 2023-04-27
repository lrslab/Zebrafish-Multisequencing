import argparse
import numpy as np
import math

parser = argparse.ArgumentParser(description="Calculate the estimated accuracy, error and quality value of each read!")
parser.add_argument("-in", "--input", type=str, metavar="", required=True, help="input the fastq file")
args = parser.parse_args()

"""
input a quality string and this function can help to 
calculate the mean error ratio. eg. %%$$$$%&'&%&&&&'((&7,
"""
def calculation_read_acc(string):
    output_list = []
    error_list = []

    for base_value in string:
        ascII = ord(base_value) -33
        error_proporation = math.pow(10, (-1)*int(ascII)/10)
        error_list.append(error_proporation)
    
    error_mean = np.mean(error_list)
    output_list.append(1-error_mean)
    output_list.append(error_mean)
    output_list.append((-10) * math.log10(error_mean))
    return(output_list)

def main(input_file):
    print("ID\testimated_acc\testimated_err\tQ_value")
    count = 1 # for getting the quality line
    with open(input_file) as File:
        for line in File:
            if count % 4 == 1:
                line = line.replace("\n", "")
                line = line.split(" ")
                ID = line[0]
                count += 1

            elif count % 4 == 0:
                line = line.replace("\n", "")
                quality = calculation_read_acc(line)
                print(ID + "\t" + str(quality[0]) + "\t" + str(quality[1]) + "\t" + str(quality[2]))
                count += 1

            else:
                count += 1

if __name__=="__main__":
    main(args.input)
