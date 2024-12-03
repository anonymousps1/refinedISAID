
import os
import re
import sys


def genkey():

    first = [[0 for i in range (64)] for j in range(64)]
    last = [[0 for i in range(64)] for j in range(64)]

    for i in range(16):
        last[4*i][2+4*i] = 1
        last[1+4*i][3+4*i] = 1
        last[2+4*i][0+4*i] = 1
        last[2+4*i][3+4*i] = 1
        last[3+4*i][0+4*i] = 1
        last[3+4*i][2+4*i] = 1
    for i in range(16):
        first[4*i][0+4*i] = 1
        first[4*i][2+4*i] = 1
        first[1+4*i][2+4*i] = 1
        first[2+4*i][1+4*i] = 1
        first[3+4*i][0+4*i] = 1
        first[3+4*i][1+4*i] = 1
    first1 = [[0 for i in range(64)] for j in range(64)]
    P = [0, 16, 32, 48, 1, 17, 33, 49, 2, 18, 34, 50, 3, 19, 35, 51, 4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55, 8, 24, 40, 56, 9, 25, 41, 57, 10, 26, 42, 58, 11, 27, 43, 59, 12, 28, 44, 60, 13, 29, 45, 61, 14, 30, 46, 62, 15, 31, 47, 63]

    for i in range(64):
        for j in range(64):
            first1[i][P[j]] = first[i][j]

    for i in range(64):
        print(first1[i])

    print()
    for i in range(64):
        print(last[i])

    
    tup1 = (4, 20, 36, 52, 5, 21, 37, 53, 6, 22, 38, 54, 7, 23, 39, 55)
    for i in range(48, 64):
        for j in range(16):
            with open("./12round_midkey/k"+str(i)+"__"+str(tup1[j])+".txt", "r")as rd:
                line = rd.readline()
                key_middle = ""
                for k in range(9):
                    key_middle += rd.readline()
                for k1 in range(4):
                    for k2 in range(4):
                        with open("./10round_allkey/k"+str((i-48)*4+k1)+"_"+str(j*4+k2)+".txt", "w") as wd:
                            for index in range(64):
                                wd.write(str(first1[(i-48)*4+k1][index]))
                            #wd.write(first1[(i-16)*4+k1])
                            wd.write("\n")
                            wd.write(key_middle)
                            for index in range(64):
                                wd.write(str(last[j*4+k2][index]))
                            wd.write("\n")
    
if __name__ == "__main__":
    genkey()
    pass
