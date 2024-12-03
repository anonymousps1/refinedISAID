import sys

#load(sys.argv[1])
#n = int(sys.argv[1])
r = int(sys.argv[1])
n = 4096


def noStar(row):
	for e in row:
		if e == -1:
			return False
	return True


def getOneIndex(row):
	for i in range(n):
		if row[i] == 1:
			return i

	print('Error')

def getmatrix():
	matrix1 = [[0 for i in range(4096)] for j in range(n)] 
	with open( str(r)+"r_matrix.txt", "r") as rd:
		for i in range(n):
			line = rd.readline()
			for j in range(4096):
				if line[j] == "*":
					matrix1[i][j] = -1
				else:
					matrix1[i][j] = int(line[j])
	return matrix1



def main():
	arr = getmatrix()
	for i in range(n):
		arr[i][i] = 1
	listNoStars = []
	unitVectorsIndices = set()

	for elem in arr:
		if noStar(elem):
			listNoStars.append(elem)

	mat = matrix(GF(2), listNoStars)
	print(mat.dimensions(), len(arr))

	while mat.rank() < n:
		oldRank = mat.rank()

		for row in mat.echelon_form():
			if row.hamming_weight() == 1:

				unitVectorsIndices.add(getOneIndex(row))

		for elem in arr:
			for index in unitVectorsIndices:
				if elem[index] == -1:
					elem[index] = 0


		listNoStars = []
		for elem in arr:
			if noStar(elem):
				listNoStars.append(elem)

		mat = matrix(GF(2), listNoStars)
		print(mat.dimensions())
		if mat.rank() == oldRank:

			print('\nRank {}'.format(mat.rank()))
			return
		else:
			oldRank = mat.rank()

	print('Rank {}'.format(mat.rank()))


if __name__ == '__main__':
	main()
