#chr1    6636590 C       1       193     37.21   [1, 191, 0, 1]  CT CA   0.01    -       -       -       -       -
#chr1    6636590 C       1       193     37.21   [1, 191, 0, 1]  CA CT   0.01    -       -       -       -       -

import sys
from itertools import izip_longest

file1 = sys.argv[1]
file2 = sys.argv[2]

fd1 = open(file1, "r")
fd2 = open(file2, "r")

for left, right in izip_longest(fd1, fd2, fillvalue="-"):
	not_equal = False

	left = left.strip()
	right = right.strip()
	fields1 = left.split("\t")
	fields2 = right.split("\t")

	if len(fields1) == len(fields2):
		if len(fields1) > 7:
			if set(fields1[7].split()) != set(fields2[7].split()): not_equal = True
			del fields1[7]
			del fields2[7]
		if fields1 != fields2: not_equal = True
	else:
		not_equal = True

	if not_equal:
		print("> " + left)
		print("< " + right)

fd1.close()
fd2.close()
