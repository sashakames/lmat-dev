import sys


amt = 2 ** int(sys.argv[1])

arr = [0] * amt

def div(a, b):
     return a / b

def mod(a, b):
     return a % b


for fn in sys.stdin:
     inf = open(fn.strip())
     for line in inf:
        if line[0] != '#':
          parts = line.split()
          val = int(parts[0])

	
          arr[mod(val, amt)] += 1
     inf.close()


for a in arr:
     print a
