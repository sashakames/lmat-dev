import sys

gid_to_tid_fn = sys.argv[1]
kpath_fn = sys.argv[2]

gid = int(sys.argv[3])
tid = int(sys.argv[4])
count = int(sys.argv[5])
parent = sys.argv[6]

gt_f = open(gid_to_tid_fn, "a")
kp_f = open(kpath_fn, "a")

ids_str = " "

for i in range(1, count+1):
    gid = gid + 1
    tid = tid + 1

    gt_f.write(str(gid)  + " " + str(tid) + "\n")
    
    kp_f.write(str(tid) + " 0 " + parent + "\n")
    kp_f.write("Homo sapiens seq. " + str(i) + "\n")

    ids_str = ids_str + " " + str(tid)

print ids_str


