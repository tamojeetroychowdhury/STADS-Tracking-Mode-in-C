f = open("sm_TM_SNT.txt")

m = f.readlines()
while 1:
    n = int(input("Enter -- "))
    if n==0:
        break
    print(m[n][0])
