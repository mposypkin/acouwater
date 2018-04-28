import os
print('Hello')
for d in range(0,1):
    for i in range(1,6):
        for j in range(1,6):
            f = 0 + i * 0.2
            g = 1 + j * 0.4
            fname = 'out_' + str(f) + '_' + str(g) + '_' + str(d)
            s = './acourosen.exe ' + str(f) + ' ' + str(g) + ' ' + str(d)
            es = 'echo ' + s
            fs = s + ' > ' + fname
            os.system(es)
            os.system(fs)
