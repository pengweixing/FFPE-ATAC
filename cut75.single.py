#################################################
#  File Name:test.py
#  Author: xingpengwei
#  Mail: xingwei421@qq.com
#  Created Time: Wed 18 Dec 2019 08:49:45 PM UTC
#################################################

import sys
import re
import gzip

f1 = gzip.open(sys.argv[1],'rb') ##input
f2 = gzip.open(sys.argv[2],'wb') ##output

count=1
while 1:
    p1_line = f1.readline()
    if not p1_line:
        break
    if count ==1:
        R1_head = p1_line.decode().strip()
    elif count ==2:
        R1_seq = p1_line.decode().strip()
    elif count ==3:
        R1_qua_h = p1_line.decode().strip()
    elif count ==4:
        R1_qu = p1_line.decode().strip()
        remain_R1 = R1_seq[0:75]
        remain_R1_qu = R1_qu[0:75]
        R1_out = R1_head+'\n'+remain_R1+'\n'+R1_qua_h+'\n'+remain_R1_qu+'\n'
        f2.write(R1_out.encode()) 
    count = count + 1
    if count == 5:
        count = 1
    else:
        count = count
f1.close()
f2.close()
