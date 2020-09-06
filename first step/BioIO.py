import Bio.SeqIO as SeqIO
import json
import os


#print(Bio.__version__)a


path = os.path.split(os.path.realpath(__file__))[0]
print(r'当前路径'+path)
files= os.listdir(path)
databases=[]
for f in files:
    if f.count('fasta')> 0 and f.count('result')==0:
       databases.append(f)
if (not databases):# or databases.count ==0 or em):
   input("未找到目标数据,按回车退出")
   os._exit(0)
print(r'目标数据库'+ str(databases))
resultcount = 0
for db in databases:
    results = []
    compare = []
    for seq_record in SeqIO.parse(path+'\\'+ db ,"fasta"):
        seqlen =len(seq_record.seq)
        if seqlen <= 100:
            continue
        s = 0
        s = seq_record.seq.count('S')
        t = 0
        t = seq_record.seq.count('T')
        p = 0
        p = seq_record.seq.count('P')
        stl = (s+t)/seqlen
        pl = p/seqlen
        if (s + t)/seqlen > 0.15: #and (p/seqlen) > 0.05:
            resultcount += 1
            line = {r'SeqID':seq_record.id,
                    r'S':s,
                    r'T':t,
               #     r'P':p,
                    r'seqlen':seqlen,
                    r"S+T/seqlen":stl
               #     r'P/seqlen/1000':pl
                    }
            print(line)
            compare.append(line);
            results.append(seq_record)
    if (None != results ):
      SeqIO.write(results,path+'\\result'+ db,"fasta")
      fileObject = open(path+ "\\compare"+ db.split('.')[0]+'.txt', 'w')
      jsObj = json.dumps(compare)
      for aline in jsObj:
        fileObject.write(aline)
      fileObject.close()
print('total seq count:'+str( resultcount))
input("操作完成,按回车退出")