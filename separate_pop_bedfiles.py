from __future__ import print_function
import sys
import gzip


popfile = open(sys.argv[1], 'r')

EUR_list = []
ASN_list = []
for line in popfile:
    if "samp" in line:
        continue
    ln = line.strip().split('\t')
    msp_ID = ln[0]
    pop = ln[1]
    #print(msp_ID,file=sys.stdout)
    if pop == 'EUR':
        EUR_list.append(msp_ID)
    elif pop == 'ASN' or pop == 'EAS':
        ASN_list.append(msp_ID)
popfile.close()

#print(ASN_list, file=sys.stdout)


bedfile = gzip.open(sys.argv[2], 'rb')

fname = str(sys.argv[2]).strip('.gz')

EUR_bedfile=gzip.open(fname+'_EUR.gz','wb')
ASN_bedfile=gzip.open(fname+'_ASN.gz','wb')
print(fname)
for line in bedfile:
#    print(line.strip())
    ln = line.strip().split('\t')
    msp_hap_chr = ln[0]
    msp_ID = msp_hap_chr.split(':')[0]
#    print(msp_ID)
    if msp_ID in EUR_list:
        EUR_bedfile.write(line)
    elif msp_ID in ASN_list:
        ASN_bedfile.write(line)

bedfile.close()
EUR_bedfile.close()
ASN_bedfile.close()
