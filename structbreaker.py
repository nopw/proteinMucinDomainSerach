import Bio.SeqIO as SeqIO
import os
from datetime import datetime
import sys
if sys.argv.count('drawimg') > 0:
    from domaintopng import domainpositonkind, drawdomain
    print('imgmodule imported')


# this is for static variable,which python do not surport.use classname.property
class domainkind(object):
    # cs complete section
    _domainkindcs = 0
    # ds data section contains more P|T|S
    _domiankindds = 1

    @property
    def kindcs(self):
        return domainkind._domainkindcs

    @property
    def kindds(self):
        return domainkind._domiankindds


class domain(object):
    def __init__(self, adomainkind, adomainseq):
        self.domainkind = adomainkind
        self.domainseq = adomainseq
        self.nextaddr = 0


def getdomainnew(sequence, domainlength, adjustmultiple, sequencelength, tspercent, position):
    # EOF return Nothing
    if position > (sequencelength - 1):
        return None
#    if position > 270:
#       position
    lposition = position
    ldomainkind = domainkind()
    # read exclude section and find the first P|T|S
    excludedomain = domain(0, '')
    while ((lposition < len(sequence)) and (sequence[lposition] != 'P') and (sequence[lposition] != 'T') and (
            sequence[lposition] != 'S')):
        excludedomain.domainseq = excludedomain.domainseq + sequence[lposition]
        lposition = lposition + 1
        excludedomain.nextaddr = lposition
    i = 0
    # read unit
    ldomain = domain(0, '')
    ldomain.domainseq = sequence[lposition: lposition + domainlength]
    domaintspercent = 0.0
    sumdomaints = 0
    for char in ldomain.domainseq:
        if char == 'T' or char == 'S':
            sumdomaints = sumdomaints + 1
    if (None == ldomain.domainseq or '' == ldomain.domainseq):
        domainpercent = 0.0
    else:
        domaintspercent = sumdomaints / len(ldomain.domainseq)
    # if percent of pts in domain less than total pts percent * adjustmultiple, this domain and domain before is cs
    if domaintspercent < tspercent * adjustmultiple:
        ldomain.domainkind = ldomainkind.kindcs
        ldomain.nextaddr = lposition + len(ldomain.domainseq)
    # domain is ds
    else:
        ldomain.domainkind = ldomainkind.kindds
        i = 1
        if (ldomain.domainseq.endswith('P') or ldomain.domainseq.endswith('T') or ldomain.domainseq.endswith('S')) and (
                lposition + len(ldomain.domainseq)) < len(sequence):

            nextword = sequence[lposition + len(ldomain.domainseq)]
            while ((nextword == 'P' or nextword == 'T' or nextword == 'S') and lposition + len(
                    ldomain.domainseq) + 1 < len(sequence)):
                ldomain.domainseq = ldomain.domainseq + nextword
                nextword = sequence[lposition + len(ldomain.domainseq)]
            if (lposition + 1) == len(sequence):
                nextword = sequence[len(sequence) - 1]
                ldomain.domainseq = ldomain.domainseq + nextword
            ldomain.nextaddr = lposition + len(ldomain.domainseq)
        else:
            tmpdomainstr = ldomain.domainseq
            while not (ldomain.domainseq.endswith('P') or ldomain.domainseq.endswith('T') or ldomain.domainseq.endswith(
                    'S')):
                endpositon = len(tmpdomainstr) - i
                ldomain.domainseq = tmpdomainstr[0:  endpositon]
                i = i + 1
            ldomain.nextaddr = position + len(ldomain.domainseq)
            if None != excludedomain and len(excludedomain.domainseq) > 0:
                ldomain.nextaddr = position + len(ldomain.domainseq) + len(excludedomain.domainseq)
    #   print(excludedomain.domainseq)
    #   print(ldomain.domainseq)
    if None != excludedomain and len(excludedomain.domainseq) > 0:
        return [excludedomain, ldomain]
    else:
        return [ldomain]


def splitseqtodomain(sequence):
    oRet = {}

    # position = 0
    adjustmultiple = 1
    defaultlength = 8
    sequencelength = len(sequence.seq)
    sumts = 0
    #ptspercent = 0.60
    sequencestr = str(sequence.seq)
    # get all pts percent in sequence
    for char in str(sequence.seq):
        if char == 'T' or char == 'S':
            sumts = sumts + 1
    tspercent = sumts / sequencelength

    domainlist = getdomainnew(sequencestr, defaultlength, adjustmultiple, sequencelength, tspercent, 0)
    domainindex = 0
    # loop
    while not domainlist is None:
        for domainobj in domainlist:
            if len(domainobj.domainseq) > 0:
                oRet[domainindex] = domainobj
                domainindex = domainindex + 1
        domainlist = getdomainnew(sequencestr, defaultlength, adjustmultiple, sequencelength, tspercent,
                                  domainobj.nextaddr)
    return oRet

def getseqtspercent(sequences):
    sumts = 0
    retValue = 0.0
    if sequences.domainseq == '':
        return 0
    sumts = sequences.domainseq.count('T') + sequences.domainseq.count('S')
    retValue = sumts / len(sequences.domainseq)
#    if retValue > 0:
#        print('tspercent' + str(retValue))
    return retValue


def getppercent(sequences):
    sump = 0
    retValue = 0.0
    if sequences.domainseq == '':
        return 0
    sump = sequences.domainseq.count('P')
    retValue = sump / len(sequences.domainseq)
    return retValue


def getseqptspercent(sequences):
    i = 0
    sumpts = 0
    retValue = 0.0
    if sequences.domainseq == '':
        return 0
    sumpts = sequences.domainseq.count('P') + sequences.domainseq.count('T') + sequences.domainseq.count('S')
    retValue = sumpts / len(sequences.domainseq)
    return retValue


def splicedomain(lastdomain, splitedsequencelist, index):
    i = 0
    newdomain = domain(0, '')
    newdomain.domainseq = lastdomain.domainseq
    newdomain.domainseq = newdomain.domainseq + splitedsequencelist[index].domainseq
    return newdomain


def mergedomainbycondition(splitedsequencelist,  minlength=40):
    tspercent = 0.40
    ppercent = 0.05
    ptspercent = 0.65
    oRet = []
    s = splitedsequencelist
    i = 0
    ldomainkind = domainkind()
    lastmerge = 0
    while (i < len(s)):
        if lastmerge == 1:
            seq = oRet[len(oRet) - 1]
            oRet.remove(oRet[len(oRet) - 1])
        else:
            seq = s[i]
#        if len(seq.domainseq) >= maxlength:
#            oRet.append(seq)
#            i = i + 1
#            lastmerge = 0
#            continue
        # if ptspercent less than min limit
        if getseqtspercent(seq) <= tspercent or getppercent(seq) <= ppercent or getseqptspercent(seq) <= ptspercent:
            # will not merge first cs
            if i == 0 and seq.domainkind == ldomainkind.kindcs:
                seq.domainkind = ldomainkind.kindcs
                oRet.append(seq)
                lastmerge = 0
                i = i + 1
                continue
            # try splice next domain, if ptspercent after splice less than min limit, this domain is cs
            if i < len(splitedsequencelist) - 1:
                newdomain = splicedomain(seq, splitedsequencelist, i)
            else:
                seq.domainkind = ldomainkind.kindcs
                oRet.append(seq)
                i = i + 1
                lastmerge = 0
                continue
            if getseqtspercent(newdomain) <= tspercent or getppercent(newdomain) <= ppercent or getseqptspercent(newdomain) <= ptspercent:
                seq.domainkind = ldomainkind.kindcs
                oRet.append(seq)
                i = i + 1
                lastmerge = 0
                continue
            else:
                if getseqtspercent(newdomain) > tspercent and getppercent(newdomain) > ppercent and getseqptspercent(newdomain) > ptspercent:
                    if len(newdomain.domainseq) >= minlength:
                        newdomain.domainkind = ldomainkind.kindds
                        oRet.append(newdomain)
                        i = i + 1
                        lastmerge = 1
                        continue
                    else:
                        seq.domainkind = ldomainkind.kindcs
                        oRet.append(seq)
                        i = i + 1
                        lastmerge = 0
                        continue
        # else try to splice ds
        else:
            if getseqtspercent(seq) > tspercent and getppercent(seq) > ppercent and getseqptspercent(seq) > ptspercent:
                if i < len(splitedsequencelist) - 1:
                    newdomain = splicedomain(seq, splitedsequencelist, i + 1)
                else:
                    seq.domainkind = ldomainkind.kindcs
                    oRet.append(seq)
                    i = i + 1
                    lastmerge = 0
                    continue
                if getseqtspercent(newdomain) <= tspercent or getppercent(newdomain) <= ppercent or getseqptspercent(newdomain) <= ptspercent:
                    if len(seq.domainseq) > minlength:
                        oRet.append(seq)
                        i = i + 1
                        lastmerge = 0
                        continue
                    else:
                        if getseqtspercent(seq) <= getseqtspercent(newdomain) and getppercent(seq) <= getppercent(newdomain) and getseqptspercent(seq)<= getseqptspercent(newdomain):# and len(newdomain.domainseq) >= minlength:
                            newdomain.domainkind = ldomainkind.kindds
                            oRet.append(newdomain)
                            i = i + 1
                            lastmerge = 1
                            continue
                        else:
                            if len(seq.domainseq) > minlength:
                                seq.domainkind = ldomainkind.kindds
                            else:
                                seq.domainkind = ldomainkind.kindcs
                            oRet.append(seq)
                            i = i + 1
                            lastmerge = 0
                            continue
                else:
                    if getseqtspercent(newdomain) > tspercent and getppercent(newdomain) > ppercent and getseqptspercent(newdomain) > ptspercent:# and len(newdomain.domainseq) >= minlength:
                        newdomain.domainkind = ldomainkind.kindds
                        oRet.append(newdomain)
                        lastmerge = 1
                        i = i + 1
                    else:
                        seq.domainkind = ldomainkind.kindcs
                        oRet.append(seq)
                        i = i + 1
                        lastmerge = 0
                        continue
    return oRet


def getmarkedds(adomain):
    i = 0
    s = adomain.domainseq
    pts = 'PTS'
    sRet = ''
    while i < len(s):
        if s[i] in pts:
            sRet = sRet + r'<span id ="MYPTS">' + s[i] + r'</span>'
        else:
            sRet = sRet + s[i]
        i = i + 1
    return sRet


def save_to_file(file_name, contents):
    fh = open(file_name, 'w')
    fh.write(contents)
    fh.close()


class summaryinfo:
    def __init__(self, sequence):
        pts='PTS'
        ts='TS'
        self.ptscount = 0
        self.tscount=0
        i = 0
        s = sequence
        while i < len(s):
            if s[i] in pts:
                self.ptscount = self.ptscount + 1
            if s[i] in ts:
                self.tscount = self.tscount + 1
            i = i + 1
        self.ptspercent = self.ptscount / (i)
        self.tspercent = self.tscount / (i)
        self.count = i


path = os.path.abspath(os.curdir) + '\\'
print(path)
dbs = []
ldmkd = domainkind()
files = os.listdir(path)
for f in files:
    if (str(f).endswith('fasta') and f.count('possible')==0):
        dbs.append(f)
if (None == dbs or dbs.count == 0):
    print('0 database found. program ends. press enter to exit')
else:
    print(dbs)
    totalseqcount = 0
for db in dbs:
    filehandle = path + db
    outputfilecontent = r'<!DOCTYPE html><html><head><style>#dsbackground{background-color:yellow;font-weight:bold;}#MYPTS{backgroud-color:#000000; color :blue;font-weight:bold;} .tooltip:hover:after { content: attr(data-tooltip);} .nowrap{white-space:nowrap;}</style></head><body>'
    num = 1
    oldid = ''
    fastaoutput=[]
    for seq in SeqIO.parse(filehandle, 'fasta'):
        #print(seq.id)
        start = datetime.now()
        dropseq = 1
        seqdiv = ''
        #skip same seq
        newid = str(seq.id).split('-')[0]
        if None != oldid:
            if oldid == newid:
                continue
        oldid = newid
        loopsummaryinfo = summaryinfo(str(seq.seq))

        ltooltipcontent = '[sum(PTS):{0}, sequence length:{1},percent of PTS in sequence:{2},percent of TS in sequence:{3}]'.format(
            loopsummaryinfo.ptscount, loopsummaryinfo.count, loopsummaryinfo.ptspercent, loopsummaryinfo.tspercent)
        seqdiv = seqdiv + r'<div><p style="color:red">No.:' + str(
            num) + ' Sequence id:' + seq.id + r'</p></div>' + ltooltipcontent + '<div></div><div><p class="nowrap">'
        spliteddomainlist = splitseqtodomain(seq)
        sssss = ''
        tmp = 0
        for tmp in range(0, len(spliteddomainlist)):
            sssss = sssss + spliteddomainlist[tmp].domainseq
        if sssss != str(seq.seq):
            print(sssss)
            print(str(seq.seq))
#        if seq.id == 'tr|A0A182G5J9|A0A182G5J9_AEDAL':
#            c = 1
        mergedomainlist = mergedomainbycondition(spliteddomainlist)
        i = 0
        positiondomainstart = 0
        positiondomainend = 1
        for i in range(0, len(mergedomainlist)):
            ldomain = mergedomainlist[i]
            positiondomainstart = positiondomainend
            positiondomainend = positiondomainstart + len(ldomain.domainseq)

            if ldomain.domainkind == ldmkd.kindcs:
                seqdiv = seqdiv + getmarkedds(ldomain)
            else:
                ssummaryinfo = summaryinfo(ldomain.domainseq)
                tooltipcontent = '[sum(PTS):{0}, domain length:{1},percent of PTS in domain:{2},domainposition[{3}:{4}]]'.format(
                    ssummaryinfo.ptscount, ssummaryinfo.count, ssummaryinfo.ptspercent, positiondomainstart,
                    positiondomainend - 1)
                seqdiv = seqdiv + r'<span data-tooltip="' + tooltipcontent + '" class="tooltip"><span id="dsbackground">' + getmarkedds(
                    ldomain) + r'</span></span>'
                dropseq = 0
        seqdiv = seqdiv + r'</p></div>'
        if dropseq == 0:
            outputfilecontent = outputfilecontent + seqdiv
            fastaoutput.append(seq)
            num = num + 1
        finish = datetime.now()
    if fastaoutput:
        SeqIO.write(fastaoutput, path + '\\possible_mucin_' + db, "fasta")
    outputfilecontent = r'<div><p>start:' + start.strftime(
        '%Y-%m-%d %X') + '</p></div>' + outputfilecontent + r'<div><p>finish:' + finish.strftime(
        '%Y-%m-%d %X') + '</p></div>' + r'</body></html>'
    save_to_file(path + db.replace('fasta', 'html'), outputfilecontent)


