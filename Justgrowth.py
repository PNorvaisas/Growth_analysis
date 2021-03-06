#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Justgrowth.py

Created by Povilas Norvaisas on 2015-02-26.
Copyright (c) 2015. All rights reserved.

"""

install=False

try:
    try:
        from pip import main as pipmain
    except:
        from pip._internal import main as pipmain
        
except ImportError, e:
    print "Module pip not found!"
    print "Please install pip manually to proceed!"
    sys.exit(1)
    
def pipinstall(package):
    pipmain(['install', package])

    
for mod in ['pip','string','math','re','csv','sys','os',
            'commands','datetime','operator','getopt','subprocess','shutil','glob',
            'types','math','copy','pyExcelerator','xlrd','xlwt','xlutils','types','warnings']:
    try:
        exec "import %(mod)s" % vars()
    except ImportError, e:
        print "Module not found %(mod)s" % vars()
        
        qstn = "Python module %(mod)s not found, should it be installed? (yes/no): " % vars()
        answ = raw_input(qstn)
        
        if answ in ['y','yes']:
            install=True
        else:
            print "Script cannot function without module %(mod)s, quitting!" % vars()
            sys.exit(1)
        
        print "\nTrying to install!"
        if install:
            pipinstall(mod)        
        #pass # module doesn't exist, deal with it.


#import warnings
#Disable simple warnings
warnings.filterwarnings("ignore")


import unicodedata

import numpy as np
import itertools as IT
import pandas as pd

from matplotlib import pyplot as plt
from scipy import interpolate as ip
from scipy import signal as sig
from scipy.optimize import curve_fit


help_message = '''
Bacterial growth data preparation

Flags:
    -h  Display this help message
    -v  Verbose output
Arguments:
    -i <files>   Input file
    -o <dir>     Directory to write the output to
Options:
    full         Return additional figures and data columns
'''

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def usage():
    print help_message % vars()
    return


optionsset='''
Options:
<--------------------------------------------->
      Files:    %(ifile)s
    Pattern:    %(pfile)s
        Out:    %(odir)s
       Full:    %(full)s
<--------------------------------------------->
    '''

def main(argv=None):
    ifile=""
    pfile=""
    msize=20
    odir='Output'
    odirn=''
    load=False
    full=False
    subst=''
    comparison=False
    
    if argv is None:
        argv = sys.argv
    try:
        try:
            opts, args = getopt.getopt(argv[1:], "hi:f:d:o:p:", ["help"])
        except getopt.error, msg:
            raise Usage(msg)


        for option, value in opts:
            if option in ("-h", "--help"):
                usage()
                return    
            if option in ("-i", "--input"):
                ifile=value
            if option in ("-p", "--pattern"):
                pfile=value
            if option in ("-o", "--out"):
                odir=value
    
    
        for argument in args:        
            if argument in ("full", "--full"):
                full = True

    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2

    #Check the input integrity
    
 

    try:
        if os.path.isfile(ifile):
            ipath, iname, itype = filename(ifile)
            info=readinfo(ifile)
            ilist=info['File'].values
            dlist=info['Pattern'].values
            udlist=list(set(dlist))
            descriptors=readdesc(udlist)
        else:
            raise Exception("No design file specified!")

    except Exception, e:
        print e
        #print "!!!-------------------------------------------!!!"
        sys.exit(1)

    print '''\n\n---------------------------\n\n'''


    print optionsset %vars()
    #----------------------------
    
#     if 'Design' in ifile:
#         print 'Reading experimental design information!'
#         info,ilist,dlist=readinfo(ifile)
        
#         udlist=list(set(dlist))
#         subs=[[]]*len(ilist)
#         descriptors=readdesc(udlist)
#     else:
#         info={}
#         print 'Reading data files!'
#         print ifile
#         ilist=genlist(ifile)
#         #print ilist
#         if pfile!='':
#             dlist=genlist(pfile)
#             udlist=list(set(dlist))
#             descriptors=readdesc(udlist)
#         else:
#             descriptors={}

    #Support absolute links?
    odir=dircheck(odir)

    data=collect(info)
    
    #Additional step to check if there are missing values?
    
    data=analyse(data,full)

    sheets=makesheets(data,descriptors,info)
    
    writesheets(sheets,odir)

    plot_comparison(data,odir,'all')

#-------------Functions------------




def runcmd(cmd):
    failure, output = commands.getstatusoutput(cmd)
    if failure:
        print '''Running failed \n %(cmd)s \n %(output)s'''.encode('utf-8') % vars();
    return failure, output

class NestedDict(dict):
    def __getitem__(self, key):         
        if key in self: return self.get(key)
        return self.setdefault(key, NestedDict())

def numerize(s):
    try:
        if s=='NAN':
            return s
        float(s)
        if float(s).is_integer():
            return int(float(s))
        elif float(s)==0:
            return float(s)
        else:
            return float(s)
    except ValueError:
        return s

def round_to(n, precission):
    #Round a number to desired precision
    correction = 0.5 if n >= 0 else -0.5
    return int(n/precission+correction)*precission
    
def filename(ifile):
    if ifile.split('.')[0]=='':
        ipat=''
        iname=''
        itype=ifile.split('.')[1]
    else:
        if "\\" in ifile.split('.')[0]:
            sep="\\"
        elif "/" in ifile.split('.')[0]:
            sep="/"
        else:
            ipat=''
            iname=ifile.split('.')[0]
            itype=ifile.split('.')[1]
            return ipat, iname, itype
        allpath=ifile.split('.')[0]    
        iname=allpath.split(sep)[-1]
        ipath=allpath.split(sep)[:-1]
        ipat='/'.join(ipath)
        itype=ifile.split('.')[1]
    return ipat, iname, itype

def dircheck(somedir):
    while True:       
        if os.path.exists(somedir):
            qstn = "Directory %(somedir)s already exists! Delete, quit, continue or provide a new name (d/q/c/<type name>): " % vars()
            answ = raw_input(qstn)
            if answ == "d":
                shutil.rmtree(somedir)
                os.makedirs(somedir)
                break
            elif answ == "q":
                sys.exit("Have a nice day!")
            elif answ == "c":
                break
            else:
                somedir=answ
                continue
        else:
            os.makedirs(somedir)
            break
    return somedir


def readxls(ifile):
    book=xlrd.open_workbook(ifile,formatting_info=False)
    data=book.sheet_by_name([nm for nm in book.sheet_names() if 'Magellan' in nm][0])
    sheet=[]    
    for r in range(0,data.nrows):
        if len(data.row_values(r))>200:
            #print set(data.row_values(r))
            sheet.append(data.row_values(r))
    return sheet


def readtext(ifile):
    f=open(ifile,'r')
    sheet=[]
    for l in f:
        row=[cell.strip() for cell in l.replace('\t',',').split(',')] #re.split(r"[,\t]+",l.replace('\t',','))
        row=[numerize(cell) for cell in row]
        sheet.append(row)
    f.close()
    return sheet

def readxls_s(ifile):
    book=xlrd.open_workbook(ifile,formatting_info=False)
    data=book.sheet_by_name(book.sheet_names()[0]) #Use first sheet
    sheet=[]
    for r in range(0,data.nrows):
        #print set(data.row_values(r))
        sheet.append(data.row_values(r))
    return sheet

def readcsv(ifile):
    f=open(ifile,'r')
    sheet=[]
    rdr=csv.reader(f, delimiter=',')
    data=[ln for ln in rdr]
    f.close()
    return data

def readdesc(udlist):
    
    ddata=[]
    for din,dfile in enumerate(udlist):
        
        
        try:
            book=xlrd.open_workbook(dfile,formatting_info=False)
        except Exception, e:
            raise Exception('Cannot open file: {}'.format(dfile))
            
        variables=book.sheet_names()
        
        varlist=[]
        for var in variables:
            sheet=book.sheet_by_name(var)
            columns=sheet.row_values(0)[1:]
            rows=sheet.col_values(0)[1:]

            #Skip header
            vardata=[]
            for r in range(1,sheet.nrows):
                row=sheet.row_values(r)
                vardata.append(row)
            
            #Prepare variable data
            varDF=pd.DataFrame(vardata,columns=['Row']+columns)
            varDFm=pd.melt(varDF,id_vars=['Row'],value_vars=columns,var_name='Col',value_name=var)
            varDFm['Well']=varDFm['Row']+varDFm['Col'].astype(int).astype(str)
            varDFm=varDFm.set_index('Well')
            varlist.append(varDFm[var])
        
        descDF=pd.concat(varlist,axis=1)
        descDF['Pattern']=dfile
        ddata.append(descDF)

    descriptors=pd.concat(ddata)
    descriptors['Well']=descriptors.index
                    
    return descriptors


def readinfo(ifile):
    print ifile
    info=NestedDict()
    ilist=[]
    dlist=[]
    odirs=[]
    ipt, inm, itp = filename(ifile)
    
    if itp in ['xlsx','xls']:
        data=readxls_s(ifile)
    elif itp=='csv':
        data=readcsv(ifile)
    headers=data[0]
    
    #Automatically find variables
    headin={ hd : headers.index(hd) for hd in headers}
    nec=['File','Pattern']

    addhead=[key for key in headin.keys() if key not in nec]
    
    miss=[hd for hd in nec if hd not in headers]
    
    if all(n in headers for n in nec):
        print 'Necessary headers found!'
    else:
        raise Exception('Missing essential headers in description file!\n{}'.format(miss))
     
    info=pd.DataFrame(data[1:],columns=headers)
    
    info=info[info['File']!='']
    
    #info=info.set_index('File')

    return info



def myround(a, decimals=1):
     return np.around(a-10**(-(decimals+5)), decimals=decimals)

# def growth(x,a,c):
#     y=x*a+c
#     return y

def Wiener(y, n):
    wi = sig.wiener(y, mysize=n)

    if sum(np.isnan(wi))>0:
        return y
        
    return wi



def growth(x,A,lam,u):
    return A/(1+np.exp((4*u/A)*(lam-x)+2))

def log_growth(x,A,lam,u):
    y=np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))
    #print x,y
    return np.log2(A)-np.log2(1+np.exp((4*u/A)*(lam-x)+2))

def lin(x,a,c):
    y=x*a+c
    return y


def setbar(x,bar):
    x2=[xi if xi>bar else bar for xi in x]
    x2=np.array(x2)
    return x2


def interp(x,y,x2):
    tck = ip.splrep(x, y, s=0)
    y2=ip.splev(x2, tck, der=0)
    return y2

def cut(x, y, a,b,equalize=False):
    x2=[]
    y2=[]
    maxy=max(y)
    maxind=y.tolist().index(maxy)
    last=0
    for xt, yt in IT.izip(enumerate(x),enumerate(y)):
        df=yt[0]-last
        #print df
        if yt[1]>a and yt[1]<=b:# and df<3
            x2.append(xt[1])
            last=yt[0]
            if yt[0]<maxind:
                y2.append(yt[1])
            else:
                if equalize:
                    y2.append(maxy)
                else:
                    y2.append(yt[1])
    y2=np.asarray(y2)
    x2=np.asarray(x2)

    return x2,y2




def getallwells():
    r = [l for l in string.ascii_uppercase[0:16]]
    c = [str(i) for i in range(1,25)]
    allwells=[]
    for rn in r:
        for cn in c:
            #print rn+cn
            allwells.append(rn+cn)

    return allwells


def time_to_sec(tstr):
    h,m,s=tstr.split(':')
    seconds=int(s)+60*int(m)+3600*int(h)
    return seconds


# def genlist(ifile):
#     #Generates list of input files, checks their existance
#     ilist=[]
#     if ',' in ifile:
#         ifiles=ifile.split(',')
#         for ifl in ifiles:
#             ilist.extend(genlist(ifl))
#     else:
#         ipath, iname, itype=filename(ifile)
#         if itype in ['xls','xlsx','asc','txt'] and os.path.isfile(ifile):
#             ilist.append(ifile)
#         # elif itype in [''] and os.path.isfile(ifile):#'txt',
#         #     ifl=open(ifile,'r')
#         #     idata=ifl.read().split('\n')
#         #     idata=[fl.strip() for fl in idata if fl!='']
#         #     for fld in idata:
#         #         ilist.extend(genlist(fld))

#         elif iname=='' and itype in ['xls','xlsx','asc','txt']:
#             if itype in ['xls','xlsx']:
#                 ffiles=glob.glob('*.%(itype)s' % vars())
#                 #print ffiles
#                 ilist.extend(ffiles)
#             elif itype=='txt':
#                 for tfile in glob.glob('*.%(itype)s' % vars()):
#                     ilist.extend(genlist(tfile))
#             else:
#                 print "Bad file type %(inp)s!" % vars()    
#     return ilist


def tableout(inp):
    #Read file to a list
    ifl=open(inp, 'r')
    idata=ifl.read()
    ifl.close()
    table=[]
    for row in idata.split('\n'):
        data=[it.strip() for it in row.split()]
        data=[numerize(it) for it in data]
        if len(data)!=0:
            table.append(data)

    return table



def collect(info):
    
    data=NestedDict()
    
    ilist=info['File'].values
    
    for index,row in info.iterrows():
        
        ifl=row['File']
        ipt, inm, itp = filename(ifl)
        print ifl
        if 'Reader' in info.columns.values:
            reader=row['Reader']
            
            if reader=='Tecan':
                
                if itp == 'xlsx':
                    sheet = readxls_s(ifl)
                    data[ifl] = collect_Tecan(sheet)
                elif itp in ['asc','txt']:
                    sheet = readtext(ifl)
                    data[ifl] = collect_Tecan(sheet)
                else:
                    raise Exception('Unknown file format: {}!'.format(itp))
                    
            elif reader=='Biotek':
                sheet=readtext(ifl)
                data[ifl]=collect_Biotek(sheet)

            else:
                raise Exception('Unknown reader: {}!'.format(reader))
                
        else:
            if itp in 'xlsx':
                sheet=readxls(ifl)
                data[ifl]=collect_Tecan(sheet)
            elif itp=='asc':
                sheet=readtext(ifl)
                data[ifl]=collect_Tecan(sheet)
            elif itp=='txt':
                sheet=readtext(ifl)
                data[ifl]=collect_Biotek(sheet)
            else:
                raise Exception('Unknown file format: {}!'.format(itp))

        data[ifl]['File']=ifl

    return data



def collect_Tecan(sheet):
    sheetdata=NestedDict()
    datarange=sheet[0].index('Well positions')

    stdrow=max([len(r) for r in sheet])
    
    sheet=[r for r in sheet if len(r)==stdrow]

    nrows=len(sheet)
    
    datarange=sheet[0].index('Well positions')
    nm_labels=[lab for lab in sheet[0] if lab not in ['Layout','Well positions','','Replicate Info']]

    
    if len(nm_labels)>1:
        starts=[sheet[0].index(lab) for lab in nm_labels]
    else:
        starts=[sheet[0].index(nm_labels[0])]
        if nm_labels[0]=='Raw data':
            nm_labels[0]='600'

    waves=[numerize(str(wlen).replace('nm','')) for wlen in nm_labels]
    #print 'Identified wavelengths: {}'.format(waves)
    #print datarange

    length=(datarange)/len(waves)

    #Extract time

    time_row=sheet[1][:length]
    if list(set([s for s in time_row if isinstance(s,float)])):
        time_t=time_row
    else:
        time_t=[int(str(t).replace('s','')) for t in time_row]

    #length=time_t.index(max(time_t))+1
    temp=[float(t.split()[0]) for t in sheet[2][:length]]

    timemax_min=int(round_to(float(time_t[-1])/60,5))
    timemax_h,timemax_remmin=divmod(timemax_min,60)

    #Time in seconds
    time=np.linspace(0,timemax_min*60,length,dtype=np.dtype(int))

    timestep=round_to(float(time_t[-1])/(length-1),1)

    alllabels=[r[datarange] for r in sheet if r[datarange] not in ['Well positions']]
    
    #Find first named well
    dstart=map(bool, alllabels).index(True)
    
    cleanlabels=alllabels[dstart:]
    cleandata=sheet[dstart+1:]
    
    #print "Non zero start {}".format(dstart)
    #print alllabels[dstart]
    #print sheet[dstart+1]
    #print alllabels
    
    labels=[l for l in cleanlabels if l!='']#Sheet
    #print labels

    plsize=len(labels)
    if plsize not in [12,48,96,384]:
        buffer=True
    else:
        buffer=False

    sheetdata['Labels']=labels
    sheetdata['Spectra']=[str(int(w))+'nm' for w in waves]
    sheetdata['Time']=time
    sheetdata['Temp']=temp
    sheetdata['Time_max']=timemax_min
    sheetdata['Time_step']=timestep
    sheetdata['Wells']=len(labels)
    sheetdata['Used wells']=plsize
    sheetdata['Buffered']=str(buffer)
    sheetdata['Figures']=[str(w)+'nm' for w in waves]
    #sheetdata['File']=inm
    print "Wavelengths: {}".format(waves)
    print "Run time {}, step {}min in {} wells\n".format(str(datetime.timedelta(minutes=timemax_min)),timestep/60, len(labels))
    
    for wave in waves:
        swave=str(int(wave))
        scol=(length)*(waves.index(wave))
        ecol=(length)*(waves.index(wave)+1)
        #print scol,ecol
        sheetvalues=[]
        for lab,well in enumerate(cleanlabels):
            if well!='':
                data_row=cleandata[lab][scol:ecol]
                sheetvalues.append(data_row)
                
        sheetdata[swave+'nm']=pd.DataFrame(sheetvalues,columns=time,index=labels)

    return sheetdata

def collect_Biotek(sheet):
    sheetdata=NestedDict()
    allwells=getallwells()

    rownames=[row[0] if len(row)!=0 else '' for row in sheet]

    time_ids=[ rin for rin,r in enumerate(rownames) if 'Time' in str(r) ]
    
    OD_ids=[ rin for rin,r in enumerate(rownames) if 'T' in str(r) and 'OD:' in str(r) ]

    time_row=sheet[time_ids[0]]
    temp_row=sheet[OD_ids[0]]

    nm_labels=[ str(sheet[rin][0].split()[1]) for rin in OD_ids]
    nm_labelsm={ rin:str(sheet[rin][0].split()[1]) for rin in OD_ids}

    waves_nmm={ rin: numerize(re.sub('OD:|GFP:','',str(sheet[rin][0].split()[1])) ) for rin in OD_ids}

    waves=[ numerize(re.sub('OD:|GFP:','',wlen)) for wlen in nm_labels if wlen!='']
    
    
    waves_nm=[str(int(w))+'nm' for w in waves]

    time_t=[time_to_sec(tval) for tval in time_row if tval not in ['Time','','0:00:00']]

    length=len(time_t)

    timestep=round_to(float(time_t[-1]-time_t[0])/(length-1),1)

    timemax_min=int((length-1)*timestep/60)
    
    timemax_h,timemax_remmin=divmod(timemax_min,60)

    time=np.linspace(0,timemax_min*60,length,dtype=np.dtype(int))

    temp=[float(t) for t in temp_row[1:] if t not in ['']]

    #print time_t[-1],timemax_h,time[-1],timestep

    alllabels=[rn if rn in allwells else '' for rn in rownames]
    labels=list(set([l for l in alllabels if l!='']))

    plsize=len(labels)
    if plsize not in [12,48,96,384]:
        buffer=True
    else:
        buffer=False

    sheetdata['Labels']=labels
    sheetdata['Spectra']=waves_nm
    sheetdata['Time']=time
    sheetdata['Temp']=temp
    sheetdata['Time_max']=timemax_min
    sheetdata['Time_step']=timestep
    sheetdata['Wells']=len(labels)
    sheetdata['Used wells']=plsize
    sheetdata['Buffered']=str(buffer)
    sheetdata['Figures']=waves_nm
    print "Wavelengths: {}".format(waves)
    print "Run time {}, step {}min in {} wells".format(str(datetime.timedelta(minutes=timemax_min)),timestep/60, len(labels))
    
    sheetvalues=[]
    wells=[]
    
    for rid,row in enumerate(sheet):
        if rid>min(OD_ids) and len(row)>0:
            well=str(row[0])
            OD_sel=max([ODid for ODid in OD_ids if rid>ODid])
            
            if well in allwells:
                data_row=[val for val in row[1:] if val not in ['']]
                sheetvalues.append(data_row)
                wells.append(well)

        if rid==len(sheet)-1 or rid in [ODid for ODid in OD_ids[1:]]:
            #print 'Collecting table with {} rows at row {}'.format(len(sheetvalues),rid)

            swave=str(waves_nmm[OD_sel])+'nm'
            sheetDF=pd.DataFrame(sheetvalues,columns=time,index=wells)

            sheetvalues=[]
            wells=[]

            sheetdata[swave]=sheetDF
    
    print "\n"
    
    return sheetdata



def loggrowth(x,y,y0=np.power(2.0,-4),thres=-5):
    #thres=-4
    #y=setbar(y,np.power(2,-5))
    #maxy=max(y)
    #miny=min(y)
    #yl=np.log2(y)
    #print l,y
    miny=min(y)
    maxy=max(y)
    gscale=maxy-miny

    #print miny, maxy
    x2,y2=cut(x, y, thres+(maxy-thres)*0.1, thres+(maxy-thres)*0.6) #0.6
    if len(y2)>0 and gscale>0.5 and len(set(y))>2:
        try:
            popt, pcov = curve_fit(lin, x2, y2)
            a=popt[0]
            c=popt[1]
            t0=(np.log2(y0)-c)/(a)
            tmax=(maxy-c)/(a)
        except TypeError:
            #print 'Curve_fit encountered an error!'
            a,c,t0,tmax=[np.nan,np.nan,np.nan,np.nan]
    else:
        a,c,t0,tmax=[np.nan,np.nan,np.nan,np.nan]

        
        # for par,nm in IT.izip([a,c,t0],['a','c','t0']):
    #     if allfit[tp]['Log'][nm]:
    #         allfit[tp]['Log'][nm]=allfit[tp]['Log'][nm]+[par]
    #     else:
    #         allfit[tp]['Log'][nm]=[par]
    
    return a,c,t0,tmax


def absgrowth(xvars,y,margin=0.01):
           
    maxg=max(y)
    
    scaleg=maxg-min(y)
    timec,growc=cut(xvars, y, 0, maxg,equalize=True)

    if scaleg>0.1 and len(timec)>10 and len(growc)>10 and len(np.unique(y))>2:
        #print plate,well
        #print len(timec), len(growc)
        #print min(growc),max(growc)
        #popt, pcov = curve_fit(growth, timec, growc,bounds=(0,np.inf),p0=[0.5, 5, 0.1],max_nfev=5000)
        #A,lam,u=popt
        try:
            popt, pcov = curve_fit(growth,xvars, y,bounds=(0,np.inf),p0=[0.5, 5, 0.1],max_nfev=5000)#,p0=[0.1,10,1]maxfev=5000
            A,lam,u=popt
            #print popt,tmaxf
        except (RuntimeError, ValueError, RuntimeWarning, UnboundLocalError) as e:
            #print e
            #print 'Curve_fit encountered an error in well {}!'.format(well)
            A,lam,u=[0,np.nan,0]

        if A>0 and lam<np.inf and u>0:
            yreducedf = growth(xvars,*popt) - max(growth(xvars,*popt))*(1-margin) #maxg#
            freducedf = ip.UnivariateSpline(xvars, yreducedf, s=0)
            if len(freducedf.roots())>0:
                tmaxf=freducedf.roots()[0]
            else:
                tmaxf=np.nan
        else:
            tmaxf=np.nan

        yreduced = growc - maxg*(1-margin)
        try:
            freduced = ip.UnivariateSpline(timec, yreduced, s=0)
            if len(freduced.roots()>0):
                tmax=freduced.roots()[0]
            else:
                tmax=np.nan
        except TypeError:
            #print 'Integration encountered an error in well {}!'.format(well)
            tmax=np.nan
    else:
        A,lam,u=[0,np.nan,0]
        tmaxf=np.nan
        tmax=np.nan

    #data[plate]['Summary']['GrowthFit'][well]=[]
    return A,lam,u,tmax


def analyse(data,full):

    window_h=2
    thres=np.power(2.0,-5)
    
    for plate in sorted(data.keys()):
        time=data[plate]['Time']
        time_h=time/3600.0
        
        msize=len(time)//10
        
        maxt=max(time_h)
        dt=time[1]-time[0]
        wells=data[plate]['Labels']
        
        #Estimate window for background subtraction
        if maxt>2:
            window=window_h*3600//dt
        else:
            window=2
            
        #print window
        
        #Needs to be checked
        waves=data[plate]['Spectra']
        #print waves
        #print wells
        
        for wave in waves:
            wv=wave.replace('nm','')
            print 'Analysing data in {}: {}'.format(plate,wave)
            #print time
            
            rawdata=data[plate][wave]#[well]

            #Change windows size
            nobagall=rawdata.apply( func=lambda row: pd.Series( setbar(row-np.mean(row[:window]) ,0.0),index=time),axis=1 )
            
            wfilt=nobagall.apply( func=lambda row: pd.Series( Wiener(row,msize), index=time) ,axis=1)
            
            
            logfilt=np.log2( wfilt )
            #logfilt=logfilt.apply(func=lambda row: pd.Series(np.where(np.isinf(row), None, row),index=time), axis=1)
            
            
            dtts=wfilt.apply(func=lambda row: pd.Series(ip.UnivariateSpline(time_h, row, s=0).derivative(1)(time_h),index=time),axis=1 )
            
            
#             logfiltdt=logfilt.apply( func=lambda row: pd.Series(ip.UnivariateSpline(time_h, row.replace(-np.inf, -5) ,s=0).derivative(1)(time_h) ,index=time),axis=1)
            
            #dtfilt=dtts.apply( lambda row: Wiener(row,msize//2),axis=1)
            
            
            data[plate][wave+'_b']=nobagall
            data[plate][wave+'_f']=wfilt
            data[plate][wave+'_log']=logfilt
            data[plate][wave+'_dt']=dtts
            #data[plate][wave+'_logdt']=logfiltdt
            
            if full:
                data[plate]['Figures']=data[plate]['Figures']+[wave+'_b',wave+'_f',wave+'_log',wave+'_dt']
            else:
                data[plate]['Figures']=data[plate]['Figures']+[wave+'_f',wave+'_log']
            
        summary=pd.DataFrame([],index=wells)
            
        for fg in data[plate]['Figures']:
            summaries=[summary]
            fgdata=data[plate][fg]
            
           
            if full:
                maxs=pd.DataFrame({ '{}_Max'.format(fg) : fgdata.apply(max,axis=1) })
                mins=pd.DataFrame({ '{}_Min'.format(fg) : fgdata.apply(min,axis=1) })
                summaries.extend([maxs,mins])
                
            if re.findall('_f$',fg):
                #print fgdata
                ints=fgdata.apply( lambda x: pd.Series({ '{}_AUC'.format(fg):ip.UnivariateSpline(time_h,x,s=0).integral(0, maxt)}),axis=1)
                
                logints=np.log2(ints.copy(deep=True)).rename(columns={'{}_AUC'.format(fg):'{}_logAUC'.format(fg)})
                
                summaries.extend([ints,logints])
                
                
                if full:
                    absgfit=fgdata.apply(lambda x: pd.Series( absgrowth(time_h,x), index=['{}_AG_A'.format(fg), '{}_AG_lamda'.format(fg), '{}_AG_u'.format(fg), '{}_AG_tmax'.format(fg)]),axis=1)
                    summaries.append(absgfit)
                
                
            if re.findall('_log$',fg):

                loggfit=fgdata.apply(lambda x: pd.Series( loggrowth(time_h,x), index=['{}_LG_a'.format(fg), '{}_LG_c'.format(fg), '{}_LG_t0'.format(fg), '{}_LG_tmax'.format(fg)]),axis=1)
                
                if full:
                    summaries.append(loggfit)
                else:
                    summaries.append(loggfit[['{}_LG_a'.format(fg),'{}_LG_c'.format(fg)]])
            
            
            summary=pd.concat(summaries, axis=1)

        data[plate]['Summary']=summary
        
    return data



def makesheets(data,descriptors,info):

    header_temp=['File','Well','Data']
    
    header=header_temp
    
    if len(info.keys())>0:
        hasinfo=[pl for pl in info.keys() if len(info[pl].keys())>0]
        if len(hasinfo)>0:
            infokeys=sorted(info[hasinfo[0]].keys())
            header+=infokeys
    
    if len(descriptors.keys())>0:
        hasdesc=[pl for pl in descriptors.keys() if len(descriptors[pl].keys())>0]
        if len(hasdesc)>0:
            desckeys=sorted(descriptors[hasdesc[0]].keys())
            header+=desckeys
            
    #allsumhead = ['590nm_f_Max', '590nm_log_Max','590nm_f_AUC', '590nm_f_logAUC'] # +\
                #['A', 'lamda', 'u', 'tmax','tmaxf'] + \
                 #['a_log', 'c_log', 't0_log', 'tmax_log'] + \

    allsummary=[]
    alldatats=[]

    for file in data.keys():
        wells=data[file]['Labels']
        
        summary = data[file]['Summary'].copy(deep=True)
        summary['File']=file
        summary['Data']='Summary'
        
        allsummary.append(summary)
        
        output = data[file]['Figures']
        for fig in output:
            datats=data[file][fig].copy(deep=True)     
            datats['File']=file
            datats['Data']=fig
            alldatats.append(datats)   
            
            
    allsummaryDF=pd.concat(allsummary)
    alldatatsDF=pd.concat(alldatats)
    
    allsummaryDF['Well']=allsummaryDF.index
    alldatatsDF['Well']=alldatatsDF.index
    
    allinfo=pd.merge(info,descriptors,on='Pattern')
    
    mainhead=['File','Pattern','Data','Well']
    
    allsummaryH=[ col for col in allsummaryDF.columns.values if col not in mainhead]
    alldatatsH=[col for col in alldatatsDF.columns.values if col not in mainhead]
    allinfoH=[col for col in allinfo.columns.values if col not in mainhead]
    
    
    allsummaryDF=pd.merge(allsummaryDF,allinfo,on=['File','Well'])
    alldatatsDF=pd.merge(alldatatsDF,allinfo,on=['File','Well'])
    
    #Reorder DF
    allsummaryDF=allsummaryDF[mainhead+allinfoH+allsummaryH]
    alldatatsDF=alldatatsDF[mainhead+allinfoH+alldatatsH]
    allinfo=allinfo[['File','Pattern','Well']+allinfoH]

    sheets={'Summary':allsummaryDF,'Timeseries':alldatatsDF,'Allinfo':allinfo}

    return sheets


def writesheets(sheets,odir):
    for sheetname in sheets.keys():
        #print fig
        oname='{}/{}.csv'.format(odir,sheetname)        
        sheet=sheets[sheetname]
        sheet=sheet.replace([np.inf, -np.inf], np.nan)
        sheet.to_csv(oname,index=False,float_format='%.4f')
        
def plot_comparison(data,dirn,figs):

    for plate in sorted(data.keys()):
        
        ipt, inm, itp = filename(plate)
        
        labels=data[plate]['Labels']
        nwells=data[plate]['Wells']
        plsize=data[plate]['Used wells']
        buffered=data[plate]['Buffered']
        fgdata=data[plate]
        time=data[plate]['Time']
        time_h=time/3600.0
        
        fgsummary=fgdata['Summary']
        
        figures_temp=data[plate]['Figures']
        
        if figs=='all':
            figures=figures_temp
            
        elif isinstance(figs, list):
            figs_ch=[f for f in figs if f in figures_temp]
            if len(figs_ch)>0:
                figures=figs_ch
            else:
                print 'Figures {} not found'.format(figs)
        elif isinstance(figs, str) or isinstance(figs, unicode):
            if figs in figures_temp:
                figures=[figs]
            else:
                figures=[fg for fg in figure_temp if figs in fg]
                if not figures:
                    print 'Figures {} not found'.format(figs)

        #print figures
        print 'File: {}'.format(plate)
        for fg in figures:
            if re.findall('_log$',fg) and '{}_LG_a'.format(fg) in fgsummary.columns.values:
                
                gfitc=fgsummary[['{}_LG_a'.format(fg),'{}_LG_c'.format(fg)]]
            else:
                gfitc=pd.DataFrame()

            print "\tPlotting {}...".format(fg)
            
            

            #Need to fix metabolites
            plot=plot_2D(inm,fg,fgdata[fg],time_h,labels,gfitc,plsize)
            plot.savefig('{}/{}.pdf'.format(dirn,inm+'_'+fg))
            plot.close()



def plot_2D(plate,fg,datac,time,labels,gfitc,plsize):
    
    if plsize==96:
        minrow=1
        maxrow=8
        mincol=1
        maxcol=12
        cols=12
        rows=8
    if plsize==60:
        minrow=1
        maxrow=6
        mincol=1
        maxcol=10
        cols=10
        rows=6
    elif plsize==12:
        minrow=1
        maxrow=3
        mincol=1
        maxcol=4
        cols=4
        rows=3
    elif plsize==384:
        minrow=1
        maxrow=16
        mincol=1
        maxcol=24
        cols=24
        rows=16
    elif plsize==240:
        minrow=1
        maxrow=14
        maxcol=22
        cols=20
        rows=12
        
    title="{} {}".format(plate,fg)
    
    fig,axes=plt.subplots(nrows=rows, ncols=cols, sharex=True, sharey=True,figsize=(11.69,8.27), dpi=100)
    fig.suptitle(title)
    plt.subplots_adjust(left=0.1, bottom=0.1, right=0.95, top=0.9, wspace=0.05, hspace=0.05)
    
    
    xmax=max(time)
    
    
    xlabel='Time, h'
    ylabel=fg
    decimals=1
    
    #print '1'
    
    totalmax=round( max( datac.max() ),1)
    totalmin=round( min( datac.min() ),1)
    
    
    #print '2'
    
    if totalmax<0:        
        totalmaxc=0
    #totalmin=min(min(data))
    
    if re.findall('_log$',fg):
        totalmax=0
        totalmin=-6
        decimals=1
        
    #print '3'
        
    if re.findall('_logdt',fg):
        decimals=1

    if re.findall('_dt$',fg):
        totalmin=0
        decimals=2

    #print '4'
    
    if plsize>12:
        ticks=3
    else:
        ticks=5
        
    

    ymin=totalmin
    fig.text(0.5, 0.04, xlabel, ha='center')
    fig.text(0.04, 0.5, ylabel, va='center', rotation='vertical')

    for v,l in IT.izip(range(plsize),labels):
        #print l
        if plsize==240:
            #Numbering inside plot arrangement differs from real numbers by buffer size
            row=string.uppercase.index(l[0])+1-2
            col=int(l.replace(l[0],''))-2
        elif plsize==60:
            row=string.uppercase.index(l[0])+1-1
            col=int(l.replace(l[0],''))-1
        else:
            row=string.uppercase.index(l[0])+1
            col=int(l.replace(l[0],''))

        if divmod(float(v)/12,1)[1]==0:
            sh_y=l

        v=v+1
    
        #x=time/3600.0
        
        yc=datac.loc[l,:]

        if re.findall('_log$',fg):
            ca,cc=gfitc.loc[l,:]

        plt.sca(axes[row-1,col-1])
        ax=axes[row-1,col-1]
        #print col,cols
        #print row,rows
        if col==cols:
            ax.yaxis.set_label_position("right")
            plt.ylabel(l[0],rotation='horizontal')
        if row==minrow:
            plt.title(str(int(l.replace(l[0],''))))

        if col>1 and row<rows-1:
            plt.setp(ax.get_yticklabels(), visible=False)
        
        #print xmax, ticks
        #print np.linspace(0, xmax, 3),['']+list(np.linspace(0, xmax, 3).astype(int)[1:])
        plt.xticks(np.linspace(0, xmax, 3),['']+list(np.linspace(0, xmax, 3).astype(int)[1:]), rotation='vertical')    
        plt.yticks(np.linspace(ymin, totalmax, ticks),['']+list(myround(np.linspace(totalmin, totalmax, ticks),decimals)[1:]))
        plt.ylim([totalmin,totalmax])
        #plt.axvline(x=3, ymin=0, ymax=1,c='green', hold=None)  u 

        #label=greek_check(metabolites[l]['Name'],12)
        plt.text(0.05, 0.9, l, fontsize=7,verticalalignment='top',transform=ax.transAxes)

        #print '{}: {} {} {}'.format(fg,len(x),len(yc),len(ye))
        plt.plot(time,yc,'r-')
        
        if re.findall('_log$',fg):
            if ca>0:
                yfitc=time*ca+cc
                plt.plot(time,yfitc,'r-',alpha=0.5)

    return plt

#----------------------------------

if __name__ == "__main__":
    sys.exit(main())

