#Coded by Ananta Acharya
#contact ananta@uga.edu if any questions
#This can be used as is or you can modify, The results are not guaranted to be correct, use with caution

import os, glob
from time import time
from scipy.misc import comb
from collections import defaultdict
def main():
    print "1. Make haplotype and count files from matches (For each SNP)"
    print "2. Define criteria for filtering SNPs (from outputs of 1 above)"
    #print "3. Change the format from A/C to binary for each SNP"
    print "3. Make haplotype and count files from matches (Haplotype format, like stacks)"
    print "4. Make haplotype and count files from stacks haplotype files and matches)"
    print "5. Make single snp from hap)"
    
    opt=raw_input()
    if opt=="1":
        
        infile=raw_input('''Where is the file:only directory \n
                     All matches file should be there''').strip().strip('"')
        sep=raw_input("seperator for genotypes (columns) ie , or tab").strip()
        if sep=="tab":
                sep="\t"
        allelesep=raw_input("seperator for alleles ie | , : ,/").strip()
        stackdepth(infile, sep, allelesep)
    elif opt=="2":
        
 
        infilesnp=raw_input("SNPs file processed with No. 1 option above").strip().strip('"')
        infilecount=raw_input("counts file processed with No. 1 option above").strip().strip('"')
        min_homo=raw_input("minimum depth to call homozygous")
        min_hetero=raw_input("minimum depth to call heterozygote (combined)")
        heteroratio=raw_input("minimum ratio to call heterozygous")
        flag=raw_input("what to use for homozygous corection? if binomial test b, if minratio r")
        outf=raw_input("output").strip().strip('"')
        minind=int(raw_input("min number of individuals").strip())
        sep=raw_input("seperator for genotypes (columns) ie , or tab").strip()
        if sep=="tab":
                sep="\t"
        allelesep=raw_input("seperator for alleles ie | , : ,/").strip()
        if flag=="r":
            homoratio=raw_input("maximum ratio to correct  heterozygous as homozygous ")
            rehapmapratio(infilesnp, infilecount,outf,min_homo, min_hetero, heteroratio, homoratio, minind, sep, allelesep)
        if flag=="b":
            binary=raw_input("binomial test coefficient (if diploid=0.005, if autotetraploid 0.0025 ")
        
            rehapmapbinary(infilesnp, infilecount,outf,min_homo, min_hetero, heteroratio, binary,minind, sep, allelesep)

    elif opt=="3":
        
        infile=raw_input('''Where is the file:only directory \n
                     All matches file should be there''').strip().strip('"')
        sep=raw_input("seperator for genotypes (columns) ie , or tab").strip()
        if sep=="tab":
                sep="\t"
        allelesep=raw_input("seperator for alleles ie | , : ,/").strip()
        stackdepthHaplo(infile, sep, allelesep)
    elif opt=="4":
        infile=raw_input('''Where is the Haplotype file from stacks? batch_x.haplotype.tsv''').strip().strip('"')
        dirr=raw_input('''Where is the file:only directory \n
                     All matches file should be there''').strip().strip('"')
        sep=raw_input("seperator for genotypes (columns) ie , or tab").strip()
        if sep=="tab":
                sep="\t"
        allelesep=raw_input("seperator for alleles ie | , : ,/").strip()
        stackdepthHaplofromstacksHaplo(infile, dirr,sep, allelesep)
    elif opt=="5":
        hapfile=raw_input("hapfile").strip().strip('"')
        useSNPpos=raw_input("Use snp positions?y/n")
        if "y" in useSNPpos:
            
            snpfile=raw_input("catalog snp file").strip().strip('"')
            snpposFlag=True
        else:
            snpfile=""
            snpposFlag=False
        outfile=raw_input("outfile").strip().strip('"')
        sep=raw_input("seperator for genotypes (columns) ie , or tab").strip()
        if sep=="tab":
                sep="\t"
        allelesep=raw_input("seperator for alleles ie | , : ,/").strip()
        genstart=int(raw_input("genotypes start from position default: 5: ").strip())
        hapToSingleSNP(hapfile, snpfile,outfile, sep, allelesep, genstart,snpposFlag)


    

def rehapmapbinary(hapmap,hapcount, outf, mindepth,mindepthhetero, minratio,r, sep, allelesep):
    mindepth=int(mindepth)-1
    mindepthhetero=int(mindepthhetero)-1
    minratio=float(minratio)
    r=float(r)
    Allpairdict=dict()
    hapfile=open(hapmap,"r")

    formatt="2"
    rehap=open(outf,"a")
    for line in hapfile:
        linelist=line.strip().split(sep)
        Allpairdict[linelist[0]]=linelist


    with open(hapcount) as hapcountf:
        line1=next(hapcountf)
        rehap.write(line1)
    
        
        writebuffer=""

        for line in hapcountf:

            linelist=line.strip().split(sep)
            ID=linelist[0]
            writebuffer+=ID
            for x in range(1,len(linelist)-1):

                counts=linelist[x].split(allelesep)

                if len(counts)==1:
                    if counts[0] in ["NA"]:
                        towrite=",NA"
                        writebuffer+=towrite
                        
                    else:
                        
                        alleledepth=int(counts[0])
                        
                        if alleledepth>mindepth:

                            towrite=Allpairdict[ID][x]
                            writebuffer+=sep+towrite

                        else:

                            towrite=",NA"
                            writebuffer+=towrite


                       
                    
                else:
                    allele1=int(counts[0])
                    allele2=int(counts[1])
                    if allele1+ allele2>mindepthhetero:
                         
                        snps=Allpairdict[ID][x]

                        
                        if minratio<=(allele1/float(allele1+allele2))<=(1-minratio):
                            writebuffer+=sep+snps
                        else:
                            snpslist=snps.split(allelesep)
                            
                            minorallele=min(allele1, allele2)
                            if minorallele==allele1:
                                minorsnp=snps[0]
                                majorsnp=snps[1]
                            else:
                                minorsnp=snps[1]
                                majorsnp=snps[0]
                            
                            p1=binomial_test((allele1+allele2), minorallele,r)

                            if p1>=0.1:
                                
                                writebuffer+=sep+majorsnp
                            elif p1<0.001:
                                
                                writebuffer+=sep+snps
                            else:
                                    
                                
                                writebuffer+=",NA"
                    else:
                        writebuffer+=",NA"
##                print(writebuffer)
##                raw_input()
                

            writebuffer+="\n"
        rehap.write(writebuffer)

def rehapmapratio(hapmap,hapcount, outf, mindepth,mindepthhetero, minratio, maxratio, minsample, sep, allelesep):
    
    mindepth=int(mindepth)-1
    mindepthhetero=int(mindepthhetero)-1
    minratio=float(minratio)
    maxratio=float(maxratio)
    Allpairdict=dict()
    hapfile=open(hapmap,"r")

    formatt="2"
    rehap=open(outf,"w")
    for line in hapfile:
        linelist=line.strip().split(sep)
        if "consesus" in linelist:
            continue
        #print str(linelist.count("NA"))
        elif linelist.count("NA")>len(linelist)-minsample-1:
            continue
        Allpairdict[linelist[0]]=linelist

    print "finished reading haplotypes"
    with open(hapcount) as hapcountf:
        line1=next(hapcountf)
        rehap.write(line1)
    
        
        writebuffer=""
        count=0
        for line in hapcountf:
            count+=1
            linelist=line.strip().split(sep)
            ID=linelist[0]
            if ID not in Allpairdict.keys():
                continue
            writebuffer+=ID
            
            for x in range(1,len(linelist)-1):

                counts=linelist[x].split(allelesep)

                if len(counts)==1:
                    if counts[0] in ["NA"]:
                        towrite=sep+"NA"
                        writebuffer+=towrite
                        
                    else:
                        
                        alleledepth=int(counts[0])
                        
                        if alleledepth>mindepth:

                            towrite=Allpairdict[ID][x]
                            writebuffer+=sep+towrite

                        else:

                            towrite=sep+"NA"
                            writebuffer+=towrite


                       
                    
                else:
                    allele1=int(counts[0])
                    allele2=int(counts[1])
                    if allele1+ allele2>mindepthhetero:
                         
                        snps=Allpairdict[ID][x]

                        
                        if minratio<=(allele1/float(allele1+allele2))<=(1-minratio):
                            writebuffer+=sep+snps
                        else:
                            snpslist=snps.split(allelesep)
                            
                            minorallele=min(allele1, allele2)
                            if minorallele==allele1:
                                minorsnp=snps[0]
                                majorsnp=snps[1]
                            else:
                                minorsnp=snps[1]
                                majorsnp=snps[0]

                            if min(allele1, allele2)/float(allele1+allele2)<maxratio:
                                
                            
                                
                                writebuffer+=sep+majorsnp
                            
                            else:
                                    
                                
                                writebuffer+=sep+"NA"
                    else:
                        writebuffer+=sep+"NA"
##                print(writebuffer)
##                raw_input()
                

            writebuffer+="\n"
            if count%1000==0:
                rehap.write(writebuffer)
                
                writebuffer=""
                print "writing " +str(count/1000)+"000 catalogs"
                
        rehap.write(writebuffer)


    

def stackdepth(infile, sep, allelesep):


    dirr=infile 

    dictsnp=dict()

    start=time()

    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()
    maxcat=0

    os.chdir(dirr)
    print "Looking in" +dirr
    allMatchFiles=glob.glob("*.matches.tsv")
    if len(allMatchFiles)<1:
        print "Match files not found"
        return
    gens=[]

    for eachfilex in allMatchFiles:
        print "parsing, " +eachfilex
        gen=eachfilex[:-12]
        gens.append(gen)

        eachfile=open(eachfilex, "r")

        for line in eachfile:
            linelist=line.strip().split()
            catID=linelist[2]
            maxcat=max(maxcat,int(catID))
            haplo=linelist[5]
            depth=linelist[6]
            haplodict[gen, catID, haplo]=depth
            if not haplo.startswith("consensus"):
                for i in range(len(haplo)):
                    snpdict[gen,catID,i,haplo[i]]=int(snpdict.get((gen,catID,i,haplo[i]),0))+int(depth)
            gencatdict[gen, catID]=1
            if haplo not in haplodict2.setdefault(catID,[]):
                
                haplodict2.setdefault(catID,[]).append(haplo)
        eachfile.close()

    writebuffer=""
    gens.sort()
    outfilesnp=open(infile+"/SNPAllmatches.snp.csv","w")
    outfilecount=open(infile+"/SNPAllmatches.count.csv","w")
    
    outfilesnp.write("Catalog,")
    outfilecount.write("Catalog,")
    
    outfilesnp.write(sep.join(gens))
    outfilesnp.write("\n")
    outfilecount.write(sep.join(gens))
    outfilecount.write("\n")
    writebuffersnp=""
    writebuffercount=""

    
    count=0
    for cats in range(1,maxcat+1):
        count+=1
        
        haploo=haplodict2.get(str(cats),0)
 #       print haploo
        
        if haploo<>0 and ("consensus" not in haploo):
            haploo=haploo[0]

    
        
            for i in range(len(haploo)):
                linebuffersnp=str(cats)+"_"+str(i)+sep
                linebuffercount=str(cats)+"_"+str(i)+sep
                for gen in gens:
                    gensnp=list()
                    gencount=list()
                    deptheach=0
                    eachgensnpdict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                        
                    
                        for nc in ["AsepTsepCsepG"]:
    
                            depthh=snpdict.get((gen,str(cats),i,nc),0)
                            if depthh>0:
                                eachgensnpdict[nc]=depthh
                                
                                deptheach+=depthh
                            if depth<=0:
                                eachgensnpdict[nc]=0
                        for key, val in eachgensnpdict.items():
                            gensnp.append(key)
                            gencount.append(str(val))
                            
                        gensnpwrite=allelesep.join(gensnp)
                        gencountwrite=allelesep.join(gencount)
                        #print str(i),gensnpwrite, gencountwrite
                        linebuffersnp+=gensnpwrite+sep
                        linebuffercount+=gencountwrite+sep
                    else:
                        linebuffersnp+="NA,"
                        linebuffercount+="NA,"
                        
                writebuffersnp+=linebuffersnp+"\n"
                writebuffercount+=linebuffercount+"\n"
                  

            if count%1000==0:
                noww=time()
                timee=str((noww-start)/60)
                
                print "%s rows in %s minutes" %(count,timee)               
                    
                outfilesnp.write(writebuffersnp)
                outfilecount.write(writebuffercount)
                writebuffersnp=""
                writebuffercount=""

        
    noww=time()
    timee=str((noww-start)/60)
    
    print "%s rows in %s minutes" %(count,timee)               
        
    outfilesnp.write(writebuffersnp)
    outfilecount.write(writebuffercount)
    outfilesnp.close()
    outfilecount.close()
    
def stackdepthSNPfromHaplo(haplofile, snpfile, dirr, sep, allelesep, genstart):
    from collections import defaultdict 
    start=time()
    with open(haplofile,"r") as hapin:
        line1list=hapin.readline().strip().split("\t")
    
        
    inds=line1list[genstart-1::]
    
    with open(haplofile,"r") as hapin:
        line1list=hapin.readline().strip().split("\t")
    os.chdir(dirr)
    print "Looking in" +dirr
    allMatchFiles=glob.glob("*.matches.tsv")
    allMatchfilesnoExt=[y[:-12] for y in allMatchFiles]
    notfound=False
    for ind in inds:
        if ind not in allMatchfilesnoExt:
            print ind+ " not found"
            notfound=True
        
    if notfound:        
        return
    hapcatlist=[]
    with open(haplofile,"r") as hapin:
        for line in hapin:
            try:
                hapcatlist.append(int(line.strip().split("\t")[0]))
            except:
                pass
    catsnpdict=defaultdict(list)
    with open(snpfile,"r") as snpin:
        for line in snpin:
            try:
                linelist=line.strip().split("\t")
                cat=int(linelist[2])
                pos=linelist[2]
                catsnpdict[cat]=catsnpdict[cat].append(pos)
            except:
                pass
    print str(len(hapcatlist))+" catalogs found"
    gens=[]
    dictsnp=dict()
    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()
    maxcat=0

    

    dictsnp=dict()

    start=time()

    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()


    gens=[]

    for eachfilex in allMatchFiles:
        print "parsing, " +eachfilex
        gen=eachfilex[:-12]
        gens.append(gen)

        eachfile=open(eachfilex, "r")

        for line in eachfile:
            linelist=line.strip().split()
            catID=linelist[2]
            maxcat=max(maxcat,int(catID))
            haplo=linelist[5]
            depth=linelist[6]
            haplodict[gen, catID, haplo]=depth
            if not haplo.startswith("consensus"):
                for i in range(len(haplo)):
                    snpdict[gen,catID,i,haplo[i]]=int(snpdict.get((gen,catID,i,haplo[i]),0))+int(depth)
            gencatdict[gen, catID]=1
            if haplo not in haplodict2.setdefault(catID,[]):
                
                haplodict2.setdefault(catID,[]).append(haplo)
        eachfile.close()

    writebuffer=""
    gens.sort()
    outfilesnp=open(infile+"/SNPAllmatches.snp.csv", "w")
    outfilecount=open(infile+"/SNPAllmatches.count.csv","w")
    
    outfilesnp.write("Catalog,")
    outfilecount.write("Catalog,")
    
    outfilesnp.write(sep.join(gens))
    outfilesnp.write("\n")
    outfilecount.write(sep.join(gens))
    outfilecount.write("\n")
    writebuffersnp=""
    writebuffercount=""

    
    count=0
    for cats in hapcatlist:
        count+=1
        
        haploo=haplodict2.get(str(cats),0)
 #       print haploo
        
        if haploo<>0 and ("consensus" not in haploo):
            haploo=haploo[0]

    
        
            for i in range(len(haploo)):
                catSNPpos=catsnpdict[cats][i]
                linebuffersnp=str(cats)+"_"+catSNPpos+sep
                linebuffercount=str(cats)+"_"+catSNPpos+sep
                for gen in gens:
                    gensnp=list()
                    gencount=list()
                    deptheach=0
                    eachgensnpdict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                        
                    
                        for nc in ["A","T","C","G"]:
    
                            depthh=snpdict.get((gen,str(cats),i,nc),0)
                            if depthh>0:
                                eachgensnpdict[nc]=depthh
                                
                                deptheach+=depthh
                            if depth<=0:
                                eachgensnpdict[nc]=0
                        for key, val in eachgensnpdict.items():
                            gensnp.append(key)
                            gencount.append(str(val))
                            
                        gensnpwrite=allelesep.join(gensnp)
                        gencountwrite=allelesep.join(gencount)
                        #print str(i),gensnpwrite, gencountwrite
                        linebuffersnp+=gensnpwrite+sep
                        linebuffercount+=gencountwrite+sep
                    else:
                        linebuffersnp+="NA"+sep
                        linebuffercount+="NA"+sep
                        
                writebuffersnp+=linebuffersnp+"\n"
                writebuffercount+=linebuffercount+"\n"
                  

            if count%1000==0:
                noww=time()
                timee=str((noww-start)/60)
                
                print "%s rows in %s minutes" %(count,timee)               
                    
                outfilesnp.write(writebuffersnp)
                outfilecount.write(writebuffercount)
                writebuffersnp=""
                writebuffercount=""

        
    noww=time()
    timee=str((noww-start)/60)
    
    print "%s rows in %s minutes" %(count,timee)               
        
    outfilesnp.write(writebuffersnp)
    outfilecount.write(writebuffercount)
    outfilesnp.close()
    outfilecount.close()
    
def hapToSingleSNP(infile, snpfile,outfile, sep, allelesep, genstart, snpposFlag):
    catsep="\t"
    genstart=genstart-1
    start=time()
    dictsnp=dict()
    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()
    gens=[]
    cats=set()
    hapfile=open(infile, "r")
    firstline=hapfile.readline().strip().split(sep)
    gens=firstline[genstart::]
    for j,line in enumerate(hapfile):
        linelist=line.strip().split(sep)
        catID=linelist[0]
        
        for i,hap in enumerate(linelist[genstart::]):
            if "consensus" in hap:
                break
            else:
                haplodict[gens[i], catID]=hap
               # print gens[i], catID, hap
               # raw_input()
                try:
                    cats.add(int(catID))
                except:
                    pass
##        if j>100:
##            break
        
    hapfile.close()
    
    print "finished Reading haplotypes"
    print str((time()-start)/60)+" minutes"
    writebuffer=""

    catsnpdict=defaultdict(list)
    with open(snpfile,"r") as snpin:
        for line in snpin:
        
            linelist=line.strip().split(catsep)
            cat=int(linelist[2])
            pos=linelist[3]
            catsnpdict[cat].append(pos)
            
    print "finished Reading snp positions"
    print str((time()-start)/60)+" minutes"
    gens.sort()
    outfilesnp=open(outfile,"w")
    outfilesnp.write("Catalog,")
   
    outfilesnp.write(sep.join(gens))
    outfilesnp.write("\n")
   
    writebuffersnp=""


    catslist=list(cats)
    catslist.sort()
    count=0
    for cat in cats:
        count+=1
        allgenhap=list()
        for gen in gens:
            
            genhaplo=haplodict.get((gen, str(cat)),"NA")
           # print genhaplo
            genhaplos=genhaplo.split(allelesep)
           # print genhaplos
            allgenhap.extend(genhaplos)        
        allgenhap=set(allgenhap)
        #print str(cat), allgenhap
       # raw_input()
        #print str(cat), allgenhap
        allgenhap.discard("NA")
        
        allgenhap.discard("-")
        try:
            nhap=len(allgenhap)
        except:
            nhap=0
       # print allgenhap, str(nhap), cat
        if nhap>1:
            allgenhap=list(allgenhap)
            snpPair=zip(*allgenhap)
            snpPair=[list(set(d)) for d in snpPair]
            catSNP=catsnpdict[cat]
            #print str(cat),catsnpdict[cat]
            #raw_input()
            
            for i,snp in enumerate(snpPair):
                try:
                    linebuffersnp=str(cat)+"_"+catSNP[i]+sep
                except:
                    #print str(cat),catsnpdict[cat]
                    linebuffersnp=str(cat)+"_"+str[i]+"unknown"+sep
                    #raw_input()
                for gen in gens:
                    genhap="NA"
                    genhapp=haplodict.get((gen, str(cat)),"NA")
                    #print genhapp
                   # raw_input()
                    if genhapp not in ["NA","-"]:
                
                        genhaploss=genhapp.split(allelesep)
                    
                        gen=allelesep.join(list(set(zip(*genhaploss)[i])))
                    
                                    
                
                        linebuffersnp+=gen+sep
                
                    else:
                        
                        linebuffersnp+="NA"+sep
               
                
                writebuffersnp+=linebuffersnp+"\n"
        
                  

        if count%1000==0:
            noww=time()
            timee=str((noww-start)/60)
            
            print "%s rows in %s minutes" %(count,timee)               
                
            outfilesnp.write(writebuffersnp)
            
            writebuffersnp=""
            

        
    noww=time()
    timee=str((noww-start)/60)
    
    print "%s rows in %s minutes" %(count,timee)               
        
    outfilesnp.write(writebuffersnp)
   
    outfilesnp.close()
    
    


def stackdepthHaplo(infile, sep, allelesep):


    dirr=infile 


    dictsnp=dict()

    start=time()

    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()
    maxcat=0

    os.chdir(dirr)
    print "Looking in" +dirr
    allMatchFiles=glob.glob("*.matches.tsv")
    if len(allMatchFiles)<1:
        print "Match files not found"
        return
    gens=[]
    #print allMatchFiles
    for eachfilex in allMatchFiles:
        print "parsing, " +eachfilex
        gen=eachfilex[:-12]
        gens.append(gen)
#        print gen
        eachfile=open(eachfilex, "r")

        for line in eachfile:
            linelist=line.strip().split()
            catID=linelist[2]
            maxcat=max(maxcat,int(catID))
            haplo=linelist[5]
            depth=linelist[6]
            haplodict[gen, catID, haplo]=depth

            gencatdict[gen, catID]=1
            if haplo not in haplodict2.setdefault(catID,[]):
                
                haplodict2.setdefault(catID,[]).append(haplo)
        eachfile.close()


    gens.sort()
    outfilehaplo=open(infile+"/SNPAllmatches.haplo.csv","w")
    outfilehaplocount=open(infile+"/SNPallmatches.haplocount.csv","w")
    outfilehaplodepth=open(infile+"/SNPAllmatches.readdepth.csv","w")

    outfilehaplo.write("Catalog,")
    outfilehaplocount.write("Catalog,")
    outfilehaplodepth.write("Catalog,")
    
    outfilehaplo.write(sep.join(gens))
    outfilehaplo.write("\n")
    outfilehaplocount.write(sep.join(gens))
    outfilehaplocount.write("\n")
    outfilehaplodepth.write(sep.join(gens))
    outfilehaplodepth.write("\n")
    writebufferhaplo=""
    writebufferhaplocount=""
    writebufferhaplodepth=""


    
    count=0
    for cats in range(1,maxcat+1):
        count+=1
        
        haploo=haplodict2.get(str(cats),0)

        

        if haploo<>0:
            if "consensus" in haploo:
                linebufferhaplo=str(cats)+sep
                linebufferhaplocount=str(cats)+sep
                linebufferhaplodepth=str(cats)+sep
                for gen in gens:
                    #genhaplo=list()
                    genhaplocount=list()
                    deptheach=0
                    #eachgenhaplodict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                      
                        depthh=haplodict.get((gen,str(cats),"consensus"),0)     
                            



                                
                        genhaplowrite="consensus"
                        genhaplocountwrite=str(depthh)
                        #print str(i),gensnpwrite, gencountwrite
                        linebufferhaplo+=genhaplowrite+sep
                        linebufferhaplocount+=genhaplocountwrite+sep
                        linebufferhaplodepth+=genhaplocountwrite+sep
                    else:
                        linebufferhaplo+="NA,"
                        linebufferhaplocount+="NA,"
                        linebufferhaplodepth+="NA,"
                writebufferhaplo+=linebufferhaplo+"\n"
                writebufferhaplocount+=linebufferhaplocount+"\n"
                writebufferhaplodepth+=linebufferhaplodepth+"\n"
            

    
        
            else:
                linebufferhaplo=str(cats)+sep
                linebufferhaplocount=str(cats)+sep
                linebufferhaplodepth=str(cats)+sep
                for gen in gens:
                    genhaplo=list()
                    genhaplocount=list()
                    deptheach=0
                    eachgenhaplodict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                        for haplotype in haploo:
                            depthh=haplodict.get((gen,str(cats),haplotype),0)     
                            
                            if depthh>0:
                                eachgenhaplodict[haplotype]=depthh

                        for key, val in eachgenhaplodict.items():
                            genhaplo.append(key)
                            genhaplocount.append(str(val))
                                
                        genhaplowrite=allelesep.join(genhaplo)
                        genhaplocountwrite=allelesep.join(genhaplocount)
                        genhaplodepthwrite=str(sum([int(y) for y in genhaplocount]))
                        #print str(i),gensnpwrite, gencountwrite
                        linebufferhaplo+=genhaplowrite+sep
                        linebufferhaplocount+=genhaplocountwrite+sep
                        linebufferhaplodepth+=genhaplodepthwrite+sep
                    else:
                        linebufferhaplo+="NA,"
                        linebufferhaplocount+="NA,"
                        linebufferhaplodepth+="NA,"
                        
                writebufferhaplo+=linebufferhaplo+"\n"
                writebufferhaplocount+=linebufferhaplocount+"\n"
                writebufferhaplodepth+=linebufferhaplodepth+"\n"
                  

            if count%1000==0:
                noww=time()
                timee=str((noww-start)/60)
                
                print "%s rows in %s minutes" %(count,timee)               
                    
                outfilehaplo.write(writebufferhaplo)
                outfilehaplocount.write(writebufferhaplocount)
                outfilehaplodepth.write(writebufferhaplodepth)
                writebufferhaplo=""
                writebufferhaplocount=""
                writebufferhaplodepth=""

        
    noww=time()
    timee=str((noww-start)/60)
    
    print "%s rows in %s minutes" %(count,timee)               
        
    outfilehaplo.write(writebufferhaplo)
    outfilehaplocount.write(writebufferhaplocount)
    outfilehaplodepth.write(writebufferhaplodepth)
    outfilehaplo.close()
    outfilehaplocount.close()
    outfilehaplodepth.close()
    


def stackdepthHaplofromstacksHaplo(haplofile, dirr, sep, allelesep):
    start=time()
    with open(haplofile,"r") as hapin:
        line1list=hapin.readline().strip().split(sep)
        
    inds=line1list[3::]
    
    
    os.chdir(dirr)
    print "Looking in" +dirr
    allMatchFiles=glob.glob("*.matches.tsv")
    allMatchfilesnoExt=[y[:-12] for y in allMatchFiles]
    notfound=False
    for ind in inds:
        if ind not in allMatchfilesnoExt:
            print ind+ " not found"
            notfound=True
        
    if notfound:        
        return
    hapcatlist=[]
    with open(haplofile,"r") as hapin:
        for line in hapin:
            try:
                hapcatlist.append(int(line.strip().split(sep)[0]))
            except:
                pass
    hapcatlist=list(set(hapcatlist))
    hapcatlist.sort()
    print str(len(hapcatlist))+" catalogs found"
    gens=[]
    dictsnp=dict()
    haplodict=dict()
    haplodict2=dict()
    gencatdict=dict()
    snpdict=dict()
    maxcat=0

    
    #print allMatchFiles
    for eachfilex in allMatchFiles:
        print "parsing, " +eachfilex
        gen=eachfilex[:-12]
        gens.append(gen)
#        print gen
        eachfile=open(eachfilex, "r")

        for line in eachfile:
            linelist=line.strip().split()
            catID=linelist[2]
            maxcat=max(maxcat,int(catID))
            haplo=linelist[5]
            depth=linelist[6]
            haplodict[gen, catID, haplo]=depth

            gencatdict[gen, catID]=1
            if haplo not in haplodict2.setdefault(catID,[]):
                
                haplodict2.setdefault(catID,[]).append(haplo)
        eachfile.close()


    gens.sort()
    outfilehaplo=open(haplofile+".haplo.txt","w")
    outfilehaplocount=open(haplofile+".haplocount.txt","w")
    outfilehaplodepth=open(haplofile+".readdepth.txt","w")

    outfilehaplo.write("Catalog"+sep)
    outfilehaplocount.write("Catalog"+sep)
    outfilehaplodepth.write("Catalog"+sep)
    
    outfilehaplo.write(sep.join(gens))
    outfilehaplo.write("\n")
    outfilehaplocount.write(sep.join(gens))
    outfilehaplocount.write("\n")
    outfilehaplodepth.write(sep.join(gens))
    outfilehaplodepth.write("\n")
    writebufferhaplo=""
    writebufferhaplocount=""
    writebufferhaplodepth=""
 
    count=0
    
    for cats in hapcatlist:
        count+=1
        
        haploo=haplodict2.get(str(cats),0)

        if haploo<>0:
            if "consensus" in haploo:
                linebufferhaplo=str(cats)+sep
                linebufferhaplocount=str(cats)+sep
                linebufferhaplodepth=str(cats)+sep
                for gen in gens:
                    #genhaplo=list()
                    genhaplocount=list()
                    deptheach=0
                    #eachgenhaplodict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                      
                        depthh=haplodict.get((gen,str(cats),"consensus"),0)     
         
                        genhaplowrite="consensus"
                        genhaplocountwrite=str(depthh)
                        #print str(i),gensnpwrite, gencountwrite
                        linebufferhaplo+=genhaplowrite+sep
                        linebufferhaplocount+=genhaplocountwrite+sep
                        linebufferhaplodepth+=genhaplocountwrite+sep
                    else:
                        linebufferhaplo+="NA"+sep
                        linebufferhaplocount+="NA"+sep
                        linebufferhaplodepth+="NA"+sep
                writebufferhaplo+=linebufferhaplo+"\n"
                writebufferhaplocount+=linebufferhaplocount+"\n"
                writebufferhaplodepth+=linebufferhaplodepth+"\n"
            

    
        
            else:
                linebufferhaplo=str(cats)+sep
                linebufferhaplocount=str(cats)+sep
                linebufferhaplodepth=str(cats)+sep
                for gen in gens:
                    genhaplo=list()
                    genhaplocount=list()
                    deptheach=0
                    eachgenhaplodict=dict()
                    flag=gencatdict.get((gen,str(cats)),0)
                    if flag==1:
                        for haplotype in haploo:
                            depthh=haplodict.get((gen,str(cats),haplotype),0)     
                            
                            if depthh>0:
                                eachgenhaplodict[haplotype]=depthh

                        for key, val in eachgenhaplodict.items():
                            genhaplo.append(key)
                            genhaplocount.append(str(val))
                                
                        genhaplowrite=allelesep.join(genhaplo)
                        genhaplocountwrite=allelesep.join(genhaplocount)
                        genhaplodepthwrite=str(sum([int(y) for y in genhaplocount]))
                        #print str(i),gensnpwrite, gencountwrite
                        linebufferhaplo+=genhaplowrite+sep
                        linebufferhaplocount+=genhaplocountwrite+sep
                        linebufferhaplodepth+=genhaplodepthwrite+sep
                    else:
                        linebufferhaplo+="NA"+sep
                        linebufferhaplocount+="NA"+sep
                        linebufferhaplodepth+="NA"+sep
                        
                writebufferhaplo+=linebufferhaplo+"\n"
                writebufferhaplocount+=linebufferhaplocount+"\n"
                writebufferhaplodepth+=linebufferhaplodepth+"\n"
                  

            if count%1000==0:
                noww=time()
                timee=str((noww-start)/60)
                
                print "%s rows in %s minutes" %(count,timee)               
                    
                outfilehaplo.write(writebufferhaplo)
                outfilehaplocount.write(writebufferhaplocount)
                outfilehaplodepth.write(writebufferhaplodepth)
                writebufferhaplo=""
                writebufferhaplocount=""
                writebufferhaplodepth=""

        
    noww=time()
    timee=str((noww-start)/60)
    
    print "%s rows in %s minutes" %(count,timee)               
        
    outfilehaplo.write(writebufferhaplo)
    outfilehaplocount.write(writebufferhaplocount)
    outfilehaplodepth.write(writebufferhaplodepth)
    outfilehaplo.close()
    outfilehaplocount.close()
    outfilehaplodepth.close()
    



def changeformat(inf, outf1, sep, allelesep):
    """Changes format from SNP A/C...to binary"""


    
    infile=open(inf,"r")

    header=next(infile)

    writebuffer=header
 
    
    for line in infile:
        linelist=line.strip().split(sep)
        catalog=linelist[0]
        snps=linelist[1::]
##        print linelist
##        print catalog
##        print snps
        for gensnp in snps:
            if allelesep in gensnp:
                
                nucls=gensnp.split(allelesep)
                break
##        print catalog, nucls
        for nucl in nucls:
            eachline=catalog+"_"+nucl
            
            for gensnp in snps:
                if "NA" in gensnp:
                    eachline+=",NA"
                elif nucl in gensnp:
                    eachline+=",1"
                else:
                    eachline+=",0"
##            print eachline
##            raw_input()
            writebuffer+=eachline+"\n"
    print "writing in %s.f01.csv" %inf
    outf=open(outf1,"a")
    outf.write(writebuffer)
    infile.close()
    outf.close()

def binomial_test(n, k,r):
    """Calculate binomial probability
    """
    r=0.0025
    #for diploid
    #p2 = comb(n, k) * 0.005**k * 0.995**(n-k)
    #return p2
    
    #for tetraploid
    p1 = comb(n, k) * r**k * (1-r)**(n-k)
    
    return p1
    

main()
