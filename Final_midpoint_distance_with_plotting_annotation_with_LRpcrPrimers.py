#!/usr/bin/python


'''
Python script to calculate the distance between subsequent RSE oligos and plot them separately for each of the 8 resequencing regions

Long Range PCR expected amplicon will be plotted as orange lines

'''
#PlotProbesPythonR

import pprint
pp = pprint.PrettyPrinter

import matplotlib.pyplot as plt
import numpy as np

#import matplotlib as plt

regions = "RegionsList.txt"

#NoPolyNoRepeats = "NoPolyNoRepeats_ProbesFelix.txt"
#FinalSelectedProbes ="FinalSelectedProbes_RegionsRemoved031614.txt"   #WAS NOT SORTED SO GAVE WRONG DISTACE BETWEEN PROBES PLOTTING

#FinalSelectedProbes ="PlotProbes06202014.txt" #old wrong file (check 07-27-2014 evernote duplicates removed, 4 Randy picked removed)

FinalSelectedProbes ="PlotProbes07272014.txt"


#ALLPOSSIBLEPROBES = "ALLPOSSIBLEPROBES_031414.txt"

#NoPolyNoRepeats = "WithPoly.txt"

#NoPolyNoRepeats ="NoPolyNoRepeats_Probes.txt"
#WithPoly = "WithPoly.txt"
#closetoRepeats = "With_or_closetoRepeats.txt"
#NoPolyNoRepeats_new ="WithPoly_close2Repeats_final.txt"

Primers = "LRpcrPrimers.txt"
Genes = "CandidateGenesReseqRegions.txt"
AllGenesReseqRegions = "AllGenesReseqRegions.txt"

f1 = open(regions, "r")
#f2 = open(WithPoly,"r")
f2 = open(FinalSelectedProbes,"r")
f3 = open(Genes, "r")
f4 = open(AllGenesReseqRegions, "r")
f5 = open(Primers, "r")


# read data in f2
f2_rows = []
for line in f2.readlines():
    line = line.strip().split('\t')   #.strip() removes all whitespace at the start and end, including spaces, tabs, newlines and carriage returns. Leaving it in doesn't do any harm, and allows your program to deal with unexpected extra whitespace inserted into the file.
    f2_rows.append(line)
    #print f2_rows [0]

# read data in f5
f5_all = []
f5_rows1 = []
f5_rows2start = []
f5_rows2stop = []
f5_rows2 = []
f5_allstart = []
f5_allstop = []
for line3 in f5.readlines():
    line3 = line3.strip().split('\t')
    #print line3 [0]
    f5_rows1.append(line3 [0])
    #f5_rows1.append(line3 [0])
    f5_rows2start.append(line3 [1])
    f5_rows2stop.append(line3 [2])
    f5_rows2.append(line3 [1])
    f5_rows2.append(line3 [2])
f5_all.append(f5_rows1)
f5_all.append(f5_rows2)

f5_allstart.append(f5_rows1)
f5_allstart.append(f5_rows2start)
f5_allstop.append(f5_rows1)
f5_allstop.append(f5_rows2stop)
   

# read data in f4
f4_all_start = []
f4_all_stop = []
f4_start = []
f4_stop = []
f4_chr = []
f4_rows = []
for line in f4.readlines():
    line = line.strip().split('\t')   #.strip() removes all whitespace at the start and end, including spaces, tabs, newlines and carriage returns. Leaving it in doesn't do any harm, and allows your program to deal with unexpected extra whitespace inserted into the file.
    f4_chr.append(line [0])
    f4_start.append(line [1])
    f4_stop.append(line [2])
#pprint.pprint(f4_chr)
#print (f4_start)
#print (f4_stop)

f4_all_start.append(f4_chr)
f4_all_start.append(f4_start)
 
f4_all_stop.append(f4_chr)
f4_all_stop.append(f4_stop)

#print f4_all_start [0]
    
# read data in f3
f3_all = []
f3_rows1 = []
f3_rows2start = []
f3_rows2stop = []
f3_rows2 = []
f3_allstart = []
f3_allstop = []
for line3 in f3.readlines():
    line3 = line3.strip().split('\t')
    #print line3 [0]
    f3_rows1.append(line3 [0])
    #f3_rows1.append(line3 [0])
    f3_rows2start.append(line3 [1])
    f3_rows2stop.append(line3 [2])
    f3_rows2.append(line3 [1])
    f3_rows2.append(line3 [2])
f3_all.append(f3_rows1)
f3_all.append(f3_rows2)

f3_allstart.append(f3_rows1)
f3_allstart.append(f3_rows2start)
f3_allstop.append(f3_rows1)
f3_allstop.append(f3_rows2stop)
 

                
f, axarr = plt.subplots(4, 2, sharey=True)

str_line_data = ""               
str_line_print = ""
distances = []
midpoints = []

lines1 = f1.readlines()
counter_plots = 0
for line in lines1:
    if line.startswith('chr'):
        cols = line.strip().split('\t')
        chr = cols [0]
        start = int(cols[1])
        stop = cols [2]
        RegionName = cols [3]
        
        ###############################
        
        mid_point = []
        mid_point.append ((float(cols[1]) + float(cols[1]))/2)
        for i in xrange(len(f2_rows)):
            if (f2_rows[i][0][0:3]) == (cols [0][0:3]) and (int(f2_rows[i][0][3:])) == (int(cols [0][3:])) and (int(cols[1])) <= (int(f2_rows[i][1])) and (int(cols[2])) >= (int(f2_rows[i][2])):
                mid_point.append ((float(f2_rows[i][1]) + float(f2_rows[i][2]))/2)
        mid_point.append ((float(cols[2]) + float(cols[2]))/2)
        midpoints.append (mid_point)
        #print str(mid_point) + "=midpoint"
        ###############################
        
        distance = [float(0)]
        for i in xrange(len(mid_point)-1):
            distance.append (mid_point [i+1]-mid_point [i])

        threshold = []
        for i in xrange(len(mid_point)):
            threshold.append (2000)
        
#############################################################################   my_data[my_data == 0.0] = np.nan
        
        allgene_coordystart = []
        #gene_coordystart = np.nan
        for i in xrange((len(f4_all_start [0]))):
            #print (f3_allstop [1] [1])
            #print cols[1]
            if ((int(f4_all_start [0] [i])) == (int(cols [0][3:]))) and (((int(cols[1])) <= (int(f4_all_start [1] [i])))) and ((int(cols[2]) >= (int(f4_all_stop [1] [i])))):
                #print (f3_allstart [1] [i])
                allgene_coordystart.append(int(f4_all_start [1] [i]))
                
               
        allgene_coordystop = []
        #gene_coordystop = np.nan
        for i in xrange((len(f4_all_start [0]))):
            #print (f3_allstop [1] [1])
            #print cols[1]
           if ((int(f4_all_start [0] [i])) == (int(cols [0][3:]))) and (((int(cols[1])) <= (int(f4_all_start [1] [i])))) and ((int(cols[2]) >= (int(f4_all_stop [1] [i])))):
                #print (f3_allstart [1] [i])
                allgene_coordystop.append(int(f4_all_stop [1] [i]))
              
        
        
        allgene_coordxstart = []
        for i in xrange((len(allgene_coordystart))):
            allgene_coordxstart.append (-6000)
            #gene_coordxstart[gene_coordxstart == 0.0] = np.nan    
            
        allgene_coordxstop = []
        for i in xrange((len(allgene_coordystop))):
            allgene_coordxstop.append (-6000)
            #gene_coordxstop[gene_coordxstop == 0.0] = np.nan
            
            
        #print len(allgene_coordystart)


#############################################################################   my_data[my_data == 0.0] = np.nan
#Candidate genes
############################################################################# 
        gene_coordystart = []
        #gene_coordystart = np.nan
        for i in xrange((len(f3_allstop [0]))):
            #print (f3_allstop [1] [1])
            #print cols[1]
            if ((int(f3_allstop [0] [i])) == (int(cols [0][3:]))) and (((int(cols[1])) <= (int(f3_allstart [1] [i])))) and ((int(cols[2]) >= (int(f3_allstop [1] [i])))):
                #print (f3_allstart [1] [i])
                gene_coordystart.append(int(f3_allstart [1] [i]))
                
               
        gene_coordystop = []
        #gene_coordystop = np.nan
        for i in xrange((len(f3_allstop [0]))):
            #print (f3_allstop [1] [1])
            #print cols[1]
            if ((int(f3_allstop [0] [i])) == (int(cols [0][3:]))) and ((int(cols[1]) <= (int(f3_allstart [1] [i])))) and ((int(cols[2]) >= (int(f3_allstop [1] [i])))):
                #print (f3_allstart [1] [i])
                gene_coordystop.append(int(f3_allstop [1] [i])) 
                
        
        
        gene_coordxstart = []
        for i in xrange((len(gene_coordystart))):
            gene_coordxstart.append (-6000)
            #gene_coordxstart[gene_coordxstart == 0.0] = np.nan    
            
        gene_coordxstop = []
        for i in xrange((len(gene_coordystop))):
            gene_coordxstop.append (-6000)
            #gene_coordxstop[gene_coordxstop == 0.0] = np.nan

#############################################################################  
#Primers
############################################################################# 
        primer_coordystart = []
        #gene_coordystart = np.nan
        for i in xrange((len(f5_allstop [0]))):
            #print (f5_allstop [1] [1])
            #print cols[1]
            if ((int(f5_allstop [0] [i])) == (int(cols [0][3:]))) and (((int(cols[1])) <= (int(f5_allstart [1] [i])))) and ((int(cols[2]) >= (int(f5_allstop [1] [i])))):
                #print (f5_allstart [1] [i])
                primer_coordystart.append(int(f5_allstart [1] [i]))
                
               
        primer_coordystop = []
        #gene_coordystop = np.nan
        for i in xrange((len(f5_allstop [0]))):
            #print (f5_allstop [1] [1])
            #print cols[1]
            if ((int(f5_allstop [0] [i])) == (int(cols [0][3:]))) and ((int(cols[1]) <= (int(f5_allstart [1] [i])))) and ((int(cols[2]) >= (int(f5_allstop [1] [i])))):
                #print (f5_allstart [1] [i])
                primer_coordystop.append(int(f5_allstop [1] [i])) 
                
        
        
        primer_coordxstart = []
        for i in xrange((len(primer_coordystart))):
            primer_coordxstart.append (40000)
            #gene_coordxstart[gene_coordxstart == 0.0] = np.nan    
            
        primer_coordxstop = []
        for i in xrange((len(primer_coordystop))):
            primer_coordxstop.append (40000)
            #gene_coordxstop[gene_coordxstop == 0.0] = np.nan

###########################################################################
            
        # Two subplots, the axes array is 1-d
#        
        axarr[counter_plots%4][counter_plots/4].scatter(mid_point, distance, color='grey', label='Probe distance')
        axarr[counter_plots%4][counter_plots/4].plot(mid_point, threshold, color= 'red', lw=1, alpha=0.4, label='Threshold(2kb)')

        
        axarr[counter_plots%4][counter_plots/4].scatter(allgene_coordystart, allgene_coordxstart, color= 'g', label='All Genes start_pos',marker = "<")
        axarr[counter_plots%4][counter_plots/4].scatter(allgene_coordystop, allgene_coordxstop, color= 'c', label='All Genes stop_pos', marker = ">")
        for i,j,k in zip(allgene_coordxstart,allgene_coordystart, allgene_coordystop ):
            axarr[counter_plots%4][counter_plots/4].plot([j,k],[i,i],lw=2, alpha=0.4, color = 'c')
        axarr[counter_plots%4][counter_plots/4].scatter(gene_coordystart, gene_coordxstart, color= 'red', label='Target Genes start_pos', marker = "<")
        axarr[counter_plots%4][counter_plots/4].scatter(gene_coordystop, gene_coordxstop, color= 'blue', label='Target Genes stop_pos', marker = ">")
 
        axarr[counter_plots%4][counter_plots/4].scatter(primer_coordystart, primer_coordxstart, color= 'red', label='LR PCR primer start_pos', marker = "<")
        axarr[counter_plots%4][counter_plots/4].scatter(primer_coordystop, primer_coordxstop, color= 'orange', label='LR PCR primer stop_pos', marker = ">")
        for x,y,z in zip(primer_coordxstart,primer_coordystart, primer_coordystop ):
            axarr[counter_plots%4][counter_plots/4].plot([y,z],[x,x],lw=2, alpha=0.4, color = 'orange') 
        
        axarr[counter_plots%4][counter_plots/4].set_title(RegionName)

        
        plt.xlabel('Midpoint')
        plt.ylabel('Probe distance')
        plt.ylim([-10000,70000])
        #plt.xlim([(int(cols[1])),(int(cols[1]) + 303000)])
        



        counter_plots += 1

plt.legend(loc='upper right', prop={'size':4})
#plt.suptitle('Distance from previous probe / (region start/stop), WithPol_Probes')
plt.suptitle('Distance from previous probe / (region start/stop), SelectedProbes, LR PCR Primers 082914')
#plt.suptitle('Distance from previous probe / (region start/stop), ALLPOSSIBLEPROBES_old')
plt.show()




f1.close()
f2.close()
f3.close()
f4.close()
f5.close()




