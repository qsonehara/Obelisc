#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Obelisc (Observational linkage scan) is a command line tools for
    1. SNP streak-based IBD mapping
    2. Runs of homozygosity mapping
    3. Visualizing both mapping results
'''
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pysnptools.snpreader import Bed
import copy
from datetime import datetime
import argparse
import LOCH_MappingTools

def replace(x):
    '''Replaces non-numeric chromosome ID'''
    if x == 'X':
        return 23
    elif x == 'Y':
        return 24
    elif x == 'XY':
        return 25
    elif x == 'MT':
        return 26
    return x

class Mapping():
    def __init__(self, prefix, case_file):
        self.prefix = prefix
        self.case_file = case_file
        self.snpreader = Bed(f"{prefix}.bed", count_A1=False)
        if self.snpreader.pos.dtype != 'int64':
            self.snpreader.pos[:,0] = np.vectorize(replace)(self.snpreader.pos[:,0])
        self.snpreader.pos[:,1] = self.snpreader.pos[:,0] * 100000000000 + self.snpreader.pos[:,2]
        self.snpdata = self.snpreader.read()
        print('SNP data loaded.')
        self.chr_list = list(set(self.snpreader.pos[:,0]))
        self.Chr = self.snpreader.pos[:,0]
        self.Position =  self.snpreader.pos[:,1]
        self.bp =  self.snpreader.pos[:,2]
        self.SNPID = self.snpreader.sid
        self.case = np.loadtxt(case_file, dtype=self.snpreader.iid.dtype)[:,:2]
        self.case_list = list(self.case)
        self.all_list = list([tuple(x) for x in self.snpreader.iid])
        self.caseset = set([tuple(x) for x in self.case])
        self.control_list = [list(x) for x in self.all_list if x not in self.caseset]
        self.numSNP = self.snpreader.sid_count
        self.numSample = len(self.all_list)
        self.numCase = len(self.case_list)
        self.numControl = len(self.control_list)
        self.case_geno = self.snpdata.val[self.snpreader.iid_to_index(self.case)]
        L = []
        for i in self.case_list:
            L.append(i[1].decode('utf-8'))
        self.case_list_print = '\n'.join(L)
        print('Case individuals are: \n')
        print(self.case_list_print)
        print('\n')

    def ibdmapping_gw(self, Windowkb, Stretchkb, numGapSNP, numMinSNP, WindowGap, out, point):
        '''Performs genome-wide IBD mapping based on SNP streak'''
        print("********************\n"
              f"IBD mapping started.\ninput file prefix: {self.prefix}\nWindowkb: {Windowkb}\nStretchkb: "
              f"{Stretchkb}\nnumGapSNP: {numGapSNP}\nnumMinSNP: {numMinSNP}\nWindowGap: {WindowGap}\noutput file prefix: {out}")
        with open(f"{out}.txt", 'w') as f:
            out = (f"Log_for_NonparametricIBDmapping\n\nInput_Genotype_File:\t{self.prefix}.bim/fam/bed\nInput_Case_File:\t{self.case_file}"
                   f"\nNo.SNP:\t{self.numSNP}\nWindow_kb:\t{Windowkb}\nStretch_kb:\t{Stretchkb}\nWindowGap_kb:\t{WindowGap}"
                   f"\nNo.InconsistentSNP:\t{numGapSNP}\nNo.MinSNP:\t{numMinSNP}\n\nNo.Samples:\t{self.numSample}\nNo.Cases:\t{self.numCase}\n"
                   f"\nCase individuals are: \n{self.case_list_print}\n\n")
            f.write(out)
            StretchLong = LOCH_MappingTools.LOCHMappingAll(self.case_geno, self.Position, Windowkb, numGapSNP, numMinSNP, Stretchkb, WindowGap)
            if point:
                PointHitFlag = LOCH_MappingTools.PointHitonStretch(StretchLong, self.Position)
            numStretchLong = len(StretchLong) if len(StretchLong[0]) else 0
            out = 'No.IBD stretch CaseOnly:\t{0}'.format(numStretchLong)
            print(out)
            out = (f"\nNo.IBD_Stretch_in_All_Cases:\t{numStretchLong}\nIBD_stretch\tChr\tStart_SNP\tEnd_SNP"
                   f"\tStart_Position(bp)\tEnd_Position(bp)\tLength(bp)\n")
            f.write(out)
            for i in range(numStretchLong):
                start = StretchLong[i][0]
                end = StretchLong[i][1]
                L = [str(i+1), str(self.Chr[start]), self.SNPID[start].decode('utf-8'), self.SNPID[end].decode('utf-8'),
                     str(self.bp[start]), str(self.bp[end]), f"{(self.bp[end] - self.bp[start])}\n"]
                out = '\t'.join(L)
                f.write(out)
            if point:
                PointHitResult = [[0] * self.numSNP for i in range(self.numControl)]
            StretchLongControl = [None] * self.numControl
            if self.numControl:
                for i in range(self.numControl):
                    case_control_list = self.case_list + [self.control_list[i]]
                    case_control_geno = self.snpdata.val[self.snpreader.iid_to_index(case_control_list)]
                    StretchLongControl[i] = LOCH_MappingTools.LOCHMappingAll(case_control_geno, self.Position, Windowkb,
                                                                            numGapSNP, numMinSNP, Stretchkb, WindowGap)
                    if point:
                        PointHitFlagControl = LOCH_MappingTools.PointHitonStretch(StretchLongControl[i], self.Position)
                        PointHitResult[i] = PointHitFlagControl
                    numStretchLongControl = len(StretchLongControl[i]) if len(StretchLongControl[i][0]) else 0
                    out = (f"\nNo.IBD_stretch_in_All_Cases_and_1_Control({self.control_list[i][1].decode('utf-8')}):\t{numStretchLongControl}"
                        f"\nIBD_stretch\tChr\tStart_SNP\tEnd_SNP\tStart_Position(bp)\tEnd_Position(bp)\tLength(bp)\n")
                    f.write(out)
                    numStretchLongControl = len(StretchLongControl[i])
                    if not StretchLongControl[i][0]:
                        numStretchLongControl = 0
                    for j in range(numStretchLongControl):
                        start = StretchLongControl[i][j][0]
                        end = StretchLongControl[i][j][1]
                        L = [str(j+1), str(self.Chr[start]), self.SNPID[start].decode('utf-8'), self.SNPID[end].decode('utf-8'),
                            str(self.bp[start]), str(self.bp[end]), f"{(self.bp[end] - self.bp[start])}\n"]
                        out = '\t'.join(L)
                        f.write(out)
            print('IBD mapping for Cases and Controls finished!!')

            if point:    
                PointHitSumControl = [0] * len(self.Position)
                for i in range(len(PointHitResult)):
                    for j in range(len(PointHitResult[0])):
                        PointHitSumControl[j] += PointHitResult[i][j]
                out = ("\n\nIndividual_Marker/Sample_IBDstatus_in_IBDregions\n\nSNP\tChr\tbp\tNo.Cases_in_IBD"
                      "\tNo.Controls_in_IBD\tAll_Cases_in_IBD(Yes:1/No:0)")
                f.write(out)
                L = []
                for i in range(len(self.control_list)):
                    L.append("\t{0}_in_IBD(Yes:1/No:0)".format(self.control_list[i][1].decode('utf-8')))
                L.append('\n')
                out = ''.join(L)
                f.write(out)
                for j in range(len(PointHitFlag)):
                    if PointHitFlag[j]:
                        L = [f"{self.SNPID[j].decode('utf-8')}\t{str(self.Chr[j])}\t{str(self.bp[j])}"
                        f"\t{str(len(self.case_list))}\t{str(PointHitSumControl[j])}\t{PointHitFlag[j]}"]
                        for i in range(len(PointHitResult)):
                            L.append(f'\t{PointHitResult[i][j]}')
                        L.append('\n')
                        out = ''.join(L)
                        f.write(out)

            case_regions = np.asarray(StretchLong)
            if self.numControl:
                ctrl_regions = [np.asarray(i) for i in StretchLongControl]
                num_ctrl_regions = len(ctrl_regions)
                if not len(ctrl_regions[0][0]):
                    num_ctrl_regions = 0
            else:
                num_ctrl_regions = 0
            num_case_regions = len(case_regions)
            if not len(case_regions[0]):
                num_case_regions = 0
            edges = []
            for i in range(num_case_regions):
                edges.append([case_regions[i][0], 0, 0])
                edges.append([case_regions[i][1], 0, 1])
            for i in range(num_ctrl_regions):
                num_ctrl_region = len(ctrl_regions[i])
                if not len(ctrl_regions[i][0]):
                    num_ctrl_region = 0
                for j in range(num_ctrl_region):
                    edges.append([ctrl_regions[i][j][0], i+1, 0])
                    edges.append([ctrl_regions[i][j][1], i+1, 1])
            state = [0] * (num_ctrl_regions+1)
            edges = sorted(edges)
            self.ibd_regions = []
            for i in range(len(edges)):
                if not edges[i][2]:
                    state[edges[i][1]] += 1
                else:
                    state[edges[i][1]] -= 1
                if i+1 == len(edges) or edges[i][0] != edges[i+1][0]:
                    self.ibd_regions.append([edges[i][0], copy.copy(state)])
            L = ["\n\nChr\tStart\tEnd\tIBD_in_Controls\tIBD_Case_specificity"]
            for i in range(len(self.ibd_regions) - 1):
                if self.ibd_regions[i][1][0]!=0:
                    prop =  round(1.0-sum(self.ibd_regions[i][1][1:])/self.numControl, 3) if self.numControl else 1.0
                    L.append(f"{self.Chr[self.ibd_regions[i][0]]}\t{self.bp[self.ibd_regions[i][0]]}"
                             f"\t{self.bp[self.ibd_regions[i+1][0]]}\t{sum(self.ibd_regions[i][1][1:])}\t{prop}")
            L.append('\nCalculation_finished_at:\t{}\n'.format(datetime.now().strftime("%Y/%m/%d %H:%M:%S")))
            out = '\n'.join(L)
            f.write(out)

    def rohmapping_gw(self, Windowkb, Stretchkb, numGapSNP, numMinSNP, WindowGap, out):
        '''Performs genome-wide runs of homozygosity mapping'''
        print("********************\n"
              f"ROH mapping started.\ninput file prefix: {self.prefix}\nWindowkb: {Windowkb}\nStretchkb: {Stretchkb}"
              f"\nnumGapSNP: {numGapSNP}\nnumMinSNP: {numMinSNP}\nWindowGap: {WindowGap}\noutput file prefix: {out}")
        with open(f"{out}.txt", 'w') as f:
            out = (f"Log_for_ROHmapping\n\nInput_File:\t{self.prefix}.ped/map/info/case\nNo.SNP:\t{self.numSNP}\nWindow_kb:\t{Windowkb}"
                   f"\nStretch_kb:\t{Stretchkb}\nWindowGap_kb:\t{WindowGap}\nNo.InconsistentSNP:\t{numGapSNP}"
                   f"\nNo.MinSNP:\t{numMinSNP}\n\nNo.Samples:\t{self.numSample}\nNo.Cases:\t{self.numCase}\n\n")
            f.write(out)

            ROHwin = LOCH_MappingTools.MakeROHonWindowMulti(self.case_geno, self.Position, Windowkb, numGapSNP, numMinSNP, Stretchkb, WindowGap)
            StretchLongArray = LOCH_MappingTools.DecideROHStretchMulti(ROHwin, self.Position, Stretchkb)
            
            for i in range(self.numCase):    
                numStretchLongCase = len(StretchLongArray[i]) if (StretchLongArray[i][0][1]) else 0
                out = (f"\n\nNo.ROH_in_1_Case({self.case_list[i][1].decode('utf-8')}):\t{numStretchLongCase}"
                       f"\nROH\tChr\tStart_SNP\tEnd_SNP\tStart_Position(bp)\tEnd_Position(bp)\tLength(bp)\n")
                f.write(out)
                for j in range(len(StretchLongArray[i])):
                    start = StretchLongArray[i][j][0]
                    end = StretchLongArray[i][j][1]
                    if end == 0:
                        continue
                    L = [str(j+1), str(self.Chr[start]), self.SNPID[start].decode('utf-8'), self.SNPID[end].decode('utf-8'), 
                         str(self.bp[start]), str(self.bp[end]), f"{(self.bp[end] - self.bp[start])}\n"]
                    out = '\t'.join(L)
                    f.write(out)

            if self.numControl:
                control_geno = self.snpdata.val[self.snpreader.iid_to_index(self.control_list)]
                ROHwinControl = LOCH_MappingTools.MakeROHonWindowMulti(control_geno, self.Position, Windowkb, numGapSNP, numMinSNP, Stretchkb, WindowGap)
                StretchLongControl = LOCH_MappingTools.DecideROHStretchMulti(ROHwinControl, self.Position, Stretchkb)

            for i in range(self.numControl):    
                numStretchLongControl = len(StretchLongControl[i]) if StretchLongControl[i][0][1] else 0
                out = (f"\n\nNo.ROH_in_1_Control({self.control_list[i][1].decode('utf-8')}):\t{numStretchLongControl}"
                       f"\nROH\tChr\tStart_SNP\tEnd_SNP\tStart_Position(bp)\tEnd_Position(bp)\tLength(bp)\n")
                f.write(out)
                for j in range(len(StretchLongControl[i])):
                    start = StretchLongControl[i][j][0]
                    end = StretchLongControl[i][j][1]
                    if end == 0:
                        continue
                    L = [str(j+1), str(self.Chr[start]), self.SNPID[start].decode('utf-8'), self.SNPID[end].decode('utf-8'), 
                         str(self.bp[start]), str(self.bp[end]), f"{(self.bp[end] - self.bp[start])}\n"]
                    out = '\t'.join(L)
                    f.write(out)    
            print('ROH detection for Cases and Controls finished!!')

            case_regions = [np.asarray(i) for i in StretchLongArray]
            edges = []
            for i in range(len(case_regions)):
                for j in range(len(case_regions[i])):
                    edges.append([case_regions[i][j][0], i, 0])
                    edges.append([case_regions[i][j][1], i, 1])
            state = [0] * (len(case_regions))
            edges = sorted(edges)
            self.roh_case_regions = []
            for i in range(len(edges)):
                if not edges[i][2]:
                    state[edges[i][1]] += 1
                else:
                    state[edges[i][1]] -= 1
                if (i+1 == len(edges) or edges[i][0] != edges[i+1][0]):
                    self.roh_case_regions.append([edges[i][0], copy.copy(state)])
            L = ["\n\nChr\tStart\tEnd\tNumber_of_ROH_in_Cases\n"]
            for i in range(len(self.roh_case_regions) - 1):
                if np.sum(self.roh_case_regions[i][1]):
                    L.append(f"{self.Chr[self.roh_case_regions[i][0]]}\t{self.bp[self.roh_case_regions[i][0]]}"
                             f"\t{self.bp[self.roh_case_regions[i+1][0]]}\t{np.sum(self.roh_case_regions[i][1])}\n")
            out = ''.join(L)
            f.write(out)

            if self.numControl:
                control_regions = [np.asarray(i) for i in StretchLongControl]
                edges = []
                for i in range(len(control_regions)):
                    for j in range(len(control_regions[i])):
                        edges.append([control_regions[i][j][0], i, 0])
                        edges.append([control_regions[i][j][1], i, 1])
                state = [0] * (len(control_regions))
                edges = sorted(edges)
                self.roh_control_regions = []
                for i in range(len(edges)):
                    if not edges[i][2]:
                        state[edges[i][1]] += 1
                    else:
                        state[edges[i][1]] -= 1
                    if (i+1 == len(edges) or edges[i][0] != edges[i+1][0]):
                        self.roh_control_regions.append([edges[i][0], copy.copy(state)])
                L = ["\n\nChr\tStart\tEnd\tNumber_of_ROH_in_Controls\n"]
                for i in range(len(self.roh_control_regions) - 1):
                    if np.sum(self.roh_control_regions[i][1]):
                        L.append(f"{self.Chr[self.roh_control_regions[i][0]]}\t{self.bp[self.roh_control_regions[i][0]]}"
                                f"\t{self.bp[self.roh_control_regions[i+1][0]]}\t{np.sum(self.roh_control_regions[i][1])}\n")
                out = ''.join(L)
                f.write(out)

            out = '\nCalculation_finished_at:\t{}\n'.format(datetime.now().strftime("%Y/%m/%d %H:%M:%S"))
            f.write(out)

    def draw_diagram(self, fig_name):
        '''Draws diagrams of mapping results'''
        chr_list = np.unique(self.Chr)
        chr_length = pd.DataFrame(self.snpdata.pos).groupby(0, as_index=False).max()
        size = len(chr_list)
        gap_length = 10 ** 7
        array = np.triu(np.ones((size,size))).T
        gap = np.full(size, gap_length)
        global_len = np.dot(array, chr_length[2] + gap)
        global_len = np.insert(global_len, 0, 0)
        global_len += gap_length
        fig = plt.figure(figsize=(15,5))
        cmap = plt.get_cmap("tab10")
        ax1 = fig.add_subplot(211)
        ax2 = fig.add_subplot(212)
        StChr, EdChr = 0, 22
        Stbp, Edbp = 0, global_len[EdChr]
        [i.set_xlim(Stbp, Edbp) for i in fig.get_axes()]
        ax1.set_ylim(-0.1, 1.0)
        ax1.set_xticklabels([]) 
        ax1.set_ylabel('IBD Case specificity')
        ax1.title.set_text('IBD mapping results based on SNP streak principle')
        ax1.spines['bottom'].set_visible(False)
        ax1.spines['top'].set_visible(False)
        ax1.tick_params(axis='both', which='both', length=0)
        ax1.axhline(0, color='k', linewidth=0.5)
        ax1.axhline(1, color='k', linewidth=0.5, ls='--')
        for i in range(len(chr_list)):
            ax1.add_patch(plt.Rectangle(xy=[global_len[i], -0.1], width=chr_length[2][i], height=0.05, color=cmap(i%10)))
        for i, reg in enumerate(self.ibd_regions):
            if reg[1][0]!=0:
                prop = 1 - sum(reg[1][1:])/self.numControl if self.numControl else 1.0
                ax1.add_patch(plt.Rectangle(xy=[(self.bp[reg[0]] + global_len[self.Chr[reg[0]] - 1]), 0],
                              width=self.bp[self.ibd_regions[i+1][0]] - self.bp[reg[0]], 
                              height=prop, color=cmap((self.Chr[reg[0]] - 1) % 10)))
                ax1.add_patch(plt.Rectangle(xy=[(self.bp[reg[0]] + global_len[self.Chr[reg[0]] - 1]), -0.04], 
                              width=max(self.bp[self.ibd_regions[i+1][0]] - self.bp[reg[0]], (Edbp - Stbp)*0.001),
                              height=0.03, color='r'))
        ax2.set_ylim(-self.numCase*0.1, self.numCase)
        ax2.set_xticks([(global_len[i] + global_len[i+1])/2 for i in range(StChr, EdChr)])
        ax2.set_xticklabels([i for i in np.unique(self.snpdata.pos[:, 0])])
        ax2.set_ylabel('ROH in Cases')
        ax2.title.set_text('Runs of homozygosity detection results')
        ax2.spines['bottom'].set_visible(False)
        ax2.spines['top'].set_visible(False)
        ax2.tick_params(axis='both', which='both', length=0)
        ax2.axhline(0, color='k', linewidth=0.5)
        ax2.axhline(self.numCase, color='k', linewidth=0.5, ls='--')
        for i in range(len(chr_list)):
            ax2.add_patch(plt.Rectangle(xy=[global_len[i], -self.numCase * 0.1], width=chr_length[2][i],
                                        height=self.numCase * 0.05, color=cmap(i%10)))
        for i, reg in enumerate(self.roh_case_regions):
            numROH = np.sum(reg[1])
            if numROH:
                ax2.add_patch(plt.Rectangle(xy=[(self.bp[reg[0]] + global_len[self.Chr[reg[0]] - 1]), 0], 
                              width=self.bp[self.roh_case_regions[i+1][0]] - self.bp[reg[0]], height=numROH, 
                              color=cmap((self.Chr[reg[0]] - 1) % 10)))
                if numROH == self.numCase:
                    ax2.add_patch(plt.Rectangle(xy=[(self.bp[reg[0]] + global_len[self.Chr[reg[0]] - 1]), -0.03 * numROH], 
                                width=max(self.bp[self.roh_case_regions[i+1][0]] - self.bp[reg[0]], (Edbp - Stbp)*0.001),
                                height=0.02 * numROH, color='r'))
        fig.align_labels()
        plt.savefig(fig_name)

def main():
    parser = argparse.ArgumentParser(description='Mapping genome-wide IBD regions based on SNP streak principle')
    parser.add_argument('prefix', type=str, help='plink BED file prefix')
    parser.add_argument('-c', '--case', dest='case_file', type=str, help='Case file')
    parser.add_argument('-w', '--ibd-win-kb', dest='ibd_windowkb', type=int, default=1500, help='IBD Windowkb') 
    parser.add_argument('-s', '--ibd-str-kb', dest='ibd_stretchkb', type=int, default=1500, help='IBD Stretchkb') 
    parser.add_argument('-i', '--ibd-inconsistent-snp', dest='ibd_numgapsnp', type=int, default=1, help='IBD numInconsistentSNP') 
    parser.add_argument('-m', '--ibd-min-snp', dest='ibd_numminsnp', type=int, default=25, help='IBD numMinSNP') 
    parser.add_argument('-g', '--ibd-win-gap', dest='ibd_windowgap', type=int, default=1000, help='IBD WindowGap') 
    parser.add_argument('-o', '--ibd-out', dest='ibd_out' , type=str, help='IBD output file prefix')
    parser.add_argument('-W', '--roh-win-kb', dest='roh_windowkb', type=int, default=1500, help='ROH Windowkb') 
    parser.add_argument('-S', '--roh-str-kb', dest='roh_stretchkb', type=int, default=1500, help='ROH Stretchkb') 
    parser.add_argument('-I', '--roh-inconsistent-snp', dest='roh_numgapsnp', type=int, default=1, help='ROH numInconsistentSNP') 
    parser.add_argument('-M', '--roh-min-snp', dest='roh_numminsnp', type=int, default=25, help='ROH numMinSNP') 
    parser.add_argument('-G', '--roh-win-gap', dest='roh_windowgap', type=int, default=1000, help='ROH WindowGap') 
    parser.add_argument('-O', '--roh-out', dest='roh_out', type=str, help='ROH output file prefix')
    parser.add_argument('-p', '--point', action='store_true', help='Show the statuses of all of the IBD markers in the cases')
    parser.add_argument('-f', '--fig', dest='fig_name', type=str, help='Output figure name')
    args = parser.parse_args()
    if args.case_file is None:
        args.case_file = f"{args.prefix}.case"
    mapping = Mapping(args.prefix, args.case_file)
    if args.ibd_out is None:
        args.ibd_out = f"{args.prefix}_IBDmapping_w{args.ibd_windowkb}_s{args.ibd_stretchkb}_i{args.ibd_numgapsnp}_m{args.ibd_numminsnp}_g{args.ibd_windowgap}"
    mapping.ibdmapping_gw(args.ibd_windowkb, args.ibd_stretchkb, args.ibd_numgapsnp, args.ibd_numminsnp, args.ibd_windowgap, args.ibd_out, args.point)
    if args.roh_out is None:
        args.roh_out = f"{args.prefix}_ROHmapping_W{args.roh_windowkb}_S{args.roh_stretchkb}_I{args.roh_numgapsnp}_M{args.roh_numminsnp}_G{args.roh_windowgap}"
    if args.fig_name is None:
        args.fig_name = (f"{args.prefix}_IBD_w{args.ibd_windowkb}_s{args.ibd_stretchkb}_i{args.ibd_numgapsnp}_m{args.ibd_numminsnp}_g{args.ibd_windowgap}"
                         f"_ROH_W{args.roh_windowkb}_S{args.roh_stretchkb}_I{args.roh_numgapsnp}_M{args.roh_numminsnp}_G{args.roh_windowgap}.pdf")
    mapping.rohmapping_gw(args.roh_windowkb, args.roh_stretchkb, args.roh_numgapsnp, args.roh_numminsnp, args.roh_windowgap, args.roh_out)
    mapping.draw_diagram(args.fig_name)

if __name__ == "__main__":
    main()