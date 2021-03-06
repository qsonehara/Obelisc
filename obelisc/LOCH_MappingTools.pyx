# distutils: language=c++
import cython
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.pair cimport pair
from libc.stdio cimport printf

cdef CheckLOCH(np.ndarray[double, ndim=2] geno):
    cdef:
        int i, j, flag0, flag2
        int[:] ret = np.ones(len(geno[0,:]), dtype=np.int32)
    for i in range(len(geno[0,:])):
        flag0 = 0
        flag2 = 0
        for j in range(len(geno[:,0])):
            if geno[j,i] == 0:
                flag0 = 1
            elif geno[j,i] == 2:
                flag2 = 1
            if flag0 and flag2:
                #print(i)
                ret[i] = 0
                break
    return ret

cdef object MakeLOCHonWindow(int[:] flagLOCHSNP, np.ndarray[long long, ndim=1] Position, int Windowkb, int numGapSNP, unsigned int numMinSNP, unsigned int WindowGap):
    cdef:
        unsigned int numSNP = len(Position)
        unsigned int WindowSt = 0
        unsigned int WindowEd = 1
        int[:] ret = np.zeros(numSNP, dtype=np.int32)
        int i, j, k
        int WindowFlag
        int numObsGapSNP = 2 - flagLOCHSNP[WindowSt] - flagLOCHSNP[WindowEd]
    while(WindowSt < numSNP - 1):
        while(Position[WindowEd] - Position[WindowSt] < Windowkb * 1000 and WindowEd < numSNP - 1):
            WindowEd += 1
            if flagLOCHSNP[WindowEd] == 0:
                numObsGapSNP += 1
            if Position[WindowEd] - Position[WindowEd - 1] > WindowGap * 1000:
                WindowSt, WindowEd = WindowEd, WindowEd + 1
                numObsGapSNP = 2 - flagLOCHSNP[WindowSt] - flagLOCHSNP[WindowEd]
                continue
        if WindowEd >= numSNP - 1:
            break
        WindowFlag = True
        if numObsGapSNP > numGapSNP:
            WindowFlag = False
        if WindowEd - WindowSt + 1 < numMinSNP:
            WindowFlag = False
        if WindowFlag:
            for k in range(WindowSt,WindowEd+1):
                ret[k] = 1
        if flagLOCHSNP[WindowSt] == 0:
            numObsGapSNP -= 1
        WindowSt += 1
    return ret

cdef object DecideLOCHStretch(int[:] flagLOCHWindow, np.ndarray[long long, ndim=1] Position, int Stretchkb):
    cdef:
        unsigned int numSNP = len(flagLOCHWindow)
        int numStretch = 0
        int OnStretch = 0
        int StretchSt = 0
        int StretchEd = 0
        unsigned int WindowGap = 1000
        int i = 0
        pair[int, int] StEd
        vector[pair[int, int]] ret
    
    for i in range(numSNP):
        if OnStretch == 0 and flagLOCHWindow[i] == 1:
            OnStretch = 1
            StEd.first = i
        elif OnStretch == 1 and (flagLOCHWindow[i] == 0 or Position[i] - Position[i-1] > WindowGap * 1000):
            OnStretch = 0
            StEd.second = i - 1
            if Position[StEd.second] - Position[StEd.first] >= Stretchkb * 1000:
                ret.push_back(StEd)
                numStretch += 1
            # if window is broken, a new window may start
            if flagLOCHWindow[i] == 1 and Position[i] - Position[i-1] > WindowGap * 1000:
                OnStretch = 1
                StEd.first = i
    # loop end processing
    if OnStretch == 1:
        StEd.second = i
        if Position[StEd.second] - Position[StEd.first] >= Stretchkb * 1000:
            ret.push_back(StEd)
            numStretch += 1
    if not np.asarray(ret).any():
        return np.array([[]], dtype=np.int32, ndmin=2)
    return np.array(ret, dtype=np.int32, ndmin=2)

cpdef object PointHitonStretch(int[:,:] Stretch, np.ndarray[long long, ndim=1] Position):
    cdef:
        int numStretch = len(Stretch)
        int numPoint = len(Position)
        int[:] ret = np.zeros(numPoint, dtype=np.int32)
        int i, j
    
    if not np.asarray(Stretch).any():
        return np.array(ret, dtype=np.int32)
    
    for i in range(numStretch):
        for j in range(numPoint):
            if Position[j] >= Position[Stretch[i][0]] and Position[j] <= Position[Stretch[i][1]]:
                ret[j] = 1
    return np.array(ret, dtype=np.int32)

cpdef object LOCHMappingAll(np.ndarray[double, ndim=2] geno, np.ndarray[long long, ndim=1] Position, int Windowkb=1000, int numGapSNP=1, unsigned int numMinSNP=25, int Stretchkb=1000, unsigned int WindowGap=1000):
    cdef:
        int[:] flagLOCHSNP = CheckLOCH(geno)
        int[:] flagLOCHWindow = MakeLOCHonWindow(flagLOCHSNP, Position, Windowkb, numGapSNP, numMinSNP, WindowGap)
        int[:,:] StretchLong = DecideLOCHStretch(flagLOCHWindow, Position, Stretchkb)
    return StretchLong

cpdef object MakeROHonWindowMulti(np.ndarray[double, ndim=2] geno, np.ndarray[long long, ndim=1] Position, int Windowkb, int numGapSNP,
                                  unsigned int numMinSNP, int Stretchkb, unsigned int WindowGap):
    cdef:
        int numSNP = len(Position)
        int WindowSt
        int WindowEd
        int[:,:] ret = np.zeros((geno.shape[0], geno.shape[1]), dtype=np.int32)
        int[:] ret_single = np.zeros(geno.shape[1], dtype=np.int32)
        int i, j, k
        int numSample = len(geno)
        int WindowFlag = 0
        int numObsGapSNP = 0
        int preEd
        double[:] geno_single
    
    for i in range(numSample):
        ret_single = np.zeros(geno.shape[1], dtype=np.int32)
        geno_single = geno[i]
        WindowSt = 0
        WindowEd = 1
        numObsGapSNP = <int>(geno_single[WindowSt] == 1) + <int>(geno_single[WindowEd] == 1)
        while(WindowSt < numSNP - 1):
            while(Position[WindowEd] - Position[WindowSt] < Windowkb * 1000 and WindowEd < numSNP - 1):
                WindowEd += 1
                if geno_single[WindowEd] == 1:
                    numObsGapSNP += 1
                if Position[WindowEd] - Position[WindowEd - 1] > WindowGap * 1000:
                    WindowSt, WindowEd = WindowEd, WindowEd + 1
                    numObsGapSNP = <int>(geno_single[WindowSt] == 1) + <int>(geno_single[WindowEd] == 1)
                    continue
            if WindowEd >= numSNP - 1:
                break
            WindowFlag = True
            if numObsGapSNP > numGapSNP:
                WindowFlag = False
            if WindowEd - WindowSt + 1 < numMinSNP:
                WindowFlag = False
            if WindowFlag:
                for k in range(WindowSt,WindowEd+1):
                    ret_single[k] = 1
            if geno_single[WindowSt] == 1:
                numObsGapSNP -= 1
            WindowSt += 1
        ret[i] = ret_single
    return ret

cpdef object DecideROHStretch(int[:] flagROHWindow, np.ndarray[long long, ndim=1] Position, int Stretchkb, unsigned int WindowGap):
    cdef:
        int numSNP = len(flagROHWindow)
        int numStretch = 0
        int OnStretch = 0
        int StretchSt = 0
        int StretchEd = 0
        int i
        pair[int, int] StEd
        vector[pair[int, int]] ret
    
    for i in range(numSNP):
        if OnStretch == 0 and flagROHWindow[i] == 1:
            OnStretch = 1
            StEd.first = i
            #printf("StEd.first: %I64d \n", Position[StEd.first])
        elif OnStretch == 1 and (flagROHWindow[i] == 0 or Position[i] - Position[i - 1] > WindowGap * 1000):
            OnStretch = 0
            StEd.second = i - 1
            #printf("StEd.second: %I64d \n", Position[StEd.second])
            if Position[StEd.second] - Position[StEd.first] >= Stretchkb * 1000:
                ret.push_back(StEd)
                #printf("Set region: %I64d - %I64d \n", Position[StEd.first], Position[StEd.second])
                numStretch += 1
            # if window is broken, a new window may start
            if flagROHWindow[i] == 1 and Position[i] - Position[i-1] > WindowGap * 1000:
                OnStretch = 1
                StEd.first = i
                #printf("StEd.first: %I64d \n", Position[StEd.first])
    # loop end processing
    if OnStretch == 1:
        StEd.second = i
        #printf("StEd.second: %I64d \n", Position[StEd.second])
        if Position[StEd.second] - Position[StEd.first] >= Stretchkb * 1000:
            ret.push_back(StEd)
            #printf("Set region: %I64d - %I64d \n", Position[StEd.first], Position[StEd.second])
            numStretch += 1
    if not np.asarray(ret).any():
        StEd.first = 0
        StEd.second = 0
        ret.push_back(StEd)
        return ret
    return ret

cpdef object DecideROHStretchMulti(int[:,:] flagROHWindow, np.ndarray[long long, ndim=1] Position, int Stretchkb, unsigned int WindowGap):
    cdef:
        int i
        int numSample = len(flagROHWindow)
        vector[vector[pair[int, int]]] ret
    for i in range(numSample):
        ret.push_back(DecideROHStretch(flagROHWindow[i], Position, Stretchkb, WindowGap))
    return ret
