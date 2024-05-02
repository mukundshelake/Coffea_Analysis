from coffea.util import load, save
import os


outputDir = 'outputs'
coffeaFile = "skimmerOutput.coffea"


# out = load("Output.coffea")
out = load(os.path.join(outputDir,coffeaFile))

Results = {}
for era in out:
    NttbarEvents = out[era]['entries']
    hist = out[era]['yMatrix']
    Results[era] = {
        "Full_y0_range" : {},
        "Y0_variations" : {}
    }
    
    for isFlow in [True, False]:
        if isFlow == False:
            featureName = "without_flow"

        if isFlow == True:
            featureName = "with_flow"

        Nuubar_highx1_Events = hist[sum, True, True, sum, sum, sum, :, :, :, :].sum(isFlow)
        Nuubar_highx2_Events = hist[sum, False, True, sum, sum, sum, :, :, :, :].sum(isFlow)
        Nubaru_highx1_Events = hist[sum, True, sum, sum, True, sum, :, :, :, :].sum(isFlow)
        Nubaru_highx2_Events = hist[sum, False, sum, sum, True, sum, :, :, :, :].sum(isFlow)
        Nddbar_highx1_Events = hist[sum, True, sum, True, sum, sum, :, :, :, :].sum(isFlow)
        Nddbar_highx2_Events = hist[sum, False, sum, True, sum, sum, :, :, :, :].sum(isFlow)
        Ndbard_highx1_Events = hist[sum, True, sum, sum, sum, True, :, :, :, :].sum(isFlow)
        Ndbard_highx2_Events = hist[sum, False, sum, sum, sum, True, :, :, :, :].sum(isFlow)

        N_highYt_Events = hist[True, sum, sum, sum, sum, sum, :, :, :, :].sum(isFlow)
        N_lowYt_Events = hist[False, sum, sum, sum, sum, sum, :, :, :, :].sum(isFlow)

        # print(N_highYt_Events, N_lowYt_Events)


        Nuu = Nuubar_highx1_Events + Nuubar_highx2_Events + Nubaru_highx1_Events + Nubaru_highx2_Events
        Ndd = Nddbar_highx1_Events + Nddbar_highx2_Events + Ndbard_highx1_Events + Ndbard_highx2_Events

        Fu = Nuu/NttbarEvents
        Fd = Ndd/NttbarEvents

        Fu_ud = Nuu/(Nuu + Ndd)
        Fd_ud = Ndd/(Nuu + Ndd)

        xqHigher_u = Nuubar_highx1_Events + Nubaru_highx2_Events
        xqbarHigher_u = Nuubar_highx2_Events + Nubaru_highx1_Events 

        xqHigher_d = Nddbar_highx1_Events + Ndbard_highx2_Events
        xqbarHigher_d = Nddbar_highx2_Events + Ndbard_highx1_Events 

        Du = (xqHigher_u - xqbarHigher_u)/(xqHigher_u + xqbarHigher_u)
        Dd = (xqHigher_d - xqbarHigher_d)/(xqHigher_d + xqbarHigher_d)

        Ac = (N_highYt_Events - N_lowYt_Events)/(N_highYt_Events + N_lowYt_Events)
        # print(Ac)

        Results[era]['Full_y0_range'][featureName] = {
                'Fu' : Fu,
                'Fd' : Fd,
                'Fu_ud' : Fu_ud,
                'Fd_ud' : Fd_ud,
                'Du' : Du,
                'Dd' : Dd,
                'FuDu' : Fu*Du,
                'FdDd' : Fd*Dd,
                'Ac' : Ac
                }
    miny0 = 0.0j
    maxy0 = 2.5j
    for i in range(24):
        y0 = (i+1)*0.1*1j

        ARange = [[miny0, y0], [miny0, y0]]
        BRange = [[y0, maxy0], [miny0, y0]]
        CRange = [[y0, maxy0], [y0, maxy0]]
        DRange = [[miny0, y0], [y0, maxy0]]

        Results[era]['Y0_variations'][str(0.1*(i+1))] = {}

        for region in ['A', 'B', 'C', 'D']:
            if region == 'A':
                Range = ARange
            if region == 'B':
                Range = BRange
            if region == 'C':
                Range = CRange
            if region == 'D':
                Range = DRange

            Nuubar_highx1_Events = hist[sum, True, True, sum, sum, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Nuubar_highx2_Events = hist[sum, False, True, sum, sum, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Nubaru_highx1_Events = hist[sum, True, sum, sum, True, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Nubaru_highx2_Events = hist[sum, False, sum, sum, True, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Nddbar_highx1_Events = hist[sum, True, sum, True, sum, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Nddbar_highx2_Events = hist[sum, False, sum, True, sum, sum, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Ndbard_highx1_Events = hist[sum, True, sum, sum, sum, True, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()
            Ndbard_highx2_Events = hist[sum, False, sum, sum, sum, True, Range[0][0]:Range[0][1], Range[1][0]:Range[1][1], sum, sum].sum()

            N_Yt_LowerthanY0_Events = hist[sum, sum, sum, sum, sum, sum, :y0, sum, sum, sum].sum()
            N_Ytbar_LowerthanY0_Events = hist[sum, sum, sum, sum, sum, sum, sum, :y0, sum, sum].sum()

            N_Yt_HigherthanY0_Events = hist[sum, sum, sum, sum, sum, sum, y0:, sum, sum, sum].sum()
            N_Ytbar_HigherthanY0_Events = hist[sum, sum, sum, sum, sum, sum, sum, y0:, sum, sum].sum()
       
            Nuu = Nuubar_highx1_Events + Nuubar_highx2_Events + Nubaru_highx1_Events + Nubaru_highx2_Events
            Ndd = Nddbar_highx1_Events + Nddbar_highx2_Events + Ndbard_highx1_Events + Ndbard_highx2_Events

            Fu = Nuu/NttbarEvents
            Fd = Ndd/NttbarEvents

            Fu_ud = Nuu/(Nuu + Ndd)
            Fd_ud = Ndd/(Nuu + Ndd)

            xqHigher_u = Nuubar_highx1_Events + Nubaru_highx2_Events
            xqbarHigher_u = Nuubar_highx2_Events + Nubaru_highx1_Events 

            xqHigher_d = Nddbar_highx1_Events + Ndbard_highx2_Events
            xqbarHigher_d = Nddbar_highx2_Events + Ndbard_highx1_Events

            Du = (xqHigher_u - xqbarHigher_u)/(xqHigher_u + xqbarHigher_u)
            Dd = (xqHigher_d - xqbarHigher_d)/(xqHigher_d + xqbarHigher_d)

            Ain = (N_Yt_LowerthanY0_Events - N_Ytbar_LowerthanY0_Events)/(N_Yt_LowerthanY0_Events + N_Ytbar_LowerthanY0_Events)
            Aout = (N_Yt_HigherthanY0_Events - N_Ytbar_HigherthanY0_Events)/(N_Yt_HigherthanY0_Events + N_Ytbar_HigherthanY0_Events)

            Results[era]['Y0_variations'][str(0.1*(i+1))][region] = {
                    'Fu' : Fu,
                    'Fd' : Fd,
                    'Fu_ud' : Fu_ud,
                    'Fd_ud' : Fd_ud,
                    'Du' : Du,
                    'Dd' : Dd,
                    'FuDu' : Fu*Du,
                    'FdDd' : Fd*Dd,
                    'Ain' : Ain,
                    'Aout': Aout
            } 

save(Results, os.path.join(outputDir,'extractorOutput.coffea'))