
import pandas as pd
import numpy as np
df = pd.read_excel(io='1-s2.0-S2405471219302649-mmc2.xlsx', sheet_name='3N2D-Non-Competitive')
df = df.set_index('ID')

interesting_topID = [58,66,74,82,2127,2135,2137,2184,2191,2189,2202,3111,3116,3118,3135,3140,3268,3270,3656,3658,3678,3695,3700,3702,3719,3724,3852,3854,4262,4276,4295,4460,4465,4482,4487,4613,4615,4629,5399,5401,5411,5424,5418,5417,5426,5529,5456,5449,5458, 5592, 5597, 5601, 5630, 5642, 5644, 5648, 5654, 5659, 5664]

robustness_type = ['Total Robustness', 'Topological Robustness','Extracellular Robustness',  'Intracellular Robustness' ]


for robustness in robustness_type:
    symmetric = df.loc[interesting_topID, robustness]
# symmetric = df.loc[interesting_topID, 'Extracellular Robustness']
# symmetric = df.loc[interesting_topID, 'Intracellular Robustness']
    symmetric = symmetric.loc[(symmetric != 0)]

    non_symmetric = df[robustness].drop(interesting_topID, axis=0)
    non_symmetric = non_symmetric.loc[(non_symmetric != 0)]
    print(robustness)
    print(np.mean(symmetric), np.mean(non_symmetric))
