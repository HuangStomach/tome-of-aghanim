import torch
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

from datasets.base import Base

class LRSSL(Base):
    inited = False
    base = './data/LRSSL/'
    params = {
        'epoch': 1000, 'lr': 8e-05, 'wd': 0,
        'sim_threshold': 0.7, 'loss_p_weight': 0.993, 'loss_d_weight': 0.994, 'loss_weight': 0.002,
        'a1': 1e-08, 'a2': 1e-08,
        
        'dropout': 0, 'graph_dropout': 0,
        'fr_dim': [4096, [1426, 1024], [681, 256], 4447],
        'fd_dim': [2048, [1426, 1024], [681, 256], 4447]
    }
    path = {
        # 'drugs': base + 'drug.txt',
        'drug_sim': base + 'drug_sim.txt',
        'drug_smiles': base + 'drug_smiles.csv',

        'drug_ecfps': base + 'drug_ecfps12.txt',
        'drug_go': base + 'drug_target_go_mat.txt', # 763*4447
        'rpi': base + 'drug_target_domain_mat.txt', # 763*1426
        'rdi': base + 'drug_dis_mat.txt', # 763*681
    }

    def drugs(self):
        return np.loadtxt(self.path['rdi'], dtype=str, delimiter='\t')[1:, 0]

    def init(self, mask_drugs=None):
        self.mask_drugs = mask_drugs
        self.rpi = self.mask(self.data('rpi', skip=True))
        self.rdi = self.mask(self.data('rdi', skip=True))
        # self.rri = self.mask(self.data('rri'))

        drug_fps = self.mask(self.data('drug_ecfps', delimiter=','))
        drug_go = self.mask(self.data('drug_go', skip=True))
        self.drug_A = self.mask(self.mask(
            self.data('drug_sim', dtype=float, delimiter=',')
        ).T)

        self.drug_x1 = torch.from_numpy(drug_fps).float().to(self.device)
        self.drug_x2 = torch.from_numpy(np.matmul(self.drug_A, self.rpi)).float().to(self.device)
        self.drug_x3 = torch.from_numpy(np.matmul(self.drug_A, self.rdi)).float().to(self.device)
        self.drug_x4 = torch.from_numpy(drug_go).float().to(self.device)
        self.drug_edge = self.edge(self.drug_A, self.params['sim_threshold'])[0]

        self.rnum = self.rpi.shape[0]
        self.pnum = self.rpi.shape[1]
        self.dnum = self.rdi.shape[1]

        self.inited = True
    
    def prepare(self):
        seqs = []
        radius = 6
        length = 4096

        drugs = np.loadtxt(self.path['drug_smiles'], delimiter=',', dtype=str, comments=None)
        for drug in drugs:
            try:
                name, smiles = drug
                mol = Chem.MolFromSmiles(smiles)
                seqs.append(AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=length).ToList())
            except Exception as e:
                print(drug, e)

        np.savetxt(self.path['drug_ecfps'], seqs, fmt='%s', delimiter=',')

    def data(self, name, dtype=int, delimiter='\t', skip=False):
        if hasattr(self, '_' + name):
            return getattr(self, '_' + name)()

        if name not in self.path: return []

        return pd.read_csv(self.path[name], header=0, index_col=0, sep=delimiter).to_numpy(dtype) \
        if skip else \
        np.loadtxt(self.path[name], dtype=dtype, delimiter=delimiter)

    def splits(self) -> list:
        return [
            [263,157,98,117,85,324,247,617,550,754,727,470,534,666,158,603,188,700,
            279,551,386,714,76,298,251,255,758,437,62,55,682,261,50,262,471,698,
            678,583,227,496,241,491,33,619,204,439,568,271,503,123,442,579,675,424,
            106,354,276,495,448,649,694,131,37,99,326,259,197,191,186,651,445,115,
            742,337,595,633],
            [208,647,235,71,546,78,516,428,192,281,549,405,566,201,417,94,48,10,
            226,741,245,475,541,351,151,18,476,733,601,711,340,730,453,660,625,650,
            214,330,462,185,222,487,193,257,381,72,715,533,539,143,220,368,699,520,
            672,147,290,703,177,53,738,333,119,596,515,183,19,504,557,478,697,458,
            118,506,29,364,223],
            [725,542,108,736,565,300,422,310,587,284,295,544,14,691,3,457,607,122,
            73,301,100,592,688,207,349,74,407,591,124,556,82,2,161,657,441,679,
            311,659,684,170,280,685,371,389,153,451,127,497,181,215,167,182,318,145,
            360,438,627,466,614,605,334,332,648,610,744,426,421,706,128,555,210,187,
            589,631,569,652,423],
            [315,275,383,409,588,391,111,656,266,419,645,629,644,609,359,377,294,636,
            517,307,511,411,292,171,571,194,105,141,611,125,373,581,746,705,492,450,
            32,305,572,623,205,8,553,498,536,384,138,670,224,573,26,51,248,40,
            403,578,180,112,563,358,17,582,570,38,213,522,528,577,398,165,584,168,
            190,342,704,745,269],
            [597,120,388,502,594,616,134,68,372,56,283,483,322,404,479,734,732,28,
            641,303,121,567,669,42,260,344,530,639,387,345,413,665,113,196,130,328,
            268,239,327,356,521,237,724,172,447,240,761,274,289,446,288,600,159,378,
            317,459,341,695,509,232,524,630,365,654,140,144,36,348,238,628,104,83,
            526,749,399,400],
            [156,743,61,86,527,309,142,231,580,39,661,692,427,467,312,396,69,397,
            64,708,264,35,200,66,277,667,720,218,731,366,347,545,252,464,518,15,
            456,137,308,233,70,287,31,465,707,149,299,135,668,225,148,564,228,599,
            278,59,683,12,425,762,30,24,443,532,543,139,674,548,723,4,416,490,
            429,126,712,375],
            [243,500,460,272,673,505,174,632,519,92,635,634,91,523,510,316,331,394,
            444,658,206,313,560,265,484,390,202,408,671,590,593,494,561,216,689,65,
            109,97,750,435,480,9,175,355,103,253,339,297,199,677,612,402,34,412,
            757,164,254,481,463,52,323,740,370,493,713,67,585,430,747,574,559,488,
            166,701,5,606],
            [622,719,110,229,89,90,702,538,96,49,540,75,562,709,293,45,410,598,
            267,80,646,0,304,155,41,753,680,270,21,178,513,343,726,418,604,415,
            382,531,486,95,748,643,759,474,325,690,338,132,6,162,244,615,306,721,
            379,7,602,249,468,575,393,653,525,477,676,336,319,499,473,449,681,739,
            621,461,385,512],
            [234,146,184,114,537,735,23,47,87,693,529,751,662,586,63,60,54,129,
            25,357,242,637,58,211,320,376,452,433,362,250,107,282,755,455,102,43,
            363,395,314,716,198,195,718,663,219,369,286,346,81,507,88,285,406,436,
            273,291,626,554,212,361,489,454,687,613,756,729,352,640,27,472,79,173,
            552,176,163,760],
            [209,329,302,547,230,380,13,432,160,717,664,133,485,514,321,414,508,620,
            353,77,655,576,350,1,642,46,203,93,558,501,431,217,420,728,11,150,
            335,469,374,179,256,737,44,686,236,84,608,624,638,246,16,696,57,440,
            367,136,392,221,618,722,752,482,434,296,535,152,20,189,22,710,154,401,
            169,101,116,258]
        ]
