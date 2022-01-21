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
        'epoch': 1000, 'lr': 1e-04, 'wd': 0,
        'sim_threshold': 0.8, 'loss_p_weight': 0.993, 'loss_d_weight': 0.994, 'loss_weight': 0.001,
        'a1': 0.000000001, 'a2': 0.000000001,
        
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

    # def splits(self) -> list:
    #     return [[ 13, 391, 377, 432, 322, 313, 608,  50, 177, 161, 387, 117, 460,
    #    230, 365, 505, 596,  39, 181, 105, 684,  95, 139, 119, 343,  46,
    #     79, 218, 380, 321, 375, 549, 382, 435,  21, 288, 499, 448, 114,
    #    597, 530, 299, 376, 338,  67, 219, 390, 709, 714, 127,  89, 733,
    #    452, 551, 357, 134, 109, 532, 635, 362, 195, 624, 615, 261, 531,
    #    254, 336,  32, 478, 472, 264, 407, 147, 332, 439,  72, 574], [437, 687, 467, 752,  83,  90, 744, 400, 262, 639, 379,  27, 447,
    #    184, 617, 564, 168, 356, 412, 273, 225, 211, 165, 651, 183,  76,
    #    523, 223, 274, 221, 590, 692, 458, 717, 192,  10, 330, 371, 197,
    #    420, 672,  38, 628, 631, 304, 611, 133, 138,  58, 470, 237, 757,
    #    132,  41, 594, 501, 535, 152, 681,  33, 490, 258,  19,  23, 372,
    #    374, 413, 746, 473,  63, 621, 438,  30, 289, 331,  55, 270], [159, 232, 193, 652, 679, 227, 298, 658, 189, 500, 285, 164, 202,
    #    194,  43, 228, 110, 661, 169, 712, 118, 440, 598, 415, 582, 335,
    #    204, 441, 384, 125, 685, 210, 563, 108, 185, 676, 737,  16, 297,
    #    431, 269, 418, 708, 200, 593, 485, 643,   4,  31,  15,  14, 481,
    #    525, 358, 236, 565, 393, 707,  73, 340, 128, 754, 308,   8, 493,
    #    436, 411, 483,  54,  86, 451, 129,  17, 399,  75, 655, 536], [111, 233, 552, 645,   6, 403, 318, 682, 622, 693, 762, 241, 281,
    #    516, 512, 610,  28, 587, 514, 729, 724, 324, 678, 345, 726, 348,
    #    580, 723, 507, 670,  74, 627, 106, 492, 494, 539, 648, 588, 520,
    #    576, 251, 715, 443, 327, 220, 632, 488, 725, 278, 302, 188, 347,
    #    528, 419, 445, 602,  56,  20, 167, 245, 208, 646, 410, 718, 326,
    #    721, 429, 498, 276, 179,  69, 115, 242, 217,  25, 389], [149, 585, 612, 703,  60, 392,  40, 126, 346, 284, 203, 508, 120,
    #     88,  81, 414, 713, 710, 625, 207, 334, 283, 287, 567,  99, 619,
    #    626, 517, 699, 300, 605, 677, 442,  97, 248, 742, 484,  18, 282,
    #    130, 654, 761, 614,  36, 454, 637, 727, 142, 354, 657,  48, 292,
    #    235, 312, 397, 166,   0, 546, 137, 231, 355, 173, 591, 145,  62,
    #    474, 427,  29, 339, 342, 408, 466, 265, 561,  94, 595], [ 84, 706, 421, 479, 170, 267, 301, 636, 290, 562, 446, 212, 151,
    #    121, 249, 541, 163, 756, 659, 370, 306, 751, 406, 463, 502, 349,
    #    518, 303, 323, 522, 461, 430, 257, 351, 405, 700, 579, 175,  37,
    #    259,  42, 540, 238,   1, 705, 720,  11, 491, 601, 743, 394, 255,
    #     71, 641, 280, 381, 697, 246, 680, 606, 277, 575, 239, 388, 496,
    #    423, 198, 350, 749,   7, 457, 275, 450, 113, 519,  22], [730, 583, 104, 630, 745, 360, 468, 455, 401, 510, 731, 416,  45,
    #    464, 543,  77, 521,  51, 741, 586, 694, 144, 172, 296,  91, 373,
    #    640, 325, 650, 623, 504, 758, 537, 215, 471, 252, 305,   3, 475,
    #    253, 341, 644, 286, 409, 569, 136, 526, 266, 616,  64, 487, 704,
    #    229, 250, 314, 459, 256, 307, 666, 489, 554, 158, 667, 329, 150,
    #    613, 247,   2, 503, 395, 124, 486, 642, 538, 589, 383], [  5,  68, 755,  65,  26, 647, 722, 534, 135,  53, 739, 609, 469,
    #    141, 396, 544, 570, 691,  85, 316, 673, 555, 716, 102, 449, 271,
    #    386, 174, 112, 740, 206, 295, 243, 344, 728, 153, 260, 367, 101,
    #    607, 361, 199, 620, 688, 568, 123, 572, 550, 353, 263, 557, 515,
    #     78, 663,  98,  52, 315, 476, 131, 696, 385, 558, 497, 402, 456,
    #    226, 702, 182,  24, 618, 753,  34, 317, 190, 671,  82], [584, 428, 732, 213, 599, 154, 234, 337, 100, 196, 191, 103, 205,
    #    698, 311, 592, 434, 146, 669, 660, 529, 224, 577,  59, 433, 186,
    #     87, 424, 156, 533, 378, 363, 291, 668, 604, 664, 581, 320, 665,
    #    178, 524,  70, 222, 352, 629, 369, 683,  12, 760, 735, 333, 542,
    #     80, 750, 495, 122, 527, 453, 547, 310, 559, 160,  47, 548, 553,
    #    140, 157, 662, 214,  49, 279, 560, 690, 465, 293, 107], [748, 417, 759, 603,  66, 319, 462, 656, 634, 426, 719, 272, 398,
    #    738, 695, 359, 747, 209, 686, 600, 422,   9,  96, 368, 513, 701,
    #     44, 328, 511, 689, 638, 143, 509, 573, 180, 364, 309, 244, 480,
    #    162, 578, 482, 736, 171, 268, 176, 294, 566, 404, 116, 649, 425,
    #    653,  93, 545,  57, 734, 477, 201, 571,  35, 633, 506, 556, 675,
    #    240,  92,  61, 216, 366, 444, 674, 187, 155, 711, 148]]

