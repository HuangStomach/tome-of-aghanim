import gc
import time
import numpy as np
import torch
import torch.nn as nn
import logging
from logging import handlers

import metric
import dataset
from models.ae import AutoEncoder
from models.loss import SONLoss, WeightMSELoss

device = 'cuda:0' if torch.cuda.is_available() else 'cpu'

def train(trainData, testData, mask, logger, tag='train'):
    SR = trainData.drug_A
    SR_test = testData.drug_A

    drug_edge, _ = trainData.edge(SR, trainData.params['sim_threshold'])
    drug_edge_test, _ = testData.edge(SR_test, trainData.params['sim_threshold'])

    eye_R = torch.eye(trainData.rnum).float().to(device)
    SR = torch.from_numpy(SR).float().to(device)

    RPI = torch.from_numpy(trainData.rpi).float().to(device)
    RDI = torch.from_numpy(trainData.rdi).float().to(device)
    RPI_test = testData.rpi[mask]
    RDI_test = testData.rdi[mask]

    print("Initialling model...")
    AE = AutoEncoder(
        4096, [trainData.pnum, 1024], [trainData.dnum, 2048],
        2048, [trainData.pnum, 1024], [trainData.dnum, 2048]
    ).to(device)
    optimizer = torch.optim.Adam(AE.parameters(), lr=trainData.params['lr'], weight_decay=trainData.params['wd'])
    son_loss = SONLoss(10)
    mse_loss_p = WeightMSELoss(trainData.params['loss_p_weight'])
    mse_loss_d = WeightMSELoss(trainData.params['loss_d_weight'])

    print("Starting {}...".format(tag))
    for epoch in range(trainData.params['epoch']):
        RPI_hat, SR_hat_1, RDI_hat, SR_hat_2 = AE(
            trainData.drug_x1, trainData.drug_x2, trainData.drug_x3,
            drug_edge
        )

        loss1 = mse_loss_p(RPI_hat, RPI) + trainData.params['a1'] * son_loss(SR_hat_1, SR, eye_R)
        loss2 = mse_loss_d(RDI_hat, RDI) + trainData.params['a2'] * son_loss(SR_hat_2, SR, eye_R)

        loss = loss1 + trainData.params['loss_weight'] * loss2
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        if (epoch + 1) % 10 == 0:
            try:
                RPI_hat_test, _, RDI_hat_test, _ = AE(
                    testData.drug_x1, testData.drug_x2, testData.drug_x3,
                    drug_edge_test
                )
                
                RPI_hat_test = RPI_hat_test.detach().cpu().numpy()[mask]
                RDI_hat_test = RDI_hat_test.detach().cpu().numpy()[mask]
                mp = testData.metric(RPI_test, RPI_hat_test)
                md = testData.metric(RDI_test, RDI_hat_test)

                info = 'Epoch: {} loss: {:.6f}, pauc: {:.6f}, paupr: {:.6f}, paupr_m: {:.6f}, dauc: {:.6f}, daupr: {:.6f}, daupr_m: {:.6f}'.format(
                    epoch, loss.item(), mp[0], mp[1], mp[2], md[0], md[1], md[2]
                )
                logger.info(info)
            except Exception as e:
                print('error', e, RPI_hat_test)

    np.savetxt('output/{}_RPI_hat.txt'.format(tag), RPI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RPI.txt'.format(tag), RPI.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI_hat.txt'.format(tag), RDI_hat.detach().cpu().numpy(), fmt='%f')
    np.savetxt('output/{}_RDI.txt'.format(tag), RDI.detach().cpu().numpy(), fmt='%f')
    torch.save(AE.state_dict(), 'output/{}_model.pt'.format(tag))

if __name__=='__main__':
    while True:
        print("[0] train")
        print("[1] metric")
        print("[2] exit")
        str_in = input("Plz select the opt: ");
        if not str_in.isdigit(): continue

        index = int(str_in)
        if index == 0:
            filename = './output/{}.log'.format(
                time.strftime("%Y%m%d_%H%M%S", time.localtime())
            )
            logger = logging.getLogger(filename)
            logger.setLevel(logging.DEBUG)

            sh = logging.StreamHandler()
            logger.addHandler(sh)
            fh = handlers.RotatingFileHandler(filename=filename)
            logger.addHandler(fh)

            trainData = dataset.Dataset()
            splits = trainData.splits()
            testData = dataset.Dataset()
            testData.init()
            splits = [
                [13, 248, 197, 152, 47, 228, 523, 695, 603, 657, 407, 515, 355, 287, 170, 281, 252, 125, 519, 351, 159, 364, 422, 65, 458, 678, 599, 538, 42, 405, 123, 385, 522, 528, 349, 644, 408, 529, 613, 670, 189, 676, 38, 102, 257, 21, 358, 576, 23, 276, 285, 681, 473, 275, 520, 388, 294, 153, 120, 48, 347, 687, 650, 593, 268, 486, 369, 524, 490, 463, 632],
                [667, 607, 567, 149, 58, 688, 361, 148, 451, 239, 601, 70, 227, 472, 677, 450, 470, 352, 204, 693, 295, 264, 625, 662, 89, 299, 580, 704, 350, 69, 585, 512, 390, 346, 203, 672, 438, 439, 641, 6, 28, 324, 202, 309, 468, 3, 645, 611, 375, 646, 318, 222, 304, 595, 471, 37, 502, 332, 436, 414, 633, 605, 258, 126, 124, 168, 393, 88, 56, 138, 649],
                [200, 372, 63, 652, 64, 527, 119, 525, 686, 340, 560, 627, 572, 271, 130, 245, 22, 703, 391, 401, 377, 409, 327, 578, 462, 201, 628, 291, 115, 10, 459, 371, 698, 133, 384, 59, 354, 480, 206, 629, 146, 392, 570, 244, 394, 289, 237, 692, 298, 684, 325, 419, 551, 230, 403, 701, 208, 671, 542, 196, 506, 251, 26, 55, 543, 273, 53, 292, 110, 399, 598],
                [185, 260, 689, 5, 83, 193, 497, 160, 420, 446, 99, 362, 11, 288, 647, 216, 91, 404, 27, 457, 270, 192, 215, 335, 653, 571, 182, 353, 430, 283, 508, 112, 444, 320, 305, 306, 111, 329, 651, 397, 338, 12, 0, 212, 71, 175, 379, 122, 433, 503, 549, 442, 421, 432, 93, 232, 658, 547, 184, 417, 272, 172, 635, 75, 428, 82, 489, 660, 144, 297, 142],
                [566, 240, 396, 475, 441, 4, 573, 568, 505, 18, 386, 415, 35, 466, 634, 666, 562, 707, 621, 221, 187, 514, 166, 214, 331, 151, 286, 277, 511, 105, 114, 225, 205, 532, 263, 697, 156, 582, 418, 213, 643, 178, 509, 150, 236, 103, 631, 606, 195, 217, 135, 323, 608, 614, 307, 17, 664, 501, 104, 233, 108, 141, 274, 43, 73, 555, 492, 499, 194, 706, 19],
                [249, 481, 2, 255, 504, 209, 176, 656, 284, 622, 638, 477, 33, 36, 174, 16, 494, 183, 51, 339, 296, 54, 162, 229, 360, 234, 541, 496, 493, 147, 705, 132, 445, 456, 545, 96, 448, 235, 669, 412, 589, 685, 565, 169, 238, 380, 673, 167, 553, 413, 482, 435, 163, 210, 46, 558, 406, 491, 654, 682, 539, 378, 207, 700, 211, 665, 68, 507, 680, 609, 303],
                [588, 161, 121, 561, 620, 1, 290, 586, 510, 134, 690, 544, 431, 311, 365, 259, 550, 366, 173, 154, 79, 427, 612, 92, 443, 302, 45, 694, 78, 247, 67, 363, 95, 343, 40, 596, 437, 129, 87, 359, 552, 661, 454, 637, 317, 66, 469, 623, 464, 181, 699, 465, 447, 410, 231, 617, 604, 107, 322, 241, 579, 98, 356, 679, 50, 518, 72, 106, 224, 190, 49],
                [618, 461, 556, 597, 84, 254, 226, 663, 316, 85, 702, 191, 591, 137, 74, 219, 581, 62, 370, 128, 269, 81, 564, 188, 220, 416, 577, 218, 86, 616, 474, 402, 639, 140, 90, 34, 467, 488, 20, 668, 675, 312, 155, 411, 374, 387, 348, 301, 453, 25, 590, 250, 113, 696, 674, 180, 530, 584, 143, 267, 485, 535, 242, 300, 630, 32, 314, 521, 569, 478, 615],
                [326, 39, 9, 118, 395, 265, 7, 546, 594, 373, 344, 602, 266, 60, 310, 333, 479, 80, 554, 116, 61, 164, 636, 171, 655, 495, 308, 315, 533, 368, 101, 282, 583, 336, 592, 426, 44, 157, 536, 127, 198, 429, 376, 94, 460, 425, 440, 642, 516, 253, 513, 452, 548, 659, 600, 557, 243, 14, 624, 15, 531, 76, 424, 575, 328, 367, 145, 640, 626, 199],
                [139, 455, 383, 179, 177, 223, 382, 313, 186, 165, 261, 334, 29, 540, 526, 256, 483, 500, 30, 648, 487, 279, 691, 398, 24, 389, 619, 262, 610, 57, 476, 100, 321, 158, 574, 559, 537, 449, 136, 319, 246, 77, 280, 52, 131, 8, 345, 400, 31, 534, 563, 517, 423, 498, 341, 293, 357, 587, 337, 683, 434, 330, 117, 41, 109, 484, 342, 381, 97, 278],
            ]
            logger.info(trainData.params)
            for i in range(10):
                trainData.init(mask_drugs=splits[i])

                train(trainData, testData, splits[i], logger, i)
                np.savetxt('output/{}_masks.txt'.format(i), splits[i], fmt='%d')

            del logger
            gc.collect()
        elif index == 1:
            metric.run()
        elif index == 2:
            quit()
        else: continue
