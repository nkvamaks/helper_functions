import re
import numpy as np


# DSIR sparse-21 position weight matrix and spectrum-21 dictionnary for motifs.
# Scores an ANTISENSE strand (19-mer duplex + 2 nt 3'-overhang).
# An accurate and interpretable model for siRNA efficacy prediction.
# BMC Bioinformatics. 2006 Nov 30; 7:520.
# Jean-Philippe Vert, Nicolas Foveau, Christian Lajaunie, Yves Vandenbrouck.
# doi: 10.1186/1471-2105-7-520
DSIR_spectrum21 = {'CC': -0.003932307, 'CG': -0.002849292, 'GG': -0.01153317, 'UA': 0.005858523, 'UC': 0.002756727,
                   'AAA': -0.007093154, 'AAC': 0.01706775, 'AAG': 0.001577039, 'ACA': 0.0163371, 'ACG': -0.01024546,
                   'ACU': 0.01057124, 'AGA': -0.002840863, 'AGC': 0.001306096, 'AGU': 0.0003593278, 'AUC': 0.006090182,
                   'AUG': -0.002284597, 'CCC': -0.02771415, 'CGA': -0.005897083, 'CGG': -0.0004771073, 'CUC': -0.01012959,
                   'GAU': -0.006603347, 'GCC': -0.02225377, 'GCG': -0.01653585, 'GGC': -0.02075344, 'GGG': -0.009166488,
                   'GUC': 0.002456783, 'GUU': -0.01122454, 'UAG': 0.01687206, 'UAU': 0.0009991506, 'UCG': 0.02439596,
                   'UCU': 0.01743714, 'UGA': 0.008405265, 'UGU': 0.004634205, 'UUG': 0.01635677, 'UUU': -0.0007407731}
# PWM: 21 positions x 4 nucleotides (A, C, G, U)
DSIR_sparse21 = np.array([[0.1034911, 0, -0.01902813, 0.1721253],
                          [0, -0.05682008, -0.05498681, 0.007134196],
                          [0.007138161, -0.01702346, 0, 0],
                          [0.01370094, -0.003464999, 0, 0],
                          [0, -0.003201931, 0, 0],
                          [-0.01757542, 0, 0, 0.008817762],
                          [0.002083753, -0.04129517, -0.002207137, 0],
                          [-0.004060712, 0.0001684693, -0.001775484, 0],
                          [-0.009667309, 0.005009017, 0, 0],
                          [0.02766501, 0, 0.009232002, -0.009436056],
                          [0, 0, 0, 0],
                          [0, 0.001366547, 0, 0],
                          [0.008966958, 0, -0.002907738, 0],
                          [0.0002572483, 0, -0.03721179, 0],
                          [-0.01430617, 0.004745051, 0, -0.009492091],
                          [0.002485228, 0, 0, -0.01307864],
                          [0, 0.006607245, -0.0007780886, 0],
                          [-0.04872041, 0.02892875, 0, -0.01385716],
                          [-0.0538126, 0.07225799, 0.05193115, 0],
                          [-0.01070109, -0.02676219, 0, 0.005234203],
                          [0, -0.05839426, 0.01453252, 0]])

# Optimized parameters for the micro effect. PWM for the SENSE strand of 19mer siRNA.
# Nucleic Acids Research, 2007, Vol. 35, No. 4, e27
# 19 positions x 4 nucleotides (A, C, G, U)
micro_katoh = np.array([[-1.555, 3.136, 3.492, -3.155],
                        [-2.237, -0.067, 0.364, -1.856],
                        [1.19, -0.808, -1.226, -0.341],
                        [-0.179, -0.136, -1.362, 0.07],
                        [-0.154, 0.698, -0.454, -2.258],
                        [3.034, -4.081, 0.387, -1.407],
                        [3.038, -2.676, 0.637, -1.44],
                        [2.341, -0.317, -0.735, -2.697],
                        [-2.028, 1.935, -0.938, 2.113],
                        [-1.455, -0.126, -0.787, 2.521],
                        [-0.938, -1.119, 0.614, -0.136],
                        [1.112, -1.279, -0.177, -1.519],
                        [2.601, -1.296, -2.095, 0.2],
                        [2.092, -1.388, 1.037, -2.258],
                        [1.797, -1.801, 1.962, -1.356],
                        [2.758, -1.491, -2.188, 2.128],
                        [0.915, -1.119, -1.09, -1.295],
                        [1.971, -0.53, -2.853, 2.482],
                        [5.31, -2.627, -2.785, 1.354]])
# Pa:Pc:Pg:Pu, SENSE, 19mer
macro_katoh = np.array([6.073, 0.304, 1.265, 5.111])


# PWM for the SENSE strand of 21mer siRNA (19mer duplex + 2nt 3'-overhang).
# Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities.
# M. Ichihara, Y. Murakumo, A. Masuda, T. Matsuura, N. Asai, M. Jijiwa, et al.
# Nucleic Acids Res 2007 Vol. 35 Issue 18 Pages e123. DOI: 10.1093/nar/gkm699
# 21 positions x 4 nucleotides (A, C, G, U)
sBiopredsi_pwm = np.array([[0.122562328943018, 0.33706889764487, -0.493682137955993, 0.034050911368105],
                           [0.199770705864379, 0.0940785052867202, -0.214309709924039, -0.0795395012270602],
                           [-0.191261671166253, 0.394463650967988, 0.471460725387953, -0.674662705189688],
                           [-0.000809566252917729, 0.054510925854472, 0.331933668008531, -0.385635027610085],
                           [0.099264124877664, -0.14751335282012, 0.0353863772920574, 0.0128628506503986],
                           [-0.0485723732993568, -0.0950757876526501, -0.0197787153749778, 0.163426876326985],
                           [0.0197995457567793, -0.0640742009529765, 0.0517871012927189, -0.0075124460965217],
                           [0.236584489707767, -0.392892634338025, 0.00508243649394087, 0.151225708136317],
                           [0.0912473059279591, -0.153846043872856, -0.115891090136529, 0.178489828081426],
                           [0.0461984054908072, -0.0726712088915075, -0.0302677450148171, 0.0567405484155174],
                           [0.0958521152433313, -0.120910870539342, -0.143234372142315, 0.168293127438326],
                           [-0.120734580166361, -0.0171956422019929, -0.160479455766972, 0.298409678135326],
                           [0.152434813120487, -0.120812260990701, 0.0349544804811741, -0.0665770326109601],
                           [0.135123334724672, -0.162962644130062, 0.0208038865888445, 0.0070354228165455],
                           [0.209519491263256, -0.0367328948253503, -0.392476935726244, 0.219690339288338],
                           [0.240908931114548, -0.0415245389715574, -0.0491277555195498, -0.150256636623441],
                           [0.116613628474749, -0.0815556466264104, -0.125498606880402, 0.0904406250320634],
                           [0.0786355957967573, -0.1054036260082, -0.137586003470488, 0.164354033681931],
                           [0.107775984404941, -0.078168890746038, -0.21841866273526, 0.188811569076357],
                           [0.466991470143507, -0.416985564546215, -0.368410735768436, 0.318404830171144],
                           [1.1624724676963, -0.865173403085021, -0.652772391206742, 0.355473326595463]])


# PWM for the ANTISENSE strand of 19mer duplex siRNA.
# Thermodynamic instability of siRNA duplex is a prerequisite for dependable prediction of siRNA activities.
# M. Ichihara, Y. Murakumo, A. Masuda, T. Matsuura, N. Asai, M. Jijiwa, et al.
# Nucleic Acids Res 2007 Vol. 35 Issue 18 Pages e123. DOI: 10.1093/nar/gkm699
# 19 positions x 4 nucleotides (A, C, G, U)
iScore_pwm = np.array([[6.74, -3.92, -6.85, 14.43],
                       [5.96, -1.55, -1.03, 7.03],
                       [5.16, -0.07, 1.54, 3.79],
                       [4.34, 1.03, 1.94, 3.1],
                       [3.19, 1.55, 1.39, 4.28],
                       [1.54, 1.61, 2.3, 4.96],
                       [4.79, -1.46, 1.73, 5.35],
                       [2.99, 2.29, 1.37, 3.76],
                       [2.49, 2.91, 1.44, 3.56],
                       [5.75, 1.19, 2.3, 1.16],
                       [3.36, 1.57, 1.28, 4.2],
                       [3.19, 2.37, 1.48, 3.38],
                       [4.15, 1.65, 0.85, 3.76],
                       [4.41, 2.57, -1.22, 4.64],
                       [2.43, 3.17, 2.06, 2.75],
                       [4.12, 2.24, 1.57, 2.47],
                       [3, 2.66, 1.12, 3.63],
                       [-1.13, 6.08, 3.11, 2.35],
                       [-4.02, 8.02, 6.39, 0.02]])


# Obtained by fit the sequence composition to scores from Horizon Discovery (Dharmacon), ANTISENSE
# PWM: 19 positions x 4 nucleotides (A, C, G, U), built in March 2025
DHARM_pwm = np.array([[ 6.57831431, -4.57608743, -3.58393026, 12.46139671],
                      [ 6.59078986, -1.03081689, -1.84931963,  7.16904   ],
                      [ 4.74823135, -0.57809292,  0.53845075,  6.17110415],
                      [ 1.16403921,  5.10222477,  0.03280701,  4.58062234],
                      [ 0.31774508,  4.54404767,  0.91744062,  5.10045996],
                      [ 0.61069343,  7.34230323,  0.6266646 ,  2.30003207],
                      [ 5.96625827, -4.43507956,  1.44984133,  7.89867329],
                      [ 3.2112489 ,  4.89904944,  0.95821142,  1.81118356],
                      [ 2.86874577,  2.60337925,  2.85449805,  2.55307025],
                      [ 7.38008918, -0.83092216,  2.5847085 ,  1.74581781],
                      [ 1.68711445,  3.63670132,  3.61322799,  1.94264957],
                      [ 0.34167108,  4.79629925,  0.90358063,  4.83814238],
                      [ 2.54581414,  2.60342   , -2.78438796,  8.51484715],
                      [ 4.62805863,  2.94828532, -3.26616529,  6.56951467],
                      [ 1.71303246,  1.17701633, -0.22683024,  8.21647477],
                      [ 2.5808466 , -0.55030554,  3.7395481 ,  5.10960417],
                      [ 3.72014753,  4.27111451, -3.63026765,  6.51869894],
                      [ 0.33960832,  4.52627313,  3.15845428,  2.8553576 ],
                      [-0.50375483,  6.22059019,  4.24604641,  0.91681156]])


def one_hot_encode(seq):
    base_to_idx = {'A': 0, 'C': 1, 'G': 2, 'U': 3}
    oh = np.zeros((len(seq), len(base_to_idx)), dtype=float)
    for i, base in enumerate(seq):
        oh[i, base_to_idx[base]] = 1.0
    return oh

def score_pwm(guide_seq: str, pwm) -> float:
    return np.sum(one_hot_encode(guide_seq) * pwm)

def count_motif_overlap(text: str, motif: str) -> int:
    count = 0
    start = 0
    while True:
        pos = text.find(motif, start)
        if pos == -1:
            break
        count += 1
        start = pos + 1
    return count


# Predict with algorythm from DSIR site http://biodev.cea.fr/DSIR/
# Scores an ANTISENSE 21-mer strand (19-mer duplex + 2 nt 3'-overhang).
# Related paper: An accurate and interpretable model for siRNA efficacy prediction
# J. P. Vert, N. Foveau, C. Lajaunie and Y. Vandenbrouck
# BMC Bioinformatics 2006 Vol. 7 Pages 520. DOI: 10.1186/1471-2105-7-520
def DSIR_score(anti_seq21):
    score_sparse = score_pwm(anti_seq21, DSIR_sparse21)
    score_spectrum = sum(count_motif_overlap(anti_seq21[:19], motif) * value for motif, value in DSIR_spectrum21.items())
    return (score_sparse + score_spectrum + 0.6938215) * 100


# Predict with algorythm from paper. Nucleic Acids Research, 2007, Vol. 35, No. 4, e27
# Scores a SENSE 19-mer.
def katoh_score(sense_seq19, pwm=micro_katoh):
    nt_count = np.array([sense_seq19.count(nt) for nt in 'ACGU'])
    macro_k = nt_count @ macro_katoh
    micro_k = score_pwm(sense_seq19, pwm)
    return macro_k + micro_k


# Predict with algorythm from paper. Nucleic Acids Res 2007 Vol. 35 Issue 18 Pages e123. DOI: 10.1093/nar/gkm699
# Scores a SENSE 21-mer strand (19-mer duplex + 2 nt 3'-overhang).
def sBiopredsi_score(sense_seq21, pwm=sBiopredsi_pwm):
    score = score_pwm(sense_seq21, pwm)
    score = 1/(1+np.exp(-score-0.204977962358907))
    adj_score = (-2.27627506037018 + 4.14301528286201 * score) * 0.168285137478663 + 0.581913382218149
    return adj_score


# Predict with algorythm from paper. Nucleic Acids Res 2007 Vol. 35 Issue 18 Pages e123. DOI: 10.1093/nar/gkm699
# Scores an ANTISENSE 19-mer strand.
def iScore_score(anti_seq19, pwm=iScore_pwm):
    return score_pwm(anti_seq19, pwm)


# Prediction of Horizon Discovery, pwm of March 2025
# Scores an ANTISENSE 19-mer strand.
def dharmacon_score(anti_seq19: str, pwm=DHARM_pwm) -> float:
    pwm_score = score_pwm(anti_seq19, pwm)
    pattern = re.compile(r'(AAAA|CCCC|UUUU|GGGG)')
    matches = list(pattern.finditer(anti_seq19))
    unique_homopolymers = {match.group(0) for match in matches}
    penalty_score = len(unique_homopolymers) * 8.5
    total_score = pwm_score - penalty_score
    return total_score